# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""ε, conductivity, and the fail-closed check on material types.

The global subpixel knob and the notebook hook live in :mod:`openem.scene.subpixel` (the functions
are re-exported here under the same names).
"""

from __future__ import annotations

import hashlib
import pathlib

import numpy as np
import tidy3d as td
from tidy3d import config

from openem import knobs
from openem.scene import poles as _poles
from openem.scene._util import COORD_KEYS, E_KEYS, yee_coords
from openem.scene.cache import atomic_write, cache_root, hash_data_arrays
# ANISOTROPIC / PEC_LIKE / is_pec_like / medium_parts are definitions that depend only on the
# material; they are implemented in poles.py and re-exported here under their old names.
from openem.scene.poles import (  # noqa: F401
    ANISOTROPIC, PEC_LIKE, is_pec_like, medium_parts, sigma_from_eps,
)
# The subpixel knob and the notebook hook live in scene/subpixel.py; this is only a re-export (the
# startup hook .ipython_openem/startup/00-openem.py calls media.install_notebook_staircase).
from openem.scene.subpixel import (  # noqa: F401
    cloud_emulation, install_notebook_staircase, local_subpixel, require_local_subpixel,
    set_local_subpixel,
)

#: On-disk cache directory for ε. Override it with ``OPENEM_EPS_CACHE``; set it to an empty string
#: to turn the cache off (``None``).
#: (``pathlib.Path("")`` used to be ``PosixPath('.')``, so an empty string actually wrote the cache
#: into the current directory.)
#:
#: Why cache at all: subpixel ε goes through ``tidy3d_extras`` as "write a temp file, call C++,
#: read it back" (`TempDir` + `get_material_yee_cpp`), not pure memory. A small box takes 9–18 s,
#: while CavityFOM has 7.9M cells and OptimizedL3 has 53M, and one scene needs 3 calls (6 when
#: lossy). This is the **only** bottleneck when sweeping every case.
_EPS_CACHE_RAW = knobs.env("EPS_CACHE")
# Unset: cache/epsilon next to the repo (the default has to be computed, so KNOBS registers None).
if _EPS_CACHE_RAW is None:
    _EPS_CACHE_RAW = str(cache_root() / "epsilon")
EPS_CACHE_DIR: pathlib.Path | None = pathlib.Path(_EPS_CACHE_RAW) if _EPS_CACHE_RAW else None

def _eps_cache_key(sim: td.Simulation, coord_key: str, freq: float | None) -> str:
    """Only put **things that affect ε** into the key.

    Structures, background medium, subpixel settings, grid boundary coordinates, component and
    frequency, plus the tidy3d version (changing version turns something already verified back into
    an unknown, so the cache must not be reused across versions). The spatial arrays of a structure
    go through :func:`openem.scene.cache.hash_data_arrays` separately, because ``json`` drops them
    silently.

    Sources and monitors are deliberately **not** included: they do not affect ε, and validation
    often changes only a monitor (checking ε for CavityFOM, for instance, means taking out the
    monitor that carries apodization).
    """
    h = hashlib.sha256()
    h.update(td.__version__.encode())
    h.update(repr(coord_key).encode())
    h.update(repr(freq).encode())
    h.update(str(bool(sim.subpixel)).encode())
    # The global knob changes the result of sim.epsilon too: with it on, even a sim with
    # subpixel=False goes through the extras C++ engine (PEC marks 432 cells instead of the 368 of
    # point sampling). Leaving it out of the key lets the two contaminate each other.
    # None must not be folded into False: tidy3d reads None as "use extras if it is there", a
    # different convention from the staircasing of an explicit False.
    h.update(repr(config.simulation.use_local_subpixel).encode())
    h.update(sim.subpixel.json().encode() if hasattr(sim.subpixel, "json") else b"-")
    h.update(sim.medium.json().encode())
    hash_data_arrays(sim.medium, h)
    for s in sim.structures:
        h.update(s.json().encode())
        hash_data_arrays(s, h)          # json drops the data, so feed the arrays separately
    for ax in "xyz":
        h.update(np.ascontiguousarray(
            np.asarray(getattr(sim.grid.boundaries, ax), dtype=np.float64)).tobytes())
    return h.hexdigest()[:32]


#: tidy3d expresses conductivity in S/µm, SI in S/m. The scalar path inverts it through
#: ``sigma_from_eps`` and comes out in SI naturally; the tensor path reads the ``conductivity``
#: field directly and has to multiply itself.
SIGMA_T3D_TO_SI = 1.0e6


def has_tensor(sim: td.Simulation) -> bool:
    """Whether the scene holds a ``FullyAnisotropicMedium`` (on a structure or as background)."""
    if type(sim.medium) is td.FullyAnisotropicMedium:
        return True
    return any(type(st.medium) is td.FullyAnisotropicMedium
               for st in media_structures(sim))


def tensor_free_sim(sim: td.Simulation) -> td.Simulation:
    """The copy used for scalar ε/σ queries: tensor media taken out of the scene, grid pinned.

    - A ``FullyAnisotropicMedium`` on a **structure** is deleted outright (the same idea as
      ``pec_free_sim``): the material of a tensor cell is expressed entirely by the entry table of
      :func:`tensor_entries` (staircased, the same convention as the reference implementation,
      since tidy3d itself does no subpixel averaging for fully anisotropic media, as its docstring
      says), and the scalar arrays on those cells are overwritten by pass-through and never reach
      the physics.
    - A tensor **background** is replaced by a lossless diagonal ``AnisotropicMedium`` (the
      diagonal of ε), so the downstream consumers of the scalar arrays (PML profile, shutoff energy
      weighting, the 1D background check of TFSF) get the diagonal approximation of the tensor
      rather than vacuum. σ is deliberately left out: all tensor σ goes down the tensor path, and
      carrying it here would instead trigger the "lossy + subpixel" fail closed.
      Note that the previous point, "the scalar arrays on those cells are overwritten by
      pass-through", does **not** hold for a tensor background: with a tensor background every cell
      in the domain is a true tensor entry (``bg=0`` in :func:`tensor_entries`), and the diagonal
      copy here only feeds the PML/shutoff/TFSF checks, never the physics.
    """
    kw: dict = {"grid_spec": td.GridSpec.from_grid(sim.grid)}
    if type(sim.medium) is td.FullyAnisotropicMedium:
        e = np.asarray(sim.medium.permittivity, dtype=np.float64)
        kw["medium"] = td.AnisotropicMedium(
            xx=td.Medium(permittivity=float(e[0, 0])),
            yy=td.Medium(permittivity=float(e[1, 1])),
            zz=td.Medium(permittivity=float(e[2, 2])))
    keep = [s for s in sim.structures
            if type(s.medium) is not td.FullyAnisotropicMedium]
    if len(keep) != len(sim.structures):
        kw["structures"] = keep
    sim2 = sim.updated_copy(**kw)
    _assert_same_grid(sim, sim2, "swapping out tensor media + from_grid")
    return sim2


def _dilate3(mask: np.ndarray, periodic: tuple[bool, bool, bool]) -> np.ndarray:
    """3×3×3 box dilation. Periodic axes wrap around, non-periodic axes stop at the edge (there are
    no cells outside the domain that would need pass-through)."""
    out = mask
    for ax in range(3):
        m = out
        lo = np.roll(m, 1, axis=ax)
        hi = np.roll(m, -1, axis=ax)
        if not periodic[ax]:
            sl_lo = [slice(None)] * 3
            sl_lo[ax] = slice(0, 1)
            lo[tuple(sl_lo)] = False
            sl_hi = [slice(None)] * 3
            sl_hi[ax] = slice(-1, None)
            hi[tuple(sl_hi)] = False
        out = m | lo | hi
    return out


def _tensor_media_index(sim: td.Simulation, structs: list) -> tuple[list, int, list[int]]:
    """``(tensor_media, bg, struct_idx)``: the list of tensor media, the index of the background
    (-2 = not a tensor), and the tensor medium index of each structure (-2 = not a tensor)."""
    tensor_media: list = []
    struct_idx: list[int] = []          # tensor medium index per structure, -2 = not a tensor
    if type(sim.medium) is td.FullyAnisotropicMedium:
        tensor_media.append(sim.medium)
        bg = 0
    else:
        bg = -2
    for st in structs:
        if type(st.medium) is td.FullyAnisotropicMedium:
            tensor_media.append(st.medium)
            struct_idx.append(len(tensor_media) - 1)
        else:
            struct_idx.append(-2)
    return tensor_media, bg, struct_idx


def _governing_masks(sim: td.Simulation, structs: list, struct_idx: list[int],
                     bg: int, shape) -> list[np.ndarray]:
    """One ``gov (nx,ny,nz) int32`` array per component: which tensor medium governs that Yee point
    (-2 = not a tensor). Structures follow tidy3d's priority: later ones override earlier ones."""
    masks = []
    for ax_name in "xyz":
        xs, ys, zs = yee_coords(sim, ax_name, meters=False)     # tidy3d's µm coordinates
        gov = np.full(shape, bg, dtype=np.int32)
        for st, si in zip(structs, struct_idx):
            inside = st.geometry.inside_meshgrid(xs, ys, zs)
            if inside.shape != shape:
                raise AssertionError(
                    f"inside_meshgrid shape {inside.shape} != grid {shape}")
            gov[inside] = si
        masks.append(gov)
    return masks


def _true_entries(c: int, gov: np.ndarray, tensor_media: list) -> tuple:
    """The true tensor entries ``(comp, cell, eps, sigma)`` of component ``c`` (σ already in
    S/m)."""
    flat = np.flatnonzero(gov >= 0)
    gv = gov.ravel()[flat]
    n = flat.size
    et = np.empty((n, 3, 3))
    st_ = np.empty((n, 3, 3))
    for mi, med in enumerate(tensor_media):
        sel = gv == mi
        if not sel.any():
            continue
        et[sel] = np.asarray(med.permittivity, dtype=np.float64)
        st_[sel] = np.asarray(med.conductivity, dtype=np.float64) * SIGMA_T3D_TO_SI
    return np.full(n, c, dtype=np.int32), flat.astype(np.int64), et, st_


def _halo_entries(c: int, gov: np.ndarray, dilated: np.ndarray,
                  eps_c: np.ndarray, sigma_c: np.ndarray | None) -> tuple:
    """The diagonal halo entries ``(comp, cell, eps, sigma)`` of component ``c``: the diagonal
    element takes the ε/σ the scalar path has in that cell, the other diagonal entries are 1."""
    halo = np.flatnonzero(dilated & (gov < 0))
    nh = halo.size
    eh = np.tile(np.eye(3), (nh, 1, 1))
    eh[:, c, c] = eps_c.ravel()[halo]
    sh = np.zeros((nh, 3, 3))
    if sigma_c is not None:
        sh[:, c, c] = sigma_c.ravel()[halo]
    return np.full(nh, c, dtype=np.int32), halo.astype(np.int64), eh, sh


def tensor_entries(sim: td.Simulation, grid, eps: dict[str, np.ndarray],
                   sigma: dict[str, np.ndarray | None]):
    """Entry table for ``FullyAnisotropicMedium`` (:class:`model.TensorEps`).

    **Staircasing**: if the Yee sample point of component ``c`` falls inside a tensor medium, that
    (comp, cell) gets one true tensor entry. This follows the same convention as the reference
    implementation, since tidy3d states explicitly that fully anisotropic media get no subpixel
    averaging (docstring in medium.py), so nothing extra is conceded here. Structures follow
    tidy3d's priority: later ones override earlier ones.

    The 4-neighbor average of a true tensor entry requires its neighbors to go through
    pass-through as well, so a ring of **diagonal halo** entries is laid around the tensor mask:
    the diagonal element takes the ε/σ the scalar path has in that cell and the other diagonal
    entries are 1, which makes the update equation of these entries degenerate verbatim into the
    scalar path (top of ``kernels/tensor.cu``).

    Returns:
        ``TensorEps``, or None when the scene has no tensor medium.
    """
    from openem.model import TensorEps

    structs = media_structures(sim)
    tensor_media, bg, struct_idx = _tensor_media_index(sim, structs)
    if not tensor_media:
        return None

    shape = grid.shape
    masks = _governing_masks(sim, structs, struct_idx, bg, shape)

    periodic = tuple(a.is_periodic for a in grid.axes)
    union = (masks[0] >= 0) | (masks[1] >= 0) | (masks[2] >= 0)
    dilated = _dilate3(union, periodic)

    comp_l, cell_l, eps_l, sig_l = [], [], [], []
    for c in range(3):
        key = E_KEYS[c]
        gov = masks[c]
        for entry in (_true_entries(c, gov, tensor_media),
                      _halo_entries(c, gov, dilated, eps[key], sigma[key])):
            comp_l.append(entry[0])
            cell_l.append(entry[1])
            eps_l.append(entry[2])
            sig_l.append(entry[3])

    return TensorEps(
        comp=np.concatenate(comp_l),
        cell=np.concatenate(cell_l),
        eps=np.concatenate(eps_l),
        sigma=np.concatenate(sig_l),
    )


#: Tolerance for deciding "pure conductivity". The σ of a pure conductor is **strictly**
#: independent of frequency (measured on BeerLambert: the values inverted at the two frequencies
#: are bitwise identical), so this can be kept extremely tight; any real dispersion goes far past
#: it.
#: The lenient verdict for "lossy + subpixel". Subpixel averaging only stops Im(ε) being
#: proportional to 1/ω in **the thin layer at the interface**, and the other cells stay just as
#: clean: measured on a dielectric sphere (2 THz, σ=50 S/m) the drift had median 0 and worst
#: 2.9e-4, with 16/405224 = 0.004% of cells over tolerance. When both the fraction of
#: over-tolerance cells and the drift stay inside these two lines, take the σ inverted at the
#: reference frequency, carry on and print the cell count; past them, still fail closed.
SIGMA_SUBPIXEL_SOFT_FRAC = 0.02
SIGMA_SUBPIXEL_SOFT_DRIFT = 1e-2

SIGMA_FREQ_RTOL = 1e-6   # loosened from 1e-9: the C++ read-back path of subpixel
# averaging introduces ~1e-8 of roundoff in the σ inversion (not real dispersion; real dispersion
# drifts by far more than 1e-6). It only lets things through, it changes no exported number.

#: Threshold for calling a cell PEC. Sampling ε on PEC, tidy3d always returns
#: ``constants.pec_val`` (-1e8, a real scalar): on the small-sphere scene, measured under both
#: conventions (subpixel on and off), that is the **only** negative value that shows up, it cannot
#: be confused with any physical ε, and there is no "averaged intermediate negative value". Half of
#: it is used as the threshold.
PEC_EPS_THRESHOLD = 0.5 * float(td.constants.pec_val)


def is_lossy_metal(medium) -> bool:
    """A medium handled as PEC that is not an actual ``PECMedium`` (``LossyMetalMedium``)."""
    return type(medium) is not td.PECMedium and is_pec_like(medium)


def has_pec(sim: td.Simulation) -> bool:
    """Whether the scene has a PEC structure (or a LossyMetal handled as PEC). A PEC background is
    rejected in :func:`refuse_unsupported`."""
    return any(is_pec_like(st.medium) for st in sim.structures)


def pec_point_sample_sim(sim: td.Simulation) -> td.Simulation:
    """The copy used for PEC mask queries: ``subpixel=False``, plain Yee point sampling.

    For PEC we only do **staircasing**: if the Yee sample point of that component falls inside PEC,
    E in that cell is identically 0 (ca=cb=0 in the solver). The mask has to be taken under the
    **point-sampling** convention. Measured on a sphere of radius 0.22 µm, the three conventions
    mark different numbers of cells (``PECConformal`` 304, ``Staircasing`` under the subpixel
    engine 432, point sampling 368), and only point sampling agrees cell by cell with "the Yee
    point is inside the geometry" (locked down by tests/test_pec.py). Staircasing instead of the
    conformal treatment the reference implementation uses by default is a **known accuracy
    concession**.

    This copy only produces the mask: ε is still queried on the original sim (the subpixel
    averaging between media must not be lost), and dt also comes from the original sim (the
    ``timestep_reduction`` of PECConformal has to be kept to stay step-for-step aligned with the
    reference implementation). subpixel does not affect the grid, which is asserted here as well
    (fail closed).
    """
    sts = list(sim.structures)
    if any(is_lossy_metal(s.medium) for s in sts):
        # LossyMetal → PECMedium: only point sampling returns pec_val, which is what enters the
        # mask. Changing a medium makes AutoGrid move the grid, so pin it to the original grid.
        sts = [s.updated_copy(medium=td.PECMedium()) if is_lossy_metal(s.medium) else s
               for s in sts]
        sim2 = sim.updated_copy(subpixel=False, structures=sts,
                                grid_spec=td.GridSpec.from_grid(sim.grid))
        _assert_same_grid(sim, sim2, "turning subpixel off + LossyMetal→PEC + from_grid")
        return sim2
    sim2 = sim.updated_copy(subpixel=False)
    _assert_same_grid(sim, sim2, "turning subpixel off")
    return sim2


def pec_free_sim(sim: td.Simulation) -> td.Simulation:
    """The copy used for ε queries: PEC structures removed, grid pinned to the original grid.

    PEC has no usable ε: on cells touching a PEC face, the subpixel engine returns a mixed value
    contaminated with pec_val (measured on a 624³ sphere, 31,000 cells per component with
    ε∈(1,2), which amounts to wrapping the sphere in a fake dielectric shell). The staircasing
    semantics are "E≡0 on masked cells, every other cell takes its ε from the material layout
    **without any PEC**", so the query goes straight to a copy with the PEC structures deleted: the
    subpixel averaging between media is kept exactly as it was and PEC contamination is zero. The
    grid is pinned with ``GridSpec.from_grid``, since deleting a structure moves the grid in an
    AutoGrid scene.

    Structure priority is not lost: the mask comes from point sampling on the **complete** sim, so
    a dielectric structure that comes after the PEC and covers it keeps those points out of the
    mask, while deleting the PEC here only exposes the medium underneath it.
    """
    keep = [s for s in sim.structures if not is_pec_like(s.medium)]
    sim2 = sim.updated_copy(structures=keep,
                            grid_spec=td.GridSpec.from_grid(sim.grid))
    _assert_same_grid(sim, sim2, "deleting PEC structures + from_grid")
    return sim2


def _assert_same_grid(a: td.Simulation, b: td.Simulation, what: str) -> None:
    for axn in "xyz":
        if not np.array_equal(np.asarray(getattr(a.grid.boundaries, axn)),
                              np.asarray(getattr(b.grid.boundaries, axn))):
            raise AssertionError(f"{what} must not move the grid, but the grid changed")


def epsilon_complex(
    sim: td.Simulation, coord_key: str, freq: float | None = None
) -> np.ndarray:
    """Complex ε_r at the Yee component positions. ``coord_key='Ex'`` samples the xx component of
    the tensor at the Ex position.

    With ``subpixel=False`` this is point sampling (the three components are identical); with
    ``subpixel=True`` plus extras it is the diagonal of the averaged tensor (the three components
    differ on interface cells). The off-diagonal elements are identically 0 for axis-aligned
    geometry, so only the diagonal is taken.

    The result is hashed over "the inputs that affect ε" and cached in :data:`EPS_CACHE_DIR` (see
    the note there).

    Args:
        freq: required for lossy or dispersive materials. Passing None means infinite frequency,
            where σ contributes 0 (``Im ε = σ/(ω ε₀) → 0``), which loses the loss entirely.
    """
    key = _eps_cache_key(sim, coord_key, freq)
    path = EPS_CACHE_DIR / f"{key}.npz" if EPS_CACHE_DIR is not None else None
    if path is not None and path.exists():
        with np.load(path) as z:
            return z["eps"]

    box = td.Box(center=(0, 0, 0), size=(td.inf, td.inf, td.inf))
    arr = sim.epsilon(box, coord_key=coord_key, freq=freq)
    out = np.asarray(arr.values, dtype=np.complex128)

    if EPS_CACHE_DIR is not None:
        # Store compressed: ε has large runs of the same value, measured at about 1/40 the size.
        # Write a temp file first, then rename atomically (cache.atomic_write).
        # (Pass a file handle, not a path: np.savez appends a suffix to a path not ending in .npz.)
        EPS_CACHE_DIR.mkdir(parents=True, exist_ok=True)

        def _write(tmp):
            with open(tmp, "wb") as f:
                np.savez_compressed(f, eps=out)
        atomic_write(path, _write)
    return out




def eps_and_sigma(
    sim: td.Simulation, coord_key: str, freq: float, freq2: float
) -> tuple[np.ndarray, np.ndarray | None]:
    """``(real part of ε_r, σ)``; σ is None when there is no loss.

    σ is inverted once at each of two frequencies and the two are required to agree: the σ of a
    pure conductor does not depend on frequency, that of a dispersive material does. This test is
    more reliable than a whitelist by medium type, but it **only holds when no subpixel averaging
    was done**: the averaged complex ε no longer satisfies ``Im ∝ 1/ω``. The lenient verdict for
    "lossy + subpixel" is in the block below (:data:`SIGMA_SUBPIXEL_SOFT_FRAC` /
    :data:`SIGMA_SUBPIXEL_SOFT_DRIFT`); past it, fail closed.

    Raises:
        NotImplementedError: the σ inverted at the two frequencies disagree, which means the
            material is dispersive.
    """
    eps = epsilon_complex(sim, coord_key, freq)
    sigma = sigma_from_eps(eps, freq)
    scale = float(np.max(np.abs(sigma)))
    if scale == 0.0:
        return np.real(eps), None

    sigma2 = sigma_from_eps(epsilon_complex(sim, coord_key, freq2), freq2)
    drift = float(np.max(np.abs(sigma - sigma2))) / scale
    if drift > SIGMA_FREQ_RTOL:
        # Triage on the evidence. **Median drift 0 with only a thin layer of cells over tolerance**
        # is the signature of subpixel averaging: the σ of those cells is an intermediate value
        # between pure material values, and since the averaged complex ε no longer satisfies
        # Im ∝ 1/ω (see the lenient verdict above), inverting σ there was never defined.
        # Real dispersion makes **every cell** of that material drift.
        # The triage matters: the fallback in nb.backend that rebuilds with subpixel=False on any
        # subpixel-related error matches on whether the message contains "subpixel", so blaming
        # dispersion would make it miss the case for nothing. (Measured on FullyAnisotropic:
        # median 0, 0.02–0.09% of cells over tolerance, and |σ| on those cells spread continuously
        # over 139..9140.)
        frac_bad = float(np.mean(np.abs(sigma - sigma2) / scale > SIGMA_FREQ_RTOL))
        if bool(sim.subpixel) and frac_bad < 0.5:
            if (frac_bad <= SIGMA_SUBPIXEL_SOFT_FRAC
                    and drift <= SIGMA_SUBPIXEL_SOFT_DRIFT):
                print(f"[openem] {coord_key}: σ drifts with frequency on "
                      f"{frac_bad * 100:.3f}% of the cells (worst {drift:.3e}), the familiar "
                      "lossy material + subpixel problem on interface cells; both the fraction "
                      f"and the magnitude are inside the lenient lines "
                      f"({SIGMA_SUBPIXEL_SOFT_FRAC:.0%} / {SIGMA_SUBPIXEL_SOFT_DRIFT:g}), so "
                      "take the σ at the reference frequency and carry on.")
                return np.real(eps), sigma
            raise NotImplementedError(
                f"{coord_key}: the inverted σ varies with frequency (relative drift "
                f"{drift:.3e} > {SIGMA_FREQ_RTOL:g}), but only {frac_bad * 100:.3f}% of the "
                "cells are over tolerance. This is **a lossy material plus subpixel averaging**: "
                "the averaged complex ε no longer satisfies Im ∝ 1/ω, so inverting σ is undefined "
                "on those boundary cells. Rebuild with subpixel=False, or wait until the "
                "convention for lossy subpixel is settled."
            )
        raise NotImplementedError(
            f"{coord_key}: the inverted σ varies with frequency ({freq / 1e12:.2f} → "
            f"{freq2 / 1e12:.2f} THz, relative drift {drift:.3e} > {SIGMA_FREQ_RTOL:g}, with "
            f"{frac_bad * 100:.1f}% of the cells over tolerance). This is a dispersive material, "
            "not a pure conductor. Dispersion has to wait for ADE."
        )
    return np.real(eps), sigma      # the two frequencies agree, σ is clean


def is_freq_independent(sim: td.Simulation) -> bool:
    """True only when every medium can be **proven** to have a frequency-independent ε; anything
    uncertain returns False.

    For a non-dispersive medium with zero conductivity, ``ε = permittivity`` does not depend on ω;
    with conductivity, ``ε(ω) = ε_r + iσ/(ωε₀)`` does. Dispersive media (Sellmeier, Lorentz,
    PoleResidue and so on) are all subclasses of ``DispersiveMedium`` and are excluded outright.

    This is only a **performance** knob (the ε monitor readout in nb.backend and the ε map
    perturbation in nb/shapegrad both use it to sample ε at a single frequency): returning False is
    merely slow, never wrong, while returning True has to genuinely hold, so any unrecognized type,
    any non-zero conductivity and any exception returns False (locked down by
    tests/test_fasteps.py). nb.backend and shapegrad used to carry one copy each (the shapegrad
    copy only recognized ``Medium``, so PEC, CustomMedium and non-dispersive anisotropic media went
    through multi-frequency sampling); they are now merged into this one.
    """
    try:
        from tidy3d.components.medium import DispersiveMedium
    except Exception:
        return False
    OK = ("Medium", "PECMedium", "CustomMedium", "AnisotropicMedium",
          "CustomAnisotropicMedium")
    try:
        meds = list(getattr(sim.scene, "mediums", None) or [])
        if not meds:
            return False
        for m in meds:
            # an anisotropic medium has to be checked component by component
            subs = list((getattr(m, "components", None) or {}).values()) or [m]
            for s in subs:
                if isinstance(s, DispersiveMedium):
                    return False
                if type(s).__name__ not in OK:
                    return False
                c = getattr(s, "conductivity", None)
                if c is None:
                    continue
                v = np.asarray(getattr(c, "values", c))
                if v.size and np.any(v != 0.0):
                    return False
        return True
    except Exception:
        return False


def monitor_freqs(sim: td.Simulation) -> list[np.ndarray]:
    """The ``freqs`` of each monitor (float64, at least 1-D), in monitor order; monitors without
    frequencies (time-domain monitors, empty arrays) are not in the list. The four places that
    collect monitor frequencies all share this one."""
    out = []
    for mon in sim.monitors:
        f = getattr(mon, "freqs", None)
        if f is None:
            continue
        arr = np.atleast_1d(np.asarray(f, dtype=np.float64))
        if arr.size:
            out.append(arr)
    return out


def _second_freq(sim: td.Simulation, freq: float) -> float:
    """The second probe frequency used to detect dispersion: the one in the monitor band furthest
    from ``freq``.

    With a single-frequency monitor, fall back to ``0.9*freq``. As long as the two frequencies
    differ, a pure conductor gives exactly the same σ at both.
    """
    for arr in monitor_freqs(sim):
        far = float(arr[int(np.argmax(np.abs(arr - freq)))])
        if abs(far - freq) > 1e-9 * freq:
            return far
    return 0.9 * freq


def reference_freq(sim: td.Simulation) -> float | None:
    """The reference frequency for sampling ε: the median of the monitor frequencies if there is
    one, otherwise the freq0 of the source."""
    for f in monitor_freqs(sim):
        return float(f[len(f) // 2])
    if sim.sources:
        f0 = getattr(sim.sources[0].source_time, "freq0", None)
        if f0:
            return float(f0)
        # For a source without freq0 (CustomSourceTime and the like): the center of its spectrum
        st = sim.sources[0].source_time
        try:
            fr = st.frequency_range(num_fwidth=4.0) if hasattr(st, "frequency_range") else None
            if fr and all(np.isfinite(fr)) and fr[1] > 0:
                return 0.5 * float(fr[0] + fr[1])
        except (AttributeError, TypeError, ValueError):
            pass
    # Last resort: midpoint of the Simulation's own frequency range (tidy3d derives that range from
    # the sources and the monitors)
    try:
        fr = sim.frequency_range
        if fr and all(np.isfinite(fr)) and fr[1] > 0:
            return 0.5 * float(fr[0] + fr[1])
    except (AttributeError, TypeError, ValueError):
        pass
    return None


def media_structures(sim: td.Simulation) -> list:
    """The structure list used when enumerating by medium type: ``Medium2D`` has already been
    replaced by its equivalent volumetric structure.

    ``sim.volumetric_structures`` is the conversion the tidy3d client ships: ``Medium2D`` →
    ``AnisotropicMediumFromMedium2D`` (the residues and conductivity of the in-plane components are
    scaled by 1/dl and weighted-averaged with the neighboring medium, the normal component takes
    the neighbor on the + side), with the geometry snapped to ``sim.grid``. Both paths of
    ``sim.epsilon`` (point sampling and extras subpixel) use these structures internally, so the
    whitelist, the weight columns and the pole table have to look at the same ones to stay
    consistent with the ε query. With no Medium2D present, ``sim.structures`` is returned
    unchanged.

    **The grid must come from the original sim**: AutoGrid computes a different grid for the
    converted structures, so ``sim`` must never be replaced wholesale by the converted version.
    """
    return [s for s in sim.volumetric_structures if not is_lossy_metal(s.medium)]


def refuse_unsupported(sim: td.Simulation) -> None:
    """Three classes of material that have to fail closed.

    1. **Time-varying / nonlinear**: ``modulation_spec`` and ``nonlinear_spec`` are separate fields
       on ``Medium`` and **do not show up in the values ``sim.epsilon`` returns** (which only give
       one instant at one linearization point), so looking at ε cannot find them and the fields
       have to be inspected explicitly. ``modulation_spec`` has a supported subset (scalar CW
       modulation of ε, see scene/modulation.py); anything outside that subset still raises
       loudly, and ``nonlinear_spec`` is always refused.
    2. **Anything that is not a plain ``Medium``**: decided by **exact type**, not ``isinstance``.
       ``LossyMetalMedium`` inherits from ``Medium`` and carries its own ``.conductivity``, so duck
       typing would let it through, while Tidy3D treats it with a surface impedance boundary
       (Leontovich SIBC) rather than a volumetric σ, and we would then run it with a completely
       wrong model. The same check catches dispersive and various special media along the way.
       ``AnisotropicMedium`` (a diagonal tensor) goes through the same whitelist **component by
       component**; ``FullyAnisotropicMedium`` (full tensor, including the antisymmetric σ of a
       gyrotropic medium) goes down the staircased entry table of :func:`tensor_entries`.
    3. **Lossy + subpixel**: see :func:`eps_and_sigma`. The averaged complex ε no longer satisfies
       ``Im ∝ 1/ω``, and there is no evidence for which convention the reference implementation
       uses for the equivalent σ in the time domain.

    ``Medium2D`` is first converted by :func:`media_structures` into its equivalent volume
    (diagonal anisotropic + PoleResidue) and only then passed through the whitelist here; an
    unrecognized type inside the equivalent volume still raises.
    """
    labels = _labels(sim)
    _refuse_nonlinear_and_modulation(labels)
    _refuse_exotic_types(labels)


def _labels(sim: td.Simulation) -> list[tuple[str, object]]:
    """``[(name, medium)]``: every structure (Medium2D already converted) plus the background."""
    labels = [(f"structures[{i}]", st.medium)
              for i, st in enumerate(media_structures(sim))]
    labels.append(("background medium", sim.medium))
    return labels


def _refuse_nonlinear_and_modulation(labels: list) -> None:
    """Time-varying: the supported subset is let through, everything else raises loudly in
    scene/modulation.py; nonlinear is always refused."""
    from openem.scene import modulation as _modulation

    for name, m in labels:
        _modulation.refuse_unsupported_spec(name, m)
    bad = [f"{name} {type(m).__name__} carries a nonlinear_spec"
           for name, m in labels for p in medium_parts(m)
           if getattr(p, "nonlinear_spec", None) is not None]
    if bad:
        raise NotImplementedError("nonlinear materials are not supported: " + "; ".join(bad[:3]))


#: The exact scalar medium types that are allowed: the non-dispersive Medium, plus the few that
#: reduce to PoleResidue (including spatial Custom dispersion) and CustomMedium. Still decided by
#: **exact type**, not isinstance: ``LossyMetalMedium`` inherits from ``Medium`` and carries its
#: own ``.conductivity``, so duck typing would let it through, while Tidy3D treats it with a
#: surface impedance boundary (Leontovich SIBC) rather than a volumetric σ.
ALLOWED_SCALAR = (td.Medium, *_poles.DISPERSIVE, *_poles.CUSTOM_DISPERSIVE,
                  *(getattr(td, n) for n in ("CustomMedium",) if hasattr(td, n)))


def _refuse_exotic_types(labels: list) -> None:
    """Pass the whitelist :data:`ALLOWED_SCALAR` by exact type. A **structure** may also be a
    PECMedium (staircased, E identically 0), but the background may not be PEC: E would be
    identically 0 over the whole domain, leaving nothing to solve."""
    exotic = []
    for name, m in labels:
        if type(m) in ANISOTROPIC:
            for ax, part in zip(("xx", "yy", "zz"), medium_parts(m)):
                if type(part) not in ALLOWED_SCALAR:
                    exotic.append(
                        f"{name} is an AnisotropicMedium whose {ax} component is "
                        f"{type(part).__name__}")
        elif type(m) is td.FullyAnisotropicMedium:
            continue  # full tensor: staircased entries (tensor_entries), structure or background
        elif is_pec_like(m) and name != "background medium":
            if type(m) is not td.PECMedium:
                print(f"[openem] {name} is a {type(m).__name__}, handled as staircased PEC "
                      "(surface impedance loss ignored, see media.PEC_LIKE)")
            continue  # PEC on a structure: staircased mask
        elif type(m) not in ALLOWED_SCALAR:
            exotic.append(f"{name} is a {type(m).__name__}")
    if exotic:
        raise NotImplementedError(
            "only td.Medium (optionally with conductivity), td.PECMedium on a structure, "
            "dispersive media that reduce to PoleResidue "
            f"({', '.join(t.__name__ for t in _poles.DISPERSIVE)}), "
            "td.FullyAnisotropicMedium (full tensor, staircased), and td.AnisotropicMedium whose "
            "diagonal components are of those types are supported: "
            + "; ".join(exotic[:3])
        )



def reference_freq_pair(sim: td.Simulation) -> tuple[float, float]:
    """``(reference frequency, the second frequency used to detect dispersion)``.

    Raises:
        ValueError: no reference frequency can be obtained, which makes it impossible to tell
            whether the material is lossy.
    """
    freq = reference_freq(sim)
    if freq is None:
        raise ValueError(
            "no reference frequency can be taken from the scene, cannot judge material loss")
    return freq, _second_freq(sim, freq)


def epsilon_and_sigma(
    sim: td.Simulation,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray | None],
           dict[str, np.ndarray] | None]:
    """``(real part of ε, σ, PEC mask)`` at the three Yee component positions.

    With PEC present the query splits in two, taking the best of each (both copies keep the grid
    pinned):

    - **ε** comes from :func:`pec_free_sim` (PEC structures deleted): the subpixel averaging
      between media is kept exactly as it was, and no pec_val or mixed value gets in;
    - **the mask** comes from :func:`pec_point_sample_sim` (point sampling on the complete sim):
      cells with ``ε ≤`` :data:`PEC_EPS_THRESHOLD` are recorded as flat indices for that component,
      and structure priority is handled by tidy3d's own sampling.

    Masked cells are refilled with 1.0 in the ε array and 0 in σ: E is identically 0 on those
    cells so their values never reach the physics, and a constant keeps the downstream consumers
    (the 1D background of TFSF, the PML coefficients, the energy weighting) clean.

    Returns:
        ``(eps, sigma, pec)``; ``pec`` is None when there is no PEC structure.

    Raises:
        ValueError: no reference frequency can be obtained, which makes it impossible to tell
            whether the material is lossy.
        NotImplementedError: a dispersive material, or the "lossy + subpixel" combination whose
            convention is not settled yet.
    """
    freq, freq2 = reference_freq_pair(sim)
    pec_here = has_pec(sim)
    sim_mask = pec_point_sample_sim(sim) if pec_here else None
    sim_eps = pec_free_sim(sim) if pec_here else sim
    eps: dict[str, np.ndarray] = {}
    sigma: dict[str, np.ndarray | None] = {}
    pec: dict[str, np.ndarray] = {}
    for i, key in enumerate(E_KEYS):
        eps[key], sigma[key] = eps_and_sigma(sim_eps, COORD_KEYS[i], freq, freq2)
        if np.any(eps[key] <= PEC_EPS_THRESHOLD):
            raise AssertionError(
                f"eps_{key} still contains pec_val cells; the query convention of the PEC-free "
                "copy is wrong")
        if sim_mask is None:
            continue
        # The mask query must **turn off** the extras global knob: with it on, even a sim with
        # subpixel=False goes through the C++ engine and what gets marked is not the point-sampling
        # convention (432 vs 368, see :func:`pec_point_sample_sim`). The knob is part of the cache
        # key, so the two conventions cannot contaminate each other.
        with local_subpixel(False):
            m_eps = np.real(epsilon_complex(sim_mask, COORD_KEYS[i], freq))
        mask = m_eps <= PEC_EPS_THRESHOLD
        eps[key] = np.where(mask, 1.0, eps[key])
        if sigma[key] is not None:
            sigma[key] = np.where(mask, 0.0, sigma[key])
        pec[key] = np.flatnonzero(mask.ravel()).astype(np.int64)
    return eps, sigma, (pec if pec_here else None)
