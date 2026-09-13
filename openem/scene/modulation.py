# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Extracting ``modulation_spec`` (time-modulated media), and failing closed.

The mathematical convention is **taken from the tidy3d client**,
``components/time_modulation.py``, not guessed:

    δε(r,t) = Re[amp_time(t) · amp_space(r)]
    amp_time = A_t · exp(i·φ_t − i·2π·f·t)
    amp_space = A_r(r) · exp(i·φ_r(r))
    ⇒ δε = A_t·A_r · cos(2πf·t − φ_t − φ_r)

Which cells count as modulated: **those whose Yee sample point is inside the geometry**, i.e. a
point-sampling convention. The basis is a note in the reference notebook: when subpixel averaging
is on it is applied to the passive part of the medium but not to the active modulation part. So
the static eps goes through subpixel averaging as usual while the modulated part does not, and the
mask is decided with ``geometry.inside`` at the Yee points.

The supported subset (the validation set contains exactly one such example; everything else fails
closed):

- the structure's medium is **exactly** ``td.Medium`` with ``conductivity == 0``;
- only ``permittivity`` is modulated (modulating ``conductivity`` is not supported);
- the time part is ``ContinuousWaveTimeModulation``;
- the space part is a scalar ``SpaceModulation`` (``SpatialDataArray`` is not supported);
- exactly one modulated structure, not covered by a later structure;
- no coexistence with background modulation, dispersion, PEC, loss, TFSF or mode sources, or
  dipoles overlapping it.
"""

from __future__ import annotations

import numpy as np
import tidy3d as td

from openem.grid import Grid, ravel_cells, unravel_cells
from openem.model import TimeModulation
from openem.scene import media, poles as _poles
from openem.scene._util import yee_coords


def _spec(medium) -> td.ModulationSpec | None:
    """The ModulationSpec that is **in effect** on a medium; both fields None counts as no
    modulation.
    """
    ms = getattr(medium, "modulation_spec", None)
    if ms is None or not (ms.permittivity is not None or ms.conductivity is not None):
        return None
    return ms


def any_modulated(sim: td.Simulation) -> bool:
    """Whether the scene has any modulation_spec in effect, background and anisotropic components
    included, for deciding which path to take.
    """
    for st in media.media_structures(sim):
        for p in _poles.medium_parts(st.medium):
            if _spec(p) is not None:
                return True
    return any(_spec(p) is not None for p in _poles.medium_parts(sim.medium))


def refuse_unsupported_spec(name: str, medium) -> None:
    """Whether a medium's modulation_spec falls inside the supported subset; raise loudly if not.

    Called from :func:`media.refuse_unsupported`; returns immediately when there is no
    modulation_spec.
    """
    if type(medium) in _poles.ANISOTROPIC:
        for ax, p in zip(("xx", "yy", "zz"), _poles.medium_parts(medium)):
            if _spec(p) is not None:
                raise NotImplementedError(
                    f"{name} is an AnisotropicMedium whose {ax} component carries a "
                    "modulation_spec; anisotropy combined with time modulation is not supported yet")
        return
    ms = _spec(medium)
    if ms is None:
        return
    if name == "background medium":
        raise NotImplementedError(
            "a modulation_spec on the background medium is not supported yet; the validation set "
            "only modulates structures")
    if type(medium) is not td.Medium and type(medium) not in _poles.DISPERSIVE:
        raise NotImplementedError(
            f"{name} is a {type(medium).__name__} carrying a modulation_spec; time modulation "
            "supports non-dispersive td.Medium plus dispersive media that reduce to PoleResidue "
            f"({', '.join(t.__name__ for t in _poles.DISPERSIVE)})")
    # Dispersive media (Sellmeier, Lorentz and so on) have no .conductivity field; their loss lives
    # in the poles and does not go through this sigma path
    if float(getattr(medium, "conductivity", 0.0)) != 0.0:
        raise NotImplementedError(
            f"{name} carries a modulation_spec and conductivity={medium.conductivity}; time "
            "modulation combined with conductivity is not supported yet")
    if ms.conductivity is not None:
        raise NotImplementedError(
            f"{name} modulates conductivity; only modulating permittivity is supported")
    stm = ms.permittivity
    if type(stm.time_modulation) is not td.ContinuousWaveTimeModulation:
        raise NotImplementedError(
            f"the time modulation of {name} is {type(stm.time_modulation).__name__}; only "
            "ContinuousWaveTimeModulation is supported")
    sm = stm.space_modulation
    if type(sm) is not td.SpaceModulation:
        raise NotImplementedError(
            f"the space modulation of {name} is {type(sm).__name__}; only SpaceModulation is supported")
    for f in ("amplitude", "phase"):
        v = getattr(sm, f)
        # Both a scalar and a SpatialDataArray are supported; the latter is linearly interpolated in
        # build at each entry's Yee coordinates. A travelling-wave modulation, for instance, has
        # its phase varying cell by cell along y.
        if not isinstance(v, (int, float)) and not hasattr(v, "interp"):
            raise NotImplementedError(
                f"{name} has space modulation {f} of type {type(v).__name__}; only a scalar or "
                "a SpatialDataArray is supported")


def _sample_space(v, base: float, X, Y, Z, flat):
    """Value of a space-modulated quantity at the selected entries: ``base + v(r)``.

    A scalar ``v`` broadcasts; a ``SpatialDataArray`` is linearly interpolated at each entry's Yee
    coordinates. The coordinates are first clamped into the data range: a Yee point in the
    modulated region can sit half a cell beyond the data coordinates, where extrapolation would
    produce NaN, and clamping at the boundary is physically just "keep using the nearest value".
    """
    if isinstance(v, (int, float)):
        return np.full(flat.size, base + float(v), dtype=np.float64)
    import xarray as xr
    pts = {}
    for key, arr in (("x", X), ("y", Y), ("z", Z)):
        c = np.asarray(v.coords[key], dtype=np.float64)
        pts[key] = xr.DataArray(
            np.clip(arr.ravel()[flat], c.min(), c.max()), dims="u")
    out = np.asarray(v.interp(**pts, method="linear"), dtype=np.float64).ravel()
    if not np.all(np.isfinite(out)):
        raise ValueError("interpolating the space modulation produced a non-finite value; the "
                         "coordinates or the data range are wrong")
    return base + out


def build(sim: td.Simulation, grid: Grid,
          eps: dict[str, np.ndarray],
          sigma: dict[str, np.ndarray | None]) -> TimeModulation | None:
    """Extract a :class:`TimeModulation` from sim; return None when there is no modulation.

    Precondition: :func:`media.refuse_unsupported` has already rejected the unsupported
    combinations, so any modulation_spec reaching here is inside the supported subset.

    Args:
        eps / sigma: the result of :func:`media.epsilon_and_sigma`. The positivity check on eps(t)
            needs the static eps, and sigma is used to reject a lossy modulated cell.
    """
    structures = media.media_structures(sim)
    modulated = [(i, st) for i, st in enumerate(structures)
                 if _spec(st.medium) is not None]
    if not modulated:
        return None
    if len(modulated) != 1:
        raise NotImplementedError(
            f"there are {len(modulated)} modulated structures; only exactly one is supported "
            "(more would first need a convention for overlapping cells and multiple frequencies)")
    idx, st = modulated[0]
    if any(v is not None for v in sigma.values()):
        # The ca = eps^n/eps^{n+1} the solver uses on a modulated cell carries no sigma term, and
        # the convention for loss in the same scene has never been verified. Keep the same line as
        # the guard in solver.run and fail closed.
        raise NotImplementedError("a time-modulated medium in the same scene as a lossy one is not "
                                  "supported yet")

    stm = st.medium.modulation_spec.permittivity
    tm, sm = stm.time_modulation, stm.space_modulation
    freq = float(tm.freq0)
    # amp and phase are computed per entry, since the space modulation may be a SpatialDataArray.
    # The time part is a scalar; the space part multiplies the amplitude and adds to the phase.
    amp_l, phase_l = [], []

    # The Yee sample coordinates must be on the same grid as the eps arrays: sim.grid's Yee points, in um
    shape = grid.shape
    nx, ny, nz = shape
    comp_l, cell_l = [], []
    for c, key in enumerate(("ex", "ey", "ez")):
        xs = tuple(yee_coords(sim, "xyz"[c], meters=False))
        if tuple(v.size for v in xs) != shape:
            raise AssertionError(
                f"the Yee coordinates of E{'xyz'[c]}, {tuple(v.size for v in xs)}, do not match the grid {shape}")
        X, Y, Z = np.meshgrid(*xs, indexing="ij")
        mask = np.asarray(st.geometry.inside(X, Y, Z), dtype=bool)
        # A point covered by a later structure does not belong to this medium (tidy3d gives the
        # later one priority). Such a point means the convention for "what is modulated" would be
        # wrong, so fail closed.
        for j in range(idx + 1, len(structures)):
            over = mask & np.asarray(structures[j].geometry.inside(X, Y, Z), dtype=bool)
            if over.any():
                raise NotImplementedError(
                    f"the modulated structure is covered by structures[{j}] at {int(over.sum())} "
                    f"E{'xyz'[c]} Yee points; the convention for modulation under coverage is undecided")
        flat = np.flatnonzero(mask.ravel())
        amp_c = float(tm.amplitude) * _sample_space(sm.amplitude, 0.0,
                                                    X, Y, Z, flat)
        phase_c = _sample_space(sm.phase, float(tm.phase), X, Y, Z, flat)
        # eps(t) = eps_s + delta must stay positive; the reference notebook notes that eps_inf must
        # be positive
        es = eps[key].ravel()[flat]
        if flat.size and float(np.min(es - np.abs(amp_c))) <= 0.0:
            raise ValueError(
                f"a modulated E{'xyz'[c]} cell has eps_s - |amp| = "
                f"{float(np.min(es - np.abs(amp_c))):.4f} <= 0, so eps(t) would become non-positive "
                "and the scheme unstable")
        comp_l.append(np.full(flat.size, c, dtype=np.int32))
        cell_l.append(flat.astype(np.int32))
        amp_l.append(amp_c)
        phase_l.append(phase_c)

    comp = np.concatenate(comp_l)
    cell = np.concatenate(cell_l)
    if comp.size == 0:
        raise ValueError("a modulation_spec is in effect but the geometry covers no Yee point at "
                         "all; the mask convention may be wrong, so stop rather than continue")
    mod = TimeModulation(
        freq=freq, comp=comp, cell=cell,
        amp=np.concatenate(amp_l),
        phase=np.concatenate(phase_l))
    mod.validate(shape)
    return mod


def check_placement(mod: TimeModulation, shape: tuple[int, int, int], srcs) -> None:
    """Fail closed on every combination of modulation and source placement whose convention is
    undecided.

    ``srcs`` is a :class:`openem.scene.sources.BuiltSources`.

    - **A plane wave injection plane crossing a modulated cell**: the incident table is fixed from
      the static eps, so the injection would silently have the wrong amplitude.
    - **TFSF or mode sources**: the 1D background table and the mode profile are both computed from
      the static medium, and the convention is undecided.
    - **A dipole stencil landing on a modulated cell**: the cb that multiplies the injection
      changes every step, and that behaviour has never been verified.
    """
    if srcs.tfsf_srcs:
        raise NotImplementedError("a time-modulated medium with a TFSF source is not supported yet")
    if srcs.mode_srcs:
        raise NotImplementedError("a time-modulated medium with a mode source is not supported yet")
    _, _, kk = unravel_cells(mod.cell, shape)
    for src in srcs.planes:
        if np.any(kk == src.plane_index):
            raise NotImplementedError(
                f"the plane wave injection plane k={src.plane_index} crosses a modulated cell; the "
                "convention is undecided")
    if srcs.dipoles:
        cells = {int(v) for v in mod.cell}
        for d in srcs.dipoles:
            for (i, j, k) in d.indices:
                if ravel_cells(int(i), int(j), int(k), shape) in cells:
                    raise NotImplementedError(
                        "a dipole stencil lands on a modulated cell; the injection coefficient "
                        "changes every step and the convention is undecided")
