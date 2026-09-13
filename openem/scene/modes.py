# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""ModeMonitor: project the fields on a plane onto guided modes.

This runs on **the same machinery** as the far-field projection (``scene/projection.py``): the
solver only has to compute the phasors on the plane, and the mode decomposition is left to the
local ``ModeSolver`` plus ``ModeSolverData.outer_dot``. Three measurements justify doing it
this way:

* the local ModeSolver's ``n_eff`` is **bitwise identical** to Tidy3D's (1e-14 at worst) across
  several examples and across many frequencies and modes, so the mode profiles can be computed
  locally without any difference in convention
* modes are normalized to **unit power** (``ModeSolverData.flux == 1``), so the absolute
  amplitude rests on something
* ``outer_dot(..., bidirectional=True)`` gives the mode amplitudes in both directions directly

So there is no physics in this file, only wrapping our phasors into a ``td.FieldData`` and
feeding them in.
"""

from __future__ import annotations

import contextlib as _ctx
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
import tidy3d as td

from openem import knobs
from openem.cpml import C_0
from openem.grid import cyclic_axes, normal_axis
from openem.model import MODE_PLANE_SUFFIX
from openem.scene._util import plane_index, sign_of, upstream_h_index
from openem.scene.cache import atomic_write, cache_root, hash_data_arrays
from openem.units import UM

if TYPE_CHECKING:
    from openem.model import ModeSource

#: Name suffix of the expanded FieldMonitor (the value lives in model.MODE_PLANE_SUFFIX, which
#: td_readout refers to as well).
PLANE_SUFFIX = MODE_PLANE_SUFFIX


def is_mode_monitor(mon) -> bool:
    """Whether this is a ``td.ModeMonitor``.

    Tested by class name rather than with ``isinstance``: ``ModeSolverMonitor`` belongs to the
    same family but means something different (it stores the mode profiles themselves, not the
    decomposed amplitudes) and will have to be handled separately later.
    """
    return type(mon).__name__ == "ModeMonitor"


def plane_monitor(mon) -> td.FieldMonitor:
    """ModeMonitor -> a ``FieldMonitor`` recording all six components on its plane.

    ``colocate=True``: the mode overlap integral needs E and H evaluated at the same point.
    Profiles at the Yee positions are the form needed for **injection**, not here.

    Planarity is not checked here; that is left to upstream, since Tidy3D itself rejects a
    volumetric ModeMonitor (pinned down by
    tests/test_modes.py::test_tidy3d_itself_rejects_volume_mode_monitor).
    """
    return td.FieldMonitor(
        center=mon.center, size=mon.size, freqs=mon.freqs,
        name=mon.name + PLANE_SUFFIX, colocate=True)


def mode_solver_data(sim: td.Simulation, mon):
    """Solve for the modes on the monitor plane.

    ``ModeSolver`` solves on **this sim**, so the medium, the grid and the frequencies all match
    the main simulation; that is the precondition for the local solver agreeing with Tidy3D's.
    """
    from tidy3d.plugins.mode import ModeSolver

    solver = ModeSolver(
        simulation=sim,
        plane=td.Box(center=mon.center, size=mon.size),
        mode_spec=mon.mode_spec,
        freqs=tuple(float(f) for f in np.atleast_1d(mon.freqs)),
        colocate=True,
    )
    # The mode basis depends only on (grid + medium + mode spec + frequency) and not on the
    # solution, yet it used to be re-solved on every injection (51 s for a single GroupDelay
    # monitor). It is cached to disk under a hash of the solver; OPENEM_MODE_CACHE=0 turns the
    # cache off.
    return _cached_solve(solver)


@_ctx.contextmanager
def _memo_eps(sim: td.Simulation):
    """In a non-dispersive scene, recomputing the subpixel permittivity at every frequency point
    of a mode solve is redundant (tidy3d_extras starts one gencoeffs subprocess per frequency at
    ~1 s each; the 201-frequency mode monitor of EffectiveIndexApproximation had still not
    finished after 3 h). When every medium is a Medium or a PECMedium,
    SubpixelSimulation._get_epsilon is memoized on (subpixel scheme, input file prefix) with the
    frequency ignored; dispersive media are left as they are."""
    meds = [sim.medium] + [st.medium for st in sim.structures]
    if any(type(m).__name__ not in ("Medium", "PECMedium") for m in meds):
        yield
        return
    try:
        from tidy3d_extras import subpixel as _sp
    except Exception:
        yield
        return
    cls = _sp.SubpixelSimulation
    orig = cls.__dict__.get("_get_epsilon")
    if orig is None:
        yield
        return
    cache: dict = {}
    fn = orig.__func__ if isinstance(orig, (classmethod, staticmethod)) else orig

    def memo(klass, freq, subpixel_scheme, input_file_name_prefix, tmp_path):
        key = (repr(subpixel_scheme), str(input_file_name_prefix))
        if key not in cache:
            cache[key] = tuple(fn(klass, freq, subpixel_scheme, input_file_name_prefix, tmp_path))
        return tuple(cache[key])

    cls._get_epsilon = classmethod(memo)
    try:
        yield
    finally:
        cls._get_epsilon = orig


#: Public name: nb.backend uses it too when solving modes locally (the old name ``_memo_eps``
#: is kept).
memo_eps = _memo_eps


def disk_cached(key: str, solve, load, default_dir=None, miss_note: str = "solve and write out"):
    """The skeleton of the mode cache (shared by :func:`_cached_solve` here and by
    ``_mode_solver_data`` / ``_mode_simulation_data`` in nb.backend; all three used to carry
    their own copy).

    ``OPENEM_MODE_CACHE=0`` calls ``solve`` directly; otherwise it looks for ``{key}.hdf5``
    under ``OPENEM_MODE_CACHE_DIR`` (falling back to ``default_dir``, then to
    ``cache_root/modes``): on a hit it calls ``load(path)``, on a miss it calls ``solve`` and
    writes the result out atomically (``cache.atomic_write``, ``{key}.<pid>.part.hdf5`` then
    renamed). **The key is computed by the caller**: :func:`mode_cache_key` on the scene side,
    ``sha1(_sim_cache_key(...))[:20]`` on the nb side. Both algorithms stay as they are, so
    existing mode caches still hit.

    The mode basis cache is on by default. It grows fast: on 2026-09-06 it piled up to 120 GB in
    a single day and filled the workspace quota, so clear cache/modes first when disk is tight.
    """
    import pathlib
    if knobs.env("MODE_CACHE") == "0":
        return solve()
    cdir = knobs.env("MODE_CACHE_DIR")
    if cdir is None:                 # unset: caller's default, else cache/modes (None in KNOBS)
        cdir = str(default_dir if default_dir is not None else cache_root() / "modes")
    cdir = pathlib.Path(cdir)
    cdir.mkdir(parents=True, exist_ok=True)
    fp = cdir / f"{key}.hdf5"
    if fp.exists():
        print(f"[openem] mode cache hit {key}")
        return load(str(fp))
    print(f"[openem] mode cache miss {key}, {miss_note}")
    data = solve()
    atomic_write(fp, lambda tmp: data.to_file(str(tmp)))     # {key}.<pid>.part.hdf5 then renamed
    return data


def mode_cache_key(solver) -> str:
    """Cache key of the scene-side ``ModeSolver``: ``sha1(solver.json ‖ spatial arrays)[:20]``.

    The json does not contain DataArray values (CustomMedium pixels and so on), so changing only
    the pixels leaves the JSON identical: AutoGrid follows the pixels while the mode basis hits
    the stale cache, and `_place` then reports that the coordinates do not match the grid (seen
    with a 4x4 pixel gradient probe). The arrays are fed into the hash separately, as in
    media._eps_cache_key.
    """
    import hashlib
    _h = hashlib.sha1(solver.json().encode())
    hash_data_arrays(solver.simulation, _h)
    return _h.hexdigest()[:20]


def _load_mode_solver_data(path: str):
    from tidy3d.plugins.mode.mode_solver import ModeSolverData
    return ModeSolverData.from_file(path)


def _cached_solve(solver):
    """``ModeSolver.solve`` with a disk cache, written out under :func:`mode_cache_key`.

    The mode basis depends only on (grid + medium + mode spec + frequency). Three places share
    it: the monitor basis, the ModeSource profile and the broadband mismatch check. Only the
    monitor path used to be cached, so every repeat run re-solved the source's modes from
    scratch (of order 30 minutes on a 1.17e8-cell sim).
    """
    def _solve():
        with _memo_eps(solver.simulation):
            return solver.solve()

    # The mode basis cache is on by default (OPENEM_MODE_CACHE=0 turns it off).
    if knobs.env("MODE_CACHE") != "0":
        try:
            return disk_cached(mode_cache_key(solver), _solve, _load_mode_solver_data)
        except Exception as e:
            print(f"[openem] mode cache unavailable ({type(e).__name__}), solving directly")
    return _solve()


def amplitudes(sim: td.Simulation, mon, data: td.FieldData) -> dict:
    """Fields on a plane -> mode amplitudes in both directions.

    Returns:
        ``{"+": (num_modes, nf), "-": (num_modes, nf)}``, complex.

    ``outer_dot`` itself **does not separate the directions**; the two settings give two
    different overlap integrals (``monitor_data.py:1447``)::

        bidirectional=True   ->  ¼∫(E₀*×H₁ + H₀*×E₁)·dS
        bidirectional=False  ->  ½∫(E₀*×H₁)·dS

    Write the field as ``a₊·(E,H) + a₋·(E,-H)`` (the backward mode flips the tangential H and
    leaves the tangential E alone) and substitute::

        ¼∫(E₀*×H₁ + H₀*×E₁) = a₊
        ½∫(E₀*×H₁)          = a₊ - a₋

    so ``a₋`` follows from the difference of the two. **There is no adjustable factor**: the
    modes are normalized to unit power (``ModeSolverData.flux == 1``) and the fields by
    :func:`projection.normalization`.
    """
    md = mode_solver_data(sim, mon)
    fwd = md.outer_dot(data, bidirectional=True)
    diff = md.outer_dot(data, bidirectional=False)
    dim = "mode_index_0" if "mode_index_0" in fwd.dims else "mode_index"

    def _np(a):
        return np.asarray(a.transpose(dim, "f").data, dtype=np.complex128)

    a_fwd = _np(fwd)
    return {"+": a_fwd, "-": a_fwd - _np(diff)}

#: Sign of the half-cell phase. The incident plane of H sits half a cell **upstream** of the
#: source plane (at ks-½ for propagation along +x), and the mode field there leads the source
#: plane by a phase of β·dl/2. The sign was **fixed by measurement**, the criterion being
#: directionality in a uniform waveguide (the backward amplitude should be ~1e-8). Getting the
#: direction wrong does not diverge, it only lifts the backward amplitude to the 1e-2 range,
#: which is a silent error.
HALF_CELL_SIGN = -1.0


def _half_cell_sign() -> float:
    """:data:`HALF_CELL_SIGN`; ``OPENEM_HALF_CELL_SIGN`` is a measurement back door, read on
    demand and warned about every time a source is built (it used to be read and printed at
    import time, which gave the import side effects)."""
    sign = knobs.env("HALF_CELL_SIGN")
    if sign is None:                 # unset (None in KNOBS: branch on set/unset, not on value)
        return HALF_CELL_SIGN
    import sys as _sys
    print("!! OPENEM_HALF_CELL_SIGN overridden by an environment variable: this is a "
          "measurement back door, the physics it produces cannot be trusted",
          file=_sys.stderr)
    return float(sign)


def _disc_on() -> bool:
    """Whether the mode source is injected with the **fully discrete** (Yee + leapfrog)
    eigenmode rather than the continuous Maxwell mode.

    ``OPENEM_MODESRC_DISC=0`` falls back to the old behavior (solve the mode at ω, use the
    analytic β for the half-cell phase). The environment variable is read every time a source is
    built, so a probe can compare on against off within a single process.
    """
    return knobs.env("MODESRC_DISC") != "0"


def _freq_eff(f: float, dt: float) -> float:
    """Equivalent frequency of the centered time difference.

    In leapfrog, ``∂/∂t`` acting on ``e^{iωt}`` gives not ``iω`` but ``i·2 sin(ω dt/2)/dt``. So
    the mode the discrete scheme really supports is the one you get by solving the continuous
    mode problem at ``ω̃ = 2 sin(ω dt/2)/dt``.
    """
    om = 2.0 * np.pi * f
    return float(np.sin(om * dt / 2.0) / (np.pi * dt))


def _k_num(beta: float, dl: float) -> float:
    """Numerical wavenumber of the centered spatial difference.

    On the Yee grid ``∂/∂x`` gives ``-i·2 sin(k Δ/2)/Δ``, so the β the transverse operator sees
    is ``2 sin(k Δ/2)/Δ``. That is the quantity ModeSolver returns, and this inverts it to
    recover the ``k`` that actually advances on the grid. The half-cell phase has to use it:
    using the analytic beta is off by a phase of ``(k*Delta)^3/48 * (...)``, which turns directly into
    backward leakage in the TF/SF.
    """
    s = beta * dl / 2.0
    if not (-1.0 < s < 1.0):
        # fewer than ~3.2 cells per wavelength: the discrete scheme has no such propagating
        # mode at all, so fall back to the analytic value
        return float(beta)
    return float(2.0 * np.arcsin(s) / dl)


def _solve_plane(sim: td.Simulation, src, ax: int) -> td.Box:
    """The plane used for the mode solve and for injection.

    The ``src.size`` the user gives usually boxes in only the main part of the mode field
    (typically ``3×w_wg × 5×t_wg`` in the validation set). ModeSolver puts PEC on the boundary
    of that box, so what it solves for is the mode **as modified by the box**: both its β and
    its profile differ from the true mode of the open waveguide by O(the residual field at the
    box edge). Used as the TF/SF incident field, the backward component of that difference is a
    floor **independent of resolution** (measured at -81.5 dB for a 1.5×1.1 µm box and -99.0 dB
    for a 2.0×1.6 µm box).

    ``OPENEM_MODESRC_PAD`` gives the number of µm to grow outward (each of the two tangential
    axes grows by that much, clipped to the simulation domain). The default 0 keeps the old
    behavior bitwise.
    """
    pad = float(knobs.env("MODESRC_PAD"))
    center, size = list(src.center), list(src.size)
    if pad > 0:
        for a in range(3):
            # leave the normal axis alone, along with an axis that already fills the domain
            # (inf) and the axis flattened out in a 2D simulation
            if a == ax or not np.isfinite(size[a]) or size[a] <= 0:
                continue
            lo = max(center[a] - size[a] / 2 - pad, sim.center[a] - sim.size[a] / 2)
            hi = min(center[a] + size[a] / 2 + pad, sim.center[a] + sim.size[a] / 2)
            center[a], size[a] = 0.5 * (lo + hi), hi - lo
    return td.Box(center=center, size=size)


def _place(prof, coords, grid, axis_kind: tuple[str, str],
           tangs: tuple[int, int] = (1, 2)) -> np.ndarray:
    """Place the 2D profile from the mode box into the whole transverse plane, zero outside the
    box.

    ``tangs`` are the axis numbers of the two tangential axes (in cyclic order), and
    ``axis_kind`` says whether the component sits on the grid lines (``"edge"``) or at the cell
    centers (``"center"``) along each of them.
    **Matching is exact, by coordinate**, with no index arithmetic: the mode box covers only
    part of the transverse plane and its boundary need not line up with the grid (measured on
    ModalSourcesMonitors, whose source center x=-2.0 lands on the grid line -1.98889). A failed
    match raises: being off by one cell is a silent physical error.
    """
    n1, n2 = grid.axes[tangs[0]].n, grid.axes[tangs[1]].n
    out = np.zeros((n1, n2), dtype=np.complex128)
    prof = np.asarray(prof)
    idx = []
    keep = []
    for a, (kind, n) in zip(tangs, zip(axis_kind, (n1, n2))):
        want = np.asarray(coords[a], dtype=np.float64) * UM
        if n == 1:
            # 2D simulation (the axis has a single cell): ModeSolver returns the coordinates of
            # the folded plane (a scene with z size 0 returns z=0 for all of them) and does not
            # distinguish Yee center from edge. The field is uniform along that axis, so it goes
            # into the single cell there is; more than one point means this is not the folded
            # axis, and it still raises.
            if want.size != 1:
                raise ValueError(
                    f"{'xyz'[a]} axis has a single cell (2D simulation), but the mode profile "
                    f"supplied {want.size} coordinate points.")
            idx.append(np.zeros(1, dtype=np.int64))
            keep.append(np.ones(1, dtype=bool))
            continue
        tol = 1e-3 * float(np.min(grid.axes[a].dl))
        have = grid.axes[a].edges[:-1] if kind == "edge" else grid.axes[a].centers
        # ModeSolver's plane grid carries one extra layer beyond each domain boundary (measured
        # on MoS2Waveguide: 826 points along y against 824 grid cells; edge-type components also
        # carry the outermost boundary line). There is no Yee component to inject at those
        # positions, so they are trimmed; points inside the range still have to match exactly.
        inside = (want >= have[0] - tol) & (want <= have[-1] + tol)
        want_in = want[inside]
        if want_in.size == 0:
            raise ValueError(f"every {'xyz'[a]} coordinate of the mode profile lies outside "
                             "the grid domain.")
        j = np.searchsorted(have, want_in)
        j = np.clip(j, 1, have.size - 1)
        j = np.where(np.abs(have[j - 1] - want_in) <= np.abs(have[np.minimum(j, have.size - 1)] - want_in),
                     j - 1, np.minimum(j, have.size - 1))
        err = float(np.max(np.abs(have[j] - want_in)))
        if err > tol:
            raise ValueError(
                f"the {'xyz'[a]} coordinates of the mode profile do not match the grid's "
                f"{kind}: largest difference {err / UM:.3e} µm > tolerance "
                f"{tol / UM:.3e} µm. Being off by one cell is a silent physical error and "
                "cannot be let through.")
        if np.unique(j).size != j.size:
            raise ValueError(
                f"two points fall on the same grid index along {'xyz'[a]} of the mode profile: "
                f"{j.size} points occupy only {np.unique(j).size} cells. "
                "np.ix_ would overwrite them silently, and that is losing data, not "
                "approximating.")
        idx.append(j)
        keep.append(inside)
    out[np.ix_(idx[0], idx[1])] = prof[np.ix_(keep[0], keep[1])]
    return out



def _bg_index(sim: td.Simulation, freq0: float, src_center) -> float:
    """Background refractive index at the source plane (needed for the transverse wavenumber
    k_t=k0·n·sinθ at oblique incidence).

    It is the square root of the real part of ε(freq0) at the source center (:func:`_eps_at`);
    the beam and oblique-incidence sources in the validation set all sit in a uniform background
    (BoundaryConditions has a background of n=2, Metasurface/Grating sit in air).
    """
    eps_c = complex(_eps_at(sim, tuple(float(c) for c in src_center), freq0))
    return float(np.sqrt(eps_c.real))


def _eps_of(med, freq0: float) -> complex:
    try:
        return complex(med.eps_model(freq0))
    except Exception:
        return complex(getattr(med, "permittivity", 1.0))


def _eps_at(sim: td.Simulation, point: tuple, freq0: float) -> complex:
    """ε(freq0) at ``point``: a later structure overrides an earlier one, and if none of them
    contains the point, sim.medium is used. (SourceNormalization and Fresnel-style examples put
    the source inside an n=3.5 slab while sim.medium is air; computing k_t and the impedance
    from sim.medium then makes every oblique-incidence reflection wrong, measured at R 0.43
    against the analytic 0.60 for s polarization at θ=15°.)"""
    x, y, z = (np.array([point[0]]), np.array([point[1]]), np.array([point[2]]))
    for st in reversed(list(sim.structures)):
        try:
            if bool(np.asarray(st.geometry.inside(x, y, z)).ravel()[0]):
                return _eps_of(st.medium, freq0)
        except Exception:
            continue
    return _eps_of(sim.medium, freq0)


#: The four tangential components used for injection,
#: ``(key, td component name prefix, (kind on the t1 axis, kind on the t2 axis))``.
#: E_c sits at the cell center along its own axis and on the grid lines along the others; H is
#: the other way round. The keys are called Ey/Ez/Hy/Hz following the habit for an x normal,
#: but they really mean E_t1/E_t2/H_t1/H_t2 (cyclic order, which for an x normal coincides).
_TANGENTIAL_COMPS = (("Ey", "E", ("center", "edge")),
                     ("Ez", "E", ("edge", "center")),
                     ("Hy", "H", ("edge", "center")),
                     ("Hz", "H", ("center", "edge")))


class _Plane(NamedTuple):
    """Orientation of the injection plane: the normal ``ax``, the tangential cyclic order
    ``(t1, t2) = ((ax+1)%3, (ax+2)%3)``, the primal edge index ``ks`` the plane sits on, and the
    propagation ``direction`` (±1).

    A dozen functions along the beam carrier path and the ModeSource path all take this same
    group; spelling it out would add four positional parameters to every signature.
    """

    ax: int
    t1: int
    t2: int
    ks: int
    direction: int

    @classmethod
    def of(cls, ax: int, ks: int, direction: int) -> "_Plane":
        t1, t2 = cyclic_axes(ax)
        return cls(ax, t1, t2, ks, direction)


class _Beam(NamedTuple):
    """The analytic beam profile object, together with the background refractive index needed to
    evaluate it (``prof_src._beam_fields_on_component_grid(..., background_n=[n_bg])``)."""

    prof_src: object
    n_bg: float


class _SolveCtx(NamedTuple):
    """Everything needed to solve modes on one plane; the main path and the broadband correction
    loop read the same one.

    ``comps`` holds the four tangential components as ``(key, td component name, kinds)``
    (:func:`_tangential_comps`), ``dl_h`` is the spacing from the source plane to the upstream H
    plane in meters, and ``half_sign`` is the sign of the half-cell phase
    (:func:`_half_cell_sign`).
    """

    sim: object
    plane: object
    src: object
    grid: object
    comps: tuple
    dl_h: float
    half_sign: float


def _tangential_comps(t1: int, t2: int) -> tuple:
    """:data:`_TANGENTIAL_COMPS` with the concrete tangential axis names filled in:
    ``(key, td component name, kinds)``."""
    return tuple((nm, eh + "xyz"[t1 if nm[1] == "y" else t2], kind)
                 for nm, eh, kind in _TANGENTIAL_COMPS)


def _beam_profiles(beam: _Beam, grid, pl: _Plane) -> tuple[dict, tuple]:
    """Evaluate the four tangential components analytically on **our own Yee coordinates**
    (``_beam_fields_on_component_grid``, with the coordinates cast to complex128 to dodge a
    tidy3d dtype bug on real coordinates), with no interpolation. The E profiles are taken on
    the source plane and the H profiles on the neighboring H plane. Returns ``(prof, comps)``,
    where ``comps`` is the component table carrying the plane coordinates."""
    prof_src, n_bg = beam
    ax, t1, t2, ks, direction = pl
    e_plane = float(grid.axes[ax].edges[ks])
    kh = upstream_h_index(ks, direction)
    h_plane = float(grid.axes[ax].centers[kh])

    def _coords(kind: str, a: int):
        axo = grid.axes[a]
        return (axo.edges[:-1] if kind == "edge" else axo.centers) / UM

    def _prof(td_name: str, kind: tuple, plane_m: float):
        co = [None, None, None]
        co[ax] = [plane_m / UM]
        co[t1] = _coords(kind[0], t1)
        co[t2] = _coords(kind[1], t2)
        arrs = {c: co["xyz".index(c)] for c in "xyz"}
        out = prof_src._beam_fields_on_component_grid(
            x=np.asarray(arrs["x"], dtype=np.complex128),
            y=np.asarray(arrs["y"], dtype=np.complex128),
            z=np.asarray(arrs["z"], dtype=np.complex128),
            background_n=np.array([n_bg]))
        a3 = np.asarray(out[td_name]).reshape(
            len(arrs["x"]), len(arrs["y"]), len(arrs["z"]))
        return np.ascontiguousarray(
            np.transpose(a3, (t1, t2, ax))[:, :, 0].astype(np.complex128))

    comps = tuple((nm, td_nm, kind, e_plane if nm[0] == "E" else h_plane)
                  for nm, td_nm, kind in _tangential_comps(t1, t2))
    prof = {nm: _prof(td_nm, kind, pl) for nm, td_nm, kind, pl in comps}
    return prof, comps


def _aperture_weights(grid, bounds, a: int, kind: str) -> np.ndarray:
    """Aperture weight of each sample point along one tangential axis: the **overlap fraction**
    between the cell that sample point belongs to and the source box (``bounds``).

    A cell entirely inside the box gives 1, entirely outside gives 0, and a cell straddling the
    box edge gives its covered fraction (the source box boundary usually falls exactly on a grid
    line, so an edge sample point at the boundary gets ½; the aperture is then symmetric at both
    ends and the discrete injected power matches the analytic aperture power, which uses the
    trapezoidal convention of ``monitor_flux_colocated``).
    """
    _lo_s, _hi_s = bounds
    axo = grid.axes[a]
    if axo.is_flat:                       # 2D flat axis: that one cell always counts in full
        return np.ones(axo.n, dtype=np.float64)
    eg = np.asarray(axo.edges, dtype=np.float64) / UM
    if kind == "center":                  # sample at cell center, own cell = [edges[i], edges[i+1]]
        c_lo, c_hi = eg[:-1], eg[1:]
    else:                                 # sample on grid line edges[i], own cell = between centers
        cen = (eg[:-1] + eg[1:]) / 2.0
        c_lo = np.concatenate([eg[:1], cen[:-1]])
        c_hi = cen
    ov = np.minimum(c_hi, float(_hi_s[a])) - np.maximum(c_lo, float(_lo_s[a]))
    return np.clip(ov / (c_hi - c_lo), 0.0, 1.0)


def _apply_aperture(prof: dict, comps: tuple, grid, src, pl: _Plane) -> dict:
    """Source aperture truncation (2026-09-07).

    A td.PlaneWave / GaussianBeam with a finite aperture only injects current **inside its own
    box**, while the analytic beam profile fills the whole source plane; not truncating amounts
    to firing the plane wave across the entire cross section. EffectiveIndexApproximation's
    size=(0, inf, 0.2) only covers the 0.2 µm thick waveguide core (out of 2.76 µm along the
    whole z axis): without truncation the actual injected power is about 12 times the aperture
    power used for normalization, and the part that lands in the substrate/cladding runs
    straight downstream, flattening the T spectrum (measured T≈10.7, barely varying with
    frequency, against 0.25 to 0.80 from Tidy3D's solver).
    The aperture weights each sample point by the **overlap fraction** between its own cell and
    the source box (see :func:`_aperture_weights`); on a 2D flat axis that one cell always counts
    in full.
    All three conventions were tried (hard cut at 19 cells / hard cut at 20 sample points /
    overlap fraction): the L2 error in T is 1.7% / 3.4% / 1.2% respectively, and 24% in R for
    all three, so the overlap-fraction version was kept.
    """
    bounds = src.bounds
    t1, t2 = pl.t1, pl.t2

    def _box_w(a: int, kind: str) -> np.ndarray:
        return _aperture_weights(grid, bounds, a, kind)

    if not all(_box_w(a, k).min() >= 1.0 for a in (t1, t2) for k in ("edge", "center")):
        for nm, _td_nm, kind, _pl in comps:
            w = np.outer(_box_w(t1, kind[0]), _box_w(t2, kind[1]))
            prof[nm] = prof[nm] * w[:prof[nm].shape[0], :prof[nm].shape[1]]
        print(f"[openem] beam source aperture truncated to its own box (overlap-fraction "
              f"weighting): {'xyz'[t1]} {_box_w(t1, 'edge').sum():.1f}/{grid.axes[t1].n} cells, "
              f"{'xyz'[t2]} {_box_w(t2, 'edge').sum():.1f}/{grid.axes[t2].n} cells")
    return prof


def _colocated_power(prof: dict, grid, src, pl: _Plane) -> float | None:
    """Power over the plane under the colocated convention (2026-09-06): the same integration
    rule as FluxMonitor and tidy3d's |amp|². Returns None when the colocated geometry refuses,
    in which case the caller falls back to the Yee in-place integral.

    The two tangential axes of prof are in cyclic order (t1, t2) while the colocated geometry
    wants them ascending (u, v); for ax=1 those are the other way round, hence the transpose.
    Units: prof is evaluated on µm coordinates while the dS of the colocated integral is in m²,
    so it has to be converted back to µm² below (the units do not cancel out in the 1/√P
    normalization).
    """
    from openem import flux as _flux
    ax, t1, t2, ks = pl.ax, pl.t1, pl.t2, pl.ks
    u, v = sorted((t1, t2))
    if (t1, t2) != (u, v):
        eu, ev, hu, hv = (np.ascontiguousarray(a.T)
                          for a in (prof["Ez"], prof["Ey"], prof["Hz"], prof["Hy"]))
    else:
        eu, ev, hu, hv = prof["Ey"], prof["Ez"], prof["Hy"], prof["Hz"]
    nu, nv = grid.axes[u].n, grid.axes[v].n
    plane = np.stack([eu[:nu, :nv], ev[:nu, :nv], hu[:nu, :nv], hv[:nu, :nv]])[:, None]
    lo_b, hi_b = src.bounds

    class _Spec:
        # Duck-typed stand-in for a FluxMonitor to feed monitor_flux_colocated: it only reads
        # these fields (name/axis/plane_index/frac/normal_dir/t_bounds), while
        # model.FluxMonitor has no frac and demands freqs/source_spectrum, so it is not
        # constructed directly.
        name, axis, plane_index, frac, normal_dir = src.name or "beam", ax, int(ks), 0.0, +1
        t_bounds = ((float(lo_b[u]) * UM, float(hi_b[u]) * UM),
                    (float(lo_b[v]) * UM, float(hi_b[v]) * UM))
    try:
        # monitor_flux_colocated's area elements are in meters (grid.edges), while the
        # downstream scale=1/√p_inj is calibrated in µm² (the old algorithm's dA=axis_dl/UM is
        # in µm²), so this has to be divided by UM². Missing that on 2026-09-07 blew the fields
        # up by 1e6 (EdgeCoupler/EffectiveIndexApproximation reported divergence,
        # BilayerSiNSiGC was silently wrong).
        return abs(float(_flux.monitor_flux_colocated(plane.astype(np.complex128), grid, _Spec())[0])) / (UM ** 2)
    except (ValueError, NotImplementedError) as _e:
        # The source fills a periodic/Bloch axis (an oblique-incidence plane wave): the
        # colocated geometry refuses to drop the edge samples. A plane wave is smooth in the
        # transverse direction and the Yee in-place convention is accurate as it is, so fall
        # back to the original algorithm (the Fresnel criterion in test_oblique relies on it).
        print(f"[openem] the colocated convention is unusable for beam normalization "
              f"({str(_e)[:60]}...), falling back to the Yee in-place integral")
        return None


def _yee_power(prof: dict, grid, pl: _Plane) -> float:
    """Power over the plane from the Yee in-place integral, ½Re∫(E_t1·H_t2* - E_t2·H_t1*)dA (E
    and H are half a cell apart, which is negligible for a smooth profile).
    The area elements follow the same rule as the flux (flux.axis_dl: a true 2D flat axis uses
    FLAT_AXIS_M, everything else the actual dl); otherwise a single-cell transverse axis
    (GratingEfficiency at y=0.01 µm) would normalize the source by 0.01 while the flux counts
    1 µm, a factor of 100."""
    from openem import flux as _flux
    dl1 = np.asarray(_flux.axis_dl(grid.axes[pl.t1])[0], dtype=np.float64) / UM   # µm
    dl2 = np.asarray(_flux.axis_dl(grid.axes[pl.t2])[0], dtype=np.float64) / UM
    dA = np.outer(dl1[:prof["Ey"].shape[0]], dl2[:prof["Ey"].shape[1]])

    def _sz(a):
        return a[:dA.shape[0], :dA.shape[1]]
    p_inj = 0.5 * float(np.real(np.sum((_sz(prof["Ey"]) * np.conj(_sz(prof["Hz"]))
                                        - _sz(prof["Ez"]) * np.conj(_sz(prof["Hy"]))) * dA)))
    return abs(p_inj)


def _injected_power(prof: dict, grid, src, pl: _Plane) -> float:
    """Power over the plane: if the colocated convention produced nothing (it is switched off,
    or the geometry refused), fall back to the Yee in-place integral."""
    from openem import flux as _flux
    p_inj = (_colocated_power(prof, grid, src, pl)
             if _flux.flux_quadrature() == "colocated" else None)
    if p_inj is None:
        p_inj = _yee_power(prof, grid, pl)
    return p_inj


def _beam_carrier(sim, src, grid, beam: _Beam, pl: _Plane) -> ModeSource:
    """Any tidy3d beam profile (GaussianBeam / PlaneWave) -> a ModeSource carrier.

    Each component is evaluated analytically on **our own Yee coordinates**
    (:func:`_beam_profiles`), with no interpolation. The tangential pair is in cyclic order
    t1/t2, the E profiles are taken on the source plane and the H profiles on the neighboring H
    plane (when the phase front is curved or tilted, evaluating analytically is more accurate
    than the half-cell approximation, so ph=1). At oblique incidence the transverse phase
    e^{i k_t·r_t} is already built into the profile, so it agrees naturally with the wave vector
    of the Bloch boundary (the same k_t).
    """
    from openem.model import ModeSource as _MS
    ax, ks, direction = pl.ax, pl.ks, pl.direction
    prof, comps = _beam_profiles(beam, grid, pl)
    prof = _apply_aperture(prof, comps, grid, src, pl)         # source aperture truncation
    # ---- The same convention as ModeSource (2026-09-04, calibrated on an oblique-incidence
    # probe) ----
    # 1) Normalize the profile to **unit power**: a mode profile is unit power to begin with,
    #    and both the downstream injection coefficient and normalize.field_norm (the mode branch
    #    only divides by D) are calibrated on that convention. A beam profile comes at unit
    #    **amplitude** (|E|=1 for a plane wave), so using it directly loses a factor
    #    2η₀/(n·A): a 5° plane wave measured a transmission of 1.9e-3 against 1.0 from Tidy3D's
    #    solver. The power is computed numerically over the plane as
    #    ½Re∫(E_t1·H_t2* - E_t2·H_t1*)dA (see _injected_power). For a plane wave this lands
    #    back exactly on Tidy3D's |E| = √(2η₀/(n·A_µm²)) convention.
    # 2) Time table ×2: the ½ convention of mode injection (the amp_e table is twice the source
    #    time series, and normalize.source_table divides by 2 again). The carrier used not to
    #    multiply by 2, which left the field twice and the power four times too large.
    p_inj = _injected_power(prof, grid, src, pl)
    if p_inj > 0:
        scale = 1.0 / np.sqrt(p_inj)          # not prof / sqrt(p_inj): the rounding differs
        prof = {k: v * scale for k, v in prof.items()}
        print(f"[openem] beam carrier profile normalized to unit power: plane power "
              f"{p_inj:.3e} → ×{scale:.3f}")
    n = np.arange(sim.num_time_steps + 1, dtype=np.float64)
    dt = float(sim.dt)
    # The ×2 applies only at **normal incidence** (the real path, where the injection kernel
    # takes Re under the ½ convention). Oblique incidence goes down the complex Bloch path,
    # where the injection kernel consumes the full complex amplitude and must not be multiplied
    # by 2 again: measured on an oblique probe, a 5° plane wave transmitted 4.0 with the ×2
    # (against 1.0 from Tidy3D's solver) and exactly 1.0 without it; a GaussianBeam at 0° gave
    # 0.612 with the ×2 against Tidy3D's 0.613.
    # The test is **whether the solver takes the complex Bloch path** (the same as
    # model.any_bloch), not angle_theta>0: an oblique beam in an all-PML domain (the k-space
    # projection example in FieldProjections) takes the real path, and treating it as oblique
    # incidence would halve the field and quarter the power (measured: upstream flux 0.250 for a
    # Gaussian beam at θ=30°).
    from openem.scene import boundaries as _bnd
    _cplx = any(k != 0.0 for k in _bnd.bloch_vector(sim))
    _fac = 1.0 if _cplx else 2.0
    if _cplx:
        # The complex Bloch path injects the full complex amplitude (the analytic signal), whose
        # +f spectrum is twice that of Re(·) on the real path, while the normalization divisor
        # (source_table takes amp_e/2) still follows the real-path convention, leaving the field
        # twice and the power four times too large (measured on an oblique probe: 4.000 for a 5°
        # plane wave; 2.514/0.630 = 3.99, i.e. exactly 4, for a 10° Gaussian beam).
        # A further ×½ on the profile cancels it. (This and the ×scale above are two separate
        # multiplications and must not be merged into v*(scale*0.5): the rounding differs.)
        prof = {k: v * 0.5 for k, v in prof.items()}
    amp_e = _fac * np.asarray(src.source_time.amp_time(n * dt), dtype=np.complex128)
    amp_h = _fac * np.asarray(src.source_time.amp_time((n + 0.5) * dt), dtype=np.complex128)
    if direction < 0:
        # The profile is the backward beam evaluated analytically at src.direction="-" (its
        # tangential H is already flipped relative to the forward beam), while the injection
        # kernel follows the ModeSource convention and takes the **forward** profile, flipping H
        # itself through sources_setup._MODE_SIGNS[-1]. Flipping twice means injecting -(the
        # forward beam), which still travels forward (measured: GaussianBeam "-" gave an
        # upstream flux of +1.045 and a downstream one of -0.045). So flip H back to the forward
        # convention; the analytic phase on the E plane and on the upstream half-cell H plane is
        # left alone.
        prof["Hy"] = -prof["Hy"]
        prof["Hz"] = -prof["Hz"]
    return _MS(plane_index=ks, direction=direction, axis=ax,
               ey_inc=prof["Ey"], ez_inc=prof["Ez"],
               hy_inc=prof["Hy"], hz_inc=prof["Hz"],
               amp_e=np.ascontiguousarray(amp_e),
               amp_h=np.ascontiguousarray(amp_h),
               n_eff=complex(beam.n_bg))


def _beam_source(sim: td.Simulation, src, grid, profile_cls, **extra) -> ModeSource:
    """The part shared by all beam-like sources (GaussianBeam / PlaneWave) on the way to a
    ModeSource carrier: the normal, the direction, the plane index, the background refractive
    index at the source center, and then building the analytic profile with ``profile_cls``
    (``extra`` holds the fields specific to that profile class, such as the waist)."""
    ax = normal_axis(src.size)
    direction = sign_of(src.direction)
    ks = plane_index(grid, ax, src.center[ax])
    f0 = float(src.source_time.freq0)
    n_bg = _bg_index(sim, f0, src.center)
    prof_src = profile_cls(
        center=src.center, size=src.size, resolution=10.0, **extra,
        freqs=[f0], direction=str(src.direction),
        pol_angle=float(src.pol_angle),
        angle_theta=float(src.angle_theta),
        angle_phi=float(getattr(src, "angle_phi", 0.0)))
    return _beam_carrier(sim, src, grid, _Beam(prof_src, n_bg),
                         _Plane.of(ax, ks, direction))


def gaussian_beam(sim: td.Simulation, src, grid) -> ModeSource:
    """``td.GaussianBeam`` -> :class:`openem.model.ModeSource` (normal or oblique incidence).

    Normal incidence uses the real solver; oblique incidence (angle_theta>0) triggers the
    complex path through the scene's Bloch boundaries. For the profile and the carrier, see
    :func:`_beam_carrier`.
    """
    return _beam_source(sim, src, grid, td.GaussianBeamProfile,
                        waist_radius=float(src.waist_radius),
                        waist_distance=float(getattr(src, "waist_distance", 0.0)))


def plane_wave_beam(sim: td.Simulation, src, grid) -> ModeSource:
    """An oblique-incidence ``td.PlaneWave`` -> a ModeSource carrier (via PlaneWaveBeamProfile).

    A normal-incidence plane wave goes through :func:`scene.sources.plane_wave` (the 1D
    auxiliary grid, which is cheaper); this only takes oblique incidence (angle_theta>0), whose
    profile carries a transverse phase and therefore has to take the complex Bloch path.
    """
    return _beam_source(sim, src, grid, td.PlaneWaveBeamProfile)

#: Frequency points whose source spectrum falls below this fraction of the peak stay out of the
#: polynomial interpolation: Lₖ is a polynomial of degree K-1 and outside the band it amplifies
#: numerical noise until it diverges.
_BB_SPEC_FLOOR = 1e-8


def _lagrange_basis(nodes: np.ndarray, f: np.ndarray) -> np.ndarray:
    """The Lagrange basis on the nodes; returns shape ``(len(nodes), len(f))``."""
    nodes = np.asarray(nodes, dtype=np.float64)
    f = np.atleast_1d(np.asarray(f, dtype=np.float64))
    out = np.ones((nodes.size, f.size), dtype=np.float64)
    for k in range(nodes.size):
        for j in range(nodes.size):
            if j != k:
                out[k] *= (f - nodes[j]) / (nodes[k] - nodes[j])
    return out


def _bb_amps(source_time, nodes, num_steps: int, dt: float):
    """Time-domain coefficients of the broadband expansion, ``aₖ(t) = IFFT[Lₖ(f)·S(f)]``.

    Returns a whole-step set and a half-step set, both of shape ``(K, num_steps+1)``. The
    half-step set comes from sampling ``amp_time((n+0.5)dt)`` directly and then filtering: the
    two sequences have the same length and the same dt, so they share one frequency grid and the
    half-cell offset is already in the samples.
    """
    nodes = np.asarray(nodes, dtype=np.float64)
    n = np.arange(num_steps + 1, dtype=np.float64)
    out = []
    for off in (0.0, 0.5):
        s = np.asarray(source_time.amp_time((n + off) * dt), dtype=np.complex128)
        spec = np.fft.fft(s)
        fr = np.fft.fftfreq(s.size, dt)
        # Both gates are needed: (1) frequency points whose spectrum is too weak stay out (both
        # the division and the interpolation are meaningless there); (2) **the range has to be
        # limited to the node interval**: the profiles were only solved on [f_min, f_max], there
        # is nothing to go on outside it, and the degree K-1 Lₖ grows exponentially there
        # (measured injected power of 2.4e8 when this gate was missing, where it should be ~1).
        # A zero correction term outside the interval means falling back to P₀.
        # Both the test and the evaluation use **|f|**: amp_time returns the analytic signal
        # (exp(-i2πf₀t)·envelope) and numpy's FFT puts its content on the **negative** frequency
        # -f₀. Masking by the positive frequency interval filters the whole content away
        # (measured: aₖ dropped to 2e-8 and Σₖaₖ≈0), while not limiting the interval at all lets
        # Lₖ diverge at -f₀ (measured injected power 2.4e8). The physical frequency of a
        # negative-frequency component is just |f|, so both places take the absolute value.
        # Outside the interval the value is **clamped to the boundary** rather than zeroed:
        # zeroing means falling back to P₀, whereas the profile of the nearest node is clearly
        # closer to the truth than the one at f0; after clamping, Lₖ takes its boundary value,
        # which is bounded and continuous and does not run away the way polynomial extrapolation
        # would.
        af = np.clip(np.abs(fr), nodes.min(), nodes.max())
        keep = np.abs(spec) > _BB_SPEC_FLOOR * np.abs(spec).max()
        lag = np.zeros((len(nodes), s.size), dtype=np.float64)
        lag[:, keep] = _lagrange_basis(nodes, af[keep])
        out.append(np.stack([np.fft.ifft(lag[k] * spec)
                             for k in range(len(nodes))]))
    return out[0], out[1]


#: Only a mismatch above this value actually builds the broadband correction terms (which costs
#: K extra mode solves); below it, single-profile injection is kept, bitwise identical to what
#: it was before the change.
BB_BUILD_MIN = 1e-3


def _check_broadband(sim: td.Simulation, src, plane: td.Box) -> float:
    """``num_freqs>1``: measure how badly using only the freq0 profile misses at the edges of
    the band, and raise if it is too large.

    Tidy3D uses ``num_freqs`` frequency points for a Chebyshev interpolation (or a pole fit) to
    approximate **how the profile varies with frequency**; we take only the single profile at
    ``freq0``. The error of that approximation is the normalized overlap deficit
    ``1 - |<E(f)|E(f0)>|``.

    **Measured at the monitor frequencies, not at the band edges**: in a linear system the
    readout at frequency f is determined only by the source content at f, and the two injections
    agree exactly at f0, so the mismatch is an error only at the frequencies actually read out
    (MoS2Waveguide: the single monitor frequency is exactly f0, so the mismatch is 1.5e-01 at
    the band edges but 0 where it is read). With no frequency-carrying monitor, the original
    criterion (the band edges) is kept.

    Measured on the validation set (bandwidths all 9-10%, monitors covering the band):
    GroupDelayCalculation 1.3e-03, StripToSlotConverters 2.5e-03, MMIMeepBenchmark 1.4e-03.
    Since 2026-08-30 **real broadband injection exists** (see the expansion at the end of
    :func:`mode_source`), so this no longer blocks anything; it only measures the mismatch and
    hands it to the caller: above :data:`BB_BUILD_MIN` the correction terms are built, below it
    the single profile is kept (which saves K mode solves and stays bitwise identical to before
    the change). The historical limit of 1e-2, which used to fail closed when exceeded, is noted
    here for the record.

    Returns:
        The largest mismatch over the evaluation frequencies.
    """
    from tidy3d.plugins.mode import ModeSolver

    fg = np.asarray(src.frequency_grid, dtype=np.float64)
    f0 = float(src.source_time.freq0)
    mi = int(src.mode_index)

    from openem.scene.media import monitor_freqs
    fmons = monitor_freqs(sim)
    eval_freqs = (np.unique(np.concatenate(fmons)) if fmons
                  else np.array([fg.min(), fg.max()]))

    def prof(fr: float) -> np.ndarray:
        # the same plane as mode_source (including the OPENEM_MODESRC_PAD growth); with pad=0
        # it is src's own box
        d = _cached_solve(ModeSolver(
            simulation=sim, plane=plane,
            mode_spec=src.mode_spec, freqs=[float(fr)], colocate=True))
        v = np.concatenate([np.asarray(getattr(d, n).isel(f=0, mode_index=mi)).ravel()
                            for n in ("Ey", "Ez")])
        return v / np.linalg.norm(v)

    p0 = prof(f0)
    worst = max(1.0 - abs(np.vdot(prof(float(fr)), p0)) if fr != f0 else 0.0
                for fr in eval_freqs)
    if worst > BB_BUILD_MIN and not len(np.atleast_1d(src.frequency_grid)) > 1:
        # The expansion is genuinely needed but there are no nodes: without a frequency_grid
        # there is no interpolation basis, so it is better to stop here than to fall back to the
        # single profile silently.
        raise NotImplementedError(
            f"num_freqs={src.num_freqs} has a mismatch of {worst:.2e} and needs broadband "
            f"injection, but src.frequency_grid only holds "
            f"{len(np.atleast_1d(src.frequency_grid))} nodes.")
    return float(worst)


def _mode_profile_at(ctx: _SolveCtx, pl: _Plane,
                     f_inj: float) -> tuple[dict, complex, float, complex]:
    """Solve the mode on ``plane`` at ``f_inj`` and place it into the whole transverse plane:
    ``(prof, n_eff, k_prop, ph)``.

    ``prof`` holds the four tangential components (with neither the snap phase nor the half-cell
    phase applied); ``k_prop`` is the wavenumber that actually advances on the grid; ``ph`` is
    the half-cell phase factor of the H plane. Shared by the main path and the broadband
    correction loop.

    The incident plane of H sits half a cell **upstream in the propagation direction**: ks-½
    forward, ks+½ backward (the way ``dl_h`` is taken already reflects this). In both directions
    the H plane is half a cell earlier than the E plane, so the phase factor **carries no
    direction**. With direction included, backward directionality drops from 1.08e-07 to
    2.07e-02 (four sets of measurements are recorded in the comments of
    sources_setup._MODE_SIGNS).
    """
    from tidy3d.plugins.mode import ModeSolver

    src, grid = ctx.src, ctx.grid
    t1, t2, ax = pl.t1, pl.t2, pl.ax
    md = _cached_solve(ModeSolver(
        simulation=ctx.sim, plane=ctx.plane,
        mode_spec=src.mode_spec, freqs=[f_inj], colocate=False))
    sel = {"f": 0, "mode_index": int(src.mode_index)}
    n_eff = complex(np.asarray(md.n_eff.isel(**sel)))
    prof = {}
    for nm, td_nm, kind in ctx.comps:
        a = getattr(md, td_nm).isel(**sel)
        co = [np.asarray(a.coords[c], dtype=np.float64) for c in "xyz"]
        # Reshape explicitly by coordinate rather than using squeeze: in a 2D simulation some
        # other axis may also be down to a single point, and squeeze would collapse (n, 1) into
        # one dimension, which no longer matches _place's two-dimensional indexing.
        # First lay it out in 3D as (x, y, z), then transpose to (t1, t2, normal) and take the
        # layer at the source plane.
        a3 = np.asarray(a.data, dtype=np.complex128).reshape(
            co[0].size, co[1].size, co[2].size)
        arr = np.transpose(a3, (t1, t2, ax))[:, :, 0]
        prof[nm] = _place(arr, [co[0], co[1], co[2]], grid, kind, tangs=(t1, t2))
    beta = 2 * np.pi * f_inj / C_0 * n_eff.real         # 1/m (the β the transverse operator sees)
    dl_h = ctx.dl_h
    k_prop = _k_num(beta, dl_h) if _disc_on() else beta   # wavenumber advancing on the grid
    ph = np.exp(ctx.half_sign * 1j * k_prop * dl_h / 2.0)
    return prof, n_eff, k_prop, ph


def _broadband_terms(ctx: _SolveCtx, pl: _Plane, dx: float,
                     prof: dict, ph: complex) -> dict:
    """Broadband correction terms. The expansion (an identity, since Σₖ Lₖ ≡ 1):
        P(f) ≈ P₀ + Σₖ Lₖ(f)·(Pₖ - P₀)
    It is written as a base term plus corrections rather than as Σₖ LₖPₖ so that the base term
    is still the f0 profile times the full waveform, which leaves normalize.py's reading of
    amp_e as the source spectrum untouched. Each node adds the source-plane snap phase with its
    own k_prop_k (``dx = x_e - x_c``), and then the half-cell phase for the H components.
    """
    sim, src, comps = ctx.sim, ctx.src, ctx.comps
    dt = float(sim.dt)
    nodes = np.asarray(src.frequency_grid, dtype=np.float64)
    stacks: dict = {nm: [] for nm, _, _ in comps}
    for fk in nodes:
        fk_inj = _freq_eff(float(fk), dt) if _disc_on() else float(fk)
        pk_all, _n_eff_k, k_prop_k, ph_k = _mode_profile_at(ctx, pl, fk_inj)
        for nm, _td_nm, _kind in comps:
            pk = pk_all[nm] * np.exp(1j * pl.direction * k_prop_k * dx)  # snap phase, node k_prop
            stacks[nm].append(pk * (ph_k if nm[0] == "H" else 1.0))
    ae_k, ah_k = _bb_amps(src.source_time, nodes, sim.num_time_steps, dt)
    base = {"Ey": prof["Ey"], "Ez": prof["Ez"],
            "Hy": prof["Hy"] * ph, "Hz": prof["Hz"] * ph}
    return {
        "bb_ey": np.ascontiguousarray(np.stack(stacks["Ey"]) - base["Ey"]),
        "bb_ez": np.ascontiguousarray(np.stack(stacks["Ez"]) - base["Ez"]),
        "bb_hy": np.ascontiguousarray(np.stack(stacks["Hy"]) - base["Hy"]),
        "bb_hz": np.ascontiguousarray(np.stack(stacks["Hz"]) - base["Hz"]),
        "bb_amp_e": np.ascontiguousarray(ae_k),
        "bb_amp_h": np.ascontiguousarray(ah_k),
    }


def mode_source(sim: td.Simulation, src, grid) -> ModeSource:
    """``td.ModeSource`` -> :class:`openem.model.ModeSource`.

    Any normal. The tangential pair is taken in cyclic order
    ``(t1, t2) = ((ax+1)%3, (ax+2)%3)``: the Yee curl keeps its form under cyclic permutation,
    so the sign table of the TF/SF correction is shared by all three normals (the x normal
    degenerates bitwise into the old implementation; directionality for the y and z normals was
    measured separately).
    """
    from openem.model import ModeSource as _MS

    ax = normal_axis(src.size)
    t1, t2 = cyclic_axes(ax)
    if src.mode_spec.bend_radius is not None:
        raise NotImplementedError(
            f"the mode_spec of ModeSource carries bend_radius={src.mode_spec.bend_radius}: "
            "the phase front of a bent mode is not a plane in Cartesian coordinates, so a "
            "TF/SF injection on a plane with a single β gives an incident field that does not "
            "satisfy the source-free Maxwell equations in that region, and the TF/SF becomes a "
            "positive feedback loop. ModesBentAngled was measured to diverge exponentially "
            "between steps 8,000 and 12,000 (the 4,800-step early termination happened to hide "
            "it)."
        )
    plane = _solve_plane(sim, src, ax)
    bb_mismatch = 0.0
    if int(src.num_freqs) != 1:
        bb_mismatch = _check_broadband(sim, src, plane)
    direction = sign_of(src.direction)
    ks = plane_index(grid, ax, src.center[ax])
    f0 = float(src.source_time.freq0)

    # **Fully discrete, consistent injection**: the mode of continuous Maxwell is not an
    # eigenmode of Yee + leapfrog, so using it as the TF/SF incident field leaves a backward
    # mode in the difference, and that is the -80 dB floor. The discrete eigenmode is obtained
    # by solving the same transverse operator (ModeSolver already uses Yee forward/backward
    # differences) at ω̃, then inverting the β it returns into k_num.
    f_inj = _freq_eff(f0, float(sim.dt)) if _disc_on() else f0
    comps = _tangential_comps(t1, t2)     # see _TANGENTIAL_COMPS for the table and center/edge
    pl = _Plane.of(ax, ks, direction)
    ctx = _SolveCtx(sim, plane, src, grid, comps,
                    dl_h=float(grid.axes[ax].dl[upstream_h_index(ks, direction)]),
                    half_sign=_half_cell_sign())
    prof, n_eff, k_prop, ph = _mode_profile_at(ctx, pl, f_inj)
    # Snap-phase compensation on the source plane (2026-09-05): Tidy3D references the mode phase
    # to **the source center the user gave**, x_c, and propagates it to the snapped grid line
    # x_e as e^{iβ·dir·(x_e-x_c)}. We used to reference x_e instead, which gave every downstream
    # field and mode amplitude an extra global phase of β·(x_e-x_c) (a gradient probe on a 2D
    # straight waveguide measured +3.4° = β·0.00442 µm against the point-by-point phase
    # difference from Tidy3D's solver; there is no difference when the source lands exactly on a
    # grid line). Amplitude and power are unaffected, but the adjoint gradient is the small real
    # part of a large complex number, so a few degrees make it 60% wrong (A_box_eps off by
    # 1.61×). Each node of the broadband correction applies this with its own β_k.
    x_e = float(grid.axes[ax].edges[ks])
    x_c = float(src.center[ax]) * UM
    snap_phase = np.exp(1j * direction * k_prop * (x_e - x_c))
    prof = {k: v * snap_phase for k, v in prof.items()}

    n = np.arange(sim.num_time_steps + 1, dtype=np.float64)
    dt = float(sim.dt)
    amp_e = np.asarray(src.source_time.amp_time(n * dt), dtype=np.complex128)
    amp_h = np.asarray(src.source_time.amp_time((n + 0.5) * dt), dtype=np.complex128)

    # ---- Broadband correction terms ----
    # They are only built (_broadband_terms) when the single-profile mismatch exceeds
    # BB_BUILD_MIN; when the mismatch is small enough that no expansion is needed, the
    # correction terms are empty and this degenerates bitwise into single-profile injection.
    bb: dict = {}
    if bb_mismatch > BB_BUILD_MIN:
        bb = _broadband_terms(ctx, pl, x_e - x_c, prof, ph)

    return _MS(
        plane_index=ks, direction=direction, axis=ax,
        ey_inc=np.ascontiguousarray(prof["Ey"]),
        ez_inc=np.ascontiguousarray(prof["Ez"]),
        hy_inc=np.ascontiguousarray(prof["Hy"] * ph),
        hz_inc=np.ascontiguousarray(prof["Hz"] * ph),
        amp_e=np.ascontiguousarray(amp_e), amp_h=np.ascontiguousarray(amp_h),
        n_eff=n_eff, **bb)
