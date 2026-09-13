# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Three kinds of monitor: flux, frequency-domain field, time-domain field."""

from __future__ import annotations

import dataclasses

import numpy as np
import tidy3d as td

from openem import flux as flux_mod
from openem import waveform
from openem.grid import Grid, normal_axis, transverse_axes
from openem.scene import modes, projection
from openem.scene._util import plane_index, sign_of
from openem.model import (
    ModeMonitorSpec,
    COMPONENTS,
    FieldMonitor,
    FieldTimeMonitor,
    FluxMonitor,
    FluxTimeMonitor,
    BoxFluxMonitor,
    PermittivityMonitor,
    ProjectionMonitor,
)
from openem.units import UM

#: Ceiling on the device buffer of one FieldMonitor. An H800 has 80 GB and the fields plus the
#: coefficients already take most of it; 24 GB covers every case in the validation set (the largest
#: one measured is far below it), so going over it is almost always a configuration mistake.
FIELD_MONITOR_MAX_BYTES = 24e9

#: Tolerance (in meters) for deciding whether a grid point falls inside the transverse range of a
#: FluxMonitor. The floating-point representations a parent and a child monitor give for the same
#: coordinate differ by about 5e-16, so 1e-12 is plenty.
COORD_TOL = 1e-12

#: Ceiling on the time-domain buffer. Going over it fails closed: truncating silently would leave a
#: stretch missing from the comparison without a word.
TIME_BUFFER_LIMIT_BYTES = 4 << 30

#: Name template for the six face monitors a 3-D FluxMonitor is split into (same shape as
#: projection.SURFACE_NAME).
BOX_FACE_NAME = "{name}__face_{axis}{sign}"


def _monitor_box(mon, grid: Grid) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    """The cell index box ``(origin, box)`` the monitor has to store, at the Yee positions.

    The two kinds of axis are handled separately:

    - **An axis with thickness**: cover the primal boundaries inside the monitor's range, with
      **2 extra cells at each end**. The margin is necessary: colocating onto boundary ``i`` needs
      cell ``i-1``, and the Yee points Tidy3D reports also reach past the monitor's nominal bounds
      (measured on CavityFOM: 2 cells past).
    - **A zero-thickness axis** (the normal of a point or plane monitor): take the **3-cell
      neighborhood** of the cell containing that point. A component may sit at the cell center or
      on a cell boundary, and 3 cells cover the linear interpolation in either case.
    """
    lo, hi = mon.bounds
    origin, box = [], []
    for ax in range(3):
        e = grid.axes[ax].edges
        n = grid.axes[ax].n
        a, b = float(lo[ax]) * UM, float(hi[ax]) * UM
        if float(mon.size[ax]) == 0.0:
            c = int(np.clip(np.searchsorted(e, a, side="right") - 1, 0, n - 1))
            i0 = max(c - 1, 0)
            i1 = min(c + 1, n - 1)
        else:
            first = int(np.searchsorted(e, a, side="left"))
            last = int(np.searchsorted(e, b, side="right")) - 1
            first, last = max(first, 0), min(max(last, first), n)
            i0 = max(first - 2, 0)
            # One extra cell at the high end: after the colocation endpoint backs off by 1e-9 it
            # still has to be covered by the index box
            i1 = min(last + 3, n - 1)
        origin.append(i0)
        box.append(max(i1 - i0 + 1, 1))
    return tuple(origin), tuple(box)                        # type: ignore[return-value]


def _apodization(mon) -> tuple[float | None, float | None, float] | None:
    """``mon.apodization`` → ``(start, end, width)`` in seconds, or ``None`` for no apodization.

    With neither start nor end given, the window is identically 1, which **is mathematically no
    apodization**, so it is normalized to ``None`` (giving width alone is the same: width only sets
    the width of the ramp, and with no ramp it does nothing).
    """
    apo = getattr(mon, "apodization", None)
    if apo is None:
        return None
    start, end, width = (getattr(apo, f, None) for f in ("start", "end", "width"))
    if start is None and end is None:
        return None
    if width is None:
        raise ValueError(
            f"the apodization of monitor {mon.name} gives start/end but no width")
    return (
        None if start is None else float(start),
        None if end is None else float(end),
        float(width),
    )


def _refuse_monitor_extras(mon, allow_interval_space: bool = False) -> None:
    """Options on a monitor that are not supported yet.

    ``allow_interval_space``: a frequency-domain FieldMonitor may carry spatial downsampling. The
    solver still accumulates at **full resolution** (what Tidy3D downsamples is a subset of the
    same Yee points, so taking the subset afterwards is exact), and the stride is recorded in
    :class:`FieldMonitor.interval_space` for the write-back. Time-domain monitors are let through
    the same way (2026-09-04, AnimationTutorial): the buffer stays at full resolution and the
    subset is taken at write-back following Tidy3D's point-selection rule
    (``nb.backend._downsample_inds``).
    """
    step = tuple(int(v) for v in getattr(mon, "interval_space", (1, 1, 1)))
    if step != (1, 1, 1) and not allow_interval_space:
        raise NotImplementedError(
            f"monitor {mon.name} has interval_space={step} (spatial downsampling), not supported "
            "yet"
        )


def _refuse_apodization(mon) -> None:
    """Monitor types that get expanded into other monitors (mode monitors) would drop the
    apodization, so they can only be refused.

    Projection monitors are not in this group: ``projection.surface_monitors`` passes the
    apodization straight into the face monitors.
    """
    if _apodization(mon) is not None:
        raise NotImplementedError(
            f"monitor {mon.name} ({type(mon).__name__}) carries an apodization; this kind of "
            "monitor does not support apodization yet")


def _time_window(sim: td.Simulation, mon) -> tuple[int, int, int, int]:
    """``(step_begin, step_end, interval, num_slots)`` of a time-domain monitor.

    The sampling window comes straight from Tidy3D's ``time_inds`` rather than being rewritten
    here.
    """
    beg, end = (int(v) for v in mon.time_inds(np.asarray(sim.tmesh)))
    interval = int(mon.interval)
    num_slots = (max(end - beg, 0) + interval - 1) // interval
    return beg, end, interval, num_slots


def _t_bounds(mon, ax: int) -> tuple[tuple[float, float], tuple[float, float]]:
    """The physical bounds of the two tangential axes in meters, ``((lo1, hi1), (lo2, hi2))``,
    needed by the colocated integration convention (flux.monitor_flux_colocated /
    flux_time_geometry)."""
    lo_b, hi_b = mon.bounds
    return tuple((float(lo_b[t]) * UM, float(hi_b[t]) * UM)
                 for t in transverse_axes(ax))


def field_time_monitor(sim: td.Simulation, mon, grid: Grid) -> FieldTimeMonitor:
    """``td.FieldTimeMonitor`` → :class:`FieldTimeMonitor`."""
    _refuse_monitor_extras(mon, allow_interval_space=True)
    beg, end, interval, num_slots = _time_window(sim, mon)
    origin, box = _monitor_box(mon, grid)
    comps = tuple(COMPONENTS.index(str(f)) for f in mon.fields)

    nbytes = len(comps) * num_slots * box[0] * box[1] * box[2] * 4
    if nbytes > TIME_BUFFER_LIMIT_BYTES:
        raise NotImplementedError(
            f"the time-domain buffer of monitor {mon.name} needs {nbytes / (1 << 30):.1f} GiB "
            f"({len(comps)} components × {num_slots} slots × {box}), over the ceiling of "
            f"{TIME_BUFFER_LIMIT_BYTES >> 30} GiB"
        )
    return FieldTimeMonitor(
        name=mon.name,
        comps=comps,
        origin=origin,
        box=box,
        step_begin=beg,
        step_end=end,
        interval=interval,
        num_slots=num_slots,
    )


def field_monitor(sim: td.Simulation, mon, grid: Grid) -> FieldMonitor:
    """``td.FieldMonitor`` → :class:`FieldMonitor`: work out the cell index box to store.

    The box keeps one cell of margin at the low end, because colocating onto a cell boundary needs
    ``F[i-1]`` and ``F[i]``.
    """
    _refuse_monitor_extras(mon, allow_interval_space=True)
    origin, box = _monitor_box(mon, grid)
    # Device memory guard: the device buffer is (6, nf, ni, nj, nk) of re+im float64. A monitor
    # with spatial downsampling is stored at full resolution (the subset is taken exactly at
    # write-back), so a large volume plus many frequency points blows up device memory here.
    # Better to show the user the arithmetic while building the scene than to let cupy throw an
    # opaque OOM.
    nf = int(np.atleast_1d(np.asarray(mon.freqs)).size)
    nbytes = 6 * nf * int(np.prod(box)) * 2 * 8
    if nbytes > FIELD_MONITOR_MAX_BYTES:
        raise NotImplementedError(
            f"the device buffer of FieldMonitor {mon.name} needs {nbytes / 1e9:.1f} GB "
            f"(box {box} × {nf} frequency points × 6 components × complex float64), over the "
            f"ceiling of {FIELD_MONITOR_MAX_BYTES / 1e9:.0f} GB. Use fewer frequency points, "
            "shrink the box, or push the interval_space downsampling all the way down into "
            "storage (not implemented yet).")
    return FieldMonitor(
        name=mon.name,
        freqs=np.asarray(mon.freqs, dtype=np.float64),
        origin=origin,
        box=box,
        apodization=_apodization(mon),
        interval_space=tuple(int(v) for v in getattr(mon, "interval_space", (1, 1, 1))),
    )



def _flux_plane(mon, grid: Grid) -> tuple[int, tuple | None]:
    """The spatial convention of a frequency-domain flux plane, ``(axis, transverse)`` (a range of
    staggered grid points).

    The transverse range is the primal grid points inside the monitor's bounds. It is left as None
    when the monitor covers the whole face (this saves a multiplication and makes "not filled in"
    and "filled in as the whole face" mean the same thing). **It has to be given when the monitor
    does not cover the full transverse extent**, otherwise the whole plane gets integrated.
    """
    ax = normal_axis(mon.size)
    lo, hi = mon.bounds
    trans, full = [], True
    for t in transverse_axes(ax):
        e = grid.axes[t].edges[:-1]
        a, b = float(lo[t]) * UM, float(hi[t]) * UM
        keep = np.flatnonzero((e >= a - COORD_TOL) & (e <= b + COORD_TOL))
        if keep.size == 0:
            raise ValueError(
                f"flux monitor {mon.name} covers no grid point on {'xyz'[t]}: range "
                f"[{a / UM:.4f}, {b / UM:.4f}] µm")
        trans.append((int(keep[0]), int(keep[-1])))
        full &= (keep.size == grid.axes[t].n)
    return ax, (None if full else (trans[0], trans[1]))


def flux_monitor(sim: td.Simulation, mon, grid: Grid) -> FluxMonitor:
    """``td.FluxMonitor`` → :class:`FluxMonitor`."""
    ax, transverse = _flux_plane(mon, grid)
    freqs = np.asarray(mon.freqs, dtype=np.float64)
    return FluxMonitor(
        name=mon.name,
        axis=ax,
        transverse=transverse,
        # Tangential physical bounds in meters: needed by the colocated integration convention
        # (flux.monitor_flux_colocated, 2026-09-06)
        t_bounds=_t_bounds(mon, ax),
        plane_index=plane_index(grid, ax, mon.center[ax]),
        normal_dir=sign_of(mon.normal_dir),
        freqs=freqs,
        source_spectrum=waveform.spectrum(
            sim.sources[0].source_time, sim.dt, sim.num_time_steps, freqs
        ),
        apodization=_apodization(mon),
    )


def flux_time_monitor(sim: td.Simulation, mon, grid: Grid) -> FluxTimeMonitor:
    """``td.FluxTimeMonitor`` → :class:`FluxTimeMonitor`.

    The spatial convention is **colocated integration** (the reference implementation hard-codes
    ``use_colocated_integration=True`` for FluxTimeMonitor): the tangential directions store the
    monitor's bound coordinates ``t_bounds`` (the face elements are truncated exactly by them, see
    ``flux.flux_time_geometry``), **not** the staggered grid-point range of :func:`_flux_plane`.
    The time window goes through Tidy3D's own ``time_inds`` (the same path as
    :func:`field_time_monitor`). The buffer is one float64 scalar per slot, so it does not need the
    :data:`TIME_BUFFER_LIMIT_BYTES` gate.

    The plane may sit **off the Yee boundaries** (none of NanobeamCavity's 6 faces is on one): it
    is interpolated to the exact coordinate under the ``colocate=True`` convention (weight
    ``frac``). When it lands all but exactly on a boundary (within COORD_TOL) it collapses to
    ``frac=0`` and takes exactly the same arithmetic as sitting on the boundary.
    """
    _refuse_monitor_extras(mon)
    if all(float(v) > 0 for v in mon.size):
        raise NotImplementedError(
            f"monitor {mon.name} is a 3-D FluxTimeMonitor (a closed box), not supported yet; in "
            "the validation set NanobeamCavity splits the box into six planes itself")
    if getattr(mon, "exclude_surfaces", None):
        raise NotImplementedError(
            f"monitor {mon.name} carries exclude_surfaces, not supported yet")
    ax = normal_axis(mon.size)
    e = grid.axes[ax].edges
    x0 = float(mon.center[ax]) * UM
    km = int(np.searchsorted(e, x0, side="right")) - 1
    if not 0 <= km < grid.axes[ax].n:
        raise ValueError(
            f"monitor {mon.name} has its plane {'xyz'[ax]}={x0 / UM:.6f} µm outside the grid")
    frac = (x0 - float(e[km])) / float(e[km + 1] - e[km])
    if x0 - float(e[km]) < COORD_TOL:
        frac = 0.0
    elif float(e[km + 1]) - x0 < COORD_TOL:
        km, frac = km + 1, 0.0
    beg, end, interval, num_slots = _time_window(sim, mon)
    out = FluxTimeMonitor(
        name=mon.name,
        axis=ax,
        plane_index=km,
        normal_dir=sign_of(mon.normal_dir),
        t_bounds=_t_bounds(mon, ax),
        frac=frac,
        step_begin=beg,
        step_end=end,
        interval=interval,
        num_slots=num_slots,
    )
    flux_mod.flux_time_geometry(grid, out)   # fail closed: verify interpolation and face range now
    return out


def permittivity_monitor(sim: td.Simulation, mon, grid: Grid) -> PermittivityMonitor:
    """``td.PermittivityMonitor`` → :class:`PermittivityMonitor`.

    **apodization is not checked**: it windows a time integral, and ε does not depend on time.
    ``interval_space`` still has to be checked, because spatial downsampling would make the number
    of output points disagree.
    """
    _refuse_monitor_extras(mon)
    origin, box = _monitor_box(mon, grid)
    return PermittivityMonitor(
        name=mon.name,
        freqs=np.asarray(mon.freqs, dtype=np.float64),
        origin=origin,
        box=box,
    )


@dataclasses.dataclass
class MonitorSet:
    """What :func:`build` produces: the monitors sorted into lists by type (the eight monitor
    fields of Scene)."""
    flux: list[FluxMonitor] = dataclasses.field(default_factory=list)
    field: list[FieldMonitor] = dataclasses.field(default_factory=list)
    field_time: list[FieldTimeMonitor] = dataclasses.field(default_factory=list)
    flux_time: list[FluxTimeMonitor] = dataclasses.field(default_factory=list)
    permittivity: list[PermittivityMonitor] = dataclasses.field(default_factory=list)
    projection: list[ProjectionMonitor] = dataclasses.field(default_factory=list)
    box_flux: list[BoxFluxMonitor] = dataclasses.field(default_factory=list)
    mode: list[ModeMonitorSpec] = dataclasses.field(default_factory=list)


def _diffraction_as_flux(mon) -> td.FluxMonitor:
    """A diffraction monitor is a plane, and underneath it takes the colocated tangential phasors
    (the same as a FluxMonitor); the near-field to diffraction-order FFT happens on the readout
    side (nb.backend._diffraction_data)."""
    return td.FluxMonitor(center=mon.center, size=mon.size, freqs=mon.freqs,
                          name=mon.name, normal_dir=mon.normal_dir)


def _box_flux(sim: td.Simulation, mon, grid: Grid) -> tuple[list[FluxMonitor], BoxFluxMonitor]:
    """A 3-D FluxMonitor is the total flux through a closed surface, split into six faces (for the
    meaning and the signs see :class:`BoxFluxMonitor`).

    exclude_surfaces (BullseyeCavityPSO, LedLEECalculation): an excluded face gets no monitor and
    does not count towards the total flux; a face monitor itself must not carry this field (Tidy3D
    only allows it on a box monitor, and updated_copy raises on it), so it is cleared explicitly.
    """
    faces: list[FluxMonitor] = []
    names = []
    excl = set(getattr(mon, "exclude_surfaces", None) or ())
    for fm, sign, ax in projection.surface_monitors(mon, sim, name_template=BOX_FACE_NAME):
        axl = "xyz"[ax]              # surface_monitors only ever returns int axes
        if f"{axl}{sign}" in excl:
            continue
        nm = fm.name
        faces.append(flux_monitor(
            sim, mon.updated_copy(center=fm.center, size=fm.size,
                                  name=nm, normal_dir=sign,
                                  exclude_surfaces=None), grid))
        names.append(nm)
    return faces, BoxFluxMonitor(name=mon.name, face_names=tuple(names))


def build(sim: td.Simulation, grid: Grid) -> MonitorSet:
    """Every monitor, dispatched by type.

    ``PermittivityMonitor`` has to be tested before ``FieldMonitor``: the two have the same data
    shape, but the former is not a subclass of the latter, so getting the order wrong would only
    miss a match rather than produce a wrong one; this is just being explicit about it.

    **A far-field projection monitor is expanded into six face FieldMonitors**: the solver only has
    to produce the phasors on those six faces, and the transform itself is left to the client in
    post-processing (``scene/projection.py``).
    """
    out = MonitorSet()
    for mon in sim.monitors:
        if modes.is_mode_monitor(mon):
            # expanded into one planar FieldMonitor; the mode decomposition is done in
            # ``scene/modes.py``
            _refuse_apodization(mon)
            fm = modes.plane_monitor(mon)
            out.field.append(field_monitor(sim, fm, grid))
            out.mode.append(ModeMonitorSpec(name=mon.name, plane_name=fm.name))
        elif projection.is_projection(mon):
            # no need to refuse apodization: surface_monitors has already passed it to every face
            # monitor
            faces = projection.surface_monitors(mon, sim)
            for fm, _, _ in faces:
                out.field.append(field_monitor(sim, fm, grid))
            out.projection.append(ProjectionMonitor(
                name=mon.name,
                surface_names=tuple(fm.name for fm, _, _ in faces),
                normal_dirs=tuple(sign for _, sign, _ in faces)))
        elif isinstance(mon, td.PermittivityMonitor):
            out.permittivity.append(permittivity_monitor(sim, mon, grid))
        elif isinstance(mon, td.FieldTimeMonitor):
            out.field_time.append(field_time_monitor(sim, mon, grid))
        elif isinstance(mon, td.FluxTimeMonitor):
            out.flux_time.append(flux_time_monitor(sim, mon, grid))
        elif type(mon).__name__ == "DiffractionMonitor":
            out.flux.append(flux_monitor(sim, _diffraction_as_flux(mon), grid))
        elif isinstance(mon, td.FieldMonitor):
            out.field.append(field_monitor(sim, mon, grid))
        elif type(mon).__name__ == "ModeSolverMonitor":
            # Mode profile monitor: it needs no FDTD field and is computed at assembly time with
            # tidy3d's local mode solver (the ModeSolverMonitor branch of nb.backend.build).
            continue
        elif isinstance(mon, td.FluxMonitor):
            if all(v > 0 for v in mon.size):
                faces, box = _box_flux(sim, mon, grid)
                out.flux.extend(faces)
                out.box_flux.append(box)
            else:
                out.flux.append(flux_monitor(sim, mon, grid))
        else:
            raise NotImplementedError(
                f"only FluxMonitor / FluxTimeMonitor / FieldMonitor / "
                f"FieldTimeMonitor / PermittivityMonitor / ModeMonitor are supported, "
                f"got {type(mon).__name__} ({mon.name})"
            )
    return out
