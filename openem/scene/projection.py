# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Far-field projection monitors: split into six near-field surface monitors, with the projection
itself left to the client.

What Tidy3D's solver does with a ``FieldProjection*Monitor`` is record the tangential field on the
**surface** of the monitor box and then perform a Stratton-Chu surface integral. The client
already implements the second step (``FieldProjector.project_fields``), so only the first is ours.
That is the same idea as taking the subpixel average from ``sim.epsilon``: **whatever can be taken
from the client is not rewritten here**.

The reference outputs store **only the projected result** (``FieldProjectionAngleData``) and not
the near-field surfaces, so the check has to run end to end: our solver produces the phasors on six
surfaces, the client projects them, and the result is compared against the reference.
"""

from __future__ import annotations

import contextlib

import numpy as np
import tidy3d as td

from openem import colocate, knobs
from openem.scene import modes
from openem.scene.normalize_td import normalization  # noqa: F401  re-export (global normalization convention)
from openem.units import UM


#: Order and naming suffix of the six surfaces, one per ``normal_dir``.
FACES = (("x", "-"), ("x", "+"), ("y", "-"), ("y", "+"), ("z", "-"), ("z", "+"))

#: Name template for a surface monitor. It carries an unlikely separator because these names share a
#: namespace with the user's own monitor names, and a collision would silently overwrite.
SURFACE_NAME = "{name}__n2f_{axis}{sign}"

#: Tolerance for matching coordinates, in um. Both sides are on the same grid, so the only
#: difference should be the floating-point representation.
COORD_ATOL = 1e-9

#: Tolerance between an index-box endpoint and a colocation endpoint, in um. Far below any grid
#: step; the finest is about 1e-3 um.
COORD_EDGE_TOL = 1e-6

#: An extra factor used **only** for the far field. Tidy3D's far-field convention measures the ``r``
#: in ``1/r`` in **um**, while our fields are in SI metres. Since ``E`` is proportional to ``1/r``,
#: the two differ by exactly 1e6.
#:
#: Vacuum against an analytic solution **cannot distinguish** this, since
#: ``eta0*k*IL/(4*pi*1 m)*1e-6`` and ``eta0*k*IL/(4*pi*1e6)`` are the same number, so it had to be
#: settled against the reference outputs: integrating the reference far field over a sphere with
#: ``r`` taken as the bare ``proj_distance`` value of 1e6 gives 1.067e+05, the same order as the
#: reference's own flux of 1.301e+05, whereas taking ``r`` in metres is off by 1e12. With the
#: factor applied, one scene came out at amplitude ratios of 0.9983 and 0.9971 with a radiation
#: pattern spread of 1.7e-02.
#:
#: **Far field only.** Folding it into :func:`normalization` makes flux come out 1e-12 too small.
FAR_FIELD_SCALE = UM


def is_projection(mon) -> bool:
    """Whether this is a far-field projection monitor, in the angle domain or in k space."""
    return isinstance(mon, td.components.monitor.AbstractFieldProjectionMonitor)


def _face_monitor(center, size, freqs, apod, name: str) -> td.FieldMonitor:
    """One near-field surface monitor.

    ``name`` is assembled by the caller from ``name_template``; the template takes three fields, and
    passing them flat would leak the string assembly down into here. ``colocate=True`` is left at
    its default because the data has already been colocated onto the plane in :func:`field_data`,
    and the two conventions have to agree.
    """
    return td.FieldMonitor(center=tuple(center), size=tuple(size),
                           freqs=list(np.asarray(freqs, dtype=float)),
                           apodization=apod, name=name)


def surface_monitors(mon, sim: td.Simulation,
                     name_template: str = SURFACE_NAME) -> list[tuple[td.FieldMonitor, str, int]]:
    """Projection box -> ``[(surface monitor, normal_dir, normal axis)]``.

    ``name_template`` is the naming template for a surface monitor (default :data:`SURFACE_NAME`;
    ``monitors.build`` passes its own ``__face_`` template when splitting a 3D FluxMonitor).

    Two surfaces are generated only for axes that have **thickness**. A zero-thickness axis means
    the monitor is itself a surface, in which case it is the only near-field surface. Tidy3D allows
    that, for instance when projecting in one direction only.

    Outward normals: the ``-`` surface faces ``-`` and the ``+`` surface faces ``+``, which is the
    ``normal_dirs`` convention that ``FieldProjector.from_near_field_monitors`` expects.

    The parent monitor's ``apodization`` must be **passed through to every surface**: the surface
    monitors are where the solver actually integrates the DFT, and without passing it down the
    apodization is silently lost.

    The transverse extent of a plane monitor has to be cut to something **finite** (one example
    writes ``size=inf``), and the convention is **the whole grid including the PML**. Tidy3D
    ``discretize``s an infinite monitor across the full grid including the PML (measured at +-6.84
    against a domain edge of +-6.6), where the decaying tail inside the PML acts as a soft window.
    Cutting at the domain edge is a hard truncation, which measured about 7% off the main lobe of
    one pattern and pushed that energy into the side lobes. The outermost cell is dropped: that half
    cell has no colocated data, and deep inside the PML the field has long since decayed to
    negligible.
    """
    apod = getattr(mon, "apodization", None)
    if apod is None:                     # a td monitor always has this field; this is a backstop
        apod = td.ApodizationSpec()
    out = []
    lo, hi = mon.bounds
    thick = [i for i in range(3) if mon.size[i] > 0]
    if not thick:
        raise NotImplementedError(f"projection monitor {mon.name} has zero thickness on all three axes")
    if len(thick) < 3:
        # A plane monitor is itself the near-field surface. The outward normal comes from the
        # monitor's own ``normal_dir``; the KSpace and Cartesian plane monitors have that field
        # while the Angle one does not, where it defaults to "+".
        ax = next(i for i in range(3) if mon.size[i] == 0)
        proj_axis = getattr(mon, "proj_axis", None)
        if proj_axis is not None and int(proj_axis) != ax:
            raise NotImplementedError(
                f"projection monitor {mon.name} has proj_axis={proj_axis}, "
                f"which disagrees with its zero-thickness axis {ax}")
        sign = str(getattr(mon, "normal_dir", "+"))
        center, size = list(mon.center), list(mon.size)
        for a in range(3):
            if a == ax:
                continue
            b = np.asarray(getattr(sim.grid.boundaries, "xyz"[a]), dtype=np.float64)
            # When cutting to the grid, pull in by a further 1e-6 um. The projector samples at the
            # corner of a surface monitor, and if that corner lands exactly on a grid boundary, the
            # largest coordinate of the surface field data (the same boundary in a different
            # floating-point representation) may be 1e-15 smaller, at which point scipy's hard-edged
            # interpolation raises "above the interpolation range".
            a_lo = max(float(lo[a]), float(b[1]) + 1e-6)
            a_hi = min(float(hi[a]), float(b[-2]) - 1e-6)
            center[a] = 0.5 * (a_lo + a_hi)
            size[a] = a_hi - a_lo
        out.append((_face_monitor(
            center, size, mon.freqs, apod,
            name_template.format(name=mon.name, axis="xyz"[ax], sign=sign)),
            sign, ax))
        return out
    # exclude_surfaces: one example excludes the z- surface of a 3D box because that surface cuts
    # through the substrate, and Tidy3D's FieldProjector requires every near-field surface to lie in
    # a uniform medium, raising in from_near_field_monitors otherwise. An excluded surface gets no
    # monitor and takes no part in the projection.
    excl = set(getattr(mon, "exclude_surfaces", None) or ())
    for ax_name, sign in FACES:
        if f"{ax_name}{sign}" in excl:
            continue
        ax = "xyz".index(ax_name)
        center = list(mon.center)
        size = list(mon.size)
        # Take the surface position from the **monitor's own bounds**, not from center +- size/2:
        # the two differ in floating-point representation (measured 0.285 against
        # 0.2849999999999995), and the projector looks up points by bounds, where a difference of
        # 5e-16 already raises "outside the interpolation range".
        center[ax] = float(hi[ax] if sign == "+" else lo[ax])
        size[ax] = 0.0
        out.append((_face_monitor(
            center, size, mon.freqs, apod,
            name_template.format(name=mon.name, axis=ax_name, sign=sign)),
            sign, ax))
    return out


def _coords(sim: td.Simulation, comp: str, origin, box) -> dict:
    """Yee coordinates of one component inside the index box, in um, keyed by ``x``/``y``/``z``."""
    g = sim.grid.yee.grid_dict[comp]
    out = {}
    for a, ax in enumerate("xyz"):
        c = np.asarray(getattr(g, ax), dtype=np.float64)
        out[ax] = c[origin[a]: origin[a] + box[a]]
    return out


def colocation_points(sim: td.Simulation, mon: td.FieldMonitor) -> dict:
    """Colocation points on a monitor plane: the plane itself along the normal, and the **primal
    boundaries** along the tangential axes.

    This is Tidy3D's ``colocate=True`` convention, and it is what ``FieldProjector`` expects: it
    treats all six components as values at the same set of points when doing the surface integral.
    Handing it data at native Yee positions makes it interpolate once more and then look for a
    point we do not have; measured, it raised "0.285 outside the interpolation range", which is
    exactly the nominal position of the surface.
    """
    nrm = next(k for k in range(3) if mon.size[k] == 0)
    # A ``size=td.inf`` monitor has bounds of +-inf, and using those directly to pick colocation
    # points would reach infinity, so such an axis falls back to the simulation domain edge.
    # **Finite bounds are used as they are**: a projection surface monitor deliberately extends into
    # the PML (see surface_monitors), and cutting at the domain edge would undo that.
    ml, mh = mon.bounds
    sl, sh = sim.bounds
    lo = [float(ml[a]) if np.isfinite(ml[a]) else float(sl[a]) for a in range(3)]
    hi = [float(mh[a]) if np.isfinite(mh[a]) else float(sh[a]) for a in range(3)]
    out = {}
    for a, ax in enumerate("xyz"):
        if a == nrm:
            c = float(mon.center[a])
            # Snap only on mode planes: snapping a far-field projection surface makes the projection
            # phase inconsistent with the monitor geometry (it broke the analytic dipole in
            # test_kspace), and the far-field side is evaluated separately.
            if knobs.env("COLOC_SNAP") == "1" and str(getattr(mon, "name", "")).endswith(modes.PLANE_SUFFIX):
                # tidy3d's convention: the data of a zero-thickness surface lands on the **nearest
                # primal boundary** (measured from sim.discretize(mon)), not on the nominal center.
                # When the nominal center falls between two boundary layers, linear interpolation
                # mixes the propagation phase and the amplitude comes out low; a straight-waveguide
                # probe measured |amp|^2 2.5% down.
                # tidy3d's rule is to take the boundary not greater than the center (floor), which
                # is colocate.snap_floor, shared by three call sites. The tolerance is 1e-12 um
                # here; nb.backend and td_readout pass 1e-15 m. Each passes its own and they are
                # deliberately not unified.
                c = colocate.snap_floor(getattr(sim.grid.boundaries, ax), c, 1e-12)
            out[ax] = np.array([c])
            continue
        b = np.asarray(getattr(sim.grid.boundaries, ax), dtype=np.float64)
        sel = b[(b >= lo[a] - COORD_ATOL) & (b <= hi[a] + COORD_ATOL)]
        # **Inset both ends slightly, while keeping the endpoints.** The projector samples at the
        # monitor corners, and the floating-point representation of a "corner" differs between the
        # parent monitor and the surface monitor (measured 0.285 against 0.2849999999999995).
        # scipy interpolates with hard edges, so even 5e-16 raises "below the interpolation range".
        # Insetting by COORD_ATOL (1e-9 um = 1e-15 m) is physically negligible, and the two-cell
        # margin of the index box easily absorbs it. When an endpoint coincides with a grid boundary
        # (an infinite monitor cut to b[1] or b[-2]), that boundary point is dropped and only the
        # inset endpoint kept: otherwise the two points differ by 1e-9, become exact duplicates once
        # field_data clamps them back into the index box, and xarray raises InvalidIndexError
        # (reindexing is only valid with uniquely valued Index objects).
        ends = np.array([float(lo[a]) - COORD_ATOL, float(hi[a]) + COORD_ATOL])
        keep = sel[(np.abs(sel - float(lo[a])) > 1e-6) & (np.abs(sel - float(hi[a])) > 1e-6)]
        out[ax] = np.unique(np.concatenate([keep, ends]))
    return out


def _periodic_periods(sim: td.Simulation) -> dict:
    """Periodic axes mapped to their period length in um; a non-periodic axis is absent from the dict.

    Bloch counts as periodic. The phase difference is not carried yet: only a k=0 whole-period wrap
    on a full-domain mode plane reaches here, and the Bloch phase exp(i*k*L) will be added when it
    is actually needed. The period of the domain is the span of ``grid.boundaries``.
    """
    out = {}
    for ax in "xyz":
        bs = getattr(sim.boundary_spec, ax)
        kinds = {type(bs.plus).__name__, type(bs.minus).__name__}
        if kinds & {"Periodic", "BlochBoundary"}:
            b = np.asarray(getattr(sim.grid.boundaries, ax), dtype=np.float64)
            out[ax] = float(b[-1] - b[0])
    return out


def _wrap_periodic(full, ax: str, period: float, want: np.ndarray):
    """Extend the data by one periodic image along a periodic axis, so the colocation points fall
    inside the interpolation range.

    On a periodic domain a cell-centered Yee component is missing its "upper boundary" sample, since
    F(hi) = F(lo) is the same point on the torus. The colocation points, however, include both
    domain edges, so the interpolation range is short by at most one cell. Appending the first
    sample at +period to the end and the last sample at -period to the front closes the gap exactly.
    Physically that is periodic continuation and introduces no new information. It only runs when
    ``want`` genuinely exceeds the current range, to avoid a pointless copy.
    """
    from tidy3d.components.data.data_array import ScalarFieldDataArray

    co = np.asarray(full.coords[ax], dtype=np.float64)
    if co.size < 2:
        return full
    pieces, coords = [full], [co]
    # Upper end: append samples from the start at +period until strictly past want.max. (want
    # already carries a +-COORD_ATOL margin; no atol is subtracted here, so the interpolation range
    # genuinely covers it.)
    j = 0
    while float(coords[-1][-1]) < want.max() and j < co.size:
        pieces.append(full.isel(**{ax: [j]}))
        coords.append(co[j: j + 1] + period)
        j += 1
    # Lower end: prepend samples from the end at -period until strictly past want.min
    j = co.size - 1
    while float(coords[0][0]) > want.min() and j >= 0:
        pieces.insert(0, full.isel(**{ax: [j]}))
        coords.insert(0, co[j: j + 1] - period)
        j -= 1
    if len(pieces) == 1:
        return full
    import xarray as xr
    merged = xr.concat(pieces, dim=ax)
    new_co = np.concatenate(coords)
    merged = merged.assign_coords({ax: new_co})
    return ScalarFieldDataArray(merged.data, coords={
        k: (new_co if k == ax else np.asarray(full.coords[k]))
        for k in full.dims})


def _clamp_target(sim: td.Simulation, mon, nm: str, ax: str, want: np.ndarray,
                  got: np.ndarray) -> tuple[np.ndarray, tuple | None]:
    """Clamp the colocation targets of one axis to the index-box range; returns
    ``(coordinates to write back into tgt[ax], dedup)``.

    The endpoint tolerance is raised to 1e-6 um: after an infinite monitor is cut to a grid
    boundary, the colocation endpoint and the index-box endpoint differ by about 1e-9, right on the
    old tolerance. An endpoint within tolerance is clamped back into the index-box range so that
    xarray.interp does not return NaN at the edge.

    A monitor can stick out of the computational domain (one coupler has its mode plane reaching
    y=30.75 while the simulation domain stops at 30.0). tidy3d gives colocation points from the
    monitor's own bounds, not truncated to the grid. If we have already reached the grid boundary on
    the short side, there is no field outside the domain anyway and clamping is correct; otherwise
    cells really are missing and it still raises.

    Several out-of-domain points clamping to the same boundary value produce duplicate coordinates,
    which xarray.interp does not accept. ``dedup`` then supplies ``(unique coordinates, inverse
    index)``: interpolate on the unique coordinates and expand back by index afterwards, with the
    final coordinates still the ones the monitor asked for (here what is written back is ``want``
    itself).
    """
    if want.min() < got.min() - COORD_EDGE_TOL or want.max() > got.max() + COORD_EDGE_TOL:
        _k = "xyz".index(ax)
        _dlo, _dhi = float(sim.bounds[0][_k]), float(sim.bounds[1][_k])
        # A gap because the colocation point ran outside the simulation domain is clamped back;
        # there is no field out there. A gap because cells are missing inside the domain is a real
        # bug.
        _lo_ok = want.min() >= got.min() - COORD_EDGE_TOL or want.min() < _dlo + COORD_EDGE_TOL
        _hi_ok = want.max() <= got.max() + COORD_EDGE_TOL or want.max() > _dhi - COORD_EDGE_TOL
        if _lo_ok and _hi_ok:
            print(f"[openem] monitor {mon.name} sticks out of the computational domain along {ax}: "
                  f"colocation wants [{want.min():.4f}, {want.max():.4f}] but the grid only reaches "
                  f"[{got.min():.4f}, {got.max():.4f}]; points outside are clamped to the boundary "
                  "value", flush=True)
        else:
            raise ValueError(
                f"{mon.name}/{nm}: the index box does not cover the colocation points along {ax}: "
                f"wants [{want.min():.9f}, {want.max():.9f}], "
                f"has [{got.min():.9f}, {got.max():.9f}]"
                f" (got={got.size} points, want={want.size} points)")
    _clipped = np.clip(want, got.min(), got.max())
    _u, _inv = np.unique(_clipped, return_inverse=True)
    if _u.size != _clipped.size:
        return want, (_u, np.asarray(_inv).ravel())
    return _clipped, None


def _interp_component(full, tgt: dict, flat: dict, dedup: dict) -> np.ndarray:
    """Interpolate one component onto the colocation points.

    A deduplicated axis is interpolated on its unique coordinates and then expanded by index; a flat
    axis is not interpolated at all but broadcast along.
    """
    interp = {a: (dedup[a][0] if a in dedup else v)
              for a, v in tgt.items() if a not in flat}
    arr2 = np.asarray(full.interp(**interp).data if interp else full.data)
    for a, (_u, _inv) in dedup.items():
        arr2 = np.take(arr2, _inv, axis="xyz".index(a))
    for a, want in flat.items():
        k = "xyz".index(a)
        if arr2.shape[k] != want.size:
            arr2 = np.repeat(arr2, want.size, axis=k)
    return arr2


def field_data(sim: td.Simulation, mon: td.FieldMonitor, spec, phasors: np.ndarray,
               scale: float = FAR_FIELD_SCALE, norm=None):
    """Our raw phasors, colocated onto the monitor plane as a ``td.FieldData``.

    ``phasors`` has shape ``(6, nf, ni, nj, nk)`` with components in the order Ex, Ey, Ez, Hx, Hy,
    Hz, **at their native Yee positions**. Two things happen here: each component is labelled with
    its own Yee coordinates, and then everything is interpolated onto
    :func:`colocation_points`. The three-cell neighbourhood the index box keeps on a zero-thickness
    axis exists for exactly this step (``monitors._monitor_box``).

    ``scale`` defaults to :data:`FAR_FIELD_SCALE` (that is, ``UM``), which is the factor the
    **far-field transform** needs and is not part of the field normalization itself; the mode
    overlap path must pass ``1.0``. The check is absolute: with a unit-power mode and a ModeSource
    of ``amplitude=1``, the downstream mode amplitude must be **1**. Carrying UM through measures
    1.0031e-06 instead.

    Note that ``tgt`` is modified in place inside the component loop, carrying the clamped
    coordinates over to the next component. That is the historical behaviour;
    :func:`_clamp_target` only computes and this function writes back, which preserves it.
    """
    from tidy3d.components.data.data_array import ScalarFieldDataArray

    names = ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz")
    freqs = np.asarray(mon.freqs, dtype=np.float64)
    tgt = colocation_points(sim, mon)
    if norm is None:
        norm = normalization(sim, freqs)
    else:
        # The caller supplied the same divisor as FieldData uses (1/field_norm), which already
        # contains the source-type factor, so FAR_FIELD_SCALE is not applied on top. See the note in
        # nb.backend._kspace_data.
        norm = np.asarray(norm, dtype=np.complex128).reshape(-1)
        scale = 1.0
    periods = _periodic_periods(sim)
    fields = {}
    for i, nm in enumerate(names):
        co = _coords(sim, nm, spec.origin, spec.box)
        arr = np.moveaxis(phasors[i], 0, -1) * norm * scale   # (nf,ni,nj,nk) -> (ni,nj,nk,nf)
        full = ScalarFieldDataArray(
            arr, coords={"x": co["x"], "y": co["y"], "z": co["z"], "f": freqs})
        # A full-domain mode plane on a periodic boundary: a cell-centered component is missing a
        # sample at the upper boundary (F(hi) = F(lo)) while the colocation points include both
        # domain edges. One periodic image closes the interpolation range.
        for ax, period in periods.items():
            if ax in tgt:
                full = _wrap_periodic(full, ax, period, tgt[ax])
        flat = {}
        dedup: dict = {}          # axes with duplicate coordinates after clamping: (unique coords, inverse index)
        for ax, want in tgt.items():
            got = np.asarray(full.coords[ax], dtype=np.float64)
            if got.size == 1:
                # **The flat axis of a 2D simulation**: only one sample fits, the field is constant
                # along it, and it does not matter which edge the sample sits on. So do not
                # interpolate; broadcast along that axis onto the colocation points.
                # (tidy3d's colocation_boundaries gives more than one target point on a flat axis;
                # measured want=2.) A genuine 3D zero-thickness axis never reaches here, because
                # monitors._monitor_box gives a three-cell neighbourhood.
                flat[ax] = want
                continue
            tgt[ax], dd = _clamp_target(sim, mon, nm, ax, want, got)
            if dd is not None:
                dedup[ax] = dd
        arr2 = _interp_component(full, tgt, flat, dedup)
        fields[nm] = ScalarFieldDataArray(arr2, coords=dict(**tgt, f=freqs))
    return td.FieldData(monitor=mon, symmetry=(0, 0, 0), symmetry_center=(0, 0, 0),
                        grid_expanded=sim.discretize(mon), **fields)


def _window_no_taper_on_zero(orig):
    """A wrapper around tidy3d's ``window_function``: ``window_size[dim] == 0`` means "no window on
    this dimension" and returns all ones.

    tidy3d's `window_parameters` fills window_size only for the two in-plane dimensions and leaves
    the normal at 0, while `_apply_window_to_currents` loops over all three of xyz. On the normal
    pass it takes 0 as the denominator, so any non-zero normal coordinate gives exp(-inf) = 0 and
    the currents over the whole surface are wiped out. That is how four windowed monitors in one
    example came back all zero. Every other case calls tidy3d's implementation unchanged (``orig``,
    captured in the closure).
    """
    def window_function(points, window_size, window_minus, window_plus, dim):
        if not window_size[dim]:
            return np.ones_like(points)
        return orig(points=points, window_size=window_size,
                    window_minus=window_minus, window_plus=window_plus, dim=dim)
    return window_function


@contextlib.contextmanager
def _patched_window(cls):
    """Within the scope, replace ``cls.window_function`` with :func:`_window_no_taper_on_zero`, and
    restore it on exit.
    """
    orig = cls.window_function
    cls.window_function = staticmethod(_window_no_taper_on_zero(orig))
    try:
        yield
    finally:
        cls.window_function = staticmethod(orig)


def project(sim: td.Simulation, proj_mon, faces: list, datas: list):
    """Run one near-to-far transform.

    Args:
        faces: the return value of :func:`surface_monitors`.
        datas: one ``td.FieldData`` per entry of ``faces``.
    """
    from tidy3d.components.field_projection import FieldProjector

    sim_data = td.SimulationData(simulation=sim.updated_copy(
        monitors=[f[0] for f in faces]), data=tuple(datas))
    fp = FieldProjector.from_near_field_monitors(
        sim_data=sim_data, near_monitors=[f[0] for f in faces],
        normal_dirs=[f[1] for f in faces],
        origin=tuple(proj_mon.custom_origin) if proj_mon.custom_origin is not None
        else None)
    # When windowing, work around tidy3d tapering the normal dimension too (see
    # _window_no_taper_on_zero)
    if getattr(proj_mon, "window_size", (0, 0)) != (0, 0):
        with _patched_window(type(proj_mon)):
            return fp.project_fields(proj_mon, verbose=False)
    return fp.project_fields(proj_mon, verbose=False)
