# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""phasor -> Poynting -> plane flux.

Any normal axis: the two transverse axes ``(t1, t2)`` are taken in ascending order,

    P(f) = 1/2 * Re Σ [ E_t1·conj(H_t2)·dA_1 - E_t2·conj(H_t1)·dA_2 ]

(for a z normal this is ``Ex·Hy* - Ey·Hx*``). At axis=1 the ascending order breaks the cyclic
order of (x, y, z), so P is **-S_y**; the negation is done in :func:`monitor_flux`.

The area weights take the dual face element at each component's own staggered position (``E_t1``
sits at (ctr_1, edge_2), so ``dA_1 = dl1_primal * dl2_dual``; ``E_t2`` is the other way round).
Modelled on ``_resolve_face_area_weights`` in FDTDX's ``poynting_flux.py``.

Takes bare arrays only and does not import solver, so the solver can call it inside the main loop
as a convergence criterion.
"""

from __future__ import annotations

import numpy as np

from openem import knobs

from openem.grid import ABSORBING, MIRROR_SIGN, Grid, transverse_axes

#: comp order of the phasor arrays (matching accumulate_dft in kernels/source_dft.cu)
EX, EY, HX, HY = 0, 1, 2, 3


#: Equivalent thickness [m] of the flat axis in a 2D simulation. In 2D, Tidy3D reports a "per µm"
#: quantity: the flat axis is **not** multiplied by the actual thickness of that cell. The
#: criterion is unambiguous: for a unit-power ModeSource injected into a 2D straight waveguide,
#: both the downstream ``|a_+|²`` and the flux should be 1; without this the flux comes out
#: exactly dz(µm) times too small (measured at dz=0.05 µm: flux/|a_+|² = 0.0500, while
#: |a_+|² = 0.9991 is itself correct).
#: This never surfaced before: the validation cases are all ratios (T = flux / reference flux),
#: where dz cancels top and bottom.
FLAT_AXIS_M = 1.0e-6


def axis_dl(ax) -> tuple[np.ndarray, np.ndarray]:
    """``(dl, dl_dual)`` of one axis, with a 2D flat axis replaced by :data:`FLAT_AXIS_M`.

    Everything that computes a face element or a transverse extent goes through it
    (:func:`area_weights`, scene/modes._yee_power, the diffraction orders in nb.backend).
    Renamed from the private ``_axis_dl`` to a public name on 2026-09-10.
    """
    if ax.is_flat:          # only sim.size==0 axes; a 1-cell size!=0 axis keeps its real dl
        one = np.array([FLAT_AXIS_M], dtype=np.float64)
        return one, one
    return ax.dl, ax.dl_dual


#: Old name, kept for the places that still refer to ``flux._axis_dl``.
_axis_dl = axis_dl


def area_weights(grid: Grid, axis: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """The face element of each of the two Poynting terms on a plane. ``axis`` is the normal
    axis.

    The flat axis of a 2D simulation uses :data:`FLAT_AXIS_M`, not the actual thickness of that
    cell.
    """
    t1, t2 = transverse_axes(axis)
    (dlx, dldx), (dly, dldy) = (axis_dl(grid.axes[t1]),
                                axis_dl(grid.axes[t2]))
    dA_x = dlx[:, None] * dldy[None, :]              # the Ex·Hy* term
    dA_y = dldx[:, None] * dly[None, :]              # the Ey·Hx* term
    return dA_x, dA_y


def plane_flux(data: np.ndarray, grid: Grid,
               transverse: tuple[tuple[int, int], tuple[int, int]] | None = None,
               axis: int = 2) -> np.ndarray:
    """Plane flux ``(nf,)`` in the ascending-transverse-axis convention:
    ``½ Re Σ [E_t1·H_t2* - E_t2·H_t1*]``.

    For axis=0/2 this is +S_x / +S_z (positive = energy flowing along the positive normal
    direction, the same as ``normal_dir="+"``); for axis=1 it is **-S_y**, and the negation is
    done in :func:`monitor_flux` (see the module header). Each cell of the two transverse axes is
    summed once, and a periodic end point is not counted twice.

    Args:
        data: ``(4, nf, n1, n2)`` complex phasors, comp order E_t1, E_t2, H_t1, H_t2.
        transverse: closed-interval indices of the two transverse axes,
            ``((j0, j1), (k0, k1))``; ``None`` = the whole plane. It **must** be given when the
            monitor does not span the full transverse extent, otherwise the whole plane gets
            integrated.
        axis: normal axis of the plane; it decides which two axes ``transverse`` refers to.
    """
    dA_x, dA_y = area_weights(grid, axis)
    if transverse is not None:
        mask = np.zeros(dA_x.shape, dtype=np.float64)
        (j0, j1), (k0, k1) = transverse
        mask[j0:j1 + 1, k0:k1 + 1] = 1.0
        dA_x = dA_x * mask
        dA_y = dA_y * mask
    ex, ey, hx, hy = data[EX], data[EY], data[HX], data[HY]
    term_x = np.real(ex * np.conj(hy)) * dA_x[None, :, :]
    term_y = np.real(ey * np.conj(hx)) * dA_y[None, :, :]
    return 0.5 * (term_x.sum(axis=(1, 2)) - term_y.sum(axis=(1, 2)))


def normalized(
    with_structure: dict,
    reference: dict,
    grid: Grid,
    incident_monitor: str,
) -> dict[str, np.ndarray]:
    """Normalize the flux to "incident power = 1", taking the incident power from a **reference
    run without the structure**.

    That way the absolute injection amplitude of the plane wave cancels out entirely and there is
    no need to reverse-engineer Tidy3D's amplitude convention.

    Args:
        incident_monitor: which monitor's reference flux to use as the incident power (the one on
            the transmission side).
    """
    p_inc = plane_flux(reference[incident_monitor].data, grid)
    if np.any(p_inc <= 0):
        raise ValueError(
            f"the incident flux of the reference run is not positive: {p_inc}. The plane wave "
            "may be injected in the wrong direction, or its amplitude is zero"
        )
    return {n: plane_flux(ph.data, grid) / p_inc for n, ph in with_structure.items()}


class _PlaneGeomSpec:
    """The minimal description :func:`flux_time_geometry` needs: the frequency-domain plane sits
    exactly on a primal boundary (``frac=0``) and the sign does not carry normal_dir (matching
    what :func:`monitor_flux` returns: +S_axis, with the cyclic-order fix for axis=1 done by
    flux_time_geometry itself)."""

    frac = 0.0
    normal_dir = +1

    def __init__(self, mon):
        self.name, self.axis = mon.name, mon.axis
        self.plane_index, self.t_bounds = mon.plane_index, mon.t_bounds


def _center_to_edge(a, g: dict, comp: str, t_ax: int):
    """Given a tangential component ``a`` ``(n1, n2)`` already taken on the primal plane, do the
    center->edge interpolation along axis ``t_ax``. The 2D and 3D versions share this code.

    A folded or wrapped axis goes through the ``idx_lo`` / ``wlo_c`` channel: the low neighbour of
    sample 0 falls outside the domain, so the mirrored cell / wrap-around cell is read instead,
    with the mirror parity already folded into element 0 of the weights (see
    flux_time_geometry). Only slicing and multiply-add, so it is **written identically for numpy
    and cupy arrays**.
    """
    fl = g.get("fold_lo") or (False, False)
    if t_ax == 0:
        if fl[0]:
            return g["wlo_c"][0][comp] * a[g["idx_lo"][0], :] + g["whi"][0] * a
        return g["wlo"][0] * a[:-1, :] + g["whi"][0] * a[1:, :]
    if fl[1]:
        return g["wlo_c"][1][comp] * a[:, g["idx_lo"][1]] + g["whi"][1] * a
    return g["wlo"][1] * a[:, :-1] + g["whi"][1] * a[:, 1:]


def _colocate_2d(a: np.ndarray, g: dict, comp: str) -> np.ndarray:
    """The 2D version of :func:`colocated_component`: ``a`` is a tangential component
    ``(n1, n2)`` already taken on the primal plane (the normal colocation is done by the plane DFT
    kernel), so only the tangential center->edge interpolation is left here; sl/wlo/whi/fold/wrap
    use the same geometry as the 3D version."""
    t_ax = 0 if comp in ("Eu", "Hv") else 1
    sl = [g["sl_edge"][0], g["sl_edge"][1]]
    sl[t_ax] = g["sl_cell"][t_ax]
    return _center_to_edge(a[sl[0], sl[1]], g, comp, t_ax)


def _flux_colocated_flat(data: np.ndarray, grid: Grid, mon, u: int, u_first: bool) -> np.ndarray:
    """Colocated convention for a 2D simulation (one transverse axis is flat): center->edge
    colocation plus trapezoidal face elements along the non-flat axis ``u``, while the single
    sample of the flat axis is used as it is, with width :data:`FLAT_AXIS_M` (the same convention
    as :func:`axis_dl` and the old quadrature).
    ``u_first``: u is the first of the ascending transverse pair (so E_u/H_v sit at cell centres
    along u and need interpolating along u, while E_v/H_u sit on the boundary along u and stay).
    Boundaries: a sample right at the edge may only be dropped at an absorbing end of the u axis;
    a symmetry/periodic/PEC end fails closed and the caller falls back to the old quadrature.
    """
    axu = grid.axes[u]
    lo, hi = mon.t_bounds[0 if u_first else 1]
    i0, i1, w = colocated_axis_weights(axu.edges, lo, hi)
    for end, need in (("lo", i0 < 1), ("hi", i1 > axu.n - 1)):
        if need and getattr(axu, f"boundary_{end}") not in ABSORBING:
            raise ValueError(
                f"monitor {mon.name}: the samples along {'xyz'[u]} sit right at the grid edge, "
                f"and the {end} end is not an absorbing boundary")
    if i0 < 1:
        w = w[1 - i0:]
        i0 = 1
    if i1 > axu.n - 1:
        w = w[:axu.n - 1 - i1]
        i1 = axu.n - 1
    eg, cen = axu.edges[i0:i1 + 1], axu.centers
    wlo = (cen[i0:i1 + 1] - eg) / (cen[i0:i1 + 1] - cen[i0 - 1:i1])
    whi = 1.0 - wlo

    def along_u(a):            # a: (n_u, 1) or (1, n_u), centres along u -> edges i0..i1
        v = a if u_first else a.T
        v = v[:, 0]
        return wlo * v[i0 - 1:i1] + whi * v[i0:i1 + 1]

    def on_edge(a):            # already on the boundary along u, take samples i0..i1
        v = a if u_first else a.T
        return v[i0:i1 + 1, 0]

    sign = -1.0 if mon.axis == 1 else 1.0
    dS = w * FLAT_AXIS_M
    nf = data.shape[1]
    out = np.empty(nf, dtype=np.float64)
    for f in range(nf):
        eu, ev, hu, hv = data[EX, f], data[EY, f], data[HX, f], data[HY, f]
        if u_first:      # u = t1: E_u and H_v sit at cell centres along u
            Eu, Hv, Ev, Hu = along_u(eu), along_u(hv), on_edge(ev), on_edge(hu)
        else:            # u = t2: E_v and H_u sit at cell centres along u
            Ev, Hu, Eu, Hv = along_u(ev), along_u(hu), on_edge(eu), on_edge(hv)
        out[f] = sign * 0.5 * float(np.real(np.sum((Eu * np.conj(Hv) - Ev * np.conj(Hu)) * dS)))
    return out


def monitor_flux_colocated(data: np.ndarray, grid: Grid, mon) -> np.ndarray:
    """Frequency-domain plane flux ``(nf,)`` in the colocated-integration convention (the same as
    tidy3d's ``FieldData.flux``): the four tangential components are colocated onto the primal
    boundary points, the trapezoidal face elements are clipped exactly at the monitor bounds, and
    ``½ Re Σ dS·(E_u·H_v* - E_v·H_u*)``. Returns the same convention as :func:`monitor_flux`
    (+S_axis).

    ``data`` has to be the **whole plane** ``(4, nf, n1, n2)`` (which is how the plane DFT kernel
    stores it), with ``mon.t_bounds`` filled in. Old results (ours.npz) satisfy the first
    condition, so re-running the scoring is enough.
    """
    if getattr(mon, "t_bounds", None) is None:
        raise ValueError(
            f"monitor {mon.name} has no t_bounds, the colocated convention is unavailable")
    t1, t2 = transverse_axes(mon.axis)
    f1, f2 = grid.axes[t1].is_flat, grid.axes[t2].is_flat
    if f1 and f2:
        raise ValueError(
            f"both transverse axes of monitor {mon.name} are flat, "
            "the colocated convention is meaningless")
    if f1 or f2:
        return _flux_colocated_flat(data, grid, mon, u=(t2 if f1 else t1), u_first=f2)
    g = flux_time_geometry(grid, _PlaneGeomSpec(mon))
    nf = data.shape[1]
    out = np.empty(nf, dtype=np.float64)
    for f in range(nf):
        eu = _colocate_2d(data[EX, f], g, "Eu")
        ev = _colocate_2d(data[EY, f], g, "Ev")
        hu = _colocate_2d(data[HX, f], g, "Hu")
        hv = _colocate_2d(data[HY, f], g, "Hv")
        out[f] = 0.5 * float(np.real(np.sum((eu * np.conj(hv) - ev * np.conj(hu)) * g["dS"])))
    return out


def flux_quadrature() -> str:
    """The current frequency-domain plane flux convention: ``"colocated"`` / ``"yee"``
    (environment variable ``OPENEM_FLUX_COLOCATED``, default 1; set it to 0 to fall back to the
    in-place Yee integration, for comparison)."""
    return "colocated" if knobs.env("FLUX_COLOCATED") == "1" else "yee"   # default on (2026-09-06)


def monitor_flux(data: np.ndarray, grid: Grid, mon) -> np.ndarray:
    """Flux over the transverse range and normal axis that
    :class:`~openem.model.FluxMonitor` declares for itself.

    :func:`flux_quadrature` picks the convention: ``colocated`` (the same as tidy3d, 2026-09-06)
    or ``yee`` (the old in-place staggered integration, which differs from |mode amplitude|² by
    7.7% at 16 cells per wavelength). A missing ``t_bounds`` (an old scene) falls back to ``yee``.

    **Use this, not** :func:`plane_flux`, **when comparing against the reference outputs**: the
    latter integrates the whole plane by default, so a monitor that does not span the full
    transverse extent would also count the power of a neighbouring waveguide.
    """
    # ``plane_flux`` computes ``Re(E_t1 H_t2*) - Re(E_t2 H_t1*)`` with (t1,t2) in **ascending**
    # order.
    # axis=0: E_y H_z* - E_z H_y* = +S_x   ok
    # axis=2: E_x H_y* - E_y H_x* = +S_z   ok
    # axis=1: E_x H_z* - E_z H_x* = **-S_y**: the ascending order breaks the cyclic order of
    #         (x, y, z), so a y normal has to be negated. **Verified**: for the closed box of
    #         PlasmonicYagiUda (an Ey dipole, with the y faces carrying 16% of the total power)
    #         the outward sum over the six faces equals the reference value 1.3e-02; with the
    #         sign flipped it is off by -31%. There is also the symmetry criterion in
    #         tests/test_boxflux.py.
    if (flux_quadrature() == "colocated" and getattr(mon, "t_bounds", None) is not None
            and data.ndim == 4 and data.shape[2:] == tuple(
                grid.axes[t].n for t in transverse_axes(mon.axis))):
        try:
            return monitor_flux_colocated(data, grid, mon)
        except (ValueError, NotImplementedError) as _e:
            # Geometric refusals such as samples right against a non-absorbing boundary (a box
            # face against the wall, a periodic or symmetry axis): fall back to the in-place Yee
            # convention. Such faces are a small minority and are usually auxiliary faces outside
            # R/T anyway (the box faces of BullseyeCavityPSO, 2026-09-07).
            print(f"[openem] monitor {mon.name} cannot use the colocated convention "
                  f"({str(_e)[:70]}...), falling back to the in-place Yee integration")
    sign = -1.0 if mon.axis == 1 else 1.0
    return sign * plane_flux(data, grid, transverse=mon.transverse, axis=mon.axis)


def colocated_axis_weights(edges: np.ndarray, lo: float,
                           hi: float) -> tuple[int, int, np.ndarray]:
    """Sample range and trapezoidal weights ``(i0, i1, w)`` of the colocated integration on one
    axis.

    The samples are the primal boundary points ``edges[i0..i1]`` (both ends included); the
    integration cell of each sample is the dual cell spanned by the midpoints of the neighbouring
    samples, clipped exactly at the monitor bounds ``[lo, hi]`` (Tidy3D's
    ``yee_areas.colocated_widths_1d``, measured to be bitwise identical to the client's
    ``_diff_area``). Only the contiguous run of non-zero weights is returned.

    Raises:
        ValueError: the monitor intersects the dual cell of no sample on that axis.
    """
    mid = 0.5 * (edges[:-1] + edges[1:])
    lo_e = np.concatenate([edges[:1], mid])       # lower cell edge of sample i (ends: the sample)
    hi_e = np.concatenate([mid, edges[-1:]])
    w = np.maximum(np.clip(hi_e, lo, hi) - np.clip(lo_e, lo, hi), 0.0)
    nz = np.flatnonzero(w > 0.0)
    if nz.size == 0:
        raise ValueError(
            f"the monitor range [{lo}, {hi}] intersects no integration cell on this axis")
    return int(nz[0]), int(nz[-1]), w[nz[0]:nz[-1] + 1]


def _bcast(p: int) -> tuple[int, int]:
    """Shape that broadcasts the 1D weights of transverse axis ``p`` out to ``(n1, n2)``."""
    return (-1, 1) if p == 0 else (1, -1)


def _normal_planes(grid: Grid, mon) -> list[dict]:
    """The 1-2 primal planes that go into the combination along the normal.

    E sits exactly on the primal plane; H is interpolated onto that plane by coordinate from the
    cell centres on both sides (two cells starting at ``kh``, weights ``wh``). With ``frac=0``
    there is only one plane, and since ``km+1`` may be out of range, the term with weight 0 is
    skipped outright.

    Raises:
        ValueError: the plane sits right at the grid edge, so the H colocation cannot reach the
            cell centres on both sides; fail closed.
    """
    ax = mon.axis
    e = grid.axes[ax].edges
    c = grid.axes[ax].centers
    km, fr = mon.plane_index, float(mon.frac)
    planes = []
    for w_p, k in ((1.0 - fr, km), (fr, km + 1)):
        if w_p == 0.0:
            continue
        kh = k - 1
        if not 0 <= kh < c.size - 1:
            raise ValueError(
                f"monitor {mon.name}: plane {'xyz'[ax]}[{k}] sits right at the grid edge, so "
                f"the H colocation cannot reach the cell centres on both sides "
                f"(1..{c.size - 1} allowed)")
        whlo = (float(c[kh + 1]) - float(e[k])) / float(c[kh + 1] - c[kh])
        planes.append({"w": w_p, "ke": int(k), "kh": kh, "wh": (whlo, 1.0 - whlo)})
    return planes


def _wrap_channel(axt, wlo: np.ndarray, i1: int, p: int) -> tuple:
    """A fully spanned pure-periodic axis: the low neighbour of sample 0 wraps round to ctr[n-1]
    (parity factor +1). Returns ``(idx_lo, wlo_c)``, sharing the same two channels as a folded
    axis."""
    shape = _bcast(p)
    return (np.concatenate([[axt.n - 1], np.arange(i1)]).astype(np.int64),
            {k: wlo.reshape(shape).copy() for k in ("Eu", "Hv", "Ev", "Hu")})


def _fold_channel(axt, wlo: np.ndarray, i1: int, p: int, tax: tuple) -> tuple:
    """Low end of a folded axis = the symmetry plane: the low neighbour of sample 0 reads folded
    ctr 0, with the mirror parity folded into element 0 of the weights (one copy per component).
    Returns ``(idx_lo, wlo_c)``."""
    shape = _bcast(p)
    t = tax[p]
    s_ax = -1.0 if axt.boundary_lo == "SymmetryPEC" else 1.0
    t1_, t2_ = tax
    names = ((("Eu", "E" + "xyz"[t1_]), ("Hv", "H" + "xyz"[t2_]))
             if p == 0 else
             (("Ev", "E" + "xyz"[t2_]), ("Hu", "H" + "xyz"[t1_])))
    wc = {}
    for lbl, nm in names:
        arr = wlo.copy()
        arr[0] *= MIRROR_SIGN[nm][t] * s_ax   # mirror parity has its single home in grid.py
        wc[lbl] = arr.reshape(shape)
    return np.concatenate([[0], np.arange(i1)]).astype(np.int64), wc


def _tangential_axis(axt, mon, p: int, tax: tuple) -> dict:
    """The colocated convention on one transverse axis: the sample run ``sl_edge``, the cell run
    ``sl_cell``, the center->edge interpolation weights ``wlo``/``whi``, the trapezoidal weights
    ``trap``, plus the ``idx_lo``/``wlo_c`` channel for the two special cases where the low
    neighbour of sample 0 falls outside the domain (``idx_lo`` None = not taken).

    Raises:
        ValueError: a sample sits right at the grid edge and that end is not an absorbing
            boundary, so dropping the sample would not be safe.
    """
    t = tax[p]
    lo, hi = mon.t_bounds[p]
    i0, i1, w = colocated_axis_weights(axt.edges, lo, hi)
    # Low end of a folded axis = the symmetry plane: sample 0 sits exactly on the plane, and the
    # low-side neighbour of the center->edge step (full-domain ctr m-1) equals +-(folded ctr 0)
    # by parity, replayed in the weight domain so that no sample is lost; its dual cell has
    # already been cut in half by lo = the symmetry plane, and the folding factor of 2 restores
    # the whole cell exactly.
    fold_lo = (i0 == 0
               and axt.boundary_lo in ("SymmetryPEC", "SymmetryPMC"))
    # A fully spanned pure-periodic axis (AndersonLocalization: the flux plane spans the periodic
    # direction): sample n and sample 0 are the same degree of freedom, so the end weight is
    # merged into sample 0 and the low neighbouring cell centre wraps round to ctr[n-1] (shifted
    # by one period in coordinate), which puts the integration measure back on a whole period.
    # This reuses the idx_lo/wlo_c channel of fold_lo (parity factor +1, so the solver's GPU
    # packing needs no change).
    # Pure Periodic only: a Bloch wrap-around sample carries an exp(i k L) phase, and with no
    # case that needs it, it is not done.
    wrap_lo = (i0 == 0 and i1 == axt.n
               and axt.boundary_lo == "Periodic"
               and axt.boundary_hi == "Periodic")
    # Samples right at the grid edge (the monitor spans or overruns the simulation domain; all
    # six faces of NanobeamCavity overrun): at an absorbing end they are dropped, since that is
    # the deepest part of the PML where the outermost E_t is pinned to 0 anyway, and only half a
    # dual cell is lost per end. At a non-absorbing end (periodic/PEC/PMC) the field is not
    # small, so fail closed.
    for end, need in (("lo", i0 < 1 and not fold_lo and not wrap_lo),
                      ("hi", i1 > axt.n - 1 and not wrap_lo)):
        if need and getattr(axt, f"boundary_{end}") not in ABSORBING:
            raise ValueError(
                f"monitor {mon.name}: the samples along {'xyz'[t]}, [{i0}, {i1}], sit right at "
                f"the grid edge (1..{axt.n - 1} allowed) and the {end} end is not an absorbing "
                "boundary, so they cannot be dropped safely")
    if wrap_lo:
        w = np.concatenate([[w[0] + w[-1]], w[1:-1]])
        i1 = axt.n - 1
    if i0 < 1 and not fold_lo and not wrap_lo:
        w = w[1 - i0:]
        i0 = 1
    if i1 > axt.n - 1:
        w = w[:axt.n - 1 - i1]
        i1 = axt.n - 1
    eg, cen = axt.edges[i0:i1 + 1], axt.centers
    if fold_lo:
        # The coordinate of the low-side cell centre of sample 0 is the mirror 2·x_m - cen[0]
        # (on a symmetric grid wlo[0] comes out exactly 0.5); the matching field value is read in
        # colocated_component through idx_lo as folded ctr 0, with the parity factor folded into
        # element 0 of wlo_c.
        cm1 = np.concatenate([[2.0 * eg[0] - cen[0]], cen[0:i1]])
        wlo = (cen[i0:i1 + 1] - eg) / (cen[i0:i1 + 1] - cm1)
    elif wrap_lo:
        # the low-side cell centre of sample 0 is ctr[n-1] shifted back by one period
        L = float(axt.edges[-1] - axt.edges[0])
        cm1 = np.concatenate([[cen[-1] - L], cen[0:i1]])
        wlo = (cen[i0:i1 + 1] - eg) / (cen[i0:i1 + 1] - cm1)
    else:
        wlo = (cen[i0:i1 + 1] - eg) / (cen[i0:i1 + 1] - cen[i0 - 1:i1])
    shape = _bcast(p)
    idx_lo = wlo_c = None
    if wrap_lo:
        idx_lo, wlo_c = _wrap_channel(axt, wlo, i1, p)
    elif fold_lo:
        idx_lo, wlo_c = _fold_channel(axt, wlo, i1, p, tax)
    return {
        "sl_edge": slice(i0, i1 + 1),
        "sl_cell": slice(0 if (fold_lo or wrap_lo) else i0 - 1, i1 + 1),
        "wlo": wlo.reshape(shape),
        "whi": (1.0 - wlo).reshape(shape),
        "trap": w,
        "idx_lo": idx_lo,
        "wlo_c": wlo_c,
    }


def flux_time_geometry(grid: Grid, mon) -> dict:
    """The complete spatial convention of the time-domain instantaneous flux (colocated
    integration).

    The reference implementation's convention (``use_colocated_integration=True``, which a
    FluxTimeMonitor always uses): on each primal plane the four tangential components are
    colocated onto the primal boundary points and integrated,
    ``S_k = sign · Σ dS · (E_t1·H̄_t2 - E_t2·H̄_t1)``; when the plane is not on a Yee boundary it
    is the **flux scalar** that is interpolated linearly,
    ``s = (1-frac)·S_km + frac·S_{km+1}``, not the fields (settled by a self-consistency check
    against the reference outputs: for a steeply decaying field the two differ at the
    frac(1-frac)(ΔE)² level, 9% on the z faces of NanobeamCavity). sign is made up of
    ``normal_dir`` and the cyclic-order fix for axis=1; there is **no ½ here** as there is in the
    frequency domain, where it is the CW time-averaging factor.

    Returns:
        dict: ``axis / tax``; ``planes``, one entry for each of the 1-2 planes,
        ``{"w": interpolation weight, "ke": index of the primal plane holding E, "kh","wh": H
        interpolated by coordinate from centers[kh],centers[kh+1] onto edges[ke]}``; for each
        transverse axis p in (0,1) the sample run ``sl_edge[p]``, the cell run ``sl_cell[p]`` and
        the center->edge interpolation weights ``wlo[p]/whi[p]`` (already in broadcast shape);
        and the signed face element ``dS`` ((n1, n2) float64).

    Raises:
        ValueError: a plane or a monitor bound is too close to the grid edge for the
            interpolation to reach the neighbouring cell; fail closed rather than extrapolate
            silently.
    """
    ax = mon.axis
    t1, t2 = transverse_axes(ax)
    g: dict = {"axis": ax, "tax": (t1, t2)}
    g["planes"] = _normal_planes(grid, mon)

    # ---- transverse: sample range, trapezoidal face element, center->edge interpolation
    # weights (per axis, see _tangential_axis) ----
    g["sl_edge"], g["sl_cell"], g["wlo"], g["whi"] = [], [], [], []
    g["fold_lo"] = [False, False]
    g["idx_lo"] = [None, None]
    g["wlo_c"] = [None, None]
    trap = []
    for p, t in enumerate(g["tax"]):
        part = _tangential_axis(grid.axes[t], mon, p, g["tax"])
        for key in ("sl_edge", "sl_cell", "wlo", "whi"):
            g[key].append(part[key])
        trap.append(part["trap"])
        if part["idx_lo"] is not None:
            g["fold_lo"][p] = True
            g["idx_lo"][p] = part["idx_lo"]
            g["wlo_c"][p] = part["wlo_c"]
    sign = float(mon.normal_dir) * (-1.0 if ax == 1 else 1.0)
    g["dS"] = sign * np.outer(trap[0], trap[1])
    return g


def colocated_component(F, g: dict, comp: str, plane: int):
    """Colocate a field component onto the edge points of primal plane ``plane``, 2D
    ``(n1, n2)``.

    ``comp in {"Eu","Ev","Hu","Hv"}``: u/v are the ascending transverse axes (t1,t2). In the Yee
    layout E_d sits at the cell centre along d and H_d on the boundary along d, so Eu/Hv need the
    center->edge step along t1 and Ev/Hu along t2. Along the normal, the tangential E is already
    on the primal plane (just take layer ``ke``), while H is interpolated from the cell centres on
    both sides with ``kh/wh``. **Written identically for numpy and cupy arrays** (only slicing and
    multiply-add); the solver and the tests share this one implementation.
    """
    t1, t2 = g["tax"]
    pe = g["planes"][plane]
    t_ax = 0 if comp in ("Eu", "Hv") else 1
    kn, (w0, w1) = (pe["ke"], (1.0, 0.0)) if comp[0] == "E" else (pe["kh"], pe["wh"])
    sl = [None, None, None]
    sl[t1] = g["sl_cell"][0] if t_ax == 0 else g["sl_edge"][0]
    sl[t2] = g["sl_cell"][1] if t_ax == 1 else g["sl_edge"][1]
    sl[g["axis"]] = kn
    a = F[tuple(sl)]
    if w1 != 0.0:                     # E is a single layer, kn+1 is not read
        sl[g["axis"]] = kn + 1
        a = w0 * a + w1 * F[tuple(sl)]
    return _center_to_edge(a, g, comp, t_ax)


def integrate_colocated(eu, ev, hu, hv, g: dict) -> float:
    """Instantaneous flux scalar of the colocated plane fields: ``Σ dS·(eu·hv - ev·hu)``
    (dS already carries the sign)."""
    return float(((eu * hv - ev * hu) * g["dS"]).sum())


def flux_time_sample(fields_e, fields_h, g: dict) -> float:
    """Host reference: given the 3D field arrays ``(E_t1, E_t2)`` / ``(H̄_t1, H̄_t2)``, produce one
    flux scalar in the convention of :func:`flux_time_geometry` (integrate per plane, combine by
    weight). The solver's GPU path is term-by-term identical to this (the split into two launches
    only touches the time averaging)."""
    s = 0.0
    for p, pe in enumerate(g["planes"]):
        s += pe["w"] * integrate_colocated(
            colocated_component(fields_e[0], g, "Eu", p),
            colocated_component(fields_e[1], g, "Ev", p),
            colocated_component(fields_h[0], g, "Hu", p),
            colocated_component(fields_h[1], g, "Hv", p), g)
    return s
