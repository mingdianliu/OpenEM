# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Symmetry folding: fold a Scene that declares ``symmetry`` down to a half, quarter or eighth
domain.

Tidy3D's solver runs on the reduced domain too (29 of the 47 validation cases have a symmetry
plane), so keeping the full domain gives away a factor of 2-8 for nothing. Folding **changes no
kernel**: a symmetry plane becomes a new kind of low-end boundary (``SymmetryPEC`` /
``SymmetryPMC``, see grid.py) and the physics is carried by the existing index tables and masks.

-- Discrete derivation (folding x, symmetry plane at the low end; y/z are isomorphic) -----------

Yee positions (top of yee.cu): along x, **edge positions** are Ey, Ez, Hx; **centre positions**
are Ex, Hy, Hz. The symmetry plane sits on an edge (the low boundary of cell m; folding keeps
cells m..n-1, and folded index i = full-domain m+i).

Mirror eigenvalues (Tidy3D convention: ``symmetry=+1`` is the PMC type, tangential E even;
``-1`` is the PEC type, tangential E odd):

    σ(Ey)=σ(Ez)=s   σ(Ex)=-s   σ(Hy)=σ(Hz)=-s   σ(Hx)=s

Which updates read past the plane? Checking yee.cu one by one:

* **The whole H step is safe**: the forward difference only reads towards the high end, and the
  update of Hx (an edge position, on the plane) contains only y/z derivatives, no read along x.
* **Only Ey/Ez on the plane** in the E step: their backward difference reads Hz/Hy at
  centre(m-1), whose mirror is ``σ(H_tangential) · H(centre m)``. Both components have the same
  eigenvalue (-s for both), so the existing pair of tables ``prv[0]=0, mpv[0]=-s`` is enough
  (the symmetry branch in grid.py).
* With s=-1 the tangential E on the plane is identically zero, which is the existing ``pec[0]=0``
  mechanism.
* Hx on the plane is correct automatically: with s=-1 its curl source (Ey/Ez on the plane) is
  zero, so from a zero initial value it stays zero; with s=+1 it evolves as the symmetric
  solution.
* ``dl_dual[0]``: on a symmetric grid the mirrored cell spacing equals ``dl[0]``, the **same**
  value a non-periodic axis already uses.

Hence the folded half-domain evolution == the full-domain evolution restricted to the upper half
(**bitwise** in a symmetric scene, as long as the ε arrays are themselves bitwise symmetric; the
A/B test verifies that rather than assuming it).

-- Scope (fail closed) -------------------------------------------------------------------------

v1 supports: eps/sigma, PML (the PML face at the low end of the folded axis is dropped), absorber
slabs, dispersion CSR, ModeSource (normal x, folding y/z only), FluxMonitor / FieldMonitor
(including clipping across the plane), ModeMonitorSpec, PointDipole (keep the upper-half stencil
points, verify **point by point** that the lower half is the σ mirror, and fail closed if it is
not).
Any other capability that is non-empty raises ``NotImplementedError``: better not to fold than to
fold wrongly.

Monitors that straddle the plane only store the upper half; expansion (mirror-rebuilding the full
plane) happens in **post-processing**: ``expand_field``.

State of play (2026-09): in production, folding is only used for speed benchmarks (the benchmark
paths all do ``sc, _ = fold_scene(sc)`` and throw the meta away), and ``expand_field`` /
``expand_flux_time`` plus the three full-domain tables of ``FoldMeta`` are consumed only by
tests/test_fold.py. The interface is kept (tests guard it); take it from here when wiring it into
acc_solve.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from openem.grid import (
    MIRROR_SIGN, YEE_ON_EDGE, Axis, Grid, ravel_cells, transverse_axes, unravel_cells,
)
from openem.model import (
    Dispersion, PointDipole, Scene,
)

#: Both the mirror-eigenvalue table σ = sign·s and the Yee position table (edge positions
#: i <-> 2m-i, with i=m on the plane mapping to itself; centre positions i <-> 2m-1-i) live in
#: grid.py: ``MIRROR_SIGN[comp][ax]`` and ``YEE_ON_EDGE[comp][ax]``.


@dataclass(frozen=True)
class FoldAxis:
    """One folded axis. ``m`` is the edge/cell index of the symmetry plane (folding keeps
    m..n-1)."""

    axis: int
    s: int          #: +1 (PMC type) / -1 (PEC type)
    m: int          #: cell index of the symmetry plane in the full domain
    n_full: int     #: number of cells in the full domain
    wrap: bool = False   #: double mirror wall (periodic axis, s=+1): folded domain has a wrap cell

    def folded_shape(self, shape_full: tuple) -> tuple:
        """Full-domain shape -> folded shape: the folded axis becomes ``n_full - m`` (plus one
        wrap-around cell when wrap is set), the other axes stay put. This is the only place the
        folded axis length is defined."""
        out = list(shape_full)
        out[self.axis] = self.n_full - self.m + (1 if self.wrap else 0)
        return tuple(out)


@dataclass(frozen=True)
class FoldMeta:
    """Everything the expansion (the mirror rebuild in post-processing) needs."""

    folds: tuple[FoldAxis, ...]
    #: Per FieldMonitor: name -> (full-domain origin, full-domain box), used by the expansion
    field_monitors_full: dict[str, tuple[tuple[int, int, int], tuple[int, int, int]]]
    #: Per FieldTimeMonitor: name -> (full-domain origin, full-domain box), same convention
    field_time_full: dict[str, tuple[tuple[int, int, int], tuple[int, int, int]]] | None = None
    #: Per FluxTimeMonitor: name -> ((folded monitor name, coefficient), ...); the full-domain
    #: series is Σ coefficient·folded series (see expand_flux_time)
    flux_time_map: dict[str, tuple] | None = None


def _check_axis_symmetric(ax: Axis, axis: int) -> int:
    """Verify that the grid is symmetric about its centre; return the cell index m of the
    symmetry plane."""
    e = np.asarray(ax.edges, dtype=np.float64)
    n = ax.n
    if n % 2 != 0:
        raise NotImplementedError(
            f"axis {axis}: the cell count {n} is odd, so the symmetry plane misses the edges")
    c = 0.5 * (e[0] + e[-1])
    resid = np.abs((e + e[::-1]) - 2.0 * c)
    scale = max(abs(e[0]), abs(e[-1]))
    if resid.max() > 1e-9 * scale:
        raise NotImplementedError(
            f"axis {axis}: the grid is not symmetric about its centre (largest residual "
            f"{resid.max():.3e}), not folding")
    m = n // 2
    if not np.isclose(e[m], c, rtol=0, atol=1e-9 * scale):
        raise NotImplementedError(f"axis {axis}: the centre is not on edge[{m}], not folding")
    return m


def _slice_nd(a: np.ndarray | None, axis: int, m: int) -> np.ndarray | None:
    """Take ``[m:]`` along ``axis`` (the upper half that folding keeps); ``None`` comes back
    unchanged."""
    if a is None:
        return None
    sl = [slice(None)] * a.ndim
    sl[axis] = slice(m, None)
    return np.ascontiguousarray(a[tuple(sl)])


def _slice_wrap(a: np.ndarray | None, axis: int, m: int) -> np.ndarray | None:
    """Double mirror wall: the upper half [m:] plus the wrap-around cell (full-domain row 0)."""
    if a is None:
        return None
    hi = [slice(None)] * a.ndim
    hi[axis] = slice(m, None)
    z0 = [slice(None)] * a.ndim
    z0[axis] = slice(0, 1)
    return np.ascontiguousarray(
        np.concatenate([a[tuple(hi)], a[tuple(z0)]], axis=axis))


def _fold_cells(cell, shape_full: tuple,
                fa: FoldAxis) -> tuple[np.ndarray, np.ndarray]:
    """The shared first half of folding an entry table: which entries are kept (``keep``) and the
    flat index of the kept cells in the folded domain (``new_cell``, computed for every entry; the
    caller then selects with keep).

    wrap (double mirror wall): entries on full-domain row 0 are kept as well and map to the last
    row of the folded domain.
    """
    axis, m = fa.axis, fa.m
    shape_new = fa.folded_shape(shape_full)
    cell = np.asarray(cell, dtype=np.int64)
    ijk = list(unravel_cells(cell, shape_full))
    coord = ijk[axis]
    keep = coord >= m
    if fa.wrap:
        keep = keep | (coord == 0)
        ijk[axis] = np.where(coord >= m, coord - m, shape_new[axis] - 1)
    else:
        ijk[axis] = coord - m
    return keep, ravel_cells(*ijk, shape_new)


def _fold_csr(ofs, keep: np.ndarray, *arrs):
    """CSR compaction: cut the entries by ``keep``, recompute the offsets, and select the
    pole-level arrays by expanding the entry-level keep out to the poles. Returns
    ``(new_ofs (int64), *the cut arrays)``; the caller settles the dtype."""
    ofs = np.asarray(ofs, dtype=np.int64)
    per = np.diff(ofs)
    pk = np.repeat(keep, per)
    new_ofs = np.zeros(int(keep.sum()) + 1, dtype=ofs.dtype)
    np.cumsum(per[keep], out=new_ofs[1:])
    return (new_ofs, *[np.ascontiguousarray(a[pk]) for a in arrs])


def _fold_dispersion(d: Dispersion, shape_full: tuple, fa: FoldAxis) -> Dispersion:
    """Drop the entries in the lower half and renumber cell to the folded flat index.

    wrap (double mirror wall): entries on full-domain row 0 are kept as well and map to the last
    row of the folded domain.
    """
    keep, new_cell = _fold_cells(d.cell, shape_full, fa)
    new_ofs, am1, b = _fold_csr(d.pole_ofs, keep, d.am1, d.b)
    return Dispersion(
        comp=np.ascontiguousarray(d.comp[keep]),
        cell=np.ascontiguousarray(new_cell[keep].astype(d.cell.dtype)),
        pole_ofs=new_ofs.astype(d.pole_ofs.dtype),
        am1=am1,
        b=b,
        g=np.ascontiguousarray(d.g[keep]),
    )


def _fold_mix(mx, shape_full: tuple, fa: FoldAxis):
    """Cut the mix entries (the kernel is cell-local with no neighbour reads, so this is pure CSR
    surgery)."""
    keep, new_cell = _fold_cells(mx.cell, shape_full, fa)
    p_ofs, pa, pb = _fold_csr(mx.p_ofs, keep, mx.pa, mx.pb)
    q_ofs, qa, qb = _fold_csr(mx.q_ofs, keep, mx.qa, mx.qb)
    return replace(
        mx,
        comp=np.ascontiguousarray(mx.comp[keep]),
        cell=np.ascontiguousarray(new_cell[keep].astype(mx.cell.dtype)),
        p_ofs=p_ofs.astype(mx.p_ofs.dtype), pa=pa, pb=pb,
        q_ofs=q_ofs.astype(mx.p_ofs.dtype), qa=qa, qb=qb,
        beta=np.ascontiguousarray(mx.beta[keep]),
        eps_inf=np.ascontiguousarray(mx.eps_inf[keep]),
        zeta_inf=np.ascontiguousarray(mx.zeta_inf[keep]),
    )


def _fold_dipole(d: PointDipole, fa: FoldAxis) -> PointDipole | None:
    """Keep the stencil points in the upper half and verify that the discarded lower-half points
    are the σ mirror of the kept ones.

    The mirror eigenvalue of the current density J is the same as that of the E component of the
    same name (it is a polar vector); a magnetic current matches the H component. ``coef``
    contains 1/dV, and on a symmetric grid the mirrored point has the same dV, so coef can be
    compared directly. Anything that does not match (an orphan point, a mismatched coefficient)
    fails closed: better not to fold at all.
    """
    axis, m = fa.axis, fa.m
    comp = ("Hx", "Hy", "Hz")[d.component] if d.magnetic else ("Ex", "Ey", "Ez")[d.component]
    sig = MIRROR_SIGN[comp][axis] * fa.s
    edge = YEE_ON_EDGE[comp][axis]
    idx = np.asarray(d.indices)
    coord = idx[:, axis]
    keep = coord >= m
    kept = {tuple(int(v) for v in row): float(c)
            for row, c in zip(idx[keep], d.coef[keep])}
    for row, c in zip(idx[~keep], np.asarray(d.coef)[~keep]):
        g = int(row[axis])
        gm = (2 * m - g) if edge else (2 * m - 1 - g)
        key = list(int(v) for v in row)
        key[axis] = gm
        want = kept.get(tuple(key))
        if want is None or not np.isclose(float(c), sig * want, rtol=1e-12, atol=0):
            raise NotImplementedError(
                f"PointDipole lower-half point {tuple(row)} (coef={c:.6e}) is not the σ mirror "
                f"of the upper half (expected {sig:+d}x{want}): source and symmetry do not "
                "match, not folding")
    if edge and sig < 0:
        on = coord == m
        if np.any(np.asarray(d.coef)[on] != 0.0):
            raise NotImplementedError(
                "an odd-parity PointDipole component has a non-zero coefficient on the symmetry "
                "plane, not folding")
    if not keep.any():
        return None
    new_idx = idx[keep].copy()
    new_idx[:, axis] -= m
    return replace(d, indices=new_idx, coef=np.ascontiguousarray(d.coef[keep]))


def _fold_axis(ax: Axis, fa: FoldAxis) -> Axis:
    """The axis itself: take the upper half of edges and replace the low end with a symmetry
    boundary; with wrap, append one more wrap-around cell (width = full-domain cell 0) and make
    the high end a symmetry wall as well."""
    m = fa.m
    if fa.wrap:
        e_new = np.concatenate(
            [ax.edges[m:],
             [ax.edges[-1] + (ax.edges[1] - ax.edges[0])]])
    else:
        e_new = ax.edges[m:]
    _sym = "SymmetryPEC" if fa.s < 0 else "SymmetryPMC"
    return Axis(
        edges=np.ascontiguousarray(e_new),
        boundary_lo=_sym,
        boundary_hi=_sym if fa.wrap else ax.boundary_hi,
    )


def _fold_pml(pml: dict, fa: FoldAxis) -> dict:
    """PML: the low-end face of the folded axis is dropped (that end is a symmetry plane now).
    The coefficients are stored per layer and placed by _pml_axis_arrays from the axis length, so
    the other faces are kept as they are."""
    axis = fa.axis
    pml = {k: c for k, c in pml.items() if k != (axis, "lo")}
    if (axis, "hi") in pml and pml[(axis, "hi")].num_layers >= fa.n_full - fa.m:
        raise NotImplementedError(
            f"the high-end PML on axis {axis} has at least as many layers as the folded axis "
            "is long")
    return pml


def _fold_absorbers(absorbers: list, fa: FoldAxis) -> list:
    """Absorber slabs: on the folded axis the low-end slab lies entirely in the discarded lower
    half; the others are shifted."""
    new_abs = []
    for sl in absorbers:
        if sl.axis == fa.axis:
            if sl.g0 == 0:
                continue        # the low-end slab is entirely in the discarded lower half
            new_abs.append(replace(sl, g0=sl.g0 - fa.m))
        else:
            new_abs.append(sl)
    return new_abs


def _fold_box(fm, fa: FoldAxis, what: str):
    """Box monitors (FieldMonitor and FieldTimeMonitor share this code): clip the box to the
    upper half and **grow it until the expansion closes**.

    The lower-half cells of a box that straddles the plane are rebuilt by mirroring, and the
    mirror index reaches as far as 2m-lo (edge positions), so the folded box has to store up to
    there or expand_field cannot read it. The extra cells stored are a harmless superset.
    """
    axis, m, n_full, wrap = fa.axis, fa.m, fa.n_full, fa.wrap
    o = list(fm.origin)
    b = list(fm.box)
    lo, hi = o[axis], o[axis] + b[axis]      # [lo, hi)
    if hi <= m:
        raise NotImplementedError(
            f"{what} {fm.name} is entirely in the lower half, not supported in v1")
    if lo < m:
        hi = max(hi, min(2 * m - lo + 1, n_full))   # grow until the expansion closes
    new_lo = max(lo - m, 0)
    o[axis] = new_lo
    b[axis] = (hi - m) - new_lo
    if wrap and lo == 0:
        b[axis] += 1        # the mirror of g=0 is the wrap cell, the box must cover it
    return replace(fm, origin=tuple(o), box=tuple(b))


def _fold_flux_time(flux_parts: dict, ax: Axis, fa: FoldAxis) -> dict:
    """FluxTimeMonitor: each monitor becomes a linear combination of several (monitor,
    coefficient) pairs.

    Normal == the folded axis: the normal component of S is odd under that mirror, so a
    lower-half plane is replaced by its mirror plane and the coefficient changes sign.
    Tangential: S is even under a tangential mirror, so clip to the upper half; when both sides
    cover the whole upper axis this folds into a single monitor times 2 (the dual cell is cut in
    half at the plane, and the factor 2 restores the whole cell). The general straddling case
    splits into two pieces.
    """
    axis, m = fa.axis, fa.m
    x_m = float(ax.edges[m])
    top = float(ax.edges[-1])
    out = {}
    for _nm, parts in flux_parts.items():
        new_parts = []
        for fm, cf in parts:
            if fm.axis == axis:
                km, fr = fm.plane_index, float(fm.frac)
                below = (km < m) if fr == 0.0 else (km + 1 <= m)
                if not below:
                    if km == m and fr == 0.0:
                        raise NotImplementedError(
                            f"the plane of FluxTimeMonitor {fm.name} sits exactly on the "
                            "symmetry plane (the normal flux is identically mirror-odd), "
                            "not folding")
                    fm = replace(fm, plane_index=km - m)
                else:
                    if fr == 0.0:
                        km2 = 2 * m - km
                    else:
                        km2, fr = 2 * m - km - 1, 1.0 - fr
                    fm = replace(fm, plane_index=km2 - m, frac=fr)
                    cf = -cf
                new_parts.append((fm, cf))
                continue
            p_t = transverse_axes(fm.axis).index(axis)
            lo, hi = fm.t_bounds[p_t]

            def _tb(fm, lo, hi, _p=p_t):
                tb = list(fm.t_bounds)
                tb[_p] = (float(lo), float(hi))
                return replace(fm, t_bounds=tuple(tb))

            if hi <= x_m:                    # entirely in the lower half: mirror it (even)
                new_parts.append((_tb(fm, 2 * x_m - hi,
                                      2 * x_m - lo), cf))
            elif lo >= x_m:                  # entirely in the upper half
                new_parts.append((fm, cf))
            elif (2 * x_m - lo == hi
                  or (hi >= top and 2 * x_m - lo >= top)):
                # lower-half mirror coincides with the upper-half sample set: one monitor x2
                new_parts.append((_tb(fm, x_m, hi), 2.0 * cf))
            else:                            # general straddling case: split into two
                new_parts.append((_tb(fm, x_m, hi), cf))
                new_parts.append((_tb(fm, x_m, 2 * x_m - lo), cf))
        out[_nm] = new_parts
    return out


def _fold_mode_sources(mode_sources: list, fa: FoldAxis) -> list:
    """Mode sources: folding the normal axis shifts the plane; folding a tangential axis clips
    the profile (the profile is (t1, t2) in cyclic order t1=(axis+1)%3, t2=(axis+2)%3, see
    model.ModeSource)."""
    axis, m = fa.axis, fa.m
    new_ms = []
    for msrc in mode_sources:
        if axis == msrc.axis:
            if msrc.plane_index < m:
                raise NotImplementedError(
                    "the mode source plane is in the lower half, not supported in v1")
            new_ms.append(replace(msrc, plane_index=msrc.plane_index - m))
        else:
            prof_ax = 0 if axis == (msrc.axis + 1) % 3 else 1
            # The broadband correction profiles carry an extra leading node axis,
            # (K, n_t1, n_t2), so the spatial axis is one higher. Forgetting to clip them would
            # silently drop the correction terms, and the injection would lose the profile
            # corrections at both edges of the band.
            bb = {f: _slice_nd(getattr(msrc, f), prof_ax + 1, m)
                  for f in ("bb_ey", "bb_ez", "bb_hy", "bb_hz")
                  if getattr(msrc, f, None) is not None}
            new_ms.append(replace(
                msrc,
                ey_inc=_slice_nd(msrc.ey_inc, prof_ax, m),
                ez_inc=_slice_nd(msrc.ez_inc, prof_ax, m),
                hy_inc=_slice_nd(msrc.hy_inc, prof_ax, m),
                hz_inc=_slice_nd(msrc.hz_inc, prof_ax, m),
                **bb,
            ))
    return new_ms


def _check_plane_waves(sources: list, fa: FoldAxis) -> None:
    """Plane-wave sources: no spatial profile, so only parity compatibility matters.

    E_pol is constant along the folded axis (an even function), so the mirror eigenvalue of that
    component must be +1:
      fold axis == pol axis (E normal, σ=-s)     => requires s=-1;
      fold axis ⊥ pol axis  (E tangential, σ=+s) => requires s=+1.
    The paired H component is then automatically compatible (it lives on the third axis, the one
    that is neither of those two).
    If the fold axis is the propagation axis the wave crosses the symmetry plane, which is not
    symmetric, so it is refused outright; the same goes for oblique incidence.
    Compatible => the source itself needs no change (the time table is 1D and the injection plane
    automatically covers the folded transverse extent).
    """
    axis, s = fa.axis, fa.s
    for src in sources:
        if src.axis == axis:
            raise NotImplementedError(
                f"the plane wave propagates along the folded axis {axis}, which is incompatible "
                "with the symmetry plane")
        if getattr(src, "angle_theta", 0.0):
            raise NotImplementedError(
                "an obliquely incident plane wave is incompatible with the symmetry plane")
        need = -1 if src.pol_axis == axis else +1
        if s != need:
            raise NotImplementedError(
                f"plane-wave polarization axis {src.pol_axis} is incompatible with the symmetry "
                f"of axis {axis}, s={s:+d} (needs {need:+d})")


def _fold_flux_monitors(flux_mons: list, fa: FoldAxis) -> list:
    """FluxMonitor: folding the normal axis shifts the plane; folding a tangential axis clips the
    transverse range (transverse None = the whole plane, which after folding is automatically the
    half plane; expansion happens in post-processing)."""
    axis, m = fa.axis, fa.m
    new_fx = []
    for fm in flux_mons:
        if fm.axis == axis:
            if fm.plane_index < m:
                raise NotImplementedError(
                    f"the plane of FluxMonitor {fm.name} is in the lower half, "
                    "not supported in v1")
            fm = replace(fm, plane_index=fm.plane_index - m)
        elif fm.transverse is not None:
            # transverse order = the two transverse axes in ascending order; find where the
            # folded axis sits in it
            ti = transverse_axes(fm.axis).index(axis)
            rng = list(fm.transverse)
            lo, hi = rng[ti]
            if hi < m:
                raise NotImplementedError(
                    f"FluxMonitor {fm.name} is entirely in the lower half, "
                    "not supported in v1")
            rng[ti] = (max(lo - m, 0), hi - m)
            fm = replace(fm, transverse=tuple(rng))
        new_fx.append(fm)
    return new_fx


def _refuse_unsupported_fold(sc: Scene) -> None:
    """Anything folding v1 does not cover fails closed (the scope is in the module docstring)."""
    if not any(sc.symmetry):
        raise ValueError("the scene has no symmetry plane (symmetry is all zeros), nothing to fold")
    for nm in ("modulation", "tensor"):
        if getattr(sc, nm) is not None:
            raise NotImplementedError(f"folding v1 does not support a non-empty {nm}")
    # projection_monitors are pure metadata (the surfaces are already expanded into
    # FieldMonitors and the projection happens in post-processing), so they ride along with the
    # clipping and expansion of field_monitors and need no separate handling.
    for nm in ("tfsf_sources", "permittivity_monitors",
               "box_flux_monitors",
               "pec_ex", "pec_ey", "pec_ez"):
        if getattr(sc, nm):
            raise NotImplementedError(f"folding v1 does not support a non-empty {nm}")
    if any(sc.bloch_k):
        raise NotImplementedError("folding v1 does not support Bloch")


def _check_double_mirror(sc: Scene, absorbers: list, mode_sources: list,
                         axis: int) -> None:
    """Double mirror wall: periodic plus a centre mirror means the wrap edge is a second mirror
    wall (period L with symmetry at 0 implies symmetry at L/2 too). The old conclusion "cannot be
    folded" only holds for a fold that keeps the periodic wrap-around, where the high end wraps
    round and reads the non-zero rows of the other half (measured on an MIM structure, the T
    spectrum shifted by 0.47). This path does not wrap: it puts a second symmetry wall at the high
    end and keeps one extra wrap-around column (= full-domain row 0) to carry the tangential E row
    on that wall.

    Three things are not done on this path yet, and they fail closed. Mirror pairs of dipoles need
    no extra handling: a mirror pair about the centre, (y, n-1-y), is exactly a mirror pair about
    the wrap wall as well (the algebraic identity 2m=n), which the centre-mirror check in
    :func:`_fold_dipole` already covers.
    """
    for sl in absorbers:
        if sl.axis == axis:
            raise NotImplementedError(
                "there is an absorber on the double-mirror-wall axis, not supported")
    if mode_sources:
        raise NotImplementedError("double mirror wall v1 does not support mode sources")
    if sc.flux_time_monitors:
        raise NotImplementedError("double mirror wall v1 does not support flux_time")


def _flux_time_names(flux_parts: dict) -> tuple[list, dict]:
    """Flatten "each full-domain monitor = several (folded monitor, coefficient) pairs" into the
    list of monitors actually to be run, plus the ``name -> ((folded name, coefficient), ...)``
    map that reconstructs them (see expand_flux_time)."""
    ftime_mons = []
    flux_time_map = {}
    for _nm, parts in flux_parts.items():
        entry = []
        for i, (fm, cf) in enumerate(parts):
            nm2 = _nm if i == 0 else f"{_nm}__fold{i}"
            ftime_mons.append(replace(fm, name=nm2))
            entry.append((nm2, float(cf)))
        flux_time_map[_nm] = tuple(entry)
    return ftime_mons, flux_time_map


def fold_scene(sc: Scene) -> tuple[Scene, FoldMeta]:
    """Fold according to ``sc.symmetry``. Unsupported combinations complain loudly (fail
    closed)."""
    _refuse_unsupported_fold(sc)

    folds: list[FoldAxis] = []
    axes = list(sc.grid.axes)
    dipoles = list(sc.dipoles)
    eps = [sc.eps_ex, sc.eps_ey, sc.eps_ez]
    sig = [sc.sigma_ex, sc.sigma_ey, sc.sigma_ez]
    disp = sc.dispersion
    mix = sc.dispersion_mix
    absorbers = list(sc.absorbers)
    mode_sources = list(sc.mode_sources)
    flux_mons = list(sc.flux_monitors)
    field_mons = list(sc.field_monitors)
    pml = dict(sc.pml)
    shape = list(sc.shape)
    fm_full = {fm.name: (tuple(fm.origin), tuple(fm.box)) for fm in field_mons}
    ft_mons = list(sc.field_time_monitors)
    ft_full = {fm.name: (tuple(fm.origin), tuple(fm.box)) for fm in ft_mons}
    # Each flux_time monitor folds into a linear combination of several (monitor, coefficient)
    # pairs
    flux_parts = {fm.name: [(fm, 1.0)] for fm in sc.flux_time_monitors}

    for axis in range(3):
        s = int(sc.symmetry[axis])
        if s == 0:
            continue
        ax = axes[axis]
        m = _check_axis_symmetric(ax, axis)
        fa = FoldAxis(axis=axis, s=s, m=m, n_full=ax.n,
                      wrap=bool(ax.is_periodic and s == +1))
        if fa.wrap:
            _check_double_mirror(sc, absorbers, mode_sources, axis)
        folds.append(fa)

        axes[axis] = _fold_axis(ax, fa)

        # ---- dense arrays ----
        _cut = (_slice_wrap if fa.wrap else _slice_nd)
        eps = [_cut(a, axis, m) for a in eps]
        sig = [_cut(a, axis, m) for a in sig]

        pml = _fold_pml(pml, fa)
        absorbers = _fold_absorbers(absorbers, fa)
        # FieldTimeMonitor: the same clipping and expansion closure as FieldMonitor
        ft_mons = [_fold_box(fm, fa, "FieldTimeMonitor") for fm in ft_mons]
        if flux_parts:
            flux_parts = _fold_flux_time(flux_parts, ax, fa)

        # ---- dispersion CSR ----
        shape_full = tuple(shape)
        shape = list(fa.folded_shape(shape_full))
        if disp is not None:
            disp = _fold_dispersion(disp, shape_full, fa)
        if mix is not None:
            mix = _fold_mix(mix, shape_full, fa)

        mode_sources = _fold_mode_sources(mode_sources, fa)
        _check_plane_waves(sc.sources, fa)

        # ---- dipoles: keep upper-half stencil points, verify the lower half is the σ mirror ----
        dipoles = [_fold_dipole(d, fa) for d in dipoles]
        dipoles = [d for d in dipoles if d is not None]

        flux_mons = _fold_flux_monitors(flux_mons, fa)
        field_mons = [_fold_box(fm, fa, "FieldMonitor") for fm in field_mons]

    ftime_mons, flux_time_map = _flux_time_names(flux_parts)

    folded = replace(
        sc,
        grid=Grid(*axes),
        eps_ex=eps[0], eps_ey=eps[1], eps_ez=eps[2],
        sigma_ex=sig[0], sigma_ey=sig[1], sigma_ez=sig[2],
        pml=pml,
        dispersion=disp,
        dispersion_mix=mix,
        absorbers=absorbers,
        dipoles=dipoles,
        field_time_monitors=ft_mons,
        mode_sources=mode_sources,
        flux_monitors=flux_mons,
        field_monitors=field_mons,
        flux_time_monitors=ftime_mons,
        symmetry=(0, 0, 0),
    )
    return folded, FoldMeta(folds=tuple(folds), field_monitors_full=fm_full,
                            field_time_full=ft_full,
                            flux_time_map=flux_time_map or None)



# ----------------------------- post-processing expansion -----------------------------

def expand_field(ph: np.ndarray, comp: str, fold: FoldAxis,
                 lo_folded: int, lo_full: int, n_full: int) -> np.ndarray:
    """Mirror one component's folded box back out to the full-domain box along the folded axis.

    Args:
        ph: ``(..., n_folded)``, with the folded axis already moved to the last dimension. The
            folded box starts at folded index ``lo_folded`` (0 = the cell holding the symmetry
            plane).
        comp: ``"Ex"`` ... ``"Hz"``.
        fold: information about the folded axis.
        lo_full: starting cell index of the full-domain box along that axis.
        n_full: length of the full-domain box along that axis.

    Mirroring: edge positions map full-domain cell g <-> 2m-g, centre positions g <-> 2m-1-g; the
    value is multiplied by σ = MIRROR_SIGN·s.
    """
    m = fold.m
    sig = MIRROR_SIGN[comp][fold.axis] * fold.s
    edge = YEE_ON_EDGE[comp][fold.axis]
    n_axis = ph.shape[-1]
    out_shape = ph.shape[:-1] + (n_full,)
    out = np.zeros(out_shape, dtype=ph.dtype)
    for t in range(n_full):
        g = lo_full + t                        # full-domain cell index
        if g >= m:                             # upper half: read straight off
            src = (g - m) - lo_folded
            if 0 <= src < n_axis:
                out[..., t] = ph[..., src]
            # past the folded box (fold_scene only grows along the folded axis, never shrinks,
            # so this branch should not occur)
        else:                                  # lower half: mirror
            gm = (2 * m - g) if edge else (2 * m - 1 - g)
            src = (gm - m) - lo_folded
            if 0 <= src < n_axis:
                out[..., t] = sig * ph[..., src]
            # else: the mirror falls off the axis (an edge position g=0 mirrors to edge n, the
            # wall ghost). In the full domain that cell is the wall behind the PML/absorber, no
            # criterion ever counts it, so fill 0
    return out


def expand_flux_time(meta: FoldMeta,
                     series: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Recombine the flux_time series of a folded run linearly back to the full-domain
    convention."""
    out = {}
    for name, parts in (meta.flux_time_map or {}).items():
        acc = None
        for nm, cf in parts:
            v = cf * series[nm]
            acc = v if acc is None else acc + v
        out[name] = acc
    return out
