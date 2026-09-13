# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Oblique-incidence TFSF: the 1D auxiliary line runs along k̂ and the box shell interpolates
its values at the projected position.

The normal-incidence path (the 16 terms in ``sources_setup._tfsf_box``) is **left untouched**
and oblique incidence comes here; both paths feed the same kernel (``kernels/tfsf.cu``) and
differ only in whether the table column position in the descriptor is an integer.

Geometry
--------
The 1D line runs along **+k̂** (its propagation direction is always +1, with the physical
direction folded into k̂), and its cell size dl₁ comes from
:func:`openem.tfsf1d.matched_dl`, which makes the numerical phase velocity along k̂ agree with
the 3D grid. The projection ``P = r·k̂`` of a 3D point r maps onto a table column:

    E column j = (P - proj0)/dl₁          (E sits on 1D edge col0+j)
    H column j = (P - proj0)/dl₁ + ½      (H sits in 1D cell col0-1+j, at position
                                           (col0+j-½)dl₁)

Correction terms
----------------
Writing (u, v, a) = ((axis+1)%3, (axis+2)%3, axis) in cyclic right-handed order,

    (∇×F)_u = +∂F_a/∂v - ∂F_v/∂a        (the others follow cyclically)

For the term ``σ·∂H_r/∂n`` of (∇×H)_m: E_m is corrected by ``-σ·ĥ_r/dld_n`` on the lo face of
n and by ``+σ·ĥ_r/dld_n`` on the hi face; for ``σ·∂E_r/∂n`` in (∇×E)_m: H_m gets
``+ch·σ·ê_r/dl_n`` on the lo cell and ``-ch·σ·ê_r/dl_n`` on the hi cell. At normal incidence
ê_a = ĥ_a = 0, so 8 of these vanish by themselves and the remaining 16 match the hand-written
table of ``_tfsf_box`` term for term (the derivation was cross-checked by working it back out
with sympy). Oblique incidence switches those 8 on, 24 terms in all.

The TF/SF interface is at ``edge box_lo[n]``: edges and cells with n >= lo are the total-field
region. So **the E correction takes the H outside the interface** (cell lo-1 / cell hi) and
**the H correction takes the E inside it** (edge lo / edge hi).
"""

from __future__ import annotations

import numpy as np

from openem.grid import cyclic_axes, transverse_axes, yee_on_edge

#: The two terms of (∇×F)_m: (r, n, σ), meaning σ·∂F_r/∂n. (u, v, a) are written as 0/1/2.
_CURL = {0: ((2, 1, +1), (1, 2, -1)),
         1: ((0, 2, +1), (2, 0, -1)),
         2: ((1, 0, +1), (0, 1, -1))}


def coord(grid, comp: int, magnetic: bool, idx) -> np.ndarray:
    """Physical coordinates (x, y, z) of component ``comp`` at 3D index ``idx``, in meters."""
    on_edge = yee_on_edge(comp, magnetic)
    out = np.empty(3, dtype=np.float64)
    for w in range(3):
        ax = grid.axes[w]
        out[w] = ax.edges[idx[w]] if on_edge[w] else ax.centers[idx[w]]
    return out


def basis(axis: int, k_hat, e_hat):
    """Components of k̂, ê and ĥ in the (u, v, a) basis. ĥ = k̂ x ê (the 1D line always
    propagates in the + direction)."""
    u, v, a = _uva_axes(axis)
    k = np.asarray(k_hat, dtype=np.float64)
    e = np.asarray(e_hat, dtype=np.float64)
    h = np.cross(k, e)

    def take(w):
        return float(w[u]), float(w[v]), float(w[a])

    return take(k), take(e), take(h)


def _uva_axes(axis: int) -> tuple[int, int, int]:
    """(u, v, a) in cyclic right-handed order: the two transverse axes come from
    grid.cyclic_axes and the third is the injection axis itself."""
    return (*cyclic_axes(axis), axis)


def face_terms(axis: int, e_uva, h_uva, box_lo, box_hi, grid, ch: float,
               open_axes: tuple[int, ...] = ()):
    """List the ingredients of the 24 kinds of face correction, without the table column
    position.

    Yields:
        dict: ``kind`` ("E"/"H"), ``comp`` (axis of the corrected component), ``inc`` (axis of
        the incident component), ``fixed_ax``, ``fixed_idx``, ``inc_idx_n`` (index of the
        incident point on axis n) and ``coef``. Terms whose coefficient is 0 are not produced
        (the 8 that vanish at normal incidence).
    """
    uva = _uva_axes(axis)
    lo = {q: box_lo[uva[q]] for q in range(3)}
    hi = {q: box_hi[uva[q]] for q in range(3)}
    for m in range(3):                       # corrected component, index within (u,v,a)
        for (r, n, sg) in _CURL[m]:
            n_ax = uva[n]
            if n_ax in open_axes:
                continue                     # box fills the domain on this axis, no face

            dld = grid.axes[n_ax].dl_dual
            dl = grid.axes[n_ax].dl
            # ---- E step: correct E_m on the total-field side, incident H_r from outside ----
            for side, fx, cf in (("lo", lo[n], -sg * h_uva[r] / dld[lo[n]]),
                                 ("hi", hi[n], +sg * h_uva[r] / dld[hi[n]])):
                if cf != 0.0:
                    yield {"kind": "E", "comp": uva[m], "inc": uva[r],
                           "fixed_ax": n_ax, "fixed_idx": fx,
                           "inc_idx_n": fx - 1 if side == "lo" else fx,
                           "coef": cf, "m": m, "n": n, "side": side}
            # ---- H step: correct H_m on the scattered-field side, incident E_r from inside ----
            for side, fx, cf in (("lo", lo[n] - 1, +ch * sg * e_uva[r] / dl[lo[n] - 1]),
                                 ("hi", hi[n], -ch * sg * e_uva[r] / dl[hi[n]])):
                if cf != 0.0:
                    yield {"kind": "H", "comp": uva[m], "inc": uva[r],
                           "fixed_ax": n_ax, "fixed_idx": fx,
                           "inc_idx_n": fx + 1 if side == "lo" else fx,
                           "coef": cf, "m": m, "n": n, "side": side}


def projection_range(grid, k_hat, box_lo, box_hi) -> tuple[float, float]:
    """Projection interval [P_min, P_max] of the box (grown by one cell) onto k̂, in meters.

    Every shell sampling position falls inside this interval, so taking the extremes of the end
    points axis by axis is enough (the projection is monotonic in each coordinate) and is
    sturdier than enumerating the sampling points.
    """
    pmin = pmax = 0.0
    for w in range(3):
        e = grid.axes[w].edges
        # on an open axis the box touches both ends of the domain, so growing it by one cell
        # would run past them: clamp back into the valid range
        a = float(k_hat[w]) * float(e[max(box_lo[w] - 1, 0)])
        b = float(k_hat[w]) * float(e[min(box_hi[w] + 1, e.size - 1)])
        pmin += min(a, b)
        pmax += max(a, b)
    return pmin, pmax


def descriptors(t, grid, strides, ch: float):
    """Descriptors of the 24 kinds of face correction (without ``launch``, which the caller
    fills in).

    The table column positions ``c0/cp/cq`` are real numbers. **Substituting normal incidence
    reproduces the integers of the hand-written table exactly** (checked by hand: the E_u@a
    face gives lo->0 and hi->nk_e; the H_v@a cell gives lo->0 and hi->nk_e-1; the E_a@u face
    gives c0=1 and cq=1), so these geometric formulas and ``_tfsf_box`` are two ways of writing
    the same thing; the other one keeps its hand-written constants only to stay bitwise
    unchanged.

    Returns:
        ``(e_corr, h_corr)``, whose elements are the same kind of dict ``_tfsf_box`` uses.
    """
    k = np.asarray(t.k_hat, dtype=np.float64)
    dl1 = float(t.dl1_step)
    _, e_uva, h_uva = basis(t.axis, t.k_hat, t.e_hat)
    e_corr: list[dict] = []
    h_corr: list[dict] = []
    for term in face_terms(t.axis, e_uva, h_uva, t.box_lo, t.box_hi, grid, ch,
                           getattr(t, "open_axes", ())):
        comp, fixed_ax = term["comp"], term["fixed_ax"]
        p_ax, q_ax = transverse_axes(fixed_ax)
        mag_corr = term["kind"] == "H"
        on_edge = yee_on_edge(comp, mag_corr)
        n_p = t.box_hi[p_ax] - t.box_lo[p_ax] + (1 if on_edge[p_ax] else 0)
        n_q = t.box_hi[q_ax] - t.box_lo[q_ax] + (1 if on_edge[q_ax] else 0)
        idx = [0, 0, 0]
        idx[fixed_ax] = term["fixed_idx"]
        idx[p_ax], idx[q_ax] = t.box_lo[p_ax], t.box_lo[q_ax]
        base = sum(idx[w] * strides[w] for w in range(3))
        inc_idx = list(idx)
        inc_idx[fixed_ax] = term["inc_idx_n"]
        # incident component: the E correction takes H (magnetic), the H correction takes E
        p0 = float(k @ coord(grid, term["inc"], not mag_corr, inc_idx))
        half = 0.5 if term["kind"] == "E" else 0.0      # the H table sits on half cells
        c0 = (p0 - float(t.proj0)) / dl1 + half
        cp = float(k[p_ax]) * float(grid.axes[p_ax].dl[t.box_lo[p_ax]]) / dl1
        cq = float(k[q_ax]) * float(grid.axes[q_ax].dl[t.box_lo[q_ax]]) / dl1
        (e_corr if term["kind"] == "E" else h_corr).append({
            "comp": comp, "coef": np.float32(term["coef"]),
            "c0": np.float64(c0), "cp": np.float64(cp), "cq": np.float64(cq),
            "base": np.int64(base),
            "sp": np.int64(strides[p_ax]), "sq": np.int64(strides[q_ax]),
            "np": np.int32(n_p), "nq": np.int32(n_q),
        })
    return e_corr, h_corr


def check_in_table(corr, row_n: int, what: str) -> None:
    """The table column positions of every correction must lie entirely inside the table: the
    out-of-range clamp in the kernel is a safety net, not permission to truncate silently. The
    four corners are enough (the position is linear in (p, q))."""
    for d in corr:
        for p in (0, int(d["np"]) - 1):
            for q in (0, int(d["nq"]) - 1):
                pos = float(d["c0"]) + float(d["cp"]) * p + float(d["cq"]) * q
                if not (0.0 <= pos <= row_n - 1):
                    raise ValueError(
                        f"{what} correction table column position {pos:.3f} is outside "
                        f"[0, {row_n - 1}]: the 1D table does not cover enough of the "
                        f"projection range (comp={int(d['comp'])})")
