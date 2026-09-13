# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Source tables: plane wave, TFSF box, mode source, dipole, time-varying medium. Turns the
sources in a ``Scene`` into the indices, coefficients and incident-field tables the kernel wants.
(Moved verbatim out of setup_tables.py, 2026-09.)
"""

from __future__ import annotations

from typing import NamedTuple

import cupy as cp
import numpy as np

from openem import cpml, tfsf_oblique, tfsf1d
from openem.device import grid_1d, grid_2d
from openem.grid import cyclic_axes, ravel_cells
from openem.model import Scene
from openem.pitch import Dims, _f32, _f64, _pitch, _repitch_cells



#: The sign table for the two TF/SF correction terms is tfsf1d.SIGNS (the 1D line and the 3D
#: injection must use the same table).

#: Signs of the four TF/SF correction terms of a mode source, in the order
#: ``(Hy<-Ez, Hz<-Ey, Ey<-Hz, Ez<-Hy)``.
#: Derived from the Yee curl equations (``μ∂Hy/∂t = ∂Ez/∂x - ...``, ``ε∂Ey/∂t = ... - ∂Hz/∂x``
#: and so on), but **the measurement decides**: the criterion is unidirectionality in a uniform
#: waveguide, and a wrong sign does not diverge, it only lifts the backward amplitude from 1e-8
#: to order 1e-1.
#: The ``-1`` entry is **not** the ``+1`` entry negated as a whole. Term by term: reversing the
#: direction swaps the inner and outer side of the TF/SF surface geometrically, so all four terms
#: flip once; but **the tangential H of a backward mode is itself negated** (E is unchanged), so
#: the two terms built from H_inc (``Ey<-Hz``, ``Ez<-Hy``) flip a second time, i.e. no net flip.
#: It used to be written as a whole-tuple negation, which injects -(forward mode), still a
#: **forward** wave. The criterion is unambiguous: for a ``direction="-"`` source in a straight
#: waveguide, ``|a_-|²`` should be ~1 upstream and ~0 downstream.
#: All four combinations were measured (sweeping at the same time whether the half-cell phase in
#: modes.py carries direction):
#:
#:   E kept    + phase without direction -> upstream 0.998245, downstream **1.08e-07**  <- chosen
#:   E kept    + phase with direction    -> upstream 1.018898, downstream 2.07e-02
#:   E flipped + phase with direction    -> upstream 0.000000, downstream 9.98e-01 (fires backwards)
#:   E flipped + phase without direction -> upstream 0.020653, downstream 1.02e+00
#:
#: 1.08e-07 is the same order as the unidirectionality on the forward side.
#: Not one mode source in the validation set uses ``direction="-"``, so this never surfaced;
#: and tidy3d's **mode adjoint source is precisely a backward one**, which is why inverse design
#: came back with the wrong gradient.
_MODE_SIGNS: dict[int, tuple[float, float, float, float]] = {
    +1: (-1.0, +1.0, +1.0, -1.0),
    -1: (+1.0, -1.0, +1.0, -1.0),
}


def _plane_source(sc: Scene) -> dict | None:
    """The handful of numbers a plane-wave source needs in the kernel: indices, coefficients and
    the two incident-field tables.

    ``None`` means the scene has no plane-wave source.
    """
    az = sc.grid.axes[2]
    # ---- plane-wave source ----
    plane_src = None
    if sc.sources:
        # Several plane-wave sources (TunableChiralMetasurface / Metalens: the Ex and Ey
        # polarizations are 90 degrees apart in phase and combine into circular polarization).
        # They must share the injection plane and the direction; group them by polarization and
        # add the incident tables inside a group (injection is linear and the phase is already in
        # each table). The first group is the main source, the rest go into "extra", and the
        # solver runs the injection kernel once per group.
        for s_ in sc.sources:
            if s_.axis != 2 or s_.pol_axis not in (0, 1):
                raise NotImplementedError(
                    "plane waves currently only support propagation along z with Ex/Ey "
                    "polarization")
        keys = {(int(s_.plane_index), int(s_.direction)) for s_ in sc.sources}
        if len(keys) != 1:
            raise NotImplementedError(
                f"multiple plane-wave sources must share one injection plane and one direction, "
                f"got {sorted(keys)}")
        groups: dict[int, list] = {}
        for s_ in sc.sources:
            groups.setdefault(int(s_.pol_axis), []).append(s_)
        built = [_plane_source_one(sc, az, grp) for _, grp in sorted(groups.items())]
        plane_src = built[0]
        if len(built) > 1:
            plane_src["extra"] = built[1:]
    return plane_src


def _plane_source_one(sc: Scene, az, grp: list) -> dict:
    """A group of plane-wave sources of one polarization -> kernel parameter dict (incident tables
    summed over the group)."""
    src = grp[0]
    ks = src.plane_index
    # The two TF/SF correction terms are mirror images of each other in the direction: the H
    # correction plane moves from ks-1 to ks and each of the two coefficients picks up a sign.
    # The signs were not derived, they were measured, against "unidirectionality in a uniform
    # lossless medium should be ~1e-8". The +z pair is the (+1, +1) verified originally.
    kh = ks - 1 if src.direction > 0 else ks
    sh, se = tfsf1d.SIGNS[src.direction]
    # Ey polarization lands on (Ey, Hx) and reuses the same kernel: relative to (Ex, Hy), the H
    # update carries E with -∂/∂z and the E update carries H with +∂/∂z, one structural sign
    # each; and the incident field has Hx_inc = -Hy_inc. The two signs on the E side cancel, so
    # the net effect is only that the H-side coefficient is negated, with the incident tables
    # reused as they are. Criterion: Ey and Ex polarization agree on unidirectionality and on
    # downstream flux on the same uniform grid (tests/test_polarization.py).
    sp = -1.0 if src.pol_axis == 1 else 1.0
    return {
        "pol": int(src.pol_axis),
        "kh": np.int32(kh),
        "ks": np.int32(ks),
        "coef_h": np.float32(sp * sh * sc.dt / cpml.MU_0 * float(1.0 / az.dl[kh])),
        "coef_e": np.float32(se * float(1.0 / az.dl_dual[ks])),
        "ex_inc": np.ascontiguousarray(sum(x.ex_inc for x in grp), dtype=np.float64),
        "hy_inc": np.ascontiguousarray(sum(x.hy_inc for x in grp), dtype=np.float64),
    }


def _tfsf_refuse_bad_placement(sc: Scene, t, n_total: int) -> None:
    """Check the placement of the TFSF box and the step count; anything out of spec fails closed.

    The box faces (including the +-1 cells the corrections need) must not reach into a PML or an
    absorber, where the update is not the free-space one; and the injection table must hold enough
    samples for ``n_total`` steps. On an open axis the box fills the domain and has no face, so it
    is skipped.
    """
    _open = set(getattr(t, "open_axes", ()))
    for a2, (lo_e, hi_e) in enumerate(zip(t.box_lo, t.box_hi)):
        if a2 in _open:
            continue        # box fills the domain on this axis, no face, nothing to touch a PML
        n = sc.grid.axes[a2].n
        c = sc.pml.get((a2, "lo"))
        if c is not None and lo_e - 2 < c.num_layers:
            raise ValueError(f"TFSF box reaches into the PML at the low end of {'xyz'[a2]} "
                             f"(face {lo_e}, PML {c.num_layers} layers)")
        c = sc.pml.get((a2, "hi"))
        if c is not None and hi_e + 2 > n - c.num_layers:
            raise ValueError(f"TFSF box reaches into the PML at the high end of {'xyz'[a2]} "
                             f"(face {hi_e}, PML {c.num_layers} layers)")
    for sl in sc.absorbers:
        if sl.axis in _open:
            continue
        lo_e, hi_e = t.box_lo[sl.axis], t.box_hi[sl.axis]
        nlay = sl.decay_int.size
        if not (hi_e + 2 <= sl.g0 or lo_e - 2 >= sl.g0 + nlay):
            raise ValueError(f"TFSF box overlaps the absorber slab on the {'xyz'[sl.axis]} axis")
    if n_total > t.ex_inc.size - 1:
        raise ValueError(
            f"asked to run {n_total} steps, but the TFSF injection table only sampled "
            f"{t.ex_inc.size} points (nominal step count {sc.num_time_steps}).")


def _tfsf_face_corrections(t, grid, strides, nk_e: int, ch: float):
    """The 16 face-correction descriptors of normal-incidence TFSF: ``(e_corr, h_corr)``.

    Which correction lands on which face, and where the coefficients come from, is the table in
    the :func:`_tfsf_box` docstring; this transcribes it term by term. The half with a zero
    coefficient (the absent polarization component) is not generated.
    """
    a = t.axis
    u, v = cyclic_axes(a)
    bu, bv = float(t.pol_u), float(t.pol_v)
    alo, ahi = t.box_lo[a], t.box_hi[a]
    ulo, uhi = t.box_lo[u], t.box_hi[u]
    vlo, vhi = t.box_lo[v], t.box_hi[v]
    dl_a, dld_a = grid.axes[a].dl, grid.axes[a].dl_dual
    dl_u, dld_u = grid.axes[u].dl, grid.axes[u].dl_dual
    dl_v, dld_v = grid.axes[v].dl, grid.axes[v].dl_dual

    def _corr(out, comp, fixed_ax, fixed_idx, p_ax, np_, q_ax, nq_,
              coef, c0, cq_follow_a):
        """One face-correction descriptor. Nothing is generated when coef is 0 (the half where
        the polarization component is absent).

        The parameters are deliberately flat rather than packed: the 16 calls below are one
        aligned table with one meaning per column, and grouping parentheses would only blur it.
        """
        if coef == 0.0:
            return
        pt = [0, 0, 0]
        pt[fixed_ax] = fixed_idx
        pt[p_ax] = (ulo if p_ax == u else (vlo if p_ax == v else alo))
        pt[q_ax] = (ulo if q_ax == u else (vlo if q_ax == v else alo))
        base = pt[0] * strides[0] + pt[1] * strides[1] + pt[2] * strides[2]
        out.append({
            "comp": comp,
            "coef": np.float32(coef),
            # Table column positions are real numbers (oblique incidence interpolates, see
            # kernels/tfsf.cu). Normal incidence always passes integer values, so w in the
            # kernel is exactly 0 and it reads back the end column itself: bitwise unchanged.
            "c0": np.float64(c0),
            "cp": np.float64(0.0),
            "cq": np.float64(1.0 if cq_follow_a else 0.0),
            "base": np.int64(base),
            "sp": np.int64(strides[p_ax]), "sq": np.int64(strides[q_ax]),
            "np": np.int32(np_), "nq": np.int32(nq_),
            "launch": grid_2d(np_, nq_),
        })

    nu_c, nu_e = uhi - ulo, uhi - ulo + 1
    nv_c, nv_e = vhi - vlo, vhi - vlo + 1
    na_c, na_e = ahi - alo, ahi - alo + 1

    e_corr: list[dict] = []      # E step, row = hy_tab[n]
    h_corr: list[dict] = []      # H step, row = ex_tab[n]
    # E_u on the a faces (E_u: u cells x v edges)
    _corr(e_corr, u, a, alo, u, nu_c, v, nv_e, +bu / dld_a[alo], 0, False)
    _corr(e_corr, u, a, ahi, u, nu_c, v, nv_e, -bu / dld_a[ahi], nk_e, False)
    # E_v on the a faces (E_v: u edges x v cells)
    _corr(e_corr, v, a, alo, u, nu_e, v, nv_c, +bv / dld_a[alo], 0, False)
    _corr(e_corr, v, a, ahi, u, nu_e, v, nv_c, -bv / dld_a[ahi], nk_e, False)
    # E_a on the u faces (E_a: v edges x a cells; the H table column follows a, column 1 = cell alo)
    _corr(e_corr, a, u, ulo, v, nv_e, a, na_c, -bu / dld_u[ulo], 1, True)
    _corr(e_corr, a, u, uhi, v, nv_e, a, na_c, +bu / dld_u[uhi], 1, True)
    # E_a on the v faces (E_a: u edges x a cells)
    _corr(e_corr, a, v, vlo, u, nu_e, a, na_c, -bv / dld_v[vlo], 1, True)
    _corr(e_corr, a, v, vhi, u, nu_e, a, na_c, +bv / dld_v[vhi], 1, True)
    # H_v on the a cells (H_v: u cells x v edges)
    _corr(h_corr, v, a, alo - 1, u, nu_c, v, nv_e, +bu * ch / dl_a[alo - 1], 0, False)
    _corr(h_corr, v, a, ahi, u, nu_c, v, nv_e, -bu * ch / dl_a[ahi], nk_e - 1, False)
    # H_u on the a cells (H_u: u edges x v cells)
    _corr(h_corr, u, a, alo - 1, u, nu_e, v, nv_c, -bv * ch / dl_a[alo - 1], 0, False)
    _corr(h_corr, u, a, ahi, u, nu_e, v, nv_c, +bv * ch / dl_a[ahi], nk_e - 1, False)
    # H_a on the v cells (H_a: u cells x a edges; E table column 0 = edge alo)
    _corr(h_corr, a, v, vlo - 1, u, nu_c, a, na_e, -bu * ch / dl_v[vlo - 1], 0, True)
    _corr(h_corr, a, v, vhi, u, nu_c, a, na_e, +bu * ch / dl_v[vhi], 0, True)
    # H_a on the u cells (H_a: v cells x a edges)
    _corr(h_corr, a, u, ulo - 1, v, nv_c, a, na_e, +bv * ch / dl_u[ulo - 1], 0, True)
    _corr(h_corr, a, u, uhi, v, nv_c, a, na_e, -bv * ch / dl_u[uhi], 0, True)
    return e_corr, h_corr


def _tfsf_box(sc: Scene, n_total: int) -> dict | None:
    """What the TFSF box needs in the kernel: two 1D incident tables (on device) and the face
    correction descriptors.

    ``None`` means the scene has no TFSF source; then the correction kernels are never launched
    and a scene without TFSF is bitwise unchanged (one code path, the same approach as the
    absorber).

    The correction formulas were derived term by term from the difference equations in
    kernels/yee.cu. Direction and polarization produce no code branch: every sign is folded into
    the coefficient of the descriptor. The criterion is the empty-box vacuum test
    (tests/test_tfsf.py): inside the box the field equals the polarization component times the 1D
    table at every point, and leakage outside the box is at the float32 rounding level.

    Write the injection axis as a and the transverse axes as ``u=(a+1)%3``, ``v=(a+2)%3`` (cyclic
    right-handed order, so (∇×H)_u = ∂H_a/∂v - ∂H_v/∂a and the rest all hold). The incident field
    is E = (β_u û + β_v v̂)·E1D and H = (β_u v̂ - β_v û)·H1D. Checking the across-box differences
    term by term gives 16 kinds of correction (those with β = 0 are not launched):

    =========  ===================  ===========================================
    corrected  face                 coefficient
    =========  ===================  ===========================================
    E_u        a=alo / a=ahi edge   (+β_u, -β_u) / dl_dual_a, H table col 0 / N
    E_v        a=alo / a=ahi edge   (+β_v, -β_v) / dl_dual_a, H table col 0 / N
    E_a        u=ulo / u=uhi edge   (-β_u, +β_u) / dl_dual_u, H table follows a
    E_a        v=vlo / v=vhi edge   (-β_v, +β_v) / dl_dual_v, H table follows a
    H_v        a cell alo-1 / ahi   (+β_u, -β_u)·ch / dl_a, E table col 0 / N-1
    H_u        a cell alo-1 / ahi   (-β_v, +β_v)·ch / dl_a, E table col 0 / N-1
    H_a        v cell vlo-1 / vhi   (-β_u, +β_u)·ch / dl_v, E table follows a
    H_a        u cell ulo-1 / uhi   (+β_v, -β_v)·ch / dl_u, E table follows a
    =========  ===================  ===========================================

    On the E step the coefficient is additionally multiplied by that component's per-cell cb.
    (N = ahi-alo+1; H table column 0 is 3D cell alo-1, E table column 0 is 3D edge alo.)
    The old z-axis injection with Ex polarization is just the special case a=2, β_u=1, β_v=0, with
    the same coefficients term by term.
    """
    if not sc.tfsf_sources:
        return None
    if len(sc.tfsf_sources) != 1:
        raise NotImplementedError(
            f"only a single TFSF source is supported, got {len(sc.tfsf_sources)}")
    t = sc.tfsf_sources[0]
    _, ny, nz = sc.shape
    _tfsf_refuse_bad_placement(sc, t, n_total)

    a = t.axis
    oblique = getattr(t, "k_hat", None) is not None
    # E column count: normal incidence uses the box columns (edges alo..ahi), oblique incidence
    # uses the span of the box projected onto k̂
    nk_e = (int(t.n_col_e) if oblique
            else t.box_hi[a] - t.box_lo[a] + 1)

    # 1D auxiliary grid, precomputed (pure numpy; the zero-at-the-ends check lives in tables)
    # amp_scale: injection amplitude factor for oblique incidence (1.0 for normal incidence, and
    # multiplying by 1.0 is exact => bitwise unchanged)
    _a = float(getattr(t, "amp_scale", 1.0))
    ex_tab, hy_tab = tfsf1d.tables(
        t.dl1, t.eps1, sc.dt, t.inj, t.direction, t.ex_inc * _a, t.hy_inc * _a,
        n_total, t.col0, nk_e, sigma1=getattr(t, "sigma1", None))

    ch = sc.dt / cpml.MU_0
    nzp = _pitch(nz)
    strides = (ny * nzp, nzp, 1)   # pitched strides
    if oblique:
        # Oblique incidence: 24 kinds of face correction, with real-valued table column positions
        # (tfsf_oblique). The 16 normal-incidence entries are hand-written constants and are
        # swapped out wholesale; both paths feed the same kernel.
        # (This used to build both sets and let the oblique one overwrite the normal-incidence
        #  set. The descriptors are pure computation, so building one fewer changes no result.)
        e_corr, h_corr = tfsf_oblique.descriptors(t, sc.grid, strides, ch)
        tfsf_oblique.check_in_table(e_corr, hy_tab.shape[1], "E step")
        tfsf_oblique.check_in_table(h_corr, ex_tab.shape[1], "H step")
        for d in e_corr + h_corr:
            d["launch"] = grid_2d(int(d["np"]), int(d["nq"]))
    else:
        e_corr, h_corr = _tfsf_face_corrections(t, sc.grid, strides, nk_e, ch)
    return {
        "ex_tab": cp.asarray(ex_tab),          # (n_total+1, nk_e) float64
        "hy_tab": cp.asarray(hy_tab),          # (n_total, nk_e+1)
        "e_corr": e_corr,
        "h_corr": h_corr,
    }


def _mode_source(sc: Scene) -> list[dict] | None:
    """What a mode source needs in the kernel: four complex field profiles, two coefficient pairs
    and the indices, per injection entry.

    Returns a **list of injection entries** (``None`` = no mode source). The injection kernel
    accumulates (``+=``), so one entry = one kernel call on one TF/SF face; several faces mean one
    call per entry and the fields add by superposition, because both the Yee update and the TF/SF
    correction are linear, so the multi-source solution is the sum of the single-source solutions.
    (The S-matrix adjoint simulation of Autograd27 injects 4-5 ModeSources at once, on different
    faces and with different waveforms.) The profiles are uploaded once, split into two float32
    arrays for real and imaginary part; per step, each entry only sends the two scalars ``amp``.

    Several sources on the same face with the same time waveform are still **merged into one
    entry** (profiles added, the correction terms being exactly isomorphic; BraggGratings is two
    parallel waveguides in one x plane that do not overlap in y), and that old path stays bitwise
    unchanged. Merging requires the profiles not to overlap, since with overlap the summed profile
    no longer matches the bookkeeping of a single amp. Combinations that cannot be merged now fall
    back to per-entry injection instead of raising.
    """
    if not sc.mode_sources:
        return None
    m = sc.mode_sources[0]
    if len(sc.mode_sources) > 1:
        mergeable = all(
            o.axis == m.axis
            and o.plane_index == m.plane_index and o.direction == m.direction
            and np.array_equal(o.amp_e, m.amp_e)
            and np.array_equal(o.amp_h, m.amp_h)
            for o in sc.mode_sources[1:])
        if not mergeable:
            # Different face or different waveform: one entry per source, the kernel accumulates
            # entry by entry (linear superposition).
            return [_mode_source_entry(sc, x) for x in sc.mode_sources]
        # Merging is only on the table for the same face and the same waveform (where the profile
        # shapes are necessarily equal); with overlap the summed profile does not match the
        # bookkeeping of a single amp, so the original overlap check is kept.
        ov = np.count_nonzero(
            np.logical_and.reduce([np.abs(x.ey_inc) > 0 for x in sc.mode_sources]))
        if ov:
            # Overlapping profiles cannot be merged into one entry (the single-amp bookkeeping
            # does not add up), but per-entry injection is always correct: the injection kernel
            # is +=, and both the Yee update and the TF/SF correction are linear.
            # (In an adjoint simulation, several diffraction-order plane waves on the same face
            #  with the same waveform are exactly this case.)
            return [_mode_source_entry(sc, x) for x in sc.mode_sources]
        from dataclasses import replace as _replace
        m = _replace(
            m,
            ey_inc=sum(x.ey_inc for x in sc.mode_sources),
            ez_inc=sum(x.ez_inc for x in sc.mode_sources),
            hy_inc=sum(x.hy_inc for x in sc.mode_sources),
            hz_inc=sum(x.hz_inc for x in sc.mode_sources))
    return [_mode_source_entry(sc, m)]


def _mode_source_entry(sc: Scene, m) -> dict:
    """Kernel parameter table of one mode-source injection entry (see :func:`_mode_source`)."""
    ax = sc.grid.axes[m.axis]
    ks = m.plane_index
    kh = ks - 1 if m.direction > 0 else ks
    # The tangential pair is in cyclic order and the Yee curl keeps its form under a cyclic
    # permutation, so one sign table serves all three normals
    sy_h, sz_h, sy_e, sz_e = _MODE_SIGNS[m.direction]
    ch = sc.dt / cpml.MU_0 * float(1.0 / ax.dl[kh])
    ce = float(1.0 / ax.dl_dual[ks])

    def _dev(a):
        return (cp.asarray(np.ascontiguousarray(a.real), dtype=cp.float32),
                cp.asarray(np.ascontiguousarray(a.imag), dtype=cp.float32))

    ey_re, ey_im = _dev(m.ey_inc)
    ez_re, ez_im = _dev(m.ez_inc)
    hy_re, hy_im = _dev(m.hy_inc)
    hz_re, hz_im = _dev(m.hz_inc)
    return {
        "axis": np.int32(m.axis),
        "ks": np.int32(ks), "kh": np.int32(kh),
        "ey_re": ey_re, "ey_im": ey_im, "ez_re": ez_re, "ez_im": ez_im,
        "hy_re": hy_re, "hy_im": hy_im, "hz_re": hz_re, "hz_im": hz_im,
        "coef_h": (np.float32(sy_h * ch), np.float32(sz_h * ch)),
        "coef_e": (np.float32(sy_e * ce), np.float32(sz_e * ce)),
        "amp_e": np.ascontiguousarray(m.amp_e, dtype=np.complex128),
        "amp_h": np.ascontiguousarray(m.amp_h, dtype=np.complex128),
        # Broadband correction terms: one set of profiles plus one coefficient per node. With
        # none, this is an empty list, the for loop downstream never runs a single pass, and the
        # result is bitwise identical to before the change.
        "bb": [
            {"ey_re": _dev(m.bb_ey[k])[0], "ey_im": _dev(m.bb_ey[k])[1],
             "ez_re": _dev(m.bb_ez[k])[0], "ez_im": _dev(m.bb_ez[k])[1],
             "hy_re": _dev(m.bb_hy[k])[0], "hy_im": _dev(m.bb_hy[k])[1],
             "hz_re": _dev(m.bb_hz[k])[0], "hz_im": _dev(m.bb_hz[k])[1],
             "amp_e": np.ascontiguousarray(m.bb_amp_e[k], dtype=np.complex128),
             "amp_h": np.ascontiguousarray(m.bb_amp_h[k], dtype=np.complex128)}
            for k in range(0 if m.bb_ey is None else len(m.bb_ey))
        ],
    }


def _dipole_table(dipoles, n_total: int, ny: int, nz: int,
                  num_time_steps: int, half: bool) -> dict | None:
    """Flatten the stencils of a group of point dipoles into 1D and move the amplitude table to
    device memory.

    Electric dipoles are injected on the E step and take the half-step amplitude table
    (``half=True``); magnetic dipoles are injected on the H step and take the integer-step table
    (``half=False``, because the midpoint of an H step is ``n·dt``, see the comments in
    inject_dipole_h).
    """
    if not dipoles:
        return None
    def tab(d):
        return d.waveform.amp_half if half else d.waveform.amp_int

    nsteps1 = min(len(tab(d)) for d in dipoles)
    # The waveform is only sampled out to the **nominal** step count. Running longer means
    # resampling: silently padding with zeros happens to be harmless for a Gaussian pulse, but it
    # would quietly switch off a continuous-wave source, so stop here instead.
    if n_total > nsteps1 - 1:
        raise ValueError(
            f"asked to run {n_total} steps, but the dipole waveform only sampled {nsteps1} points"
            f" (nominal step count {num_time_steps}). Raise run_time or lower num_steps.")
    flat, comp, coefs, src_of = [], [], [], []
    for s, d in enumerate(dipoles):
        for (i, j, kk), cf in zip(d.indices, d.coef):
            flat.append(ravel_cells(int(i), int(j), int(kk), (ny, nz)))
            comp.append(d.component)
            coefs.append(cf)
            src_of.append(s)
    return {
        "n": np.int32(len(flat)),
        "nsteps1": np.int32(nsteps1),
        "flat": cp.asarray(np.asarray(flat, dtype=np.int32)),
        "comp": cp.asarray(np.asarray(comp, dtype=np.int32)),
        "coef": _f32(np.asarray(coefs)),
        "src_of": cp.asarray(np.asarray(src_of, dtype=np.int32)),
        "amp": _f64(np.stack([tab(d)[:nsteps1] for d in dipoles])),
    }


def _modulation_setup(sc: Scene, ny: int, nz: int, nzp: int) -> dict | None:
    """Time-varying medium: ``modulate_coeffs`` rewrites ca/cb of the modulated cells every step.

    The D form is folded into ca/cb (``kernels/modulation.cu``), so update_e itself needs no
    change. A scene without modulation returns ``None``, the kernel is never launched once and
    the result is bitwise unchanged (one path, the same approach as the absorber).
    """
    if not sc.any_modulation:
        return None
    mo = sc.modulation
    mo.validate(sc.shape)
    # PEC / mix / tensor are still refused: each of them rewrites the same ca/cb, and the
    # convention for combining them is not settled.
    # **Dispersion and loss are supported now**: the kernel uses the joint update form
    #   ca = (ε_inf^n - s)/(ε_inf^{n+1} + s + G),  cb = (dt/ε₀)/(...)
    # with ε_inf taken at step n in the numerator and at step n+1 in the denominator; s and G are
    # passed in per entry.
    if sc.any_pec or sc.any_dispersion_mix or sc.any_tensor:
        raise NotImplementedError(
            "a time-varying medium together with PEC/mix/tensor in one scene is not supported "
            "yet (ca/cb would be rewritten by two paths)")
    eps_maps = (sc.eps_ex, sc.eps_ey, sc.eps_ez)
    es = np.empty(mo.n_entry, dtype=np.float64)
    for c in range(3):
        m = mo.comp == c
        es[m] = eps_maps[c].ravel()[mo.cell[m]]
    # Per-entry s = σ dt/(2ε₀) and the implicit dispersion coupling G; lossless and non-dispersive
    # means all zeros, and the kernel's ca=(e0-0)/(e1+0+0) falls back to the original form exactly.
    _sd = np.zeros(mo.n_entry, dtype=np.float64)
    if sc.any_loss:
        sig_maps = (sc.sigma_ex, sc.sigma_ey, sc.sigma_ez)
        for c in range(3):
            m = mo.comp == c
            _sd[m] = (sig_maps[c].ravel()[mo.cell[m]]
                      * sc.dt / (2.0 * cpml.EPSILON_0))
    _gg = np.zeros(mo.n_entry, dtype=np.float64)
    if sc.any_dispersion:
        g_maps = sc.dispersion.g_dense(sc.shape)
        for c in range(3):
            m = mo.comp == c
            _gg[m] = g_maps[c].ravel()[mo.cell[m]]
    if float(np.min(es - np.abs(mo.amp))) <= 0.0:
        raise ValueError(
            f"modulated cells: min of ε_s - |amp| is {float(np.min(es - np.abs(mo.amp))):.4f}"
            " <= 0, so ε(t) would go non-positive, which is unstable")
    return {
        "n": np.int32(mo.n_entry),
        "amp": _f32(mo.amp),
        "cph": _f32(np.cos(mo.phase)),
        "sph": _f32(np.sin(mo.phase)),
        "eps": _f32(es),
        "sd": _f32(_sd),
        "gg": _f32(_gg),
        "comp": cp.asarray(np.ascontiguousarray(mo.comp, dtype=np.int32)),
        "cell": cp.asarray(_repitch_cells(mo.cell, ny, nz, nzp)),
        "om": 2.0 * np.pi * mo.freq,
        "dte": np.float32(sc.dt / cpml.EPSILON_0),
        "launch": grid_1d(mo.n_entry),
    }


def _dipole_table_c(dipoles, n_total: int, ny: int, nz: int,
                    num_time_steps: int, half: bool) -> dict | None:
    """Complex variant of :func:`_dipole_table`: the amplitude table becomes the full
    ``exp(iφ)·g(t)``.

    The stencil geometry (flat/comp/coef/src_of) is literally the same code as the real path.
    """
    if not dipoles:
        return None
    base = _dipole_table(dipoles, n_total, ny, nz, num_time_steps, half=half)
    tabs = []
    for d in dipoles:
        w = d.waveform
        t = w.amp_half_complex if half else w.amp_int_complex
        if t is None:
            raise ValueError(
                "the electric dipole waveform has no amp_half_complex (older npz files only "
                "stored the real part): the complex-field path injects the full complex "
                "amplitude, so re-export scene.npz")
        tabs.append(np.asarray(t, dtype=np.complex128))
    n1 = int(base["nsteps1"])
    stack = np.stack([t[:n1] for t in tabs])
    del base["amp"]
    base["amp_re"] = _f64(stack.real)
    base["amp_im"] = _f64(stack.imag)
    return base


class ModeLaunch(NamedTuple):
    """Return value of :func:`_mode_source_launch`: launch geometry and amplitude tables of one
    injection entry."""

    grid: tuple
    block: tuple
    dims: tuple          #: (axis, n1, n2, w2), the kernel's profile geometry
    amp_e: cp.ndarray    #: (n, 2) float64, integer-step amp (used by the H-side injection)
    amp_h: cp.ndarray    #: (n, 2) float64, half-step amp (used by the E-side injection)
    bb: list             #: ``(amp_e, amp_h)`` of each broadband correction term, empty if none


def _mode_source_launch(mode_src, nz, nzp) -> ModeLaunch:
    """Mode source: pad the field planes along z up to the pitch, then work out the launch
    geometry and the amplitude tables.

    The eight field planes are padded in place inside ``mode_src``. The amplitude table is
    (n, 2) float64; the (float) conversion in the kernel rounds the same way as the old
    ``np.float32(a.real)``.
    """
    def _pad(a):
        if nzp == nz or int(mode_src["axis"]) != 0:
            return a
        o = cp.zeros((a.shape[0], nzp), a.dtype)
        o[:, :nz] = a
        return o

    keys8 = ("ey_re", "ey_im", "ez_re", "ez_im",
             "hy_re", "hy_im", "hz_re", "hz_im")
    for key in keys8:
        mode_src[key] = _pad(mode_src[key])
    for _b in mode_src.get("bb", ()):           # broadband correction terms are padded too
        for key in keys8:
            _b[key] = _pad(_b[key])
    n1, w2 = (int(v) for v in mode_src["ey_re"].shape)
    n2 = w2 if int(mode_src["axis"]) != 0 else nz
    g2ms, b2ms = grid_2d(n1, w2)
    dims = (mode_src["axis"], np.int32(n1), np.int32(n2), np.int32(w2))
    amp = [_f64(np.ascontiguousarray(np.stack([a.real, a.imag], axis=1)))
           for a in (np.asarray(mode_src["amp_e"]),
                     np.asarray(mode_src["amp_h"]))]
    bb_amp = [(_f64(np.ascontiguousarray(
                   np.stack([b["amp_e"].real, b["amp_e"].imag], axis=1))),
               _f64(np.ascontiguousarray(
                   np.stack([b["amp_h"].real, b["amp_h"].imag], axis=1))))
              for b in mode_src.get("bb", ())]
    return ModeLaunch(g2ms, b2ms, dims, amp[0], amp[1], bb_amp)


def _mode_launch_tables(sc: Scene, d: Dims, broadband_ok: bool = True) -> dict:
    """Mode source tables plus one set of launch parameters per injection entry (several sources =
    per-entry injection; the kernel accumulates, so they superpose naturally).

    For normal x the second dimension of the profile is z, padded to (ny, nzp) (the kernel's
    o = o1*w2+o2); for the other normals the second dimension is not z and the row width is that
    axis' own n_t2.
    ``broadband_ok=False`` (the Bloch path): the complex injection kernel does not handle
    broadband mode sources yet, so any correction term fails closed. Returns
    ``{"mode_srcs", "ms_launch"}``, both None when there is no mode source.
    """
    mode_srcs = _mode_source(sc)
    ms_launch = None
    if mode_srcs is not None:
        ms_launch = [_mode_source_launch(m, d.nz, d.nzp) for m in mode_srcs]
        if not broadband_ok and any(launch.bb for launch in ms_launch):
            raise NotImplementedError(
                "the complex injection kernel of the Bloch path does not handle broadband mode "
                f"sources yet ({sum(len(launch.bb) for launch in ms_launch)} correction terms). "
                "Silently injecting only the base term would drop the profile corrections at both "
                "edges of the band, so stop instead.")
    return {"mode_srcs": mode_srcs, "ms_launch": ms_launch}


def _dipole_tables(sc: Scene, n_total: int, d: Dims,
                   complex_amp: bool = False) -> dict:
    """One table each for electric dipoles (injected on the E step, half-step table) and magnetic
    dipoles (injected on the H step, integer-step table), with the matching 1D launch geometry;
    when one class of dipole is absent, both of its entries are None.
    ``complex_amp=True`` goes through :func:`_dipole_table_c` (the full complex amplitude of the
    Bloch path). Returns ``{"dip", "dip_h", "g_dip", "g_dip_h"}``."""
    build = _dipole_table_c if complex_amp else _dipole_table
    ny, nzp = d.ny, d.nzp
    dip = build([dp for dp in sc.dipoles if not dp.magnetic],
                n_total, ny, nzp, sc.num_time_steps, half=True)
    dip_h = build([dp for dp in sc.dipoles if dp.magnetic],
                  n_total, ny, nzp, sc.num_time_steps, half=False)
    return {"dip": dip, "dip_h": dip_h,
            "g_dip": grid_1d(int(dip["n"]), block=64) if dip else None,
            "g_dip_h": grid_1d(int(dip_h["n"]), block=64) if dip_h else None}


def _modulation_phase_table(mod, n_total: int, dt: float):
    """Per-step cos/sin table of a time-varying medium,
    ``(n_total, 4) = [cos w0, sin w0, cos w1, sin w1]``, computed on the host in float64 and cast
    down to f32: bitwise identical to the old per-step np.float32(np.cos(...)).
    Returns None when there is no modulation (``mod is None``)."""
    if mod is None:
        return None
    _t0 = mod["om"] * (np.arange(n_total, dtype=np.float64) * dt)
    _t1 = mod["om"] * ((np.arange(n_total, dtype=np.float64) + 1.0) * dt)
    return _f32(np.ascontiguousarray(np.stack(
        [np.cos(_t0), np.sin(_t0), np.cos(_t1), np.sin(_t1)], axis=1)))
