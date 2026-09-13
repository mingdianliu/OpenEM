# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Equivalence criteria for symmetry folding (openem/fold.py).

The central claim: **stepping the folded half domain == stepping the full domain restricted to
the upper half, bitwise identical**, under a single precondition: the input itself is bitwise
symmetric. The tests build such an input themselves (upper half random, lower half mirrored
according to each component's edge/ctr position and the mirror eigenvalue σ), with no dependence
on tidy3d.

The feature set matches RingResonator: absorbers on all three axes + dispersion everywhere + an
x-normal mode source + a FieldMonitor spanning the symmetry plane. Both s=+1 and s=−1 are tested.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

pytest.importorskip("cupy")

from openem import fold, solver, waveform
from openem.device import Kernels
from openem.fold import expand_field
from openem.grid import MIRROR_SIGN as _SIG, YEE_ON_EDGE as _EDGE
from openem.grid import Axis, Grid
from openem.model import (
    AbsorberSlab, Dispersion, FieldMonitor, FluxTimeMonitor, ModeSource,
    PointDipole, Scene,
)

NX, NY, NZ = 20, 18, 24          # z is the fold axis, NZ has to be even
M = NZ // 2
DT = 2.0e-17
STEPS = 260


def _mirror_z(a: np.ndarray, comp: str, s: int) -> np.ndarray:
    """Mirror the upper-z values into the lower half, following the component's edge/ctr position
    and σ.

    An odd-parity (σ<0) edge-located component has to be zero on the symmetry-plane row (g=M). A
    physical field profile is like that anyway (tangential E on a PEC symmetry plane, say), but a
    random profile has to be zeroed by hand, otherwise the odd symmetry of the full-domain
    reference is broken from the start.
    """
    sig = _SIG[comp][2] * s
    edge = _EDGE[comp][2]
    out = a.copy()
    if edge and sig < 0:
        out[..., M] = 0.0
    for g in range(M):
        gm = (2 * M - g) if edge else (2 * M - 1 - g)
        if gm >= NZ:            # edge g=0 mirrors to 2M, a cell past the end that nobody uses
            continue
        out[..., g] = sig * out[..., gm]
    return out


def _mirror_eps_z(a: np.ndarray, comp: str) -> np.ndarray:
    """Mirroring ε/σ carries no sign (the material is a scalar), only the position mapping."""
    edge = _EDGE[comp][2]
    out = a.copy()
    for g in range(M):
        gm = (2 * M - g) if edge else (2 * M - 1 - g)
        if gm >= NZ:
            continue
        out[..., g] = out[..., gm]
    return out


def _make_scene(s: int, rng: np.ndarray | None = None) -> Scene:
    r = np.random.default_rng(7)

    # ---- grid: uniform on all three axes, symmetric in z ----
    def ax(n):
        return Axis(edges=np.linspace(0.0, n * 25e-9, n + 1),
                    boundary_lo="Absorber", boundary_hi="Absorber")
    grid = Grid(ax(NX), ax(NY), ax(NZ))

    # ---- ε: upper half random, lower half mirrored; conductivity σ unused (left as None) ----
    eps = {}
    for comp in ("Ex", "Ey", "Ez"):
        a = 1.0 + 1.5 * r.random((NX, NY, NZ))
        eps[comp] = _mirror_eps_z(np.asarray(a, np.float64), comp)

    # ---- absorbers: both ends of x/y plus both ends of z; the low-z slab mirrors the high-z one
    # The low-z slab has to be the **exact mirror** of the high-z one, and the two staggered
    # positions map differently:
    #   ctr (decay_half):  ctr g ↔ ctr n-1-g   → hi ctr {19..23} mirrors to {0..4}
    #   edge (decay_int):  edge g ↔ edge n-g   → hi edge {19..23} mirrors to {1..5}
    # That is, the mirrored edge set and the mirrored ctr set are **offset by one cell**, and a
    # lo slab with the same number of layers cannot cover both (brute-forcing all four profile
    # combinations came out asymmetric every time). The fix: one extra layer in the lo slab
    # (nlay = L+1, covering cells 0..5):
    #   edge 0    = dead cell on the PEC wall (tangential E is zeroed by pec and Hz on the wall
    #               does not affect any live cell), so any value will do;
    #   edge 1..5 = prof_i reversed (mirrors hi edge 23..19);
    #   ctr 0..4  = prof_h reversed, ctr 5 = 1.0 (the identity: its mirror, ctr 18, is not in
    #               the hi slab).
    L = 5
    prof_i = np.linspace(0.995, 0.9, L)
    prof_h = np.linspace(0.993, 0.88, L)
    lo_int = np.empty(L + 1)
    lo_int[0] = prof_i[-1]              # dead cell
    lo_int[1:] = prof_i[::-1]           # lo edge g = hi edge 24-g
    lo_half = np.empty(L + 1)
    lo_half[:L] = prof_h[::-1]          # lo ctr g = hi ctr 23-g
    lo_half[L] = 1.0                    # ctr 5 is the identity
    slabs = []
    for axis, n in ((0, NX), (1, NY)):
        slabs.append(AbsorberSlab(axis=axis, g0=0,
                                  decay_int=prof_i[::-1].copy(),
                                  decay_half=prof_h[::-1].copy()))
        slabs.append(AbsorberSlab(axis=axis, g0=n - L,
                                  decay_int=prof_i.copy(), decay_half=prof_h.copy()))
    slabs.append(AbsorberSlab(axis=2, g0=NZ - L,
                              decay_int=prof_i.copy(), decay_half=prof_h.copy()))
    slabs.append(AbsorberSlab(axis=2, g0=0,
                              decay_int=lo_int, decay_half=lo_half))

    # ---- dispersion: everywhere, one entry per (comp, cell), 1 pole; the weights are random in
    # the upper half and mirrored in the lower ----
    comp_l, cell_l, b_l = [], [], []
    bw = {c: _mirror_eps_z(0.05 * r.random((NX, NY, NZ)) + 0.01, c)
          for c in ("Ex", "Ey", "Ez")}
    for c_i, c in enumerate(("Ex", "Ey", "Ez")):
        idx = np.arange(NX * NY * NZ)
        comp_l.append(np.full(idx.size, c_i, np.int32))
        cell_l.append(idx.astype(np.int32))
        b_l.append(bw[c].ravel())
    comp = np.concatenate(comp_l)
    cell = np.concatenate(cell_l)
    q = -1.0e13 + 2.0e15j
    a_pole = (1 + q * DT / 2) / (1 - q * DT / 2)
    b = (np.concatenate(b_l) * (DT / 2) / (1 - q * DT / 2)).astype(np.complex128)
    disp = Dispersion(
        comp=comp, cell=cell,
        pole_ofs=np.arange(comp.size + 1, dtype=np.int32),
        am1=np.full(comp.size, a_pole - 1.0, np.complex128),
        b=b,
        g=np.ascontiguousarray(2.0 * b.real),
    )

    # ---- mode source: x-normal plane, with a profile satisfying the σ symmetry ----
    n_t = STEPS + 2
    t = np.arange(n_t) * DT
    amp = np.exp(-((t - 60 * DT) / (18 * DT)) ** 2) * np.exp(2j * np.pi * 2e14 * t)
    ey = r.random((NY, NZ)) + 1j * r.random((NY, NZ))
    ez = r.random((NY, NZ)) + 1j * r.random((NY, NZ))
    hy = r.random((NY, NZ)) + 1j * r.random((NY, NZ))
    hz = r.random((NY, NZ)) + 1j * r.random((NY, NZ))
    ey = _mirror_z(ey, "Ey", s); ez = _mirror_z(ez, "Ez", s)
    hy = _mirror_z(hy, "Hy", s); hz = _mirror_z(hz, "Hz", s)
    # The z=0 row of an edge-located profile sits on the low wall, and its mirror partner
    # (edge 2M) is the implicit ghost ≡ 0 on the high wall. A physical mode field is 0 on a
    # PEC-terminated wall anyway; leaving random values there lets update_h pick them up in the
    # window "after injection, before the next update_e zeroes them with pec", which makes the
    # full-domain reference asymmetric (located with an ablation matrix).
    ey[:, 0] = 0.0; hz[:, 0] = 0.0
    msrc = ModeSource(plane_index=8, direction=+1,
                      ey_inc=ey, ez_inc=ez, hy_inc=hy, hz_inc=hz,
                      amp_e=amp, amp_h=np.exp(
                          -((t + 0.5 * DT - 60 * DT) / (18 * DT)) ** 2
                      ) * np.exp(2j * np.pi * 2e14 * (t + 0.5 * DT)),
                      n_eff=1.5 + 0j)

    # ---- dipole: a mirror pair of stencil points (Ey, edge position) ----
    # σ(Ey) along z = +s: lower-half coefficient = +s × upper-half coefficient
    from openem import waveform as wf_mod
    n_t = STEPS + 2
    tt = np.arange(n_t) * DT
    wamp = (np.exp(-((tt - 50 * DT) / (15 * DT)) ** 2)).astype(np.float64)
    wf = wf_mod.Waveform(amp_int=wamp, amp_half=wamp.copy(),
                         amp_int_complex=wamp.astype(np.complex128), dt=DT)
    g_up = M + 3                      # upper-half edge point
    g_dn = 2 * M - g_up               # mirrored edge
    dip = PointDipole(
        component=1,                  # Ey
        indices=np.array([[9, 8, g_up], [9, 8, g_dn]], dtype=np.int32),
        coef=np.array([2.5e3, s * 2.5e3], dtype=np.float64),
        waveform=wf,
    )

    # ---- FieldMonitor spanning the symmetry plane ----
    # The box has to be "closed under expansion": the mirror of every lower-half point must land
    # inside the range the folded scene still stores. The edge mirror is g↔2M−g, and the smallest
    # g is 4, which mirrors to 20, so the high end of the box has to be ≥ 21 (to include 20).
    fmon = FieldMonitor(name="span", freqs=np.array([2.0e14, 2.2e14]),
                        origin=(6, 5, 4), box=(3, 4, 17))

    # ---- FluxTimeMonitor: the three folding shapes ----
    H = 25e-9
    ftm = [
        # tangential, spans the plane, covers both z sides -> a single monitor ×2
        FluxTimeMonitor(name="fx_full", axis=0, plane_index=10, normal_dir=1,
                        step_begin=0, step_end=STEPS, interval=1,
                        num_slots=STEPS,
                        t_bounds=((2.0 * H, 16.0 * H), (-H, (NZ + 1) * H)),
                        frac=0.3),
        # tangential, spans the plane, asymmetric -> split into two pieces and recombine
        FluxTimeMonitor(name="fx_split", axis=0, plane_index=7, normal_dir=1,
                        step_begin=0, step_end=STEPS, interval=1,
                        num_slots=STEPS,
                        t_bounds=((2.0 * H, 16.0 * H), (4.3 * H, 17.6 * H)),
                        frac=0.0),
        # normal = fold axis, plane in the lower half -> mirror the plane and flip the flux sign
        FluxTimeMonitor(name="fz_below", axis=2, plane_index=4, normal_dir=1,
                        step_begin=0, step_end=STEPS, interval=1,
                        num_slots=STEPS,
                        t_bounds=((2.0 * H, 16.0 * H), (3.0 * H, 15.0 * H)),
                        frac=0.4),
    ]

    return Scene(
        grid=grid, dt=DT, num_time_steps=STEPS, shutoff=0.0,
        eps_ex=eps["Ex"], eps_ey=eps["Ey"], eps_ez=eps["Ez"],
        pml={},
        dispersion=disp,
        symmetry=(0, 0, s),
        mode_sources=[msrc],
        dipoles=[dip],
        field_monitors=[fmon],
        flux_time_monitors=ftm,
        absorbers=slabs,
    )


@pytest.fixture(scope="module")
def kernels():
    return Kernels()


@pytest.mark.parametrize("s", [+1, -1])
def test_fold_bitwise_equals_full_upper_half(kernels, s):
    sc = _make_scene(s)
    folded, meta = fold.fold_scene(sc)
    assert folded.shape == (NX, NY, NZ - M)
    assert len(folded.absorbers) == 5           # the low-z slab is dropped
    assert folded.dispersion.n_entry == sc.dispersion.n_entry // 2

    r_full = solver.run(sc, num_steps=STEPS, use_shutoff=False,
                        kernels=kernels, verbose=False, return_fields=True)
    r_fold = solver.run(folded, num_steps=STEPS, use_shutoff=False,
                        kernels=kernels, verbose=False, return_fields=True)

    # ---- fields bitwise: folded domain == upper half of the full domain ----
    for name in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz"):
        up = r_full.fields[name][:, :, M:]
        fd = r_fold.fields[name]
        assert up.shape == fd.shape
        assert np.array_equal(up, fd), (
            f"s={s} {name}: folded vs full upper half not bitwise, "
            f"max|Δ|={np.abs(up - fd).max():.3e}, "
            f"max|full|={np.abs(up).max():.3e}")

    # ---- monitors: mirror-expanded folded phasors == full-domain phasors, bitwise ----
    (fold_ax,) = meta.folds
    o_full, b_full = meta.field_monitors_full["span"]
    ph_full = r_full.field_phasors["span"]      # (6, nf, ni, nj, nk)
    ph_fold = r_fold.field_phasors["span"]
    for c_i, comp in enumerate(("Ex", "Ey", "Ez", "Hx", "Hy", "Hz")):
        got = expand_field(ph_fold[c_i], comp, fold_ax,
                           lo_folded=0, lo_full=o_full[2], n_full=b_full[2])
        want = ph_full[c_i]
        assert np.array_equal(want, got), (
            f"s={s} monitor {comp}: expansion not bitwise, "
            f"max|Δ|={np.abs(want - got).max():.3e}")


def test_fold_fail_closed(kernels):
    sc = _make_scene(+1)
    with pytest.raises(NotImplementedError):
        fold.fold_scene(replace(sc, bloch_k=(0.1, 0, 0)))
    with pytest.raises(ValueError):
        fold.fold_scene(replace(sc, symmetry=(0, 0, 0)))


@pytest.mark.parametrize("s", [+1, -1])
def test_fold_flux_time(kernels, s):
    """Folding flux_time: a normal mirror flips the sign, a tangential one is ×2 or split in two,
    and the linear recombination == the full domain.

    The flux is an f64 sum over f32 fields that are themselves bitwise identical, and folding only
    changes the association order of the sum and the trapezoid split, so the deviation is at the
    1e-13 relative level: a tight allclose rather than bitwise.
    """
    sc = _make_scene(s)
    folded, meta = fold.fold_scene(sc)
    assert set(meta.flux_time_map) == {"fx_full", "fx_split", "fz_below"}
    assert len(meta.flux_time_map["fx_split"]) == 2      # split in two
    assert meta.flux_time_map["fx_full"][0][1] == 2.0    # ×2
    assert meta.flux_time_map["fz_below"][0][1] == -1.0  # sign flip

    r_full = solver.run(sc, num_steps=STEPS, use_shutoff=False,
                        kernels=kernels)
    r_fold = solver.run(folded, num_steps=STEPS, use_shutoff=False,
                        kernels=kernels)
    rec = fold.expand_flux_time(meta, r_fold.flux_time)
    for name, want in r_full.flux_time.items():
        got = rec[name]
        assert want.shape == got.shape, name
        scale = float(np.max(np.abs(want))) or 1.0
        np.testing.assert_allclose(got, want, rtol=5e-9,
                                   atol=1e-12 * scale, err_msg=name)


# ───────────────────── two mirror walls (periodic axis + s=+1) ─────────────────────

MY = NY // 2


def _mirror_eps_y_wrap(a: np.ndarray, comp: str) -> np.ndarray:
    """Mirror ε about the central y edge; edge row g=0 is its own mirror on the wrap-around wall,
    so it is kept."""
    edge = _EDGE[comp][1]
    out = a.copy()
    for g in range(MY):
        gm = (2 * MY - g) if edge else (2 * MY - 1 - g)
        if gm >= NY:
            continue                     # edge g=0 ↔ edge NY ≡ 0 (self-mirroring)
        out[:, g, :] = out[:, gm, :]
    return out


def _make_scene_wrap() -> Scene:
    """Periodic in y with s=+1: eps mirrored in y, dispersion everywhere mirrored in y, and a
    mirror pair of Ex dipoles."""
    r = np.random.default_rng(11)

    def ax(n, periodic=False):
        b = "Periodic" if periodic else "Absorber"
        return Axis(edges=np.linspace(0.0, n * 25e-9, n + 1),
                    boundary_lo=b, boundary_hi=b)
    grid = Grid(ax(NX), ax(NY, periodic=True), ax(NZ))

    eps = {}
    for comp in ("Ex", "Ey", "Ez"):
        a = 1.0 + 1.5 * r.random((NX, NY, NZ))
        eps[comp] = _mirror_eps_y_wrap(np.asarray(a, np.float64), comp)

    L = 5
    prof_i = np.linspace(0.995, 0.9, L)
    prof_h = np.linspace(0.993, 0.88, L)
    slabs = []
    for axis, n in ((0, NX), (2, NZ)):
        slabs.append(AbsorberSlab(axis=axis, g0=0,
                                  decay_int=prof_i[::-1].copy(),
                                  decay_half=prof_h[::-1].copy()))
        slabs.append(AbsorberSlab(axis=axis, g0=n - L,
                                  decay_int=prof_i.copy(),
                                  decay_half=prof_h.copy()))

    # one-pole dispersion everywhere, weights mirrored in y
    comp_l, cell_l, b_l = [], [], []
    bw = {c: _mirror_eps_y_wrap(0.05 * r.random((NX, NY, NZ)) + 0.01, c)
          for c in ("Ex", "Ey", "Ez")}
    for c_i, c in enumerate(("Ex", "Ey", "Ez")):
        idx = np.arange(NX * NY * NZ)
        comp_l.append(np.full(idx.size, c_i, np.int32))
        cell_l.append(idx.astype(np.int32))
        b_l.append(bw[c].ravel())
    comp = np.concatenate(comp_l)
    cell = np.concatenate(cell_l)
    q = -1.0e13 + 2.0e15j
    a_pole = (1 + q * DT / 2) / (1 - q * DT / 2)
    b = (np.concatenate(b_l) * (DT / 2) / (1 - q * DT / 2)
         ).astype(np.complex128)
    disp = Dispersion(
        comp=comp, cell=cell,
        pole_ofs=np.arange(comp.size + 1, dtype=np.int32),
        am1=np.full(comp.size, a_pole - 1.0, np.complex128),
        b=b, g=np.ascontiguousarray(2.0 * b.real),
    )

    # a mirror pair of Ex dipoles (Ex is ctr-located along y, σ = +s = +1)
    n_t = STEPS + 2
    tt = np.arange(n_t) * DT
    wamp = (np.exp(-((tt - 50 * DT) / (15 * DT)) ** 2)).astype(np.float64)
    wf = waveform.Waveform(amp_int=wamp, amp_half=wamp.copy(),
                           amp_int_complex=wamp.astype(np.complex128), dt=DT)
    g_up = MY + 3
    g_dn = 2 * MY - g_up        # Ex is edge-located along y
    dip = PointDipole(
        component=0,
        indices=np.array([[9, g_up, 12], [9, g_dn, 12]], dtype=np.int32),
        coef=np.array([2.5e3, 2.5e3], dtype=np.float64),
        waveform=wf,
    )

    fmon = FieldMonitor(name="span", freqs=np.array([2.0e14]),
                        origin=(6, 0, 4), box=(3, NY, 5))

    return Scene(
        grid=grid, dt=DT, num_time_steps=STEPS, shutoff=0.0,
        eps_ex=eps["Ex"], eps_ey=eps["Ey"], eps_ez=eps["Ez"],
        pml={},
        dispersion=disp,
        symmetry=(0, 1, 0),
        dipoles=[dip],
        field_monitors=[fmon],
        absorbers=slabs,
    )


def test_fold_wrap_bitwise(kernels):
    """Two-mirror-wall theorem: the folded domain == full-domain rows [m..n−1] plus row 0, bitwise."""
    sc = _make_scene_wrap()
    res_full = solver.run(sc, num_steps=STEPS, use_shutoff=False,
                          kernels=kernels, verbose=False, return_fields=True)
    folded, meta = fold.fold_scene(sc)
    fa = meta.folds[0]
    assert fa.wrap and fa.axis == 1
    m = fa.m
    res_f = solver.run(folded, num_steps=STEPS, use_shutoff=False,
                       kernels=kernels, verbose=False, return_fields=True)
    for c in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz"):
        full = res_full.fields[c]
        want = np.concatenate([full[:, m:, :], full[:, 0:1, :]], axis=1)
        got = res_f.fields[c]
        assert got.shape == want.shape, (c, got.shape, want.shape)
        if not np.array_equal(want.view(np.uint32), got.view(np.uint32)):
            d = np.abs(want - got)
            bad = np.unravel_index(np.argmax(d), d.shape)
            raise AssertionError(f"{c} not bitwise: max|Δ|={d.max():.3e} at {bad}")


# ───────────────────── folding mix (slanted-interface dispersion) ─────────────────────

def _make_mix(sc: Scene, s: int):
    """Add a band of mix entries, mirror-symmetric in z, to the scene from _make_scene."""
    from openem.model import DispersionMix
    r = np.random.default_rng(23)
    cells, comps = [], []
    beta_v, ei_v, zi_v = [], [], []
    # one band at z ∈ [M+2, M+6) in the upper half plus its mirror (ctr mapped per component)
    for c_i, c in enumerate(("Ex", "Ey", "Ez")):
        edge = _EDGE[c][2]
        for i in range(6, 10):
            for j in range(5, 9):
                for g_up in range(M + 2, M + 6):
                    g_dn = (2 * M - g_up) if edge else (2 * M - 1 - g_up)
                    b_ = 0.2 + 0.6 * r.random()
                    e_ = 2.0 + r.random()
                    z_ = 1.0 / e_ + 0.05 * r.random()
                    for g in (g_up, g_dn):
                        comps.append(c_i)
                        cells.append((i * NY + j) * NZ + g)
                        beta_v.append(b_)
                        ei_v.append(e_)
                        zi_v.append(z_)
    n_e = len(cells)
    q = -2.0e13 + 1.5e15j
    a_pole = (1 + q * DT / 2) / (1 - q * DT / 2)
    rb = 0.03 * r.random(n_e) + 0.01
    # a mirror pair must share the pole weight: entries are generated in pairs (2k, 2k+1) = (up, dn)
    rb = np.repeat(rb[::2], 2)[:n_e]
    pb = (rb * (DT / 2) / (1 - q * DT / 2)).astype(np.complex128)
    mx = DispersionMix(
        comp=np.asarray(comps, np.int32),
        cell=np.asarray(cells, np.int32),
        p_ofs=np.arange(n_e + 1, dtype=np.int32),
        pa=np.full(n_e, a_pole - 1.0, np.complex128),
        pb=pb,
        q_ofs=np.arange(n_e + 1, dtype=np.int32),
        qa=np.full(n_e, a_pole - 1.0, np.complex128),
        qb=(pb * 0.5).astype(np.complex128),
        beta=np.asarray(beta_v, np.float64),
        eps_inf=np.asarray(ei_v, np.float64),
        zeta_inf=np.asarray(zi_v, np.float64),
    )
    return replace(sc, dispersion_mix=mx)


@pytest.mark.parametrize("s", [+1, -1])
def test_fold_mix_bitwise(kernels, s):
    # the same (comp,cell) cannot be in both dispersion and mix (the model checks this), and
    # folding of dispersion already has its own theorem test, so this one focuses on mix
    sc = _make_mix(replace(_make_scene(s), dispersion=None), s)
    res_full = solver.run(sc, num_steps=STEPS, use_shutoff=False,
                          kernels=kernels, verbose=False, return_fields=True)
    folded, meta = fold.fold_scene(sc)
    m = meta.folds[0].m
    res_f = solver.run(folded, num_steps=STEPS, use_shutoff=False,
                       kernels=kernels, verbose=False, return_fields=True)
    for c in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz"):
        want = res_full.fields[c][:, :, m:]
        got = res_f.fields[c]
        assert got.shape == want.shape
        if not np.array_equal(want.view(np.uint32), got.view(np.uint32)):
            d = np.abs(want - got)
            raise AssertionError(f"{c} not bitwise: max|Δ|={d.max():.3e}")
