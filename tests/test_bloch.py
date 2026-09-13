# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Criteria for the Bloch (k≠0) complex-field path.

The first two of the three layers of criteria are locked here (seconds on a GPU):

1. **k=0 equivalence**: the solution of the complex path = the real solution from injecting Re,
   plus i times the real solution from injecting Im. At k=0 every update operator is real, so
   this is an identity rather than an approximation (in floating point the two differ at the
   level of FMA contraction).
2. **The k=0.5 antiperiodic shortcut**: exp(iπ) = −1 is still real, so the k=0.5 complex solution
   on a domain of length L has to equal the real solution on a 2L domain with "the original
   dipole plus a sign-flipped dipole shifted by L", restricted to the first half. This holds
   discretely, cell by cell, since both sides use the same Yee differences.
3. **1D discrete Bloch dispersion**: on a uniform vacuum chain the allowed modes are
   k_m = 2π(m+bv)/L, with frequencies given by the discrete dispersion relation
   sin(ωΔt/2) = S·|sin(k_m·Δx/2)| (S=cΔt/Δx). That has a closed form and moves with bv, so it
   verifies that the **magnitude** of the wrap-around phase is used correctly. The **sign** of
   the phase is settled by the third layer, which compares the complex series of sim_1..11
   against the reference output.

The third layer lives outside this file.
"""

import numpy as np
import os
import pytest

cp = pytest.importorskip("cupy")

from openem import cpml, serialize, solver
from openem.grid import Axis, Grid
from openem.model import (FieldTimeMonitor, FluxMonitor, FluxTimeMonitor,
                             PointDipole, Scene)
from openem.waveform import Waveform

C0 = cpml.C_0


# ------------------------------------------------------------ building blocks

def _axis(n: int, dl: float, kind: str) -> Axis:
    edges = np.arange(n + 1, dtype=np.float64) * dl
    return Axis(edges, kind, kind)


def _wave(num_steps: int, dt: float, f0: float, fwidth: float, phase: float,
          real_only: bool = False) -> Waveform:
    """A complex Gaussian pulse g(t)·exp(iφ), the same shape as GaussianPulse (none of the
    criteria depend on the exact spectrum)."""
    n = np.arange(num_steps + 1, dtype=np.float64)
    sig = 1.0 / (2.0 * np.pi * fwidth)
    t0 = 5.0 * sig

    def g(t):
        env = np.exp(-0.5 * ((t - t0) / sig) ** 2)
        if real_only:
            return (env * np.cos(2 * np.pi * f0 * (t - t0))).astype(np.complex128)
        return env * np.exp(-2j * np.pi * f0 * (t - t0)) * np.exp(1j * phase)

    a_int, a_half = g(n * dt), g((n + 0.5) * dt)
    return Waveform(amp_int=a_int.real, amp_half=a_half.real,
                    amp_int_complex=a_int, dt=dt, amp_half_complex=a_half)


def _dip(comp: int, ijk, dV: float, wf: Waveform, magnetic: bool) -> PointDipole:
    return PointDipole(component=comp, indices=np.array([ijk], dtype=np.int64),
                       coef=np.array([1.0 / dV]), waveform=wf, magnetic=magnetic)


def _scene(nx=10, ny=8, nz=32, dl=5e-8, n_pml=8, steps=400,
           bloch_k=(0.0, 0.0, 0.0), phases=(0.7, 2.1), real_only=False,
           x_kind="Periodic", monitors=True) -> Scene:
    """A small vacuum scene, periodic (or Bloch) in x/y with PML along z: one magnetic Hz dipole
    and one electric Ex dipole."""
    dt = 0.9 / (C0 * np.sqrt(3.0) / dl)
    ax = _axis(nx, dl, x_kind)
    ay = _axis(ny, dl, "Periodic")
    az = _axis(nz, dl, "PML")
    grid = Grid(ax, ay, az)
    eps = np.ones(grid.shape, dtype=np.float64)
    params = cpml.PMLParams(num_layers=n_pml)
    pml = {(2, s): cpml.build(params, dl, s, dt) for s in ("lo", "hi")}
    f0, fw = 1.5e14, 0.8e14
    dV = dl ** 3
    dips = [
        _dip(2, (nx // 2, ny // 2, nz // 2), dV,
             _wave(steps, dt, f0, fw, phases[0], real_only), magnetic=True),
        _dip(0, (nx // 3, ny // 3, nz // 2 - 2), dV,
             _wave(steps, dt, f0, fw, phases[1], real_only), magnetic=False),
    ]
    tmons = []
    if monitors:
        tmons = [FieldTimeMonitor(
            name="t0", comps=(0, 5), origin=(1, 1, nz // 2 - 3), box=(3, 3, 3),
            step_begin=0, step_end=steps + 1, interval=1, num_slots=steps + 1)]
    return Scene(grid=grid, dt=dt, num_time_steps=steps, shutoff=0.0,
                 eps_ex=eps, eps_ey=eps.copy(), eps_ez=eps.copy(), pml=pml,
                 dipoles=dips, field_time_monitors=tmons, bloch_k=bloch_k)


@pytest.fixture(scope="module")
def kernels():
    from openem.device import Kernels
    return Kernels()


# ------------------------------------------------------------ compilation self-check

def test_kernel_count_and_complex_variants(kernels):
    """A baseline of 23 kernels plus 5 complex variants = 28; every name has to compile.

    23 = the original 22 + dft_phase_table − absorb_disp_hist (folded into dispersion_post)
    + dispersion_step (post+pre fused).
    """
    names = kernels.names
    # 37 + 4 (update_e_ade_lean / update_e_lut_shell / update_h_lean / update_h_shell:
    # the interior/boundary-shell split, measured)
    assert len(names) == 59, names  # +4: sample_time_batch(P27), dft_*_sub(P33),
                                    #     absorb_batch(P36)
                                    # +2: update_he_fused / update_e_lut_edge
                                    #     (P41 H/E fusion, off by default, measured)
                                    # +2: inject_mode_e_c / inject_mode_h_c
                                    #     (complex injection on a single oblique TF/SF face)
                                    # +3: dispersion_step[_dec[_we]]_p2
                                    #     (P48 2-pole bucket, float2 fast path)
    assert "absorb_batch" in names
    assert "dft_phase_table" in names
    assert "dispersion_step" in names
    for _n in ("update_e_ade_lean", "update_e_lut_shell",
               "update_h_lean", "update_h_shell",
               "update_he_fused", "update_e_lut_edge", "update_e_ade",
               "dispersion_step_dec_we"):
        assert _n in names, _n
    assert "step_advance" in names       # device-side step counter
    assert "absorb_disp_hist" not in names
    for n in ("update_h_c", "update_e_c", "inject_dipole_c",
              "inject_dipole_h_c", "sample_time_box_c"):
        assert n in names


# ------------------------------------------------------------ mask phases

def test_index_tables_bloch_phase():
    a = _axis(8, 1e-7, "Periodic")
    t0 = a.index_tables()
    assert t0["mnx"].dtype == np.float64 and np.all(t0["mnx"] == 1.0)
    k = 0.3
    t = a.index_tables(bloch_k=k)
    assert t["mnx"].dtype == np.complex128
    assert np.allclose(t["mnx"][:-1], 1.0) and np.allclose(t["mpv"][1:], 1.0)
    # u(x+L) = u(x)·exp(+i2πk): a forward wrap multiplies by e^{+i2πk}, a backward one by e^{−i2πk}
    assert np.isclose(t["mnx"][-1], np.exp(2j * np.pi * k))
    assert np.isclose(t["mpv"][0], np.exp(-2j * np.pi * k))
    with pytest.raises(ValueError):
        _axis(8, 1e-7, "PML").index_tables(bloch_k=0.5)


# ------------------------------------------------------------ serialize v9

def test_serialize_roundtrip_bloch(tmp_path):
    sc = _scene(bloch_k=(0.125, 0.25, 0.0))
    p = tmp_path / "scene.npz"
    serialize.save(sc, p)
    back = serialize.load(p)
    assert back.bloch_k == (0.125, 0.25, 0.0)
    for d0, d1 in zip(sc.dipoles, back.dipoles):
        np.testing.assert_array_equal(d0.waveform.amp_int_complex,
                                      d1.waveform.amp_int_complex)
        np.testing.assert_array_equal(d0.waveform.amp_half_complex,
                                      d1.waveform.amp_half_complex)


# ------------------------------------------------------------ criterion 1: k=0 equivalence

def test_k0_complex_equals_real_re_plus_im(kernels):
    """At k=0 the complex path = the real path (Re injection) + i times the real path (Im)."""
    from dataclasses import replace

    sc = _scene(steps=400)
    res_c = solver.run(sc, kernels=kernels, verbose=False, use_shutoff=False,
                       return_fields=True, force_complex=True)
    res_re = solver.run(sc, kernels=kernels, verbose=False, use_shutoff=False,
                        return_fields=True)

    def _im_dip(d):
        w = d.waveform
        wi = Waveform(amp_int=w.amp_int_complex.imag,
                      amp_half=w.amp_half_complex.imag,
                      amp_int_complex=w.amp_int_complex.imag.astype(np.complex128),
                      dt=w.dt,
                      amp_half_complex=w.amp_half_complex.imag.astype(np.complex128))
        return replace(d, waveform=wi)

    sc_im = replace(sc, dipoles=[_im_dip(d) for d in sc.dipoles])
    res_im = solver.run(sc_im, kernels=kernels, verbose=False, use_shutoff=False,
                        return_fields=True)

    for name in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz"):
        c, r, i = res_c.fields[name], res_re.fields[name], res_im.fields[name]
        scale = max(float(np.max(np.abs(r))), float(np.max(np.abs(i))), 1e-30)
        dre = float(np.max(np.abs(c.real - r))) / scale
        dim = float(np.max(np.abs(c.imag - i))) / scale
        assert dre < 1e-6 and dim < 1e-6, (name, dre, dim)

    # the same identity for the time-domain monitor (including the timing of the two-launch H
    # average)
    tc = res_c.time_samples["t0"]
    tr, ti = res_re.time_samples["t0"], res_im.time_samples["t0"]
    scale = max(float(np.max(np.abs(tr))), 1e-30)
    assert float(np.max(np.abs(tc.real - tr))) / scale < 1e-6
    assert float(np.max(np.abs(tc.imag - ti))) / scale < 1e-6


# ------------------------------------------------------------ criterion 2a: k=0.5 antiperiodic

def test_k05_equals_antiperiodic_doubled_domain(kernels):
    """The k=0.5 complex solution on an L domain == the real solution on a 2L domain with "the
    source plus a sign-flipped source shifted by L", restricted to the first half."""
    from dataclasses import replace

    nx, steps = 16, 600
    sc = _scene(nx=nx, ny=4, steps=steps, real_only=True,
                x_kind="BlochBoundary", bloch_k=(0.5, 0.0, 0.0), monitors=False)
    res_c = solver.run(sc, kernels=kernels, verbose=False, use_shutoff=False,
                       return_fields=True)

    # the 2L domain **reuses** the L-domain dipoles as they are (bitwise identical positions and
    # waveforms) and adds the sign-flipped image shifted by L
    sc2 = _scene(nx=2 * nx, ny=4, steps=steps, real_only=True, monitors=False)
    shifted = []
    for d in sc.dipoles:
        idx = d.indices.copy()
        idx[:, 0] += nx
        shifted.append(replace(d, indices=idx, coef=-d.coef))
    sc2 = replace(sc2, dipoles=list(sc.dipoles) + shifted)
    res_d = solver.run(sc2, kernels=kernels, verbose=False, use_shutoff=False,
                       return_fields=True)

    for name in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz"):
        c = res_c.fields[name]
        d = res_d.fields[name]
        scale = max(float(np.max(np.abs(d))), 1e-30)
        # the 2L solution has to be antiperiodic on its own (a self-check on the construction,
        # not the criterion itself)
        anti = float(np.max(np.abs(d[:nx] + d[nx:]))) / scale
        diff = float(np.max(np.abs(c.real - d[:nx]))) / scale
        leak = float(np.max(np.abs(c.imag))) / scale
        assert anti < 1e-6, (name, anti)
        assert diff < 1e-5, (name, diff)
        assert leak < 1e-6, (name, leak)


# ------------------------------------------------------------ criterion 2b: 1D discrete dispersion

def test_1d_discrete_bloch_dispersion(kernels):
    """The spectral peaks of a uniform vacuum chain have to land on the discrete Bloch dispersion
    ω(k_m), within ±1 bin, with zero parameters in the closed form.

    At bv=0.25, k_m = 2π(m+0.25)/L, and the four lowest modes are tens of bins away from the
    bv=0 (purely periodic) modes at 0, 187.5 THz and so on; get the magnitude of the phase wrong
    and the peaks move wholesale.
    """
    nx, dl, bv, steps = 32, 5e-8, 0.25, 16384
    dt = 0.5 * dl / C0
    grid = Grid(_axis(nx, dl, "BlochBoundary"),
                _axis(1, dl, "Periodic"), _axis(1, dl, "Periodic"))
    eps = np.ones(grid.shape, dtype=np.float64)
    wf = _wave(steps, dt, f0=1.6e14, fwidth=1.2e14, phase=0.4)
    dip = _dip(2, (5, 0, 0), dl ** 3, wf, magnetic=False)
    mon = FieldTimeMonitor(name="p", comps=(2,), origin=(17, 0, 0), box=(1, 1, 1),
                           step_begin=0, step_end=steps + 1, interval=1,
                           num_slots=steps + 1)
    sc = Scene(grid=grid, dt=dt, num_time_steps=steps, shutoff=0.0,
               eps_ex=eps, eps_ey=eps.copy(), eps_ez=eps.copy(), pml={},
               dipoles=[dip], field_time_monitors=[mon],
               bloch_k=(bv, 0.0, 0.0))
    res = solver.run(sc, kernels=kernels, verbose=False, use_shutoff=False)
    series = res.time_samples["p"][0, :, 0, 0, 0].astype(np.complex128)

    n = series.size
    F = np.abs(np.fft.fft(series * np.hanning(n)))
    freqs = np.fft.fftfreq(n, dt)
    df = 1.0 / (n * dt)

    L = nx * dl
    S = C0 * dt / dl

    def f_disc(m):
        km = 2.0 * np.pi * (m + bv) / L
        return np.arcsin(S * abs(np.sin(km * dl / 2.0))) / (np.pi * dt)

    targets = sorted(f_disc(m) for m in (0, -1, 1, -2))
    fmax = float(F.max())
    for fa in targets:
        win = np.flatnonzero(np.abs(np.abs(freqs) - fa) < 5 * df)
        j = win[int(np.argmax(F[win]))]
        assert F[j] > 0.05 * fmax, (fa, F[j] / fmax)
        # parabolic refinement to sub-bin resolution
        y0, y1, y2 = F[j - 1], F[j], F[j + 1]
        delta = 0.5 * float(y0 - y2) / float(y0 - 2 * y1 + y2)
        f_meas = abs(freqs[j]) + delta * df * np.sign(freqs[j])
        assert abs(abs(f_meas) - fa) <= df, (fa, f_meas, df)
    # no energy is allowed at the purely periodic (bv=0) modes: if the phase were dropped and the
    # run silently went ahead at k=0, a peak would show up right here
    f_periodic = f_disc(1.0 - bv)  # m+bv=1 -> the first purely periodic mode
    win = np.flatnonzero(np.abs(np.abs(freqs) - f_periodic) < 3 * df)
    assert float(F[win].max()) < 0.02 * fmax


# ------------------------------------------------------------ fail closed

def test_bloch_refuses_outside_subset(kernels):
    from dataclasses import replace

    sc = _scene(bloch_k=(0.25, 0.0, 0.0), x_kind="BlochBoundary")
    # FluxMonitor / FieldMonitor are supported (complex colocated flux / frequency-domain
    # accumulation), so it just has to run through
    ok = replace(sc, flux_monitors=[FluxMonitor(
        name="f", axis=2, plane_index=5, normal_dir=1,
        freqs=np.array([2e14]), source_spectrum=np.array([1.0 + 0j]))])
    res = solver.run(ok, kernels=kernels, verbose=False)
    assert "f" in res.phasors and res.phasors["f"].data.shape[0] == 4
    # FluxTimeMonitor still has no complex variant, so it has to fail closed
    ftm = replace(sc, flux_time_monitors=[FluxTimeMonitor(
        name="ft", axis=2, plane_index=5, normal_dir=1, step_begin=0, step_end=10,
        interval=1, num_slots=10, t_bounds=((0.0, 1e-6), (0.0, 1e-6)), frac=0.0)])
    with pytest.raises(NotImplementedError, match="FluxTimeMonitor"):
        solver.run(ftm, kernels=kernels, verbose=False)
    with pytest.raises(NotImplementedError, match="dipole"):
        solver.run(replace(sc, dipoles=[]), kernels=kernels, verbose=False)


def test_old_npz_electric_dipole_refused_in_complex_path(tmp_path, kernels):
    """An electric dipole read back from an old npz (no complex tables) has to be rejected loudly
    on the complex path, not run silently on the real part."""
    from dataclasses import replace

    sc = _scene(bloch_k=(0.25, 0.0, 0.0), x_kind="BlochBoundary")
    stripped = []
    for d in sc.dipoles:
        w = d.waveform
        stripped.append(replace(d, waveform=Waveform(
            amp_int=w.amp_int, amp_half=w.amp_half,
            amp_int_complex=w.amp_int.astype(np.complex128), dt=w.dt)))
    with pytest.raises(ValueError, match="amp_half_complex"):
        solver.run(replace(sc, dipoles=stripped), kernels=kernels, verbose=False)


# ------------------------------------------------- extraction side (needs tidy3d + the archive)

def test_bandstructure_sim1_extracts_bloch_k():
    td = pytest.importorskip("tidy3d")
    from openem import scene as scene_mod

    TABLEA = os.environ.get("OPENEM_TABLEA", "/nonexistent")  # external reference archive
    sc = scene_mod.from_file(f"{TABLEA}/Bandstructure/sims/sim_1/simulation.json")
    assert sc.bloch_k == (0.125, 0.0, 0.0)
    assert sc.any_bloch
    assert len(sc.dipoles) == 7 and all(d.magnetic for d in sc.dipoles)
    assert all(d.waveform.amp_half_complex is not None for d in sc.dipoles)
