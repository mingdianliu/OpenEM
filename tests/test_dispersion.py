# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Dispersion ADE.

**The criterion is exact equivalence, not "about right"**: at a single frequency a dispersive
medium is physically indistinguishable from a lossy non-dispersive medium with the same complex
ε (in a linear medium every frequency propagates independently). So once the reference scene's
``(ε, σ)`` has been solved back out of the **discrete** transfer function, the phasors of the two
runs have to be bitwise identical: zero tunable parameters, the only error left is float32
rounding.

The discrete transfer function. In the z domain the two updates give

    dispersive: [(ε_∞ + χ_disc(ω)) (z−1) + s (z+1)] Ê = (Δt/ε₀) Ĉ z^{1/2}
    lossy:      [ ε_B              (z−1) + s_B(z+1)] Ê = the same right-hand side

with ``z = exp(−iωΔt)`` and ``s = σΔt/(2ε₀)``. Setting ``σ=0`` in the dispersive scene leaves two
real unknowns ``(ε_B, σ_B)`` and one complex equation, which is exactly solvable.

The poles come from the Tidy3D material library rather than being written by hand: hand-written
poles easily land on the **gain** side (Im ε > 0, which under Tidy3D's convention means
amplification), the σ_B solved back out is then negative and the reference run diverges outright.
Fits from the library are passive by construction. Silicon at 1.5 µm has ε≈12.1, very little loss,
and a non-zero real part of the residue; that last property is shared by gold, silver and silicon,
and it is the term that drops a real-valued second-order ADE to first-order accuracy.
"""

from __future__ import annotations

import os
from dataclasses import replace

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import cpml
from openem import scene as scene_mod
from openem.model import Dispersion

FREQ0 = 2.0e14                      # 1.5 µm

#: Dipole amplitude. In Tidy3D ``amplitude=1`` is a **current moment** of 1 A·µm; at dl=0.05 µm
#: the current density on the cell is ``1/dV ≈ 4e15`` and the fields climb to 1e9, still far from
#: the float32 ceiling but pointless. Only the relative difference between the two scenes matters
#: here, so the absolute amplitude can be as small as we like.
AMPLITUDE = 1e-16


def _material():
    """Five-pole fit for silicon. Returns ``(eps_inf, [(a, c), ...])``."""
    med = td.material_library["cSi"]["Green2008"]
    pr = med.pole_residue if hasattr(med, "pole_residue") else med
    return float(pr.eps_inf), [(complex(a), complex(c)) for a, c in pr.poles]


def eps_analytic(freq: float) -> complex:
    """``ε(ω) = ε_∞ − Σ [c/(jω+a) + c*/(jω+a*)]`` (Tidy3D medium.py:3214)."""
    eps_inf, poles = _material()
    w = 2 * np.pi * freq
    out = complex(eps_inf)
    for a, c in poles:
        out -= c / (1j * w + a) + np.conj(c) / (1j * w + np.conj(a))
    return out


def trap_coeffs(dt: float) -> tuple[np.ndarray, np.ndarray]:
    """``(A, B)`` per pole; trapezoid rule on ``dP/dt = qP + rE``, with ``q=a, r=c``."""
    _, poles = _material()
    a_arr, b_arr = [], []
    for q, r in poles:
        den = 1.0 - q * dt / 2.0
        a_arr.append((1.0 + q * dt / 2.0) / den)
        b_arr.append((r * dt / 2.0) / den)
    return np.array(a_arr), np.array(b_arr)


def chi_disc(freq: float, dt: float) -> complex:
    """Sum of the **discrete** susceptibility over all conjugate pairs.

    The recurrence for a conjugate pole is identically the conjugate of the original one (real E
    drive, conjugated coefficients, both starting from zero), so ``2Re(P)`` is exactly the sum of
    the two recurrences for that pole and its conjugate.
    """
    _, poles = _material()
    z = np.exp(-1j * 2 * np.pi * freq * dt)

    def one(q: complex, r: complex) -> complex:
        return r * dt * (z + 1) / 2 / (z * (1 - q * dt / 2) - (1 + q * dt / 2))

    return sum(one(q, r) + one(np.conj(q), np.conj(r)) for q, r in poles)


def equivalent_lossy(freq: float, dt: float) -> tuple[float, float]:
    """Solve for the ``(ε_B, σ_B)`` **bitwise equivalent** to the dispersive update at ``freq``."""
    eps_inf, _ = _material()
    z = np.exp(-1j * 2 * np.pi * freq * dt)
    target = (eps_inf + chi_disc(freq, dt)) * (z - 1.0)
    m = np.array([[(z - 1).real, (z + 1).real],
                  [(z - 1).imag, (z + 1).imag]])
    eps_b, s_b = np.linalg.solve(m, np.array([target.real, target.imag]))
    return float(eps_b), float(s_b * 2.0 * cpml.EPSILON_0 / dt)


def _simulation() -> "td.Simulation":
    """Small box, PML on all six sides, one point dipole, one volume FieldMonitor.

    A point dipole rather than a plane wave: the dipole's injection coefficient contains only the
    geometric volume, not the medium's η or the half-cell delay, so the excitation is **bitwise
    identical** in the two scenes and the measured difference comes purely from the medium update.
    """
    eps_inf, _ = _material()
    return td.Simulation(
        size=(1.2, 1.2, 1.2),
        grid_spec=td.GridSpec.uniform(dl=0.05),
        medium=td.Medium(permittivity=eps_inf),
        sources=[td.PointDipole(
            center=(0.0, 0.0, 0.0), polarization="Ex",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.4 * FREQ0,
                                         amplitude=AMPLITUDE))],
        monitors=[td.FieldMonitor(
            center=(0.15, 0.1, 0.05), size=(0.3, 0.2, 0.1),
            freqs=[FREQ0], name="probe")],
        run_time=2.5e-12,      # covers the 20000 steps of the longest test
        boundary_spec=td.BoundarySpec.all_sides(td.PML()),
    )


def bulk_dispersion(sc, a: np.ndarray, b: np.ndarray) -> Dispersion:
    """Spread the same set of poles over all three E components of the whole domain."""
    ncell = int(np.prod(sc.shape))
    npole = a.size
    n = 3 * ncell
    return Dispersion(
        comp=np.repeat(np.arange(3, dtype=np.int32), ncell),
        cell=np.tile(np.arange(ncell, dtype=np.int32), 3),
        pole_ofs=(np.arange(n + 1) * npole).astype(np.int32),
        am1=np.tile(a - 1.0, n),        # stores A−1, see kernels/dispersion.cu
        b=np.tile(b, n),
        g=np.full(n, 2.0 * float(np.sum(b.real))),
    )


def test_pole_mapping_matches_tidy3d():
    """The ``q=a, r=c`` mapping lines up with Tidy3D's ``eps_model``."""
    freqs = np.linspace(1.0e14, 5.0e14, 17)
    med = td.material_library["cSi"]["Green2008"]
    ref = np.array([complex(v) for v in med.eps_model(freqs)])
    ours = np.array([eps_analytic(float(f)) for f in freqs])
    np.testing.assert_allclose(ours, ref, rtol=1e-12)


def test_material_is_passive():
    """The material has to be passive.

    Tidy3D's convention is ``ε = ε' + i σ/(ω ε₀)`` (``medium.py:786``; our ``sigma_from_eps`` is
    written to match it), so **Im ε ≥ 0 is the lossy side** and a negative value means gain. Gain
    turns the equivalent σ_B negative and makes the reference run diverge; the failure would then
    have nothing to do with the ADE, and this assertion catches it first.
    """
    for f in np.linspace(1.0e14, 6.0e14, 21):
        assert eps_analytic(float(f)).imag >= 0, f"gain at {f / 1e12:.0f} THz"


def test_discrete_chi_converges_second_order():
    """The discrete ε has to converge to the analytic ε at **second order**.

    This locks the choice of scheme: switching back to the real-valued second-order form (what
    FDTDX uses today) drops the order to 1 whenever ``Re(c)≠0``, and this assertion then fails.
    Measured order 2.00.
    """
    eps_inf, _ = _material()
    freqs = np.linspace(1.0e14, 4.0e14, 9)
    ana = np.array([eps_analytic(float(f)) for f in freqs])
    errs = []
    for dt in (8e-18, 4e-18, 2e-18):
        num = np.array([eps_inf + chi_disc(float(f), dt) for f in freqs])
        errs.append(float(np.max(np.abs(num - ana) / np.abs(ana))))
    order = np.log2(errs[0] / errs[-1]) / np.log2(8e-18 / 2e-18)
    assert 1.9 < order < 2.1, f"convergence order {order:.3f}, expected 2; errors {errs}"


def test_validate_rejects_duplicate_entries():
    """A duplicate ``(comp, cell)`` has to raise: otherwise hist is overwritten and the
    physics goes silently wrong."""
    d = Dispersion(
        comp=np.array([0, 0], np.int32), cell=np.array([5, 5], np.int32),
        pole_ofs=np.array([0, 1, 2], np.int32),
        am1=np.ones(2, complex), b=np.ones(2, complex), g=np.ones(2))
    with pytest.raises(ValueError, match="duplicate"):
        d.validate((4, 4, 4))


def test_validate_rejects_out_of_range_cell():
    d = Dispersion(
        comp=np.array([0], np.int32), cell=np.array([9999], np.int32),
        pole_ofs=np.array([0, 1], np.int32),
        am1=np.ones(1, complex), b=np.ones(1, complex), g=np.ones(1))
    with pytest.raises(ValueError, match="out of range"):
        d.validate((4, 4, 4))


def test_no_dispersion_is_bit_identical():
    """``dispersion=None`` and "poles everywhere, residue 0" must be bitwise identical.

    This locks "there is only one code path": with no dispersion the new hist and G terms have to
    be exactly the identity, with no drift at the 1e-7 level.
    """
    pytest.importorskip("cupy")
    from openem import solver
    from openem.device import Kernels

    sc = scene_mod.from_simulation(_simulation())
    k = Kernels()
    zero = bulk_dispersion(sc, np.ones(1, complex), np.zeros(1, complex))
    r0 = solver.run(sc, num_steps=300, use_shutoff=False, kernels=k, verbose=False)
    r1 = solver.run(replace(sc, dispersion=zero), num_steps=300,
                    use_shutoff=False, kernels=k, verbose=False)
    np.testing.assert_array_equal(r0.field_phasors["probe"],
                                  r1.field_phasors["probe"])


def test_dispersive_equals_equivalent_lossy_medium():
    """**Decisive criterion**: the dispersive run and the equivalent lossy run agree at ``FREQ0``.

    The equivalent ``(ε_B, σ_B)`` is solved back out of the discrete transfer function, with no
    tunable parameter anywhere. The σ path has already been verified independently, so this test
    amounts to calibrating the ADE against it.

    The criterion is not a fixed threshold but that **the error drops as the DFT window gets
    longer**. The two scenes are equivalent only at ``FREQ0``; at other frequencies their
    time-domain transients differ, a finite-window DFT cannot remove that, so a short window
    necessarily leaves a residual. Measured: 4.0e-05 at 2000 steps, 2.7e-06 at 40000 steps (the
    reference medium's own decay time is about 31700 steps). If the scheme itself were wrong the
    error would be systematic and would **not** decay with the window, which is much harder to
    fake than any threshold.
    """
    pytest.importorskip("cupy")
    from openem import solver
    from openem.device import Kernels

    sc = scene_mod.from_simulation(_simulation())
    eps_inf, _ = _material()
    a, b = trap_coeffs(sc.dt)
    eps_b, sigma_b = equivalent_lossy(FREQ0, sc.dt)
    assert eps_b > 0 and sigma_b >= 0, f"unphysical equivalent medium: ε_B={eps_b}, σ_B={sigma_b}"

    ones = np.ones(sc.shape)
    ref = replace(sc, eps_ex=ones * eps_b, eps_ey=ones * eps_b, eps_ez=ones * eps_b,
                  sigma_ex=ones * sigma_b, sigma_ey=ones * sigma_b,
                  sigma_ez=ones * sigma_b)
    dis = replace(sc, eps_ex=ones * eps_inf, eps_ey=ones * eps_inf,
                  eps_ez=ones * eps_inf, dispersion=bulk_dispersion(sc, a, b))

    k = Kernels()
    errs = {}
    for n in (2000, 20000):
        r_ref = solver.run(ref, num_steps=n, use_shutoff=False, kernels=k, verbose=False)
        r_dis = solver.run(dis, num_steps=n, use_shutoff=False, kernels=k, verbose=False)
        p_ref, p_dis = r_ref.field_phasors["probe"], r_dis.field_phasors["probe"]
        scale = float(np.max(np.abs(p_ref)))
        assert np.isfinite(scale) and scale > 0, f"{n}-step reference run invalid (diverged or zero)"
        assert np.all(np.isfinite(p_dis)), f"{n}-step dispersive run diverged"
        errs[n] = float(np.max(np.abs(p_dis - p_ref))) / scale

    detail = (f"\n  ε(f0) analytic = {eps_analytic(FREQ0):.6f}"
              f"\n  ε(f0) discrete = {eps_inf + chi_disc(FREQ0, sc.dt):.6f}"
              f"\n  equivalent (ε_B, σ_B) = ({eps_b:.6f}, {sigma_b:.4e})"
              f"\n  errors {errs}")
    assert errs[2000] < 1e-4, f"short-window error {errs[2000]:.3e} is already too large{detail}"
    assert errs[20000] < errs[2000] / 3, (
        f"a 10x longer window only took the error from {errs[2000]:.3e} down to {errs[20000]:.3e}, "
        f"not like a transient residual; the scheme may be systematically wrong{detail}")


def test_dispersion_actually_changes_the_result():
    """Negative criterion: the poles have to actually do something.

    Otherwise the previous test can pass because neither side took effect, which is the worst
    kind of false pass. Silicon has ε_∞ = 1.0 while ε(1.5 µm)≈12.1, so the difference should be
    of order 1.
    """
    pytest.importorskip("cupy")
    from openem import solver
    from openem.device import Kernels

    sc = scene_mod.from_simulation(_simulation())
    eps_inf, _ = _material()
    a, b = trap_coeffs(sc.dt)
    ones = np.ones(sc.shape)
    base = replace(sc, eps_ex=ones * eps_inf, eps_ey=ones * eps_inf,
                   eps_ez=ones * eps_inf)
    k = Kernels()
    r0 = solver.run(base, num_steps=800, use_shutoff=False, kernels=k, verbose=False)
    r1 = solver.run(replace(base, dispersion=bulk_dispersion(sc, a, b)),
                    num_steps=800, use_shutoff=False, kernels=k, verbose=False)
    p0, p1 = r0.field_phasors["probe"], r1.field_phasors["probe"]
    rel = float(np.max(np.abs(p1 - p0))) / float(np.max(np.abs(p0)))
    assert rel > 0.1, f"adding poles changed the result by only {rel:.3e}; the poles did nothing"


TABLEA = os.environ.get("OPENEM_TABLEA", "/nonexistent")  # external reference archive; skipped if unset
#: The only ``subpixel=False`` dispersive case in the validation set: it has no mixed cells, so
#: this test checks pole extraction on its own with nothing from subpixel mixed in. It covers all
#: three kinds at once: a Drude model (including the DC pole at a=0), a ``Medium`` carrying a
#: conductivity, and vacuum.
DRUDE_CASE = f"{TABLEA}/DielectricMetasurfaceAbsorber/Tidy3D/simulation.json"


def test_extracted_poles_reproduce_tidy3d_epsilon(monkeypatch):
    """**Zero free parameters**: recompute ε(ω) from the ``Dispersion`` that was built and match
    ``sim.epsilon`` cell by cell.

    In the frequency domain a conjugate pair contributes ``t(A,B) + t(Ā,B̄)`` with
    ``t = B(1+z)/(z−A)``; taking ``2Re(P)`` in the time domain corresponds to this exactly
    (**taking only ``Re(t)`` is wrong**, it gives an identically zero imaginary part).

    This test locks three things at once, all of which fail silently: the pole mapping, the
    subpixel weights, and **the units of σ** (Tidy3D's ``Medium.conductivity`` is in S/µm, so
    using it directly is a factor 1e6 too small; the one time that happened it showed up as a
    perfectly correct real part with the imaginary part missing entirely). Largest relative
    difference measured: 2.9e-10.

    Poles inside the ψ-PML carry a depth-graded damping by default (boundaries.PML_POLE_DAMP,
    which **deliberately** departs from the physical ε inside the layer, by ~1e-2 deep in); it is
    switched off explicitly here. What this test locks is the fidelity of the extraction itself;
    the gating and the magnitude of the in-layer damping are locked by test_pml_pole_damp.py.
    """
    import os

    monkeypatch.setenv("OPENEM_PML_POLE_DAMP", "0")

    if not os.path.exists(DRUDE_CASE):
        pytest.skip("needs the external archive (OPENEM_TABLEA)")
    from openem.scene.media import epsilon_complex

    sim = td.Simulation.from_file(DRUDE_CASE)
    sc = scene_mod.from_file(DRUDE_CASE)
    d = sc.dispersion
    assert d is not None and d.n_entry > 0, "this case should have dispersive cells"
    ncell = int(np.prod(sc.shape))
    ofs = d.pole_ofs.astype(np.int64)
    A, B = d.am1 + 1.0, d.b
    freqs = np.asarray(sim.monitors[0].freqs, dtype=np.float64)

    worst = 0.0
    for c, key in enumerate(("Ex", "Ey", "Ez")):
        eps_inf = (sc.eps_ex, sc.eps_ey, sc.eps_ez)[c].ravel()
        sig = (sc.sigma_ex, sc.sigma_ey, sc.sigma_ez)[c]
        sig = np.zeros(ncell) if sig is None else sig.ravel()
        m = np.flatnonzero(d.comp == c)
        lo, hi = ofs[m], ofs[m + 1]
        assert np.all(hi > lo), "an entry has zero poles"
        for f in freqs[[0, len(freqs) // 2, -1]]:
            w = 2 * np.pi * float(f)
            z = np.exp(-1j * w * sc.dt)
            t = B * (1 + z) / (z - A) + np.conj(B) * (1 + z) / (z - np.conj(A))
            # reduceat only over this component's own pole range: reduceat's last segment always
            # runs to the end of the array, so without the slice it would also accumulate the
            # poles of the components that follow
            seg = np.add.reduceat(t[int(lo[0]):int(hi[-1])], lo - lo[0])
            chi = np.zeros(ncell, dtype=np.complex128)
            chi[d.cell[m]] = seg
            got = eps_inf + chi + 1j * sig / (w * cpml.EPSILON_0)
            ref = np.asarray(epsilon_complex(sim, key, float(f))).ravel()
            rel = np.abs(got - ref) / np.maximum(np.abs(ref), 1.0)
            worst = max(worst, float(rel.max()))
    assert worst < 1e-8, f"largest per-cell relative difference in the ε reconstruction: {worst:.3e}"


# ---------------------------------------------------------------- D form

def zeta_disc(freq: float, dt: float, zinf: float, am1, b) -> complex:
    """The **discrete** ζ of the D form: ``ζ_∞ + Σ B_m(1+z)/(z−A_m)``, summed over all roots."""
    z = np.exp(-1j * 2 * np.pi * freq * dt)
    A = am1 + 1.0
    return zinf + complex(np.sum(b * (1 + z) / (z - A)))


def equivalent_lossy_d(freq: float, dt: float, zinf: float, am1, b):
    """Solve for the ``(ε_B, σ_B)`` equivalent to the D form at ``freq``.

    In the z domain the D form is ``E(z−1)/ζ_disc = (Δt/ε₀)Ĉz^{1/2}`` and the lossy form is
    ``[ε_B(z−1)+s_B(z+1)]E = the same right-hand side``, so
    ``ε_B(z−1) + s_B(z+1) = (z−1)/ζ_disc``.
    """
    z = np.exp(-1j * 2 * np.pi * freq * dt)
    target = (z - 1.0) / zeta_disc(freq, dt, zinf, am1, b)
    m = np.array([[(z - 1).real, (z + 1).real],
                  [(z - 1).imag, (z + 1).imag]])
    eps_b, s_b = np.linalg.solve(m, np.array([target.real, target.imag]))
    return float(eps_b), float(s_b * 2.0 * cpml.EPSILON_0 / dt)


def bulk_dispersion_d(sc, zinf: float, am1: np.ndarray, b: np.ndarray):
    """**Pure harmonic branch**: the same set of ζ poles everywhere, ``β = 1``, arithmetic branch
    empty.

    This used to be a separate ``DispersionD``, but at β=1 the mixed form is term for term the
    same thing (``k1=0`` → denominator 1 → ``d_h = d_tot``, ``E = k2·d_tot + SQ``), and keeping a
    second code path only means more code to maintain and to verify. So we construct that limit
    instead.
    """
    from openem.model import DispersionMix

    ncell = int(np.prod(sc.shape))
    n = 3 * ncell
    z = np.zeros(0, dtype=np.complex128)
    return DispersionMix(
        comp=np.repeat(np.arange(3, dtype=np.int32), ncell),
        cell=np.tile(np.arange(ncell, dtype=np.int32), 3),
        p_ofs=np.zeros(n + 1, dtype=np.int32), pa=z, pb=z,
        q_ofs=(np.arange(n + 1) * am1.size).astype(np.int32),
        qa=np.tile(am1, n), qb=np.tile(b, n),
        beta=np.ones(n), eps_inf=np.ones(n),
        zeta_inf=np.full(n, zinf, dtype=np.float64),
    )


def _zeta_tables(dt: float):
    from openem.scene import dispersion as dsp

    eps_inf, poles = _material()
    zinf, p, r = dsp.zeta_poles(eps_inf, poles)
    am1, b = dsp.trap_from_poles(p, r, dt)
    return zinf, am1, b


def test_zeta_poles_reproduce_inverse_epsilon():
    """The pole expansion of ``1/ε`` has to reproduce ``1/eps_model``, with every pole in the
    left half-plane.

    The variable has to be ``s = −iω``, the one the already verified ``q=a, r=c`` mapping lives
    in. Using ``+iω`` moves the poles into the right half-plane and makes a stable scheme look
    unstable; we walked into exactly that.
    """
    from openem.scene import dispersion as dsp

    eps_inf, poles = _material()
    zinf, p, r = dsp.zeta_poles(eps_inf, poles)
    assert np.all(p.real < 0), f"pole in the right half-plane: max Re = {p.real.max():.3e}"
    freqs = np.linspace(1.0e14, 5.0e14, 17)
    s = -1j * 2 * np.pi * freqs
    got = zinf + sum(r[m] / (s - p[m]) for m in range(p.size))
    med = td.material_library["cSi"]["Green2008"]
    ref = 1.0 / np.array([complex(v) for v in med.eps_model(freqs)])
    np.testing.assert_allclose(got, ref, rtol=1e-9)


def test_zeta_rejects_right_half_plane_poles():
    """A pole in the right half-plane has to raise: the recurrence would diverge, so it must not
    just run quietly."""
    from openem.scene import dispersion as dsp

    # a negative ε_∞ puts a root of the numerator of ε(s) in the right half-plane
    with pytest.raises(ValueError, match="right half plane"):
        dsp.zeta_poles(-2.0, _material()[1])


# MgF2_Horiba from MIMResonator: lossless Sellmeier-type single pole (a and c both imaginary).
_MGF2_EPS_INF = 1.0
_MGF2_POLES = [(-2.5358092974503356e16j, 1.1398462792039258e16j)]


def test_zeta_accepts_lossless_on_axis_poles():
    """For a lossless material the ζ poles sit **exactly** on the imaginary axis (the longitudinal
    zeros of ε, ±iω_L), and building the table must not reject them.

    This once shut MIMResonator out entirely, even though it was measured to have zero mixed
    cells and never uses the ζ table at all. The criterion: after being projected back onto the
    imaginary axis the poles still have to reproduce 1/eps_model (the analytic zero is
    ω_L² = ω₀² + 2c₀ω₀, zero free parameters).
    """
    from openem.scene import dispersion as dsp

    zinf, p, r = dsp.zeta_poles(_MGF2_EPS_INF, _MGF2_POLES)
    assert np.all(p.real == 0.0), "ζ poles of a lossless material should land on the imaginary axis"
    w0, c0 = 2.5358092974503356e16, 1.1398462792039258e16
    wl = np.sqrt(w0 ** 2 + 2 * c0 * w0)
    np.testing.assert_allclose(np.sort(np.abs(p.imag)), [wl, wl], rtol=1e-12)
    freqs = np.linspace(2.0e14, 8.0e14, 13)
    s = -1j * 2 * np.pi * freqs
    got = zinf + sum(r[m] / (s - p[m]) for m in range(p.size))
    ref = 1.0 / (_MGF2_EPS_INF + 2 * c0 * w0 / (s ** 2 + w0 ** 2))
    np.testing.assert_allclose(got, ref, rtol=1e-9)


def test_mix_cell_accepts_marginal_zeta_pole_with_bounded_coeffs():
    """A marginally stable ζ pole landing on a slanted-interface cell is **accepted**, and its
    coefficients stay bounded in fp32.

    This used to fail closed, on the grounds that long-run boundedness in fp32 had never been
    verified. It was verified on 2026-08-31 (tests/test_zeta_stability.py: free oscillation for
    50000 steps, written the way the kernel writes it, amplitude ratio within ±3e-3), so it is
    now let through. This test locks two things: **no exception any more**, and the recurrence
    coefficients produced satisfy ``|A| ≤ 1 + tolerance``; a coefficient with |A| clearly above
    1 really would make a long run diverge.
    """
    from openem.scene import dispersion as dsp

    dt = 1e-17
    kept = [[], _MGF2_POLES]
    eps_inf_of = np.array([1.0, _MGF2_EPS_INF])
    sigma_of = np.array([0.0, 0.0])
    zeta_of = [(1.0, np.zeros(0, complex), np.zeros(0, complex)),
               dsp.zeta_poles(_MGF2_EPS_INF, _MGF2_POLES)]
    mixp = {"cells": np.array([0]), "f": np.array([0.5]),
            "beta": np.array([0.5]), "mi": np.array([1]), "mj": np.array([0])}
    out = dsp._mix_entries(0, mixp, kept, sigma_of, eps_inf_of, zeta_of, dt)

    qa = np.asarray(out["qa"])
    assert qa.size, "the harmonic branch should produce recurrence coefficients for the ζ poles"
    # what is stored is A−1; the kernel uses complex64
    mag = np.abs(1.0 + qa.astype(np.complex64).astype(np.complex128))
    assert float(mag.max()) <= 1.0 + 1e-6, (
        f"max |A| is {mag.max():.9f}, above 1; a long run would diverge")


def test_harmonic_branch_equals_equivalent_lossy_medium():
    """**Decisive criterion for the pure harmonic branch**: it agrees with the equivalent lossy
    run, and the error drops as the DFT window grows.

    Structurally identical to the E-form test
    (:func:`test_dispersive_equals_equivalent_lossy_medium`), except that the equivalent medium
    comes from inverting the D form's own discrete transfer function.
    """
    pytest.importorskip("cupy")
    from openem import solver
    from openem.device import Kernels

    sc = scene_mod.from_simulation(_simulation())
    zinf, am1, b = _zeta_tables(sc.dt)
    eps_b, sigma_b = equivalent_lossy_d(FREQ0, sc.dt, zinf, am1, b)
    assert eps_b > 0 and sigma_b >= 0, f"unphysical equivalent medium: ε_B={eps_b}, σ_B={sigma_b}"

    ones = np.ones(sc.shape)
    ref = replace(sc, eps_ex=ones * eps_b, eps_ey=ones * eps_b, eps_ez=ones * eps_b,
                  sigma_ex=ones * sigma_b, sigma_ey=ones * sigma_b,
                  sigma_ez=ones * sigma_b)
    dis = replace(sc, dispersion_mix=bulk_dispersion_d(sc, zinf, am1, b))

    k = Kernels()
    errs = {}
    for n in (2000, 20000):
        r_ref = solver.run(ref, num_steps=n, use_shutoff=False, kernels=k, verbose=False)
        r_dis = solver.run(dis, num_steps=n, use_shutoff=False, kernels=k, verbose=False)
        p_ref, p_dis = r_ref.field_phasors["probe"], r_dis.field_phasors["probe"]
        scale = float(np.max(np.abs(p_ref)))
        assert np.isfinite(scale) and scale > 0, f"{n}-step reference run invalid"
        assert np.all(np.isfinite(p_dis)), f"{n}-step D-form run diverged"
        errs[n] = float(np.max(np.abs(p_dis - p_ref))) / scale

    detail = (f"\n  ζ(f0) discrete = {zeta_disc(FREQ0, sc.dt, zinf, am1, b):.6e}"
              f"\n  1/ε(f0) analytic = {1 / eps_analytic(FREQ0):.6e}"
              f"\n  equivalent (ε_B, σ_B) = ({eps_b:.6f}, {sigma_b:.4e})"
              f"\n  errors {errs}")
    assert errs[2000] < 1e-3, f"short-window error {errs[2000]:.3e} is already too large{detail}"
    assert errs[20000] < errs[2000] / 3, (
        f"a 10x longer window only took the error from {errs[2000]:.3e} down to {errs[20000]:.3e}, "
        f"not like a transient residual{detail}")


def test_harmonic_and_arithmetic_branches_agree_on_same_material():
    """The same material down both code paths has to give the same answer, to second order.

    The E form discretizes the poles of ε, the D form those of 1/ε: two different
    discretizations, so they should agree only to O(Δt²) and must not be bitwise identical. What
    this criterion locks is that **they describe the same material**: the difference has to be
    far smaller than the difference switching to another material makes.
    """
    pytest.importorskip("cupy")
    from openem import solver
    from openem.device import Kernels

    sc = scene_mod.from_simulation(_simulation())
    eps_inf, _ = _material()
    a, b_e = trap_coeffs(sc.dt)
    zinf, am1, b_d = _zeta_tables(sc.dt)
    ones = np.ones(sc.shape)
    e_form = replace(sc, eps_ex=ones * eps_inf, eps_ey=ones * eps_inf,
                     eps_ez=ones * eps_inf,
                     dispersion=bulk_dispersion(sc, a, b_e))
    d_form = replace(sc, dispersion_mix=bulk_dispersion_d(sc, zinf, am1, b_d))
    # reference difference: treat ε_∞ as the whole material, i.e. "another material"
    other = replace(sc, eps_ex=ones * eps_inf, eps_ey=ones * eps_inf,
                    eps_ez=ones * eps_inf)

    k = Kernels()
    ph = {}
    for tag, s_ in (("e", e_form), ("d", d_form), ("other", other)):
        r = solver.run(s_, num_steps=4000, use_shutoff=False, kernels=k, verbose=False)
        ph[tag] = r.field_phasors["probe"]
        assert np.all(np.isfinite(ph[tag])), f"{tag} diverged"
    scale = float(np.max(np.abs(ph["e"])))
    same = float(np.max(np.abs(ph["d"] - ph["e"]))) / scale
    diff = float(np.max(np.abs(ph["other"] - ph["e"]))) / scale
    assert same < 0.02, f"the two paths differ by {same:.3e}, not like the same material"
    assert same < diff / 20, (
        f"the two paths differ by {same:.3e}, while another material differs by only "
        f"{diff:.3e}; the criterion has no discriminating power")


# ---------------------------------------------------------------- mixed form

def bulk_mix(sc, beta: float, pa, pb, qa, qb, eps_inf: float, zeta_inf: float):
    """The same coefficients over all three components of the whole domain, with constant β."""
    from openem.model import DispersionMix

    ncell = int(np.prod(sc.shape))
    n = 3 * ncell
    return DispersionMix(
        comp=np.repeat(np.arange(3, dtype=np.int32), ncell),
        cell=np.tile(np.arange(ncell, dtype=np.int32), 3),
        p_ofs=(np.arange(n + 1) * pa.size).astype(np.int32),
        pa=np.tile(pa, n), pb=np.tile(pb, n),
        q_ofs=(np.arange(n + 1) * qa.size).astype(np.int32),
        qa=np.tile(qa, n), qb=np.tile(qb, n),
        beta=np.full(n, float(beta)),
        eps_inf=np.full(n, float(eps_inf)),
        zeta_inf=np.full(n, float(zeta_inf)),
    )


#: The second material in the bulk tests (non-dispersive). **There have to be two materials**:
#: with only one, the arithmetic and the harmonic mean both equal the material itself, β does
#: nothing at all, and the discriminating-power criterion becomes vacuous (the first version of
#: this test fell into exactly that).
EPS2 = 2.25
FILL = 0.5


def _mix_pieces(sc):
    """The six things the mixed form needs when cSi and ε=2.25 are mixed at f=0.5.

    Arithmetic branch: ``ε_∞ = f·ε_∞1 + (1−f)·ε2``, residues times f (the non-dispersive material
    has no poles). Harmonic branch: ``ζ_∞ = f/ε_∞1 + (1−f)/ε2``, the ζ residues likewise times f.
    Measured at 1.5 µm the two branches give 7.19 and 3.80, which is far enough apart.
    """
    from openem.scene import dispersion as dsp

    eps_inf1, poles = _material()
    f = FILL
    pa_full, pb_full = trap_coeffs(sc.dt)
    pa = pa_full - 1.0                              # stores A−1
    pb = pb_full * f                                # residues scaled by the fill fraction
    eps_inf = f * eps_inf1 + (1 - f) * EPS2

    z1, p1, r1 = dsp.zeta_poles(eps_inf1, poles)
    qa, qb = dsp.trap_from_poles(p1, r1 * f, sc.dt)
    zeta_inf = f * z1 + (1 - f) / EPS2
    return eps_inf, zeta_inf, pa, pb, qa, qb


def test_mix_beta_zero_equals_e_form():
    """**β=0 has to reduce exactly to the E form.**

    The two are algebraically identical: the mixed form accumulates ``d_tot`` and then solves
    ``E=(d_tot−SP)/(ε_∞+G_P)``, while ``d_tot^{n+1} = ε_∞E^{n+1} + Σ2Re(P^{n+1})`` follows by
    induction (both start at 0). So the only difference should be float32 rounding. The E-form
    path has already been verified independently, so this test hangs the mixed form off it.
    """
    pytest.importorskip("cupy")
    from openem import solver
    from openem.device import Kernels

    sc = scene_mod.from_simulation(_simulation())
    eps_inf, zinf, pa, pb, qa, qb = _mix_pieces(sc)
    ones = np.ones(sc.shape)
    e_form = replace(sc, eps_ex=ones * eps_inf, eps_ey=ones * eps_inf,
                     eps_ez=ones * eps_inf,
                     dispersion=bulk_dispersion(sc, pa + 1.0, pb))
    mixed = replace(sc, dispersion_mix=bulk_mix(sc, 0.0, pa, pb, qa, qb,
                                                eps_inf, zinf))
    k = Kernels()
    p1 = solver.run(e_form, num_steps=3000, use_shutoff=False, kernels=k,
                    verbose=False).field_phasors["probe"]
    p2 = solver.run(mixed, num_steps=3000, use_shutoff=False, kernels=k,
                    verbose=False).field_phasors["probe"]
    assert np.all(np.isfinite(p2)), "the mixed form at β=0 diverged"
    err = float(np.max(np.abs(p2 - p1))) / float(np.max(np.abs(p1)))
    assert err < 1e-4, f"β=0 differs from the E form by {err:.3e}, only float32 rounding should remain"


def test_mix_middle_beta_lies_between_and_is_stable():
    """An intermediate β has to land between the two limits, and must not diverge.

    This is not an accuracy criterion (there is no independent reference for an intermediate β),
    it is a **discriminating-power** criterion: if the mixed form were not sensitive to β at all,
    the two tests above would still each pass while the whole thing was wrong.
    """
    pytest.importorskip("cupy")
    from openem import solver
    from openem.device import Kernels

    sc = scene_mod.from_simulation(_simulation())
    eps_inf, zinf, pa, pb, qa, qb = _mix_pieces(sc)
    k = Kernels()
    out = {}
    for beta in (0.0, 0.35, 0.7, 1.0):
        sc2 = replace(sc, dispersion_mix=bulk_mix(sc, beta, pa, pb, qa, qb,
                                                  eps_inf, zinf))
        r = solver.run(sc2, num_steps=3000, use_shutoff=False, kernels=k, verbose=False)
        p = r.field_phasors["probe"]
        assert np.all(np.isfinite(p)), f"β={beta} diverged"
        out[beta] = p
    span = float(np.max(np.abs(out[1.0] - out[0.0])))
    scale = float(np.max(np.abs(out[0.0])))
    assert span / scale > 0.05, (
        f"going from β=0 to β=1 changed the result by only {span / scale:.3e}; β has no effect "
        f"and the criterion has no discriminating power")
    # monotone crossing: each middle setting should be closer to either end than the two ends
    # are to each other
    for beta in (0.35, 0.7):
        d0 = float(np.max(np.abs(out[beta] - out[0.0])))
        d1 = float(np.max(np.abs(out[beta] - out[1.0])))
        assert d0 < span * 1.05 and d1 < span * 1.05, (
            f"β={beta} falls outside the interval between β=0 and β=1 "
            f"({d0:.3e}, {d1:.3e} vs {span:.3e})")


def test_mix_validate_rejects_beta_out_of_range():
    """β outside [0,1] has to raise: it is the square of the normal projection, so out of range
    means the fit is not trustworthy."""
    from openem.model import DispersionMix

    z = np.zeros(0, complex)
    mm = DispersionMix(
        comp=np.array([0], np.int32), cell=np.array([3], np.int32),
        p_ofs=np.array([0, 0], np.int32), pa=z, pb=z,
        q_ofs=np.array([0, 0], np.int32), qa=z, qb=z,
        beta=np.array([1.4]), eps_inf=np.array([2.0]), zeta_inf=np.array([0.5]))
    with pytest.raises(ValueError, match=r"\[0,1\]"):
        mm.validate((4, 4, 4))


def test_material_weights_are_nonnegative():
    """Lock "material weights are non-negative": a pole with a negative weight is a reversed
    oscillator, that is, a narrow-band gain medium.

    Unconstrained least squares leaves noise weights of ±1e-6..1e-4 on materials that are not
    present, and the negative half of those has gain in the narrow k band where the pole crosses
    the grid dispersion branch: the SiO2 background cells of RingResonator carried cSi poles at
    the −3e-6 level and overflowed to NaN at step 159,200. ``weights.decompose`` clamps the noisy
    negative weights to zero, and what is locked here is the product: within any pole family, the
    projection of every b onto the family's reference direction has to be non-negative.
    """
    sio2 = td.PoleResidue(eps_inf=1.5385442336875639, poles=[
        (-11504139.374277674 - 1.595196740783775e16j,
         7507685.43042605 + 4535416182817100.0j),
        (-249390.0121716501 - 172284807034385.44j,
         46272.59577676453 + 99705881297053.99j)])
    csi = td.PoleResidue(eps_inf=1.0, poles=[
        (-1.7473849958109988 - 6409829457220535.0j,
         0.06947645444424029 + 3.4268436708700284e16j)])
    f0 = 1.9e14
    sim = td.Simulation(
        size=(0.6, 0.6, 0.6), grid_spec=td.GridSpec.uniform(dl=0.025),
        medium=sio2,
        structures=[td.Structure(
            geometry=td.Sphere(center=(0.0, 0.0, 0.0), radius=0.17), medium=csi)],
        sources=[td.PointDipole(
            center=(0.0, 0.0, 0.21), polarization="Ex",
            source_time=td.GaussianPulse(freq0=f0, fwidth=0.2 * f0))],
        monitors=[td.FieldMonitor(center=(0.0, 0.0, 0.0), size=(0.1, 0.1, 0.1),
                                  freqs=[0.9 * f0, f0, 1.1 * f0], name="p")],
        run_time=1e-13,
        boundary_spec=td.BoundarySpec.all_sides(td.PML()))
    sc = scene_mod.from_simulation(sim)
    d = sc.dispersion
    assert d is not None and d.n_entry > 0
    fam = np.stack([np.round(d.am1.real, 15), np.round(d.am1.imag, 15)], axis=1)
    checked = 0
    for f in np.unique(fam, axis=0):
        m = (fam == f).all(axis=1)
        bm = d.b[m]
        ref = bm[np.abs(bm).argmax()]
        w = (bm.real * ref.real + bm.imag * ref.imag) / (abs(ref) ** 2)
        assert w.min() >= 0.0, (
            f"pole family am1={f} has a negative projected weight {w.min():.3e} "
            f"({int((w < 0).sum())} entries)")
        checked += 1
    assert checked >= 2, f"expected at least the two SiO2/cSi pole families, only saw {checked}"


def test_mix_sigma_beta_zero_matches_semi_implicit():
    """**The s=0 integrator for σ (mix path) ≡ semi-implicit σ (ordinary path).**

    Discretizing D = ε_B·E + ∫(σ/ε₀)E dt with the trapezoid rule is algebraically identical to
    semi-implicit σ, so the only difference should be float32 rounding. This pins down the
    residue convention of the step in dispersion._mix_entries that turns σ back into an s=0 pole
    (c = σ/(2ε₀), with the kernel taking 2Re when it recombines).
    """
    pytest.importorskip("cupy")
    from openem import cpml, solver
    from openem.device import Kernels
    from openem.scene import dispersion as dsp

    sc = scene_mod.from_simulation(_simulation())
    eps_b, sigma_b = 4.0, 40.0
    ones = np.ones(sc.shape)
    ref = replace(sc, eps_ex=ones * eps_b, eps_ey=ones * eps_b, eps_ez=ones * eps_b,
                  sigma_ex=ones * sigma_b, sigma_ey=ones * sigma_b,
                  sigma_ez=ones * sigma_b)

    a_dc, b_dc = trap_coeffs_sigma(sigma_b, sc.dt)
    z0, pp, rr = dsp.zeta_poles(eps_b, [], sigma_b)
    qa, qb = dsp.trap_from_poles(pp, rr, sc.dt)
    mixed = replace(sc, eps_ex=ones * eps_b, eps_ey=ones * eps_b, eps_ez=ones * eps_b,
                    dispersion_mix=bulk_mix(sc, 0.0, a_dc, b_dc, qa, qb,
                                            eps_b, z0))
    k = Kernels()
    p1 = solver.run(ref, num_steps=3000, use_shutoff=False, kernels=k,
                    verbose=False).field_phasors["probe"]
    p2 = solver.run(mixed, num_steps=3000, use_shutoff=False, kernels=k,
                    verbose=False).field_phasors["probe"]
    assert np.all(np.isfinite(p2)), "the σ integrator diverged"
    err = float(np.max(np.abs(p2 - p1))) / float(np.max(np.abs(p1)))
    assert err < 1e-4, f"integrator σ and semi-implicit σ differ by {err:.3e}, expected float32 rounding"


def trap_coeffs_sigma(sigma: float, dt: float):
    """(A−1, B) for the s=0 integrator pole of σ: the same formula as dispersion._mix_entries.

    Note that this file has its own single-argument material helper also called ``trap_coeffs``;
    the one wanted here is ``scene.poles.trap_coeffs``.
    """
    from openem import cpml
    from openem.scene.poles import trap_coeffs as _tc
    a, b = _tc([(0.0, sigma / (2.0 * cpml.EPSILON_0))], dt)
    return a, b
