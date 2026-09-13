# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Time-modulated media (modulation_spec).

Three main criteria:

* **Extraction**: amp = A_t·A_r, phase = φ_t + φ_r, freq = freq0, and the mask is "the Yee point
  lies inside the geometry" (point sampling; subpixel acts on the static ε only, as the notebook
  Notes say).
* **Degeneracy at zero amplitude**: a modulated scene with amplitude=0 agrees numerically with a
  scene carrying no modulation_spec (the kernel still runs every step, but the ca/cb it writes
  back differ from the static values only by float32 rounding).
* **First-order sidebands**: for an index-matched (ε_s=1) weakly modulated slab the transmitted
  sidebands match the perturbative closed form ``|E±/E0| = k±·δ·S/4``, where S is the discrete
  structure factor of the modulated column; zero fitted parameters.
"""

from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")
cp = pytest.importorskip("cupy")

from openem import flux as flux_mod
from openem import scene as scene_mod
from openem import solver
from openem.device import Kernels

C0 = 2.99792458e8
FREQ0 = C0 / 1e-6                       # λ0 = 1 µm
FM = 0.1 * FREQ0                        # modulation frequency


@pytest.fixture(scope="module")
def kernels():
    return Kernels()


def _slab_sim(medium, *, dl=0.05, lz=3.0, fwidth=FREQ0 / 10,
              mon_freqs=(FREQ0,), run_time=1.5e-12):
    """Slab scene, periodic in x/y with PML along z: a plane wave plus a transmission flux
    monitor."""
    return td.Simulation(
        size=(4 * dl, 4 * dl, lz),
        grid_spec=td.GridSpec.uniform(dl=dl),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, 0), size=(td.inf, td.inf, 1.0)),
            medium=medium)],
        sources=[td.PlaneWave(
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=fwidth),
            size=(td.inf, td.inf, 0), center=(0, 0, -lz / 2 + 0.4),
            direction="+")],
        monitors=[td.FluxMonitor(
            center=(0, 0, lz / 2 - 0.4), size=(td.inf, td.inf, 0),
            freqs=list(mon_freqs), name="flux")],
        run_time=run_time,
        boundary_spec=td.BoundarySpec.pml(x=False, y=False, z=True),
        normalize_index=None,
    )


def _modulated(eps=2.0, amp_t=1.0, ph_t=0.0, amp_r=0.4, ph_r=0.0, freq=FM):
    spec = td.ModulationSpec(permittivity=td.SpaceTimeModulation(
        time_modulation=td.ContinuousWaveTimeModulation(
            freq0=freq, amplitude=amp_t, phase=ph_t),
        space_modulation=td.SpaceModulation(amplitude=amp_r, phase=ph_r)))
    return td.Medium(permittivity=eps, modulation_spec=spec)


# ------------------------------------------------------------ extraction

def test_extraction_amp_phase_product():
    """amp is the product of the time and space parts and phase is their sum; the convention is
    copied from the tidy3d client."""
    sim = _slab_sim(_modulated(amp_t=0.5, ph_t=0.3, amp_r=0.2, ph_r=0.4))
    sc = scene_mod.from_simulation(sim)
    mo = sc.modulation
    assert mo is not None and mo.n_entry > 0
    np.testing.assert_allclose(mo.amp, 0.5 * 0.2, rtol=1e-12)
    np.testing.assert_allclose(mo.phase, 0.3 + 0.4, rtol=1e-12)
    assert mo.freq == FM
    # mask: Ex sits on the z edges and the slab [-0.5, 0.5] has grid-aligned faces -> 21 z planes
    # including both ends
    nz = sc.shape[2]
    kk = np.unique(mo.cell[mo.comp == 0] % nz)
    z = sc.grid.axes[2].edges[kk]
    assert abs(z.min() + 0.5e-6) < 1e-12 and abs(z.max() - 0.5e-6) < 1e-12


def test_unmodulated_scene_has_none():
    """An unmodulated scene has modulation None, so the solver never launches the modulation
    kernel."""
    sc = scene_mod.from_simulation(_slab_sim(td.Medium(permittivity=2.0)))
    assert sc.modulation is None and not sc.any_modulation


# ------------------------------------------------------------ fail closed

def test_spatial_data_array_phase_varies_per_cell():
    """A spatially varying modulation phase (supported since 2026-08-30).

    The travelling-wave modulation of TimeModulationTutorial has exactly this shape: a scalar
    amplitude with a phase that varies cell by cell along one axis. The criterion is not "it does
    not raise" but that **the phase really does vary with position, entry by entry**: build used
    to write ``np.full(comp.size, phase)``, and filling in a constant would not raise either.
    """
    ys = np.linspace(-0.4, 0.4, 9)
    ph = td.SpatialDataArray(
        np.linspace(-2.4, 2.4, 9).reshape(1, 9, 1),
        coords=dict(x=[0.0], y=list(ys), z=[0.0]))
    spec = td.ModulationSpec(permittivity=td.SpaceTimeModulation(
        time_modulation=td.ContinuousWaveTimeModulation(freq0=FM, amplitude=1),
        space_modulation=td.SpaceModulation(amplitude=0.4, phase=ph)))
    sc = scene_mod.from_simulation(
        _slab_sim(td.Medium(permittivity=2.0, modulation_spec=spec)))
    assert sc.modulation is not None
    phase = np.asarray(sc.modulation.phase)
    assert phase.size > 1
    assert np.ptp(phase) > 1.0, (
        f"the phase should vary with position, measured spread is only {np.ptp(phase):.3f} rad; "
        "it has probably been filled with a constant again")


def test_refuse_conductivity_modulation():
    stm = td.SpaceTimeModulation(
        time_modulation=td.ContinuousWaveTimeModulation(freq0=FM, amplitude=0.5),
        space_modulation=td.SpaceModulation(amplitude=1.0))
    med = td.Medium(permittivity=2.0, conductivity=100.0, allow_gain=True,
                    modulation_spec=td.ModulationSpec(conductivity=stm))
    with pytest.raises(NotImplementedError):
        scene_mod.from_simulation(_slab_sim(med))


def test_dispersive_modulated_scene_builds():
    """Dispersion plus time modulation: this works at the Scene level (2026-08-30).

    The modulation rides on the ε_∞ of the dispersive medium: the first item dispersion.build
    returns is ε_∞, ca/cb are recomputed every step by kmod, and the ADE recurrence coefficients
    depend only on the pole parameters and dt, so they are unaffected.
    """
    med = td.Sellmeier.from_dispersion(n=1.5, freq=FREQ0, dn_dwvl=-0.1)
    med = med.updated_copy(modulation_spec=td.ModulationSpec(
        permittivity=td.SpaceTimeModulation(
            time_modulation=td.ContinuousWaveTimeModulation(freq0=FM, amplitude=1),
            space_modulation=td.SpaceModulation(amplitude=0.1))))
    sc = scene_mod.from_simulation(_slab_sim(med))
    assert sc.any_dispersion and sc.modulation is not None


def test_dispersive_modulated_zero_amplitude_is_identity():
    """Physical criterion for the combined update: **a zero-amplitude modulation has to be the
    identity**.

    The combined formula for dispersion plus modulation (2026-08-30):

        ca = (ε_inf^n − s)/(ε_inf^{n+1} + s + G),  cb = (dt/ε₀)/(ε_inf^{n+1} + s + G)

    The numerator takes ε_inf at step n and the denominator at step n+1. In the static case the
    two are equal and it falls back to pure dispersion, (ε−s)/(ε+s+G); at s=G=0 it falls back to
    pure modulation, ε^n/ε^{n+1}.

    So adding an **amplitude=0** modulation to a dispersive scene has to leave the result exactly
    as it was. The criterion is not circular (the formula is not rewritten to check itself) and it
    runs through the whole combined path: building the tables, sampling sd/gg, and the kernel. Get
    any single term of the formula wrong and zero amplitude stops being the identity.
    """
    from openem import solver
    from openem.device import Kernels

    med = td.Sellmeier.from_dispersion(n=1.5, freq=FREQ0, dn_dwvl=-0.1)
    mod0 = med.updated_copy(modulation_spec=td.ModulationSpec(
        permittivity=td.SpaceTimeModulation(
            time_modulation=td.ContinuousWaveTimeModulation(freq0=FM, amplitude=0.0),
            space_modulation=td.SpaceModulation(amplitude=0.0))))

    k = Kernels()
    out = []
    for m in (med, mod0):
        sc = scene_mod.from_simulation(_slab_sim(m))
        res = solver.run(sc, num_steps=400, use_shutoff=False,
                         kernels=k, verbose=False)
        out.append(np.concatenate(
            [np.asarray(v.data).ravel() for v in res.phasors.values()]))
    a, b = out
    scale = max(float(np.abs(a).max()), 1e-30)
    dev = float(np.abs(a - b).max()) / scale
    assert dev < 1e-6, (
        f"a zero-amplitude modulation changed the result (relative deviation {dev:.3e}); "
        "the combined ca/cb formula is wrong")


def test_dispersive_modulated_tables_carry_s_and_g():
    """The s and G the combined formula needs really are picked up per entry; all zeros would
    mean they are not wired up."""
    from openem import setup_tables
    med = td.Sellmeier.from_dispersion(n=1.5, freq=FREQ0, dn_dwvl=-0.1)
    med = med.updated_copy(modulation_spec=td.ModulationSpec(
        permittivity=td.SpaceTimeModulation(
            time_modulation=td.ContinuousWaveTimeModulation(freq0=FM, amplitude=1),
            space_modulation=td.SpaceModulation(amplitude=0.1))))
    sc = scene_mod.from_simulation(_slab_sim(med))
    mod = setup_tables._modulation_setup(sc, sc.shape[1], sc.shape[2], sc.shape[2])
    assert mod is not None
    gg = np.asarray(mod["gg"].get() if hasattr(mod["gg"], "get") else mod["gg"])
    assert np.abs(gg).max() > 0.0, "G is all zeros in a dispersive scene; it was not read from g_dense"


def test_refuse_lossy_scene_with_modulation():
    sim = _slab_sim(_modulated())
    lossy = td.Structure(
        geometry=td.Box(center=(0, 0, 1.0), size=(td.inf, td.inf, 0.2)),
        medium=td.Medium(permittivity=2.0, conductivity=50.0))
    sim = sim.updated_copy(structures=list(sim.structures) + [lossy],
                           subpixel=False)
    with pytest.raises(NotImplementedError, match="lossy"):
        scene_mod.from_simulation(sim)


# ------------------------------------------------------------ degeneracy at zero amplitude

def test_zero_amplitude_matches_unmodulated(kernels):
    """amplitude=0: the modulation kernel still rewrites ca/cb every step, but the values differ
    from the static coefficients only by float32 rounding (a float64 division on the host then
    narrowed, vs a float32 division in the kernel). Tolerance |ΔE| ≤ 1e-5·max|E|."""
    n_steps = 1200
    sim_a = _slab_sim(_modulated(amp_r=0.0))
    sim_b = _slab_sim(td.Medium(permittivity=2.0))
    sc_a = scene_mod.from_simulation(sim_a)
    sc_b = scene_mod.from_simulation(sim_b)
    assert sc_a.any_modulation and not sc_b.any_modulation
    ra = solver.run(sc_a, num_steps=n_steps, use_shutoff=False,
                    kernels=kernels, verbose=False, return_fields=True)
    rb = solver.run(sc_b, num_steps=n_steps, use_shutoff=False,
                    kernels=kernels, verbose=False, return_fields=True)
    for name in ("Ex", "Hy"):
        a, b = ra.fields[name], rb.fields[name]
        scale = float(np.max(np.abs(b)))
        assert scale > 0
        assert float(np.max(np.abs(a - b))) <= 1e-5 * scale, name
    fa = flux_mod.plane_flux(ra.phasors["flux"].data, sc_a.grid)
    fb = flux_mod.plane_flux(rb.phasors["flux"].data, sc_b.grid)
    np.testing.assert_allclose(fa, fb, rtol=1e-4)


# ------------------------------------------------------------ first-order sidebands (analytic)

def test_first_order_sideband(kernels):
    """An index-matched weakly modulated slab: the perturbative closed form for the transmitted
    sideband amplitude, with zero fitted parameters.

    Inside the slab ε(z,t) = 1 + δ·cos(Ωt) and the background is ε=1, so the static problem is
    pure vacuum and t₀=1. First-order Born (1D Helmholtz, G = e^{ik|z−z'|}/(2ik)) gives

        |E±(ω0±Ω)| / |E0(ω0)| = (k±·δ/4)·S,
        S = |Σ_i dz_i·e^{i(k0−k±)z_i}|   (summed over the modulated Ex planes)

    The flux ratio is the square of the amplitude ratio. At δ=0.01 the O(δ²) correction is
    negligible, and the 5% tolerance covers numerical dispersion (λ/40) and DFT truncation."""
    delta = 0.01
    freqs = (FREQ0 - FM, FREQ0, FREQ0 + FM)
    fwidth = FREQ0 / 50                 # Ω/fwidth = 5 -> source leakage at the sidebands ~1e-11
    sim = _slab_sim(_modulated(eps=1.0, amp_r=delta),
                    dl=0.025, lz=4.0, fwidth=fwidth,
                    mon_freqs=freqs, run_time=10 / fwidth)
    sc = scene_mod.from_simulation(sim)
    res = solver.run(sc, use_shutoff=True, kernels=kernels, verbose=False)
    fl = flux_mod.plane_flux(res.phasors["flux"].data, sc.grid)
    ratio_lo = float(fl[0] / fl[1])
    ratio_hi = float(fl[2] / fl[1])

    # discrete structure factor: the modulated Ex planes (z edges) times the dual spacing
    mo = sc.modulation
    nz = sc.shape[2]
    kk = np.unique(mo.cell[mo.comp == 0] % nz)
    z = sc.grid.axes[2].edges[kk]
    dz = sc.grid.axes[2].dl_dual[kk]
    k0 = 2 * np.pi * FREQ0 / C0
    for ratio, f in ((ratio_lo, freqs[0]), (ratio_hi, freqs[2])):
        kpm = 2 * np.pi * f / C0
        s = abs(np.sum(dz * np.exp(1j * (k0 - kpm) * z)))
        pred = (kpm * delta / 4.0 * s) ** 2
        assert ratio == pytest.approx(pred, rel=0.05), (
            f"f={f:.3e}: measured {ratio:.4e} vs analytic {pred:.4e}")
    # the upper/lower sideband ratio is (k+/k−)²·(S+/S−)²; each one already matches its own
    # analytic value above, and this adds one more check on top
    assert ratio_hi > ratio_lo
