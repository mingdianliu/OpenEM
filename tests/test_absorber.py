# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Adiabatic absorber (``td.Absorber``): folding σ in, reflectance calibration, and the stability
of "dispersion + absorber".

How it works (corrected 2026-09-07): an Absorber is a graded lossy layer with **conductivity only,
no magnetic loss**, backed by PEC. σ uses the same normalization as the PML and is folded into the
E update coefficients (``scene/boundaries.fold_absorber_sigma``), each cell with its own ε. There
is no separate damping kernel.

Three things are locked here:

* σ really is folded into the layer band (nonzero inside, zero outside, monotonic outwards) and
  there is no AbsorberSlab;
* normal-incidence reflectance in vacuum lands near a **1D transfer-matrix model built from the
  same σ definition**: too low means it has turned back into a "matched" absorber (about 10x
  stronger than the reference implementation), too high means σ was lost or the wave hit the PEC;
* "absorber + dispersion" stays bounded and converges (the old ψ recursion diverged exponentially
  here).
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")
pytest.importorskip("cupy")

from openem import cpml
from openem import flux as flux_mod
from openem import scene as scene_mod
from openem import setup_tables, solver
from openem.device import Kernels

F0 = 5.0e14
DL = 0.03          # λ/20: σ ∝ 1/dl, at λ/60 the layer band becomes a wall


def _absorber_z(nl: int = 40, sigma_max: float = 6.4) -> "td.Boundary":
    par = td.AbsorberParams(sigma_order=3, sigma_min=0.0, sigma_max=sigma_max)
    return td.Boundary(plus=td.Absorber(num_layers=nl, parameters=par),
                       minus=td.Absorber(num_layers=nl, parameters=par))


def _vacuum_sim(nl: int) -> "td.Simulation":
    return td.Simulation(
        size=(4 * DL, 4 * DL, 2.0), grid_spec=td.GridSpec.uniform(dl=DL),
        sources=[td.PlaneWave(center=(0, 0, 0.0), size=(td.inf, td.inf, 0),
            direction="-", source_time=td.GaussianPulse(freq0=F0, fwidth=0.2 * F0))],
        monitors=[td.FluxMonitor(center=(0, 0, 0.4), size=(td.inf, td.inf, 0),
                                 freqs=[F0], name="back")],
        run_time=6e-13,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.periodic(), y=td.Boundary.periodic(), z=_absorber_z(nl)))


def _tmm_reflection_db(nl: int, dl_m: float, order: int = 3,
                       sigma_max: float = 6.4, sub: int = 8) -> float:
    """1D transfer matrix: an E-only graded-conductivity layer in vacuum backed by PEC, normal
    incidence reflectance (dB).

    The σ definition is the same as ``cpml.profile`` (avg_speed=1). Both of the reference
    implementation's AbsorbingBoundaryReflection curves, over layer count and over resolution,
    were matched with this model.
    """
    w = 2 * np.pi * F0
    k0 = w / cpml.C_0
    d = nl * dl_m
    m = nl * sub
    s = (np.arange(m) + 0.5) / m
    sig = sigma_max / (cpml.ETA_0 * dl_m) * s ** order
    eps = 1.0 + 1j * sig / (w * cpml.EPSILON_0)
    nn = np.sqrt(eps)
    eta = 1.0 / nn
    z = 0.0 + 0j                      # PEC
    h = d / m
    for j in range(m - 1, -1, -1):
        t = np.tan(k0 * nn[j] * h)
        z = eta[j] * (z - 1j * eta[j] * t) / (eta[j] - 1j * z * t)
    return float(10 * np.log10(abs((z - 1.0) / (z + 1.0)) ** 2))


def _reflection_db(sc, kernels) -> float:
    src = sc.sources[0]
    dn = int(np.clip(src.plane_index + src.direction * 12, 10, sc.shape[2] - 11))
    probe = replace(sc.flux_monitors[0], name="_dn", plane_index=dn, normal_dir=1)
    res = solver.run(replace(sc, flux_monitors=list(sc.flux_monitors) + [probe]),
                     use_shutoff=False, kernels=kernels, verbose=False)
    pin = float(np.abs(flux_mod.plane_flux(res.phasors["_dn"].data, sc.grid))[0])
    back = float(np.abs(flux_mod.plane_flux(res.phasors["back"].data, sc.grid))[0])
    return float(10 * np.log10(back / pin))


def test_absorber_sigma_folded_into_layers():
    """σ lives only inside the layer band, grows monotonically outwards, and is present on all
    three components; no AbsorberSlab and no ψ coefficients."""
    sc = scene_mod.from_simulation(_vacuum_sim(40))
    assert not sc.absorbers and not sc.pml
    nz = sc.shape[2]
    for arr in (sc.sigma_ex, sc.sigma_ey, sc.sigma_ez):
        assert arr is not None
        col = arr[0, 0, :]
        assert np.all(col[40:nz - 40] == 0.0), "outside the layers there must be no loss"
        lo, hi = col[:40], col[nz - 40:]
        assert lo[0] > 0 and np.all(np.diff(lo) <= 0), "low end rises outwards (index decreasing)"
        assert hi[-1] > 0 and np.all(np.diff(hi) >= 0), "high end rises outwards (index increasing)"
    # Magnitude: outermost whole-cell tangential component = σ_max/(η₀·dl) (vacuum avg_speed=1,
    # step=1; dl is the grid's actual step, tidy3d nudges a uniform grid to fit the domain length)
    expect = 6.4 / (cpml.ETA_0 * float(sc.grid.axes[2].dl[0]))
    assert abs(sc.sigma_ex[0, 0, 0] / expect - 1.0) < 1e-9


def test_absorber_reflection_matches_1d_model():
    """Normal incidence in vacuum: reflectance tracks the 1D model (40 layers ≈ -34 dB), and
    80 layers is another >= 20 dB lower.

    A matched absorber (the old implementation) lands around -300 dB here, which is exactly the
    "10x stronger than the reference implementation" failure mode, so the lower bound is locked
    as well.
    """
    k = Kernels()
    got = {}
    for nl in (40, 80):
        got[nl] = _reflection_db(scene_mod.from_simulation(_vacuum_sim(nl)), k)
    model40 = _tmm_reflection_db(40, DL * 1e-6)
    assert abs(got[40] - model40) < 6.0, (
        f"40 layers reflect {got[40]:.1f} dB, 1D model says {model40:.1f} dB")
    assert got[80] < got[40] - 20.0, (
        f"80 layers {got[80]:.1f} dB should be clearly below 40 layers {got[40]:.1f} dB")
    assert got[80] > -120.0, "absurdly low: back to a matched absorber?"


def test_absorber_with_dispersion_is_stable():
    """An absorber plus a dispersive medium must stay bounded and converge; the old ψ route
    diverged exponentially here.

    The medium is the library's 5-pole cSi (including poles with Re(c) != 0, exactly the kind that
    triggered the mismatch at the time), filling the whole domain, with absorbers at both z ends.
    Verdict: double the step count and the maximum of the accumulated phasor must grow by less than
    1.5x (divergence gives 10x per 8000 steps; on convergence the DFT sum tends to a constant).
    """
    med = td.material_library["cSi"]["Green2008"]
    f0 = 2.0e14
    sim = td.Simulation(
        size=(0.4, 0.4, 1.6), grid_spec=td.GridSpec.uniform(dl=0.025),
        medium=med,
        sources=[td.PointDipole(center=(0, 0, 0), polarization="Ex",
            source_time=td.GaussianPulse(freq0=f0, fwidth=0.2 * f0,
                                         amplitude=1e-16))],
        monitors=[td.FieldMonitor(center=(0, 0, 0.3), size=(0.2, 0.2, 0),
                                  freqs=[f0], name="p")],
        run_time=2e-12,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.periodic(), y=td.Boundary.periodic(), z=_absorber_z(20)))
    sc = scene_mod.from_simulation(sim)
    assert not sc.absorbers and sc.dispersion is not None and sc.dispersion.n_entry > 0
    assert sc.sigma_ex is not None and (sc.sigma_ex[0, 0, :20] > 0).all()
    assert setup_tables._absorber_disp_decay(sc) is None, "no slab means no P damping path"

    k = Kernels()
    v = {}
    for ns in (4000, 8000):
        r = solver.run(sc, num_steps=ns, use_shutoff=False, kernels=k, verbose=False)
        v[ns] = float(np.abs(r.field_phasors["p"]).max())
    growth = v[8000] / max(v[4000], 1e-300)
    assert growth < 1.5, (
        f"steps 4000 -> 8000, phasor maximum {v[4000]:.3e} -> {v[8000]:.3e}"
        f" ({growth:.2f}x): diverging again")
