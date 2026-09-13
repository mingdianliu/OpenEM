# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Magnetic dipoles (polarization=Hx/Hy/Hz): the H dual of inject_dipole.

The criterion has no free parameters: radiated power in vacuum against the analytic formula. The
electric-dipole chain was already checked against the analytic far field to 0.2% (with IL =
1 A·µm), and by duality J<->M, ε₀<->μ₀ the magnetic side is

    P_m(ω) = ε₀ ω² |K|² / (12π c),   K = 1 V·µm = 1e-6 V·m

(the electric side is P_e = μ₀ ω² |S|²/(12π c); dividing the two gives P_m/P_e = ε₀/μ₀, which is
independent of geometry, so that ratio doubles as a duality self-check). The amplitude convention
"amplitude=1 = a magnetic current moment of 1 V·µm" follows from duality, and was confirmed end to
end by comparing Bandstructure sim_0's FieldTimeMonitor against the reference output.

Note: power is |·|², so it is insensitive to overall phase; this **cannot** tell whether the
amplitude is taken at the whole step or the half step. That half is covered by the step-by-step
time-domain comparison of sim_0, which is phase sensitive.
"""

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

import pytest

# This module imports openem.solver at the top level, which imports cupy at its own top level: on a
# CPU-only machine that would crash during collection.
pytest.importorskip("cupy", reason="needs CUDA and cupy")

from openem import cpml
from openem import scene as scene_mod

C0 = 1.0 / np.sqrt(cpml.EPSILON_0 * cpml.MU_0)
UM = 1e-6
FREQS = [150e12, 200e12, 250e12]


def _sim(pol: str, symmetry=(0, 0, 0), center=(0.012, 0.007, 0.0)):
    return td.Simulation(
        size=(3.0, 3.0, 3.0),
        grid_spec=td.GridSpec.uniform(dl=0.05),
        run_time=2.0e-13,
        symmetry=symmetry,
        sources=[td.PointDipole(
            center=center, polarization=pol,
            source_time=td.GaussianPulse(freq0=200e12, fwidth=50e12))],
        monitors=[td.FluxMonitor(center=(0, 0, 0), size=(1.0, 1.0, 1.0),
                                 freqs=FREQS, name="box")],
        boundary_spec=td.BoundarySpec.all_sides(boundary=td.PML()),
    )


# ---------------------------------------------------------------- extraction (CPU)

def test_magnetic_dipole_accepted_and_flagged():
    sc = scene_mod.from_simulation(_sim("Hz"))
    d = sc.dipoles[0]
    assert d.magnetic and d.component == 2


def test_magnetic_stencil_weights_sum_to_one():
    """Restored against the **H dual volume**, the weights sum to 1. Using the E volume by mistake
    would show up on a non-uniform grid; on this uniform grid what is checked is the chain itself,
    that a volume is divided out and multiplied back in."""
    sc = scene_mod.from_simulation(_sim("Hz"))
    d = sc.dipoles[0]
    vol = sc.grid.yee_h_dual_volumes()[d.component]
    w = np.array([c * float(vol[i, j, k]) / scene_mod.UM
                  for (i, j, k), c in zip(d.indices, d.coef)])
    np.testing.assert_allclose(w.sum(), 1.0, rtol=1e-12)
    assert np.all(w > 0)


def test_hz_on_pmc_plane_doubles():
    """An Hz source on a symmetry=+1 (PMC) z plane: H has eigenvalue +1 along the reflection axis,
    so the image lands back on the original position and the source doubles (the same rule carried
    over to the pseudovector table). Bandstructure sim_0 is exactly this configuration."""
    plain = scene_mod.from_simulation(_sim("Hz", center=(0.012, 0.007, 0.0)))
    mirrored = scene_mod.from_simulation(
        _sim("Hz", symmetry=(0, 0, 1), center=(0.012, 0.007, 0.0)))
    v = plain.grid.yee_h_dual_volumes()[2]
    tot = lambda sc: sum(c * float(v[i, j, k])
                         for (i, j, k), c in zip(sc.dipoles[0].indices,
                                                 sc.dipoles[0].coef))
    np.testing.assert_allclose(tot(mirrored), 2.0 * tot(plain), rtol=1e-12)


def test_ez_on_pmc_plane_still_contradictory():
    """Control case: **Ez** at the same position (a polar vector, eigenvalue -1 along the axis) is
    contradictory input on a PMC plane, and the existing fail-closed behaviour must not be
    disturbed by the change to the pseudovector table."""
    with pytest.raises(NotImplementedError, match="symmetry plane"):
        scene_mod.from_simulation(_sim("Ez", symmetry=(0, 0, 1)))


# ---------------------------------------------------------------- solving (GPU, seconds)

@pytest.fixture(scope="module")
def powers():
    from openem import flux as flux_mod, solver
    from openem.device import Kernels
    from openem.scene.projection import normalization

    k = Kernels()
    out = {}
    for pol in ("Hz", "Ez"):
        sim = _sim(pol)
        sc = scene_mod.from_simulation(sim)
        res = solver.run(sc, kernels=k, verbose=False)
        box = sc.box_flux_monitors[0]
        faces = {m.name: m for m in sc.flux_monitors}
        m0 = faces[box.face_names[0]]
        norm2 = np.abs(normalization(sim, m0.freqs)) ** 2
        total = np.zeros(m0.freqs.size)
        for nm in box.face_names:
            m = faces[nm]
            total += flux_mod.monitor_flux(res.phasors[nm].data, sc.grid, m) \
                     * m.normal_dir * norm2
        out[pol] = total
    return out


def test_magnetic_dipole_power_matches_analytic(powers):
    """An absolute criterion plus "the error matches the electric side point by point".

    For this configuration (30 cells per wavelength, flux box 0.5 µm from the source) the measured
    discretization floor is +5.0 to 5.3%, and halving dl brings it down to +2.8%, converging towards
    the analytic value. **The already validated electric side is off by the same +4.9 to 5.2% under
    exactly the same settings**, so those 5% belong to the test configuration, not to the magnetic
    path. The real parameter-free criterion is the second assert: the error of the magnetic path
    must agree frequency point by frequency point with the electric path that was validated against
    the reference output (measured difference <= 1.3e-3). Any amplitude, volume or injection error
    on the H side breaks it.
    """
    w = 2.0 * np.pi * np.asarray(FREQS)
    ref_m = cpml.EPSILON_0 * w ** 2 * UM ** 2 / (12.0 * np.pi * C0)
    ref_e = cpml.MU_0 * w ** 2 * UM ** 2 / (12.0 * np.pi * C0)
    rel_m = powers["Hz"] / ref_m - 1.0
    rel_e = powers["Ez"] / ref_e - 1.0
    assert np.abs(rel_m).max() < 8e-2, f"P_m differs from analytic by {rel_m}"
    assert np.abs(rel_m - rel_e).max() < 5e-3,         f"magnetic path departs from the electric baseline: Hz {rel_m} vs Ez {rel_e}"


def test_duality_power_ratio(powers):
    """P_m / P_e = ε₀/μ₀: geometry and discretization errors cancel in the ratio (measured
    difference ~6e-4)."""
    ratio = powers["Hz"] / powers["Ez"]
    ref = cpml.EPSILON_0 / cpml.MU_0
    rel = np.abs(ratio / ref - 1.0)
    assert rel.max() < 1e-2, f"duality ratio off by {rel}"
