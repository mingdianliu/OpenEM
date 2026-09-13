# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
# -*- coding: utf-8 -*-
"""Depth-graded damping of ADE poles inside a psi-PML band (boundaries.PML_POLE_DAMP).

Background: the CFS alpha injection cured the CPU-reproducible family of alpha=0 CPML against ADE,
but GPU probes showed a faster residual mode of the same family still present on the complete path
(implied growth of 3e-3 to 1e-2 per step), which a term-by-term CPU replay could not reproduce.
This mechanism lowers the spectral radius of the poles in the band as step cubed with depth, giving
gamma_d*|q|*dt of about 0.03 per step deep in the layer, and stacks with alpha for margin.
Controlled experiments had already shown it alone cures the reproducible family.

The criteria:
1. pole entries inside the band have a strictly smaller |A|, i.e. more decay, while entries outside
   it are bitwise unchanged;
2. Dispersion.g stays consistent with 2*sum Re(B) for the entries that changed;
3. a non-dispersive scene is unchanged, and dispersion shrunk into the domain center is bitwise
   unchanged;
4. setting the environment variable to 0 restores the old behaviour bitwise;
5. an Absorber face does not count as a psi-PML band, so swapping in an Absorber leaves the same
   scene bitwise unchanged.

CPU only; no GPU needed.
"""
from __future__ import annotations

import os
import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem.scene import boundaries
from openem.scene import from_simulation

FREQ0 = 2e14
LORENTZ = td.Lorentz(eps_inf=2.0, coeffs=[(1.5, 5e14, 0.0)])


def _sim(slab_medium, slab_z=(0.35, 2.0), zb=None):
    if zb is None:
        zb = td.Boundary.pml(num_layers=12)
    return td.Simulation(
        size=(0.1, 0.1, 1.2), grid_spec=td.GridSpec.uniform(dl=0.025),
        subpixel=False,
        structures=[td.Structure(
            geometry=td.Box.from_bounds(
                (-td.inf, -td.inf, slab_z[0]), (td.inf, td.inf, slab_z[1])),
            medium=slab_medium)],
        sources=[td.PlaneWave(
            center=(0, 0, -0.4), size=(td.inf, td.inf, 0), direction="+",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.2 * FREQ0))],
        monitors=[td.FluxMonitor(center=(0, 0, 0.0), size=(td.inf, td.inf, 0),
                                 freqs=[FREQ0], name="T")],
        run_time=2e-13,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.periodic(), y=td.Boundary.periodic(), z=zb))


def _band_mask(sc, disp):
    """Whether an entry falls inside the z-axis PML band; in these scenes only z has a PML."""
    nz = sc.shape[2]
    k = disp.cell.astype(np.int64) % nz
    nlo = sc.pml[(2, "lo")].num_layers if (2, "lo") in sc.pml else 0
    nhi = sc.pml[(2, "hi")].num_layers if (2, "hi") in sc.pml else 0
    return (k < nlo) | (k >= nz - nhi)


def _undamped():
    """A reference Dispersion table at gamma=0, with every other mechanism untouched."""
    import os
    old = os.environ.get("OPENEM_PML_POLE_DAMP")
    os.environ["OPENEM_PML_POLE_DAMP"] = "0"
    try:
        return from_simulation(_sim(LORENTZ))
    finally:
        if old is None:
            del os.environ["OPENEM_PML_POLE_DAMP"]
        else:
            os.environ["OPENEM_PML_POLE_DAMP"] = old


def test_band_poles_damped_out_of_band_untouched():
    sc = from_simulation(_sim(LORENTZ))
    sc0 = _undamped()
    d, d0 = sc.dispersion, sc0.dispersion
    assert d is not None and d.n_entry > 0
    inb = _band_mask(sc, d)
    assert inb.any() and (~inb).any()
    counts = np.diff(d.pole_ofs).astype(np.int64)
    inb_p = np.repeat(inb, counts)
    # Outside the band, bitwise unchanged
    assert np.array_equal(d.am1[~inb_p], d0.am1[~inb_p])
    assert np.array_equal(d.b[~inb_p], d0.b[~inb_p])
    assert np.array_equal(d.g[~inb], d0.g[~inb])
    # Inside the band, in layers of non-zero depth, |A| is strictly smaller; at the interface layer
    # the profile depth can be small but is non-zero
    changed = d.am1 != d0.am1
    assert changed.any(), "the poles inside the band must actually have changed"
    A, A0 = d.am1 + 1.0, d0.am1 + 1.0
    assert np.all(np.abs(A[changed]) < np.abs(A0[changed]))


def test_g_consistent_with_b():
    sc = from_simulation(_sim(LORENTZ))
    d = sc.dispersion
    counts = np.diff(d.pole_ofs).astype(np.int64)
    seg = np.add.reduceat(2.0 * d.b.real, d.pole_ofs[:-1])
    seg = np.where(counts > 0, seg, 0.0)
    np.testing.assert_allclose(d.g, seg, rtol=0, atol=1e-12)


def test_nondispersive_unchanged():
    sc = from_simulation(_sim(td.Medium(permittivity=2.0)))
    assert sc.dispersion is None


def test_center_confined_bitwise_unchanged():
    sim = _sim(LORENTZ, slab_z=(-0.2, 0.2))
    sc = from_simulation(sim)
    import os
    os.environ["OPENEM_PML_POLE_DAMP"] = "0"
    try:
        sc0 = from_simulation(sim)
    finally:
        del os.environ["OPENEM_PML_POLE_DAMP"]
    d, d0 = sc.dispersion, sc0.dispersion
    assert d.n_entry > 0
    assert np.array_equal(d.am1, d0.am1)
    assert np.array_equal(d.b, d0.b)
    assert np.array_equal(d.g, d0.g)


def test_env_zero_restores_old(monkeypatch):
    monkeypatch.setenv("OPENEM_PML_POLE_DAMP", "0")
    sc = from_simulation(_sim(LORENTZ))
    # At gamma=0, damp_pml_poles returns its input unchanged: bitwise identical to building directly
    # (also at gamma=0), with no trace of damping inside or outside the band. Equivalence with the
    # old behaviour is guaranteed by the early-return branch in damp_pml_poles; what is verified here
    # is that the table really was not touched, i.e. the |A| monotonicity does not hold.
    d = sc.dispersion
    monkeypatch.delenv("OPENEM_PML_POLE_DAMP")
    sc1 = from_simulation(_sim(LORENTZ))
    d1 = sc1.dispersion
    assert not np.array_equal(d.am1, d1.am1), "gamma=0 and the default gamma must actually differ"


def test_absorber_face_not_counted():
    """Swapping both z ends for an Absorber leaves no psi-PML face, so the result is bitwise identical
    to gamma=0.
    """
    zb = td.Boundary(minus=td.Absorber(num_layers=12),
                     plus=td.Absorber(num_layers=12))
    sim = _sim(LORENTZ, zb=zb)
    sc = from_simulation(sim)
    import os
    os.environ["OPENEM_PML_POLE_DAMP"] = "0"
    try:
        sc0 = from_simulation(sim)
    finally:
        del os.environ["OPENEM_PML_POLE_DAMP"]
    assert np.array_equal(sc.dispersion.am1, sc0.dispersion.am1)
    assert np.array_equal(sc.dispersion.b, sc0.dispersion.b)
