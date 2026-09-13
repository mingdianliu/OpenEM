# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""The O(dl^2) floor of subpixel averaging, measured against the analytic Fresnel solution.

Two things are locked down:
1. With the interface at a **cell center**, the worst position, the reflectance error at dl=0.025
   is of order 1e-03. That is the accuracy floor of every flux example, and on a small quantity
   like R it divides out into a 5% relative error.
2. Halving the cell size drops the error to a quarter, i.e. **second order**. Falling to first
   order means the arithmetic/harmonic mixing of subpixel averaging (the beta in
   ``scene/weights.py``) has degenerated.

This test is the explanation for the "5.3% error" on one grating example: that is not a bug but an
absolute power error of 2e-03 at dl=0.025 landing on an R of about 0.06.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")
pytest.importorskip("cupy")

from openem import flux as flux_mod
from openem import scene as scene_mod
from openem import solver
from openem.device import Kernels

N = 1.5
FREQ0 = 3.6e14
R_EXACT = ((1 - N) / (1 + N)) ** 2      # 0.04, exact


def _reflectance(frac: float, dl: float, kernels) -> float:
    """Put the interface at position ``frac`` within a cell and measure the normalized reflectance.

    ``frac`` is measured from the grid line as 0. The half space is a box 10 um thick, far beyond
    the simulation domain, so only **one** interface, the top surface, takes part.
    """
    z0 = -1.5 + round((-0.4 + 1.5) / dl) * dl + frac * dl
    sim = td.Simulation(
        size=(4 * dl, 4 * dl, 3.0),
        grid_spec=td.GridSpec.uniform(dl=dl),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, z0 - 5.0), size=(td.inf, td.inf, 10.0)),
            medium=td.Medium(permittivity=N * N))],
        sources=[td.PlaneWave(
            center=(0, 0, 0.6), size=(td.inf, td.inf, 0), direction="-",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.1 * FREQ0))],
        monitors=[td.FluxMonitor(center=(0, 0, 0.9), size=(td.inf, td.inf, 0),
                                 freqs=[FREQ0], name="R")],
        run_time=8e-13,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.periodic(), y=td.Boundary.periodic(),
            z=td.Boundary.pml(num_layers=int(round(12 * 0.025 / dl)))))
    sc = scene_mod.from_simulation(sim)
    src = sc.sources[0]
    dn = int(np.clip(src.plane_index + src.direction * 16, 14, sc.shape[2] - 15))
    probe = replace(sc.flux_monitors[0], name="_dn", plane_index=dn, normal_dir=1)
    # Incident power: the same scene with no structures, so the absolute amplitude of the plane wave
    # cancels entirely
    ref = solver.run(scene_mod.lossless(scene_mod.homogeneous(
        replace(sc, flux_monitors=list(sc.flux_monitors) + [probe]))),
        kernels=kernels, verbose=False)
    pin = float(np.abs(flux_mod.plane_flux(ref.phasors["_dn"].data, sc.grid))[0])
    res = solver.run(sc, kernels=kernels, verbose=False)
    return abs(float(flux_mod.plane_flux(res.phasors["R"].data, sc.grid)[0]) / pin)


def test_subpixel_reflectance_is_second_order():
    k = Kernels()
    e_coarse = _reflectance(0.5, 0.025, k) - R_EXACT
    e_fine = _reflectance(0.5, 0.0125, k) - R_EXACT
    # The floor itself: of order 1e-03 on the coarse grid
    assert 5e-4 < abs(e_coarse) < 3e-3, f"the coarse-grid error {e_coarse:+.3e} is not of the expected magnitude"
    # Second order: halving the cell size drops it to a quarter
    order = np.log2(abs(e_coarse) / abs(e_fine))
    assert 1.7 < order < 2.3, (
        f"convergence order {order:.2f}, should be 2. Coarse {e_coarse:+.3e}, fine {e_fine:+.3e}")


def test_quarter_cell_is_the_sweet_spot():
    """With the interface a quarter of a cell in, the error is nearly zero, a direct consequence of
    the mixing weight beta.

    This guards against a silent degeneration such as beta being replaced by a constant: the quarter
    position would then no longer be special and its error would be of the same order as at the cell
    center.
    """
    k = Kernels()
    e_quarter = abs(_reflectance(0.25, 0.025, k) - R_EXACT)
    e_center = abs(_reflectance(0.5, 0.025, k) - R_EXACT)
    assert e_quarter < 0.25 * e_center, (
        f"the quarter position {e_quarter:.3e} is not clearly better than the cell center {e_center:.3e}")
