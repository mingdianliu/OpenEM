# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Residual reflection from the PML: below 1e-08 at normal incidence.

This locks the calibration of sigma and kappa. A weakened PML fails **silently**: the field does
not diverge and the result is not obviously wrong, there is simply a little more reflection, which
then contaminates every flux criterion.

**Normal incidence only.** Grazing incidence, i.e. a periodic structure above its Rayleigh anomaly,
is a different matter: there our 12 layers are worth about 24 of Tidy3D's.
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

F0, DL = 5.0e14, 0.01


def _pml_reflection(num_layers: int, kernels) -> float:
    """Vacuum, no structures, one-way TF/SF injection; the backward flux measured upstream of the
    source is the PML reflection.
    """
    sim = td.Simulation(
        size=(4 * DL, 4 * DL, 2.0), grid_spec=td.GridSpec.uniform(dl=DL),
        sources=[td.PlaneWave(
            center=(0, 0, 0.0), size=(td.inf, td.inf, 0), direction="-",
            source_time=td.GaussianPulse(freq0=F0, fwidth=0.2 * F0))],
        monitors=[td.FluxMonitor(center=(0, 0, 0.4), size=(td.inf, td.inf, 0),
                                 freqs=[F0], name="back")],
        run_time=6e-13,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.periodic(), y=td.Boundary.periodic(),
            z=td.Boundary.pml(num_layers=num_layers)))
    sc = scene_mod.from_simulation(sim)
    src = sc.sources[0]
    dn = int(np.clip(src.plane_index + src.direction * 12, 10, sc.shape[2] - 11))
    probe = replace(sc.flux_monitors[0], name="_dn", plane_index=dn, normal_dir=1)
    res = solver.run(replace(sc, flux_monitors=list(sc.flux_monitors) + [probe]),
                     kernels=kernels, verbose=False)
    pin = float(np.abs(flux_mod.plane_flux(res.phasors["_dn"].data, sc.grid))[0])
    back = float(np.abs(flux_mod.plane_flux(res.phasors["back"].data, sc.grid))[0])
    return back / pin


def test_pml_reflection_is_negligible_at_normal_incidence():
    """12 layers, the most common configuration in the validation set, measured 3.8e-09."""
    r = _pml_reflection(12, Kernels())
    assert r < 1e-8, f"residual reflection of a 12-layer PML is {r:.3e}, should be < 1e-8"


def test_pml_reflection_improves_with_layers():
    """More layers must monotonically improve it, which guards against the sigma grading being
    written independently of the layer count.
    """
    k = Kernels()
    r6, r12 = _pml_reflection(6, k), _pml_reflection(12, k)
    assert r12 < r6 / 100, f"6 layers {r6:.3e} to 12 layers {r12:.3e} is less than two orders of improvement"


def _sim_z(minus) -> "td.Simulation":
    """Vacuum with TF/SF injecting towards -z; the low-z boundary is given by ``minus``."""
    return td.Simulation(
        size=(4 * DL, 4 * DL, 2.0), grid_spec=td.GridSpec.uniform(dl=DL),
        sources=[td.PlaneWave(
            center=(0, 0, 0.0), size=(td.inf, td.inf, 0), direction="-",
            source_time=td.GaussianPulse(freq0=F0, fwidth=0.2 * F0))],
        monitors=[td.FluxMonitor(center=(0, 0, 0.4), size=(td.inf, td.inf, 0),
                                 freqs=[F0], name="back")],
        run_time=6e-13,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.periodic(), y=td.Boundary.periodic(),
            z=td.Boundary(minus=minus, plus=td.PML())))


def test_pec_wall_reflects_fully():
    """A PEC wall reflects everything: the reflected flux measured upstream equals the incident flux
    of the reference run.

    The normalization cannot use the PEC run's own downstream probe: under total reflection the
    downstream field is a standing wave whose net flux is about 0. The reference run, with a PML at
    low z, measures the injected power at that same downstream position.
    Mixing the axis, a PEC wall at low z and a PML at high z, also exercises the path through
    refuse_unsupported that admits such a combination.
    """
    k = Kernels()
    ref = scene_mod.from_simulation(_sim_z(td.PML()))
    src = ref.sources[0]
    dn = int(np.clip(src.plane_index + src.direction * 12, 10, ref.shape[2] - 11))
    probe = replace(ref.flux_monitors[0], name="_dn", plane_index=dn, normal_dir=1)
    res_ref = solver.run(replace(ref, flux_monitors=list(ref.flux_monitors) + [probe]),
                         kernels=k, verbose=False)
    pin = float(np.abs(flux_mod.plane_flux(res_ref.phasors["_dn"].data, ref.grid))[0])

    sc = scene_mod.from_simulation(_sim_z(td.PECBoundary()))
    assert sc.grid.axes[2].boundary_lo == "PECBoundary"
    res = solver.run(sc, kernels=k, verbose=False)
    back = float(np.abs(flux_mod.plane_flux(res.phasors["back"].data, sc.grid))[0])
    r = back / pin
    assert 0.98 < r < 1.02, f"PEC wall reflection |r|^2 = {r:.4f}, should be 1"
