# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Closed-box flux: signs and symmetry across all six faces, the two y-normal ones included.

A point dipole in vacuum with a closed box around it. Two criteria with no tunable parameters:

* **The outward flux through every face must be positive.** The source is inside the box, so
  energy only leaves. Get the Poynting sign on a y normal wrong (the axis=1 negation in
  ``flux.monitor_flux``) and the two y faces go negative, which is exactly how it once failed
  silently.
* **Symmetry.** With the dipole at the box center and Ez polarization, the field is symmetric
  under x to -x, y to -y and z to -z, so the three opposite pairs of faces must match.

The end-to-end absolute check lives elsewhere; this file only locks the signs and the symmetry,
which are cheap and are precisely the parts that break without saying anything.
"""

from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")
pytest.importorskip("cupy")

from openem import flux as flux_mod
from openem import scene as scene_mod
from openem import solver
from openem.device import Kernels

FREQ0 = 3.0e14


def _closed_box_fluxes():
    sim = td.Simulation(
        size=(1.6, 1.6, 1.6), grid_spec=td.GridSpec.uniform(dl=0.025),
        sources=[td.PointDipole(
            center=(0, 0, 0), polarization="Ez",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.2 * FREQ0))],
        monitors=[td.FluxMonitor(center=(0, 0, 0), size=(0.8, 0.8, 0.8),
                                 freqs=[FREQ0], name="box")],
        run_time=1.5e-13,
        boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=10)))
    sc = scene_mod.from_simulation(sim)
    res = solver.run(sc, use_shutoff=False, kernels=Kernels(), verbose=False)
    faces = {m.name: m for m in sc.flux_monitors}
    box = sc.box_flux_monitors[0]
    out = {}
    for nm in box.face_names:
        m = faces[nm]
        out[nm[-2:]] = float(
            flux_mod.monitor_flux(res.phasors[nm].data, sc.grid, m)[0] * m.normal_dir)
    return out


def test_all_faces_radiate_outward():
    out = _closed_box_fluxes()
    for face, v in out.items():
        assert v > 0, f"outward flux through face {face} is {v:.3e} <= 0, but the source is " \
                      f"inside the box so energy can only leave. A negative y face means the " \
                      f"axis=1 Poynting sign has been lost again. All faces: {out}"


def test_opposite_faces_are_symmetric():
    out = _closed_box_fluxes()
    for a, b in (("x-", "x+"), ("y-", "y+"), ("z-", "z+")):
        r = out[a] / out[b]
        assert abs(r - 1) < 0.02, f"{a}/{b} = {r:.4f}; with the dipole at the box center these " \
                                  f"should be symmetric. All faces: {out}"
