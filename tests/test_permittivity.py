# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""PermittivityMonitor, and the out-of-range guard on colocation."""
import os
import json

import h5py
import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import colocate
from openem import scene as scene_mod
from openem.grid import Axis, Grid

TABLEA = os.environ.get("OPENEM_TABLEA", "/nonexistent")  # external reference archive; without
#                                                          OPENEM_TABLEA these cases skip
CAV = f"{TABLEA}/CavityFOM/Tidy3D/simulation.json"
PAIR = {"eps_xx": "Ex", "eps_yy": "Ey", "eps_zz": "Ez"}


@pytest.fixture(scope="module")
def cav():
    """Keep only the PermittivityMonitor: field_vol carries apodization, which is not supported yet."""
    sim = td.Simulation.from_file(CAV)
    return sim, scene_mod.from_simulation(sim.updated_copy(
        monitors=[m for m in sim.monitors if isinstance(m, td.PermittivityMonitor)]))


def test_permittivity_matches_archive_to_machine_precision(cav):
    """Comparing eps needs no solver run: a PermittivityMonitor stores values at their native Yee
    positions (``monitor.py:1270``), the same convention as our internal ``eps_ex/ey/ez``.

    This also locks the subpixel-averaged eps, the Yee coordinate alignment and the reduced-domain
    coordinate mapping.
    """
    sim, sc = cav
    mon = sc.permittivity_monitors[0]
    with h5py.File(f"{TABLEA}/CavityFOM/Tidy3D/data.hdf5", "r") as f:
        raw = f["JSON_STRING"][()]
        js = json.loads(raw.decode() if isinstance(raw, bytes) else raw)
        i = next(n for n, it in enumerate(js["data"])
                 if it.get("type") == "PermittivityData")
        g = f[f"data/{i}"]
        for (name, comp), arr in zip(PAIR.items(), sc.eps_box(mon)):
            ref = np.real(np.squeeze(
                np.asarray(g[f"{name}/__xarray_dataarray_variable__"])))
            dst = [np.asarray(g[f"{name}/{a}"], float) * scene_mod.UM for a in "xyz"]
            src = colocate.sample_coords(sc.grid, mon.origin, mon.box, comp)
            got = colocate.interp_to(arr[None], src, dst)[0]
            rel = np.max(np.abs(got - ref)) / np.max(np.abs(ref))
            assert rel < 1e-12, f"{name} relative difference {rel:.3e}"


def test_permittivity_ignores_apodization(cav):
    """eps does not depend on time, so apodization must not cause a refusal; interval_space is still
    checked.
    """
    sim, _ = cav
    mon = next(m for m in sim.monitors if isinstance(m, td.PermittivityMonitor))
    bad = sim.updated_copy(monitors=[mon.updated_copy(interval_space=(2, 1, 1))])
    with pytest.raises(NotImplementedError, match="interval_space"):
        scene_mod.from_simulation(bad)


def test_colocate_refuses_out_of_range():
    """Going out of range by **more than one sample spacing** must raise rather than clamp silently.

    Within half a cell it is allowed: the outermost primal boundary that ``colocate=True`` asks for
    sits half a cell outside the outermost cell center, and there is no sample to interpolate across
    that half cell, so the nearest value is all there is. Telling those two apart is the whole point
    of this test: one scene was short by two full cells, got silently clamped, and reported a
    relative difference of 3.2.
    """
    e = np.linspace(0.0, 4.0, 5)
    g = Grid(*(Axis(e, "PML", "PML") for _ in range(3)))
    coords = colocate.sample_coords(g, (0, 0, 0), (4, 4, 4), "Ez")
    vals = np.zeros((1, 4, 4, 4))
    step = float(coords[2][1] - coords[2][0])

    # Half a cell: allowed, take the nearest endpoint
    near = [coords[0], coords[1], np.array([coords[2][-1] + 0.5 * step])]
    colocate.interp_to(vals, coords, near)

    # Two and a half cells: raises
    far = [coords[0], coords[1], np.array([coords[2][-1] + 2.5 * step])]
    with pytest.raises(ValueError, match="outside the sampled range"):
        colocate.interp_to(vals, coords, far)


def test_monitor_box_has_margin_on_both_sides(cav):
    """The Yee points Tidy3D reports extend past the monitor's nominal bounds, so both ends need a
    margin.
    """
    sim, sc = cav
    mon = sc.permittivity_monitors[0]
    tdm = next(m for m in sim.monitors if m.name == mon.name)
    lo, hi = tdm.bounds
    for ax in range(3):
        src = colocate.sample_coords(sc.grid, mon.origin, mon.box, "Ex")[ax]
        assert src[0] / scene_mod.UM < float(lo[ax]), f"{'xyz'[ax]}: no margin at the low end"
        assert src[-1] / scene_mod.UM > float(hi[ax]), f"{'xyz'[ax]}: no margin at the high end"
