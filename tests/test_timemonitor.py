# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""FieldTimeMonitor and colocation.

Every criterion for the time convention comes from the reference outputs and from the Tidy3D
source.
"""

import os
import json

import h5py
import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import colocate
from openem import scene as scene_mod
from openem.grid import Axis, Grid
from openem.model import COMPONENTS

TABLEA = os.environ.get("OPENEM_TABLEA", "/nonexistent")  # external reference archive; without
#                                                          OPENEM_TABLEA these cases skip
BC = f"{TABLEA}/BoundaryConditions/Tidy3D/simulation.json"
RF = f"{TABLEA}/ResonanceFinder/Tidy3D/simulation.json"


def _archive_t(case: str, kind: str, comp: str) -> np.ndarray:
    with h5py.File(f"{TABLEA}/{case}/Tidy3D/data.hdf5", "r") as f:
        raw = f["JSON_STRING"][()]
        js = json.loads(raw.decode() if isinstance(raw, bytes) else raw)
        i = next(n for n, it in enumerate(js["data"]) if it.get("type") == kind)
        return np.asarray(f[f"data/{i}/{comp}/t"], dtype=np.float64)


@pytest.fixture(scope="module")
def bc():
    return scene_mod.from_file(BC)


# ------------------------------------------------------------ sampling window

def test_step_window_from_tidy3d_time_inds(bc):
    """The window uses Tidy3D's ``time_inds`` (``monitor.py:331``) directly rather than being
    reimplemented.
    """
    sim = td.Simulation.from_file(BC)
    mon = bc.field_time_monitors[0]
    beg, end = (int(v) for v in sim.monitors[0].time_inds(np.asarray(sim.tmesh)))
    assert (mon.step_begin, mon.step_end) == (beg, end)
    assert mon.interval == int(sim.monitors[0].interval) == 50


def test_sample_times_match_archive(bc):
    """The reference ``t`` is ``m*dt``, with ``m`` starting at ``tind_beg`` and stepping by
    ``interval``.
    """
    sim = td.Simulation.from_file(BC)
    mon = bc.field_time_monitors[0]
    t = _archive_t("BoundaryConditions", "FieldTimeData", "Ex")
    want = mon.step_begin + np.arange(t.size) * mon.interval
    np.testing.assert_allclose(t / sim.dt, want, atol=1e-6)


def test_resonancefinder_window_starts_late():
    """One scene starts recording only at step 9736, so ``start`` has to be respected; otherwise the
    source transient of the first ten thousand steps mixes in, which is exactly what a resonance
    finder is trying to avoid.
    """
    sc = scene_mod.from_file(RF)
    mon = sc.field_time_monitors[0]
    assert (mon.step_begin, mon.step_end, mon.interval) == (9736, 97354, 1)
    assert mon.num_slots == 87618
    t = _archive_t("ResonanceFinder", "FieldTimeData", "Ey")
    assert t.size == mon.num_slots


def test_only_requested_components(bc):
    """``fields`` is a subset: one scene wants only Ey, and nothing else is allocated or sampled."""
    rf = scene_mod.from_file(RF)
    assert [COMPONENTS[c] for c in rf.field_time_monitors[0].comps] == ["Ey"]
    assert len(bc.field_time_monitors[0].comps) == 6


def test_slot_of_respects_window_and_interval(bc):
    mon = bc.field_time_monitors[0]
    assert mon.slot_of(0) == 0
    assert mon.slot_of(50) == 1
    assert mon.slot_of(49) is None            # not in the sampled set
    assert mon.slot_of(-1) is None            # before start
    assert mon.slot_of(mon.step_end) is None  # after stop


def test_refuses_interval_space():
    """Spatial thinning is not implemented; ignoring it silently would make the output point count
    disagree with the reference.
    """
    sim = td.Simulation.from_file(BC)
    bad = sim.updated_copy(
        monitors=[sim.monitors[0].updated_copy(interval_space=(2, 1, 1))]
    )
    with pytest.raises(NotImplementedError, match="interval_space"):
        scene_mod.from_simulation(bad)


def test_time_buffer_guard_fails_closed():
    """Exceeding the buffer ceiling must raise rather than truncate silently; a truncation would
    leave a silent gap in the comparison.
    """
    sim = td.Simulation.from_file(BC)
    huge = sim.updated_copy(
        monitors=[sim.monitors[0].updated_copy(interval=1, size=(td.inf, td.inf, td.inf))]
    )
    with pytest.raises(NotImplementedError, match="time-domain buffer"):
        scene_mod.from_simulation(huge)


# ------------------------------------------------------------ colocation

def test_monitor_box_covers_three_cells_on_zero_axis(bc):
    """A zero-thickness axis takes a three-cell neighbourhood: a component may sit at a cell center
    or on a cell boundary, and three cells cover both.
    """
    mon = bc.field_time_monitors[0]
    assert mon.box[1] == 3


def test_at_point_interp_reproduces_linear_field():
    """On a linear field, interpolation must reproduce the value exactly; snapping does not."""
    e = np.linspace(0.0, 4.0, 5)
    g = Grid(*(Axis(e, "PML", "PML") for _ in range(3)))
    coords = colocate.sample_coords(g, (0, 0, 0), (4, 4, 4), "Ez")   # ctr_z
    # Build a field that is linear along x
    vals = np.zeros((1, 4, 4, 4))
    vals[0] = coords[0][:, None, None] * np.ones((1, 4, 4))
    x = 1.37
    got = colocate.at_point(vals, coords, (x, coords[1][1], coords[2][1]), "interp")
    np.testing.assert_allclose(got[0], x, rtol=1e-12)
    snapped = colocate.at_point(vals, coords, (x, coords[1][1], coords[2][1]), "snap")
    assert abs(snapped[0] - x) > 0.05, "snapping lands on the nearest grid point, not the exact value"


def test_interp_to_is_identity_on_matching_coords():
    e = np.linspace(0.0, 4.0, 5)
    g = Grid(*(Axis(e, "Periodic", "Periodic") for _ in range(3)))
    coords = colocate.sample_coords(g, (0, 0, 0), (4, 4, 4), "Hx")
    vals = np.arange(64, dtype=float).reshape(1, 4, 4, 4)
    out = colocate.interp_to(vals, coords, coords)
    np.testing.assert_allclose(out, vals, rtol=1e-12)


# ------------------------------------------------------------ end to end

def test_boundaryconditions_time_series_matches_archive(bc):
    """Step-by-step comparison against a reference run, locking three conventions:

    - ``t = m*dt`` with ``m`` starting at ``tind_beg``
    - **H is time-averaged onto the whole step**; without the average the Hy residual is 1.83e-1,
      and with it 5.9e-3
    - time-domain data is not normalized in frequency, so the two sides differ only by the 1e6 of
      V/um against V/m
    """
    pytest.importorskip("cupy")
    from openem import solver

    mon = bc.field_time_monitors[0]
    with h5py.File(f"{TABLEA}/BoundaryConditions/Tidy3D/data.hdf5", "r") as f:
        raw = f["JSON_STRING"][()]
        js = json.loads(raw.decode() if isinstance(raw, bytes) else raw)
        i = next(n for n, it in enumerate(js["data"]) if it.get("type") == "FieldTimeData")
        g = f[f"data/{i}"]
        ref = {c: np.asarray(g[f"{c}/__xarray_dataarray_variable__"], float)
               for c in COMPONENTS}
        dst = [np.asarray(g[f"Ex/{a}"], float) * scene_mod.UM for a in "xyz"]

    res = solver.run(bc, use_shutoff=True, verbose=False)
    buf = res.time_samples[mon.name]
    n = min(buf.shape[1], ref["Ex"].shape[-1])
    assert n >= 20, f"only {n} overlapping samples, too few"

    gmax = max(np.max(np.abs(ref[c][..., :n])) for c in COMPONENTS)
    worst = 0.0
    for ci, cidx in enumerate(mon.comps):
        c = COMPONENTS[cidx]
        coords = colocate.sample_coords(bc.grid, mon.origin, mon.box, c)
        ours = colocate.interp_to(buf[ci, :n], coords, dst) / 1e6
        a = np.moveaxis(ref[c][..., :n], -1, 0)
        if np.max(np.abs(a)) <= 1e-4 * gmax:      # noise in the reference
            continue
        worst = max(worst, np.max(np.abs(ours - a)) / gmax)
    assert worst < 1e-2, f"time-domain residual {worst:.3e} exceeds 1e-2"
