# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""How ModeMonitor is wired up.

The real physical criterion needs a full FDTD run plus the ModeSolver, which is too slow to keep in
the tests. What is locked here are three cheap things that break easily and quietly:

* a ModeMonitor expands into one planar ``FieldMonitor``, with a matching name
* that FieldMonitor must have ``colocate=True``: the mode overlap integral needs E and H sampled at
  the same point, whereas the in-place Yee profile is the form needed for **injection**
* a volumetric ModeMonitor is rejected by **Tidy3D itself**, so we do not check it again
"""

from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import scene as scene_mod
from openem.scene import modes

FREQ0 = 1.934e14


def _sim(mode_size=(0.0, 1.2, 1.0)) -> "td.Simulation":
    return td.Simulation(
        size=(2.0, 1.6, 1.4), grid_spec=td.GridSpec.uniform(dl=0.05),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, 0), size=(td.inf, 0.5, 0.22)),
            medium=td.Medium(permittivity=12.1))],
        sources=[td.PointDipole(
            center=(-0.6, 0, 0), polarization="Ey",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.1 * FREQ0))],
        monitors=[td.ModeMonitor(
            center=(0.5, 0, 0), size=mode_size, freqs=[FREQ0],
            mode_spec=td.ModeSpec(num_modes=2), name="m")],
        run_time=2e-13,
        boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=8)))


def test_mode_monitor_expands_to_field_monitor():
    sc = scene_mod.from_simulation(_sim())
    assert [m.name for m in sc.mode_monitors] == ["m"]
    spec = sc.mode_monitors[0]
    assert spec.plane_name == "m" + modes.PLANE_SUFFIX
    assert spec.plane_name in {m.name for m in sc.field_monitors}


def test_plane_monitor_is_colocated():
    """``colocate=False`` would put E and H at different points, which breaks the overlap
    integral."""
    mon = _sim().monitors[0]
    fm = modes.plane_monitor(mon)
    assert fm.colocate is True
    assert tuple(fm.center) == tuple(mon.center)
    assert tuple(fm.size) == tuple(mon.size)
    np.testing.assert_array_equal(np.asarray(fm.freqs), np.asarray(mon.freqs))


def test_tidy3d_itself_rejects_volume_mode_monitor():
    """plane_monitor does not check planarity, it relies on this guarantee. If it ever goes away we
    have to add the check ourselves."""
    mon = _sim().monitors[0]
    with pytest.raises(Exception, match="planar"):
        mon.updated_copy(size=(0.3, 1.2, 1.0))


def test_mode_solver_monitor_is_not_treated_as_mode_monitor():
    """A ``ModeSolverMonitor`` stores field profiles, not amplitudes, so it must not take the same
    path."""
    mon = _sim().monitors[0]
    assert modes.is_mode_monitor(mon)
    solver_mon = td.ModeSolverMonitor(
        center=mon.center, size=mon.size, freqs=mon.freqs,
        mode_spec=mon.mode_spec, name="ms")
    assert not modes.is_mode_monitor(solver_mon)


def _periodic_mode_sim():
    """A fully periodic transverse axis plus a single-plane mode monitor: Autograd20's geometry in
    miniature.

    y is periodic and the mode plane, size=(0, inf, inf), covers the whole transverse domain. A
    cell-centred component is one sample short at the upper y boundary (on a periodic domain
    F(hi)=F(lo)), while the colocation points include both domain edges. Without a periodic
    extension they are not covered and projection.field_data reports that the index box does not
    cover the colocation points.
    """
    return td.Simulation(
        size=(3.0, 0.5, 2.0), grid_spec=td.GridSpec.uniform(dl=0.025),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, 0), size=(td.inf, td.inf, 0.22)),
            medium=td.Medium(permittivity=12.1))],
        sources=[td.PointDipole(
            center=(-1.0, 0, 0), polarization="Ey",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.1 * FREQ0))],
        monitors=[td.ModeMonitor(
            center=(0.8, 0, 0), size=(0.0, td.inf, td.inf), freqs=[FREQ0],
            mode_spec=td.ModeSpec(num_modes=1), name="m")],
        run_time=2e-13,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.pml(num_layers=8), y=td.Boundary.periodic(),
            z=td.Boundary.pml(num_layers=8)))


def test_periodic_mode_plane_field_data_covers_colocation():
    """A full-domain mode plane on a periodic transverse axis: synthetic phasors must assemble onto
    the colocation points without a coverage complaint."""
    from openem.scene import projection

    sim = _periodic_mode_sim()
    sc = scene_mod.from_simulation(sim)
    mm = sim.monitors[0]
    pm = modes.plane_monitor(mm)
    spec = {m.name: m for m in sc.field_monitors}[pm.name]
    freqs = np.atleast_1d(np.asarray(mm.freqs, float))
    rng = np.random.RandomState(0)
    ph = (rng.randn(6, freqs.size, *spec.box)
          + 1j * rng.randn(6, freqs.size, *spec.box)).astype(np.complex128)
    fd = projection.field_data(sim, pm, spec, ph, scale=1.0)   # must not raise
    assert not np.isnan(np.asarray(fd.Ex.data)).any()
    # The colocated y axis must reach the domain edges at ±0.25, once the periodic extension fills
    # it in
    yc = np.asarray(fd.Ex.y)
    assert yc.max() >= 0.25 - 1e-6 and yc.min() <= -0.25 + 1e-6


def test_periodic_wrap_matches_at_domain_edges():
    """The periodic extension has to be physically right: F(y=-L/2) and F(y=+L/2) must be equal,
    since they are the same point on the torus."""
    from openem.scene import projection
    from tidy3d.components.data.data_array import ScalarFieldDataArray

    sim = _periodic_mode_sim()
    sc = scene_mod.from_simulation(sim)
    mm = sim.monitors[0]
    pm = modes.plane_monitor(mm)
    spec = {m.name: m for m in sc.field_monitors}[pm.name]
    freqs = np.atleast_1d(np.asarray(mm.freqs, float))
    period = 0.5
    co = projection._coords(sim, "Ey", spec.origin, spec.box)
    prof = np.exp(2j * np.pi * np.asarray(co["y"]) / period)   # strictly periodic
    ph = np.zeros((6, freqs.size, *spec.box), np.complex128)
    ph[1] = prof[None, None, :, None]
    fd = projection.field_data(sim, pm, spec, ph, scale=1.0)
    yt = np.asarray(fd.Ey.y)
    yax = list(fd.Ey.dims).index("y")
    lo = fd.Ey.data.take(int(np.argmin(abs(yt + 0.25))), axis=yax).ravel()[0]
    hi = fd.Ey.data.take(int(np.argmin(abs(yt - 0.25))), axis=yax).ravel()[0]
    assert abs(lo - hi) < 1e-9
