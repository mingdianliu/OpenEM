# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""PointDipole / FieldMonitor / the general three-axis CPML.

The criteria come from Tidy3D's source and from StartHere's stored reference output, not from
expected values we computed ourselves.
"""

import numpy as np
import os
import pytest

td = pytest.importorskip("tidy3d")

from openem import scene as scene_mod
from openem import waveform
from openem.grid import Axis

TABLEA = os.environ.get("OPENEM_TABLEA", "/nonexistent")  # reference outputs; skip if unset
START = f"{TABLEA}/StartHere/Tidy3D/simulation.json"
BIO = f"{TABLEA}/BiosensorGrating/Tidy3D/simulation.json"


@pytest.fixture(scope="module")
def start():
    return scene_mod.from_file(START)


# ---------------------------------------------------------------- boundary topology

def test_index_tables_periodic_wraps():
    a = Axis(np.linspace(0.0, 4.0, 5), "Periodic", "Periodic")
    t = a.index_tables()
    assert list(t["nxt"]) == [1, 2, 3, 0]
    assert list(t["prv"]) == [3, 0, 1, 2]
    assert np.all(t["mnx"] == 1) and np.all(t["mpv"] == 1) and np.all(t["pec"] == 1)


def test_index_tables_absorbing_masks_outside():
    """An absorbing axis: the mask for out-of-domain neighbours is 0, and at i=0 the tangential E is
    zeroed by pec."""
    a = Axis(np.linspace(0.0, 4.0, 5), "PML", "PML")
    t = a.index_tables()
    assert list(t["nxt"]) == [1, 2, 3, 3]      # clamped at the high end
    assert list(t["prv"]) == [0, 0, 1, 2]      # clamped at the low end
    assert list(t["mnx"]) == [1, 1, 1, 0]      # outermost high cell has no i+1
    assert list(t["mpv"]) == [0, 1, 1, 1]      # outermost low cell has no i-1
    assert list(t["pec"]) == [0, 1, 1, 1]      # the PML is backed by PEC


def test_all_pml_scene_now_accepted(start):
    """PML on all six sides used to be run silently as periodic by the kernel; it must build
    normally now."""
    assert [a.boundary_lo for a in start.grid.axes] == ["PML"] * 3
    assert sorted(start.pml) == [(0, "hi"), (0, "lo"), (1, "hi"),
                                (1, "lo"), (2, "hi"), (2, "lo")]


def test_accepts_pec_and_pmc_walls():
    """A PEC wall uses the same masks as an absorbing end; a PMC wall is a magnetic mirror
    (tangential H odd, E even). Both are accepted.

    The low PMC wall is exact (mpv[0]=-1, reflection phase 180 degrees opposite to PEC); under
    drop-last the high wall is a first-order approximation (|R|=1, phase off by about 6 degrees),
    see test_pml_reflection.py.
    """
    sim = td.Simulation.from_file(START)
    for bnd, name in ((td.Boundary.pec(), "PECBoundary"),
                      (td.Boundary.pmc(), "PMCBoundary")):
        sc = scene_mod.from_simulation(sim.updated_copy(
            boundary_spec=sim.boundary_spec.updated_copy(x=bnd)))
        assert sc.grid.axes[0].boundary_lo == name


def test_accepts_pml_and_absorber_on_same_axis():
    """Both ends absorbing is enough; the type names need not match, since the index tables and
    masks are the same."""
    sim = td.Simulation.from_file(START)
    ok = sim.updated_copy(
        boundary_spec=sim.boundary_spec.updated_copy(
            x=td.Boundary(minus=td.PML(), plus=td.Absorber())
        )
    )
    sc = scene_mod.from_simulation(ok)
    assert sc.grid.axes[0].boundary_lo == "PML"
    assert sc.grid.axes[0].boundary_hi == "Absorber"


# ---------------------------------------------------------------- PointDipole

def test_dipole_stencil_weights_sum_to_one(start):
    """The multilinear weights must sum to 1, otherwise the current moment is not the one Tidy3D
    declares."""
    d = start.dipoles[0]
    vol = start.grid.yee_dual_volumes()[d.component]
    w = np.array([c * float(vol[i, j, k]) / scene_mod.UM
                  for (i, j, k), c in zip(d.indices, d.coef)])
    np.testing.assert_allclose(w.sum(), 1.0, rtol=1e-12)
    assert np.all(w > 0)


def test_dipole_coef_carries_um_to_m_conversion(start):
    """``coef = w * UM / dV``: amplitude=1 is a current moment of 1 A·µm, in SI 1e-6 A·m.

    Missing this 1e-6 makes the absolute field 1e6 times too large, which is how it was found.
    """
    d = start.dipoles[0]
    vol = start.grid.yee_dual_volumes()[d.component]
    for (i, j, k), c in zip(d.indices, d.coef):
        w = c * float(vol[i, j, k]) / scene_mod.UM
        assert 0.0 < w <= 1.0, "the recovered weight must lie in (0, 1]"


def test_dipole_component_matches_polarization(start):
    sim = td.Simulation.from_file(START)
    assert str(sim.sources[0].polarization) == "Ey"
    assert start.dipoles[0].component == 1


def test_magnetic_dipole_now_accepted():
    """Magnetic dipoles are supported (Bandstructure's Hz): they take the H injection path, see
    test_dipole_h.py."""
    sim = td.Simulation.from_file(START)
    sc = scene_mod.from_simulation(sim.updated_copy(
        sources=[sim.sources[0].updated_copy(polarization="Hy")]
    ))
    assert sc.dipoles[0].magnetic and sc.dipoles[0].component == 1


def test_no_plane_wave_source_is_fine(start):
    """A pure dipole scene has no plane-wave source and no FluxMonitor."""
    assert start.sources == []
    assert start.flux_monitors == []
    assert len(start.dipoles) == 1


# ---------------------------------------------------------------- FieldMonitor

def test_field_monitor_box_covers_cells_needed_for_colocation(start):
    """The box has to cover the cells colocation needs, without hard-coding a cell count.

    The reference output has a single z point, at z=0, which is a **cell centre** (147 is odd). The
    components sitting at the cell centre (Ez/Hx/Hy) only need ctr_z[73]; those sitting on the cell
    boundary (Ex/Ey/Hz) need edges[73] and edges[74], that is, cells 73 and 74.
    """
    m = start.field_monitors[0]
    assert m.box[:2] == start.shape[:2]        # x/y are inf, so they cover the whole axis
    lo, n = m.origin[2], m.box[2]
    assert lo <= 73 and lo + n >= 75, f"origin={lo} box={n} does not cover cells 73/74"


def test_archive_coords_are_our_primal_boundaries():
    """The colocated coordinates in the reference output are exactly the primal cell boundaries,
    matching point by point on x and y."""
    import h5py, json
    sc = scene_mod.from_file(START)
    with h5py.File(f"{TABLEA}/StartHere/Tidy3D/data.hdf5", "r") as f:
        raw = f["JSON_STRING"][()]
        js = json.loads(raw.decode() if isinstance(raw, bytes) else raw)
        i = next(n for n, it in enumerate(js["data"]) if it.get("type") == "FieldData")
        cx = np.asarray(f[f"data/{i}/Ex/x"], dtype=np.float64)
        cy = np.asarray(f[f"data/{i}/Ex/y"], dtype=np.float64)
    np.testing.assert_allclose(cx, sc.grid.axes[0].edges / scene_mod.UM, atol=1e-12)
    np.testing.assert_allclose(cy, sc.grid.axes[1].edges / scene_mod.UM, atol=1e-12)


# ---------------------------------------------------------------- normalization convention

def test_spectrum_is_dt_dft_over_sqrt_two_pi():
    """Tidy3D's ``spectrum`` carries ``dt/√(2π)`` (``components/time.py:112``).

    The absolute-amplitude comparison for StartHere relies on this to derive the constant rather
    than fit it, so the convention has to be locked: change it and that parameter-free criterion
    no longer holds.
    """
    sim = td.Simulation.from_file(START)
    st = sim.sources[0].source_time
    freqs = np.array([4.0e14])
    times = np.arange(sim.num_time_steps + 1, dtype=np.float64) * sim.dt
    S = np.asarray(st.spectrum(times, freqs, sim.dt))[0]

    amps = np.real(np.asarray(st.amp_time(times)))
    dft = np.sum(amps * np.exp(2j * np.pi * freqs[0] * times))
    np.testing.assert_allclose(S, sim.dt * dft / np.sqrt(2 * np.pi), rtol=1e-9)


# ---------------------------------------------------------------- the general kernel

def test_all_pml_scene_actually_absorbs(start):
    """PML on all six sides must actually absorb.

    The old xyper kernel treated x/y as periodic, so energy wrapped around and after 2000 steps it
    only decayed to 7.7e-2. This criterion exists so that cannot happen silently again.
    """
    pytest.importorskip("cupy")
    from openem import solver
    res = solver.run(start, num_steps=2400, use_shutoff=True, verbose=False)
    ratio = res.decay_trace[-1][1]
    assert ratio < 1e-2, f"PML on six sides but decay only reached {ratio:.3e}; is x/y absorbing?"
