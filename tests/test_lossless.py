# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""``lossless`` has to rebuild the TF/SF incident tables.

If σ is merely cleared without rebuilding, the reference run injects into a lossless medium with an
H that was built for the lossy one, and the launched amplitude is off by a factor
``gain x |n_c|/n``. On BeerLambert that factor measures 1.00353 and goes straight into
``T = flux_lossy/flux_lossless``, so T comes out 0.64% low, disguised as "the optical path is half
a cell too long".

These tests are pure CPU, no GPU needed.
"""

from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import scene as scene_mod
from openem.model import lossless
from openem.scene import sources as src_mod

#: Use BeerLambert's medium and frequency directly. Rounding them to "nice" values, ε=9 / σ=4e-3,
#: trips the "σ drifts with frequency" fail-closed check in media.py (measured drift 8.3e-08 >
#: 1e-09), which is a different criterion and should not be mixed into this test.
FREQ0 = 1.9986163866666666e14
PERMITTIVITY = 8.996437927137574
CONDUCTIVITY = 0.003981628091989607          # S/µm


def _lossy_plane_wave_sim() -> "td.Simulation":
    """A uniform lossy medium plus a plane wave, the same class of scene as BeerLambert."""
    return td.Simulation(
        size=(0.1, 0.1, 1.2), grid_spec=td.GridSpec.uniform(dl=0.025),
        # There are no structures, so subpixel has nothing to do; but it would break the "σ drifts
        # with frequency" discriminator in media.eps_and_sigma (see that function's docstring), so
        # turning it off is cleaner.
        subpixel=False,
        medium=td.Medium(permittivity=PERMITTIVITY, conductivity=CONDUCTIVITY),
        sources=[td.PlaneWave(
            center=(0, 0, -0.4), size=(td.inf, td.inf, 0), direction="+",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.2 * FREQ0))],
        monitors=[td.FluxMonitor(center=(0, 0, 0.3), size=(td.inf, td.inf, 0),
                                 freqs=[FREQ0], name="T")],
        run_time=2e-13,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.periodic(), y=td.Boundary.periodic(),
            z=td.Boundary.pml(num_layers=12)))


def test_lossless_rebuilds_incident_table():
    sc = scene_mod.from_simulation(_lossy_plane_wave_sim())
    assert sc.sigma_ex is not None, "this scene has to be lossy, otherwise there is nothing to test"
    got = lossless(sc).sources[0]

    # It really must change: unchanged means only σ was cleared
    assert not np.array_equal(got.hy_inc, sc.sources[0].hy_inc), (
        "lossless() did not rebuild hy_inc")
    # And it must equal exactly what rebuilding the table with σ=0 gives
    kw = dict(sc.sources[0].incident)
    kw["sigma"] = 0.0
    _, hy_ref = src_mod._incident_tables(**kw)
    np.testing.assert_array_equal(got.hy_inc, hy_ref)
    # ex_inc does not depend on σ and must not move
    np.testing.assert_array_equal(got.ex_inc, sc.sources[0].ex_inc)


def test_lossless_fails_closed_without_incident_args():
    """Better to raise than to proceed without the original arguments: silently using the wrong
    table is impossible to track down."""
    from dataclasses import replace
    sc = scene_mod.from_simulation(_lossy_plane_wave_sim())
    broken = replace(sc, sources=[replace(sc.sources[0], incident=None)])
    with pytest.raises(ValueError, match="incident"):
        lossless(broken)


def test_lossless_is_noop_on_lossless_scene():
    """An already lossless scene must not be touched, which keeps it to a single code path."""
    sim = _lossy_plane_wave_sim().updated_copy(medium=td.Medium(permittivity=PERMITTIVITY))
    sc = scene_mod.from_simulation(sim)
    got = lossless(sc)
    np.testing.assert_array_equal(got.sources[0].hy_inc, sc.sources[0].hy_inc)
