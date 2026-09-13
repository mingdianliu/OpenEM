# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""``CustomFieldSource`` -> equivalent surface currents (J = n̂×H, M = -n̂×E, Love's equivalence
principle).

The second stage of Autograd29SourceGradients uses a ``CustomFieldSource`` as its forward source,
injecting the tangential E/H fields that the first stage projected onto a plane. OpenEM did not
recognise this type before. What is locked here is the conversion: the component mapping, the
signs, and the dispatch that sends electric currents to E* and magnetic currents to H*. These tests
only build the scene, they do not run a simulation.
"""

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import scene as scene_mod
from openem.scene import sources as S

FREQ = td.C_0 / 0.94
N = 6
_ax = np.linspace(-0.3, 0.3, N)
_st = td.GaussianPulse(freq0=FREQ, fwidth=FREQ / 10)


def _sfa(v):
    coords = dict(x=_ax, y=_ax, z=np.array([0.0]), f=np.array([FREQ]))
    return td.ScalarFieldDataArray(np.full((N, N, 1, 1), v, dtype=complex),
                                   coords=coords)


def _cfs(**field_components):
    ds = td.FieldDataset(**field_components)
    return td.CustomFieldSource(center=(0, 0, 0), size=(0.4, 0.4, 0),
                                source_time=_st, field_dataset=ds)


def _sim(src):
    return td.Simulation(
        size=(0.8, 0.8, 1.2), run_time=2e-14,
        grid_spec=td.GridSpec.uniform(dl=0.05), sources=[src],
        monitors=[td.FieldMonitor(center=(0, 0, 0.4), size=(td.inf, td.inf, 0),
                                  freqs=[FREQ], name="f")],
        boundary_spec=td.BoundarySpec.all_sides(td.PML()))


def test_custom_field_source_builds_scene():
    """Given tangential E and H the scene builds, landing as a set of dipoles."""
    sc = scene_mod.from_simulation(_sim(_cfs(Ex=_sfa(1.0), Ey=_sfa(0.5),
                                             Hx=_sfa(0.2), Hy=_sfa(1.0))))
    assert len(sc.dipoles) > 0
    assert not sc.sources and not sc.mode_sources


def test_equivalent_currents_component_mapping():
    """Injection along z: J = n̂×H gives Jx=-Hy, Jy=Hx; M = -n̂×E gives Mx=Ey, My=-Ex.

    With only Ex and Hy given, the result should be one electric-current component (Jx, from -Hy,
    not magnetic) and one magnetic-current component (My, from -Ex, magnetic). OpenEM puts the
    electric current into a dipole with component=x and the magnetic current into one with
    component=y.
    """
    from openem.scene import boundaries
    sim = _sim(_cfs(Ex=_sfa(2.0), Hy=_sfa(3.0)))
    g = boundaries.build_grid(sim)
    dips = S.custom_field(sim, sim.sources[0], g)
    j = [d for d in dips if not d.magnetic]      # equivalent electric current
    m = [d for d in dips if d.magnetic]          # equivalent magnetic current
    assert j and m, "both the equivalent electric and magnetic currents must be present"
    assert all(d.component == 0 for d in j), "Jx = -Hy must land on the x component"
    assert all(d.component == 1 for d in m), "My = -Ex must land on the y component"
    assert np.any(np.concatenate([d.coef for d in j]) != 0)
    assert np.any(np.concatenate([d.coef for d in m]) != 0)
