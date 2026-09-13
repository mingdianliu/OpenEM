# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Non-finite values are refused at the assembly entry point (nb.backend._refuse_nonfinite_ours).

Background: a parameter-sweep batch once produced NaN mode amplitudes inside a solve job. The NaN
passed through np.angle, np.unwrap and np.interp and only blew up four cells later as
td.Box(size=nan). The solve job's stdout was lost along with a failed notebook conversion, so the
print sentinel at the solver exit could not point at the source. Assembly and scene building were
then shown to be innocent (11 of 11 clean on synthetic finite data), so the refusal moved forward
to the build entry: when the data is bad, name the monitor on the spot and fail closed.

The checks run entirely on CPU and never touch the GPU:
* finite data: build returns a SimulationData as before, bitwise unchanged;
* data with NaN: RuntimeError whose message names fld:xxx;
* OPENEM_NB_ALLOW_NONFINITE=1: downgraded to a warning, assembly continues.
"""

from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem.nb import backend as nbmod
from openem import scene as scene_mod

FREQ0 = 2.0e14


def _tiny_sim() -> "td.Simulation":
    return td.Simulation(
        size=(1.0, 1.0, 1.0), grid_spec=td.GridSpec.uniform(dl=0.1),
        structures=[],
        sources=[td.PointDipole(
            center=(0, 0, 0), polarization="Ez",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.2 * FREQ0))],
        monitors=[td.FieldMonitor(center=(0, 0, 0.2), size=(0.6, 0.6, 0),
                                  freqs=[FREQ0], name="field")],
        run_time=1e-13,
        boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=6)))


@pytest.fixture(scope="module")
def scene_and_ours():
    sim = _tiny_sim()
    sc = scene_mod.from_simulation(sim)
    spec = {m.name: m for m in sc.field_monitors}["field"]
    rng = np.random.default_rng(7)
    ph = (rng.standard_normal((6, 1, *spec.box))
          + 1j * rng.standard_normal((6, 1, *spec.box))).astype(np.complex128)
    return sim, sc, ph


def test_finite_data_builds(scene_and_ours):
    sim, sc, ph = scene_and_ours
    sd = nbmod.build(sim, sc, {"fld:field": ph})
    arr = np.asarray(sd["field"].Ez.values)
    assert arr.size and np.all(np.isfinite(arr))


def test_nan_data_fails_closed(scene_and_ours):
    sim, sc, ph = scene_and_ours
    bad = ph.copy()
    bad[2, 0, 1, 1, 0] = np.nan + 1j * np.nan
    with pytest.raises(RuntimeError, match=r"fld:field"):
        nbmod.build(sim, sc, {"fld:field": bad})


def test_env_override_downgrades(scene_and_ours, monkeypatch, capsys):
    sim, sc, ph = scene_and_ours
    bad = ph.copy()
    bad[0, 0, 0, 0, 0] = np.inf
    monkeypatch.setenv("OPENEM_NB_ALLOW_NONFINITE", "1")
    sd = nbmod.build(sim, sc, {"fld:field": bad})
    assert sd is not None
    out = capsys.readouterr().out
    assert "fld:field" in out and "WARNING" in out
