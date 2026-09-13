# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Injection from several independent ModeSource planes, the shape of the S-matrix adjoint sim in
Autograd27.

The adjoint sim turns every ModeMonitor the objective touches into a ModeSource, on different
planes, with different time waveforms and all amplitudes nonzero, injected at the same time. The
fix: ``_mode_source`` returns a list of injection entries and the solver calls the accumulating
injection kernel once per entry. The Yee update and the TF/SF correction are linear, so the
multi-source solution is the sum of the single-source solutions. The criterion locks exactly that
**superposition** (two sources run together == two single-source runs added), plus the grouping
done while building the tables: different planes no longer raise NotImplementedError, and the same
plane with the same waveform still merges into one entry, leaving the old BraggGratings path
bitwise unchanged.

The scene is deliberately small (40x32x28 cells, a few hundred steps), a few seconds on a GPU.
"""

from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

import pytest

# This module imports openem.solver at the top level, which imports cupy at its own top level: on a
# CPU-only machine that would crash during collection.
pytest.importorskip("cupy", reason="needs CUDA and cupy")

from openem import scene as scene_mod
from openem.setup_tables import _mode_source

FREQ0 = 1.934e14


def _wg_sim(sources) -> "td.Simulation":
    """A silicon strip waveguide along x, with PML all around."""
    return td.Simulation(
        size=(2.0, 1.6, 1.4), grid_spec=td.GridSpec.uniform(dl=0.05),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, 0), size=(td.inf, 0.5, 0.22)),
            medium=td.Medium(permittivity=12.1))],
        sources=sources,
        monitors=[],
        run_time=2e-13,
        boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=8)))


def _mode_src(x: float, direction: str, freq0: float = FREQ0) -> "td.ModeSource":
    return td.ModeSource(
        center=(x, 0, 0), size=(0, 1.2, 1.0),
        source_time=td.GaussianPulse(freq0=freq0, fwidth=0.1 * freq0),
        direction=direction, mode_spec=td.ModeSpec(num_modes=1), mode_index=0)


def _fields(sim, num_steps):
    from openem import solver
    from openem.device import Kernels

    sc = scene_mod.from_simulation(sim)
    res = solver.run(sc, num_steps=num_steps, use_shutoff=False,
                     kernels=Kernels(), verbose=False, return_fields=True)
    return res.fields


def test_mode_source_entries_split_by_plane():
    """Different planes and different waveforms: one entry per source, no NotImplementedError."""
    src_a = _mode_src(-0.6, "+")
    src_b = _mode_src(+0.6, "-", freq0=1.05 * FREQ0)   # plane, direction and waveform all differ
    sc = scene_mod.from_simulation(_wg_sim([src_a, src_b]))
    assert len(sc.mode_sources) == 2
    entries = _mode_source(sc)
    assert isinstance(entries, list) and len(entries) == 2
    ks = {(int(e["axis"]), int(e["ks"])) for e in entries}
    assert len(ks) == 2, "the two entries must land on different injection planes"


def test_mode_source_same_plane_still_merges():
    """The same plane with the same waveform (the BraggGratings shape) still merges into one entry,
    so the old path is unchanged."""
    wg = dict(size=(td.inf, 0.3, 0.22), medium=td.Medium(permittivity=12.1))
    sim = td.Simulation(
        size=(2.0, 2.4, 1.4), grid_spec=td.GridSpec.uniform(dl=0.05),
        structures=[
            td.Structure(geometry=td.Box(center=(0, -0.6, 0), size=wg["size"]),
                         medium=wg["medium"]),
            td.Structure(geometry=td.Box(center=(0, +0.6, 0), size=wg["size"]),
                         medium=wg["medium"])],
        sources=[
            td.ModeSource(center=(-0.6, -0.6, 0), size=(0, 0.9, 1.0),
                          source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.1 * FREQ0),
                          direction="+", mode_spec=td.ModeSpec(num_modes=1), mode_index=0),
            td.ModeSource(center=(-0.6, +0.6, 0), size=(0, 0.9, 1.0),
                          source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.1 * FREQ0),
                          direction="+", mode_spec=td.ModeSpec(num_modes=1), mode_index=0)],
        monitors=[], run_time=2e-13,
        boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=8)))
    sc = scene_mod.from_simulation(sim)
    assert len(sc.mode_sources) == 2
    entries = _mode_source(sc)
    assert isinstance(entries, list) and len(entries) == 1, (
        "same plane, same waveform must merge")


def test_two_plane_injection_superposes():
    """Superposition: the fields from running both sources together == the two single-source runs
    added (linear, within float32 tolerance)."""
    src_a = _mode_src(-0.6, "+")
    src_b = _mode_src(+0.6, "-")
    n_steps = 400

    f_a = _fields(_wg_sim([src_a]), n_steps)
    f_b = _fields(_wg_sim([src_b]), n_steps)
    f_ab = _fields(_wg_sim([src_a, src_b]), n_steps)

    scale = max(float(np.abs(f_a["Ey"]).max()), float(np.abs(f_b["Ey"]).max()))
    assert scale > 0, "a single-source field must not be zero: was the mode injected at all?"
    worst = 0.0
    for c in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz"):
        worst = max(worst, float(np.abs(f_ab[c] - (f_a[c] + f_b[c])).max()))
    rel = worst / scale
    assert rel < 1e-4, f"superposition broken: max|F_ab-(F_a+F_b)|/peak = {rel:.3e}"
