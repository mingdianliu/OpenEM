# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""ModeSource with a y or z normal: serialization round-trip and backward compatibility.

The physics of the mode profile and the injection is verified elsewhere, where directionality
along all three axes came out identical to the last bit. This file only guards the serialization
contract.
"""
from __future__ import annotations

import dataclasses

import numpy as np
import tidy3d as td

from openem import scene as scene_mod
from openem import serialize
from openem.model import ModeSource

FREQ0 = 2e14


def _sc():
    sim = td.Simulation(
        size=(1.0, 1.0, 1.0), grid_spec=td.GridSpec.uniform(dl=0.1),
        run_time=1e-14, subpixel=False,
        sources=[td.PointDipole(center=(0, 0, 0), polarization="Ez",
                                source_time=td.GaussianPulse(
                                    freq0=FREQ0, fwidth=FREQ0 / 10))],
        monitors=[td.FieldMonitor(center=(0, 0, 0), size=(0, 0, 0),
                                  freqs=[FREQ0], name="f")],
        boundary_spec=td.BoundarySpec.all_sides(td.Periodic()))
    return scene_mod.from_simulation(sim)


def _msrc(axis: int) -> ModeSource:
    r = np.random.default_rng(7)
    f = lambda: (r.random((8, 6)) + 1j * r.random((8, 6)))
    t = np.arange(12, dtype=np.float64)
    return ModeSource(plane_index=3, direction=-1, axis=axis,
                      ey_inc=f(), ez_inc=f(), hy_inc=f(), hz_inc=f(),
                      amp_e=np.exp(1j * t), amp_h=np.exp(1j * (t + 0.5)),
                      n_eff=2.4 + 1e-4j)


def test_serialize_roundtrip_mode_source_axis(tmp_path):
    sc = dataclasses.replace(_sc(), mode_sources=[_msrc(axis=1)])
    back = serialize.load(serialize.save(sc, tmp_path / "s.npz"))
    b, m = back.mode_sources[0], sc.mode_sources[0]
    assert (b.axis, b.plane_index, b.direction) == (1, 3, -1)
    for f in ("ey_inc", "ez_inc", "hy_inc", "hz_inc", "amp_e", "amp_h"):
        np.testing.assert_array_equal(getattr(b, f), getattr(m, f))
    assert b.n_eff == m.n_eff


def test_serialize_old_two_field_meta_means_x_axis(tmp_path):
    """An npz saved before the change, whose meta holds only [plane_index, direction], reads back
    as axis=0.
    """
    sc = dataclasses.replace(_sc(), mode_sources=[_msrc(axis=0)])
    path = serialize.save(sc, tmp_path / "s.npz")
    d = dict(np.load(path, allow_pickle=False))
    d["msrc0_meta"] = d["msrc0_meta"][:2]
    np.savez(path, **d)
    back = serialize.load(path)
    assert back.mode_sources[0].axis == 0
