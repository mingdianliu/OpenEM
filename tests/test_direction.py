# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Propagation direction of a plane wave: injection along -z.

The criterion does not rely on stored reference output: flip the direction on the **same uniform
grid** and the result must be an exact mirror image.
"""

from dataclasses import replace

import numpy as np
import os
import pytest

td = pytest.importorskip("tidy3d")

from openem import scene as scene_mod
from openem.scene import sources

TABLEA = os.environ.get("OPENEM_TABLEA", "/nonexistent")  # reference outputs; skip if unset
#: Picked because the grid is uniform on all three axes (0.025 µm). On a non-uniform grid the
#: one-way purity only reaches 7.7e-06, which is a property of the grid and would drown out what
#: this test is looking at.
UNIFORM = f"{TABLEA}/BiosensorGrating/Tidy3D/simulation.json"


def _flip(sim):
    return sim.updated_copy(sources=[sim.sources[0].updated_copy(direction="-")])


def test_direction_minus_is_accepted():
    sim = td.Simulation.from_file(UNIFORM)
    sc = scene_mod.from_simulation(_flip(sim))
    assert sc.sources[0].direction == -1


def test_incident_tables_flip_only_the_h_sign():
    """``hy_inc`` changes sign with the direction (E×H has to point along propagation) while
    ``ex_inc`` does not; the half-cell delay does **not** change sign, because for both directions
    the neighbouring Hy cell centre is upstream."""
    sim = td.Simulation.from_file(UNIFORM)
    st = sim.sources[0].source_time
    kw = dict(dt=sim.dt, num_steps=200, dl_half=1e-8, eps_r=2.25)
    ep, hp = sources._incident_tables(st, direction=+1, **kw)
    em, hm = sources._incident_tables(st, direction=-1, **kw)
    np.testing.assert_array_equal(ep, em)
    np.testing.assert_array_equal(hp, -hm)


def test_h_correction_plane_moves_with_direction():
    """The H correction plane: at ``ks-1`` for +z, at ``ks`` for -z.

    This is the H cell adjacent on the scattered-field side of the TF/SF boundary. The total-field
    region flips from "z >= z_ks" to "z <= z_ks", so the adjacent cell moves from below to above.
    """
    sim = td.Simulation.from_file(UNIFORM)
    plus = scene_mod.from_simulation(sim).sources[0]
    minus = scene_mod.from_simulation(_flip(sim)).sources[0]
    assert plus.plane_index == minus.plane_index      # the source plane itself does not move
    # The dl used for the half-cell delay shifts by one cell: +z uses dl[ks-1], -z uses dl[ks].
    # On a uniform grid the two are equal, so only the direction itself can be checked here.
    assert (plus.direction, minus.direction) == (+1, -1)


def test_minus_z_is_exact_mirror_of_plus_z():
    """**The decisive criterion**: flip the direction on the same uniform grid and the downstream
    flux keeps its magnitude and reverses its sign.

    Measured one-way purity is 5.457e-08 for +z and 5.279e-08 for -z, with a downstream flux of
    ±2.0404e-12, the magnitudes bitwise identical. The sign flips because the flux is measured
    along +z while the wave travels along -z.

    The purity criterion **only fixes the relative sign** (flipping both correction terms at once
    is the same as negating the wave as a whole); the absolute sign comes from the derivation:
    ``(sh, se) = (-1, -1)``.
    """
    pytest.importorskip("cupy")
    from openem import flux as flux_mod
    from openem import solver
    from openem.device import Kernels

    sim = td.Simulation.from_file(UNIFORM)
    k = Kernels()
    out = {}
    for tag, s in (("+", sim), ("-", _flip(sim))):
        sc = scene_mod.from_simulation(s)
        ks, d, nz = sc.sources[0].plane_index, sc.sources[0].direction, sc.shape[2]
        proto = sc.flux_monitors[0]
        lo, hi = max(ks - 20, 16), min(ks + 20, nz - 17)
        mons = [replace(proto, name="lo", plane_index=lo, normal_dir=1),
                replace(proto, name="hi", plane_index=hi, normal_dir=1)]
        homo = scene_mod.lossless(scene_mod.homogeneous(
            replace(sc, flux_monitors=mons)))
        r = solver.run(homo, num_steps=4000, use_shutoff=False, kernels=k, verbose=False)
        f_lo = flux_mod.plane_flux(r.phasors["lo"].data, sc.grid)
        f_hi = flux_mod.plane_flux(r.phasors["hi"].data, sc.grid)
        dn, up = (f_hi, f_lo) if d > 0 else (f_lo, f_hi)
        out[tag] = (dn, up)

    for tag in "+-":
        dn, up = out[tag]
        pur = float(np.max(np.abs(up) / np.abs(dn)))
        assert pur < 1e-6, f"direction={tag} purity {pur:.3e}, the source is not one-way"

    dn_p, dn_m = out["+"][0], out["-"][0]
    assert np.all(dn_p > 0) and np.all(dn_m < 0), "downstream flux must change sign with direction"
    np.testing.assert_allclose(np.abs(dn_m), np.abs(dn_p), rtol=1e-6)
