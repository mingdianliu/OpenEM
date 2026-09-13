# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Plane wave with ``pol_angle=+-pi/2`` (Ey polarization): equivalence against Ex polarization.

No free parameters. In isotropic vacuum on the same uniform grid, switching the polarization from
Ex to Ey must leave both the directionality and the downstream flux unchanged, because the
discretization is symmetric between the two transverse axes.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import os
import pytest

td = pytest.importorskip("tidy3d")

from openem import scene as scene_mod

TABLEA = os.environ.get("OPENEM_TABLEA", "/nonexistent")  # external reference archive; without
#                                                          OPENEM_TABLEA these cases skip
#: The same example as test_direction.py: uniform 0.025 um grid on all three axes, so the
#: directionality check is not drowned by the 7.7e-06 floor a non-uniform grid imposes.
UNIFORM = f"{TABLEA}/BiosensorGrating/Tidy3D/simulation.json"


def _pol(sim, angle):
    return sim.updated_copy(sources=[sim.sources[0].updated_copy(pol_angle=angle)])


def test_pol_angle_maps_to_second_transverse_axis():
    sim = td.Simulation.from_file(UNIFORM)
    for angle, axis in ((0.0, 0), (np.pi / 2, 1), (-np.pi / 2, 1)):
        sc = scene_mod.from_simulation(_pol(sim, angle))
        assert sc.sources[0].pol_axis == axis, (angle, axis)


def test_other_pol_angles_fail_closed():
    sim = td.Simulation.from_file(UNIFORM)
    with pytest.raises(NotImplementedError):
        scene_mod.from_simulation(_pol(sim, 0.3))


def test_ey_incident_tables_reused():
    """Ey polarization reuses the same pair of incident tables; the sign difference is folded
    entirely into the kernel coefficients.
    """
    sim = td.Simulation.from_file(UNIFORM)
    sx = scene_mod.from_simulation(_pol(sim, 0.0)).sources[0]
    sy = scene_mod.from_simulation(_pol(sim, np.pi / 2)).sources[0]
    np.testing.assert_array_equal(sx.ex_inc, sy.ex_inc)
    np.testing.assert_array_equal(sx.hy_inc, sy.hy_inc)


def test_ey_matches_ex_in_vacuum():
    """The decisive check: in uniform vacuum, Ex and Ey polarization give the same directionality
    and the same downstream flux.

    Measured over 4000 steps with one flux plane upstream and one downstream: the directionality
    of both polarizations is of order 1e-07 and the downstream fluxes agree to better than 1e-6
    relative. The discretization is symmetric between the two transverse axes, and the sign
    derivation introduces no asymmetric injection error.
    """
    pytest.importorskip("cupy")
    from openem import flux as flux_mod
    from openem import solver
    from openem.device import Kernels

    sim = td.Simulation.from_file(UNIFORM)
    k = Kernels()
    down = {}
    for tag, angle in (("Ex", 0.0), ("Ey", np.pi / 2)):
        sc = scene_mod.from_simulation(_pol(sim, angle))
        ks, d, nz = sc.sources[0].plane_index, sc.sources[0].direction, sc.shape[2]
        proto = sc.flux_monitors[0]
        lo, hi = max(ks - 20, 16), min(ks + 20, nz - 17)
        mons = [replace(proto, name="lo", plane_index=lo, normal_dir=1),
                replace(proto, name="hi", plane_index=hi, normal_dir=1)]
        homo = scene_mod.lossless(scene_mod.homogeneous(
            replace(sc, flux_monitors=mons)))
        r = solver.run(homo, num_steps=4000, use_shutoff=False, kernels=k,
                       verbose=False)
        f_lo = flux_mod.plane_flux(r.phasors["lo"].data, sc.grid)
        f_hi = flux_mod.plane_flux(r.phasors["hi"].data, sc.grid)
        dn, up = (f_hi, f_lo) if d > 0 else (f_lo, f_hi)
        pur = float(np.max(np.abs(up) / np.abs(dn)))
        assert pur < 1e-6, f"{tag} polarization directionality {pur:.3e}"
        down[tag] = dn

    np.testing.assert_allclose(down["Ey"], down["Ex"], rtol=1e-6)
