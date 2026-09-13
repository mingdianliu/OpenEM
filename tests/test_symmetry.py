# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Mirror sources across a symmetry plane.

The sign table is checked directly against Tidy3D's ``field_symmetry_eigenvalue``, so if that
changes, this goes red.
"""

import numpy as np
import os
import pytest

td = pytest.importorskip("tidy3d")
from tidy3d.components.data.em_fields import field_symmetry_eigenvalue

from openem import scene as scene_mod
from openem.scene import sources

TABLEA = os.environ.get("OPENEM_TABLEA", "/nonexistent")  # external reference archive; without
#                                                          OPENEM_TABLEA these cases skip
RF = f"{TABLEA}/ResonanceFinder/Tidy3D/simulation.json"
START = f"{TABLEA}/StartHere/Tidy3D/simulation.json"


def _weights(sc, d):
    """Recover the multilinear weights, mirror signs included, from coef."""
    vol = sc.grid.yee_dual_volumes()[d.component]
    return np.array([c * float(vol[i, j, k]) / scene_mod.UM
                     for (i, j, k), c in zip(d.indices, d.coef)])


def test_reflection_sign_matches_tidy3d_table():
    """The sign is ``symmetry[dim] * field_symmetry_eigenvalue`` (``em_fields.py:20-25`` and
    ``monitor_data.py:450-452``), checked term by term against Tidy3D's own function.
    """
    for comp in range(3):
        for dim in range(3):
            for sym in (1, -1):
                want = sym * field_symmetry_eigenvalue("E" + "xyz"[comp], dim)
                assert sources.reflection_sign(comp, dim, sym) == want


def test_reflection_sign_physics():
    """Beside a PEC plane (symmetry=-1), the image of a **tangential** dipole takes -1, since
    tangential E vanishes on the plane, and a **normal** one takes +1.

    Two physical self-checks that do not depend on the source code.
    """
    # The plane normal is x
    assert sources.reflection_sign(1, 0, -1) == -1        # tangential (y polarized)
    assert sources.reflection_sign(2, 0, -1) == -1        # tangential (z polarized)
    assert sources.reflection_sign(0, 0, -1) == +1        # normal (x polarized)


def test_no_symmetry_leaves_stencil_alone():
    """A scene with no symmetry adds no images at all."""
    sc = scene_mod.from_file(START)
    w = _weights(sc, sc.dipoles[0])
    assert w.size == 4
    np.testing.assert_allclose(w.sum(), 1.0, rtol=1e-12)


def test_generic_position_gets_distinct_images():
    """symmetry=(1,-1,1), with the dipole at a general position in x and y and sitting on the z plane.

    x and y each generate one image at a different position while the z image lands back on the
    original, so the total weight is 2^3 = 8.
    """
    sc = scene_mod.from_file(RF)
    for d in sc.dipoles:
        w = _weights(sc, d)
        assert w.size == 32, "8 images times a 4-point stencil"
        np.testing.assert_allclose(w.sum(), 8.0, rtol=1e-12)
        assert np.all(w > 0), "with Ey and symmetry=(1,-1,1) every mirror sign should be +1"


def test_on_plane_source_is_doubled():
    """A source sitting on a symmetry plane is doubled by its own image.

    Measured: skipping the z image made the field come out exactly 0.5 times Tidy3D's (the ratio
    across 7 monitors ran 0.4988 to 0.5011), and adding it back brought that to about 1.00.
    """
    sim = td.Simulation.from_file(RF)
    ctr = tuple(float(v) for v in sim.center)
    ext = tuple(abs(float(v)) for v in sim.size)
    # Ey, with symmetry only in z and the source sitting exactly on the z plane
    imgs = sources.mirror_images((0.3, 0.4, ctr[2]), 1, (0, 0, 1), ctr, ext)
    assert len(imgs) == 2
    assert all(abs(p[2] - ctr[2]) < 1e-15 for p, _ in imgs), "both images are on the plane"
    assert sum(w for _, w in imgs) == 2.0


def test_on_plane_with_negative_sign_fails_closed():
    """Sitting on the plane while requiring odd symmetry means the source must be identically zero:
    contradictory input, so it raises rather than silently producing a zero field.
    """
    sim = td.Simulation.from_file(RF)
    ctr = tuple(float(v) for v in sim.center)
    ext = tuple(abs(float(v)) for v in sim.size)
    with pytest.raises(NotImplementedError, match="symmetry plane"):
        # Ey sits on the y plane with symmetry[y]=+1, so the sign is (+1)(-1) = -1
        sources.mirror_images((0.3, ctr[1], 0.5), 1, (0, 1, 0), ctr, ext)


def test_images_stay_inside_domain():
    """An image must still lie inside the domain, or _linear_weights raises about a point outside
    the coordinate range.
    """
    sc = scene_mod.from_file(RF)
    for ax in range(3):
        e = sc.grid.axes[ax].edges
        for d in sc.dipoles:
            assert d.indices[:, ax].min() >= 0
            assert d.indices[:, ax].max() < sc.grid.axes[ax].n
