# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""The grid and the ghost count.

The criterion comes from BiosensorGrating's log from the reference implementation:

    Simulation domain Nx, Ny, Nz: [22, 16, 140]
    Number of computational grid points: 4.9984e+04

22 x 16 x (140+2) = 49,984, so the ghost rule holds here independently.
"""

import numpy as np
import pytest

from openem.grid import Axis, Grid, ghost_count, nearest_index


def _uniform(n, dl=2.5e-8, lo="PML", hi="PML"):
    return Axis(np.arange(n + 1) * dl, lo, hi)


def test_ghost_by_boundary():
    assert ghost_count("PML") == 1
    assert ghost_count("StablePML") == 1
    assert ghost_count("Absorber") == 1
    assert ghost_count("Periodic") == 0
    assert ghost_count("BlochBoundary") == 0
    # A periodic axis that also carries symmetry still needs one ghost layer
    assert ghost_count("Periodic", has_symmetry=True) == 1
    # An end that is both PML and a symmetry plane only gets the layer once
    assert ghost_count("PML", has_symmetry=True) == 1


def test_unknown_boundary_refuses():
    """An unknown boundary type has to raise. Missing a ghost layer gives no symptom at all, it
    just shifts the cell count convention."""
    with pytest.raises(KeyError, match="unknown boundary type"):
        ghost_count("SomethingNew")


def test_biosensorgrating_cloud_cell_count():
    """Reproduce the 4.9984e+04 from the reference implementation's log."""
    grid = Grid(
        x=_uniform(22, lo="Periodic", hi="Periodic"),
        y=_uniform(16, lo="Periodic", hi="Periodic"),
        z=_uniform(140, lo="PML", hi="PML"),
    )
    assert grid.shape == (22, 16, 140)
    assert grid.num_cells == 49_280
    assert grid.shape_with_ghost == (22, 16, 142)
    assert grid.num_cells_cloud == 49_984


def test_periodic_axis_has_no_ghost():
    ax = _uniform(22, lo="Periodic", hi="Periodic")
    assert ax.ghost_lo == 0 and ax.ghost_hi == 0
    assert ax.n_with_ghost == 22
    assert ax.is_periodic


def test_primal_and_dual_spacing_uniform():
    """On a uniform grid the dual spacing equals the primal spacing."""
    dl = 2.5e-8
    ax = _uniform(10, dl)
    np.testing.assert_allclose(ax.dl, dl)
    np.testing.assert_allclose(ax.dl_dual, dl)


def test_dual_spacing_nonuniform():
    """dl_dual[i] = (dl[i] + dl[i-1])/2; at an absorbing end the ghost at i=0 extends with the same
    cell spacing."""
    edges = np.array([0.0, 1.0, 3.0, 6.0, 10.0])  # dl = 1,2,3,4
    ax = Axis(edges, "PML", "PML")
    np.testing.assert_allclose(ax.dl, [1, 2, 3, 4])
    # i=0: the ghost takes dl[0]=1 -> (1+1)/2 = 1
    np.testing.assert_allclose(ax.dl_dual, [1.0, 1.5, 2.5, 3.5])


def test_dual_spacing_periodic_wraps():
    """On a periodic axis, dl[-1] at i=0 wraps around and takes dl[n-1]."""
    edges = np.array([0.0, 1.0, 3.0, 6.0, 10.0])  # dl = 1,2,3,4
    ax = Axis(edges, "Periodic", "Periodic")
    np.testing.assert_allclose(ax.dl_dual, [2.5, 1.5, 2.5, 3.5])


def test_cell_volumes():
    grid = Grid(_uniform(2, 1.0), _uniform(3, 2.0), _uniform(4, 3.0))
    v = grid.cell_volumes()
    assert v.shape == (2, 3, 4)
    np.testing.assert_allclose(v, 6.0)


def test_nearest_index():
    edges = np.array([0.0, 1.0, 2.0, 3.0])
    assert nearest_index(edges, 0.0) == 0
    assert nearest_index(edges, 0.5) == 0
    assert nearest_index(edges, 1.0) == 1
    assert nearest_index(edges, 2.9) == 2
    assert nearest_index(edges, 3.0) == 2  # the right endpoint belongs to the last cell


def test_nearest_index_out_of_range_refuses():
    """Outside the domain it has to raise: silent clamping would make a monitor record on the wrong
    plane."""
    edges = np.array([0.0, 1.0, 2.0])
    with pytest.raises(ValueError, match="outside the grid"):
        nearest_index(edges, -0.1)
    with pytest.raises(ValueError, match="outside the grid"):
        nearest_index(edges, 2.1)


def test_non_monotonic_edges_refuse():
    with pytest.raises(ValueError, match="strictly increasing"):
        Axis(np.array([0.0, 2.0, 1.0]), "PML", "PML")


def test_yee_dual_volumes_uniform_degenerate():
    """On a uniform grid the dual volumes of all three components are equal: the degenerate,
    runtime-uniform case of the OpenEM contract."""
    grid = Grid(_uniform(4, 1.0), _uniform(5, 1.0), _uniform(6, 1.0))
    vx, vy, vz = grid.yee_dual_volumes()
    np.testing.assert_allclose(vx, 1.0)
    np.testing.assert_allclose(vy, 1.0)
    np.testing.assert_allclose(vz, 1.0)


def test_yee_dual_volumes_nonuniform_differ():
    """On a non-uniform grid the three must differ, otherwise the shutoff weighting is wrong."""
    ex = np.array([0.0, 1.0, 3.0, 6.0])          # dl = 1, 2, 3
    grid = Grid(Axis(ex, "PML", "PML"), _uniform(2, 1.0), _uniform(2, 1.0))
    vx, vy, vz = grid.yee_dual_volumes()
    # Along x: Ex uses primal (1,2,3), Ey/Ez use dual (1,1.5,2.5)
    np.testing.assert_allclose(vx[:, 0, 0], [1.0, 2.0, 3.0])
    np.testing.assert_allclose(vy[:, 0, 0], [1.0, 1.5, 2.5])
    np.testing.assert_allclose(vz[:, 0, 0], [1.0, 1.5, 2.5])
    assert not np.allclose(vx, vy)
