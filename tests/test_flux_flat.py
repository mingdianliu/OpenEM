# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""How the colocated flux convention behaves in 2D (one flat tangential axis) and when the geometry
fails (2026-09-07).

1. Flat axis: the old plane_flux convention counts the width of a flat axis as FLAT_AXIS_M, and the
   colocated convention has to use the same rule. Multiply the real thickness of that flat cell by
   10 and neither convention may change its result; on a smooth field the two differ only by
   O(h²). EffectiveIndexApproximation used to come out as 3e12 precisely because the real
   thickness was used.
2. Flush against a non-absorbing boundary: the geometry fails closed, monitor_flux falls back to
   the old convention, and the result is bitwise identical to plane_flux (the box faces of
   BullseyeCavityPSO).
"""
from __future__ import annotations
import numpy as np
import pytest

from openem import flux as flux_mod
from openem.grid import Axis, Grid
from openem.model import FluxMonitor


def _axis(edges, lo="PML", hi="PML"):
    return Axis(edges=np.asarray(edges, dtype=np.float64), boundary_lo=lo, boundary_hi=hi)


def _grid2d(thick, n=48, L=1.2e-6):
    e = np.linspace(-L / 2, L / 2, n + 1)
    return Grid(x=_axis(e), y=_axis(e), z=_axis(np.array([-thick / 2, thick / 2]), "Periodic", "Periodic"))


def _plane_gauss(grid, w=0.25e-6, eta=376.73):
    """Plane with x normal, tangential (y, z), z flat: E_y=g(y), H_z=g(y)/η, each at its own Yee
    position."""
    ay = grid.axes[1]
    cy, ey = 0.5 * (ay.edges[:-1] + ay.edges[1:]), ay.edges[:-1]
    g = lambda y: np.exp(-(y ** 2) / w ** 2)
    eu = g(cy)[:, None]        # E_t1 = E_y: at the cell centre along y
    hv = (g(ey) / eta)[:, None]   # H_t2 = H_z: on the cell boundary along y
    zero = np.zeros_like(eu)
    return np.stack([eu, zero, zero, hv])[:, None].astype(np.complex128)


@pytest.mark.parametrize("thick", [1e-8, 1e-7])
def test_flat_axis_uses_flat_axis_m(thick, monkeypatch):
    monkeypatch.setenv("OPENEM_FLUX_COLOCATED", "1")
    grid = _grid2d(thick)
    mon = FluxMonitor(name="m", axis=0, plane_index=24, normal_dir=+1, freqs=np.array([2e14]),
                      source_spectrum=np.array([1 + 0j]), transverse=None,
                      t_bounds=((-0.4e-6, 0.4e-6), (-thick, thick)))
    plane = _plane_gauss(grid)
    got = flux_mod.monitor_flux(plane, grid, mon)[0]
    e = grid.axes[1].edges[:-1]
    keep = np.flatnonzero((e >= -0.4e-6 - 1e-12) & (e <= 0.4e-6 + 1e-12))
    legacy = flux_mod.plane_flux(plane, grid, transverse=((keep[0], keep[-1]), (0, 0)), axis=0)[0]
    # Both conventions count the flat axis as FLAT_AXIS_M, so the magnitudes agree; on a smooth
    # field they differ by O(h²)
    assert abs(got - legacy) / legacy < 5e-3, (got, legacy)
    return got


def test_flat_axis_result_independent_of_thickness(monkeypatch):
    monkeypatch.setenv("OPENEM_FLUX_COLOCATED", "1")
    vals = []
    for thick in (1e-8, 1e-7, 1e-6):
        grid = _grid2d(thick)
        mon = FluxMonitor(name="m", axis=0, plane_index=24, normal_dir=+1, freqs=np.array([2e14]),
                          source_spectrum=np.array([1 + 0j]), transverse=None,
                          t_bounds=((-0.4e-6, 0.4e-6), (-thick, thick)))
        vals.append(flux_mod.monitor_flux(_plane_gauss(grid), grid, mon)[0])
    assert np.allclose(vals, vals[0], rtol=1e-12), vals


def test_both_flat_axes_fall_back_to_legacy_bitwise(monkeypatch, capsys):
    """A 1D simulation (DistributedBraggReflectorCavity, grid [1,1,241]): both tangential axes are
    flat, the colocated convention is meaningless, so it must fall back to the old convention and be
    bitwise identical. Before the fix a run came out as 6e-6, off by a factor of 1e4."""
    monkeypatch.setenv("OPENEM_FLUX_COLOCATED", "1")
    th = 1e-8
    ez = np.linspace(-0.8e-6, 0.8e-6, 242)
    grid = Grid(x=_axis(np.array([-th / 2, th / 2]), "Periodic", "Periodic"),
                y=_axis(np.array([-th / 2, th / 2]), "Periodic", "Periodic"), z=_axis(ez))
    mon = FluxMonitor(name="R", axis=2, plane_index=30, normal_dir=-1, freqs=np.array([2e14, 3e14]),
                      source_spectrum=np.array([1 + 0j, 1 + 0j]), transverse=None,
                      t_bounds=((-th, th), (-th, th)))
    rng = np.random.default_rng(9)
    plane = rng.standard_normal((4, 2, 1, 1)) + 1j * rng.standard_normal((4, 2, 1, 1))
    got = flux_mod.monitor_flux(plane, grid, mon)
    assert np.array_equal(got, flux_mod.plane_flux(plane, grid, transverse=None, axis=2))
    assert "in-place Yee integration" in capsys.readouterr().out


def test_wall_falls_back_to_legacy_bitwise(monkeypatch, capsys):
    monkeypatch.setenv("OPENEM_FLUX_COLOCATED", "1")
    rng = np.random.default_rng(5)
    n = 14
    e = np.linspace(-0.7e-6, 0.7e-6, n + 1)
    grid = Grid(x=_axis(e), y=_axis(e, "PECBoundary", "PML"), z=_axis(e))
    mon = FluxMonitor(name="face", axis=0, plane_index=7, normal_dir=+1, freqs=np.array([2e14]),
                      source_spectrum=np.array([1 + 0j]), transverse=None,
                      t_bounds=((-1e-6, 1e-6), (-0.3e-6, 0.3e-6)))     # the y end is flush with a PEC wall
    plane = rng.standard_normal((4, 1, n, n)) + 1j * rng.standard_normal((4, 1, n, n))
    got = flux_mod.monitor_flux(plane, grid, mon)
    assert np.array_equal(got, flux_mod.plane_flux(plane, grid, transverse=None, axis=0))
    assert "in-place Yee integration" in capsys.readouterr().out
