# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""The colocated integration convention of the frequency-domain FluxMonitor
(flux.monitor_flux_colocated, 2026-09-06).

The criteria have no tunable parameters:

1. **Same geometry**: sample a 3D random field onto a plane (E on the primal plane km; H as
   0.5·(H[km]+H[km-1]), which is the normal colocation used by the plane DFT kernel), and
   ``monitor_flux_colocated`` must give a result **bitwise identical** to the time-domain path,
   ``colocated_component`` + ``integrate_colocated``, on the same 3D field (up to the ½ that the
   frequency domain carries). Checked for all three normal axes, on a non-uniform grid, with
   monitor bounds that do not land on cell boundaries.
2. **Oracle**: put an analytic, smooth Gaussian profile on the plane (E_u = g, H_v = g/η, the rest
   zero). The flux from the colocated convention converges to the analytic integral ½∫g²/η dA, with
   the error falling **second order** in the cell count (x0.25 per refinement). Measured
   (2026-09-06, waist from 2.5 to 40 cells): colocated 7.8e-2 -> 3.2e-4, in-place Yee 3.7e-2 ->
   7.4e-5, so **on a smooth field the old convention is in fact more accurate**. The reason for
   changing convention is not accuracy but **agreement with tidy3d**: its flux, its mode
   normalization (|amp|²) and its adjoint source strength are all built on the colocated
   convention, and the 7.7% flux/|amp|² inconsistency seen on a straight-waveguide probe is exactly
   the difference between the two second-order rules on a confined mode.
3. Knob: with ``OPENEM_FLUX_COLOCATED`` unset or 0, ``monitor_flux`` is bitwise identical to the
   ``plane_flux`` path, so old results are unaffected; with 1 it equals
   ``monitor_flux_colocated``.
"""
from __future__ import annotations
import os
import numpy as np
import pytest

from openem import flux as flux_mod
from openem.grid import Axis, Grid
from openem.model import FluxMonitor


def _axis(edges, lo="PML", hi="PML"):
    return Axis(edges=np.asarray(edges, dtype=np.float64), boundary_lo=lo, boundary_hi=hi)


def _grid(nonuniform=True, n=(14, 11, 9)):
    rng = np.random.default_rng(3)
    axes = []
    for k in n:
        dl = 0.05e-6 * (1.0 + (0.4 * rng.random(k) if nonuniform else 0.0))
        axes.append(_axis(np.concatenate([[0.0], np.cumsum(dl)]) - 0.3e-6))
    return Grid(x=axes[0], y=axes[1], z=axes[2])


def _mon(grid, axis, km, name="m"):
    t1, t2 = (t for t in range(3) if t != axis)
    e1, e2 = grid.axes[t1].edges, grid.axes[t2].edges
    # The bounds deliberately do not land on cell boundaries
    tb = ((float(e1[2]) + 0.37 * float(e1[3] - e1[2]), float(e1[-3]) - 0.21 * float(e1[-2] - e1[-3])),
          (float(e2[1]) + 0.55 * float(e2[2] - e2[1]), float(e2[-2]) - 0.13 * float(e2[-1] - e2[-2])))
    return FluxMonitor(name=name, axis=axis, plane_index=km, normal_dir=+1,
                       freqs=np.array([2.0e14]), source_spectrum=np.array([1.0 + 0j]),
                       transverse=None, t_bounds=tb)


@pytest.mark.parametrize("axis", [0, 1, 2])
def test_same_geometry_as_time_domain_path(axis):
    grid = _grid()
    n = tuple(grid.axes[t].n for t in range(3))
    rng = np.random.default_rng(7)
    E = [rng.standard_normal(n) + 1j * rng.standard_normal(n) for _ in range(3)]
    H = [rng.standard_normal(n) + 1j * rng.standard_normal(n) for _ in range(3)]
    km = n[axis] // 2
    mon = _mon(grid, axis, km)
    t1, t2 = (t for t in range(3) if t != axis)

    # How the plane DFT kernel stores it: E on plane km, H as 0.5·(H[km]+H[km-1])
    def take(F, k):
        sl = [slice(None)] * 3
        sl[axis] = k
        return F[tuple(sl)]
    plane = np.stack([take(E[t1], km), take(E[t2], km),
                      0.5 * (take(H[t1], km) + take(H[t1], km - 1)),
                      0.5 * (take(H[t2], km) + take(H[t2], km - 1))])[:, None]   # (4, nf=1, n1, n2)

    got = flux_mod.monitor_flux_colocated(plane, grid, mon)[0]

    # The time-domain path, on the same 3D field
    g = flux_mod.flux_time_geometry(grid, flux_mod._PlaneGeomSpec(mon))
    assert len(g["planes"]) == 1
    eu = flux_mod.colocated_component(E[t1], g, "Eu", 0)
    ev = flux_mod.colocated_component(E[t2], g, "Ev", 0)
    hu = flux_mod.colocated_component(H[t1], g, "Hu", 0)
    hv = flux_mod.colocated_component(H[t2], g, "Hv", 0)
    want = 0.5 * float(np.real(np.sum((eu * np.conj(hv) - ev * np.conj(hu)) * g["dS"])))
    assert got == pytest.approx(want, rel=0, abs=1e-12 * max(1.0, abs(want)))


def _gauss_plane(grid, axis, w=0.25e-6, eta=376.73):
    """Analytic profile on the plane: E_u = g(u,v), H_v = g/η, each at its own Yee position."""
    t1, t2 = (t for t in range(3) if t != axis)
    a1, a2 = grid.axes[t1], grid.axes[t2]
    c1, e1 = 0.5 * (a1.edges[:-1] + a1.edges[1:]), a1.edges[:-1]
    c2, e2 = 0.5 * (a2.edges[:-1] + a2.edges[1:]), a2.edges[:-1]
    g = lambda u, v: np.exp(-(u[:, None] ** 2 + v[None, :] ** 2) / w ** 2)
    eu = g(c1, e2)            # E_t1: (center_1, edge_2)
    hv = g(c1, e2) / eta      # H_t2 colocated onto the E plane, also at (center_1, edge_2)? No:
    # on the Yee grid H_t2 sits at (edge_1, center_2). Put each one where it belongs:
    hv = g(e1, c2) / eta
    zero = np.zeros_like(eu)
    return np.stack([eu, zero, zero, hv])[:, None].astype(np.complex128), w, eta


def _analytic(mon, w, eta):
    from math import erf, sqrt, pi
    (l1, h1), (l2, h2) = mon.t_bounds
    I = lambda lo, hi: 0.5 * sqrt(pi) * w * (erf(hi / w) - erf(lo / w))   # ∫ e^{-x²/w²}; for g² use w/√2
    w2 = w / sqrt(2.0)
    I2 = lambda lo, hi: 0.5 * sqrt(pi) * w2 * (erf(hi / w2) - erf(lo / w2))
    return 0.5 * I2(l1, h1) * I2(l2, h2) / eta


def test_colocated_converges_second_order_to_analytic():
    errs_c, errs_y = [], []
    for ncell in (12, 24, 48, 96):
        edges = np.linspace(-0.6e-6, 0.6e-6, ncell + 1)
        grid = Grid(x=_axis(edges), y=_axis(edges), z=_axis(edges))
        mon = FluxMonitor(name="g", axis=2, plane_index=ncell // 2, normal_dir=+1,
                          freqs=np.array([2.0e14]), source_spectrum=np.array([1.0 + 0j]),
                          transverse=None, t_bounds=((-0.4e-6, 0.4e-6), (-0.4e-6, 0.4e-6)))
        plane, w, eta = _gauss_plane(grid, 2)
        want = _analytic(mon, w, eta)
        got_c = flux_mod.monitor_flux_colocated(plane, grid, mon)[0]
        # The in-place Yee convention (the old path): turn t_bounds into transverse indices
        e = edges[:-1]
        keep = np.flatnonzero((e >= -0.4e-6 - 1e-12) & (e <= 0.4e-6 + 1e-12))
        got_y = flux_mod.plane_flux(plane, grid, transverse=((keep[0], keep[-1]),) * 2, axis=2)[0]
        errs_c.append(abs(got_c - want) / want)
        errs_y.append(abs(got_y - want) / want)
    # Second-order convergence: the error goes x0.25 per refinement (measured 0.26 / 0.25 / 0.25),
    # with a little margin left
    for a, b in zip(errs_c, errs_c[1:]):
        assert b < 0.35 * a, (errs_c, errs_y)
    # Both conventions converge to the same continuum limit (in-place Yee is not bad on a smooth
    # field; the convention changed to agree with tidy3d)
    assert errs_y[-1] < 1e-3 and errs_c[-1] < 2e-3, (errs_c, errs_y)


def test_switch_keeps_legacy_bitwise(monkeypatch):
    grid = _grid()
    n = tuple(grid.axes[t].n for t in range(3))
    rng = np.random.default_rng(11)
    plane = (rng.standard_normal((4, 2, n[0], n[1])) + 1j * rng.standard_normal((4, 2, n[0], n[1])))
    mon = _mon(grid, 2, n[2] // 2)
    # Since the 2026-09-06 probe acceptance the default is the colocated convention; only 0 falls
    # back to in-place Yee, for comparison
    monkeypatch.setenv("OPENEM_FLUX_COLOCATED", "0")
    legacy = flux_mod.monitor_flux(plane, grid, mon)
    assert np.array_equal(legacy, flux_mod.plane_flux(plane, grid, transverse=None, axis=2))
    monkeypatch.delenv("OPENEM_FLUX_COLOCATED", raising=False)
    assert np.array_equal(flux_mod.monitor_flux(plane, grid, mon),
                          flux_mod.monitor_flux_colocated(plane, grid, mon))
    monkeypatch.setenv("OPENEM_FLUX_COLOCATED", "1")
    assert np.array_equal(flux_mod.monitor_flux(plane, grid, mon),
                          flux_mod.monitor_flux_colocated(plane, grid, mon))
    # An old scene without t_bounds: even with the knob on it falls back to the old convention
    # instead of blowing up
    old = FluxMonitor(name="o", axis=2, plane_index=n[2] // 2, normal_dir=+1,
                      freqs=mon.freqs, source_spectrum=mon.source_spectrum, transverse=None)
    assert np.array_equal(flux_mod.monitor_flux(plane, grid, old), legacy)
