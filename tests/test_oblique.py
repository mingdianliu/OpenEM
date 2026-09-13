# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Obliquely incident plane wave on the Bloch complex path: R against Fresnel at a single
interface, and R+T conservation.

The machinery under test is complex single-face TF/SF injection (inject_mode_*_c), Bloch phase
matching, and the complex-field frequency-domain FieldMonitor (accumulate_dft_box_c).
"""
from __future__ import annotations
import numpy as np
import pytest

td = pytest.importorskip("tidy3d")
pytest.importorskip("cupy")

from openem import scene as scene_mod, solver
from openem.device import Kernels

F0, N1, N2, THETA, DL = 2e14, 1.0, 1.5, np.pi / 6, 0.02


def _build(interface: bool):
    lam = td.C_0 / F0
    Lx = 20 * DL
    kt = 2 * np.pi / lam * N1 * np.sin(THETA)
    bx = kt * Lx / (2 * np.pi)
    structs = ([td.Structure(geometry=td.Box(center=(0, 0, 1.0), size=(td.inf, td.inf, 2.0)),
                             medium=td.Medium(permittivity=N2**2))] if interface else [])
    src = td.PlaneWave(center=(0, 0, -0.8), size=(td.inf, td.inf, 0), direction="+",
                       source_time=td.GaussianPulse(freq0=F0, fwidth=0.1 * F0),
                       angle_theta=THETA, angle_phi=0.0, pol_angle=np.pi / 2)
    mons = [td.FieldMonitor(center=(0, 0, z), size=(td.inf, td.inf, 0), freqs=[F0], name=nm)
            for nm, z in (("up", -1.3), ("mid", -0.3), ("dn", 1.6))]
    per = td.Boundary(minus=td.BlochBoundary(bloch_vec=bx), plus=td.BlochBoundary(bloch_vec=bx))
    return td.Simulation(size=(Lx, Lx, 4.0), grid_spec=td.GridSpec.uniform(dl=DL),
                         structures=structs, sources=[src], monitors=mons, run_time=8e-13,
                         medium=td.Medium(permittivity=N1**2),
                         boundary_spec=td.BoundarySpec(x=per, y=td.Boundary.periodic(),
                             z=td.Boundary(minus=td.PML(), plus=td.PML())))


def _sz(res, nm):
    ph = res.field_phasors[nm]
    Ex, Ey, Hx, Hy = ph[0, 0], ph[1, 0], ph[3, 0], ph[4, 0]
    return float(np.sum(0.5 * np.real(Ex * np.conj(Hy) - Ey * np.conj(Hx))))


def test_oblique_planewave_matches_fresnel():
    k = Kernels()
    ref = solver.run(scene_mod.from_simulation(_build(False)), kernels=k, verbose=False)
    P_inc = abs(_sz(ref, "mid"))
    assert abs(_sz(ref, "up")) / P_inc < 1e-4, "TF/SF directionality is broken: there is field in the reflected region"
    it = solver.run(scene_mod.from_simulation(_build(True)), kernels=k, verbose=False)
    R = abs(_sz(it, "up")) / P_inc
    T = abs(_sz(it, "dn")) / P_inc
    sint2 = N1 * np.sin(THETA) / N2
    r_s = (N1 * np.cos(THETA) - N2 * np.sqrt(1 - sint2**2)) / \
          (N1 * np.cos(THETA) + N2 * np.sqrt(1 - sint2**2))
    assert abs(R - r_s**2) < 0.01, f"R={R:.4f} vs Fresnel {r_s**2:.4f}"
    assert abs(R + T - 1) < 0.03, f"energy is not conserved, R+T={R + T:.4f}"
