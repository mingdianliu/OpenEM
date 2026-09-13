#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Plane wave at normal incidence on a dielectric slab: read the transmitted and reflected flux
and compare against the analytic thin-film formula.

    python3 examples/01_slab_planewave.py

Scene: a lossless slab 0.3 um thick with eps=4 (n=2), sitting in vacuum. Periodic boundaries in
x and y, PML on both z faces. The plane wave enters from the -z side at normal incidence, with a
FluxMonitor on each side of the slab and a frequency-domain field monitor on the xz plane between
them.

What to expect: T+R is very close to 1 at all 11 frequency points (lossless medium, so energy is
conserved), and T agrees with the analytic single-layer transmittance to about 1% relative. At the
center frequency, lambda = 1 um, T is about 0.84.
"""
import numpy as np
import tidy3d as td

import openem

C0 = 2.99792458e8          # speed of light in vacuum, m/s
LAM0 = 1.0                 # center wavelength, um
N_SLAB, D_SLAB = 2.0, 0.3  # slab refractive index and thickness in um


def analytic_transmission(freqs_hz):
    """Analytic transmittance of an air/dielectric/air single layer at normal incidence.

    The layer is lossless, so T = 1 - |r|^2.
    """
    lam_um = C0 / np.asarray(freqs_hz) * 1e6
    r1 = (1.0 - N_SLAB) / (1.0 + N_SLAB)          # front interface; the back one is -r1
    phase = np.exp(2j * (2 * np.pi * N_SLAB * D_SLAB / lam_um))
    r = (r1 - r1 * phase) / (1.0 - r1 * r1 * phase)
    return 1.0 - np.abs(r) ** 2


def main():
    openem.install()                       # from here on, td.web.run lands on the local GPU

    f0 = C0 / (LAM0 * 1e-6)
    freqs = np.linspace(0.85 * f0, 1.15 * f0, 11)
    periodic = td.Boundary(plus=td.Periodic(), minus=td.Periodic())

    sim = td.Simulation(
        size=(0.6, 0.6, 3.0),
        grid_spec=td.GridSpec.uniform(dl=0.02),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, 0), size=(td.inf, td.inf, D_SLAB)),
            medium=td.Medium(permittivity=N_SLAB ** 2))],
        sources=[td.PlaneWave(
            center=(0, 0, -1.0), size=(td.inf, td.inf, 0), direction="+",
            source_time=td.GaussianPulse(freq0=f0, fwidth=0.2 * f0))],
        monitors=[
            td.FluxMonitor(center=(0, 0, 1.0), size=(td.inf, td.inf, 0),
                           freqs=list(freqs), name="T"),
            td.FluxMonitor(center=(0, 0, -1.2), size=(td.inf, td.inf, 0),
                           freqs=list(freqs), name="R"),
            td.FieldMonitor(center=(0, 0, 0), size=(td.inf, 0, td.inf),
                            freqs=[f0], name="fxz"),
        ],
        run_time=4e-13,
        boundary_spec=td.BoundarySpec(x=periodic, y=periodic, z=td.Boundary.pml()))

    sim_data = td.web.run(sim, task_name="slab_planewave")

    T = np.asarray(sim_data["T"].flux.values, dtype=float)
    R = -np.asarray(sim_data["R"].flux.values, dtype=float)   # reflected wave travels -z, so its flux is negative
    T_ref = analytic_transmission(freqs)
    i0 = len(freqs) // 2                                      # the midpoint is exactly f0

    print()
    print(f"cells              : {np.prod(sim.grid.num_cells)}  {tuple(sim.grid.num_cells)}")
    print(f"T(lambda=1um)      : {T[i0]:.4f}   analytic {T_ref[i0]:.4f}")
    print(f"R(lambda=1um)      : {R[i0]:.4f}   analytic {1 - T_ref[i0]:.4f}")
    print(f"T+R over all freqs : min {(T + R).min():.4f}  max {(T + R).max():.4f}")
    print(f"|T-analytic|/analytic : max {np.abs(T - T_ref).max() / T_ref.max():.2e}")
    Ex = np.asarray(sim_data["fxz"].Ex.values)
    print(f"xz plane max|Ex|   : {np.abs(Ex).max():.4f}   array shape {Ex.shape}")


if __name__ == "__main__":
    main()
