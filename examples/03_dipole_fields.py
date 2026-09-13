#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Point dipole in a uniform medium: box flux against the analytic radiated power, plus a
frequency-domain field slice for the peak position and the symmetry.

    python3 examples/03_dipole_fields.py

Scene: a uniform medium of n=1.5 with no structures at all, an Ez-polarized point dipole at the
origin, PML on all six faces. Two concentric box flux monitors (0.4 and 0.8 um on a side) enclose
the dipole, and one frequency-domain field monitor sits on the z=0 plane.

Why this is a real check: Tidy3D's PointDipole is a fixed current-density source, so amplitude=1
means a current moment of 1 A.um. Monitor data is normalized by the source spectrum by default,
so the box flux reads out as "radiated power per unit squared current moment", which is exactly the
Hertzian dipole result P = n.mu0.omega^2.(1 A.um)^2 / (12.pi.c). There is no free parameter here.

What to expect: the outer box flux agrees with the analytic value to about 1e-4 relative at all
three frequency points, and the two boxes agree with each other to about 2e-3. In a uniform
lossless medium the flux does not depend on the box size, which makes that second comparison a
conservation check independent of any amplitude convention. The |Ez| peak lands on the dipole cell,
and |Ez| on the z=0 plane is exactly symmetric under swapping x and y, since the discretization is
identical along the two axes.
"""
import numpy as np
import tidy3d as td

import openem

C0 = 2.99792458e8            # speed of light in vacuum, m/s
MU0 = 4.0e-7 * np.pi         # vacuum permeability, H/m
UM = 1e-6                    # a current moment of 1 A.um expressed in A.m
N_BG = 1.5                   # refractive index of the background medium
FREQS = [1.7e14, 2.0e14, 2.3e14]      # centered on 200 THz, vacuum wavelength 1.5 um


def analytic_power(freqs_hz):
    """Total radiated power of a Hertzian dipole in a uniform medium, for a current moment of 1 A.um."""
    w = 2 * np.pi * np.asarray(freqs_hz)
    return N_BG * MU0 * w ** 2 * UM ** 2 / (12 * np.pi * C0)


def row(values):
    return "  ".join(f"{v:+.4e}" for v in values)


def main():
    openem.install()                       # from here on, td.web.run lands on the local GPU

    sim = td.Simulation(
        size=(2.0, 2.0, 2.0),
        grid_spec=td.GridSpec.uniform(dl=0.02),
        medium=td.Medium(permittivity=N_BG ** 2),
        sources=[td.PointDipole(
            center=(0, 0, 0), polarization="Ez",
            source_time=td.GaussianPulse(freq0=FREQS[1], fwidth=0.25 * FREQS[1]))],
        monitors=[
            td.FluxMonitor(center=(0, 0, 0), size=(0.4, 0.4, 0.4), freqs=FREQS, name="inner"),
            td.FluxMonitor(center=(0, 0, 0), size=(0.8, 0.8, 0.8), freqs=FREQS, name="outer"),
            td.FieldMonitor(center=(0, 0, 0), size=(1.0, 1.0, 0), freqs=[FREQS[1]], name="fxy"),
        ],
        run_time=2e-13,
        boundary_spec=td.BoundarySpec.all_sides(td.PML()))

    sim_data = td.web.run(sim, task_name="dipole_fields")

    p_ref = analytic_power(FREQS)
    p_in = np.asarray(sim_data["inner"].flux.values, dtype=float)
    p_out = np.asarray(sim_data["outer"].flux.values, dtype=float)
    ez = np.abs(np.asarray(sim_data["fxy"].Ez.values)).squeeze()
    x = np.asarray(sim_data["fxy"].Ez.coords["x"].values)
    y = np.asarray(sim_data["fxy"].Ez.coords["y"].values)
    ix, iy = np.unravel_index(int(np.argmax(ez)), ez.shape)

    print()
    print(f"cells                  : {np.prod(sim.grid.num_cells)}  {tuple(sim.grid.num_cells)}")
    print(f"frequencies (THz)      : " + "  ".join(f"{f / 1e12:.0f}" for f in FREQS))
    print(f"analytic radiated power: {row(p_ref)}")
    print(f"outer box flux         : {row(p_out)}")
    print(f"outer/analytic - 1     : {row(p_out / p_ref - 1)}")
    print(f"inner/analytic - 1     : {row(p_in / p_ref - 1)}")
    print(f"inner/outer - 1        : {row(p_in / p_out - 1)}")
    print(f"|Ez| peak at           : x={x[ix]:+.3f} y={y[iy]:+.3f} um (dipole is at the origin)")
    print(f"z=0 plane x-y symmetry : max|A - A.T| / max A = {np.abs(ez - ez.T).max() / ez.max():.1e}")


if __name__ == "__main__":
    main()
