#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Silicon strip waveguide: a mode source injects the fundamental mode, a mode monitor reads the
transmission, and a mode solver monitor reads the effective index.

    python3 examples/02_waveguide_mode.py

Scene: a strip waveguide of eps=12.1, 0.5 um wide and 0.22 um thick, in vacuum with PML on all six
faces. A ModeSource at x = -0.6 um injects mode_index=0. At x = +0.5 um sit a ModeMonitor (2 modes,
3 frequency points) and a FluxMonitor, with a ModeSolverMonitor between them.

What to expect: the waveguide is straight, lossless and non-scattering, so the forward |amp|^2 of
the fundamental mode and the total flux both come out very close to 1. The higher-order mode
(mode_index=1) should be several orders of magnitude smaller. The fundamental n_eff is about 2.4,
between the cladding index of 1 and the core index sqrt(12.1) = 3.48.
"""
import numpy as np
import tidy3d as td

import openem

F0 = 1.934e14                 # about 1.55 um


def main():
    openem.install()

    mode_spec = td.ModeSpec(num_modes=2)
    plane = dict(size=(0, 1.8, 1.4))
    sim = td.Simulation(
        size=(2.0, 2.4, 2.0),
        grid_spec=td.GridSpec.uniform(dl=0.05),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, 0), size=(td.inf, 0.5, 0.22)),
            medium=td.Medium(permittivity=12.1))],
        sources=[td.ModeSource(
            center=(-0.6, 0, 0), **plane, direction="+", mode_index=0,
            mode_spec=td.ModeSpec(num_modes=1),
            source_time=td.GaussianPulse(freq0=F0, fwidth=0.1 * F0))],
        monitors=[
            td.ModeMonitor(center=(0.5, 0, 0), **plane, mode_spec=mode_spec,
                           freqs=[0.98 * F0, F0, 1.02 * F0], name="m"),
            td.FluxMonitor(center=(0.5, 0, 0), **plane, freqs=[F0], name="fl"),
            td.ModeSolverMonitor(center=(0.0, 0, 0), **plane, mode_spec=mode_spec,
                                 freqs=[F0], name="ms"),
        ],
        run_time=2e-13,
        boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=8)))

    sim_data = td.web.run(sim, task_name="waveguide_mode")

    amps = sim_data["m"].amps.sel(direction="+")
    p_fwd = np.abs(np.asarray(amps.values)) ** 2          # shape (f, mode_index)
    flux = float(np.real(sim_data["fl"].flux.values[0]))
    n_eff = np.real(np.asarray(sim_data["ms"].n_eff.values)).ravel()

    print()
    print(f"cells                : {np.prod(sim.grid.num_cells)}  {tuple(sim.grid.num_cells)}")
    print(f"n_eff (f0)           : mode0 {n_eff[0]:.4f}   mode1 {n_eff[1]:.4f}")
    print(f"|amp+|^2 mode0       : " + "  ".join(f"{v:.4f}" for v in p_fwd[:, 0]))
    print(f"|amp+|^2 mode1       : " + "  ".join(f"{v:.3e}" for v in p_fwd[:, 1]))
    print(f"total flux (f0)      : {flux:.4f}")
    print(f"|amp0|^2/flux (f0)   : {p_fwd[1, 0] / flux:.4f}")


if __name__ == "__main__":
    main()
