#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Adjoint gradient against a central difference: the derivative of a downstream mode amplitude
with respect to the length of a loading block.

    python3 examples/04_autograd_gradient.py

Scene: a silicon strip waveguide of eps=12.1, 0.5 um wide and 0.22 um thick, with a loading block
of the same material suspended 0.1 um above it (length t, width 0.4 um, thickness 0.2 um) and PML
on all six faces. An upstream ModeSource injects the fundamental mode and a downstream ModeMonitor
reads its forward amplitude. The objective is J(t) = |amp+(f0)|^2 and t is the only differentiable
parameter: a longer block means a worse mode match, so J drops.

The check: dJ/dt from tidy3d's autograd (one forward solve plus one adjoint solve, both on the
local GPU) against a central difference of the same objective, (J(t+h) - J(t-h)) / (2h). Four
solves in total, about two minutes from cold. Three things that will bite you:

1. The objective uses a mode amplitude, not a flux. The adjoint source of a flux objective is a
   CustomCurrentSource that expands into several thousand point dipoles; batched injection
   accumulates with atomicAdd, float addition is not associative, and the gradient stops being
   bitwise reproducible. A mode objective gets a ModeSource adjoint, which is reproducible.
2. openem.install() does not pin the step count (steps=None): the forward solve stops on its own
   criterion and the adjoint copies that count. Pin it too low and the phase has not converged, so
   J itself is noisy and the central difference falls apart first.
3. The block is separated from the waveguide by 0.1 um of air and is deliberately a different width.
   Butt them together, with coincident faces and the same material, and the same code returns an
   adjoint gradient only 0.65 times the true slope. A shape derivative is a continuous surface
   integral over an interface, which requires that interface to separate exactly two things: inside
   the block, and background. The difference step is a full cell (h = dl), which steps over the
   S-shaped response that subpixel averaging has within one cell.

What to expect: J is about 0.933, the adjoint gradient about -1.66, the central difference about
-1.69, a relative difference of roughly 2%.
"""
import time
import autograd
import autograd.numpy as anp
import tidy3d as td

import openem

F0 = 1.934e14          # about 1.55 um
DL = 0.025             # grid step, um
T0 = 0.40              # operating point for the block length, um
H_FD = DL              # central-difference step: one full cell


def make_sim(t):
    """t is the length of the loading block along x, in um.

    When t carries an autograd tracer, so does the structure built from it.
    """
    si = td.Medium(permittivity=12.1)
    plane = dict(size=(0, 1.2, 1.2), mode_spec=td.ModeSpec(num_modes=1))
    return td.Simulation(
        size=(2.0, 1.8, 1.8), grid_spec=td.GridSpec.uniform(dl=DL),
        structures=[
            td.Structure(geometry=td.Box(center=(0, 0, 0), size=(td.inf, 0.5, 0.22)), medium=si),
            td.Structure(geometry=td.Box(center=(0, 0, 0.31), size=(t, 0.4, 0.2)), medium=si)],
        sources=[td.ModeSource(
            center=(-0.6, 0, 0), **plane, direction="+", mode_index=0,
            source_time=td.GaussianPulse(freq0=F0, fwidth=0.1 * F0))],
        monitors=[td.ModeMonitor(center=(0.6, 0, 0), **plane, freqs=[F0], name="m")],
        run_time=3e-13, boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=8)))


def objective(t):
    """J(t): squared magnitude of the downstream fundamental mode's forward amplitude."""
    sim_data = td.web.run(make_sim(t), task_name="autograd_gradient", verbose=False)
    amp = sim_data["m"].amps.sel(direction="+", f=F0, mode_index=0)
    return anp.sum(anp.abs(amp.values) ** 2)


def main():
    openem.install()                       # from here on, td.web.run and autograd both use the local GPU
    mark = time.time()
    J, grad = autograd.value_and_grad(objective)(T0)
    t_adjoint = time.time() - mark
    mark = time.time()
    j_plus, j_minus = float(objective(T0 + H_FD)), float(objective(T0 - H_FD))
    t_fd = time.time() - mark
    fd = (j_plus - j_minus) / (2 * H_FD)

    print()
    print(f"block length t        : {T0} um    grid dl = {DL} um, difference step h = {H_FD} um")
    print(f"J(t) = |amp+|^2       : {float(J):.6f}")
    print(f"J(t+h) / J(t-h)       : {j_plus:.6f} / {j_minus:.6f}")
    print(f"adjoint dJ/dt         : {float(grad):+.4f}    forward + adjoint, two solves, {t_adjoint:.0f}s")
    print(f"central difference    : {fd:+.4f}    two more forward solves, {t_fd:.0f}s")
    print(f"relative difference   : {abs(float(grad) - fd) / abs(fd):.2%}")


if __name__ == "__main__":
    main()
