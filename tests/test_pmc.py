# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""PMCBoundary, a perfect magnetic conductor wall: a Neumann condition on tangential E.

An analytic criterion with no free parameters. Send a plane wave into the wall in vacuum and
measure the standing wave of the single-frequency phasor:

* PEC wall: tangential E has a **node** on the wall;
* PMC wall: tangential E has an **antinode** on the wall.

The last sample point of the monitor is not on the wall, so the criterion is written as the
analytic value at its distance from the wall:

===========  ==============================  ==================
wall            analytic \\|Ex\\|/max at the last point   effective wall
=============  =====================================  ==========================
PEC high wall   ``sin(2*pi/N)``                        nominal position (node n)
PMC low wall    ``1.0`` (the first point is on it)     nominal position (node 0)
PMC high wall   ``cos(pi/N)``                          **nominal position, pulled in by dl/2**
===========  ==============================  ==================

``N`` is the number of cells per wavelength. The last row is not a typo but a measured conclusion,
confirmed to three decimal places at both N=40 and N=80, with "the wall is at its nominal position"
ruled out:

**A PMC high wall is pulled in by half a cell from its nominal position.** The root cause is
drop-last: tangential E on the high wall is a free degree of freedom yet is not stored, so all that
can be done is to force the difference in the H update to zero, which puts dE_t/dn = 0 at the last
H node z_{n-1/2} rather than at the wall, z_n. The low wall does not have this problem:
``mpv[0] = -1`` makes H(-1/2) = -H(+1/2), so H = 0 lands exactly on node 0.

Fixing it needs tangential E stored at node n, one extra column on that axis, or a special case in
the kernel. Until then it is an **O(dl) geometric error**: a PEC-plus-PMC cavity is effectively
dl/2 short. This test pins the current state down, so it goes red when someone does fix it, at
which point the expected value on the PMC high wall row becomes ``cos(2*pi/N)``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

import pytest

# This module imports openem.solver at top level, and that imports cupy at top level, so a CPU-only
# machine would crash during collection.
pytest.importorskip("cupy", reason="needs CUDA and cupy")

from openem import scene as scene_mod, solver
from openem.device import Kernels

F0 = 2.0e14
LAM = 2.99792458e8 / F0 * 1e6        # µm
LZ = 3.0


def _sim(wall_lo, wall_hi, n_per_lam, mon_lo, mon_hi, src_z, direction):
    dl = LAM / n_per_lam
    return td.Simulation(
        size=(4 * dl, 4 * dl, LZ), run_time=4e-13,
        grid_spec=td.GridSpec.uniform(dl=dl),
        sources=[td.PlaneWave(
            center=(0, 0, src_z), size=(td.inf, td.inf, 0),
            direction=direction, pol_angle=0.0,
            source_time=td.GaussianPulse(freq0=F0, fwidth=0.1 * F0))],
        monitors=[td.FieldMonitor(
            center=(0, 0, 0.5 * (mon_lo + mon_hi)),
            size=(0, 0, mon_hi - mon_lo), freqs=[F0],
            fields=["Ex"], name="line")],
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.periodic(), y=td.Boundary.periodic(),
            z=td.Boundary(minus=wall_lo, plus=wall_hi)),
        subpixel=False)


def _profile(sim):
    sc = scene_mod.from_simulation(sim)
    res = solver.run(sc, num_steps=None, use_shutoff=True,
                     kernels=Kernels(), verbose=False)
    return np.abs(np.asarray(res.field_phasors["line"])[0, 0]).ravel()


@pytest.mark.parametrize("n_per_lam", [40, 80])
def test_pec_high_wall_node_at_nominal_position(n_per_lam):
    """PEC high wall: the last point sits one cell inside the wall, so |Ex| should be sin(2*pi/N)."""
    ex = _profile(_sim(td.PML(), td.PECBoundary(), n_per_lam,
                       0.0, LZ / 2, -1.2, "+"))
    got = ex[-1] / ex.max()
    assert got == pytest.approx(math.sin(2 * math.pi / n_per_lam), abs=3e-3)


@pytest.mark.parametrize("n_per_lam", [40, 80])
def test_pmc_high_wall_antinode_half_cell_inside(n_per_lam):
    """PMC high wall: an antinode rather than a node, but with the effective wall pulled in by half a
    cell.

    Both assertions are needed: the first proves the boundary condition type is right, since
    otherwise it would be a node as with PEC, and the second pins down that half cell while
    explicitly ruling out "the wall is at its nominal position".
    """
    ex = _profile(_sim(td.PML(), td.PMCBoundary(), n_per_lam,
                       0.0, LZ / 2, -1.2, "+"))
    got = ex[-1] / ex.max()
    assert got > 0.9, f"a PMC wall should carry an antinode, measured {got:.4f}"
    # Tolerance 5e-4: it measures to within 1e-4, while the two hypotheses differ by 2.3e-3 at N=80.
    # A loose 3e-3 tolerance would leave them indistinguishable at high resolution, which is exactly
    # what the first version of this test got wrong.
    assert got == pytest.approx(math.cos(math.pi / n_per_lam), abs=5e-4)
    assert abs(got - math.cos(2 * math.pi / n_per_lam)) > 5e-4, (
        "it now matches the analytic value for a wall at its nominal position. Has that half cell on "
        "the high wall been fixed? If so, change this test's expectation to cos(2*pi/N) and update "
        "the comment in boundaries.py")


@pytest.mark.parametrize("n_per_lam", [40, 80])
def test_pmc_low_wall_antinode_at_nominal_position(n_per_lam):
    """A PMC low wall has no such half cell: ``mpv[0]=-1`` puts H=0 exactly on node 0.

    The first sample point of the low-end monitor sits **exactly on the wall**, unlike the high end
    which is one cell short, so the expectation here is the antinode itself, 1.0, rather than
    cos(2*pi/N).
    """
    ex = _profile(_sim(td.PMCBoundary(), td.PML(), n_per_lam,
                       -LZ / 2, 0.0, 1.2, "-"))
    got = ex[0] / ex.max()
    assert got == pytest.approx(1.0, abs=2e-3), (
        f"a PMC low wall should carry the antinode peak, measured {got:.4f}")
