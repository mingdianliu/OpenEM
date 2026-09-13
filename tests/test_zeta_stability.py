# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
# -*- coding: utf-8 -*-
"""Long-run boundedness in float32 of a marginally stable zeta pole, the +-i*omega_L of a lossless
material.

`dispersion.py` used to reject such poles on the grounds that it had never been verified. This
pins the conclusion down: stepped the way the kernel does it, ``P += (A-1)*P``, entirely in
complex64 and freely oscillating, it does not diverge over a long run.

The key design is that ``trap_coeffs`` stores **A-1** rather than A: at small omega*dt, A is about
1 + i*omega*dt, and storing A-1 puts the rounding error at the scale of omega*dt instead of 1.
This test also locks that design down: change it to store A and the rounding of |A|-1 jumps to
about 6e-8, which a long run will not survive.
"""
from __future__ import annotations

import numpy as np
import pytest

from openem.scene.poles import trap_coeffs

#: Step count. The longest run in the validation set is 248,863 steps; this takes half of that,
#: tightens the bound proportionally, and finishes in a few seconds.
N_STEPS = 50_000

#: Bound on the amplitude ratio. Measured, the worst case over 250k steps was 1.0016, at the extreme
#: omega*dt = 1; 50k steps corresponds to about 1.0003, leaving one digit of margin.
BOUND = 3e-3


@pytest.mark.parametrize("wdt", [0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0])
def test_marginally_stable_pole_bounded_in_fp32(wdt):
    """A pole with Re=0: after 50,000 fp32 recursion steps the amplitude does not run away."""
    am1, _b = trap_coeffs([(1j * wdt, 1.0 + 0j)], 1.0)   # dt=1, so q*dt = i*omega*dt
    assert abs(abs(1.0 + complex(am1[0])) - 1.0) < 1e-12, "in fp64, |A| should be exactly 1"

    m = np.complex64(am1[0])
    p = np.complex64(1.0 + 0j)
    for _ in range(N_STEPS):
        p = np.complex64(p + m * p)
    ratio = abs(complex(p))
    assert abs(ratio - 1.0) < BOUND, (
        f"omega*dt={wdt}: after {N_STEPS} steps the amplitude ratio is {ratio:.6f}, outside +-{BOUND}")


def test_storing_a_minus_one_is_what_makes_it_work():
    """Storing A-1 rounds far better than storing A. That is why the test above holds; do not change
    it back.
    """
    wdt = 0.05
    am1, _ = trap_coeffs([(1j * wdt, 1.0 + 0j)], 1.0)
    a64 = 1.0 + complex(am1[0])
    err_am1 = abs(1.0 + complex(np.complex64(am1[0]))) - 1.0   # storing A-1
    err_a = abs(complex(np.complex64(a64))) - 1.0              # storing A
    assert abs(err_am1) < abs(err_a) / 10, (
        f"the |A| error when storing A-1, {err_am1:.3e}, should be far below that of storing A, {err_a:.3e}")
