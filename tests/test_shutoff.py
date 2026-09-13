# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Early-termination criteria.

The focus is :class:`PhasorConvergence`, the only one of the three guaranteed to fire; see the
module docstring of openem/shutoff.py for why. A guided-mode resonance device has a Q high enough
that "wait for the field to decay" never holds within the nominal step count.
"""

import numpy as np
import pytest

from openem import shutoff


def test_field_decay_ratio_is_vs_running_peak():
    """The denominator is the maximum over **time**, not the current value. One of the four easy
    mistakes listed in docs/architecture.md.
    """
    fd = shutoff.FieldDecay(threshold=1e-6)
    assert not fd.update(400, 1.0)
    assert fd.ratio == 1.0
    assert not fd.update(800, 0.5)      # the peak is still 1.0
    assert fd.ratio == 0.5
    assert not fd.update(1200, 2.0)     # the peak rises to 2.0
    assert fd.ratio == 1.0
    assert fd.update(1600, 1e-7)        # 1e-7/2.0 < 1e-6
    assert fd.deepest == (1600, 5e-8)


def test_field_decay_phase_lottery_is_reproduced():
    """An instantaneous criterion can miss a deep trough inside an oscillation, which is exactly how
    the reference solver once stopped by luck.
    """
    fd = shutoff.FieldDecay(threshold=1e-6)
    # The envelope is 1e-3, but each sample lands at a different phase
    for i, m in enumerate([1.0, 1e-3, 8e-4, 1.2e-3, 5e-4, 9e-4]):
        assert not fd.update(400 * (i + 1), m)
    # Under the same envelope, hitting one deep trough fires it
    assert fd.update(2800, 1e-7)


def test_energy_decay_is_smooth():
    ed = shutoff.EnergyDecay(threshold=1e-6)
    assert not ed.update(400, 1.0)
    assert not ed.update(800, 1e-3)
    assert ed.update(1200, 1e-7)


def test_phasor_convergence_needs_patience():
    """One small change is not convergence; it takes patience consecutive ones, which prevents a
    false trigger at an inflection point.
    """
    pc = shutoff.PhasorConvergence(tol=1e-4, patience=3, window_steps=400)
    tau = 400.0  # tau equals the window length, so the extrapolation factor is 1 and only the
                 # patience logic is under test
    assert not pc.update(400, np.array([1.0, 2.0]), tau)          # the first call only records a baseline
    assert not pc.update(800, np.array([1.0, 2.0 + 1e-5]), tau)   # streak 1
    assert pc.streak == 1
    assert not pc.update(1200, np.array([1.0, 2.0 + 2e-5]), tau)  # streak 2
    assert not pc.update(1600, np.array([1.0, 2.1]), tau)         # a large change resets the count
    assert pc.streak == 0
    for _ in range(3):
        fired = pc.update(2000, np.array([1.0, 2.1]), tau)
    assert fired


def test_phasor_convergence_relative_not_absolute():
    """The relative change is scaled by the largest component, so an overall rescaling does not
    affect the verdict.
    """
    pc = shutoff.PhasorConvergence(tol=1e-4, patience=1, window_steps=400)
    pc.update(400, np.array([1e-12, 2e-12]), 400.0)
    pc.update(800, np.array([1e-12, 2e-12 * (1 + 1e-6)]), 400.0)
    assert pc.rel < 1e-4


def test_phasor_extrapolation_scales_with_tau():
    """The extrapolated change is the per-window change times tau over the window length. A large tau
    means a large remaining change and must not fire.

    This is the heart of the criterion: 1.2e-3 per window looks small, but at high Q it accumulates
    over hundreds of windows. Measured, going from step 68,400 to step 209,802 drifted T by about
    1e-2.
    """
    pc = shutoff.PhasorConvergence(tol=1e-3, patience=1, window_steps=400)
    pc.update(400, np.array([1.0]), None)
    # A per-window change of 1e-4, but tau is 100 window lengths, so the remainder is 1e-2, far above tol
    pc.update(800, np.array([1.0 + 1e-4]), 40000.0)
    # The relative change is scaled by the largest component: 1e-4 / 1.0001
    assert pc.per_window == pytest.approx(1e-4, rel=1e-3)
    assert pc.rel == pytest.approx(pc.per_window * 100, rel=1e-9)   # tau over window length = 100
    assert pc.streak == 0


def test_phasor_falls_back_without_tau():
    """Without a tau, i.e. when the decay is too slow to fit, it falls back to the per-window change."""
    pc = shutoff.PhasorConvergence(tol=1e-3, patience=1)
    pc.update(400, np.array([1.0]), None)
    pc.update(800, np.array([1.0 + 1e-4]), float("inf"))
    assert pc.rel == pc.per_window


def test_energy_decay_tau_fit():
    """Fitting tau: given a known exponential decay, it should recover the time constant."""
    ed = shutoff.EnergyDecay(threshold=0.0)
    tau_true = 5000.0
    for i in range(40):
        step = 400 * (i + 1)
        ed.update(step, np.exp(-step / tau_true))
    assert ed.tau == pytest.approx(tau_true, rel=1e-6)


def test_energy_decay_tau_inf_when_not_decaying():
    """With no decay, or with growth, it returns inf, meaning "far from converged" rather than
    "converged".
    """
    ed = shutoff.EnergyDecay(threshold=0.0)
    for i in range(40):
        ed.update(400 * (i + 1), 1.0)
    assert ed.tau == float("inf")


def test_policy_reports_which_fired():
    pol = shutoff.Policy(
        field_decay=shutoff.FieldDecay(threshold=1e-6),
        phasor=shutoff.PhasorConvergence(tol=1e-4, patience=1),
    )
    assert not pol.check(400, 1.0, 1.0, np.array([1.0]))
    assert pol.check(800, 1e-9, 1.0, np.array([1.0]))
    assert "field_decay" in pol.reason
    assert "phasor_converged" in pol.reason


def test_default_policy_has_all_three():
    pol = shutoff.default_policy(1e-6)
    assert pol.field_decay is not None
    assert pol.energy_decay is not None
    assert pol.phasor is not None
    assert pol.interval == 400
    assert pol.phasor.tol == 1e-3        # target tolerance
    assert pol.phasor.window_steps == 400
