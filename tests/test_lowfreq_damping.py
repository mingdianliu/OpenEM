# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Extra damping for low-frequency lossless poles (the dt branch of poles.stabilize_lossless).

The contract: by default (env unset, no dt passed, or the constant at 0) nothing changes bitwise.
With the env set, only poles that are both lossless and satisfy |a|*dt < LOSSLESS_LOWFREQ_ADT_MAX
move, and after the injection |A| < 1. Lossy poles, non-low-frequency poles and residues are
always left bitwise untouched.
"""
import numpy as np
import pytest
import tidy3d as td

from openem.scene.poles import (
    LOSSLESS_AXIS_RTOL, LOSSLESS_LOWFREQ_ADT_MAX, pole_residue,
    split_dc_pole, stabilize_lossless, trap_coeffs)

DT = 4.8366e-17     # dt of the probe grid


def _kept(med):
    return split_dc_pole(pole_residue(med), DT)[0]


SIO2 = td.material_library["SiO2"]["Palik_Lossless"]
CSI = td.material_library["cSi"]["Palik_Lossless"]
AU = td.material_library["Au"]["JohnsonChristy1972"]


def test_default_bitwise_identical():
    """With the env unset, the new signature (taking dt) matches the old one bitwise."""
    for med in (SIO2, CSI, AU):
        kept = _kept(med)
        assert stabilize_lossless(kept, dt=DT) == stabilize_lossless(kept)


def test_no_dt_disables_lowfreq(monkeypatch):
    """Env set but no dt passed: the low-frequency branch stays off, since the test needs |a|*dt."""
    monkeypatch.setenv("OPENEM_LOSSLESS_LOWFREQ_DAMPING_REL", "0.3")
    kept = _kept(SIO2)
    assert stabilize_lossless(kept) == stabilize_lossless(kept, dt=None)
    monkeypatch.delenv("OPENEM_LOSSLESS_LOWFREQ_DAMPING_REL")
    ref = stabilize_lossless(kept)
    monkeypatch.setenv("OPENEM_LOSSLESS_LOWFREQ_DAMPING_REL", "0.3")
    assert stabilize_lossless(kept) == ref


def test_only_lowfreq_lossless_poles_move(monkeypatch):
    monkeypatch.setenv("OPENEM_LOSSLESS_LOWFREQ_DAMPING_REL", "0.3")
    for med in (SIO2, CSI, AU):
        kept = _kept(med)
        base = stabilize_lossless(kept, dt=None)     # uniform gamma_rel only
        out = stabilize_lossless(kept, dt=DT)
        n_moved = 0
        for (qb, rb), (qo, ro), (qk, rk) in zip(base, out, kept):
            assert ro == rb                          # residues never move
            lossless = abs(qk.real) <= LOSSLESS_AXIS_RTOL * abs(qk)
            low = abs(qk) * DT < LOSSLESS_LOWFREQ_ADT_MAX
            if lossless and low:
                assert qo == complex(-0.3 * abs(qk), qk.imag)
                n_moved += 1
            else:
                assert qo == qb                      # everything else stays bitwise identical
        if med is SIO2:
            assert n_moved == 1                      # exactly the low-frequency pole
        if med is AU:
            assert n_moved == 0                      # a lossy material leaves the whole table alone


def test_damped_pole_strictly_stable(monkeypatch):
    """After the injection the trapezoidal |A| < 1, and decays more than a uniform gamma_rel=1e-6."""
    monkeypatch.setenv("OPENEM_LOSSLESS_LOWFREQ_DAMPING_REL", "0.3")
    kept = _kept(SIO2)
    am1, _ = trap_coeffs(stabilize_lossless(kept, dt=DT), DT)
    am1_0, _ = trap_coeffs(stabilize_lossless(kept), DT)
    a = np.abs(am1 + 1.0)
    a0 = np.abs(am1_0 + 1.0)
    assert np.all(a < 1.0)
    assert np.min(a0 - a) >= 0.0 and np.max(a0 - a) > 1e-4


def test_env_zero_is_off(monkeypatch):
    monkeypatch.setenv("OPENEM_LOSSLESS_LOWFREQ_DAMPING_REL", "0")
    kept = _kept(SIO2)
    assert stabilize_lossless(kept, dt=DT) == stabilize_lossless(kept)
