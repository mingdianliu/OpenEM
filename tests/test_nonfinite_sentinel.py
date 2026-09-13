# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Divergence sentinel at the solve-data exit (solver._warn_nonfinite).

It warns and does not block: NaN or Inf get named at the exit before they escape, so the source is
visible even when the job log is crowded. On healthy data it prints nothing and behaviour is
bitwise unchanged. A "1 validation error for Box" seen deep inside a notebook was exactly this:
a NaN escaping silently and surfacing much later through amps, np.angle, np.interp and finally
td.Box(size=nan).
"""

from __future__ import annotations

import numpy as np

import pytest

# This module imports openem.solver at top level, and that imports cupy at top level, so a CPU-only
# machine would crash during collection.
pytest.importorskip("cupy", reason="needs CUDA and cupy")

from openem.solver import _warn_nonfinite


def test_finite_data_is_silent(capsys):
    _warn_nonfinite({"m": np.ones((2, 3), complex)}, {"f": np.zeros(4)})
    assert capsys.readouterr().out == ""


def test_nan_and_inf_are_reported_not_raised(capsys):
    bad = np.ones(5, complex)
    bad[1] = complex("nan")
    worse = np.ones((2, 2))
    worse[0, 0] = np.inf
    _warn_nonfinite({"amps": bad}, {"vol": worse})   # must not raise
    out = capsys.readouterr().out
    assert "mon:amps" in out and "1/5" in out
    assert "fld:vol" in out and "1/4" in out


def test_data_attr_and_empty_are_handled(capsys):
    class _P:  # stands in for an object in res.phasors that carries .data
        data = np.array([np.nan + 0j])

    _warn_nonfinite({"p": _P()}, {"e": np.zeros(0)})
    out = capsys.readouterr().out
    assert "mon:p" in out and "1/1" in out
