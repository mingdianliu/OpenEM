# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
# -*- coding: utf-8 -*-
"""Unit tests for the epsilon(omega) fast-path decision.

This is a **performance** switch, but deciding it wrongly changes results silently: a dispersive
medium treated as non-dispersive collapses epsilon(omega) to a single-frequency value and the
geometry gradients fall apart with it. So the decision is pinned down:
* non-dispersive with zero conductivity: True, one frequency point suffices
* dispersive (Sellmeier, PoleResidue, Lorentz and so on): False
* with conductivity: False, since epsilon(omega) = eps_r + i*sigma/(omega*eps0) depends on omega
* an unrecognized type: False, preferring slow over wrong
"""
from __future__ import annotations

import numpy as np
import pytest
import tidy3d as td

from openem.nb import backend as cnb


def _sim(medium):
    return td.Simulation(
        size=(2, 2, 2),
        grid_spec=td.GridSpec.uniform(dl=0.2),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, 0), size=(1, 1, 1)), medium=medium)],
        sources=[td.PointDipole(
            center=(0, 0, -0.8), polarization="Ey",
            source_time=td.GaussianPulse(freq0=2e14, fwidth=2e13))],
        monitors=[], run_time=1e-13,
        boundary_spec=td.BoundarySpec.all_sides(td.Periodic()))


def test_plain_medium_is_freq_independent():
    assert cnb._eps_is_freq_independent(_sim(td.Medium(permittivity=4.0)))


def test_conductive_medium_is_not():
    """With conductivity, epsilon(omega) = eps_r + i*sigma/(omega*eps0) and depends on frequency."""
    assert not cnb._eps_is_freq_independent(
        _sim(td.Medium(permittivity=4.0, conductivity=0.3)))


@pytest.mark.parametrize("med", [
    td.Sellmeier(coeffs=[(1.0, 0.1)]),
    td.Lorentz(eps_inf=1.0, coeffs=[(1.0, 2e14, 1e13)]),
    td.material_library["Si3N4"]["Philipp1973Sellmeier"],
])
def test_dispersive_media_are_not(med):
    assert not cnb._eps_is_freq_independent(_sim(med))


def test_custom_medium_without_conductivity_is_freq_independent():
    """A CustomMedium is not itself dispersive; inverse-design regions are typically of this kind."""
    n = 5
    c = np.linspace(-0.4, 0.4, n)
    da = td.SpatialDataArray(
        np.full((n, n, n), 3.0), coords=dict(x=c, y=c, z=c))
    assert cnb._eps_is_freq_independent(_sim(td.CustomMedium(permittivity=da)))


def test_unknown_type_falls_back_to_slow_path():
    """An unrecognized type returns False: slow, but correct."""

    class _Weird:
        pass

    class _FakeScene:
        mediums = [_Weird()]

    class _FakeSim:
        scene = _FakeScene()

    assert not cnb._eps_is_freq_independent(_FakeSim())
