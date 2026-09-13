# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""The CPML profile lined up against Tidy3D.

What makes this file worth having: **it does not restate our own formula, it calls tidy3d's own
``create_sfactor_f`` / ``create_sfactor_b`` / ``s_value`` and compares point by point.**
The moment Tidy3D changes a convention, this goes red.
"""

import numpy as np
import pytest

from openem import cpml

# The whole value of this file is comparing point by point against tidy3d's own implementation, so
# without tidy3d the file is skipped entirely: a GPU job container may not have tidy3d installed,
# and there only the tests that do not depend on it run.
pytest.importorskip("tidy3d", reason="this file cross-checks against tidy3d's source")

from tidy3d.components.mode.derivatives import create_sfactor_b, create_sfactor_f
from tidy3d.constants import C_0 as TD_C_0
from tidy3d.constants import EPSILON_0 as TD_EPS_0
from tidy3d.constants import ETA_0 as TD_ETA_0
from tidy3d.constants import MU_0 as TD_MU_0

#: tidy3d works in µm, so lengths differ by 1e6
UM_PER_M = 1e6


def test_constants_match_tidy3d():
    """Converted, our SI constants must be exactly equal to tidy3d's µm-based ones.

    This locks the decision to work in SI metres internally: change a constant and it goes red at
    once.
    """
    assert cpml.ETA_0 == TD_ETA_0, "impedance is unit-system independent, must be bitwise equal"
    assert cpml.EPSILON_0 == pytest.approx(TD_EPS_0 * UM_PER_M, rel=1e-15)
    assert cpml.MU_0 == pytest.approx(TD_MU_0 * UM_PER_M, rel=1e-15)
    assert cpml.C_0 == pytest.approx(TD_C_0 / UM_PER_M, rel=1e-15)


def _our_sfactor(params, N, n_pml, dl_um, omega, avg_speed, at):
    """Assemble an array from our profile with the same length as tidy3d's create_sfactor_*."""
    s = np.ones(N, dtype=np.complex128)
    dl_m = dl_um / UM_PER_M
    for side, sl in (("lo", slice(0, n_pml)), ("hi", slice(N - n_pml, N))):
        sigma, kappa, _alpha = cpml.profile(params, dl_m, side, at, avg_speed)
        # tidy3d's s = kappa + 1j*sigma/(omega*eps0), where sigma and eps0 are both in µm units
        sigma_um = sigma / UM_PER_M
        s[sl] = kappa + 1j * sigma_um / (omega * TD_EPS_0)
    return s


@pytest.mark.parametrize("n_pml", [8, 12, 16])
@pytest.mark.parametrize("dl_um", [0.025, 0.01, 0.1])
def test_sfactor_matches_tidy3d(n_pml, dl_um):
    """The s-factor must agree point by point at both the E position (backward) and the H position
    (forward).

    **Mind what this test actually proves.** ``create_sfactor_b/f`` comes from
    ``tidy3d/components/mode/derivatives.py``, which is the PML of the **mode solver**, not of the
    FDTD solver. The FDTD PML lives in the hosted solver, is invisible to the client, and **cannot
    be checked**.

    So what is locked here is "our profile formula == Tidy3D's mode-solver profile formula". It
    does **not** license the claim that FDTD uses the same normalization. Measurement says the
    opposite: at grazing incidence σ has to be scaled up 4x to reproduce the hosted FDTD result.
    """
    N = 4 * n_pml
    omega = 2 * np.pi * 3.612e14
    avg_speed = 0.75  # BiosensorGrating background n=1.333
    # Use tidy3d s_value's own defaults, otherwise the comparison means nothing
    params = cpml.PMLParams(
        num_layers=n_pml,
        sigma_max=2.0,
        sigma_order=3,
        kappa_min=1.0,
        kappa_max=3.0,
        kappa_order=3,
        alpha_max=0.0,
    )
    dls = np.full(N, dl_um)

    theirs_b = create_sfactor_b(omega, dls, N, n_pml, True, (avg_speed, avg_speed))
    ours_b = _our_sfactor(params, N, n_pml, dl_um, omega, avg_speed, "E")
    np.testing.assert_allclose(ours_b, theirs_b, rtol=1e-12, atol=0)

    theirs_f = create_sfactor_f(omega, dls, N, n_pml, True, (avg_speed, avg_speed))
    ours_f = _our_sfactor(params, N, n_pml, dl_um, omega, avg_speed, "H")
    np.testing.assert_allclose(ours_f, theirs_f, rtol=1e-12, atol=0)


def test_kappa_endpoints_asymmetric():
    """Record Tidy3D's min/max asymmetry: not a bug, but a product of the backward derivative's
    staggering.

    min side, outermost step is exactly 1.0 -> κ is exactly kappa_max
    max side, innermost step is 0           -> κ = kappa_min, σ = 0, that cell does nothing
    """
    p = cpml.PMLParams(num_layers=12)
    s_lo = cpml.steps(12, "lo", "E")
    s_hi = cpml.steps(12, "hi", "E")
    assert s_lo[0] == 1.0
    assert s_hi[0] == 0.0
    assert s_hi[-1] == pytest.approx(11 / 12)

    _, kappa_lo, _ = cpml.profile(p, 2.5e-8, "lo", "E")
    _, kappa_hi, _ = cpml.profile(p, 2.5e-8, "hi", "E")
    assert kappa_lo[0] == pytest.approx(p.kappa_max)
    assert kappa_hi[0] == pytest.approx(p.kappa_min)


def test_kappa_grading_starts_at_one():
    """κ has to grow smoothly from 1.0 at the inner side: any jump at the interface reflects."""
    p = cpml.PMLParams(num_layers=12, kappa_min=1.0, kappa_max=3.0, kappa_order=3)
    _, kappa, _ = cpml.profile(p, 2.5e-8, "lo", "E")
    inner_to_outer = kappa[::-1]  # on the lo side increasing global index = outside to inside
    assert inner_to_outer[0] == pytest.approx(1.0 + 2.0 * (1 / 12) ** 3)
    assert np.all(np.diff(inner_to_outer) > 0), "κ must increase monotonically"
    assert inner_to_outer[-1] == pytest.approx(3.0)


def test_recursion_coeffs_sigma_zero_gives_identity():
    """In a layer with σ=0 the ψ recursion must be the identity (b=1, a=0), never nan."""
    sigma = np.array([0.0, 1e4])
    kappa = np.array([1.0, 2.0])
    alpha = np.zeros(2)
    a, b = cpml.recursion_coeffs(sigma, kappa, alpha, dt=4.77e-17)
    assert a[0] == 0.0 and b[0] == 1.0
    assert np.isfinite(a).all() and np.isfinite(b).all()
    assert 0.0 < b[1] < 1.0, "in a lossy layer b must be in (0,1), otherwise ψ diverges"


def test_recursion_coeffs_alpha_zero_closed_form():
    """At α=0, a = (b-1)/κ: the algebra of FDTDX pml.py:124 simplified for α=0."""
    sigma = np.array([1e3, 1e5])
    kappa = np.array([1.5, 3.0])
    a, b = cpml.recursion_coeffs(sigma, kappa, np.zeros(2), dt=4.77e-17)
    np.testing.assert_allclose(a, (b - 1.0) / kappa, rtol=1e-13)


def test_sigma_min_nonzero_refuses():
    """In Tidy3D's s_value, σ always starts from 0; σ_min != 0 has no basis, so fail closed."""
    p = cpml.PMLParams(num_layers=12, sigma_min=0.1)
    with pytest.raises(NotImplementedError, match="sigma_min"):
        cpml.profile(p, 2.5e-8, "lo", "E")


def test_build_shapes():
    p = cpml.PMLParams(num_layers=12)
    c = cpml.build(p, dl=2.5e-8, side="lo", dt=4.766437173827506e-17, avg_speed=0.75)
    assert c.num_layers == 12
    for arr in (c.a_E, c.b_E, c.inv_kappa_E, c.a_H, c.b_H, c.inv_kappa_H):
        assert arr.shape == (12,)
        assert np.isfinite(arr).all()
    # inv_kappa must lie in [1/kappa_max, 1/kappa_min] = [1/3, 1]
    assert c.inv_kappa_E.min() >= 1 / 3 - 1e-12
    assert c.inv_kappa_E.max() <= 1.0 + 1e-12
