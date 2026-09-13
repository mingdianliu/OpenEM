# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
# -*- coding: utf-8 -*-
"""Injecting a CFS alpha when the PML band contains ADE pole cells (boundaries.PML_DISP_ALPHA).

Background: a CPML psi recursion with alpha=0 diverges at the **scheme level** over a long run
against ADE poles in the same cells. A term-by-term whole-domain CPU replay measured identical
fp32 and fp64 trajectories growing exponentially at about 6.7e-4 per step, with the mode inside the
PML band in the low-frequency region around omega*dt = 0.0084. Removing the poles from the band,
or injecting a CFS alpha normalized like sigma, made it entirely stable. The fix: inject
alpha = PML_DISP_ALPHA * avg_speed/(eta0*dl) * step into a psi-PML face only when its layer band
genuinely contains pole cells and the user's alpha_max is 0.

The criteria:
1. a face with poles in its band gets coefficients bitwise equal to _cfs_coeffs, and they really do
   decay more than the originals;
2. a face of the same scene with no poles in its band is bitwise identical to cpml.build;
3. a non-dispersive scene leaves every face bitwise unchanged;
4. a dispersive structure shrunk into the domain center, out of the PML band, leaves every face
   bitwise unchanged;
5. setting the environment variable to 0 restores the old behaviour bitwise;
6. a face where the user gave an explicit CFS alpha is copied through unchanged, with no injection.

CPU only; no GPU needed.
"""
from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import cpml
from openem.scene import boundaries
from openem.scene import from_simulation

FREQ0 = 2e14
#: A lossless Lorentz model with its pole at 500 THz and the working band at 200 THz, where eps is
#: real and lossless.
LORENTZ = td.Lorentz(eps_inf=2.0, coeffs=[(1.5, 5e14, 0.0)])


def _sim(slab_medium, slab_z=(0.35, 2.0), pml=None):
    """Air background, a plane wave along z, and a slab that either does or does not cross the z PML
    band.
    """
    zb = pml if pml is not None else td.Boundary.pml(num_layers=12)
    return td.Simulation(
        size=(0.1, 0.1, 1.2), grid_spec=td.GridSpec.uniform(dl=0.025),
        subpixel=False,
        structures=[td.Structure(
            geometry=td.Box.from_bounds(
                (-td.inf, -td.inf, slab_z[0]), (td.inf, td.inf, slab_z[1])),
            medium=slab_medium)],
        sources=[td.PlaneWave(
            center=(0, 0, -0.4), size=(td.inf, td.inf, 0), direction="+",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.2 * FREQ0))],
        monitors=[td.FluxMonitor(center=(0, 0, 0.0), size=(td.inf, td.inf, 0),
                                 freqs=[FREQ0], name="T")],
        run_time=2e-13,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.periodic(), y=td.Boundary.periodic(), z=zb))


def _plain_and_cfs(sim, sc, side):
    """The expected coefficients for that face without and with injection, taking parameters by the
    same route as build_pml.
    """
    params = boundaries._pml_params(sim, 2, side)
    axis = sc.grid.axes[2]
    dl = float(axis.dl[0] if side == "lo" else axis.dl[-1])
    speed = boundaries._avg_speed(sc.eps_ez, 2, side, params.num_layers)
    plain = cpml.build(params, dl=dl, side=side, dt=float(sc.dt), avg_speed=speed)
    cfs = boundaries._cfs_coeffs(params, dl, side, float(sc.dt), speed,
                                 boundaries.PML_DISP_ALPHA)
    return plain, cfs


def _same(c1, c2) -> bool:
    return all(np.array_equal(getattr(c1, f), getattr(c2, f))
               for f in ("a_E", "b_E", "inv_kappa_E", "a_H", "b_H", "inv_kappa_H"))


def test_gated_face_gets_cfs_alpha():
    """The slab crosses the z-hi PML band, so that face gets a CFS alpha and really does decay more."""
    sim = _sim(LORENTZ)
    sc = from_simulation(sim)
    assert sc.dispersion is not None and sc.dispersion.n_entry > 0
    plain, cfs = _plain_and_cfs(sim, sc, "hi")
    got = sc.pml[(2, "hi")]
    assert _same(got, cfs), "a face with poles in its band should receive the CFS coefficients"
    assert not _same(got, plain), "injection must actually change something, or the fix did nothing"
    # alpha only adds decay: b is at most the original in every layer, and strictly smaller in at
    # least one
    assert np.all(got.b_E <= plain.b_E) and np.any(got.b_E < plain.b_E)
    # sigma and kappa are untouched
    assert np.array_equal(got.inv_kappa_E, plain.inv_kappa_E)


def test_pole_free_face_unchanged():
    """The z-lo face of the same scene, whose band is all air, is bitwise unchanged."""
    sim = _sim(LORENTZ)
    sc = from_simulation(sim)
    plain, _cfs = _plain_and_cfs(sim, sc, "lo")
    assert _same(sc.pml[(2, "lo")], plain)


def test_nondispersive_scene_unchanged():
    """A non-dispersive scene, the same geometry with a static eps, leaves both faces bitwise
    unchanged.
    """
    sim = _sim(td.Medium(permittivity=2.0))
    sc = from_simulation(sim)
    assert sc.dispersion is None
    for side in ("lo", "hi"):
        plain, _cfs = _plain_and_cfs(sim, sc, side)
        assert _same(sc.pml[(2, side)], plain)


def test_center_confined_dispersive_unchanged():
    """A dispersive slab shrunk into the domain center, entering no PML band, leaves both faces
    bitwise unchanged.
    """
    sim = _sim(LORENTZ, slab_z=(-0.2, 0.2))
    sc = from_simulation(sim)
    assert sc.dispersion is not None and sc.dispersion.n_entry > 0
    for side in ("lo", "hi"):
        plain, _cfs = _plain_and_cfs(sim, sc, side)
        assert _same(sc.pml[(2, side)], plain)


def test_env_override_zero_restores_old(monkeypatch):
    """OPENEM_PML_DISP_ALPHA=0 is bitwise identical to the old behaviour."""
    monkeypatch.setenv("OPENEM_PML_DISP_ALPHA", "0")
    sim = _sim(LORENTZ)
    sc = from_simulation(sim)
    plain, _cfs = _plain_and_cfs(sim, sc, "hi")
    assert _same(sc.pml[(2, "hi")], plain)


def test_user_cfs_alpha_untouched():
    """A face where the user gave a non-zero alpha_max copies the tidy3d parameters through
    unchanged, with no injection.
    """
    user = td.Boundary.pml(num_layers=12)
    user = td.Boundary(
        minus=td.PML(num_layers=12, parameters=td.PMLParams(
            sigma_order=3, sigma_min=0.0, sigma_max=1.5,
            kappa_order=3, kappa_min=1.0, kappa_max=3.0,
            alpha_order=1, alpha_min=0.0, alpha_max=0.9)),
        plus=td.PML(num_layers=12, parameters=td.PMLParams(
            sigma_order=3, sigma_min=0.0, sigma_max=1.5,
            kappa_order=3, kappa_min=1.0, kappa_max=3.0,
            alpha_order=1, alpha_min=0.0, alpha_max=0.9)))
    sim = _sim(LORENTZ, pml=user)
    sc = from_simulation(sim)
    plain, _cfs = _plain_and_cfs(sim, sc, "hi")   # plain already carries the user's raw alpha
    assert _same(sc.pml[(2, "hi")], plain)
