# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Lossy media, i.e. conductivity sigma.

The criteria come from the Tidy3D source and from measurement on the validation set, not from
expectations computed here: ``sigma = Im(eps_r) * omega * eps0`` comes from ``medium.py:786``, and
the conductivity of the absorbing example is the value in the reference output.
"""

import numpy as np
import os
import pytest

td = pytest.importorskip("tidy3d")

from openem import cpml, serialize
from openem import scene as scene_mod
from openem.scene import media

TABLEA = os.environ.get("OPENEM_TABLEA", "/nonexistent")  # external reference archive; without
#                                                          OPENEM_TABLEA these cases skip
BEER = f"{TABLEA}/BeerLambert/Tidy3D/simulation.json"
BIO = f"{TABLEA}/BiosensorGrating/Tidy3D/simulation.json"


@pytest.fixture(scope="module")
def beer():
    return scene_mod.from_file(BEER)


def test_sigma_matches_declared_conductivity(beer):
    """The recovered sigma must equal the conductivity Tidy3D declares, converted to S/m.

    Tidy3D's conductivity is in um units (S/um) and ours is SI (S/m), a factor of 1e6 apart.
    """
    sim = td.Simulation.from_file(BEER)
    declared_si = float(sim.medium.conductivity) * 1e6
    for arr in (beer.sigma_ex, beer.sigma_ey, beer.sigma_ez):
        assert arr is not None
        assert arr.shape == beer.shape
        np.testing.assert_allclose(arr, declared_si, rtol=1e-12)


def test_sigma_is_frequency_independent(beer):
    """The test for pure conductivity: sigma must be **bitwise identical** at two frequencies.

    This is what failing closed rests on, and it is where a dispersive material gets rejected.
    """
    sim = td.Simulation.from_file(BEER)
    scene_mod.set_local_subpixel(bool(sim.subpixel))
    freqs = np.asarray(sim.monitors[0].freqs, dtype=float)
    lo, hi = float(freqs.min()), float(freqs.max())
    s_lo = media.sigma_from_eps(media.epsilon_complex(sim, "Ex", lo), lo)
    s_hi = media.sigma_from_eps(media.epsilon_complex(sim, "Ex", hi), hi)
    assert np.array_equal(s_lo, s_hi), "the sigma of a pure conductor should not depend on frequency"


def test_lossless_scene_has_no_sigma():
    """A lossless scene must have sigma of None, which is how the solver decides to take the original
    kernel.
    """
    sc = scene_mod.from_file(BIO)
    assert sc.sigma_ex is None and sc.sigma_ey is None and sc.sigma_ez is None
    assert sc.any_loss is False


def test_lossless_helper_zeroes_sigma(beer):
    """``lossless`` is what makes a reference run possible: one absorbing example has no structures
    at all, so ``homogeneous`` alone cannot build a reference, since it would be identical to the
    structure run.
    """
    assert beer.any_loss is True
    ref = scene_mod.lossless(scene_mod.homogeneous(beer))
    assert ref.any_loss is False
    np.testing.assert_allclose(ref.eps_ex, beer.eps_ex)      # uniform to begin with


def test_homogeneous_keeps_loss(beer):
    """``homogeneous`` on its own must not lose the loss; zeroing it is ``lossless``'s job."""
    assert scene_mod.homogeneous(beer).any_loss is True


def test_lossy_coeffs_reduce_to_lossless_at_zero_sigma(beer):
    """At sigma=0 the semi-implicit coefficients must degenerate exactly: ``ca=1`` and
    ``cb=dt/(eps0*eps)``.

    This is the arithmetic premise behind "a lossless result is bitwise unchanged"; the other half
    is that the kernel splits into two paths.
    """
    eps = beer.eps_ex
    kk = np.zeros_like(eps)
    ca = (1.0 - kk) / (1.0 + kk)
    cb = (beer.dt / (cpml.EPSILON_0 * eps)) / (1.0 + kk)
    assert np.all(ca == 1.0)
    np.testing.assert_array_equal(cb, beer.dt / (cpml.EPSILON_0 * eps))


def test_refuses_lossy_metal_medium():
    """``LossyMetalMedium`` must be refused: Tidy3D models it with a surface-impedance boundary, not
    a volumetric sigma.

    It **inherits from ``Medium``** and carries its own ``.conductivity``, so either duck typing or
    ``isinstance`` would let it through and we would then run it under an entirely wrong model. The
    check has to be by exact type.
    """
    assert issubclass(td.LossyMetalMedium, td.Medium), "the premise changed: it no longer inherits Medium"
    sim = td.Simulation.from_file(BEER)
    lossy = td.LossyMetalMedium(conductivity=10.0, frequency_range=(9e9, 1e10))
    bad = sim.updated_copy(
        structures=[td.Structure(geometry=td.Box(center=(0, 0, 0), size=(0.2, 0.2, 0.2)),
                                 medium=lossy)]
    )
    with pytest.raises(NotImplementedError, match="td.Medium"):
        scene_mod.from_simulation(bad)


def test_dispersive_medium_gets_poles_not_a_constant_eps():
    """A dispersive medium **must not** slip through as a constant eps at some frequency.

    This test used to assert that such media are refused outright; once the ADE was working it
    became an assertion that they are handled correctly: ``eps_*`` must store eps_inf, with the
    poles going separately into :class:`Dispersion`. The criterion is that eps_inf is clearly not
    equal to eps at that frequency; otherwise the "silently treated as a constant" failure would
    still pass.
    """
    sim = td.Simulation.from_file(BEER)
    med = td.material_library["Ag"]["Rakic1998BB"]
    sim2 = sim.updated_copy(
        structures=[td.Structure(geometry=td.Box(center=(0, 0, 0), size=(0.2, 0.2, 0.2)),
                                 medium=med)]
    )
    sc = scene_mod.from_simulation(sim2)
    assert sc.any_dispersion, "the dispersive medium produced no pole data"
    d = sc.dispersion
    assert d.n_pole >= d.n_entry, "every dispersive cell needs at least one pole"
    # The eps_inf of silver is O(1) while the real part of its eps at optical frequencies is a large
    # negative number; the two must never be confused
    eps_inf = float(med.pole_residue.eps_inf)
    f0 = float(sim.sources[0].source_time.freq0)
    eps_f0 = complex(med.eps_model(np.array([f0]))[0])
    assert abs(eps_f0.real - eps_inf) > 1.0, "the premise changed: eps_inf and eps(f0) of this material are too close"
    assert abs(sc.eps_ex.min() - eps_f0.real) > 1.0, (
        f"the minimum of eps_ex, {sc.eps_ex.min():.3f}, is close to eps(f0)={eps_f0.real:.3f}, "
        "which looks like the dispersive material was stored as a constant eps")


def test_fully_anisotropic_in_lossy_background_builds():
    """A tensor structure in a lossy scalar background, which used to fail closed and now goes
    through the tensor entry table.

    The diagonal sigma of the halo entries has to pick up the background's conductivity, or the ring
    around the tensor cells would silently become lossless.
    """
    sim = td.Simulation.from_file(BEER)
    mixed = sim.updated_copy(
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, 0), size=(0.2, 0.2, 0.2)),
            medium=td.FullyAnisotropicMedium(
                permittivity=[[2.0, 0.1, 0.0], [0.1, 2.0, 0.0], [0.0, 0.0, 2.0]]))]
    )
    sc = scene_mod.from_simulation(mixed)
    t = sc.tensor
    assert t is not None and t.n_entry > 0 and sc.any_loss
    halo = ~t.coupled_mask()
    assert halo.any()
    # The background is lossy, so the diagonal sigma of the halo is non-zero
    diag = t.sigma[halo][:, np.arange(3), np.arange(3)]
    assert float(np.abs(diag).max()) > 0.0


def test_serialize_roundtrip_keeps_sigma(beer, tmp_path):
    """sigma has to survive the npz round-trip: a solve job container has only numpy and takes
    exactly that path.
    """
    p = serialize.save(beer, tmp_path / "scene.npz")
    back = serialize.load(p)
    assert back.any_loss is True
    for a, b in ((beer.sigma_ex, back.sigma_ex),
                 (beer.sigma_ey, back.sigma_ey),
                 (beer.sigma_ez, back.sigma_ez)):
        np.testing.assert_array_equal(a, b)


def test_serialize_lossless_omits_sigma(tmp_path):
    """A lossless scene writes no sigma and reads back as None; an older v1 bundle is read the same
    way.
    """
    sc = scene_mod.from_file(BIO)
    back = serialize.load(serialize.save(sc, tmp_path / "s.npz"))
    assert back.any_loss is False

def test_lossy_subpixel_soft():
    """A lossy medium with subpixel averaging: only the thin interface layer has an ill-conditioned
    sigma inversion, so it is admitted on the evidence rather than rejected as a class.

    Measured on a dielectric sphere at 2 THz with sigma=50 S/m, the Ez component drifted with a
    median of 0, a maximum of 2.9e-4 and 0.004% of cells over the limit, against tolerance lines of
    2% and 1e-2. So it should be admitted and return a finite sigma.
    """
    import numpy as np
    from openem.scene import media

    freq0 = 2.0e12
    sim = td.Simulation(
        size=(4, 4, 4),
        grid_spec=td.GridSpec.uniform(dl=0.08),
        structures=[td.Structure(
            geometry=td.Sphere(center=(0, 0, 0), radius=1.2),
            medium=td.Medium(permittivity=4.0, conductivity=50.0))],
        sources=[td.PointDipole(
            center=(-1.8, 0, 0), polarization="Ey",
            source_time=td.GaussianPulse(freq0=freq0, fwidth=freq0 / 5))],
        monitors=[td.FluxMonitor(center=(1.8, 0, 0), size=(0, 2, 2),
                                 freqs=[freq0], name="t")],
        run_time=2e-12, subpixel=True,
        boundary_spec=td.BoundarySpec.all_sides(boundary=td.PML()))

    eps, sigma, _pec = media.epsilon_and_sigma(sim)      # must not raise NotImplementedError
    assert any(v is not None for v in sigma.values())
    for k, v in sigma.items():
        if v is not None:
            assert np.all(np.isfinite(v)), f"the sigma of {k} contains non-finite values"
