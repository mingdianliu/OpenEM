# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Diagonal-tensor anisotropic media (td.AnisotropicMedium).

**The criterion is bitwise identity, with zero tunable parameters**: an Ex-polarized plane wave
travelling along z has Ey = Ez = Hx = Hz ≡ 0 on the Yee grid at all times (the corresponding curl
components vanish term by term), so on its way through an
``AnisotropicMedium(xx=a, yy=b, zz=c)`` slab it can only see ε_xx, and the run has to be
**bitwise identical** to one through a ``Medium(permittivity=a)`` slab. The same holds for Ey
polarization and ε_yy. This also locks "the component mapping is not written backwards": get it
backwards and the whole ε is wrong, not just the rounding.

Three degrees of freedom have to be pinned down, otherwise the two runs are not comparable in the
first place and the test would not be measuring anisotropy at all:

* **dt**: Tidy3D takes dt from the smallest ``n_cfl`` in the scene. Putting the medium in a slab
  and leaving the background as vacuum makes vacuum set dt for both runs.
* **PML coefficients**: ``build_pml`` uses ``eps_ez`` to set the absorption profile, and the two
  scenes have different ε_zz. Making every boundary Periodic removes this term (one-way plane-wave
  injection does not depend on an absorbing boundary).
* **The source plane**: it has to be in vacuum. Tidy3D itself rejects a plane-wave source cutting
  through an anisotropic medium, and our fail-closed check rejects a source plane that lands
  inside a dispersive medium.
"""

from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import scene as scene_mod

FREQ0 = 2.0e14
EPS_XX, EPS_YY, EPS_ZZ = 2.0, 3.0, 4.0


def _slab_sim(medium, normal: str) -> "td.Simulation":
    """A slab in vacuum with normal ``normal`` (its faces aligned to the dl=0.05 grid lines), with
    a plane wave injected along that axis.

    With pol_angle=0 the E field lies along the first transverse axis: normal='z' -> Ex
    polarization, normal='x' -> Ey polarization.
    """
    ax = "xyz".index(normal)
    size = [0.4, 0.4, 0.4]
    size[ax] = 2.4
    slab_center = [0.0, 0.0, 0.0]
    slab_center[ax] = 0.3
    slab_size = [td.inf, td.inf, td.inf]
    slab_size[ax] = 0.6
    src_center = [0.0, 0.0, 0.0]
    src_center[ax] = -0.8
    src_size = [td.inf, td.inf, td.inf]
    src_size[ax] = 0.0
    mon_size = [0.2, 0.2, 0.2]
    mon_size[ax] = 0.4
    mon_center = [0.05, 0.0, 0.0]
    mon_center[ax] = 0.3
    return td.Simulation(
        size=tuple(size),
        grid_spec=td.GridSpec.uniform(dl=0.05),
        medium=td.Medium(permittivity=1.0),
        structures=[td.Structure(
            geometry=td.Box(center=tuple(slab_center), size=tuple(slab_size)),
            medium=medium)],
        sources=[td.PlaneWave(
            center=tuple(src_center), size=tuple(src_size), direction="+",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.3 * FREQ0))],
        monitors=[td.FieldMonitor(
            center=tuple(mon_center), size=tuple(mon_size),
            freqs=[FREQ0], name="probe")],
        run_time=2.0e-13,
        boundary_spec=td.BoundarySpec.all_sides(td.Periodic()),
        subpixel=False,
    )


def _run(sim, num_steps: int):
    from openem import solver
    from openem.device import Kernels

    sc = scene_mod.from_simulation(sim)
    return sc, solver.run(sc, num_steps=num_steps, use_shutoff=False,
                          kernels=Kernels(), verbose=False)


def _aniso():
    return td.AnisotropicMedium(xx=td.Medium(permittivity=EPS_XX),
                                yy=td.Medium(permittivity=EPS_YY),
                                zz=td.Medium(permittivity=EPS_ZZ))


def test_static_eps_maps_components():
    """After the scene is built, eps_ex/ey/ez on the slab cells have to equal xx/yy/zz, which
    locks all three mappings at once.

    The plane-wave tests only cover xx and yy (pol_angle=0 never reaches Ez polarization), so zz
    is locked here.
    """
    sc = scene_mod.from_simulation(_slab_sim(_aniso(), "z"))
    for arr, want in ((sc.eps_ex, EPS_XX), (sc.eps_ey, EPS_YY), (sc.eps_ez, EPS_ZZ)):
        assert sorted(np.unique(arr)) == [1.0, want]


def test_static_plane_wave_bit_identical():
    """An Ex plane wave through the anisotropic slab == the run with an isotropic
    ``Medium(ε_xx)``.

    The solver supports only one kind of plane wave, "along z, Ex polarization" (solver.py fails
    closed otherwise), so at run level only xx can be locked directly; yy/zz are locked by the
    blindness test below combined with the array-mapping test above.
    """
    pytest.importorskip("cupy")
    sc_a, r_a = _run(_slab_sim(_aniso(), "z"), num_steps=2000)
    sc_i, r_i = _run(_slab_sim(td.Medium(permittivity=EPS_XX), "z"), num_steps=2000)
    assert sc_a.dt == sc_i.dt, "the two runs have different dt, a bitwise comparison is meaningless"
    p_a, p_i = r_a.field_phasors["probe"], r_i.field_phasors["probe"]
    assert float(np.max(np.abs(p_i))) > 0, "the reference run is all zeros, the criterion is vacuous"
    np.testing.assert_array_equal(p_a, p_i)


def test_static_plane_wave_blind_to_yy_zz():
    """Swap yy/zz and the Ex plane-wave result has to stay bitwise unchanged: ε_yy/ε_zz do not
    leak into the Ex path."""
    pytest.importorskip("cupy")
    other = td.AnisotropicMedium(xx=td.Medium(permittivity=EPS_XX),
                                 yy=td.Medium(permittivity=7.0),
                                 zz=td.Medium(permittivity=9.0))
    sc_a, r_a = _run(_slab_sim(_aniso(), "z"), num_steps=2000)
    sc_b, r_b = _run(_slab_sim(other, "z"), num_steps=2000)
    assert sc_a.dt == sc_b.dt
    np.testing.assert_array_equal(r_a.field_phasors["probe"],
                                  r_b.field_phasors["probe"])


def test_uniaxial_dispersive_bit_identical():
    """Uniaxial dispersion (xx a PoleResidue, yy/zz constants) == the isotropic run with the same
    PoleResidue.

    An Ex plane wave only drives the pole recurrence on the xx path; whether yy/zz are constants
    or poles, their fields stay identically zero. So Ex/Hy have to be bitwise identical between
    the two runs, which locks the assembly of the "(medium, axis) -> pole table" mapping: hang the
    poles on the wrong component and the whole dispersive response is wrong.
    """
    pytest.importorskip("cupy")
    csi = td.material_library["cSi"]["Green2008"]
    pr = csi.pole_residue if hasattr(csi, "pole_residue") else csi
    uniax = td.AnisotropicMedium(xx=pr, yy=td.Medium(permittivity=2.0),
                                 zz=td.Medium(permittivity=3.0))
    sc_a, r_a = _run(_slab_sim(uniax, "z"), num_steps=2000)
    sc_i, r_i = _run(_slab_sim(pr, "z"), num_steps=2000)
    assert sc_a.dt == sc_i.dt
    # first match at the assembly level: the dispersion table of comp 0 is bitwise identical
    d_a, d_i = sc_a.dispersion, sc_i.dispersion
    m_a, m_i = d_a.comp == 0, d_i.comp == 0
    np.testing.assert_array_equal(d_a.cell[m_a], d_i.cell[m_i])
    np.testing.assert_array_equal(d_a.am1[np.repeat(m_a, np.diff(d_a.pole_ofs))],
                                  d_i.am1[np.repeat(m_i, np.diff(d_i.pole_ofs))])
    # the anisotropic run has poles only on comp 0 (yy/zz are constants)
    assert set(np.unique(d_a.comp)) == {0}
    p_a, p_i = r_a.field_phasors["probe"], r_i.field_phasors["probe"]
    assert float(np.max(np.abs(p_i))) > 0, "the reference run is all zeros, the criterion is vacuous"
    assert np.all(np.isfinite(p_a)), "the anisotropic run diverged"
    np.testing.assert_array_equal(p_a, p_i)


def test_fully_anisotropic_builds_tensor_entries():
    """FullyAnisotropicMedium now goes through the tensor entry table (it used to fail closed).

    Structural assertions: the entry table is non-empty, all three components are present, and a
    diagonal halo is laid down around the true tensor cells. The numerical criteria live in
    tests/test_tensor.py (diagonal degeneracy / rotation superposition identity / Faraday
    rotation).
    """
    r = np.array([[1, 0, 0], [0, 0.8, -0.6], [0, 0.6, 0.8]])
    perm = r @ np.diag([2.0, 3.0, 4.0]) @ r.T
    full = td.FullyAnisotropicMedium(permittivity=perm.tolist())
    sim = td.Simulation(
        size=(0.4, 0.4, 0.4),
        grid_spec=td.GridSpec.uniform(dl=0.05),
        medium=td.Medium(permittivity=1.0),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, 0), size=(0.2, 0.2, 0.1)), medium=full)],
        sources=[td.PointDipole(
            center=(0, 0, 0.15), polarization="Ez",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.3 * FREQ0))],
        monitors=[td.FieldMonitor(center=(0, 0, 0), size=(0.2, 0.2, 0),
                                  freqs=[FREQ0], name="probe")],
        run_time=2.0e-13,
        boundary_spec=td.BoundarySpec.all_sides(td.Periodic()),
        subpixel=False,
    )
    sc = scene_mod.from_simulation(sim)
    t = sc.tensor
    assert t is not None and t.n_entry > 0
    assert set(np.unique(t.comp)) == {0, 1, 2}
    cpl = t.coupled_mask()
    assert cpl.any() and (~cpl).any()          # both the true tensor cells and the halo are there
    # the ε of a true tensor cell is the tensor that was passed in; the halo is purely diagonal
    np.testing.assert_allclose(t.eps[cpl][0], perm)
    off = ~np.eye(3, dtype=bool)
    assert np.abs(t.eps[~cpl][:, off]).max() == 0.0
    # on a tensor cell the scalar ε array holds the background (vacuum): the pass-through
    # overwrites it, it never enters the physics
    assert set(np.round(np.unique(sc.eps_ex), 9)) == {1.0}


# ---------------------------------------------------------------------------
# Two assembly problems the hBN scene exposed: collinear constant-medium columns, and
# serialization of a scene containing a dipole, a field monitor and mix entries
# ---------------------------------------------------------------------------

def _lorentz():
    # the real in-plane parameters of hBN (NanostructuredBoronNitride)
    return td.Lorentz(eps_inf=3.0,
                      coeffs=((1.0958877712259605, 41821047891000.0, 29979245800.0),))


def _hbn_like_sim(subpixel: bool) -> "td.Simulation":
    """Vacuum + a dispersive slab + a ``Medium(12)`` slab: the design matrix then holds two
    lossless constant columns at once.

    The dispersive slab is 0.55 thick, so its faces land **between** grid lines; with
    subpixel=True the normal (Ez) cells are a harmonic mean, which necessarily produces mix
    entries.
    """
    freq0 = 42.42e12
    return td.Simulation(
        size=(0.4, 0.4, 2.4),
        grid_spec=td.GridSpec.uniform(dl=0.1),
        medium=td.Medium(permittivity=1.0),
        structures=[
            td.Structure(geometry=td.Box(center=(0, 0, 0.0), size=(td.inf, td.inf, 0.55)),
                         medium=_lorentz()),
            td.Structure(geometry=td.Box(center=(0, 0, -0.7), size=(td.inf, td.inf, 0.6)),
                         medium=td.Medium(permittivity=12.0)),
        ],
        sources=[td.PointDipole(
            center=(0, 0, 0.9), polarization="Ez",
            source_time=td.GaussianPulse(freq0=freq0, fwidth=0.1 * freq0))],
        # The whole band has to sit **below** the resonance (41.8 THz): there Re ε > 0 and tidy3d
        # uses PolarizedAveraging, whereas above the resonance Re ε < 0 is treated as a metal and
        # switched to Staircasing, so ε_eff(f) on an interface cell jumps between the two
        # conventions and the mix model cannot fit it (and should not).
        # It also must not be too narrow: in a narrow band an arithmetic combination fits the
        # harmonic mixing to within 1e-4 and mix becomes undetectable.
        monitors=[td.FieldMonitor(
            center=(0, 0, 0.6), size=(0.2, 0.2, 0),
            freqs=[0.80 * freq0, 0.85 * freq0, 0.90 * freq0], name="probe")],
        run_time=2.0e-13,
        boundary_spec=td.BoundarySpec.all_sides(td.Periodic()),
        subpixel=subpixel,
    )


def test_two_constant_media_plus_dispersive_builds():
    """The vacuum and Medium(12) columns are exactly collinear (one is a multiple of the other):
    this structural degeneracy has to be recognized and merged, not rejected as "the materials are
    linearly dependent". The criterion has no parameters: compare the ε arrays cell by cell."""
    sc = scene_mod.from_simulation(_hbn_like_sim(subpixel=False))
    for arr in (sc.eps_ex, sc.eps_ey, sc.eps_ez):
        assert set(np.round(np.unique(arr), 9)) == {1.0, 3.0, 12.0}
    # poles can only hang on the Lorentz cells (eps_inf=3.0), and all three components have them
    d = sc.dispersion
    assert d is not None and d.n_entry > 0
    assert set(np.unique(d.comp)) == {0, 1, 2}
    for c in range(3):
        arr = (sc.eps_ex, sc.eps_ey, sc.eps_ez)[c]
        on = np.zeros(arr.size, dtype=bool)
        on[d.cell[d.comp == c]] = True
        # ε is reconstructed by least squares and carries noise at the 1e-10 level, so round to 9
        # digits, the same convention as the unique() above
        np.testing.assert_array_equal(on.reshape(arr.shape),
                                      np.round(arr, 9) == 3.0)


def test_serialize_roundtrip_dipole_fieldmonitor_mix(tmp_path):
    """v4 serialization: dipole, field monitor and dispersion_mix all have to come back unchanged,
    bitwise."""
    from openem import serialize

    # slab faces land between grid lines -> the normal component cells are a harmonic mean ->
    # there are necessarily mix entries
    sc = scene_mod.from_simulation(_hbn_like_sim(subpixel=True))
    assert sc.dipoles and sc.field_monitors
    assert sc.dispersion_mix is not None and sc.dispersion_mix.n_entry > 0, \
        "the scene has no mix entries, the round-trip test is vacuous"
    back = serialize.load(serialize.save(sc, tmp_path / "s.npz"))
    np.testing.assert_array_equal(back.eps_ex, sc.eps_ex)
    dp, bp = sc.dipoles[0], back.dipoles[0]
    assert bp.component == dp.component
    np.testing.assert_array_equal(bp.indices, dp.indices)
    np.testing.assert_array_equal(bp.coef, dp.coef)
    np.testing.assert_array_equal(bp.waveform.amp_int, dp.waveform.amp_int)
    np.testing.assert_array_equal(bp.waveform.amp_half, dp.waveform.amp_half)
    fm, bm = sc.field_monitors[0], back.field_monitors[0]
    assert (bm.name, bm.origin, bm.box) == (fm.name, fm.origin, fm.box)
    np.testing.assert_array_equal(bm.freqs, fm.freqs)
    mm, bmix = sc.dispersion_mix, back.dispersion_mix
    for f in ("comp", "cell", "p_ofs", "pa", "pb", "q_ofs", "qa", "qb",
              "beta", "eps_inf", "zeta_inf"):
        np.testing.assert_array_equal(getattr(bmix, f), getattr(mm, f))
    dd, bd = sc.dispersion, back.dispersion
    for f in ("comp", "cell", "pole_ofs", "am1", "b", "g"):
        np.testing.assert_array_equal(getattr(bd, f), getattr(dd, f))
