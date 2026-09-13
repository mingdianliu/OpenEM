# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Two-dimensional sheet materials (td.Medium2D, the MoS2Waveguide family).

How they are supported: **reuse the tidy3d client's own conversion**,
``sim.volumetric_structures`` (Medium2D -> ``AnisotropicMediumFromMedium2D``: the in-plane
response is scaled by 1/dl and weight-averaged with the neighboring medium, and the geometry is
snapped to the grid), then take the existing diagonal-anisotropic + PoleResidue path.
``sim.epsilon`` uses that same conversion internally, so the ε queries agree automatically; all we
have to do is look at the converted structures in the three places that enumerate media by type
(the allowlist, the weight columns, any_dispersive).

**The criterion is bitwise identity with zero tunable parameters**: compare scene A, with the
Medium2D, against scene B, where the equivalent volumetric structures from
``sim.volumetric_structures`` are placed by hand. B contains no Medium2D and takes the path
test_anisotropic locks. The ε arrays, the dispersion tables and the run results of the two must be
bitwise identical. What this locks is the wiring (all three enumerations use the same conversion);
the physics is backed by the existing verification of the anisotropic path plus tidy3d's own
conversion.

The grid is uniform: AutoGrid produces different grids for the original and the converted
structures (measured), which is why the implementation **always takes the grid from the original
sim**. A uniform grid makes A and B share a grid automatically, which is what makes the comparison
valid at all.
"""

from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import scene as scene_mod

FREQ0 = 2.0e14
DL = 0.05


def _lorentz_sheet(thickness: float = 0.01) -> "td.Medium2D":
    """A 2D material sheet with a single Lorentz pole (a minimal MoS2: the exciton resonance)."""
    med3d = td.Lorentz(eps_inf=1.0,
                       coeffs=((8.0, 1.05 * FREQ0, 0.05 * FREQ0),))
    return td.Medium2D.from_dispersive_medium(med3d, thickness=thickness)


def _sheet_sim(structures: list, run_time: float = 2.0e-13,
               z_pml: bool = False) -> "td.Simulation":
    """A sheet in vacuum with normal z (z=0.3 lands exactly on a dl=0.05 grid line), with an Ex
    plane wave at normal incidence.

    The bitwise criteria use Periodic everywhere, which removes PML as a variable; the
    transmission criterion needs ``z_pml=True``, because with periodic boundaries the pulse wraps
    around and re-enters the monitor, which makes a single-frequency flux ratio meaningless.
    """
    zb = td.Boundary(plus=td.PML(), minus=td.PML()) if z_pml else td.Boundary(
        plus=td.Periodic(), minus=td.Periodic())
    return td.Simulation(
        size=(0.4, 0.4, 2.4),
        grid_spec=td.GridSpec.uniform(dl=DL),
        medium=td.Medium(permittivity=1.0),
        structures=structures,
        sources=[td.PlaneWave(
            center=(0, 0, -0.8), size=(td.inf, td.inf, 0.0), direction="+",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.3 * FREQ0))],
        monitors=[
            td.FieldMonitor(center=(0.05, 0, 0.6), size=(0.2, 0.2, 0.4),
                            freqs=[FREQ0], name="probe"),
            td.FluxMonitor(center=(0, 0, 0.7), size=(td.inf, td.inf, 0),
                           freqs=[FREQ0], name="T"),
        ],
        run_time=run_time,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary(plus=td.Periodic(), minus=td.Periodic()),
            y=td.Boundary(plus=td.Periodic(), minus=td.Periodic()),
            z=zb),
        subpixel=False,
    )


def _sheet_structure() -> "td.Structure":
    return td.Structure(
        geometry=td.Box(center=(0, 0, 0.3), size=(td.inf, td.inf, 0)),
        medium=_lorentz_sheet())


def _converted_twin(sim: "td.Simulation") -> "td.Simulation":
    """Scene B: the equivalent volumetric structures from tidy3d's conversion, placed by hand,
    with no Medium2D.

    On a uniform grid the grid is the same before and after the conversion (the snapping is the
    identity), and that is asserted here: if the grid changed, a bitwise comparison of A and B
    would be meaningless.
    """
    twin = sim.updated_copy(structures=list(sim.volumetric_structures))
    for ax in "xyz":
        np.testing.assert_array_equal(
            np.asarray(getattr(sim.grid.boundaries, ax)),
            np.asarray(getattr(twin.grid.boundaries, ax)))
    assert not any(isinstance(s.medium, td.Medium2D) for s in twin.structures)
    return twin


def test_sheet_build_bit_identical_to_converted():
    """The Scenes built from A (Medium2D) and B (hand-placed equivalent volume) have to be
    bitwise identical.

    Covers: the allowlist letting it through, the weight columns, the (medium, axis) pole table,
    and the normal component (ε_zz) taking the background.
    """
    sim = _sheet_sim([_sheet_structure()])
    sc_a = scene_mod.from_simulation(sim)
    sc_b = scene_mod.from_simulation(_converted_twin(sim))
    for f in ("eps_ex", "eps_ey", "eps_ez"):
        np.testing.assert_array_equal(getattr(sc_a, f), getattr(sc_b, f))
    d_a, d_b = sc_a.dispersion, sc_b.dispersion
    assert d_a is not None and d_a.n_entry > 0
    for f in ("comp", "cell", "pole_ofs", "am1", "b", "g"):
        np.testing.assert_array_equal(getattr(d_a, f), getattr(d_b, f))
    # the sheet's normal component (Ez) has to be the background: the zz of the equivalent volume
    # takes the medium on the + side, which is vacuum here
    assert set(np.unique(sc_a.eps_ez)) == {1.0}
    # the in-plane components (Ex/Ey) are not vacuum on the cells holding the sheet
    assert set(np.unique(sc_a.eps_ex)) != {1.0}
    # dispersion hangs only on the in-plane components
    assert set(np.unique(d_a.comp)) <= {0, 1}


def test_sheet_run_bit_identical_to_converted():
    """Run scenes A and B for 2000 steps: the probe phasors are bitwise identical."""
    pytest.importorskip("cupy")
    from openem import solver
    from openem.device import Kernels

    sim = _sheet_sim([_sheet_structure()])
    sc_a = scene_mod.from_simulation(sim)
    sc_b = scene_mod.from_simulation(_converted_twin(sim))
    assert sc_a.dt == sc_b.dt, "the two runs have different dt, a bitwise comparison is meaningless"
    r_a = solver.run(sc_a, num_steps=2000, use_shutoff=False,
                     kernels=Kernels(), verbose=False)
    r_b = solver.run(sc_b, num_steps=2000, use_shutoff=False,
                     kernels=Kernels(), verbose=False)
    p_a, p_b = r_a.field_phasors["probe"], r_b.field_phasors["probe"]
    assert float(np.max(np.abs(p_b))) > 0, "the reference run is all zeros, the criterion is vacuous"
    assert np.all(np.isfinite(p_a)), "the Medium2D run diverged"
    np.testing.assert_array_equal(p_a, p_b)


def test_sheet_transmission_matches_thin_film():
    """Physical check: the transmittance matches the **Fresnel formula generalized to a
    zero-thickness conducting sheet**, with no fitted parameters.

    At normal incidence with vacuum on both sides, a sheet with surface conductivity
    σ_s(ω) = −iωε₀χ_2D(ω) gives ``t = 1/(1 + η₀σ_s/2)`` and ``T = |t|²``. χ_2D is taken straight
    from the conversion result (χ_2D = (ε_eq − ε_bg)·dl), the formula is not re-derived.

    Measured convergence to that analytic value is second order in dl: at dl = 0.05/0.025/0.0125
    the absolute deviations are 5.28e-3 / 1.34e-3 / 3.47e-4 (ratio ≈ 3.9, i.e. O(dl²)). The test
    runs dl=0.05 (2.4% relative deviation) with a tolerance of 5e-2: the convergence order already
    settles that the deviation is discretization error rather than a wrong model, and a wrongly
    wired conversion or ADE would be off by orders of magnitude.
    """
    pytest.importorskip("cupy")
    from openem import solver
    from openem.device import Kernels

    from openem import flux as flux_mod

    def _t_flux(structures):
        sim = _sheet_sim(structures, run_time=6.0e-13, z_pml=True)
        sc = scene_mod.from_simulation(sim)
        res = solver.run(sc, num_steps=int(sim.num_time_steps),
                         use_shutoff=False, kernels=Kernels(), verbose=False)
        return sim, float(flux_mod.plane_flux(res.phasors["T"].data, sc.grid)[0])

    sim, f_sheet = _t_flux([_sheet_structure()])
    _, f_ref = _t_flux([])   # vacuum reference: the normalization cancels exactly in the ratio
    t_num = f_sheet / f_ref

    # χ_2D straight from the conversion: the ε_eq(ω) of the equivalent cell smeared back into a
    # surface response
    equiv = list(sim.volumetric_structures)[0].medium
    eps_w = complex(equiv.xx.eps_model(FREQ0))
    chi_2d = (eps_w - 1.0) * DL * 1e-6            # SI [m]
    eps0, eta0 = 8.8541878128e-12, 376.730313668
    sigma_s = -1j * (2 * np.pi * FREQ0) * eps0 * chi_2d
    t_an = abs(1.0 / (1.0 + eta0 * sigma_s / 2.0)) ** 2
    assert abs(t_num - t_an) / t_an < 5e-2, (t_num, t_an)


def test_2d_sheet_waveguide_mini():
    """A miniature MoS2Waveguide: **a 2D simulation (nz=1) with a y-normal sheet and a
    ModeSource**.

    It locks three things that are specific to the real case and that no other test touches:
    * the sheet's normal is y (the equivalent volume has xx/zz in-plane and yy as background),
      not the z of the tests above;
    * the collapsed axis at nz=1: the z coordinate of the mode profile is the plane position (0),
      not a Yee center;
    * the ModeSolver's plane grid extends one cell past the domain, and after the out-of-range
      points are trimmed the match is exact point by point.

    The sheet uses a Lorentz model far from resonance (the resonance is at 2×FREQ0), so at FREQ0
    it is nearly real with ε_eq>1, giving a bound sheet-waveguide mode. The criteria: the scene
    builds, n_eff>1 (it really guides), and a short run stays finite.
    """
    pytest.importorskip("cupy")
    from openem import solver
    from openem.device import Kernels

    sheet = td.Medium2D.from_dispersive_medium(
        td.Lorentz(eps_inf=1.0, coeffs=((30.0, 2.0 * FREQ0, 0.01 * FREQ0),)),
        thickness=0.02)
    sim = td.Simulation(
        size=(6.0, 4.0, 0.0),
        grid_spec=td.GridSpec.uniform(dl=DL),
        medium=td.Medium(permittivity=1.0),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, 0), size=(td.inf, 0, td.inf)),
            medium=sheet)],
        sources=[td.ModeSource(
            center=(-2.0, 0, 0), size=(0, td.inf, td.inf), direction="+",
            mode_spec=td.ModeSpec(num_modes=1), mode_index=0,
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.1 * FREQ0))],
        monitors=[td.FieldMonitor(center=(1.5, 0, 0), size=(0.5, 1.0, 0),
                                  freqs=[FREQ0], name="probe")],
        run_time=4.0e-13,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary(plus=td.PML(), minus=td.PML()),
            y=td.Boundary(plus=td.PML(), minus=td.PML()),
            z=td.Boundary(plus=td.Periodic(), minus=td.Periodic())),
        subpixel=True,
    )
    sc = scene_mod.from_simulation(sim)
    assert sc.grid.shape[2] == 1
    assert len(sc.mode_sources) == 1
    assert sc.mode_sources[0].n_eff.real > 1.0, "the sheet mode is not bound, the test scene is void"
    assert sc.dispersion is not None and sc.dispersion.n_entry > 0
    # in-plane is x/z; the normal direction y must carry no dispersion
    assert 1 not in set(np.unique(sc.dispersion.comp))
    r = solver.run(sc, num_steps=int(sim.num_time_steps), use_shutoff=False,
                   kernels=Kernels(), verbose=False)
    p = r.field_phasors["probe"]
    assert float(np.max(np.abs(p))) > 0, "the run is all zeros, the criterion is vacuous"
    assert np.all(np.isfinite(p)), "the 2D sheet waveguide run diverged"


def test_sheet_normal_component_stays_background():
    """Replace the sheet's in-plane response with a completely different value and ε_zz is still
    all background: nothing leaks into the normal direction."""
    sim = _sheet_sim([td.Structure(
        geometry=td.Box(center=(0, 0, 0.3), size=(td.inf, td.inf, 0)),
        medium=_lorentz_sheet(thickness=0.05))])
    sc = scene_mod.from_simulation(sim)
    assert set(np.unique(sc.eps_ez)) == {1.0}
