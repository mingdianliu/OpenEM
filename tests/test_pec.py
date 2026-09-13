# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""PECMedium, a perfect electric conductor: a staircased implementation where E is identically 0 on
the masked cells.

Two criteria with no free parameters:

1. A PEC slab in vacuum reflects everything: ``|R| = 1`` to within 1e-3, and the field behind the
   slab is 0 (below 1e-6 relative, though it should be exactly 0, since once the slab spans the
   cross section no difference path can get around it).
2. A scene with no PEC is bitwise unchanged: an empty mask modifies no coefficient, and two runs of
   the same scene are bitwise identical.

The mask convention (tidy3d's eps sampling returns ``pec_val = -1e8`` for PEC) is locked by a
separate test: the mask must be exactly "the Yee sample point of that component lies inside the PEC
geometry".
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import scene as scene_mod
from openem.scene import media

F0 = 2.0e14
DL = 0.05
#: z extent of the slab; 0.35 lands exactly on a primal boundary of the 0.05 grid
SLAB_LO, SLAB_HI = 0.35, 0.65


def _slab_sim(pec: bool) -> "td.Simulation":
    """Vacuum plus a slab spanning the cross section; with ``pec=False`` the slab becomes vacuum, as
    the reference scene.
    """
    med = td.PECMedium() if pec else td.Medium(permittivity=1.0)
    zc, zs = 0.5 * (SLAB_LO + SLAB_HI), SLAB_HI - SLAB_LO
    freqs = [0.8 * F0, F0, 1.2 * F0]
    return td.Simulation(
        size=(8 * DL, 8 * DL, 4.0), run_time=1.5e-13,
        grid_spec=td.GridSpec.uniform(dl=DL),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, zc), size=(td.inf, td.inf, zs)),
            medium=med)],
        sources=[td.PlaneWave(
            center=(0, 0, -0.9), size=(td.inf, td.inf, 0), direction="+",
            pol_angle=0.0,
            source_time=td.GaussianPulse(freq0=F0, fwidth=0.25 * F0))],
        monitors=[
            td.FluxMonitor(center=(0, 0, -1.3), size=(td.inf, td.inf, 0),
                           freqs=freqs, name="back"),
            td.FluxMonitor(center=(0, 0, -0.4), size=(td.inf, td.inf, 0),
                           freqs=freqs, name="net"),
            td.FluxMonitor(center=(0, 0, 1.2), size=(td.inf, td.inf, 0),
                           freqs=freqs, name="trans"),
        ],
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.periodic(), y=td.Boundary.periodic(),
            z=td.Boundary.pml()),
        subpixel=False)


def test_pec_scene_builds_with_mask_and_sane_eps():
    """A PEC structure is accepted, the eps arrays are cleaned (no pec_val remains), and the mask is
    non-empty on all three components.
    """
    sc = scene_mod.from_simulation(_slab_sim(pec=True))
    assert sc.any_pec
    for eps in (sc.eps_ex, sc.eps_ey, sc.eps_ez):
        np.testing.assert_array_equal(eps, 1.0)     # vacuum plus the cleaned slab
    for idx in (sc.pec_ex, sc.pec_ey, sc.pec_ez):
        assert idx.size > 0
    # The slab is 6 cells thick and spans the cross section: Ez, sampled at cell centers, has one
    # layer fewer than Ex and Ey, which are sampled on edges
    nx, ny, _ = sc.shape
    assert sc.pec_ex.size % (nx * ny) == 0


def test_pec_slab_reflects_everything_and_blocks_everything():
    """Criterion 1: |R| = 1 to within 1e-3, and everything behind the slab, both the monitor and the
    final field, is 0.
    """
    pytest.importorskip("cupy")
    from openem import flux as flux_mod
    from openem import solver
    from openem.device import Kernels

    sc = scene_mod.from_simulation(_slab_sim(pec=True))
    res = solver.run(sc, use_shutoff=False, kernels=Kernels(), verbose=False,
                     return_fields=True)

    back = flux_mod.plane_flux(res.phasors["back"].data, sc.grid)   # reflection, travelling -z
    net = flux_mod.plane_flux(res.phasors["net"].data, sc.grid)     # incident minus reflected
    trans = flux_mod.plane_flux(res.phasors["trans"].data, sc.grid)
    p_in = net + np.abs(back)              # back travels -z, so take its absolute value
    r_amp = np.sqrt(np.abs(back) / p_in)
    assert np.all(np.abs(r_amp - 1.0) < 1e-3), f"|R| = {r_amp}"
    # Flux behind the slab: the field is identically 0, so the DFT is too, exactly. Far stronger than
    # the 1e-6 criterion.
    assert np.all(np.abs(trans) <= 1e-12 * p_in), f"transmission {trans / p_in}"

    # Final field: E is exactly 0 on the PEC cells, and every component is exactly 0 in the region
    # behind the slab, with a two-cell margin
    for comp, idx in (("Ex", sc.pec_ex), ("Ey", sc.pec_ey), ("Ez", sc.pec_ez)):
        vals = res.fields[comp].ravel()[idx]
        assert np.all(vals == 0.0), f"{comp} is non-zero on a PEC cell: max {np.abs(vals).max()}"
    k_hi = int(np.searchsorted(sc.grid.z.edges, SLAB_HI * 1e-6 * 1.0000001)) + 2
    peak = max(float(np.abs(res.fields[c]).max()) for c in res.fields)
    behind = max(float(np.abs(res.fields[c][:, :, k_hi:]).max()) for c in res.fields)
    assert behind <= 1e-6 * peak, f"residual behind the slab {behind:.3e} against a peak of {peak:.3e}"


def test_no_pec_scene_is_bitwise_unchanged():
    """Criterion 2: a scene with no PEC (a None mask) is bitwise identical to one with empty mask
    arrays, i.e. a single code path.
    """
    pytest.importorskip("cupy")
    from openem import solver
    from openem.device import Kernels

    sc = scene_mod.from_simulation(_slab_sim(pec=False))
    assert not sc.any_pec and sc.pec_ex is None
    empty = np.zeros(0, dtype=np.int64)
    sc2 = replace(sc, pec_ex=empty, pec_ey=empty.copy(), pec_ez=empty.copy())
    k = Kernels()
    r1 = solver.run(sc, use_shutoff=False, kernels=k, verbose=False)
    r2 = solver.run(sc2, use_shutoff=False, kernels=k, verbose=False)
    for name in r1.phasors:
        np.testing.assert_array_equal(r1.phasors[name].data, r2.phasors[name].data)


def test_pec_mask_is_yee_point_sampling():
    """The mask equals "the Yee sample point of that component lies inside the PEC geometry",
    exactly, component by component.
    """
    r = 0.2371                        # deliberately a radius that never ties with a sample point
    sim = td.Simulation(
        size=(1, 1, 1), run_time=1e-13, grid_spec=td.GridSpec.uniform(dl=DL),
        structures=[td.Structure(geometry=td.Sphere(center=(0, 0, 0), radius=r),
                                 medium=td.PECMedium())],
        sources=[td.PointDipole(center=(0, 0, 0.4), polarization="Ex",
                                source_time=td.GaussianPulse(freq0=F0,
                                                             fwidth=0.25 * F0))],
        boundary_spec=td.BoundarySpec.all_sides(boundary=td.Periodic()),
        subpixel=False)
    sc = scene_mod.from_simulation(sim)
    for comp, idx in (("Ex", sc.pec_ex), ("Ey", sc.pec_ey), ("Ez", sc.pec_ez)):
        g = sim.grid.yee.grid_dict[comp]
        x, y, z = (np.asarray(getattr(g, ax), dtype=np.float64)[: n]
                   for ax, n in zip("xyz", sc.shape))
        inside = (x[:, None, None] ** 2 + y[None, :, None] ** 2
                  + z[None, None, :] ** 2) < r ** 2
        np.testing.assert_array_equal(idx, np.flatnonzero(inside.ravel()),
                                      err_msg=f"the {comp} mask disagrees with the point-sampling convention")


def test_pec_with_dispersive_medium_fails_closed():
    sim = _slab_sim(pec=True)
    lorentz = td.Lorentz(eps_inf=2.0, coeffs=[(1.0, 3e14, 1e13)])
    sim = sim.updated_copy(structures=list(sim.structures) + [
        td.Structure(geometry=td.Box(center=(0, 0, -1.7), size=(td.inf, td.inf, 0.1)),
                     medium=lorentz)])
    with pytest.raises(NotImplementedError, match="dispersive"):
        scene_mod.from_simulation(sim)


def test_pec_background_fails_closed():
    sim = _slab_sim(pec=False).updated_copy(medium=td.PECMedium())
    with pytest.raises(NotImplementedError, match="background"):
        media.refuse_unsupported(sim)


def test_serialize_roundtrip_keeps_pec_and_projection():
    from openem import serialize
    from openem.model import ProjectionMonitor
    import tempfile, pathlib

    sc = scene_mod.from_simulation(_slab_sim(pec=True))
    sc = replace(sc, projection_monitors=[ProjectionMonitor(
        name="ff", surface_names=("ff__n2f_x-", "ff__n2f_x+"),
        normal_dirs=("-", "+"))])
    with tempfile.TemporaryDirectory() as d:
        p = pathlib.Path(d) / "scene.npz"
        serialize.save(sc, p)
        back = serialize.load(p)
    for f in ("pec_ex", "pec_ey", "pec_ez"):
        np.testing.assert_array_equal(getattr(back, f), getattr(sc, f))
    assert back.projection_monitors == sc.projection_monitors
