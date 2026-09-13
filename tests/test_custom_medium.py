# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""The ``CustomMedium`` family: spatially varying ε and dispersion coefficients must not be lost
in serialization.

Background (2026-08-30): ``openem_run`` used to write the sim with ``to_file(simulation.json)``
and read it back before building the Scene. tidy3d's **JSON serialization silently drops** the
spatial arrays of CustomMedium / TriangleMesh, so the coefficients of ``CustomLorentz`` were gone,
no poles could be fitted, and ``dispersion.build`` raised "the scene has no non-DC poles". The
whole CustomMediumTutorial notebook was stuck on this, while the solver itself had supported it
all along.

This test locks two things:

1. all three Custom media (spatially varying coefficients included) build a Scene;
2. **a JSON round trip really does lose data**. If tidy3d ever fixes JSON this test turns red,
   and only then may the hdf5 in ``nb.backend`` go back to json.
"""

from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import scene as scene_mod

F0 = 2.0e14
N = 8
_ax = np.linspace(-0.2, 0.2, N)
_one = np.ones((N, N, N))
_vary = 1.0 + 0.3 * np.linspace(0, 1, N * N * N).reshape(N, N, N)


def _sda(v):
    return td.SpatialDataArray(v, coords=dict(x=_ax, y=_ax, z=_ax))


def _lorentz(de):
    return td.CustomLorentz(
        eps_inf=_sda(_one * 2.0),
        coeffs=[(_sda(de), _sda(_one * 3.0e14), _sda(_one * 1.0e13))])


def _sim(med):
    return td.Simulation(
        size=(0.8, 0.8, 0.8), run_time=2e-14,
        grid_spec=td.GridSpec.uniform(dl=0.05),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, 0), size=(0.4, 0.4, 0.4)), medium=med)],
        sources=[td.PointDipole(
            center=(0, 0, -0.3), polarization="Ex",
            source_time=td.GaussianPulse(freq0=F0, fwidth=0.3 * F0))],
        monitors=[td.FluxMonitor(center=(0, 0, 0.3), size=(td.inf, td.inf, 0),
                                 freqs=[F0], name="flux")],
        boundary_spec=td.BoundarySpec.all_sides(td.PML()), subpixel=False)


CASES = {
    "CustomMedium": lambda: td.CustomMedium(permittivity=_sda(2.0 + _vary)),
    "CustomLorentz per-cell coefficients": lambda: _lorentz(_vary),
    # the shape of the real CustomMediumTutorial case: only zz is dispersive
    "CustomAnisotropic zz dispersion": lambda: td.CustomAnisotropicMedium(
        xx=td.CustomMedium(permittivity=_sda(_one * 2.0)),
        yy=td.CustomMedium(permittivity=_sda(_one * 2.0)),
        zz=_lorentz(_vary)),
}


#: The two per-cell-coefficient cases, supported since 2026-08-31 by reading the medium's own
#: coefficients directly.
PERCELL = {"CustomLorentz per-cell coefficients", "CustomAnisotropic zz dispersion"}


@pytest.mark.parametrize("name", list(CASES))
def test_custom_media_build_from_object(name):
    """Build the tables straight from the in-memory sim object, which is the path openem_run
    takes now."""
    sc = scene_mod.from_simulation(_sim(CASES[name]()))
    assert sc.grid.num_cells > 0


@pytest.mark.parametrize("name", sorted(PERCELL))
def test_percell_weights_are_geometric_fills(name):
    """With per-cell coefficients the weights have to go back to their literal meaning, the
    geometric fill fraction.

    This is the whole argument: without it, a two-column fit has a residual at **machine
    precision** (2.6e-15) while asking for a weight of −0.25. A residual criterion cannot catch
    that at all, and a negative weight is a gain medium. The sim here has ``subpixel=False`` and
    the box is aligned with the grid, so every cell is either pure or empty and the weight can
    only be 0 or 1.
    """
    from openem.scene import weights as W

    sim = _sim(CASES[name]())
    ck = "Ez" if "zz" in name else "Ex"
    meds = W.media_columns(sim, "xyz".index(ck[1]))
    freqs = W.fit_frequencies(sim, len(meds))
    w, mix, harm, _ = W.decompose_full(sim, ck, freqs, meds)

    assert mix is None and not harm.any(), "an aligned grid with subpixel=False has no interface cells"
    assert float(w.min()) >= 0.0, f"negative geometric fill fraction {w.min():.3e}"
    near = np.minimum(np.abs(w), np.abs(w - 1.0))
    assert float(near.max()) < 1e-6, (
        f"the fill fraction of a pure cell should be 0 or 1, largest deviation {near.max():.3e}")


def test_percell_pole_scale_matches_the_coefficient_array():
    """The per-cell residue strength s_c has to equal de(cell)/de(reference cell) from the
    coefficient array.

    This is an **independent** criterion: s_c is measured off ε(f) while de is what the medium
    stores itself, so the two agreeing is what shows the direct-read path really reads the right
    thing.
    """
    from openem.scene import weights as W

    sim = _sim(CASES["CustomLorentz per-cell coefficients"]())
    meds = W.media_columns(sim, 0)
    freqs = W.fit_frequencies(sim, len(meds))
    w, _mix, _harm, ps = W.decompose_full(sim, "Ex", freqs, meds)
    assert ps is not None, "the per-cell path has to produce pole_scale"

    pc = W.percell_indices(meds)
    assert len(pc) == 1
    inside = w.reshape(len(meds), -1)[pc[0]] > 0.5
    got = ps[pc[0]][inside]

    # independently interpolate the de array onto the same Yee coordinates (through the
    # **coefficient array**, not through ε)
    med = meds[pc[0]]
    de = med.coeffs[0][0]
    co = sim.grid["Ex"]
    de_grid = np.asarray(
        de.interp(x=np.asarray(co.x), y=np.asarray(co.y), z=np.asarray(co.z),
                  method="nearest").values, dtype=np.float64).ravel()[inside]

    # the criterion is **independent of the reference cell**: normalize each side by its own
    # minimum, then compare point by point
    a = got / got.min()
    b = de_grid / de_grid.min()
    assert a.shape == b.shape and a.size > 1, "too few cells picked up, the criterion does not hold"
    assert float(np.abs(a - b).max()) < 1e-6, (
        f"s_c does not match de in the coefficient array, largest pointwise diff {np.abs(a - b).max():.3e}"
        f" (s_c ranges {a.min():.6f}~{a.max():.6f}, de ranges {b.min():.6f}~{b.max():.6f})")


def test_json_roundtrip_still_loses_custom_data(tmp_path):
    """A JSON round trip still loses data. If this turns red tidy3d has fixed it, and only then
    can we go back to json."""
    sim = _sim(CASES["CustomLorentz per-cell coefficients"]())
    f = tmp_path / "sim.json"
    sim.to_file(str(f))
    with pytest.raises(Exception):
        scene_mod.from_simulation(td.Simulation.from_file(str(f)))


def test_hdf5_roundtrip_preserves_custom_data(tmp_path):
    """An hdf5 round trip loses nothing, which is why openem_run archives with it."""
    sim = _sim(CASES["CustomLorentz per-cell coefficients"]())
    f = tmp_path / "sim.hdf5"
    sim.to_file(str(f))
    # hdf5 really does **not** lose data: the tables built after a round trip are **identical** to
    # the ones built from the in-memory object. (Before 2026-08-31 this could only be asserted
    # indirectly, as "the error changed to a different one", because per-cell coefficients could
    # not build a table at all back then. Now they can, so the criterion compares directly.)
    a = scene_mod.from_simulation(td.Simulation.from_file(str(f)))
    b = scene_mod.from_simulation(sim)
    np.testing.assert_allclose(a.eps_ex, b.eps_ex, rtol=1e-12)
    np.testing.assert_array_equal(a.dispersion.cell, b.dispersion.cell)
    np.testing.assert_allclose(a.dispersion.b, b.dispersion.b, rtol=1e-12)


F0_DISP = 2.0e14
_N = 8
_XY = np.linspace(-0.4, 0.4, _N)
_Z = np.linspace(-0.25, 0.25, _N)
_C = dict(x=_XY, y=_XY, z=_Z)


def _sda3(v):
    return td.ScalarFieldDataArray(v, coords=_C)


def _polresidue(eps_inf_arr, residue_arr):
    """A CustomPoleResidue with two conjugate poles whose per-cell residue is residue_arr."""
    poles = [(_sda3(np.full((_N, _N, _N), a)), _sda3(residue_arr))
             for a in (-1e14 - 3e14j, -1e14 + 3e14j)]
    return td.CustomPoleResidue(
        eps_inf=_sda3(eps_inf_arr), poles=poles, interp_method="linear")


def _disp_sim(med, source):
    """The design region fills the whole domain in x/y (so the source plane is uniform) and is
    bounded only in z."""
    return td.Simulation(
        size=(0.8, 0.8, 0.9), run_time=2e-14,
        grid_spec=td.GridSpec.uniform(dl=0.05),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, 0), size=(td.inf, td.inf, 0.5)),
            medium=med)],
        sources=[source],
        monitors=[td.FluxMonitor(center=(0, 0, 0.4), size=(td.inf, td.inf, 0),
                                 freqs=[F0_DISP], name="flux")],
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.periodic(), y=td.Boundary.periodic(),
            z=td.Boundary.pml()), subpixel=False)


def test_zero_strength_dispersion_leaves_no_dispersive_cells():
    """A Custom dispersive medium whose per-cell residue is 0 everywhere **must not** produce
    dispersive cells.

    When the per-cell strength s_c = 0 the residue b = base·w·s_c ≡ 0, that cell contributes
    nothing to ε(ω) and is equivalent to a non-dispersive background. Recording it as a dispersive
    cell anyway makes ``sources.plane_wave`` reject a plane-wave source plane that lands on it, on
    the mistaken grounds that it sits inside a dispersive medium; that was the real failure of the
    normalization run of Autograd15 with an all-zero density. Here the source plane (z=0) lands
    inside the design region.
    """
    med = _polresidue(np.ones((_N, _N, _N)), np.zeros((_N, _N, _N)))  # eps_inf=1, residue 0
    src = td.PlaneWave(
        center=(0, 0, 0.0), size=(td.inf, td.inf, 0), direction="+",
        source_time=td.GaussianPulse(freq0=F0_DISP, fwidth=0.3 * F0_DISP),
        pol_angle=0.0)
    sc = scene_mod.from_simulation(_disp_sim(med, src))    # must not raise "in a dispersive medium"
    ncell = 0 if sc.dispersion is None else int(sc.dispersion.cell.size)
    assert ncell == 0, f"a zero-strength medium should have no dispersive cells, measured {ncell}"
    assert float(sc.eps_ez.min()) == pytest.approx(1.0, abs=1e-9)


def test_partial_strength_keeps_only_nonzero_cells():
    """Half the residues zero, half non-zero: only the non-zero half of z becomes dispersive."""
    de = np.zeros((_N, _N, _N))
    de[:, :, _N // 2:] = 1.0                               # only the upper half of z has a residue
    med = _polresidue(np.ones((_N, _N, _N)), 3e14 * de)    # eps_inf is constant 1 (no coupling)
    src = td.PointDipole(
        center=(0, 0, -0.4), polarization="Ex",
        source_time=td.GaussianPulse(freq0=F0_DISP, fwidth=0.3 * F0_DISP))
    sc = scene_mod.from_simulation(_disp_sim(med, src))
    assert sc.dispersion is not None
    ncell = int(sc.dispersion.cell.size)
    total = int(np.prod(sc.shape))
    assert 0 < ncell < total, f"only some cells should be dispersive, measured {ncell}/{total}"
