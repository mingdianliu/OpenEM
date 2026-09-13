# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""End to end: BiosensorGrating's simulation.json -> Scene.

Every criterion comes from measured values in the validation set and from the reference solver's
log, not from expected values we worked out ourselves.
"""

import os
import pathlib

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import scene as scene_mod

TABLEA = os.environ.get("OPENEM_TABLEA", "/nonexistent")  # external reference archive; skipped if unset
SIM_JSON = f"{TABLEA}/BiosensorGrating/Tidy3D/simulation.json"


@pytest.fixture(scope="module")
def sc():
    return scene_mod.from_file(SIM_JSON)


def test_tidy3d_version_is_2_12():
    """OpenEM has to run in the same environment as the validation set (tidy3d 2.12.0).

    Every conclusion in the validation set was measured on 2.12.0; changing the version turns
    what has been verified back into an unknown. The FDTDX venv has 2.10.0 and cannot be used to
    run OpenEM.
    """
    assert td.__version__.startswith("2.12"), (
        f"tidy3d is {td.__version__}, OpenEM requires 2.12.x: "
        "use /usr/bin/python3, not the FDTDX venv"
    )


def test_grid_matches_cloud_log(sc):
    """Both cell counts have to match the numbers in the reference solver's log."""
    assert sc.shape == (22, 16, 140)          # Simulation domain Nx, Ny, Nz
    assert sc.grid.num_cells == 49_280
    assert sc.grid.num_cells_cloud == 49_984  # Number of computational grid points


def test_dt_bit_exact(sc):
    """dt has to be bitwise equal to sim.dt: it comes from Tidy3D, we do not compute it."""
    assert sc.dt == 4.766437173827506e-17
    assert sc.num_time_steps == 209_802
    assert sc.shutoff == 1e-6


def test_units_are_meters(sc):
    """Everything internal is in SI metres. BiosensorGrating spans 0.55 µm in x."""
    x = sc.grid.x.edges
    assert x[0] == pytest.approx(-0.275e-6)
    assert x[-1] == pytest.approx(0.275e-6)
    np.testing.assert_allclose(sc.grid.x.dl, 0.025e-6, rtol=1e-12)


def test_grid_is_uniform_all_axes(sc):
    """All three axes of this case are a uniform 0.025 µm, which means a non-uniform Yee grid
    cannot be exercised on it.

    The z grid of 140 is 116 physical cells plus 12+12 PML, not a non-uniform grid. This one sits
    on the list of things that remain unverified.
    """
    for ax in sc.grid.axes:
        np.testing.assert_allclose(ax.dl, 0.025e-6, rtol=1e-12)


def test_boundaries(sc):
    assert sc.grid.x.is_periodic and sc.grid.y.is_periodic
    assert not sc.grid.z.is_periodic
    # only the two z ends are PML
    assert set(sc.pml) == {(2, "lo"), (2, "hi")}
    for c in sc.pml.values():
        assert c.num_layers == 12


def test_epsilon_values(sc):
    """ε is point-sampled from sim.epsilon, so there should be exactly three media."""
    for eps in (sc.eps_ex, sc.eps_ey, sc.eps_ez):
        assert eps.shape == (22, 16, 140)
        vals = np.unique(np.round(eps, 6))
        assert set(vals) <= {1.776889, 2.25, 4.2025}, f"unexpected ε values: {vals}"
        assert 1.776889 in vals  # the background (water, n=1.333)


def test_source(sc):
    assert len(sc.sources) == 1
    s = sc.sources[0]
    assert s.axis == 2 and s.direction == +1
    assert s.angle_theta == 0.0
    assert s.num_freqs == 3          # broadband injection, chebyshev
    assert s.pol_axis == 0           # pol_angle=0 -> E along the first transverse axis
    # The source plane snaps to the grid edge at z = -1.1 µm.
    # Note this must not be phrased as "which cell does it fall into": edges[26] is slightly
    # larger than -1.1e-6 in floating point, so that criterion would land on 25, off by a whole
    # cell (see the docstring of grid.nearest_edge_index).
    assert s.plane_index == 26
    assert abs(sc.grid.z.edges[s.plane_index] - (-1.1e-6)) < 1e-18
    # the waveform is sampled over the full nominal step count (an early shutoff must not
    # truncate the normalization waveform)
    assert s.waveform.amp_int.size == 209_803
    assert s.waveform.amp_half.size == 209_803


def test_waveform_is_real_and_peaks_inside(sc):
    w = sc.sources[0].waveform
    assert np.isrealobj(w.amp_int)
    assert np.isfinite(w.amp_int).all()
    # GaussianPulse offset=5.0 -> the peak is very early, not at the end
    assert 0 < w.peak_index < w.num_steps // 10


def test_flux_monitors(sc):
    assert [m.name for m in sc.flux_monitors] == ["T", "R"]
    for m in sc.flux_monitors:
        assert m.axis == 2
        assert m.normal_dir == +1
        assert m.freqs.shape == (5,)
        assert m.source_spectrum.shape == (5,)
        assert np.all(np.abs(m.source_spectrum) > 0), "a zero source spectrum blows up the normalization"
    t, r = sc.flux_monitors
    assert (t.plane_index, r.plane_index) == (116, 20)
    assert abs(sc.grid.z.edges[t.plane_index] - 1.15e-6) < 1e-18
    assert abs(sc.grid.z.edges[r.plane_index] - (-1.25e-6)) < 1e-18


def test_reference_flux_sums_to_one():
    """The reference data of the validation set carries its own check: the structure is lossless,
    so T + |R| = 1.

    This does not depend on our code; it confirms the convention of the reference data. The
    incident power is normalized to 1, so the absolute injection amplitude of the plane wave is
    not an unknown, which is what unblocked that question.
    """
    import h5py

    if not os.path.exists(f"{TABLEA}/BiosensorGrating/Tidy3D/data.hdf5"):
        pytest.skip("external archive missing (OPENEM_TABLEA must point at an archive with data.hdf5)")
    with h5py.File(f"{TABLEA}/BiosensorGrating/Tidy3D/data.hdf5", "r") as f:
        t = np.asarray(f["data/0/flux/__xarray_dataarray_variable__"])
        r = np.asarray(f["data/1/flux/__xarray_dataarray_variable__"])
    total = t + np.abs(r)
    np.testing.assert_allclose(total, 1.0, atol=2e-3)


def test_subpixel_supported_with_extras(sc):
    """subpixel=True is **supported** now: with tidy3d-extras installed the averaged tensor is
    available.

    Since 2026-08-12 the earlier assumption has been overturned: it is not that the client cannot
    produce the averaged ε, it is that local_subpixel from tidy3d-extras is needed.
    """
    pytest.importorskip("tidy3d_extras", reason="tidy3d-extras is needed to turn subpixel on")
    sim = td.Simulation.from_file(SIM_JSON).updated_copy(subpixel=td.SubpixelSpec())
    sub = scene_mod.from_simulation(sim)
    # anisotropy is direct evidence that subpixel took effect: eps_xx != eps_zz on the same cell
    aniso = float(np.max(np.abs(sub.eps_ex - sub.eps_ez)))
    assert aniso > 0, "subpixel is on but there is no anisotropy; it likely fell back to point sampling"
    # point sampling gives only 3 pure values; averaging produces intermediate ones
    assert len(np.unique(np.round(sub.eps_ex, 6))) > 3
    # Control: the same case with subpixel off has to give pure values and zero anisotropy.
    # This also locks "having extras installed does not contaminate a subpixel=False scene",
    # which is what tidy3d's use_local_subpixel=None default would do.
    off = scene_mod.from_file(SIM_JSON)
    assert float(np.max(np.abs(off.eps_ex - off.eps_ez))) == 0.0
    assert len(np.unique(np.round(off.eps_ex, 6))) == 3


def test_subpixel_fails_closed_without_extras():
    """Without tidy3d-extras, subpixel=True has to raise and **must not** fall back silently to
    point sampling.

    This one matters: tidy3d defaults to ``use_local_subpixel=None`` and does not even print a
    warning when extras is unavailable (``quiet=(preference is None)`` in ``packaging.py:327``).
    A silent fallback means carrying the wrong ε without noticing.
    """
    import subprocess
    import sys

    if not os.path.exists(SIM_JSON):
        pytest.skip(f"archive missing: {SIM_JSON}")
    code = (
        "import tidy3d as td\n"
        "from openem import scene\n"
        f"sim = td.Simulation.from_file({SIM_JSON!r}).updated_copy(subpixel=td.SubpixelSpec())\n"
        "try:\n"
        "    scene.from_simulation(sim)\n"
        "except NotImplementedError as e:\n"
        "    assert 'tidy3d-extras' in str(e), str(e)\n"
        "    print('OK')\n"
        "else:\n"
        "    raise SystemExit('no extras, yet no error: it may have fallen back to point sampling')\n"
    )
    # PYTHONPATH keeps only the repo root plus nb_pylibs (tidy3d lives there; the job environment
    # has no system-wide tidy3d) and leaves extras out
    root = str(pathlib.Path(__file__).resolve().parent.parent)
    env = {**os.environ, "PYTHONPATH": os.pathsep.join([root, str(pathlib.Path(root).parent / "nb_pylibs")])}
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)
    assert r.returncode == 0, f"stdout={r.stdout}\nstderr={r.stderr[-600:]}"


def test_time_modulated_media_not_silently_ignored():
    """Time-modulated media: either extract them or fail loudly, never ignore them silently.

    History: ``modulation_spec`` is a separate field on ``Medium`` and **does not show up in the
    values of sim.epsilon**, so TimeModulationTutorial was once silently ignored and its
    "comparison complete" verdict was false; after that everything failed closed. Since
    2026-08-19 the scalar CW subset is supported (scene/modulation.py), and this test now locks
    the new invariant: the Scene that gets built has to carry a non-empty modulation, because an
    empty one means we are back to ignoring it silently.
    """
    sim = td.Simulation.from_file(f"{TABLEA}/TimeModulationTutorial/Tidy3D/simulation.json")
    sc = scene_mod.from_simulation(sim)
    assert sc.any_modulation and sc.modulation.n_entry > 0
    assert sc.modulation.freq == 29979245800000.0


#: A case whose only blocker was symmetry to begin with (PlaneWave + Medium + symmetry=(-1,1,0))
SYM_CASE = "AllDielectricStructuralColor"


def _sym_scene():
    """Build the Scene for SYM_CASE, working around a blocker irrelevant to these tests.

    Its source is ``direction='-'`` (a backward plane wave, not implemented yet), but these two
    tests are about **the grid and the cell-count convention under symmetry**, which has nothing
    to do with which way the source travels. Flipping the direction to '+' to build the scene
    keeps the 176,610 that was checked against the reference solver's log.
    """
    sim = td.Simulation.from_file(f"{TABLEA}/{SYM_CASE}/Tidy3D/simulation.json")
    src = sim.sources[0]
    return scene_mod.from_simulation(
        sim.updated_copy(sources=[src.updated_copy(direction="+")]))


def test_symmetry_now_accepted_full_domain():
    """Symmetry boundaries are accepted now: we **run the full domain** and do not reduce it.

    Why: the client's sim.grid is the full-domain grid (measured on CavityFOM with
    symmetry=(1,-1,1): all three axes span the full domain and the symmetry planes sit in the
    middle); only the reference solver reduces the domain internally. Reduction is a performance
    optimization, not a prerequisite for correctness.
    """
    sc = _sym_scene()
    assert sc.symmetry == (-1, 1, 0)
    # full domain: the cell count is exactly the one the client reports, not a quarter of it
    assert sc.shape == (54, 54, 208)


def test_cloud_cell_count_with_symmetry():
    """Cell count in the reference solver's convention when symmetry is present: ∏(N/2 + 2).

    AllDielectricStructuralColor: (54/2+2)(54/2+2)(208+2) = 29*29*210 = 176,610. The reference
    solver's log reports 1.7661e+05, an exact match (the log carries 5 significant digits).

    On a symmetry axis the reference solver computes only half the domain and adds one ghost
    layer at each end (one on the symmetry plane, one on the outer boundary); a non-symmetric
    axis follows the original rule (here z has PML at both ends, +2).
    """
    sc = _sym_scene()
    assert sc.num_cells_cloud == 29 * 29 * 210 == 176_610
    # what we allocate is the full domain, 4x the reference solver's count (two symmetry axes)
    assert sc.grid.num_cells > 3 * sc.num_cells_cloud


def test_cloud_cell_count_without_symmetry(sc):
    """Without symmetry it reduces to "+1 for every end with no neighbour", which is the 49,984
    of BiosensorGrating."""
    assert sc.symmetry == (0, 0, 0)
    assert sc.num_cells_cloud == sc.grid.num_cells_cloud == 49_984
