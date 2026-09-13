# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""The k-space far-field projection monitor (FieldProjectionKSpaceMonitor).

Three layers of criteria:

* **Expansion** (pure CPU): apodization reaches every face monitor unchanged; a plane monitor takes
  its outward normal from its own ``normal_dir``; a ``proj_axis`` that disagrees with the
  zero-thickness axis fails closed.
* **Vacuum dipole vs the analytic far field** (GPU, seconds): the ``|Eθ|`` of a closed-box KSpace
  monitor along ``(ux, uy=0)`` matches the ideal dipole formula, in the same scene and at the same
  order of tolerance as the angle-domain path, which reached 0.19%.
* **KSpace vs Angle agreement** (within one run): over the same set of directions
  (``ux = sinθ``) the two monitors give far fields that agree with each other; the two paths differ
  only in how the observation grid is parameterized.
"""

from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import scene as scene_mod
from openem.scene import projection as pj

FREQ0 = 2.0e14
ETA0, C0 = 376.73031366857, 299792458.0


# ------------------------------------------------------------ expansion (pure CPU)

def _kspace_plane(**kw):
    args = dict(center=(0, 0, 0.21), size=(td.inf, td.inf, 0), freqs=[FREQ0],
                name="n2f", ux=[-0.5, 0.0, 0.5], uy=[-0.5, 0.0, 0.5],
                proj_axis=2, proj_distance=1e6)
    args.update(kw)
    return td.FieldProjectionKSpaceMonitor(**args)


def _sim(monitors=()):
    return td.Simulation(
        size=(1.2, 1.2, 1.2), grid_spec=td.GridSpec.uniform(dl=0.04),
        sources=[td.PointDipole(
            center=(0, 0, 0), polarization="Ez",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.2 * FREQ0))],
        monitors=list(monitors), run_time=8e-14,
        boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=10)))


def test_is_projection_recognizes_kspace():
    assert pj.is_projection(_kspace_plane())


def test_plane_face_carries_apodization_normal_dir_and_clips_inf():
    apod = td.ApodizationSpec(start=4e-13, width=2e-13)
    sim = _sim()
    mon = _kspace_plane(apodization=apod, normal_dir="-")
    faces = pj.surface_monitors(mon, sim)
    assert len(faces) == 1
    fm, sign, ax = faces[0]
    assert (sign, ax) == ("-", 2)
    assert fm.apodization.start == 4e-13 and fm.apodization.width == 2e-13
    # A size=inf plane is clipped to a finite one: **the whole grid including the PML** minus the
    # outermost cell. Tidy3D discretizes an inf monitor over the full grid including the PML, and
    # the PML tail acts as a soft window for the projection integral
    for a, axn in enumerate("xy"):
        b = np.asarray(getattr(sim.grid.boundaries, axn))
        # After clipping to the grid, pull in by another 1e-6 µm (see projection.surface_monitors:
        # when a corner lands exactly on a grid line the projector's interpolation returns NaN)
        assert fm.bounds[0][a] == pytest.approx(b[1] + 1e-6, abs=1e-9)
        assert fm.bounds[1][a] == pytest.approx(b[-2] - 1e-6, abs=1e-9)
        assert fm.size[a] > sim.size[a]          # it really does reach into the PML
    assert fm.size[2] == 0.0


def test_box_faces_carry_apodization():
    apod = td.ApodizationSpec(start=4e-13, width=2e-13)
    mon = td.FieldProjectionKSpaceMonitor(
        center=(0, 0, 0), size=(1, 1, 1), freqs=[FREQ0], name="n2f",
        ux=[0.0], uy=[0.0], proj_axis=2, proj_distance=1e6, apodization=apod)
    faces = pj.surface_monitors(mon, _sim())
    assert len(faces) == 6
    for fm, _, _ in faces:
        assert fm.apodization.start == 4e-13 and fm.apodization.width == 2e-13


def test_proj_axis_mismatch_fails_closed():
    mon = td.FieldProjectionKSpaceMonitor(
        center=(0, 0.21, 0), size=(td.inf, 0, td.inf), freqs=[FREQ0], name="n2f",
        ux=[0.0], uy=[0.0], proj_axis=2, proj_distance=1e6)
    with pytest.raises(NotImplementedError, match="proj_axis"):
        pj.surface_monitors(mon, _sim())


def test_scene_build_accepts_apodized_kspace_monitor():
    """This used to raise NotImplementedError outright, which is where OptimizedL3 got stuck."""
    sim = _sim([_kspace_plane(center=(0, 0, 0.3),
                              apodization=td.ApodizationSpec(start=4e-13,
                                                             width=2e-13))])
    sc = scene_mod.from_simulation(sim)
    assert len(sc.projection_monitors) == 1
    (fm,) = sc.field_monitors
    assert fm.apodization == (4e-13, None, 2e-13)


def test_serialize_roundtrip_keeps_projection_monitor(tmp_path):
    """A solver job rebuilds the Scene from an npz: the projection monitor, its name and its
    normals, must come back unchanged."""
    from openem import serialize
    sim = _sim([_kspace_plane(center=(0, 0, 0.3),
                              apodization=td.ApodizationSpec(start=4e-13,
                                                             width=2e-13))])
    sc = scene_mod.from_simulation(sim)
    back = serialize.load(serialize.save(sc, tmp_path / "s.npz"))
    bp, p = back.projection_monitors[0], sc.projection_monitors[0]
    assert (bp.name, bp.surface_names, bp.normal_dirs) == \
        (p.name, p.surface_names, p.normal_dirs)
    assert back.field_monitors[0].apodization == sc.field_monitors[0].apodization


# ------------------------------------------- vacuum dipole (GPU, whole scene in seconds)

#: Sampling directions: ux = sinθ (uy = 0, the φ = 0 half plane). θ=0 is avoided, since Eθ goes to
#: 0 there and a relative difference would be meaningless
UX = np.array([0.1, 0.3, 0.5, 0.7])


@pytest.fixture(scope="module")
def vacuum_far():
    """One run, feeding the same near field to both the KSpace and the Angle monitor."""
    pytest.importorskip("cupy")
    from openem import solver
    from openem.device import Kernels

    th = np.arcsin(UX)
    amp = 1e-12
    kmon = td.FieldProjectionKSpaceMonitor(
        center=(0, 0, 0), size=(1.0, 1.0, 1.0), freqs=[FREQ0], name="k",
        ux=list(UX), uy=[0.0], proj_axis=2, proj_distance=1e6,
        far_field_approx=True)
    amon = td.FieldProjectionAngleMonitor(
        center=(0, 0, 0), size=(1.0, 1.0, 1.0), freqs=[FREQ0], name="a",
        theta=list(th), phi=[0.0], proj_distance=1e6, far_field_approx=True)
    sim = td.Simulation(
        size=(1.6, 1.6, 1.6), grid_spec=td.GridSpec.uniform(dl=0.03),
        medium=td.Medium(permittivity=1.0),
        sources=[td.PointDipole(
            center=(0, 0, 0), polarization="Ez",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.5 * FREQ0,
                                         amplitude=amp))],
        monitors=[kmon, amon], run_time=3e-13,
        boundary_spec=td.BoundarySpec.all_sides(td.PML()))
    sc = scene_mod.from_simulation(sim)
    res = solver.run(sc, kernels=Kernels(), verbose=False)

    out = {}
    for pm in sc.projection_monitors:
        proj_mon = next(m for m in sim.monitors if m.name == pm.name)
        faces = pj.surface_monitors(proj_mon, sim)
        specs = {m.name: m for m in sc.field_monitors}
        datas = [pj.field_data(sim, fm, specs[fm.name], res.field_phasors[fm.name])
                 for fm, _, _ in faces]
        out[pm.name] = pj.project(sim, proj_mon, faces, datas)
    return out


def test_kspace_matches_analytic_dipole(vacuum_far):
    """|Eθ(ux)| against the ideal dipole formula. The angle domain reached 0.19% before, so the
    tolerance is set at the same order, 1%."""
    r = 1e6 * 1e-6                       # proj_distance µm -> m
    k = 2 * np.pi * FREQ0 / C0
    il = 1e-6                            # amplitude 1e-12 is already divided out by the source
                                         # spectrum normalization; current moment 1 A·µm
    sinth = UX                           # with uy=0, sinθ = ux
    ana = ETA0 * k * il * sinth / (4 * np.pi * r) * 1e-6   # far field r in the µm convention
    ours = np.abs(np.asarray(vacuum_far["k"].Etheta).squeeze())
    rel = np.abs(ours - ana) / ana
    assert rel.max() < 1e-2, f"KSpace vs analytic: {rel}"


def test_kspace_matches_angle_monitor(vacuum_far):
    """Over the same set of directions the KSpace and Angle paths agree (same near field, same
    client)."""
    a = np.abs(np.asarray(vacuum_far["k"].Etheta).squeeze())
    b = np.abs(np.asarray(vacuum_far["a"].Etheta).squeeze())
    rel = np.abs(a - b) / b
    assert rel.max() < 1e-3, f"KSpace vs Angle: {rel}"
