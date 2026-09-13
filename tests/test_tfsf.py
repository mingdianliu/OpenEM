# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Criteria for the TFSF box. The core one is the **empty-box vacuum criterion** (the lesson
learned: a sign or phase error does not diverge, it only leaks silently, so the criteria come
first):

1. the field inside the box equals, point by point, the incident table from the 1D auxiliary
   grid (same discrete operators, the only difference is float32 rounding);
2. outside the box (the scattered-field region) it is ≈ 0, and whatever leakage is measured gets
   reported;
3. with an infinite substrate crossing the box (the MultipoleExpansion layout), 1 and 2 still
   hold;
4. a scene without TFSF never once touches the new kernels (a single code path);
5. empty-box leakage at oblique incidence (1D line along k̂ plus interpolation over real-valued
   table columns); unsupported backgrounds still fail closed;
6. scene.npz serialization round-trips bitwise.

The scenes are deliberately small (64³ cells, about 1300 steps), so one test takes seconds on an
H800.
"""

from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import scene as scene_mod
from openem import tfsf1d

#: Empty-box leakage threshold. The measured magnitudes are recorded in the assertion messages;
#: the threshold is about 10x the measured value, which locks regressions without being sensitive
#: to floating-point noise.
LEAK_TOL = 1e-4
#: Threshold on the relative deviation between the in-box field and the 1D table (float32
#: propagation vs a float64 table).
MATCH_TOL = 1e-4


def _make_sim(direction: str = "+", substrate_top: float | None = None,
              angle_theta: float = 0.0, pol_angle: float = 0.0,
              injection_axis: int = 2, angle_phi: float = 0.0) -> td.Simulation:
    """A centered TFSF box in vacuum (or with an infinite substrate), uniform grid of 0.02 µm."""
    freq0 = 6.0e14
    structures = []
    if substrate_top is not None:
        structures.append(td.Structure(
            geometry=td.Box.from_bounds((-100, -100, -100), (100, 100, substrate_top)),
            medium=td.Medium(permittivity=2.25)))
    src = td.TFSF(
        center=(0, 0, 0), size=(0.4, 0.4, 0.4),
        source_time=td.GaussianPulse(freq0=freq0, fwidth=freq0 / 6),
        injection_axis=injection_axis, direction=direction,
        angle_theta=angle_theta, pol_angle=pol_angle, angle_phi=angle_phi)
    return td.Simulation(
        size=(0.8, 0.8, 0.8),
        grid_spec=td.GridSpec.uniform(dl=0.02),
        sources=[src],
        structures=structures,
        run_time=5.0e-14,
        boundary_spec=td.BoundarySpec.all_sides(boundary=td.PML()),
    )


def _run(sim, num_steps, kernels=None):
    from openem import solver
    from openem.device import Kernels

    sc = scene_mod.from_simulation(sim)
    res = solver.run(sc, num_steps=num_steps, use_shutoff=False,
                     kernels=kernels or Kernels(), verbose=False, return_fields=True)
    return sc, res


def _leak(sc, res) -> float:
    """Largest |E| in the scattered-field region (between the box grown by 2 cells and the PML
    pulled in by 2 cells), divided by the incident peak."""
    t = sc.tfsf_sources[0]
    ilo, jlo, klo = t.box_lo
    ihi, jhi, khi = t.box_hi
    nx, ny, nz = sc.shape
    pml = {a: (sc.pml[(a, "lo")].num_layers, sc.pml[(a, "hi")].num_layers)
           for a in range(3)}
    mask = np.zeros(sc.shape, dtype=bool)
    mask[pml[0][0] + 2:nx - pml[0][1] - 2,
         pml[1][0] + 2:ny - pml[1][1] - 2,
         pml[2][0] + 2:nz - pml[2][1] - 2] = True
    mask[ilo - 2:ihi + 2, jlo - 2:jhi + 2, klo - 2:khi + 2] = False
    peak = max(float(np.abs(t.ex_inc).max()), 1e-300)
    return max(float(np.abs(res.fields[c][mask]).max())
               for c in ("Ex", "Ey", "Ez")) / peak


def _tables(sc, n_steps):
    t = sc.tfsf_sources[0]
    return tfsf1d.tables(t.dl1, t.eps1, sc.dt, t.inj, t.direction,
                         t.ex_inc, t.hy_inc, n_steps,
                         t.col0, t.box_hi[t.axis] - t.box_lo[t.axis] + 1)


def _check_box(sim):
    """While the pulse is in the box: the in-box column == polarization component × the 1D table.
    Both while the pulse is in the box and after the run: leakage outside the box.

    A deviation is taken for each of the two transverse E components (the component with β = 0
    has to be zero point by point, which the same formula covers).
    """
    sc = scene_mod.from_simulation(sim)
    ex_tab, _ = _tables(sc, sc.num_time_steps)
    t = sc.tfsf_sources[0]
    a = t.axis
    u, v = (a + 1) % 3, (a + 2) % 3
    mid = (t.box_lo[a] + t.box_hi[a]) // 2 - t.box_lo[a]
    n_mid = int(np.argmax(np.abs(ex_tab[:, mid])))     # the step where the pulse peak is at the box center

    _, res = _run(sim, n_mid)
    scale = float(np.abs(ex_tab).max())
    mismatch = 0.0
    for comp, beta in ((u, t.pol_u), (v, t.pol_v)):
        idx = [0, 0, 0]
        idx[u] = (t.box_lo[u] + t.box_hi[u]) // 2
        idx[v] = (t.box_lo[v] + t.box_hi[v]) // 2
        idx[a] = slice(t.box_lo[a], t.box_hi[a] + 1)
        col = res.fields[("Ex", "Ey", "Ez")[comp]][tuple(idx)]
        mismatch = max(mismatch,
                       float(np.abs(col - beta * ex_tab[n_mid]).max()) / scale)
    leak_mid = _leak(sc, res)

    _, res_end = _run(sim, sc.num_time_steps)          # after the pulse has left
    leak_end = _leak(sc, res_end)
    return mismatch, leak_mid, leak_end


@pytest.mark.parametrize("direction", ["+", "-"])
def test_empty_box_vacuum(direction):
    """Empty-box vacuum criterion, checked in both directions.

    Measured (H800, 2026-08-18): in-box vs 1D-table deviation 2.98e-07, leakage 2.53e-07 (pulse
    in the box) / 5.47e-08 (after the run), which is float32 rounding level and the same order
    as the 5.5e-08 one-way-ness of a single TF/SF face.
    """
    pytest.importorskip("cupy")
    mismatch, leak_mid, leak_end = _check_box(_make_sim(direction=direction))
    assert mismatch < MATCH_TOL, f"in-box field differs from the 1D incident table by {mismatch:.3e}"
    assert leak_mid < LEAK_TOL, f"leakage while the pulse is in the box {leak_mid:.3e}"
    assert leak_end < LEAK_TOL, f"leakage after the run {leak_end:.3e}"


def test_empty_box_layered_substrate():
    """An infinite substrate crossing the box (the MultipoleExpansion layout): the incident field
    then contains interface reflection and transmission, and once the 1D auxiliary grid carries
    the ε(z) profile the empty-box criterion still has to hold.

    Measured: deviation 3.60e-07, leakage 2.18e-07 / 5.64e-08, the same order as the empty box in
    vacuum.
    """
    pytest.importorskip("cupy")
    mismatch, leak_mid, leak_end = _check_box(_make_sim(substrate_top=-0.06))
    assert mismatch < MATCH_TOL, f"in-box field differs from the 1D incident table by {mismatch:.3e}"
    assert leak_mid < LEAK_TOL, f"leakage while the pulse is in the box {leak_mid:.3e}"
    assert leak_end < LEAK_TOL, f"leakage after the run {leak_end:.3e}"


def test_no_tfsf_never_touches_new_kernels():
    """A scene without TFSF never calls the new kernels: direct evidence of a single code path.

    How: replace the four tfsf kernels with stubs that blow up on the first call, then run an
    ordinary plane-wave scene.
    """
    pytest.importorskip("cupy")
    from openem.device import Kernels

    k = Kernels()

    def _boom(*a, **kw):
        raise AssertionError("a scene without TFSF called a TFSF kernel")

    for name in list(k._fn):
        if name.startswith("tfsf_"):
            k._fn[name] = _boom

    sim = _make_sim()
    pw = td.PlaneWave(
        center=(0, 0, -0.3), size=(td.inf, td.inf, 0),
        source_time=sim.sources[0].source_time, direction="+")
    sc, res = _run(sim.updated_copy(sources=[pw]), 50, kernels=k)
    assert not sc.tfsf_sources
    assert res.steps_run == 50


#: Empty-box leakage threshold at oblique incidence. The 1D line matches the discrete dispersion
#: only at freq0, so a broadband pulse leaves a residual: measured 2.4e-05 at 25.71° and 1.0e-07
#: at 45°, and the threshold is about 10x the measured value.
OBLIQUE_LEAK_TOL = 3e-4


def test_matched_dl_analytic():
    """The matched cell size of the 1D line agrees with the analytic ``dl₁ = dl·√(Σ k̂ᵢ⁴)``.

    Expanding the two discrete dispersion relations to O(k²dl²) gives this formula (the
    fourth-order term is ``k²dl²Σk̂ᵢ⁴/12`` in 3D and ``k²dl₁²/12`` in 1D). Along a coordinate
    axis, a face diagonal and a body diagonal, Σk̂ᵢ⁴ is 1, 1/2 and 1/3 respectively, so the
    bisection should land **exactly** on dl, dl/√2 and dl/√3 (in those three cases the two
    formulas have the same shape and the result is frequency-independent).
    """
    from openem import cpml, tfsf1d
    dl = 5e-9
    dt = 0.99 * dl / (cpml.C_0 * np.sqrt(3))
    for k in ((1, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1),
              (np.cos(np.pi / 7), np.sin(np.pi / 7), 0)):
        k = np.asarray(k, float) / np.linalg.norm(k)
        got = tfsf1d.matched_dl(k, (dl, dl, dl), dt, 2.99792458e14)
        want = dl * np.sqrt(float(np.sum(k ** 4)))
        assert abs(got - want) < 1e-4 * dl, f"k̂={k}: {got:.6e} vs {want:.6e}"


@pytest.mark.parametrize("kw,tol", [
    (dict(angle_theta=np.pi / 4), LEAK_TOL),                 # face diagonal: the match is exact
    (dict(angle_theta=np.pi / 7), OBLIQUE_LEAK_TOL),
    (dict(angle_theta=np.pi / 7, angle_phi=np.pi / 6), OBLIQUE_LEAK_TOL),
    (dict(angle_theta=0.3, injection_axis=0), OBLIQUE_LEAK_TOL),
])
def test_oblique_empty_box(kw, tol):
    """Empty-box criterion at oblique incidence: leakage outside the box.

    At 45° (the face diagonal) ``dl₁ = dl/√2`` is exact independently of frequency and the
    leakage falls back to the **same order as normal incidence** (measured 1.04e-07 vs
    1.23e-07), which is independent corroboration of the whole geometry and the 24 coefficients.
    The other angles are matched only at freq0, with a residual measured at ~2e-05.
    """
    pytest.importorskip("cupy")
    sc, res = _run(_make_sim(**kw), 400)
    leak = _leak(sc, res)
    assert leak < tol, f"{kw}: leakage {leak:.3e}"


def test_oblique_fails_closed_on_unsupported():
    """Oblique incidence currently handles only a uniform lossless background; anything else
    still fails closed instead of guessing."""
    with pytest.raises(NotImplementedError):     # substrate through the box -> non-uniform background
        scene_mod.from_simulation(
            _make_sim(angle_theta=0.3, substrate_top=-0.1))


def test_oblique_serialize_roundtrip():
    """The oblique geometry (k̂/ê/dl₁/proj0/table column count) is bitwise identical after a
    round trip through npz."""
    import tempfile
    from openem import serialize

    sc = scene_mod.from_simulation(_make_sim(angle_theta=np.pi / 7))
    with tempfile.TemporaryDirectory() as d:
        back = serialize.load(serialize.save(sc, f"{d}/scene.npz"))
    a, b = sc.tfsf_sources[0], back.tfsf_sources[0]
    assert a.k_hat == b.k_hat and a.e_hat == b.e_hat
    assert (a.dl1_step, a.proj0, a.n_col_e) == (b.dl1_step, b.proj0, b.n_col_e)


@pytest.mark.parametrize("kw", [
    dict(injection_axis=0, angle_phi=np.pi / 3),   # the FullyAnisotropic layout
    dict(injection_axis=1, direction="-", pol_angle=0.4),
    dict(injection_axis=2, pol_angle=0.5),
])
def test_empty_box_rotated_pol_injection(kw):
    """Empty-box criterion for an arbitrary injection axis with in-plane rotated polarization.

    Inside the box each of the two transverse E components equals β·(1D table), including the one
    with β=0, which is zero point by point; the leakage outside is float32 rounding level.
    Measured (H800, 2026-08-19): across the three layouts, deviation ≤3.32e-07 and leakage
    ≤2.17e-07 (pulse in the box) / ≤5.4e-08 (after the run), the same order as the older
    criterion for Ex polarization on the z axis.
    """
    pytest.importorskip("cupy")
    mismatch, leak_mid, leak_end = _check_box(_make_sim(**kw))
    assert mismatch < MATCH_TOL, f"in-box field differs from β×(1D table) by {mismatch:.3e}"
    assert leak_mid < LEAK_TOL, f"leakage while the pulse is in the box {leak_mid:.3e}"
    assert leak_end < LEAK_TOL, f"leakage after the run {leak_end:.3e}"


def test_structure_crossing_box_face_fails_closed():
    """A structure cutting into a box face (transverse non-uniformity) has to raise, not silently
    hand back a wrong "incident field"."""
    sim = _make_sim()
    blk = td.Structure(
        geometry=td.Box(center=(0.2, 0, 0), size=(0.3, 0.1, 0.1)),  # straddles the x+ face
        medium=td.Medium(permittivity=4.0))
    with pytest.raises(NotImplementedError):
        scene_mod.from_simulation(sim.updated_copy(structures=[blk]))


def test_runway_guard_fires():
    """Too short a 1D runway has to raise: this is the backstop for trusting the incident table."""
    n1 = 200
    dl1 = np.full(n1, 2e-8)
    eps1 = np.ones(n1)
    dt = 3.8e-17
    n_steps = 2000                      # long enough for the signal to reach the end
    ex = np.ones(n_steps + 1)
    hy = np.ones(n_steps + 1) / 376.73
    with pytest.raises(ValueError):
        tfsf1d.tables(dl1, eps1, dt, 100, +1, ex, hy, n_steps, 110, 20)


def test_serialize_roundtrip():
    """TFSF sources and monitors are bitwise identical after a scene.npz save/load (job
    containers only read npz)."""
    import tempfile
    from openem import serialize

    sim = _make_sim()
    sc = scene_mod.from_simulation(sim)
    with tempfile.TemporaryDirectory() as d:
        p = serialize.save(sc, f"{d}/scene.npz")
        back = serialize.load(p)
    assert len(back.tfsf_sources) == 1
    a, b = sc.tfsf_sources[0], back.tfsf_sources[0]
    assert (a.direction, a.box_lo, a.box_hi, a.inj, a.col0) == \
        (b.direction, b.box_lo, b.box_hi, b.inj, b.col0)
    for f in ("dl1", "eps1", "ex_inc", "hy_inc"):
        np.testing.assert_array_equal(getattr(a, f), getattr(b, f))
    np.testing.assert_array_equal(a.waveform.amp_int, b.waveform.amp_int)
    np.testing.assert_array_equal(a.waveform.amp_half, b.waveform.amp_half)

# ---------------------------------------------------------------- transversely infinite box


def _make_open_sim(angle_theta: float = np.pi / 4, angle_phi: float = np.pi / 4,
                   size_z: float = 0.4) -> td.Simulation:
    """The tfsf2 layout from TFSF.ipynb: ``size=(inf, inf, 0.4)`` with Bloch boundaries
    transversely and PML on the injection axis.

    The box fills the computational domain in x and y, so those four faces do not exist: only the
    two z faces need a TF/SF correction.
    """
    freq0 = 6.0e14
    src = td.TFSF(
        center=(0, 0, 0), size=(td.inf, td.inf, size_z),
        source_time=td.GaussianPulse(freq0=freq0, fwidth=freq0 / 6),
        injection_axis=2, direction="+",
        angle_theta=angle_theta, angle_phi=angle_phi)
    size = (0.8, 0.8, 1.2)
    return td.Simulation(
        size=size, grid_spec=td.GridSpec.uniform(dl=0.02),
        sources=[src], structures=[], run_time=5.0e-14,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.bloch_from_source(source=src, domain_size=size[0], axis=0),
            y=td.Boundary.bloch_from_source(source=src, domain_size=size[1], axis=1),
            z=td.Boundary.pml()))


def test_open_box_axes_detected():
    """Transversely infinite -> recorded as an open axis, box indices spanning the whole axis."""
    sc = scene_mod.from_simulation(_make_open_sim())
    t = sc.tfsf_sources[0]
    assert t.open_axes == (0, 1)
    assert t.box_lo[0] == 0 and t.box_lo[1] == 0
    assert t.box_hi[0] == sc.grid.axes[0].edges.size - 1
    assert t.box_hi[1] == sc.grid.axes[1].edges.size - 1
    assert 0 < t.box_lo[2] < t.box_hi[2] < sc.shape[2]      # z is still a finite box


def test_open_box_only_injection_faces():
    """An open axis has no faces: of the 24 face corrections only the two on the injection axis
    are left."""
    from openem import cpml, tfsf_oblique

    sc = scene_mod.from_simulation(_make_open_sim())
    t = sc.tfsf_sources[0]
    _, e_uva, h_uva = tfsf_oblique.basis(t.axis, t.k_hat, t.e_hat)
    ch = sc.dt / cpml.MU_0
    axes = {term["fixed_ax"] for term in tfsf_oblique.face_terms(
        t.axis, e_uva, h_uva, t.box_lo, t.box_hi, sc.grid, ch, t.open_axes)}
    assert axes == {t.axis}


def test_open_box_serialize_roundtrip():
    """open_axes has to survive a round trip through npz."""
    import tempfile
    from openem import serialize

    sc = scene_mod.from_simulation(_make_open_sim())
    with tempfile.TemporaryDirectory() as d:
        back = serialize.load(serialize.save(sc, f"{d}/scene.npz"))
    assert back.tfsf_sources[0].open_axes == sc.tfsf_sources[0].open_axes


def test_open_box_normal_incidence_fails_closed():
    """The normal-incidence path uses the 16 hand-written terms in setup_tables and has not been
    split up per face yet: say so clearly, then block it."""
    with pytest.raises(NotImplementedError, match="normal-incidence TFSF box has infinite size"):
        scene_mod.from_simulation(_make_open_sim(angle_theta=0.0, angle_phi=0.0))


def test_open_box_infinite_injection_axis_fails_closed():
    """If the injection axis itself is infinite there is no total-field/scattered-field interface.

    tidy3d also rejects such a source when the Simulation is constructed ("must not touch or
    cross the simulation boundary along its injection axis"), so this tests ``_tfsf_box_edges``
    directly: a scene read back from a file bypasses tidy3d's validation, so we need this gate
    of our own.
    """
    from types import SimpleNamespace
    from openem.scene.sources import _tfsf_box_edges

    sc = scene_mod.from_simulation(_make_open_sim())
    src = SimpleNamespace(
        bounds=((-np.inf, -0.2, -np.inf), (np.inf, 0.2, np.inf)),
        injection_axis=2)
    with pytest.raises(NotImplementedError, match="injection axis"):
        _tfsf_box_edges(src, sc.grid)
