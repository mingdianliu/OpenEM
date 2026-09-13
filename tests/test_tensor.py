# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Full-tensor ε/σ (FullyAnisotropicMedium) and UniformCurrentSource.

Every criterion has zero tunable parameters:

1. **Diagonal degeneracy**: a diagonal FullyAnisotropicMedium and an AnisotropicMedium build
   scenes with the same ε cell by cell and the same dt, so the fields from the two update paths
   (the tensor trio vs the scalar update_e) differ only by float32 rounding.
2. **Analytic amplitude of a current sheet**: in vacuum an infinite current sheet ``K = amp A/µm``
   has the plane-wave solution ``|E(ω)| = η₀|K(ω)|/2``, which locks the amplitude convention of
   UniformCurrentSource.
3. **Rotation superposition identity**: in a transversely uniform scene a diagonal medium rotated
   about the propagation axis decouples **exactly** in the rotated basis (4-neighbor averaging is
   the identity on a uniform field), so the tensor run has to equal a trigonometric combination of
   two scalar runs. This locks the coefficients and the signs of the off-diagonal coupling, once
   for each of the three axes.
4. **Faraday rotation**: the polarization rotation angle of an antisymmetric (gyrotropic) σ has a
   closed form ``ρ = π/λ·(n₋ − n₊)``, compared point by point against the rotation factor applied
   to an unperturbed reference run.

Every tolerance was measured first and written down afterwards: the identity criteria sit at
float32 rounding level, the analytic ones at discretization level (values given in each test).
"""

from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from dataclasses import replace

import pytest

# This module imports openem.solver at the top level, and that imports cupy at its own top level,
# so on a CPU-only machine collection itself would crash.
pytest.importorskip("cupy", reason="needs CUDA and cupy")

from openem import scene as scene_mod
from openem import serialize, solver
from openem.model import TensorEps

FREQ0 = 2.0e14


@pytest.fixture(scope="module")
def kernels():
    from openem.device import Kernels
    return Kernels()


def _run(sc, steps, kernels):
    return solver.run(sc, num_steps=steps, use_shutoff=False,
                      kernels=kernels, verbose=False)


# ---------------------------------------------------------------------------
# 1 · diagonal degeneracy
# ---------------------------------------------------------------------------

def _slab_pair():
    """The same geometry: a slab of diagonal FullyAnisotropicMedium vs AnisotropicMedium."""
    diag_e = (2.0, 3.0, 4.0)
    diag_s = (1.0e-4, 2.0e-4, 3.0e-4)     # S/µm; only a lossy medium exercises the semi-implicit half
    full = td.FullyAnisotropicMedium(
        permittivity=np.diag(diag_e).tolist(),
        conductivity=np.diag(diag_s).tolist())
    aniso = td.AnisotropicMedium(
        xx=td.Medium(permittivity=diag_e[0], conductivity=diag_s[0]),
        yy=td.Medium(permittivity=diag_e[1], conductivity=diag_s[1]),
        zz=td.Medium(permittivity=diag_e[2], conductivity=diag_s[2]))
    def sim_of(medium):
        return td.Simulation(
            size=(0.4, 0.4, 2.4),
            grid_spec=td.GridSpec.uniform(dl=0.05),
            medium=td.Medium(permittivity=1.0),
            structures=[td.Structure(
                geometry=td.Box(center=(0, 0, 0.3), size=(td.inf, td.inf, 0.6)),
                medium=medium)],
            sources=[td.PlaneWave(
                center=(0, 0, -0.8), size=(td.inf, td.inf, 0), direction="+",
                source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.3 * FREQ0))],
            monitors=[td.FieldMonitor(center=(0, 0, 0.3), size=(0.2, 0.2, 0.4),
                                      freqs=[FREQ0], name="probe")],
            run_time=2.0e-13,
            boundary_spec=td.BoundarySpec.all_sides(td.Periodic()),
            subpixel=False)
    return sim_of(full), sim_of(aniso)


def test_diagonal_degenerate_matches_anisotropic(kernels):
    """A diagonal tensor == the existing AnisotropicMedium path.

    The two sims have the same dt (locked by an assertion) and both extract ε/σ by staircased
    point sampling, so the only difference is which update path runs. Tolerance 2e-6·|E|max;
    measured 3.2e-7 (float32 rounding: M1/M2 are rows of an inverted 3×3 while the scalar path
    divides directly, so the last digit differs).
    """
    sim_full, sim_aniso = _slab_pair()
    assert sim_full.dt == sim_aniso.dt
    sc_f = scene_mod.from_simulation(sim_full)
    sc_a = scene_mod.from_simulation(sim_aniso)
    assert sc_f.any_tensor and not sc_a.any_tensor
    # in the scalar arrays a tensor cell holds "the background with the structure removed"
    # (vacuum): the pass-through overwrites it, it never enters the physics
    r_f = _run(sc_f, 1500, kernels)
    r_a = _run(sc_a, 1500, kernels)
    a = r_a.field_phasors["probe"]
    b = r_f.field_phasors["probe"]
    scale = float(np.max(np.abs(a)))
    assert scale > 0
    diff = float(np.max(np.abs(a - b))) / scale
    print(f"[measure] diagonal degenerate residual = {diff:.3e}")
    assert diff < 2.0e-6, f"diagonal degeneracy residual {diff:.3e}"


# ---------------------------------------------------------------------------
# 2 · analytic amplitude of a current sheet (the amplitude convention of UniformCurrentSource)
# ---------------------------------------------------------------------------

def test_sheet_current_analytic_amplitude(kernels):
    """A current sheet in vacuum: ``|Ex(ω)| = η₀·|K(ω)|/2``, with K = amp A/µm.

    The phase is locked as well: ``Ex(ω, z) = −(η₀/2)K(ω)e^{ik|z−z₀|}``, where k is the ``k̃`` of
    the **discrete** dispersion relation. All that is left in the residual is then the half-cell
    filtering from spreading the sheet over two planes (cos(kΔ/2) ~ 1e-4), DFT truncation and
    residual PML reflection. Tolerance 1e-2; measured 5.6e-3.
    """
    from openem import cpml
    sim = td.Simulation(
        size=(0.2, 0.2, 6.0), grid_spec=td.GridSpec.uniform(dl=0.05),
        medium=td.Medium(permittivity=1.0),
        sources=[td.UniformCurrentSource(
            center=(0, 0, -2.0), size=(td.inf, td.inf, 0), polarization="Ex",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.1 * FREQ0),
            current_amplitude_definition="density")],
        monitors=[td.FieldMonitor(center=(0, 0, 0), size=(0, 0, td.inf),
                                  freqs=[FREQ0], name="line")],
        run_time=6.0e-13,
        boundary_spec=td.BoundarySpec(x=td.Boundary.periodic(),
                                      y=td.Boundary.periodic(),
                                      z=td.Boundary.pml()))
    sc = scene_mod.from_simulation(sim)
    # sheet exactly on an Ex z edge -> one plane, 16 points; between two edges -> two planes, 32
    assert len(sc.dipoles) == 1 and sc.dipoles[0].indices.shape[0] in (16, 32)
    n_steps = sc.num_time_steps
    res = _run(sc, n_steps, kernels)

    d = sc.dipoles[0]
    dt = sc.dt
    n = np.arange(n_steps)
    om = 2 * np.pi * FREQ0
    # We inject with amp_half (the half-step instants) and the DFT takes E at (n+1)dt: the same
    # convention on both sides.
    # K(ω) = ∫K(t)e^{iωt}dt, K = amp·1e6 A/m (a current sheet in the density convention, see the
    # derivation in sources.uniform_current); the phasor buffer carries no dt, so put it back in
    # before comparing.
    k_spec = 1.0e6 * dt * np.sum(
        d.waveform.amp_half[:n_steps] * np.exp(1j * om * (n + 0.5) * dt))
    # discrete dispersion: sin(k̃ dz/2) = (dz/(c dt))·sin(ω dt/2) (vacuum)
    dz = float(sc.grid.axes[2].dl[0])
    k_num = 2.0 / dz * np.arcsin(np.clip(
        dz / (cpml.C_0 * dt) * np.sin(om * dt / 2.0), -1.0, 1.0))

    mon = sc.field_monitors[0]
    ours = dt * res.field_phasors[mon.name][0, 0, 0, 0, :]   # Ex along z, dt put back in
    zc = sc.grid.axes[2].edges[mon.origin[2]:mon.origin[2] + mon.box[2]]
    z0 = -2.0e-6
    pred = -(cpml.ETA_0 / 2.0) * k_spec * np.exp(1j * k_num * np.abs(zc - z0))
    # compare only the middle stretch, away from the PML and the source
    n_pml = 12 * 2
    sel = (zc > z0 + 0.5e-6) & (zc < zc[-n_pml - 1] - 0.5e-6)
    rel = np.abs(ours[sel] - pred[sel]) / (cpml.ETA_0 / 2.0 * abs(k_spec))
    print(f"[measure] sheet analytic residual = {float(rel.max()):.3e}")
    assert float(rel.max()) < 1.0e-2, f"current sheet analytic residual {rel.max():.3e}"


# ---------------------------------------------------------------------------
# 3 · rotation superposition identity (off-diagonal coupling coefficients and signs, once per axis)
# ---------------------------------------------------------------------------

def _rot(axis: int, th: float) -> np.ndarray:
    c, s = np.cos(th), np.sin(th)
    if axis == 2:
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    if axis == 0:
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def _uniform_tensor_scene(prop_axis: int, eps_main: tuple, sig_main: tuple,
                          th: float):
    """A uniform tensor background rotated by ``th`` about the propagation axis ``prop_axis``,
    plus a current sheet.

    The source polarization is the first transverse axis t1 and the monitor line runs along the
    propagation axis. Returns (sc, t1, t2, mon_name).
    """
    t1, t2 = [a for a in range(3) if a != prop_axis]
    r = _rot(prop_axis, th)
    perm = r @ np.diag(eps_main) @ r.T
    cond = r @ np.diag(sig_main) @ r.T
    size = [0.2, 0.2, 0.2]
    size[prop_axis] = 6.0
    src_center = [0.0, 0.0, 0.0]
    src_center[prop_axis] = -2.0
    src_size = [td.inf, td.inf, td.inf]
    src_size[prop_axis] = 0.0
    mon_size = [0.0, 0.0, 0.0]
    mon_size[prop_axis] = td.inf
    bnd = {a: td.Boundary.periodic() for a in "xyz"}
    bnd["xyz"[prop_axis]] = td.Boundary.pml()
    sim = td.Simulation(
        size=tuple(size), grid_spec=td.GridSpec.uniform(dl=0.05),
        medium=td.FullyAnisotropicMedium(permittivity=perm.tolist(),
                                         conductivity=cond.tolist()),
        sources=[td.UniformCurrentSource(
            center=tuple(src_center), size=tuple(src_size),
            polarization="E" + "xyz"[t1],
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.1 * FREQ0),
            current_amplitude_definition="density")],
        monitors=[td.FieldMonitor(center=(0, 0, 0), size=tuple(mon_size),
                                  freqs=[FREQ0], name="line")],
        run_time=4.0e-13,
        boundary_spec=td.BoundarySpec(**bnd))
    return scene_mod.from_simulation(sim), t1, t2


def _scalar_ref(sc, eps_val: float, sig_val: float):
    """A scalar reference built from the same Scene: ε/σ flattened to constants, tensor path off.

    dt, the waveform tables, the PML coefficients and the source are all shared unchanged, so the
    identity holds **exactly** (item 3 of the module docstring) and is not disturbed by tidy3d
    picking a different dt for a different medium.
    """
    const = np.full(sc.eps_ex.shape, eps_val, dtype=np.float64)
    kw = dict(tensor=None, eps_ex=const, eps_ey=const.copy(), eps_ez=const.copy())
    if sig_val:
        s = np.full(sc.eps_ex.shape, sig_val, dtype=np.float64)
        kw.update(sigma_ex=s, sigma_ey=s.copy(), sigma_ez=s.copy())
    return replace(sc, **kw)


@pytest.mark.parametrize("prop_axis", [2, 0, 1])
def test_rotation_superposition_identity(prop_axis, kernels):
    """A diagonal medium rotated by θ about the propagation axis: the tensor run == a
    trigonometric combination of two scalar runs.

    Transverse uniformity ⇒ 4-neighbor averaging is the identity ⇒ exact decoupling in the
    rotated basis:

        E_t1 = cos²θ·E(ε₁) + sin²θ·E(ε₂)
        E_t2 = cosθ·sinθ·(E(ε₁) − E(ε₂))

    The z axis also carries a diagonal σ (rotated as well), and the three axes lock the signs of
    the xy/yz/xz coupling pairs. Tolerance 2e-5·|E|max; measured ≤ 7.6e-7 on all three axes
    (float32 rounding plus the difference between an inverted row and a direct division).
    """
    th = np.deg2rad(30.0)
    eps_main = (2.0, 3.0, 4.0)
    sig_main = (1.0e-4, 2.0e-4, 0.0) if prop_axis == 2 else (0.0, 0.0, 0.0)
    sc, t1, t2 = _uniform_tensor_scene(prop_axis, eps_main, sig_main, th)
    assert sc.any_tensor
    steps = min(sc.num_time_steps, 2500)
    res_t = _run(sc, steps, kernels)
    e_main = [eps_main[t1], eps_main[t2]]
    s_main = [sig_main[t1], sig_main[t2]]
    r1 = _run(_scalar_ref(sc, e_main[0], s_main[0] * 1.0e6), steps, kernels)
    r2 = _run(_scalar_ref(sc, e_main[1], s_main[1] * 1.0e6), steps, kernels)

    ph_t = res_t.field_phasors["line"]
    ph_1 = r1.field_phasors["line"]
    ph_2 = r2.field_phasors["line"]
    # Components in (t1, t2) of the rotated principal axis u = R·t̂1. Do **not** simply write
    # (cosθ, sinθ): in the (x, z) ordering R_y has the opposite handedness (the cyclic order is
    # (z, x)), so the sign flips.
    u = _rot(prop_axis, th)[:, t1]
    c, s = float(u[t1]), float(u[t2])
    want_t1 = c * c * ph_1[t1] + s * s * ph_2[t1]
    want_t2 = c * s * (ph_1[t1] - ph_2[t1])
    scale = float(np.max(np.abs(ph_1[t1])))
    assert scale > 0
    d1 = float(np.max(np.abs(ph_t[t1] - want_t1))) / scale
    d2 = float(np.max(np.abs(ph_t[t2] - want_t2))) / scale
    print(f"[measure] superposition axis={prop_axis}: d1={d1:.3e} d2={d2:.3e}")
    assert d1 < 2.0e-5 and d2 < 2.0e-5, f"superposition identity residual {d1:.3e} / {d2:.3e}"


# ---------------------------------------------------------------------------
# 4 · Faraday rotation (gyrotropic: antisymmetric σ)
# ---------------------------------------------------------------------------

def test_faraday_rotation_analytic(kernels):
    """Polarization rotation in a magneto-optic medium, compared point by point against the
    rotation factor applied to an unperturbed reference run.

    The closed form from the Gyrotropic notebook: ``ρ = π/λ·(n₋ − n₊)`` with
    ``n± = √(ε₀ᵉ ± g)``; after a distance d, ``Ex → Ex·cos(ρd)`` and ``Ey → −Ex·sin(ρ|d|)``.
    The formula is for a continuous medium, so the tolerance is at discretization level: with
    dl=0.05 (about 21 cells per wavelength inside the medium) over 10 µm of propagation the
    measured values are dx=3.0e-2 / dy=2.4e-2 (the phase error accumulates with distance), so the
    bound is 5e-2.
    """
    eps0 = 2.0
    g = 0.1
    wvl = 1.4989622900000002    # C_0/FREQ0, µm
    cond = np.array([[0, g, 0], [-g, 0, 0], [0, 0, 0]]) * 2 * np.pi * FREQ0 \
        * td.EPSILON_0
    sim = td.Simulation(
        size=(0.2, 0.2, 10.0), grid_spec=td.GridSpec.uniform(dl=0.05),
        medium=td.FullyAnisotropicMedium(
            permittivity=np.diag([eps0] * 3).tolist(), conductivity=cond.tolist()),
        sources=[td.UniformCurrentSource(
            center=(0, 0, -4.0), size=(td.inf, td.inf, 0), polarization="Ex",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.1 * FREQ0),
            current_amplitude_definition="density")],
        monitors=[td.FieldMonitor(center=(0, 0, 0), size=(0, 0, td.inf),
                                  freqs=[FREQ0], name="line")],
        run_time=8.0e-13,
        boundary_spec=td.BoundarySpec(x=td.Boundary.periodic(),
                                      y=td.Boundary.periodic(),
                                      z=td.Boundary.pml()))
    sc = scene_mod.from_simulation(sim)
    steps = sc.num_time_steps
    res_g = _run(sc, steps, kernels)
    res_0 = _run(_scalar_ref(sc, eps0, 0.0), steps, kernels)

    mon = sc.field_monitors[0]
    ex_g = res_g.field_phasors["line"][0, 0, 0, 0, :]
    ey_g = res_g.field_phasors["line"][1, 0, 0, 0, :]
    ex_0 = res_0.field_phasors["line"][0, 0, 0, 0, :]

    rho = np.pi / (wvl * 1e-6) * (np.sqrt(eps0 - g) - np.sqrt(eps0 + g))
    zc = sc.grid.axes[2].edges[mon.origin[2]:mon.origin[2] + mon.box[2]]
    dist = zc - (-4.0e-6)
    th_x = ex_0 * np.cos(rho * dist)
    th_y = -ex_0 * np.sin(rho * np.abs(dist))
    n_pml = 12 * 2
    sel = (zc > -3.5e-6) & (zc < zc[-n_pml - 1] - 0.5e-6)
    scale = float(np.max(np.abs(ex_0[sel])))
    dx = float(np.max(np.abs(ex_g[sel] - th_x[sel]))) / scale
    dy = float(np.max(np.abs(ey_g[sel] - th_y[sel]))) / scale
    print(f"[measure] faraday residual dx={dx:.3e} dy={dy:.3e}")
    assert dx < 5.0e-2 and dy < 5.0e-2, f"Faraday rotation residual {dx:.3e} / {dy:.3e}"
    # the rotation really happened (rather than both sides degenerating to 0)
    assert float(np.max(np.abs(th_y[sel]))) / scale > 0.2


# ---------------------------------------------------------------------------
# 5 · structural: serialization, closure validation, fail-closed combinations
# ---------------------------------------------------------------------------

def test_serialize_roundtrip_tensor(tmp_path):
    sc, _, _ = _uniform_tensor_scene(2, (2.0, 3.0, 4.0), (0.0, 0.0, 0.0),
                                     np.deg2rad(30.0))
    p = serialize.save(sc, tmp_path / "scene.npz")
    back = serialize.load(p)
    assert back.any_tensor
    np.testing.assert_array_equal(back.tensor.comp, sc.tensor.comp)
    np.testing.assert_array_equal(back.tensor.cell, sc.tensor.cell)
    np.testing.assert_array_equal(back.tensor.eps, sc.tensor.eps)
    np.testing.assert_array_equal(back.tensor.sigma, sc.tensor.sigma)


def test_tensor_closure_validation():
    """A broken closure has to stop at validate: a missing entry among the 4-neighbor stencil
    points ⇒ ValueError."""
    from openem.grid import Axis, Grid
    edges = np.linspace(0.0, 4e-7, 5)
    grid = Grid(*(Axis(edges, "Periodic", "Periodic") for _ in range(3)))
    eps = np.tile(np.array([[2.0, 0.5, 0], [0.5, 3.0, 0], [0, 0, 4.0]]), (1, 1, 1))
    t = TensorEps(comp=np.array([0]), cell=np.array([0]),
                  eps=eps, sigma=np.zeros((1, 3, 3)))
    tabs = [a.index_tables() for a in grid.axes]
    with pytest.raises(ValueError, match="closure"):
        t.validate(grid.shape, tabs)


def test_tensor_plus_dispersive_fails_closed():
    sim, _ = _slab_pair()
    sim = sim.updated_copy(structures=list(sim.structures) + [td.Structure(
        geometry=td.Box(center=(0, 0, -0.6), size=(td.inf, td.inf, 0.2)),
        medium=td.Lorentz(eps_inf=2.0, coeffs=((1.0, 3e14, 1e13),)))])
    with pytest.raises(NotImplementedError, match="dispersive"):
        scene_mod.from_simulation(sim)


def test_tensor_plus_pec_fails_closed():
    sim, _ = _slab_pair()
    sim = sim.updated_copy(structures=list(sim.structures) + [td.Structure(
        geometry=td.Box(center=(0, 0, -0.6), size=(0.1, 0.1, 0.1)),
        medium=td.PECMedium())])
    with pytest.raises(NotImplementedError, match="PEC"):
        scene_mod.from_simulation(sim)
