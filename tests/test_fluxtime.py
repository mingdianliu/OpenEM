# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""FluxTimeMonitor: instantaneous flux in the time domain (colocated-integration convention).
Every criterion has zero tunable parameters.

1. Converter convention: km+frac along the normal, the monitor's boundary coordinates in the
   tangential directions, and a time window equal to Tidy3D's ``time_inds``.
2. **Integration oracle**: on random colocated fields the area elements and the integral agree
   bitwise with the Tidy3D client's ``FieldTimeData.flux`` (all three normal axes, a non-uniform
   grid, monitor bounds off the grid points).
3. **Colocation oracle**: on affine fields the interpolation in ``colocated_component`` is exact
   (linear interpolation has no error on an affine function).
4. Step-by-step agreement: the GPU series ≡ samples from a FieldTimeMonitor on the same plane put
   through the same host colocation/integration.
5. Parseval: ``Σ_m s(m) = 4·dt·∫₀^∞ P(f) df``, with the frequency side computing the phasor flux
   from the **same** colocation geometry. That checks only the time mechanism (the two launches,
   H̄, the DFT); the spatial convention is checked by 2 and 3.
6. A single code path: removing the FluxTimeMonitor leaves the frequency-domain phasors bitwise
   unchanged.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import flux as flux_mod
from openem import scene as scene_mod
from openem import serialize
from openem.grid import Axis, Grid
from openem.model import COMPONENTS, FluxTimeMonitor
from openem.units import UM

FREQ0 = 3.0e14
FWIDTH = 0.6e14
#: Frequency grid for Parseval: the source power spectrum is down to e^{-20} at ±4.5·fwidth, and
#: the grid is dense enough for the trapezoid rule.
FREQS = np.linspace(FREQ0 - 4.5 * FWIDTH, FREQ0 + 4.5 * FWIDTH, 181)
PLANE_Z = 0.3
#: Transverse geometry of the flux plane in the small scene: the bounds deliberately miss the
#: 0.025 grid points.
FT_CENTER = (0.013, -0.017)
FT_SIZE = (0.27, 0.29)


def _sim(monitors) -> td.Simulation:
    """Small vacuum plane-wave scene: periodic in x/y, propagating along z, monitors near
    z=PLANE_Z."""
    per = td.Boundary(plus=td.Periodic(), minus=td.Periodic())
    return td.Simulation(
        size=(0.4, 0.4, 1.6),
        grid_spec=td.GridSpec.uniform(dl=0.025),
        sources=[td.PlaneWave(
            center=(0, 0, -0.55), size=(td.inf, td.inf, 0), direction="+",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=FWIDTH))],
        monitors=monitors,
        run_time=1.5e-13,
        boundary_spec=td.BoundarySpec(
            x=per, y=per, z=td.Boundary(plus=td.PML(), minus=td.PML())),
    )


def _ftm(**kw):
    args = dict(center=(*FT_CENTER, PLANE_Z), size=(*FT_SIZE, 0),
                name="ft", interval=1)
    args.update(kw)
    return td.FluxTimeMonitor(**args)


def _yee_coords(grid: Grid, comp: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Yee coordinates of a component: E_d sits at the cell center along d and on the boundary
    otherwise; H_d is the other way round."""
    kind, d = comp[0], "xyz".index(comp[1])
    out = []
    for i, ax in enumerate(grid.axes):
        on_center = (i == d) if kind == "E" else (i != d)
        out.append(ax.centers if on_center else ax.edges[:-1])
    return tuple(out)


# ------------------------------------------------------------ converter

def test_converter_geometry():
    """km+frac along the normal agrees with the plane of FluxMonitor; the tangential directions
    store the monitor's boundary coordinates in meters."""
    fx = td.FluxMonitor(center=(*FT_CENTER, PLANE_Z), size=(*FT_SIZE, 0),
                        freqs=[FREQ0], name="fx")
    sim = _sim([_ftm(), fx])
    sc = scene_mod.from_simulation(sim)
    m = sc.flux_time_monitors[0]
    assert (m.axis, m.plane_index) == (sc.flux_monitors[0].axis,
                                       sc.flux_monitors[0].plane_index)
    assert m.normal_dir == +1 and m.frac == 0.0      # the plane lands exactly on a boundary
    lo, hi = sim.monitors[0].bounds
    for p, t in enumerate((0, 1)):
        np.testing.assert_allclose(
            m.t_bounds[p], (lo[t] * UM, hi[t] * UM), rtol=1e-12)
    beg, end = (int(v) for v in sim.monitors[0].time_inds(np.asarray(sim.tmesh)))
    assert (m.step_begin, m.step_end, m.interval) == (beg, end, 1)
    assert m.num_slots == end - beg


def test_converter_off_grid_plane_gets_frac():
    """A plane that is not on a Yee boundary gets a normal interpolation weight (all 6 planes of
    NanobeamCavity are like this).

    With dl=0.025, z=0.3125 sits exactly between the edges 0.3 and 0.325 -> frac=0.5.
    """
    on = scene_mod.from_simulation(_sim([_ftm()])).flux_time_monitors[0]
    m = scene_mod.from_simulation(
        _sim([_ftm(center=(*FT_CENTER, 0.3125))])).flux_time_monitors[0]
    assert m.plane_index == on.plane_index          # the same lower boundary
    np.testing.assert_allclose(m.frac, 0.5, atol=1e-9)


def test_converter_respects_start_interval():
    part = _ftm(start=5e-14, interval=2)
    sim = _sim([part])
    m = scene_mod.from_simulation(sim).flux_time_monitors[0]
    beg, end = (int(v) for v in sim.monitors[0].time_inds(np.asarray(sim.tmesh)))
    assert beg > 0 and (m.step_begin, m.interval) == (beg, 2)
    assert m.num_slots == (end - beg + 1) // 2


def test_refuses_box():
    """A 3-D monitor (a closed box) is not supported yet. ``interval_space≠1`` is blocked by
    tidy3d's own schema (the field is Literal[1]), so our gate never sees it for this type."""
    with pytest.raises(NotImplementedError, match="3-D FluxTimeMonitor"):
        scene_mod.from_simulation(_sim([
            _ftm(center=(0, 0, 0), size=(0.2, 0.2, 0.2))]))


def test_full_periodic_plane_wraps():
    """Filling a purely periodic axis: the samples wrap around (AndersonLocalization) and the
    integration measure comes back to the full period."""
    sc = scene_mod.from_simulation(_sim([
        _ftm(center=(0, 0, PLANE_Z), size=(td.inf, td.inf, 0))]))
    g = flux_mod.flux_time_geometry(sc.grid, sc.flux_time_monitors[0])
    for p in (0, 1):
        assert g["fold_lo"][p] and g["idx_lo"][p][0] == sc.grid.axes[p].n - 1
        # parity factor +1: wlo_c gives all four components the same weights
        ws = list(g["wlo_c"][p].values())
        assert all(np.array_equal(w, ws[0]) for w in ws[1:])
    # the area elements cover the full-period cross section 0.4×0.4 µm², no cell more or less
    np.testing.assert_allclose(float(np.abs(g["dS"]).sum()) / UM ** 2,
                               0.4 * 0.4, rtol=1e-12)


def test_periodic_wrap_translation_invariance():
    """Shift the whole field by one cell in a periodic domain: the wrapped flux has to follow the
    shift bitwise (discrete translation symmetry)."""
    rng = np.random.default_rng(3)
    edges = [np.linspace(0.0, 0.4 * UM, 17), np.linspace(0.0, 0.4 * UM, 17),
             np.linspace(0.0, 1.6 * UM, 65)]
    grid = Grid(Axis(edges=edges[0], boundary_lo="Periodic", boundary_hi="Periodic"),
                Axis(edges=edges[1], boundary_lo="Periodic", boundary_hi="Periodic"),
                Axis(edges=edges[2], boundary_lo="PML", boundary_hi="PML"))
    mon = FluxTimeMonitor(
        name="wrap", axis=2, plane_index=30, normal_dir=+1,
        step_begin=0, step_end=1, interval=1, num_slots=1,
        t_bounds=((0.0, 0.4 * UM), (0.0, 0.4 * UM)), frac=0.3)
    g = flux_mod.flux_time_geometry(grid, mon)
    shape = (16, 16, 64)
    f = {c: rng.standard_normal(shape) for c in ("Ex", "Ey", "Hx", "Hy")}
    s0 = flux_mod.flux_time_sample((f["Ex"], f["Ey"]), (f["Hx"], f["Hy"]), g)
    for axis_t in (0, 1):
        fr = {c: np.roll(v, 1, axis=axis_t) for c, v in f.items()}
        s1 = flux_mod.flux_time_sample((fr["Ex"], fr["Ey"]), (fr["Hx"], fr["Hy"]), g)
        # roll changes the floating-point summation order, so the floor is ~1e-12; 1e-10 is still
        # far below any physical effect
        np.testing.assert_allclose(s1, s0, rtol=1e-10)


def test_serialize_roundtrip(tmp_path):
    """Since v4 both FluxTimeMonitor and FieldTimeMonitor go into the npz (the solver job needs
    them)."""
    sim = _sim([
        _ftm(start=5e-14),
        td.FieldTimeMonitor(center=(0, 0, PLANE_Z), size=(0, 0, 0),
                            fields=["Ey"], name="pt"),
    ])
    sc = scene_mod.from_simulation(sim)
    p = serialize.save(sc, tmp_path / "scene.npz")
    back = serialize.load(p)
    assert back.flux_time_monitors == sc.flux_time_monitors
    assert back.field_time_monitors == sc.field_time_monitors


# ------------------------------------------------------------ oracle: against the Tidy3D client

def test_integration_matches_tidy3d_client():
    """Random colocated fields: the area elements and the integral agree with
    ``FieldTimeData.flux`` bitwise (≤1e-12).

    What this locks is the spatial convention of ``use_colocated_integration=True``: the sample
    set (tangential primal boundary points), the trapezoidal area elements truncated exactly at
    the monitor bounds, and the sign flip of H for a y normal. Non-uniform grid (auto grid plus a
    dielectric block), all three normal axes, monitor bounds off the grid points.
    """
    from tidy3d.components.data.data_array import ScalarFieldTimeDataArray

    rng = np.random.default_rng(0)
    sim = td.Simulation(
        size=(2.0, 2.4, 2.8),
        grid_spec=td.GridSpec.auto(min_steps_per_wvl=8, wavelength=1.0),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, 0), size=(0.8, 0.9, 1.0)),
            medium=td.Medium(permittivity=4.0))],
        sources=[td.PointDipole(center=(0, 0, 0), polarization="Ey",
                                source_time=td.GaussianPulse(freq0=3e14, fwidth=1e14))],
        run_time=1e-13,
        boundary_spec=td.BoundarySpec.all_sides(td.PML()),
    )
    grid = Grid(*(Axis(edges=np.asarray(sim.grid.boundaries.to_dict[d]) * UM,
                       boundary_lo="PML", boundary_hi="PML") for d in "xyz"))
    t = np.array([0.0, 1e-15, 2e-15])

    for ax in range(3):
        size = [1.1, 1.3, 1.5]
        size[ax] = 0.0
        center = [0.037, -0.051, 0.062]     # deliberately off the grid points
        t1, t2 = (i for i in range(3) if i != ax)
        comps = [f"E{'xyz'[t1]}", f"E{'xyz'[t2]}", f"H{'xyz'[t1]}", f"H{'xyz'[t2]}"]
        tmon = td.FieldTimeMonitor(center=center, size=size, name="ft",
                                   fields=comps, colocate=True, interval=1)
        gexp = sim.discretize_monitor(tmon)
        coords = {d: np.asarray(gexp.boundaries.to_dict[d])[:-1] for d in "xyz"}
        shape = tuple(coords[d].size for d in "xyz") + (t.size,)
        das = {c: ScalarFieldTimeDataArray(
            rng.standard_normal(shape), coords={**coords, "t": t}) for c in comps}
        data = td.FieldTimeData(monitor=tmon, grid_expanded=gexp, **das)
        ref_flux = np.asarray(data.flux)
        ref_dS = data._diff_area.to_numpy()

        # ---- our side: a FluxTimeMonitor with the same geometry + flux_time_geometry ----
        e = grid.axes[ax].edges
        x0 = center[ax] * UM
        km = int(np.searchsorted(e, x0, side="right")) - 1
        mon = FluxTimeMonitor(
            name="ft", axis=ax, plane_index=km, normal_dir=+1,
            step_begin=0, step_end=1, interval=1, num_slots=1,
            t_bounds=tuple((float(center[d] - size[d] / 2) * UM,
                            float(center[d] + size[d] / 2) * UM) for d in (t1, t2)),
            frac=(x0 - float(e[km])) / float(e[km + 1] - e[km]))
        g = flux_mod.flux_time_geometry(grid, mon)

        # sample sets line up: our edge run has to be a contiguous sub-run of the client's
        # tangential coordinates, and the area elements outside it are 0
        offs, ns = [], []
        for p, tt in enumerate((t1, t2)):
            ours = grid.axes[tt].edges[g["sl_edge"][p]]
            theirs = np.asarray(coords["xyz"[tt]]) * UM
            off = int(np.argmin(np.abs(theirs - ours[0])))
            np.testing.assert_allclose(theirs[off:off + ours.size], ours, rtol=1e-12)
            offs.append(off)
            ns.append(ours.size)
        mask = np.ones_like(ref_dS, dtype=bool)
        mask[offs[0]:offs[0] + ns[0], offs[1]:offs[1] + ns[1]] = False
        assert float(np.abs(ref_dS[mask]).max(initial=0.0)) == 0.0
        sub = (slice(offs[0], offs[0] + ns[0]), slice(offs[1], offs[1] + ns[1]))
        np.testing.assert_allclose(np.abs(g["dS"]) / UM ** 2, ref_dS[sub],
                                   rtol=1e-12, atol=0.0)

        vals = {c: np.squeeze(das[c].to_numpy(), axis=ax)[sub] for c in comps}
        ours_flux = np.array([
            flux_mod.integrate_colocated(vals[comps[0]][..., i], vals[comps[1]][..., i],
                                         vals[comps[2]][..., i], vals[comps[3]][..., i], g)
            for i in range(t.size)]) / UM ** 2
        np.testing.assert_allclose(ours_flux, ref_flux, rtol=1e-12)


def test_colocation_exact_on_affine_fields():
    """Colocation is exact on an affine field F=a+bx+cy+dz (linear interpolation has no error),
    which locks the indices and the weights."""
    rng = np.random.default_rng(1)
    edges = [np.concatenate([[0.0], np.cumsum(rng.uniform(0.5, 1.5, n))])
             for n in (12, 14, 16)]
    grid = Grid(*(Axis(edges=e, boundary_lo="PML", boundary_hi="PML")
                  for e in edges))
    for ax in range(3):
        t1, t2 = (i for i in range(3) if i != ax)
        e = grid.axes[ax].edges
        x0 = 0.6 * e[4] + 0.4 * e[5]
        km = int(np.searchsorted(e, x0, side="right")) - 1
        mon = FluxTimeMonitor(
            name="af", axis=ax, plane_index=km, normal_dir=+1,
            step_begin=0, step_end=1, interval=1, num_slots=1,
            t_bounds=tuple((float(grid.axes[t].edges[2]) + 0.3,
                            float(grid.axes[t].edges[-3]) - 0.3) for t in (t1, t2)),
            frac=(x0 - float(e[km])) / float(e[km + 1] - e[km]))
        g = flux_mod.flux_time_geometry(grid, mon)
        assert [pe["ke"] for pe in g["planes"]] == [km, km + 1]
        np.testing.assert_allclose(
            [pe["w"] for pe in g["planes"]], [1 - mon.frac, mon.frac], rtol=1e-12)
        eg1 = grid.axes[t1].edges[g["sl_edge"][0]]
        eg2 = grid.axes[t2].edges[g["sl_edge"][1]]
        for key, comp in (("Eu", f"E{'xyz'[t1]}"), ("Ev", f"E{'xyz'[t2]}"),
                          ("Hu", f"H{'xyz'[t1]}"), ("Hv", f"H{'xyz'[t2]}")):
            coef = rng.standard_normal(4)
            cs = _yee_coords(grid, comp)
            F = (coef[0] + coef[1] * cs[0][:, None, None]
                 + coef[2] * cs[1][None, :, None] + coef[3] * cs[2][None, None, :])
            for p, pe in enumerate(g["planes"]):
                pos = [None, None, None]
                pos[t1], pos[t2] = eg1[:, None], eg2[None, :]
                pos[ax] = float(e[pe["ke"]])
                expect = (coef[0] + coef[1] * pos[0] + coef[2] * pos[1]
                          + coef[3] * pos[2])
                got = flux_mod.colocated_component(F, g, key, p)
                np.testing.assert_allclose(got, expect, rtol=1e-12)


# ------------------------------------------------------------ GPU: three criteria share one run

@pytest.fixture(scope="module")
def runs():
    pytest.importorskip("cupy")
    from openem import solver
    from openem.device import Kernels

    # fld covers the whole plane plus a 3-cell neighborhood along the normal and feeds the host
    # reference for both ft (a boundary plane) and ft2 (an interpolated plane); ft and ft2
    # themselves are finite faces whose transverse bounds miss the grid points
    fld = td.FieldTimeMonitor(center=(0, 0, 0.3125), size=(td.inf, td.inf, 0),
                              fields=["Ex", "Ey", "Hx", "Hy"], name="fld")
    fx = td.FluxMonitor(center=(*FT_CENTER, PLANE_Z), size=(*FT_SIZE, 0),
                        freqs=list(FREQS), name="fx")
    ft2 = _ftm(center=(*FT_CENTER, 0.3125), name="ft2")    # off the boundary, frac=0.5
    k = Kernels()
    sc = scene_mod.from_simulation(_sim([_ftm(), ft2, fld, fx]))
    res = solver.run(sc, use_shutoff=False, kernels=k, verbose=False)
    # the same scene with the FluxTimeMonitor removed: the single-code-path criterion
    sc_no = replace(sc, flux_time_monitors=[])
    res_no = solver.run(sc_no, use_shutoff=False, kernels=k, verbose=False)
    return sc, res, res_no


def test_series_matches_fieldtime_flux(runs):
    """The GPU series ≡ samples from a FieldTimeMonitor on the same plane put through the same
    host colocation/integration, step by step.

    The H of a FieldTimeMonitor is already time-averaged onto the whole step, and on the host
    each slot goes through ``flux.colocated_component`` (with the normal index shifted into fld's
    box) plus ``integrate_colocated``. That is mathematically identical to the GPU's two-launch
    reduction, so the only differences are float32 rounding and summation order. Both a plane on
    a boundary (ft) and one off it (ft2) are covered.
    """
    sc, res, _ = runs
    fm = sc.field_time_monitors[0]
    buf = res.time_samples[fm.name]                      # (4, nslots, ni, nj, nk)
    c = {COMPONENTS[ci]: i for i, ci in enumerate(fm.comps)}
    oz = fm.origin[2]

    for ftm in sc.flux_time_monitors:
        g = flux_mod.flux_time_geometry(sc.grid, ftm)
        gl = dict(g)
        gl["planes"] = [{**pe, "ke": pe["ke"] - oz, "kh": pe["kh"] - oz}
                        for pe in g["planes"]]
        assert all(0 <= pe["kh"] and pe["ke"] < fm.box[2] for pe in gl["planes"])
        s = res.flux_time[ftm.name]
        n = min(s.size, buf.shape[1])
        assert n > 1000
        s_host = np.array([
            flux_mod.flux_time_sample(
                (buf[c["Ex"], i], buf[c["Ey"], i]),
                (buf[c["Hx"], i], buf[c["Hy"], i]), gl)
            for i in range(n)])
        err = float(np.max(np.abs(s[:n] - s_host)) / np.max(np.abs(s_host)))
        assert err < 1e-5, f"{ftm.name} (frac={ftm.frac}) differs from the FieldTime path by {err:.3e}"


def test_parseval_against_flux_monitor(runs):
    """Σ_m s(m) = 4·dt·∫₀^∞ P(f) df (derivation in model.FluxTimeMonitor).

    On the frequency side P(f) is computed straight from the phasor plane using the **same**
    colocation geometry as the time side (fx and ft share the plane and the bounds; the H of a
    phasor is already colocated onto the E plane along the normal, which is equivalent to our H
    interpolation on a frac=0 plane, the grid being uniform and x0 sitting exactly between two
    cell centers). The tolerance is not fitted: the deviation comes from the O((ωΔt)²)
    time-averaging error of H̄ and from the trapezoid rule, and with ω₀Δt ≈ 0.014 it is of order
    1e-3.
    """
    sc, res, _ = runs
    ftm = next(m for m in sc.flux_time_monitors if m.frac == 0.0)
    fx = sc.flux_monitors[0]
    assert fx.plane_index == ftm.plane_index
    g = flux_mod.flux_time_geometry(sc.grid, ftm)
    ph = res.phasors[fx.name].data                       # (4, nf, n1, n2)

    def tang(a, c2e):
        """On a phasor plane (already colocated along the normal), only the tangential
        center→edge step plus trimming of the edge samples."""
        if c2e == 0:
            b = a[:, g["sl_cell"][0], g["sl_edge"][1]]
            return g["wlo"][0][None] * b[:, :-1, :] + g["whi"][0][None] * b[:, 1:, :]
        b = a[:, g["sl_edge"][0], g["sl_cell"][1]]
        return g["wlo"][1][None] * b[:, :, :-1] + g["whi"][1][None] * b[:, :, 1:]

    eu, ev = tang(ph[flux_mod.EX], 0), tang(ph[flux_mod.EY], 1)
    hu, hv = tang(ph[flux_mod.HX], 1), tang(ph[flux_mod.HY], 0)
    p = 0.5 * np.real(np.einsum(
        "fuv,uv->f", eu * np.conj(hv) - ev * np.conj(hu), g["dS"]))
    lhs = float(np.sum(res.flux_time[ftm.name]))
    # trapezoid rule written out by hand: np.trapz is deprecated in numpy 2.x and the numpy
    # version differs between environments
    rhs = 4.0 * sc.dt * float(np.sum((p[1:] + p[:-1]) * np.diff(fx.freqs)) / 2.0)
    assert lhs > 0 and rhs > 0
    rel = abs(lhs - rhs) / rhs
    assert rel < 2e-2, f"Parseval deviation {rel:.3e} (lhs={lhs:.6e}, rhs={rhs:.6e})"


def test_scene_without_ftm_is_bitwise_unchanged(runs):
    """A scene without this monitor is bitwise unchanged: the frequency-domain phasors have to be
    exactly the same."""
    _, res, res_no = runs
    a = res.phasors["fx"].data
    b = res_no.phasors["fx"].data
    assert np.array_equal(a, b)
