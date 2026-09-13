# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Adjoint gradients with respect to geometry parameters (center, size, radius, vertices and so on
of Box, Cylinder, PolySlab and GeometryGroup) computed by perturbing the eps map, rather than by the
continuous surface integral in the tidy3d client.

## Why

The tidy3d client with ``local_gradient=True`` computes a geometry-parameter gradient from the
**continuous** shape-derivative formula: sample E and D on the surface of the body and integrate
delta_eps * E_par E_par - delta(1/eps) * D_perp D_perp over it. That assumes the field is well
resolved near the interface. When a feature is only 2 to 4 cells across it departs badly from the
**true slope of the discrete objective**, and not even consistently in one direction. Measured
against our own central differences:

| scene | resolution | parameter | continuous integral / central difference |
|---|---|---|---|
| circular hole in a 3D silicon slab | r = 2.25 cells | radius | 0.09 to 0.23 |
| the same | r = 4.5 cells | radius | 0.87 to 0.97 |
| first step of a photonic crystal case | r about 3.2 cells | radius | 1.63 |
| a 2D box face at a cell center / between cells / on a cell edge | - | gap | 2.49 / 1.04 / 0.71 |

The last row shows something else: after subpixel averaging, eps is **not linear in the face
position within a cell** (the eps against fill fraction that extras produces looks like an S
curve), so the response of the discrete objective to "move the face by 0.3 cells" depends on where
in the cell the face sits. The continuous formula gives the cell-averaged slope, independent of
position. Both are "right", but what an optimizer needs is the slope of the discrete objective.

## How

The discrete chain rule: dJ/dp = sum over nodes of Re[E_fwd . E_adj] * w_node * d(eps_node)/dp.
The first factor is exactly the volume VJP that tidy3d uses for CustomMedium pixels, already
verified against the reference. d(eps)/dp comes from central-differencing the subpixel eps map
(``sim.epsilon`` through tidy3d-extras, the same source as the eps fed to the solver) with respect
to the parameter, so no extra FDTD run is needed. Members of one GeometryGroup whose local boxes do
not overlap are perturbed together and their eps taken in one query: the 780 parameters of one
photonic crystal case need only 3 groups x 2 sides x 3 components = 18 eps queries.

## Usage

``apply`` installs the patches (autograd_hook.install calls it at the end) and ``revert`` removes
them. The mode comes from ``OPENEM_SHAPEGRAD``:

- ``mirror`` (default): does exactly one thing, namely record a zero gradient for any sub-geometry
  lying entirely on the mirror side of a symmetry plane. Both the solver (which folds) and
  ``sim.epsilon`` (through extras) see only the simulated half and mirror the other one, so those
  parameters never enter the physics at all. The tidy3d client, however, computes mirrored values
  for them from the unfolded field, which made one gradient norm larger than the reference by
  sqrt(2); the reference solver computes those components as zero. Every other geometry gradient
  still takes tidy3d's continuous surface integral.
- ``eps``: on top of ``mirror``, geometry-parameter gradients switch to the eps-map perturbation in
  this file, i.e. the discrete chain rule. The cost is a number of extra subpixel eps queries per
  body: about 18 per step in one case and about 170 in another.
- ``0``: install nothing.

``OPENEM_SHAPEGRAD_H`` overrides the eps perturbation step, in um; the default is a quarter of the
smallest cell.
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np

from openem import knobs

_STATE: dict = {"on": False, "orig_cmp": None, "orig_proc": None, "sim": None, "index": None}

#: How many extra cells the local box keeps beyond the geometry bounding box; it has to cover every
#: Yee node whose control volume touches the perturbed interface
MARGIN_CELLS = 1.5
COORD_TOL = 1e-9
_DEBUG = knobs.env("SHAPEGRAD_DEBUG") == "1"


#: The perturbation step as a fraction of the smallest cell. The subpixel eps against fill fraction
#: is a steep S curve with plateaus, so too small a step (1e-3 um) only touches the few nodes that
#: happen to flip. A quarter of a cell is the same scale as an optimizer step and as the central
#: differences used elsewhere, 0.004 to 0.006 um.
H_FRACTION = 0.25


def _h_default(dl_min: float) -> float:
    env = knobs.env("SHAPEGRAD_H")
    return float(env) if env else H_FRACTION * dl_min


# ---- Reading and writing geometry parameters, following tidy3d autograd's path convention ----
def _get(geom, sub):
    key, *rest = sub
    val = getattr(geom, key)
    if key == "geometries":
        return _get(val[rest[0]], tuple(rest[1:]))
    if rest:
        return val[rest[0]]
    return val


def _with(geom, sub, value):
    """Return a copy of the geometry with the parameter at path ``sub`` replaced by ``value``."""
    key, *rest = sub
    if key == "geometries":
        k = rest[0]
        geos = list(geom.geometries)
        geos[k] = _with(geos[k], tuple(rest[1:]), value)
        return geom.updated_copy(geometries=geos)
    cur = getattr(geom, key)
    if rest:                               # center[i] / size[i] / slab_bounds[i]
        lst = list(cur)
        lst[rest[0]] = value
        return geom.updated_copy(**{key: tuple(lst)})
    return geom.updated_copy(**{key: value})


def _sub_geometry(geom, sub):
    """The (sub-)geometry being perturbed, used to compute the local box."""
    key, *rest = sub
    if key == "geometries":
        return _sub_geometry(geom.geometries[rest[0]], tuple(rest[1:]))
    return geom


def _scalar_tasks(geom, sub, h):
    """Expand one path into scalar perturbation tasks, [(elem, geom_plus, geom_minus)].

    ``elem`` is the index within an array parameter such as vertices, and None for a scalar
    parameter.
    """
    val = _get(geom, sub)
    arr = np.asarray(val, dtype=float)
    if arr.ndim == 0:
        v = float(arr)
        return [(None, _with(geom, sub, v + h), _with(geom, sub, v - h))]
    tasks = []
    for idx in np.ndindex(arr.shape):
        up, dn = arr.copy(), arr.copy()
        up[idx] += h
        dn[idx] -= h
        tasks.append((idx, _with(geom, sub, up), _with(geom, sub, dn)))
    return tasks


def _bounds_union(a, b):
    lo = np.minimum(np.asarray(a[0], float), np.asarray(b[0], float))
    hi = np.maximum(np.asarray(a[1], float), np.asarray(b[1], float))
    return lo, hi


def _overlap(b1, b2) -> bool:
    return bool(np.all(b1[0] <= b2[1]) and np.all(b2[0] <= b1[1]))


def _on_mirror_side(bounds_hi, center, sym_axes) -> bool:
    """Whether the geometry lies entirely on the mirror side of a symmetry plane, i.e. its bounding
    box upper edge is below the symmetry center on any symmetry axis.
    """
    return any(float(bounds_hi[d]) < center[d] - COORD_TOL for d in sym_axes)


class _MonCtx(NamedTuple):
    """Sampling context on the monitor side: per-component node coordinates ``coords`` and volume
    weights ``weights``, the adjoint field map ``E_der``, the monitor frequencies ``freqs``, and the
    frequencies ``eps_freqs`` at which eps is taken (with ``flat=True`` the medium is
    frequency-independent and only the first frequency is used).

    The per-task summation (:func:`_task_contribution`) only reads it, so it travels as one group.
    """

    coords: dict
    weights: dict
    E_der: dict
    freqs: np.ndarray
    eps_freqs: np.ndarray
    flat: bool


class _BoxClamp(NamedTuple):
    """How far the local box is grown and where it is clamped: ``margin`` is the outward growth,
    ``mon_lo``/``mon_hi`` the extent of the monitor box, and ``sym_axes``/``center`` the symmetry
    planes, since only the simulated half is counted.
    """

    margin: np.ndarray
    mon_lo: np.ndarray
    mon_hi: np.ndarray
    sym_axes: list
    center: np.ndarray


# ---- Core ----
def _monitor_geometry(info):
    """Coordinates and volume weights of the monitor box:
    ``(coords, weights, mon_lo, mon_hi, dl_min, dl_max)``.
    """
    from tidy3d.components.autograd.derivative_utils import compute_spatial_weights

    E_der = info.E_der_map                       # {Ex,Ey,Ez: (x,y,z,f)}
    coords = {c: [np.asarray(E_der[c].coords[a].values, float) for a in "xyz"] for c in ("Ex", "Ey", "Ez")}
    dl_min = min(float(np.min(np.diff(coords[c][d]))) for c in coords for d in range(3) if len(coords[c][d]) > 1)
    weights = {}
    for c in ("Ex", "Ey", "Ez"):
        w = compute_spatial_weights(E_der[c])
        w3 = np.ones(tuple(len(x) for x in coords[c]))
        if np.ndim(w.values) > 0:
            wv = np.asarray(w.values, float)
            shp = [1, 1, 1]
            for d, a in enumerate("xyz"):
                if a in w.dims:
                    shp[d] = len(coords[c][d])
            w3 = w3 * wv.reshape(shp)
        weights[c] = w3
    mon_lo = np.array([min(float(coords[c][d][0]) for c in coords) for d in range(3)])
    mon_hi = np.array([max(float(coords[c][d][-1]) for c in coords) for d in range(3)])
    dl_max = np.array([max(float(np.max(np.diff(coords[c][d]))) if len(coords[c][d]) > 1 else 0.0
                       for c in coords) for d in range(3)])
    return coords, weights, mon_lo, mon_hi, dl_min, dl_max


def _expand_tasks(geom, info, h, clamp: _BoxClamp):
    """Expand info.paths into scalar perturbation tasks
    ``[(path, elem, gp, gm, box_lo, box_hi)]`` and compute their local boxes.

    A sub-geometry lying entirely on the mirror side becomes no task at all and records a zero
    gradient (the returned ``zero_out``).

    On symmetry: sim.epsilon through extras, like the solver, computes only the simulated half and
    mirrors the other. Measured with symmetry=(0,-1,0), perturbing a hole at y<0 changed eps not at
    all, while perturbing one at y>0 changed both sides. So a sub-geometry entirely on the mirror
    side has zero gradient, since it never enters the physics, and every other task sums only over
    nodes on the simulated side. Each body therefore reports a "single hole" value, matching both
    the tidy3d client, which computes per hole after unfolding the field, and the reference solver.
    """
    sym_axes, center = clamp.sym_axes, clamp.center
    zero_out: dict = {}
    tasks = []      # (path, elem, gp, gm, box_lo, box_hi)
    for path in info.paths:
        sub = tuple(path[1:])
        for elem, gp, gm in _scalar_tasks(geom, sub, h):
            sg_p, sg_m = _sub_geometry(gp, sub), _sub_geometry(gm, sub)
            lo, hi = _bounds_union(sg_p.bounds, sg_m.bounds)
            if sym_axes and _on_mirror_side(hi, center, sym_axes):
                if elem is None:
                    zero_out[path] = 0.0
                else:
                    zero_out.setdefault(path, np.zeros(np.asarray(_get(geom, sub), float).shape))
                continue
            lo = np.maximum(lo - clamp.margin, clamp.mon_lo)
            hi = np.minimum(hi + clamp.margin, clamp.mon_hi)
            for d in sym_axes:
                lo[d] = max(lo[d], center[d] - COORD_TOL)      # count only the simulated half
            # (There used to be another clamp before subtracting margin. Clamp, subtract margin,
            #  clamp is the same as subtract margin, clamp: worked through by hand for both lo<c and
            #  lo>=c, and max does not round.)
            tasks.append((path, elem, gp, gm, lo, hi))
    return tasks, zero_out


def _batch_nonoverlapping(tasks) -> list[list[int]]:
    """Group tasks whose local boxes do not overlap into one batch; a batch costs only
    2 sides x number of components eps queries.
    """
    batches: list[list[int]] = []
    for ti, t in enumerate(tasks):
        for b in batches:
            if not any(_overlap((t[4], t[5]), (tasks[j][4], tasks[j][5])) for j in b):
                b.append(ti)
                break
        else:
            batches.append([ti])
    return batches


def _delta_eps(sims, box, comps, eps_freqs, h) -> dict:
    """``(comp, f) -> (array of delta_eps/2h, coordinates on all three axes)``, the central
    difference of the subpixel eps maps of the plus and minus perturbations.

    ``sims`` is the pair of perturbed Simulations, ``(sim_p, sim_m)``.
    """
    sim_p, sim_m = sims
    deps = {}
    for c in comps:
        for f in eps_freqs:
            ep = sim_p.epsilon(box, coord_key=c, freq=float(f))
            em = sim_m.epsilon(box, coord_key=c, freq=float(f))
            d = (np.asarray(ep.values, complex) - np.asarray(em.values, complex)) / (2 * h)
            deps[(c, float(f))] = (d, [np.asarray(ep.coords[a].values, float) for a in "xyz"])
    return deps


def _task_contribution(task, deps, mon: _MonCtx) -> float:
    """dJ/dp for one scalar task: sum over nodes of Re[E_fwd . E_adj] * w * d(eps)/dp, restricted to
    that task's own local box, with the eps nodes matched to the monitor nodes.

    Accumulation order is components outer, frequencies inner.
    """
    coords, weights, E_der = mon.coords, mon.weights, mon.E_der
    path, elem, _, _, tlo, thi = task
    val = 0.0
    for c in ("Ex", "Ey", "Ez"):
        for fi, f in enumerate(mon.freqs):
            d, dc = deps[(c, float(mon.eps_freqs[0] if mon.flat else f))]
            sel = []
            ok = True
            for a in range(3):
                if len(coords[c][a]) == 1 and len(dc[a]) == 1:
                    # A flat 2D axis: the single coordinate the monitor gives and the one
                    # sim.epsilon gives may differ (Ez at a cell center against one snapped to 0),
                    # so pair them directly
                    sel.append((np.array([0]), np.array([0])))
                    continue
                m = (dc[a] >= tlo[a] - COORD_TOL) & (dc[a] <= thi[a] + COORD_TOL)
                if not np.any(m):
                    ok = False
                    break
                j = np.searchsorted(coords[c][a], dc[a][m])
                j = np.clip(j, 0, len(coords[c][a]) - 1)
                jm = np.clip(j - 1, 0, len(coords[c][a]) - 1)
                near = np.where(np.abs(coords[c][a][jm] - dc[a][m]) < np.abs(coords[c][a][j] - dc[a][m]), jm, j)
                good = np.abs(coords[c][a][near] - dc[a][m]) <= 1e-6
                sel.append((np.flatnonzero(m)[good], near[good]))
            if not ok:
                continue
            dsub = d[np.ix_(sel[0][0], sel[1][0], sel[2][0])]
            esub = np.asarray(E_der[c].values)[np.ix_(sel[0][1], sel[1][1], sel[2][1])][..., fi]
            wsub = weights[c][np.ix_(sel[0][1], sel[1][1], sel[2][1])]
            contrib = float(np.real(np.sum(dsub * esub * wsub)))
            if _DEBUG:
                print(f"[shapegrad] path={path} elem={elem} {c} f={f:.4e} nodes={dsub.shape} "
                      f"max|dε|={np.abs(dsub).max():.3e} max|E_der|={np.abs(esub).max():.3e} "
                      f"w[{wsub.min():.2e},{wsub.max():.2e}] contrib={contrib:+.4e}", flush=True)
            val += contrib
    if _DEBUG:
        print(f"[shapegrad] task path={path} elem={elem} -> {val:+.6e}", flush=True)
    return val


def _numerical_geometry_vjp(structure, info) -> dict:
    import tidy3d as td

    sim = _STATE["sim"]
    idx = _STATE["index"]
    if sim is None or idx is None:
        raise RuntimeError("shapegrad: the original Simulation is unavailable, so the backward patch "
                           "did not take effect")
    geom = structure.geometry
    E_der = info.E_der_map                       # {Ex,Ey,Ez: (x,y,z,f)}
    freqs = np.atleast_1d(np.asarray(getattr(info.frequencies, "values", info.frequencies), float))
    coords, weights, mon_lo, mon_hi, dl_min, dl_max = _monitor_geometry(info)
    h = _h_default(dl_min)
    if _DEBUG:
        print(f"[shapegrad] dl_min={dl_min:.4g} h={h:.4g}", flush=True)
    margin = MARGIN_CELLS * dl_max + h
    # When the medium is frequency-independent, eps is taken at one frequency only. The test is
    # shared with nb.backend through scene/media.is_freq_independent. This used to recognize only
    # Medium, so PEC, CustomMedium and non-dispersive anisotropic media all went through multi-
    # frequency sampling: the same eps values, merely slower.
    from openem.scene import media as _media
    flat = _media.is_freq_independent(sim)
    eps_freqs = freqs[:1] if flat else freqs

    sym = tuple(int(v) for v in getattr(sim, "symmetry", (0, 0, 0)))
    center = np.asarray(sim.center, float)
    sym_axes = [d for d in range(3) if sym[d] != 0]

    # 1) expand into scalar tasks and compute their local boxes (mirror-side ones record 0 directly);
    # 2) batch together the tasks whose local boxes do not overlap
    mon = _MonCtx(coords, weights, E_der, freqs, eps_freqs, flat)
    clamp = _BoxClamp(margin, mon_lo, mon_hi, sym_axes, center)
    tasks, out = _expand_tasks(geom, info, h, clamp)
    batches = _batch_nonoverlapping(tasks)

    def _apply_all(sign, members):
        g = geom
        for ti in members:
            path = tasks[ti][0]
            sub = tuple(path[1:])
            src = tasks[ti][2] if sign > 0 else tasks[ti][3]
            g = _with(g, sub, _get(src, sub))
        structures = list(sim.structures)
        structures[idx] = structure.updated_copy(geometry=g)
        return sim.updated_copy(structures=structures, grid_spec=td.GridSpec.from_grid(sim.grid))

    # The eps queries must use the same subpixel convention as the solver; the global extras switch
    # may have been turned off by the mode solve or another path (see media.set_local_subpixel)
    from openem.scene import subpixel as _subpixel
    from tidy3d import config as _config
    saved_flag = _config.simulation.use_local_subpixel
    with _subpixel.local_subpixel(bool(sim.subpixel)):
        if _DEBUG:
            print(f"[shapegrad] sim.subpixel={bool(sim.subpixel)} use_local_subpixel {saved_flag} -> "
                  f"{_config.simulation.use_local_subpixel}", flush=True)
            # Whether the eps_data actually fed to tidy3d in the job is subpixel averaged: count the
            # mixed cells
            try:
                for c in ("eps_xx", "eps_yy", "eps_zz"):
                    ev = np.real(np.asarray(info.eps_data[c].values))
                    u = np.unique(np.round(ev, 4))
                    mixed = int(np.sum((ev > u.min() + 1e-3) & (ev < u.max() - 1e-3)))
                    print(f"[shapegrad] eps_data {c}: n={ev.size} mixed={mixed} distinct={len(u)} min={u.min():.3f} max={u.max():.3f}", flush=True)
            except Exception as _e:
                print("[shapegrad] eps_data check failed", repr(_e)[:80], flush=True)
        for b in batches:
            lo = np.min([tasks[ti][4] for ti in b], axis=0)
            hi = np.max([tasks[ti][5] for ti in b], axis=0)
            box = td.Box.from_bounds(rmin=tuple(lo), rmax=tuple(hi))
            sim_p, sim_m = _apply_all(+1, b), _apply_all(-1, b)
            if _DEBUG:
                print(f"[shapegrad] batch of {len(b)} tasks, box {np.round(lo, 4)}..{np.round(hi, 4)}, "
                      f"freqs={eps_freqs.tolist()}", flush=True)
            deps = _delta_eps((sim_p, sim_m), box, ("Ex", "Ey", "Ez"), eps_freqs, h)
            for ti in b:
                path, elem = tasks[ti][0], tasks[ti][1]
                val = _task_contribution(tasks[ti], deps, mon)
                if elem is None:
                    out[path] = val
                else:
                    arr = out.setdefault(path, np.zeros(np.asarray(_get(geom, tuple(path[1:])), float).shape))
                    arr[elem] += val
    return out


# ---- Patches ----
def _mode() -> str:
    return knobs.env("SHAPEGRAD").strip().lower()


def _mirror_side_paths(structure, info) -> list:
    """The geometry paths in info.paths that lie entirely on the mirror side of a symmetry plane.

    Those parameters never enter the physics, so their gradient should be 0.
    """
    sim = _STATE["sim"]
    if sim is None:
        return []
    sym = tuple(int(v) for v in getattr(sim, "symmetry", (0, 0, 0)))
    sym_axes = [d for d in range(3) if sym[d] != 0]
    if not sym_axes:
        return []
    center = np.asarray(sim.center, float)
    out = []
    for path in info.paths:
        if path[0] != "geometry":
            continue
        sg = _sub_geometry(structure.geometry, tuple(path[1:]))
        if _on_mirror_side(np.asarray(sg.bounds[1], float), center, sym_axes):
            out.append(path)
    return out


def _zero_like(structure, path):
    val = np.asarray(_get(structure.geometry, tuple(path[1:])), float)
    return 0.0 if val.ndim == 0 else np.zeros(val.shape)


def apply() -> bool:
    mode = _mode()
    if _STATE["on"] or mode in ("0", "off", "false", ""):
        return False
    if mode not in ("mirror", "eps"):
        raise ValueError(f"OPENEM_SHAPEGRAD={mode!r}: only mirror, eps and 0 are recognized")
    from tidy3d.components.structure import Structure
    from tidy3d.web.api.autograd import backward as _bw

    orig_cmp = Structure._compute_derivatives
    orig_proc = _bw._process_structure_gradients

    def _proc(sim_data_adj, sim_data_orig, sim_data_fwd, structure_index, *a, **kw):
        prev = (_STATE["sim"], _STATE["index"])
        _STATE["sim"], _STATE["index"] = sim_data_orig.simulation, int(structure_index)
        try:
            return orig_proc(sim_data_adj, sim_data_orig, sim_data_fwd, structure_index, *a, **kw)
        finally:
            _STATE["sim"], _STATE["index"] = prev

    def _cmp(self, derivative_info, vjp_fns=None):
        mirror = set(_mirror_side_paths(self, derivative_info))
        out = {p: _zero_like(self, p) for p in mirror}
        geo = [p for p in derivative_info.paths if p[0] == "geometry" and p not in mirror]
        rest = [p for p in derivative_info.paths if p[0] != "geometry"]
        if mode == "eps":
            if rest:
                out.update(orig_cmp(self, derivative_info.updated_copy(paths=rest, deep=False), vjp_fns=vjp_fns))
            if geo:
                out.update(_numerical_geometry_vjp(self, derivative_info.updated_copy(paths=geo, deep=False)))
        else:
            keep = rest + geo
            if keep:
                out.update(orig_cmp(self, derivative_info.updated_copy(paths=keep, deep=False), vjp_fns=vjp_fns))
        return out

    _STATE.update(orig_cmp=orig_cmp, orig_proc=orig_proc, on=True, mode=mode)
    _bw._process_structure_gradients = _proc
    Structure._compute_derivatives = _cmp
    return True


def revert() -> None:
    if not _STATE["on"]:
        return
    from tidy3d.components.structure import Structure
    from tidy3d.web.api.autograd import backward as _bw
    Structure._compute_derivatives = _STATE["orig_cmp"]
    _bw._process_structure_gradients = _STATE["orig_proc"]
    _STATE["on"] = False
