# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Wire tidy3d's autograd inverse design through to OpenEM, with **no change to the notebook**.

tidy3d's ``local_gradient`` mode adds its own monitors to the simulation (a ``FieldMonitor`` plus a
``PermittivityMonitor`` over the design region), builds its own adjoint sources, and assembles the
gradient locally. **Execution leaves through only two entry points**:

    forward   ``web/api/autograd/engine.py::_run_tidy3d(sim, task_name, **kw)``
              -> ``(SimulationData, task_id)``
    adjoint   ``_run_async_tidy3d(sims: dict, **kw)``
           → ``(BatchData, {task_name: task_id})``

Replace those two and both the forward and the adjoint solve of an inverse design run here, while
the gradient assembly stays tidy3d's, so its VJPs need no rewriting.

**Two things to watch out for**, both learned the hard way:

1. ``config.adjoint.local_gradient`` defaults to **False**, in which case the gradient is computed
   on the server, never passes through these two entry points, and **really does go to the cloud
   and cost money**. So :func:`install` sets it to True explicitly and replaces ``Job.run`` and
   ``Batch.run`` with gates that raise: a cloud path that slipped through should stop rather than
   bill silently.
2. The call site is ``hooks._run_tidy3d``, which captured its own reference at import time, so
   patching ``engine._run_tidy3d`` alone does nothing. All three modules have to be patched.

The adjoint source is a ``td.CustomCurrentSource`` (point-by-point complex-amplitude current),
supported by ``scene/sources.py::custom_current``.
"""
from __future__ import annotations

import contextlib
import re
import time
from types import SimpleNamespace

import numpy as np

from openem import scene as scene_mod
from openem import solver
from openem.device import Kernels

#: One record per solve, ``(label, cells, steps, seconds)``, for diagnostics.
CALLS: list[tuple] = []

_STATE: dict = {"kernels": None, "steps": None, "max_cells": None}

#: Steps each forward solve actually ran, ``{task_name: steps}``; the adjoint copies them
_FWD_STEPS: dict[str, int] = {}
_ADJ_RE = re.compile(r"^(.*)_adjoint_\d+$")


def _result_to_ours(res) -> dict:
    """``solver.Result`` to the dict that ``nb.backend.build`` consumes, matching
    nb.solve_worker.save_result.
    """
    arrs: dict = {}
    for n, p in res.phasors.items():
        arrs[f"mon:{n}"] = p.data
        arrs[f"monf:{n}"] = np.asarray(p.freqs)
    for attr, tag in (("field_phasors", "fld"), ("time_samples", "tim"),
                      ("flux_time", "ftm")):
        for n, a in (getattr(res, attr, {}) or {}).items():
            arrs[f"{tag}:{n}"] = a
    import json as _json
    from openem.results import result_meta
    arrs["__meta"] = np.array(_json.dumps(result_meta(res)))   # same as in ours.npz: a 0-d string array, read with .item
    return arrs


def _dispatch_non_fdtd(simulation, tag):
    """Dispatch containers and mode solvers; every entry point (Job.run, Batch.run, web.run) ends up
    here. A plain FDTD simulation returns None and is handed back to :func:`_solve`.

    - list/tuple/dict: solved one at a time (some notebooks hand a whole sequence to web.run);
    - ModeSolver / ModeSimulation: not FDTD, so handed to tidy3d's local mode solver (one notebook
      submits mode solves in bulk through Batch). OpenEM does not do mode solving itself; this step
      only connects a notebook's incidental mode-solve need to Tidy3D's own solver.
    """
    from openem.nb import backend as nb

    if isinstance(simulation, dict):
        return _BatchShim({k: _solve(s, f"{tag}:{k}") for k, s in simulation.items()})
    if isinstance(simulation, (list, tuple)):
        return [_solve(s, f"{tag}[{i}]") for i, s in enumerate(simulation)]
    if type(simulation).__name__.endswith("EMESimulation"):
        # EME (eigenmode expansion) is a different solver; OpenEM is pure FDTD. Without this guard
        # it walks all the way into from_simulation and raises "'EMESimulation' object has no
        # attribute 'dt'", which says nothing about what is actually missing.
        raise NotImplementedError(
            "the EME (eigenmode expansion) solver is out of scope for OpenEM, which is pure FDTD. "
            "A notebook containing EME can only have its FDTD part compared here.")
    if type(simulation).__name__.endswith("ModeSimulation"):
        # This branch used to require hasattr(sim, "data"); a ModeSimulation has no such attribute
        # and slipped through, raising "'ModeSimulation' object has no attribute
        # 'frequency_range'". It has its own readout function, so use that.
        return nb._mode_simulation_data(simulation, False)
    if type(simulation).__name__.endswith("ModeSolver"):
        return nb._mode_solver_data(simulation, False)
    return None


def _build_scene_with_subpixel_fallback(simulation):
    """``(simulation, scene)``: build the scene, falling back to subpixel=False when the subpixel
    interface cells cannot be handled.
    """
    # Strip autograd tracer values (ArrayBox). On the normal path tidy3d has already stripped them
    # before _run_tidy3d, but the higher-level wrappers (tidy3d.plugins.autograd.optimize, or a
    # direct web.run) hand a traced Simulation straight through: in one case
    # GaussianBeam.angle_theta was still an ArrayBox and float(...) inside from_simulation raised
    # TypeError. to_static returns a sim with no tracers unchanged (it short-circuits on
    # _has_tracers), so the normal path is unaffected.
    to_static = getattr(simulation, "to_static", None)
    if callable(to_static):
        simulation = to_static()
    try:
        sc = scene_mod.from_simulation(simulation)
    except NotImplementedError as _exc:
        # The same backstop as nb.backend._openem_run_inner: when our builder cannot handle the
        # interface cells that subpixel averaging produced, rebuild with subpixel=False. That is
        # staircasing, not an equivalent, and how much it differs is answered by the comparison
        # against the reference notebook. Some scenes take the Batch route and never reach the copy
        # of this in openem_run.
        if "subpixel" not in str(_exc):
            raise
        print(f"[openem] our builder does not support subpixel with this medium combination "
              f"({str(_exc)[:48]}...); rebuilding with subpixel=False")
        simulation = simulation.updated_copy(subpixel=False)
        sc = scene_mod.from_simulation(simulation)
    return simulation, sc


def _ladder_log(tag, last, suffix, stair):
    """Log text for the divergence ladder on this path; the copy in openem_run appends the task name
    (see nb._ladder_log).
    """
    print(f"[openem] {tag}: the solve diverged ({str(last).splitlines()[0][:60]}...); retrying with "
          f"fallback {suffix}"
          f"{' (dispersive interface cells staircased)' if stair else ''}", flush=True)


@scene_mod.media.cloud_emulation()
def _solve(simulation, tag, steps=None):
    from openem.nb import backend as nb

    other = _dispatch_non_fdtd(simulation, tag)
    if other is not None:
        return other

    t0 = time.time()
    simulation, sc = _build_scene_with_subpixel_fallback(simulation)
    cells = int(np.prod(sc.shape))
    cap = _STATE["max_cells"]
    if cap is not None and cells > cap:
        raise RuntimeError(
            f"this simulation has {cells:,} cells, above the install(max_cells={cap:,}) ceiling; "
            "raise the ceiling, or lower the resolution in the notebook")
    if _STATE["kernels"] is None:
        _STATE["kernels"] = Kernels()
    n = _STATE["steps"] if steps is None else steps
    res = solver.run(sc, num_steps=n, use_shutoff=n is None,
                     kernels=_STATE["kernels"], verbose=False)
    def _rerun(sim2, _suffix):
        """One rung of the ladder: rebuild the scene, rerun, reassemble. sc and res are kept for the
        CALLS bookkeeping below.
        """
        nonlocal sc, res
        sc = scene_mod.from_simulation(sim2)
        res = solver.run(sc, num_steps=n, use_shutoff=n is None,
                         kernels=_STATE["kernels"], verbose=False)
        return nb.build(sim2, sc, _result_to_ours(res))

    try:
        sd = nb.build(simulation, sc, _result_to_ours(res))
    except RuntimeError as _e:
        # The same ladder as openem_run (nb.run_with_divergence_ladder / DIVERGENCE_LADDER):
        # Courant 0.8, then 0.6, then staircasing plus 0.6. If every rung diverges, failed_data
        # marks it failed by zeroing the monitor data, and the optimizer discards that design.
        sd = nb.run_with_divergence_ladder(simulation, _rerun, tag, first_exc=_e, log=_ladder_log)
    CALLS.append((tag, cells, int(res.steps_run), round(time.time() - t0, 2)))
    return sd


def _run_one(simulation, task_name, **kw):
    sd = _solve(simulation, str(task_name))
    _FWD_STEPS[str(task_name)] = CALLS[-1][2]      # record how many steps the forward solve ran
    return sd, "openem-local"


class _BatchShim(dict):
    """A minimal stand-in for ``BatchData``: downstream code only looks up a ``SimulationData`` by
    task_name.
    """

    @property
    def task_ids(self) -> dict:
        return {k: "openem-local" for k in self}

    @property
    def task_paths(self) -> dict:
        # Some notebooks display this as a log or archive path; a local solve writes nothing to disk
        return {k: "openem-local" for k in self}


def _adjoint_steps(task_name: str) -> int | None:
    """Map an adjoint task name, ``<forward name>_adjoint_<n>``, to the step count of that forward
    solve.
    """
    m = _ADJ_RE.match(str(task_name))
    return _FWD_STEPS.get(m.group(1)) if m else None


def _run_many(simulations, **kw):
    # **An adjoint solve must run the same number of steps as its forward solve.** The adjoint is
    # the transpose of the same discrete operator, and letting each terminate early on its own makes
    # the result not the gradient of that objective. Measured on one inverse-design case, with the
    # forward at 5600 steps and the adjoint at only 4000, the objective got worse the more it was
    # optimized and the gradient norm collapsed quickly.
    out = _BatchShim({k: _solve(s, f"batch:{k}", steps=_adjoint_steps(k))
                      for k, s in simulations.items()})
    return out, dict(out.task_ids)


def _has_tracers(sim) -> bool:
    """Whether this Simulation carries any autograd-traced field, i.e. an ArrayBox in a structure or
    a source.
    """
    try:
        return bool(sim._strip_traced_fields())
    except Exception:
        return False


@scene_mod.media.cloud_emulation()
def _run_plain(simulation, task_name="plain", **kw):
    """``web.run`` lands on OpenEM, along one of two routes.

    1. **A single Simulation carrying tracers**, i.e. the notebook calls ``web.run(sim)``
       **directly** under ``ag.value_and_grad``, ``ag.grad`` or ``plugins.autograd.optimize``.
       Here it **must be handed back to tidy3d's native autograd ``run``**, which is wrapped in
       ``@primitive``: forward, it extracts ``.data`` of the ``SimulationData`` into a traced field
       map; backward, it builds the adjoint and computes the gradient itself. The native ``run``
       still calls ``_run_tidy3d`` internally, which :func:`install` has replaced with a local
       solve, so nothing goes to the cloud or over the network.

       This used to send **every** ``web.run`` straight to ``_solve`` and return a plain
       ``SimulationData``, bypassing that primitive so the output carried no tracer. The objective
       then read **plain numpy** from ``sim_data[...].amps.data`` or ``mon.flux.data`` and the
       gradient broke right there; in one case ``aux._value`` raised ``AttributeError`` precisely
       because aux had received a ``numpy.float64``.

    2. **No tracer, or a dict or list**, i.e. a normal forward run or a multi-objective batch:
       solve locally one at a time, returning a BatchData stand-in for a dict and a list in order
       for a list.
    """
    # A ModeSolver is not an FDTD simulation; several notebooks submit a remote mode solve with
    # web.run(simulation=mode_solver). Send it to tidy3d's local mode solver, the same branch as in
    # nb.backend.run. Otherwise it is treated as a source-free Simulation and the scene build dies
    # on "no reference frequency available".
    if simulation is None:
        simulation = kw.get("simulation")
    # Dispatch by type name only, no longer by hasattr(.., "data"): ModeSolver.data is a property
    # that triggers a local solve, and when that raised AttributeError internally inside a job,
    # hasattr went False and the whole dispatch stopped working.
    if type(simulation).__name__.endswith("ModeSolver"):
        from openem.nb import backend as _nb
        return _nb._mode_solver_data(simulation, kw.get("verbose", True))
    if type(simulation).__name__.endswith("ModeSimulation"):
        from openem.nb import backend as _nb
        return _nb._mode_simulation_data(simulation, kw.get("verbose", True))
    ag_run = _STATE.get("autograd_run")
    if (ag_run is not None
            and not isinstance(simulation, (dict, list, tuple))
            and _has_tracers(simulation)):
        return ag_run(simulation, task_name=task_name,
                      verbose=kw.get("verbose", False))
    # The same holds for a dict or list: with tracers present it must be handed back to the native
    # run_async, or the gradient of a batched objective breaks here. One multi-objective case had a
    # batched gradient of identically zero while the same objective run one at a time was non-zero.
    # The _run_async_tidy3d that the native run_async calls has already been replaced by install
    # with a local solve, so nothing goes to the cloud.
    ag_run_async = _STATE.get("autograd_run_async")
    if (ag_run_async is not None
            and isinstance(simulation, (dict, list, tuple))
            and any(_has_tracers(s) for s in (
                simulation.values() if isinstance(simulation, dict)
                else simulation))):
        return ag_run_async(simulation, verbose=kw.get("verbose", False))
    if hasattr(simulation, "sim_dict"):
        # The Smatrix plugin's ComponentModeler: hand it back to tidy3d's _run_local, which uses
        # web.Batch internally, and Batch.run has already been replaced with a local solve.
        from tidy3d.plugins.smatrix.run import _run_local
        return _run_local(simulation, verbose=False)
    # A single Simulation and a dict or list container both go to _solve, which dispatches containers
    return _solve(simulation, f"plain:{task_name}")



def _noop(self, *a, **kw):
    return None


def _batch_get_info(self, *a, **kw):
    sims = getattr(self, "simulations", None) or {}
    return {nm: SimpleNamespace(taskName=nm, taskId=f"openem-local-{nm}", status="success")
            for nm in sims}


def _job_get_info(self, *a, **kw):
    return SimpleNamespace(taskName=getattr(self, "task_name", "openem"),
                           taskId="openem-local", status="success")


def _job_run(self, path=None, **kw):
    """``web.Job(...).run`` -> a local solve; several notebooks use this entry point."""
    return _solve(self.simulation, "job:" + str(self.task_name))


_BATCH_CACHE = {}


def _batch_run(self, path_dir=None, **kw):
    """``web.Batch(...).run`` -> solve locally one at a time, returning a BatchData stand-in.

    Results are cached under the sorted tuple of simulation names so a later ``Batch.load`` hits
    them. Notebooks often run and then load again from disk, but a local solve writes nothing to
    disk and a load would go to the cloud's get_info and get a 404.
    """
    res = _solve(dict(self.simulations), "batch")
    shim = res if isinstance(res, _BatchShim) else _BatchShim(dict(res))
    _BATCH_CACHE[tuple(sorted(self.simulations))] = shim
    _BATCH_CACHE["__last__"] = shim
    return shim


def _batch_load(self, path_dir=None, **kw):
    """``web.Batch(...).load``, including after a rebuild, hits the local results of the previous run."""
    key = tuple(sorted(getattr(self, "simulations", {}) or {}))
    return _BATCH_CACHE.get(key) or _BATCH_CACHE.get("__last__") or _BatchShim({})


def _patch_autograd_entrypoints() -> None:
    """Point tidy3d autograd's two execution entry points at a local solve, keeping the native
    ``run`` and ``run_async``.
    """
    from tidy3d import config
    from tidy3d.web.api.autograd import engine, hooks, strategy
    from tidy3d.web.api.autograd.autograd import run as _autograd_run
    from tidy3d.web.api.autograd.autograd import run_async as _autograd_run_async

    # Keep tidy3d's native autograd ``run``, the one wrapped in primitive: a ``web.run(sim)``
    # carrying tracers is handed back to it, which extracts the traced field map forward and builds
    # the adjoint backward. It still calls the ``_run_tidy3d`` we replace below, so the solve stays
    # local and nothing goes to the cloud.
    _STATE["autograd_run"] = _autograd_run
    _STATE["autograd_run_async"] = _autograd_run_async
    config.adjoint.local_gradient = True
    # The call site is in hooks, which captured its own reference at import time, so all three
    # modules have to be patched
    for mod in (engine, hooks, strategy):
        if hasattr(mod, "_run_tidy3d"):
            mod._run_tidy3d = _run_one
        if hasattr(mod, "_run_async_tidy3d"):
            mod._run_async_tidy3d = _run_many


def _patch_containers(Batch, Job) -> None:
    """Replace run, load, start, monitor and get_info of ``web.Job`` and ``web.Batch`` with local
    stand-ins.
    """
    Job.run = _job_run
    Batch.run = _batch_run
    Batch.load = _batch_load
    # One notebook reads .taskName and .taskId from batch.get_info.values. The rest of the web
    # container interface is a no-op under a local solve (start, monitor), and load is the same as
    # run.
    Batch.get_info = _batch_get_info
    Batch.start = _noop
    Batch.monitor = _noop
    Job.load = _job_run
    Job.start = _noop
    Job.monitor = _noop
    Job.get_info = _job_get_info


def _patch_costs(webapi, Batch, Job) -> None:
    """Billing and query interfaces: a local solve has no task_id and costs nothing, so these all
    short-circuit to placeholder values.

    Otherwise real_cost and estimate_cost go through SimulationTask.get and hit a real 404 in the
    cloud. The stand-in is the one in nb.backend._zero_cost: a local solve is not billed, and it
    returns a non-zero placeholder of 1e-3.
    """
    from openem.nb.backend import _zero_cost
    Batch.real_cost = _zero_cost
    Batch.estimate_cost = _zero_cost
    Job.real_cost = _zero_cost
    Job.estimate_cost = _zero_cost
    webapi.estimate_cost = _zero_cost
    # ``web.estimate_cost(job.task_id)``: task_id is a lazily uploaded attribute, and merely
    # touching it goes to the cloud.
    with contextlib.suppress(Exception):
        Job.task_id = property(lambda self: "openem-local")
    try:
        import tidy3d.web as _w
        _w.real_cost = _zero_cost
        _w.estimate_cost = _zero_cost
        _w.run = _run_plain
    except Exception:
        pass


def install(steps: int | None = None, max_cells: int | None = 80_000_000) -> None:
    """Install the hooks. After this, ``web.run`` and ``ag.grad`` both land on OpenEM.

    Args:
        steps: a fixed step count; ``None`` means the scene's nominal count with early termination
            enabled. **A forward and its adjoint must run the same number of steps**, since the
            adjoint is the transpose of the same discrete operator and a mismatch makes the result
            not the gradient of that objective. Early termination here is judged per scene
            independently, so pass ``steps`` explicitly when a strict comparison is needed.
        max_cells: per-simulation cell ceiling; exceeding it raises, which stops one large example
            in a notebook from filling the GPU. ``None`` means no limit.
    """
    from tidy3d.web.api import webapi
    from tidy3d.web.api.container import Batch, Job

    _STATE["steps"] = steps
    _STATE["max_cells"] = max_cells
    _patch_autograd_entrypoints()
    webapi.run = _run_plain
    _patch_containers(Batch, Job)
    _patch_costs(webapi, Batch, Job)
    # Geometry gradient patches (openem/nb/shapegrad.py): zero the gradient of a mirror-side body;
    # with OPENEM_SHAPEGRAD=eps, geometry gradients switch to the eps-map perturbation; =0 disables
    from openem.nb import shapegrad as _sg
    if _sg.apply():
        print(f"[openem] geometry gradient patch shapegrad.apply (mode {_sg._STATE['mode']})", flush=True)

