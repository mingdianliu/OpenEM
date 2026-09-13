# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""OpenEM's ``tidy3d.web.run``: solve a Tidy3D ``Simulation`` on the local GPU and return a real
``td.SimulationData``.

    import openem
    openem.install                      # from here on, td.web.run and autograd run locally
    sim_data = td.web.run(sim, task_name="demo")
    # or call it directly, without installing the hook
    sim_data = openem.run(sim)

The flow: ``td.Simulation`` -> ``scene.from_simulation`` (grid, geometry and materials taken
directly from what the Tidy3D client already computed) -> ``solver.run`` (with early termination)
-> ``normalize``, a normalization with no fitting -> assembled into ``FluxData``, ``FieldData``,
``ModeData`` and the rest. Downstream APIs such as ``sim_data.plot_field`` and
``sim_data["mon"].Ex`` work unchanged.

Environment variables, all optional:
- ``OPENEM_WORK_DIR``   root of the per-solve working directory (default ``./openem_runs``); a
  second run of the same scene hits the cache and does not recompute.
- ``OPENEM_SUBPROCESS`` set to 1 to run each solve in a subprocess, isolating device memory in a
  long notebook; in-process by default.
- ``OPENEM_RESULTS_DIR`` / ``OPENEM_RESULT_CASE``  optional: reuse an already-computed
  ``simdata.hdf5`` by content hash, for regression comparisons.
"""
from __future__ import annotations

import pathlib
import os
import re
import time

import numpy as np


import tidy3d as td

#: Root of the per-solve working directory, by default openem_runs/ under the current directory.
#: scene.npz and ours.npz land here, and an ours.npz newer than its scene.npz is reused without
#: recomputing.
WORK_DIR = pathlib.Path(os.environ.get("OPENEM_WORK_DIR", "openem_runs"))
#: Optional directory of reusable results: for a regression comparison, look up an
#: already-computed simdata.hdf5 by content hash. Unset means no lookup.
_r = os.environ.get("OPENEM_RESULTS_DIR")
RESULTS_DIR = pathlib.Path(_r) if _r else None
#: Root of the disk cache for mode solves and the like, by default .openem_cache/ under the current
#: directory.
CACHE_DIR = pathlib.Path(os.environ.get("OPENEM_CACHE_DIR", ".openem_cache"))

from openem import normalize, serialize
from openem import flux as flux_mod
from openem import grid as grid_mod
from openem.grid import YEE_ON_EDGE
from openem.model import COMPONENTS
from openem import scene as scene_mod

_LAST_QUEUE: dict = {}




_RESULT_CASE = None


def _zero_cost(*_a, **_k):
    """Stand-in for ``web.estimate_cost`` and ``real_cost``: a local solve is not billed, and the
    notebook only prints the number.

    This is the only copy; autograd_hook._patch_costs installs this same function on Job, Batch and
    webapi.
    """
    print("[openem] estimate_cost/real_cost: the simulation runs on a local GPU, so cloud billing "
          "does not apply")
    return 1e-3   # non-zero placeholder, so a cost ratio or division in a notebook does not hit ZeroDivision


def _shim_web_billing() -> None:
    """Stub out the billing and account queries in tidy3d.web; the simulation runs locally and there
    is no cloud bill.

    Only the **queries** are stubbed (estimate_cost, real_cost). run and upload are left alone: they
    are replaced by openem_run and openem_Job, so a genuine submission call in a notebook fails
    visibly rather than being swallowed.
    """
    try:
        import tidy3d.web as _w
    except Exception:
        return
    _w.estimate_cost = _zero_cost
    if hasattr(_w, "real_cost"):
        _w.real_cost = _zero_cost



def _shim_design_batch() -> None:
    """Replace the batch executor of tidy3d.plugins.design with one openem_run per simulation.

    Inside DesignSpace.run(fn=...), `_run_batch(batch, ...)` calls `batch.run` and goes straight to
    the cloud. It contains neither "web." nor "Batch" as text, so a textual substitution cannot
    intercept it, and one notebook raised a WebError once the account balance expired. After this
    replacement each sim goes through openem_run, and a plain dict is returned, which is what the
    plugin consumes by key (post_input_loader(key, value) in design.py).

    Note that openem_run is referenced through the **module global** here, so that when a capture
    replaces nb.backend.openem_run with a stub, this shim follows the stub too.
    """
    try:
        from tidy3d.plugins.design.design import DesignSpace as _DS
    except Exception:
        return

    class _LocalBatchData(dict):
        """A thin local shell around BatchData: the plugin wants task_ids and task_paths as dicts
        keyed by task name.
        """

        @property
        def task_ids(self):
            return {k: f"openem-local-{k}" for k in self}

        @property
        def task_paths(self):
            return {k: "" for k in self}

    def _local_batch(batch, path_dir=None, priority=None):
        sims = getattr(batch, "simulations", {}) or {}
        print(f"[openem] DesignSpace sweep of {len(sims)} simulations: solved locally, one at a "
              "time, nothing goes to the cloud")
        return _LocalBatchData(
            {nm: openem_run(sim, task_name=nm) for nm, sim in sims.items()})

    _DS._run_batch = staticmethod(_local_batch)


def _shim_web_batch() -> None:
    """Replace tidy3d.web.Batch with a local stand-in.

    Some notebooks call `td.web.Batch(...).run` directly, bypassing web.run, where a textual
    substitution cannot intercept it.
    """
    try:
        import tidy3d.web as _w
    except Exception:
        return

    class _LocalBatch:
        def __init__(self, simulations=None, **kw):
            self.simulations = dict(simulations or {})

        def run(self, path_dir=None, priority=None, **kw):
            print(f"[openem] Batch: {len(self.simulations)} simulations solved locally, one at a "
                  "time, nothing goes to the cloud")
            from tidy3d.plugins.design.design import DesignSpace as _DS  # noqa: F401
            data = {nm: openem_run(sim, task_name=nm)
                    for nm, sim in self.simulations.items()}
            out = dict(data)
            if path_dir:
                # Tidy3D's Batch.run writes batch.hdf5 into path_dir, and notebooks often read it
                # back with web.Batch.from_file(path_dir/batch.hdf5).load(path_dir); one died on
                # FileNotFoundError. Our to_file writes a json index under the same name, which both
                # from_file and load accept.
                try:
                    self.to_file(os.path.join(path_dir, "batch.hdf5"))
                except Exception as _e:
                    print(f"[openem] failed to write the Batch index batch.hdf5: {str(_e)[:80]}")

            class _BD(dict):
                @property
                def task_ids(self):
                    return {k: f"openem-local-{k}" for k in self}

                @property
                def task_paths(self):
                    return {k: "" for k in self}

            return _BD(out)

        def upload(self, *a, **k):
            pass

        def estimate_cost(self, *a, **k):
            return 1e-3

        # One notebook reads .taskName and .taskId from batch.get_info.values; the rest are the
        # usual no-op parts of the web.Batch interface.
        def get_info(self, *a, **k):
            from types import SimpleNamespace
            return {nm: SimpleNamespace(taskName=nm, taskId=f"openem-local-{nm}",
                                        status="success")
                    for nm in (self.simulations or {})}

        def start(self, *a, **k):
            pass

        def monitor(self, *a, **k):
            pass

        def load(self, path_dir=None, **k):
            return self.run(path_dir=path_dir)

        def real_cost(self, *a, **k):
            return 1e-3

        def to_file(self, path, *a, **k):
            import json as _json
            os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
            base = os.path.abspath(path)
            paths = {}
            for nm, sim in (self.simulations or {}).items():
                p = f"{base}.{re.sub(r'[^A-Za-z0-9_.-]', '_', nm)}.openem_sim.json"
                sim.to_file(p); paths[nm] = p
            with open(path, "w") as f:
                _json.dump({"sims": paths}, f)

        @classmethod
        def from_file(cls, path, *a, **k):
            import json as _json
            d = _json.load(open(path))
            sims = {nm: td.Simulation.from_file(p) for nm, p in d.get("sims", {}).items()}
            return cls(simulations=sims)

    _w.Batch = _LocalBatch




def _with_courant(sim, c: float):
    """Set the courant of a Simulation, or of every Simulation inside a dict or list, to ``c``;
    anything that is not a Simulation is returned unchanged.
    """
    if isinstance(sim, dict):
        return {k: _with_courant(v, c) for k, v in sim.items()}
    if isinstance(sim, (list, tuple)):
        return type(sim)(_with_courant(v, c) for v in sim)
    try:
        return sim.updated_copy(courant=c) if type(sim).__name__ == "Simulation" else sim
    except Exception:
        return sim


def _is_divergence(exc: BaseException) -> bool:
    m = str(exc)
    return isinstance(exc, DivergenceError) or "solve output contains NaN/Inf" in m or "looks divergent" in m


#: Fallbacks after a divergence, retried in order; the first rung that does not diverge is used, and
#: if all fail, failed_data marks it failed. Each entry is ``(task name suffix, courant,
#: staircase?)``.
#:
#: Why the order lowers Courant first and only staircases last: on one optimization member,
#: measured values were 1.593 dB as given, 1.594 at Courant 0.8 and 1.592 at Courant 0.6, while
#: staircasing gave NaN and staircasing plus 0.6 gave 4.095 dB at |amp| 0.75. Staircasing a 0.43 um
#: wide waveguide onto a 0.029 um grid manufactures 2.5 dB of false loss, which is not the same
#: physical problem. The old order (staircase, then staircase plus 0.6) replaced the values of three
#: members of generation 0 with staircased ones. The divergence itself, a scheme-level instability
#: between the silicon poles and the psi-PML in the deepest layer of the +x PML, is suppressed by
#: shrinking dt alone. Staircasing stays as the last rung, as a backstop for scenes where the
#: dispersive interface cells are themselves unstable.
DIVERGENCE_LADDER = (("c08", 0.8, False), ("c06", 0.6, False), ("stair-c06", 0.6, True))


def divergence_retries(sim):
    """Yield ``(suffix, sim with adjusted courant, staircase?)`` in the order of
    :data:`DIVERGENCE_LADDER`.

    A rung is skipped when the user already set courant below it, since rerunning at the same dt
    proves nothing.
    """
    cur = getattr(sim, "courant", None)
    for suffix, c, stair in DIVERGENCE_LADDER:
        if cur is not None and cur <= c + 1e-9 and not stair:
            continue
        yield suffix, (_with_courant(sim, c) if cur is None or cur > c + 1e-9 else sim), stair


def _ladder_log(tag, last, suffix, stair):
    """Log text for the ladder in openem_run; autograd_hook._solve passes its own wording."""
    print(f"[openem] the solve diverged ({str(last).splitlines()[0][:60]}...); retrying with "
          f"fallback {suffix}"
          f"{' (OPENEM_DISP_STAIRCASE=1, dispersive interface cells staircased)' if stair else ''}"
          f", task name suffixed -{suffix}", flush=True)


def run_with_divergence_ladder(sim, run_once, tag, first_exc=None, log=_ladder_log):
    """The divergence fallback ladder, shared by openem_run and autograd_hook._solve.

    ``run_once(sim, suffix)`` runs once and returns the data. The first attempt (``suffix=""``)
    happens here; a caller that has already run and hit a divergence passes that exception in as
    ``first_exc``. (autograd_hook's first solve happens after the scene is built and is tied to its
    CALLS bookkeeping, which makes wrapping it in run_once awkward.)

    Anything that is not a divergence, or a divergence with OPENEM_DISP_STAIRCASE=1 already on, is
    re-raised unchanged. Otherwise the rungs of :data:`DIVERGENCE_LADDER` are retried in order:
    Courant 0.8, then 0.6, then staircasing plus 0.6. Frequency-domain results do not depend on dt,
    and the source-spectrum normalization follows automatically. If every rung diverges,
    :func:`failed_data` marks it failed by zeroing the monitor data, and the optimizer discards that
    design. The staircasing rung works by setting and clearing the ``OPENEM_DISP_STAIRCASE``
    environment variable. ``log(tag, last, suffix, stair)`` only prints, and each call site keeps
    its own wording.
    """
    if first_exc is None:
        try:
            return run_once(sim, "")
        except RuntimeError as _e:
            first_exc = _e
    if not _is_divergence(first_exc) or os.environ.get("OPENEM_DISP_STAIRCASE") == "1":
        raise first_exc
    data, last = None, first_exc
    for suffix, sim2, stair in divergence_retries(sim):
        log(tag, last, suffix, stair)
        if stair:
            os.environ["OPENEM_DISP_STAIRCASE"] = "1"
        try:
            data = run_once(sim2, suffix)
            break
        except RuntimeError as _e2:
            if not _is_divergence(_e2):
                raise
            last = _e2
        finally:
            os.environ.pop("OPENEM_DISP_STAIRCASE", None)
    if data is None:
        data = failed_data(last, tag)
    return data


def _traced_run(sim, task_name, verbose):
    """A simulation carrying autograd tracers (ArrayBox) is handed back to tidy3d's native autograd
    ``run``; without tracers this returns None.

    When a notebook calls ``web.run`` under ``ag.grad`` or ``plugins.autograd.optimize``, the entry
    point has been rewritten to ``openem_run``, which bypasses the dispatch in autograd_hook's
    ``_run_plain``, goes straight to ``_openem_run_inner`` and returns a plain SimulationData with
    **no tracer**. autograd then reports "Output seems independent of input" and the gradient is
    identically zero; one case ran 40 steps with the gradient norm at zero throughout and J never
    moving. tidy3d's native ``run`` is a ``@primitive``: forward it extracts the traced field map,
    backward it builds the adjoint source and assembles the gradient. Its internal ``_run_tidy3d``
    has already been replaced by ``autograd_hook.install`` with a local solve, so nothing goes to
    the cloud. Both the forward and the adjoint have been through ``to_static`` by the time they
    reach ``_run_tidy3d``, so they carry no tracer when they come back here and there is no
    recursion.
    """
    from openem.nb import autograd_hook as _ah
    if isinstance(sim, dict):
        sims = list(sim.values())
    elif isinstance(sim, (list, tuple)):
        sims = list(sim)
    else:
        sims = [sim]
    if not any(_ah._has_tracers(s) for s in sims if s is not None):
        return None
    if _ah._STATE.get("autograd_run") is None:
        # The hook is not installed, i.e. we did not come in through the notebook kernel startup
        # script, so install it now with the same parameters
        _s = os.environ.get("OPENEM_NB_STEPS")
        _ah.install(steps=int(_s) if _s else None,
                    max_cells=int(os.environ.get("OPENEM_NB_MAX_CELLS", "80000000")))
    if verbose:
        print(f"[openem] {task_name}: the simulation carries autograd tracers, handing it back to "
              "tidy3d's native autograd run (still solved locally)")
    return _ah._run_plain(sim, task_name=task_name, verbose=verbose)



def _openem_run_inner(sim: "td.Simulation" = None, task_name: str = "openem",
                      verbose: bool = True, **kw):
    """OpenEM's stand-in for ``web.run``, solving on the local GPU.

    Each route goes its own way: a sequence of simulations recurses one at a time; a ModeSolver or
    ModeSimulation goes to tidy3d's local mode solver; an already-computed simdata is read
    directly; everything else builds a scene and solves.
    """
    # Some notebooks pass web.run the simulation through the simulation= keyword
    if sim is None:
        sim = kw.get("simulation")
    if sim is None:
        raise TypeError("openem.run needs sim or simulation=")
    kw2 = {k: v for k, v in kw.items() if k != "simulation"}
    # A sequence of simulations: each takes the same route, and the list is returned in order.
    # Every level of a nested list appends its index to the task name; otherwise different
    # simulations share a name and write into the same working directory.
    if isinstance(sim, (list, tuple)):
        return [run(s, task_name=f"{task_name}-{i}", verbose=verbose, **kw2) for i, s in enumerate(sim)]
    # Batch style, {task name: sim}: some notebooks hand a dict straight to web.run. Returns a thin
    # shell dict carrying task_ids and task_paths, matching what _LocalBatch.run returns.
    if isinstance(sim, dict):
        out = {nm: run(s, task_name=nm, verbose=verbose, **kw2) for nm, s in sim.items()}

        class _BD(dict):
            @property
            def task_ids(self): return {k: f"openem-local-{k}" for k in self}
            @property
            def task_paths(self): return {k: "" for k in self}
        return _BD(out)
    # Match on the **suffix** of the name: the waveguide plugin passes a wrapper or subclass, and an
    # exact match would miss it.
    if type(sim).__name__.endswith("EMESimulation"):
        raise NotImplementedError(
            "the EME (eigenmode expansion) solver is out of scope for OpenEM, which is pure FDTD. "
            "A notebook containing EME can only have its FDTD part compared here.")
    if type(sim).__name__.endswith("ModeSolver"):
        return _mode_solver_data(sim, verbose)
    if type(sim).__name__.endswith("ModeSimulation"):
        return _mode_simulation_data(sim, verbose)
    if hasattr(sim, "sim_dict") and type(sim).__name__.endswith("ComponentModeler"):
        # The smatrix plugin: hand it back to tidy3d's own _run_local, which runs port by port
        # through web.Batch, and Batch.run has already been replaced by the hook with a local solve,
        # so each port simulation still goes through our solver.
        if verbose:
            print("[openem] ComponentModeler: solving locally, one port at a time")
        from tidy3d.plugins.smatrix.run import _run_local
        return _run_local(sim, verbose=bool(verbose))
    cached = _cached_simdata(sim, task_name, verbose)
    if cached is not None:
        return cached
    # The working directory carries an optional notebook tag (OPENEM_NB_TAG): when several run in
    # parallel their task names are often all sim_0.
    _tag = re.sub(r"[^a-z0-9-]", "-", os.environ.get("OPENEM_NB_TAG", "").lower())[:14]
    label = "openem-nb-" + (_tag + "-" if _tag else "") + re.sub(r"[^a-z0-9-]", "-", task_name.lower())[:24]
    work = WORK_DIR / label
    work.mkdir(parents=True, exist_ok=True)
    # Build the tables from the in-memory sim rather than round-tripping through a file: tidy3d's
    # JSON serialization silently drops the spatial data of CustomMedium and TriangleMesh. The hdf5
    # written to disk is only a record; downstream reads scene.npz.
    simf = work / "simulation.hdf5"
    sim.to_file(str(simf))
    try:
        sc = scene_mod.from_simulation(sim)
    except NotImplementedError as exc:
        if "subpixel" not in str(exc):
            raise
        if verbose:
            print(f"[openem] the builder does not support subpixel with this medium combination "
                  f"({str(exc)[:40]}...); rebuilding with subpixel=False, which is physically "
                  f"equivalent against a uniform background")
        sim = sim.updated_copy(subpixel=False)
        sim.to_file(str(simf))
        sc = scene_mod.from_simulation(sim)
    serialize.save(sc, work / "scene.npz")
    oursf = work / "ours.npz"
    # Cache: an ours.npz newer than its scene.npz is reused directly
    if oursf.exists() and oursf.stat().st_mtime >= (work / "scene.npz").stat().st_mtime:
        if verbose:
            print(f"[openem] cache hit on {oursf.name}, skipping the solve")
    else:
        if verbose:
            print(f"[openem] scene {sc.shape} = {sc.grid.num_cells:,} cells, nominal steps "
                  f"{sc.num_time_steps:,}; solving...")
        from openem.nb import solve_worker
        solve_worker.solve_dir(work, work, subprocess_=os.environ.get("OPENEM_SUBPROCESS") == "1")
        if not oursf.exists():
            raise RuntimeError(f"the solve produced no {oursf}")
    ours = np.load(oursf)
    if verbose:
        import json as _json
        meta = _json.loads(ours["__meta"].item())
        print(f"[openem] done: {meta['steps_run']:,} steps, early stop {meta.get('stop') or 'none'}, "
              f"GPU {meta['wall']:.1f}s")
    # Assembly is left to build, which covers flux, field, mode, time-domain and projection monitors
    # alike.
    return build(sim, sc, ours)


@scene_mod.media.cloud_emulation()
def run(sim: "td.Simulation" = None, task_name: str = "openem",
                 verbose: bool = True, **kw):
    """The web.run stand-in entry point: goes through _openem_run_inner, and when the notebook passes
    path=, writes the SimulationData to disk, since some notebooks read that file back later.
    """
    path = kw.pop("path", None)
    traced = _traced_run(sim if sim is not None else kw.get("simulation"), task_name, verbose)
    if traced is not None:
        return traced
    # The divergence fallbacks (Courant 0.8, then 0.6, then staircasing plus 0.6, and failed_data if
    # all diverge) live in run_with_divergence_ladder
    data = run_with_divergence_ladder(
        sim,
        lambda s, suffix: _openem_run_inner(
            s, task_name=f"{task_name}-{suffix}" if suffix else task_name, verbose=verbose, **kw),
        task_name)
    if path and not isinstance(data, (list, tuple)):
        try:
            import pathlib as _pl
            _pl.Path(path).parent.mkdir(parents=True, exist_ok=True)
            data.to_file(str(path))
        except Exception as _e:
            print(f"[openem] failed to write the archive at path (the return value is unaffected): {_e}")
    return data


#: Former name, kept for the ``openem_run(...)`` spelling in earlier notebooks.
openem_run = run

_shim_web_billing()
_shim_design_batch()
_shim_web_batch()
def _mode_simulation_data(msim, verbose):
    """A ``td.ModeSimulation``, which some notebooks submit as a remote mode solve through
    ``web.run(simulation=td.ModeSimulation.from_simulation(...))``.

    It goes to tidy3d's own ``run_local``, with the same disk cache as ModeSolver. It is not an FDTD
    simulation (it has no source) and entering the solver would only die on "no reference frequency
    available".
    """
    if verbose:
        print("[openem] this is a ModeSimulation rather than an FDTD simulation; using tidy3d's "
              "local mode solver")
    import hashlib as _hl
    from tidy3d.components.mode.data.sim_data import ModeSimulationData
    # The cache skeleton is shared with scene/modes (modes.disk_cached); the key algorithm is
    # unchanged (an msim- prefix plus sha1(_sim_cache_key)[:20])
    key = "msim-" + _hl.sha1(_sim_cache_key(msim).encode()).hexdigest()[:20]
    return scene_mod.modes.disk_cached(key, msim.run_local, ModeSimulationData.from_file,
                                       default_dir=CACHE_DIR / "modes", miss_note="solved locally and written to disk")


def _mode_solver_data(sim, verbose):
    """A ModeSolver goes to tidy3d's local mode solver; it is not FDTD and never reaches our solver.

    It carries a disk cache keyed by the content hash of the solver, reusable across processes, with
    the same skeleton as scene/modes.disk_cached: a cold solve takes twenty-odd minutes and a hit
    takes twenty seconds. A failure on the cache path is not swallowed, since that amounts to
    re-solving every time, which once pushed a notebook past its cell timeout with no visible
    reason.
    """
    if verbose:
        print("[openem] this is a ModeSolver rather than an FDTD simulation; using tidy3d's local "
              "mode solver")
    try:
        import hashlib as _hl
        from tidy3d.plugins.mode.mode_solver import ModeSolverData
        # The cache skeleton is shared with scene/modes (modes.disk_cached); the key algorithm is
        # unchanged (sha1(_sim_cache_key)[:20])
        key = _hl.sha1(_sim_cache_key(sim).encode()).hexdigest()[:20]
        return scene_mod.modes.disk_cached(key, lambda: sim.data, ModeSolverData.from_file,
                                           default_dir=CACHE_DIR / "modes", miss_note="solved locally and written to disk")
    except Exception as exc:
        print(f"[openem] the mode cache is unavailable ({type(exc).__name__}: "
              f"{str(exc)[:80]}); solving directly")
    return sim.data


def _cached_simdata(sim, task_name, verbose):
    """Look for an already-computed simdata; None means none was found and the caller should solve.

    Dispatch order: by the content hash of the simulation (the most reliable for a notebook with
    several simulations, since task names and variable names can repeat), then a task_name
    subdirectory, then a single simdata.hdf5 (the single-simulation case, whose behaviour is
    unchanged).

    The latter two are fallback paths and may point at **another** simulation's result. So when a
    by_hash directory exists, meaning several simulations, **and this simulation's entry is
    missing, it raises** rather than falling back. A single-simulation case has no by_hash directory
    and behaves as before. A fallback result additionally has its monitor names checked for
    completeness. Quietly returning another simulation's result means the figures and numbers are
    already wrong while the notebook says nothing.
    """
    if not _RESULT_CASE or RESULTS_DIR is None:
        return None
    import hashlib as _hl
    h = _hl.sha1(_sim_cache_key(sim).encode()).hexdigest()[:16]
    byh = RESULTS_DIR / _RESULT_CASE / "by_hash"
    exact = byh / h / "simdata.hdf5"
    if byh.is_dir() and not exact.exists():
        # A multi-simulation case, which is what a by_hash directory means, **must stop** when it
        # cannot reach the exact entry. With duplicate monitor names the two fallback paths below
        # catch nothing: in one notebook all three simulations named their monitors flux, field_yz
        # and field_xy, one of them received another's fields, and the whole "Reference" curve in
        # the figure was rotated by 60 degrees while the notebook said nothing.
        note = byh / h / "UNSUPPORTED.txt"
        why = (f"it was refused at table-build time: {note.read_text(encoding='utf-8').strip()[:300]}"
               if note.exists() else
               f"there is no simdata.hdf5 under by_hash/{h}/")
        raise KeyError(
            f"{_RESULT_CASE} has several simulations, but this one (content hash {h}) has no "
            f"result.\n"
            f"{why}\n"
            f"A fallback path would return **another** simulation's result, which cannot be "
            f"detected when monitor names collide and would corrupt the figures silently, so this "
            f"stops here instead.")
    for sdf in (exact,
                RESULTS_DIR / _RESULT_CASE / task_name / "simdata.hdf5",
                RESULTS_DIR / _RESULT_CASE / "simdata.hdf5"):
        if not sdf.exists():
            continue
        if verbose:
            print(f"[openem] reading {sdf}"
                  f" (reusing the coordinates with OpenEM values, nothing is recomputed)")
        import tidy3d as _td
        sd = _td.SimulationData.from_file(str(sdf))
        want = {m.name for m in sim.monitors}
        have = set()
        for d in sd.data:
            nm = getattr(getattr(d, "monitor", None), "name", None)
            if nm:
                have.add(nm)
        miss = want - have
        if miss and not os.environ.get("OPENEM_NB_ALLOW_PARTIAL"):
            raise KeyError(
                f"this simulation needs {sorted(want)}, but {sdf} only has "
                f"{sorted(have)} and is missing {sorted(miss)}.\n"
                f"The content hash of this simulation is {h}, so there is no result under "
                f"by_hash/{h}/ and it fell back to a different simdata.\n"
                f"(To use an incomplete result deliberately, set OPENEM_NB_ALLOW_PARTIAL=1.)")
        return sd
    return None


def _data_array_digest(model) -> str:
    """Walk a tidy3d model and feed the values and coordinates of every xarray.DataArray (CustomMedium
    pixels, Custom source datasets, ModeSpec and so on) into a sha1; an empty string when there are
    no arrays.

    tidy3d's sim.json stores those arrays as hdf5 references and keeps them out of the JSON, so
    changing only the pixels leaves the JSON identical. One 75-iteration optimization turned out to
    have genuinely computed its forward and adjoint solves once each, with the other 148 all hitting
    the first iteration's cache.
    """
    import hashlib as _hl
    try:
        import xarray as _xr
    except Exception:
        return ""
    h = _hl.sha1(); seen = [0]

    def walk(o, depth=0):
        if depth > 12:
            return
        if isinstance(o, _xr.DataArray):
            seen[0] += 1
            h.update(np.ascontiguousarray(np.asarray(o.values)).tobytes())
            for k in sorted(o.coords):
                h.update(str(k).encode()); h.update(np.ascontiguousarray(np.asarray(o.coords[k].values)).tobytes())
            return
        if isinstance(o, dict):
            for k in sorted(o, key=str):
                walk(o[k], depth + 1)
            return
        if isinstance(o, (list, tuple)):
            for v in o:
                walk(v, depth + 1)
            return
        d = getattr(o, "__dict__", None)
        if d and hasattr(o, "__fields__") or (d and type(o).__module__.startswith("tidy3d")):
            for k in sorted(d, key=str):
                if not str(k).startswith("_"):
                    walk(d[k], depth + 1)

    walk(model)
    return f"|data:{h.hexdigest()}" if seen[0] else ""


def _sim_cache_key(sim, extra: str = "") -> str:
    """Ingredients of the simulation cache key: json | extra | a digest of the data arrays. With no
    arrays it is exactly the old key, so existing caches stay valid.
    """
    return sim.json() + extra + _data_array_digest(sim)


def _load_simdata(path: pathlib.Path):
    """Read the simdata.hdf5 a job wrote back, retrying with backoff since a network filesystem may
    expose it before the write completes.
    """
    for _try in range(8):
        try:
            return td.SimulationData.from_file(str(path))
        except Exception as _e:
            if _try == 7:
                raise RuntimeError(f"cannot read the job result {path} ({type(_e).__name__}: {str(_e)[:120]})") from _e
            time.sleep(5 + 5 * _try)


# ---- Native SimulationData assembler, structurally identical to Tidy3D data.hdf5 ----
from tidy3d.components.data.data_array import (
    ScalarFieldDataArray, FluxDataArray, FluxTimeDataArray,
    ModeAmpsDataArray, ModeIndexDataArray)

# Both the component order and the Yee staggering table come from their single source on the solver
# side: model.COMPONENTS and grid.YEE_ON_EDGE (True means that axis lands on a primal boundary,
# which is what the 1 in the old _YEE table meant).


def _comp_coords(grid, comp, colocate=False):
    """Target (x, y, z) coordinates of this component, in um.

    ``colocate=False``: each component at its own native Yee position, the same as Tidy3D's
    semantics.
    ``colocate=True``: Tidy3D's output convention, where all six components are interpolated onto
    ``boundaries[:-1]`` (its ``colocation_boundaries`` rule) and a zero-thickness axis snaps to the
    center. **They must be uniform**: otherwise the six components have different shapes and a
    cross-component operation in a notebook fails to broadcast.
    """
    out = []
    for ax in range(3):
        b = np.asarray(grid.boundaries.to_list[ax], dtype=float)
        c = np.asarray(grid.centers.to_list[ax], dtype=float)
        if colocate:
            a = c if b.size <= 2 else b[:-1]     # zero-thickness axis (boundaries are just [z,z]) -> center
        else:
            # Tidy3D stores native Yee positions with the boundary dimension **dropping its last
            # entry** (the drop-last rule in colocation_boundaries; measured on a reference output as
            # Ex = (70 centers, 70 boundaries, 70 boundaries)). A zero-thickness axis, where b is
            # just [z,z], stays a single point.
            if YEE_ON_EDGE[comp][ax]:
                a = b if b.size <= 2 else b[:-1]
            else:
                a = c
        out.append(np.unique(a))                 # already monotonic; deduplicates a degenerate axis
    return out


def _field_data_td(sim, mon, ours, nm, freqs):
    """Direct copy into FieldData: ``tdf:<name>|<component>`` is already in Tidy3D's convention,
    (x, y, z, f).

    A scene where the reference solver did no source normalization has to have its source spectrum
    multiplied back; see :func:`_denorm_spectrum`.
    """
    grid = sim.discretize_monitor(mon)
    s = _denorm_spectrum(sim, freqs)
    comps = {}
    for comp in COMPONENTS:
        k = f"tdf:{nm}|{comp}"
        if k not in ours:
            continue
        co = {ax: ours[f"tdc:{nm}|{comp}:{ax}"] for ax in "xyz"}
        vals = ours[k] if s is None else ours[k] * s
        comps[comp] = ScalarFieldDataArray(
            vals, coords=dict(x=co["x"], y=co["y"], z=co["z"], f=freqs))
    return td.FieldData(monitor=mon, grid_expanded=grid, **comps)


def _snap_thin_axis_target(t_m, srcax, edges_m):
    """On a zero-thickness axis, where the target is a single point, move the interpolation target to
    the nearest primal boundary, i.e. the E layer.

    This matches tidy3d's "take the layer, label it with the center" behaviour; the caller keeps the
    label. OPENEM_COLOC_SNAP=0 disables it.
    """
    import os as _os
    # Its own knob, OPENEM_FIELDMON_SNAP, default 0: across the real comparison cases it improved 5
    # field rows and worsened 13, so it was rolled back. Snapping on mode planes goes through a
    # separate knob, OPENEM_COLOC_SNAP, and is unaffected.
    if t_m.size != 1 or _os.environ.get("OPENEM_FIELDMON_SNAP", "0") != "1":
        return t_m
    e = np.asarray(edges_m, dtype=np.float64)
    if e.size < 2:
        return t_m
    # tidy3d's rule: take the boundary that is not greater than the center (ind_min in
    # discretize_inds is the last boundary <= center), falling back to the nearest one when the
    # center is below every boundary. An earlier version used "nearest", which picked a different
    # layer than tidy3d whenever the offset exceeded half a cell.
    from openem import colocate as _colocate
    snapped = _colocate.snap_floor(e, float(t_m[0]), 1e-15)      # tolerance 1e-15 m (scene/projection uses 1e-12 um)
    lo, hi = float(np.min(srcax)), float(np.max(srcax))
    return np.array([min(max(snapped, lo), hi)])


def _field_data(sim, sc, mon, raw, freqs, norm):
    """FieldData: raw (6, nf, ni, nj, nk) turned into six ScalarFieldDataArray."""
    grid = sim.discretize_monitor(mon)
    from openem import colocate
    comps = {}
    _colo = bool(getattr(mon, "colocate", False))
    for ci, comp in enumerate(COMPONENTS):
        xc, yc, zc = _apply_interval_space(
            mon, *_comp_coords(grid, comp, colocate=_colo))
        # OpenEM stores its own native Yee positions; interpolate onto this monitor's component
        # coordinates, in metres
        src = colocate.sample_coords(sc.grid, mon_origin(sc, mon),
                                     mon_box(sc, mon), comp)
        # Clip the target coordinates into the sampled range; the monitor grid may extend half a cell
        # beyond our storage box
        tgt = []
        for ai, (s_um, srcax) in enumerate(zip((xc, yc, zc), src)):
            s_m = s_um * 1e-6
            lo, hi = float(srcax.min()), float(srcax.max())
            t_m = np.clip(s_m, lo, hi)
            # Zero-thickness axis: snap the interpolation target to the nearest E layer (see
            # _snap_thin_axis_target); the output coordinate keeps its nominal value.
            # Snap only when colocate=True: with colocate=False, H sits at cell centers and must not
            # be interpolated onto the E layer.
            tgt.append(_snap_thin_axis_target(t_m, srcax, sc.grid.axes[ai].edges) if _colo else t_m)
        vals = colocate.interp_to(raw[ci], src, tgt)
        vals = (vals / norm[:, None, None, None]).transpose(1, 2, 3, 0)
        comps[comp] = ScalarFieldDataArray(
            vals, coords=dict(x=xc, y=yc, z=zc, f=freqs))
    return td.FieldData(monitor=mon, grid_expanded=grid, **comps)


def _adjoint_coefficient(sim) -> complex | None:
    """The source coefficient to multiply back in an **adjoint simulation**; ``None`` when this is
    not one.

    tidy3d's autograd packs the complex coefficient ``dJ/d(amps)`` into the adjoint source's
    ``source_time.amplitude`` and ``phase`` and sets ``post_norm`` to **1.0**, meaning it expects
    the returned data to **carry** that coefficient.

    Our normalization divisor, however, is the source's own injection table
    (``normalize.source_table``), which divides it straight back out. The consequence is a mode
    objective whose gradient is off by 1/|A|: on one inverse-design case the ratio of the autograd
    sum to the finite difference measured **-31,256** (|A| = 6.10e-05, 1/|A| = 16,391, with the
    2.90 rad phase explaining the sign). Multiplying the coefficient back gives **+1.35**, the same
    range as the flux-objective comparison at 1.14 and 0.81.

    The test is **whether ``post_norm`` is a ``FreqDataArray``**, which is how tidy3d itself decides
    (``_testing/synthetic_monitor_data.py``: ``is_adjoint_sim = isinstance(simulation.post_norm,
    td.FreqDataArray)``). In an ordinary scene it is a float. ``simulation_type`` **cannot** be
    used: on the local-gradient path it is always ``"tidy3d"``. Nor can "the amplitude is not 1",
    since an ordinary scene is free to have a source amplitude other than 1, in which case the
    reference solver really does divide it out and the behaviour must not change.
    """
    # A probe against the reference solver: in an **ordinary** simulation with normalize_index=0 and
    # a source of amplitude=2.5 and phase=0.7, the normalized mode amplitude comes back as
    # 2.5 at +0.7 rad. Normalization divides by the spectrum of a **unit-amplitude** source and the
    # coefficient survives. Our divisor is the injection table, which includes the coefficient, so
    # whenever normalization happened at all (normalize_index is not None), the coefficient of the
    # normalize_index source has to be multiplied back, adjoint simulation or not. An adjoint
    # simulation (post_norm a FreqDataArray, normalize_sim=True) is only a special case of that
    # rule; multiplying only for it left the coefficient missing in ordinary simulations whose
    # source amplitude was not 1.
    if getattr(sim, "normalize_index", 0) is None:
        return None            # raw fields: the full source spectrum that _denorm_spectrum multiplies
                               # back already includes the coefficient
    if not sim.sources:
        return None
    ni = int(getattr(sim, "normalize_index", 0) or 0)
    ni = ni if 0 <= ni < len(sim.sources) else 0
    st = sim.sources[ni].source_time
    a = float(getattr(st, "amplitude", 1.0))
    p = float(getattr(st, "phase", 0.0))
    if abs(a - 1.0) < 1e-12 and abs(p) < 1e-12:
        return None
    return complex(a * np.exp(1j * p))


def _denorm_spectrum(sim, freqs):
    """The complex factor to multiply back at readout; ``None`` when none is needed.

    Two cases, which can stack:

    * the reference solver did no source normalization (``normalize_index is None``): multiply back
      the source spectrum ``S(f)``; see :func:`openem.normalize.td_source_spectrum` for the rule and
      its derivation.
    * **an adjoint simulation**: multiply back the source coefficient; see
      :func:`_adjoint_coefficient`.
    """
    out = None
    if getattr(sim, "normalize_index", 0) is None and sim.sources:
        out = normalize.td_source_spectrum(sim.sources[0].source_time, freqs,
                                           sim.dt, sim.num_time_steps)
    c = _adjoint_coefficient(sim)
    if c is not None:
        out = c if out is None else out * c
    return out


def _box_flux_data(sim, sc, mon, ours, box, faces):
    """Net outward power through a closed box: each face's flux times that face's outward normal
    sign, summed.

    A 3D FluxMonitor is split into six faces on the solver side (BoxFluxMonitor in
    scene/monitors.py), and they are added back here. The convention matches the one verified
    against the reference outputs.
    """
    m0 = faces[box.face_names[0]]
    fs = np.asarray(m0.freqs, float)
    vals = np.zeros(fs.size)
    for fn in box.face_names:
        f = faces[fn]
        vals = vals + flux_mod.monitor_flux(
            ours[f"mon:{fn}"], sc.grid, f) * f.normal_dir
    vals = vals / normalize.flux_norm(sc, fs)
    s = _denorm_spectrum(sim, fs)
    if s is not None:
        vals = vals * np.abs(s) ** 2
    return td.FluxData(monitor=mon,
                       flux=FluxDataArray(vals, coords=dict(f=fs)))


#: The (monitor name, order) pairs already warned about; each is reported once per process
_RAYLEIGH_WARNED: set = set()


def _rayleigh_check(sim, sc, spec, tol=0.02):
    """A plane flux covering a whole periodic face: check whether a diffraction order sits right on a
    Rayleigh cutoff.

    Along a periodic direction, order (m1, m2) has normal wavenumber kx^2 = (n*k0)^2 - k_t^2. A
    grazing order with kx approaching 0 carries power proportional to 1/kx through the Poynting
    surface integral, and the continuous solution **diverges** at kx = 0. The "total flux through
    the whole periodic face" then has no finite value at all, and what any discrete implementation
    reads is a regularized value set jointly by numerical dispersion, run_time and the reflection of
    grazing waves off the PML, which **does not converge with the grid**. One notebook's Lz was
    exactly 5.0000 lambda, putting orders (0, +-5) precisely at grazing.

    This only prints a warning and changes no number. ``OPENEM_RAYLEIGH_WARN=0`` disables it. The
    test has no free parameter: the orders follow from the period length and the cutoff from n*k0.
    """
    if os.environ.get("OPENEM_RAYLEIGH_WARN", "1") != "1":
        return
    ax = int(spec.axis)
    tax = grid_mod.transverse_axes(ax)
    Ls, bvs = [], []
    for p, t in enumerate(tax):
        axt = sc.grid.axes[t]
        if axt.boundary_lo not in ("Periodic", "BlochBoundary"):
            return                                  # not a periodic axis, so there are no orders
        if axt.is_flat:
            return
        lo, hi = (spec.t_bounds[p] if getattr(spec, "t_bounds", None) is not None
                  else (float(axt.edges[0]), float(axt.edges[-1])))
        L = float(axt.edges[-1] - axt.edges[0])
        if (hi - lo) < 0.999 * L:                   # does not cover the periodic face, so the order
                                                    # decomposition does not hold
            return
        Ls.append(L)
        bv = 0.0
        try:
            b = sim.boundary_spec.to_list[t][0]
            bv = float(getattr(b, "bloch_vec", 0.0) or 0.0)
        except Exception:
            bv = 0.0
        bvs.append(bv)
    try:
        eps = complex(sim.medium.eps_model(float(np.max(spec.freqs)))).real
    except Exception:
        eps = 1.0
    n = float(np.sqrt(max(eps, 1e-12)))
    for f in np.atleast_1d(np.asarray(spec.freqs, float)):
        k0 = 2 * np.pi * n * f / td.C_0 * 1e6       # 1/m, since grid.edges is in metres
        M = [int(k0 * L / (2 * np.pi)) + 1 for L in Ls]
        worst, wm = None, 1e9
        for m1 in range(-M[0], M[0] + 1):
            for m2 in range(-M[1], M[1] + 1):
                kt2 = (2 * np.pi * (m1 + bvs[0]) / Ls[0]) ** 2 + \
                      (2 * np.pi * (m2 + bvs[1]) / Ls[1]) ** 2
                r = abs(k0 ** 2 - kt2) / k0 ** 2
                if r < wm:
                    wm, worst = r, (m1, m2)
        if wm < tol:
            key = (spec.name, worst)
            if key in _RAYLEIGH_WARNED:
                continue
            _RAYLEIGH_WARNED.add(key)
            print(f"[openem] monitor {spec.name} covers a whole periodic face, and diffraction order "
                  f"{worst} sits on a Rayleigh cutoff (|kx^2|/k0^2 = {wm:.2e} <= {tol}, "
                  f"period {Ls[0]*1e6*f/td.C_0:.4f} x {Ls[1]*1e6*f/td.C_0:.4f} wavelengths). "
                  "That order carries power proportional to 1/kx, the continuous solution diverges "
                  "there, and **the total flux through the whole face does not converge with the "
                  "grid**; the same is true of the reference solver. For a convergent quantity, use "
                  "the power in a diffraction order instead, or move the period off the cutoff. "
                  "(OPENEM_RAYLEIGH_WARN=0 disables this message.)")


def _plane_flux_data(sim, sc, mon, ours, spec):
    """Plane flux: divided by the analytic incident power when there is a plane wave source, and by
    the source spectrum otherwise.

    A scene where the reference solver did no source normalization has to have its source spectrum
    multiplied back; see :func:`_denorm_spectrum`.
    """
    fs = np.asarray(spec.freqs, float)
    # monitor_flux always returns +S_axis (see the comment on flux.monitor_flux), so normal_dir has
    # to be multiplied in here: Tidy3D's FluxMonitor reports the power flowing out on the normal_dir
    # side. The closed-box path (_box_flux_data) already multiplied per face; the plane path did not.
    # Across the whole comparison set, all three plane flux monitors with normal_dir="-" had their
    # entire spectrum sign-flipped, and two of them came out with an error of exactly 2.00.
    try:
        _rayleigh_check(sim, sc, spec)
    except Exception as _e:            # diagnostics only; never affects the numbers
        pass
    vals = flux_mod.monitor_flux(ours[f"mon:{mon.name}"], sc.grid, spec) * float(spec.normal_dir)
    if getattr(sc, "sources", None):
        vals = vals / normalize.plane_flux_divisor(sc, spec, fs)
    else:
        vals = vals / normalize.flux_norm(sc, fs)
    s = _denorm_spectrum(sim, fs)
    if s is not None:
        vals = vals * np.abs(s) ** 2
    return td.FluxData(monitor=mon,
                       flux=FluxDataArray(vals, coords=dict(f=fs)))


def _norm_source_time(sim, sc):
    """Pick the matching Tidy3D source under the same priority as normalize.source_table (TFSF, then
    plane wave, then mode source or Gaussian beam, then current source) and take its source_time;
    fall back to the normalize_index source when nothing matches.
    """
    order = []
    if getattr(sc, "tfsf_sources", None):
        order.append((td.TFSF,))
    if getattr(sc, "sources", None):
        order.append((td.PlaneWave,))
    if getattr(sc, "mode_sources", None):
        order.append((td.ModeSource, td.GaussianBeam, td.PlaneWave))
    if getattr(sc, "dipoles", None):
        order.append((td.PointDipole, td.UniformCurrentSource, td.CustomCurrentSource, td.CustomFieldSource))
    for types in order:
        for src in sim.sources:
            if isinstance(src, types):
                return src.source_time
    return sim.sources[sim.normalize_index or 0].source_time


def _time_scale(sim, sc) -> float:
    """The time-domain conversion factor lambda, where E_ours(t) = lambda * E_td3d(t); see
    normalize.time_scale.
    """
    return normalize.time_scale(sc, _norm_source_time(sim, sc), sim.dt, sim.num_time_steps)


def _flux_time_data(sim, sc, mon, ours):
    """Instantaneous time-domain flux: the raw value divided by lambda^2 * UM^2.

    Time-domain quantities do not go through the frequency-domain source-spectrum normalization, but
    both E and H carry a factor of lambda (normalize.time_scale), and our area element dS is in SI
    m^2 while the reference uses um^2. The time axis is the sample step times dt. Returns None when
    the scene has no matching spec.
    """
    spec = {m.name: m for m in
            (getattr(sc, "flux_time_monitors", []) or [])}.get(mon.name)
    if spec is None:
        return None
    from openem.units import UM as _UM
    lam = _time_scale(sim, sc)
    arr = np.asarray(ours[f"ftm:{mon.name}"], float) / (lam * lam * _UM * _UM)
    ts = (spec.step_begin + np.arange(arr.size) * spec.interval) * sc.dt
    return td.FluxTimeData(monitor=mon,
                           flux=FluxTimeDataArray(arr, coords=dict(t=ts)))


def _refuse_nonfinite_ours(ours) -> None:
    """Refuse NaN and Inf in the solve output at the **assembly entry point**, failing closed.

    The solver exit carries a sentinel that only warns (solver._warn_nonfinite), but when a notebook
    runs inside a job the cell's stdout is lost along with a failed notebook conversion, which makes
    that print useless. A NaN in some batched mode amplitudes once passed through np.angle,
    np.unwrap and np.interp exactly that way and only blew up four cells later as a ValidationError
    on td.Box(size=nan), with no way to identify the source. Assembly and scene building were then
    shown to be innocent (with synthetic finite phasors, the amplitudes of all 11 swept simulations
    were finite, and a full array scan of the Scene was finite too), so a non-finite value can only
    come from a divergence in the GPU time stepping. This names the bad monitor on the spot and
    stops. To assemble with NaN deliberately, for instance just to look at a field plot, set
    OPENEM_NB_ALLOW_NONFINITE=1 and it degrades to a printed warning.

    Cost: one summation probe per array, since NaN and Inf propagate into a sum, the same as the
    solver sentinel. On healthy data the behaviour is bitwise unchanged.
    """
    keys = getattr(ours, "files", None)
    if keys is None:
        keys = list(ours.keys())
    bad = []
    for k in keys:
        if not str(k).startswith(("mon:", "fld:", "tim:", "ftm:", "tdf:")):
            continue
        a = np.asarray(ours[k])
        if a.dtype.kind not in "fc" or not a.size:
            continue
        s = complex(np.sum(a, dtype=np.complex128))
        if not (np.isfinite(s.real) and np.isfinite(s.imag)):
            n = int(np.count_nonzero(~np.isfinite(a)))
            bad.append(f"{k} ({n}/{a.size} non-finite values)")
    if not bad:
        return
    msg = ("solve output contains NaN/Inf: " + ", ".join(bad) + ".\n"
           "The field most likely diverged during GPU time stepping; assembly and scene building "
           "have been ruled out.\n"
           "This fails closed at the assembly entry point so that a NaN does not propagate through "
           "the notebook post-processing and blow up far from its source (amps, np.unwrap, "
           "np.interp, td.Box(size=nan)).\n"
           "To assemble with NaN deliberately, set OPENEM_NB_ALLOW_NONFINITE=1.")
    if os.environ.get("OPENEM_NB_ALLOW_NONFINITE"):
        print("[openem] WARNING (allowed by OPENEM_NB_ALLOW_NONFINITE=1): " + msg,
              flush=True)
        return
    # Raise DivergenceError (a RuntimeError subclass) so that in-job assembly writes
    # simdata.DIVERGENCE and the client goes straight into the staircasing and Courant retries. One
    # case used to write simdata.FAIL instead and reported only "scene build failed".
    raise DivergenceError(msg)


def build(sim, sc, ours, provenance=None):
    """Assemble a complete SimulationData. ``ours`` is the solve output from np.load(...)."""
    try:
        _refuse_nonfinite_ours(ours)
        _refuse_energy_growth(ours)
    except DivergenceError as _e:
        raise _e.attach(sim, sc, ours)
    fmon = {m.name: m for m in sc.field_monitors}
    xmon = {m.name: m for m in sc.flux_monitors}
    data = []
    # The denominator for diffraction efficiency is the **incident power**.
    #
    # It used to be "the sum of the flux over all diffraction surfaces", i.e. energy conservation for
    # a lossless grating, which only holds when **diffraction monitors sit on both sides**. With one
    # on the transmission side only, sum |amps|^2 is **identically 1**, an identity manufactured by
    # construction: one case measured T = 0.99999999823 against a transfer-matrix truth of 0.786.
    # Another passed only because it had monitors on both the R and T sides.
    #
    # For a plane wave scene the incident power has an analytic value (the same divisor used
    # elsewhere), so use it. Only without a plane wave source does it fall back to energy
    # conservation.
    _diff_pinc: dict = {}
    _diff_names = [m.name for m in sim.monitors
                   if type(m).__name__ == "DiffractionMonitor"
                   and f"mon:{m.name}" in ours]
    if _diff_names:
        from openem import flux as _flux_mod
        _xm = {m.name: m for m in (getattr(sc, "flux_monitors", None) or [])}
        _cons = 0.0
        for _nm in _diff_names:
            # monitor_flux must be used, since it carries the monitor's own normal axis and
            # transverse extent; plane_flux cannot be called bare, because it defaults to axis=2 and
            # a diffraction monitor with an x or y normal would get the area element of the wrong
            # axis. The stored phasor is (4, nf, n1, n2) where n1 and n2 are that monitor's two
            # transverse axes, while the area element would be computed as (nx, ny) for a z normal,
            # and the broadcast fails outright (measured, (1,58,57) against (1,81,58)).
            _fm = _xm.get(_nm)
            _flx = (_flux_mod.monitor_flux(np.asarray(ours[f"mon:{_nm}"]),
                                           sc.grid, _fm)
                    if _fm is not None
                    else _flux_mod.plane_flux(np.asarray(ours[f"mon:{_nm}"]),
                                              sc.grid))
            _cons += abs(float(np.real(np.atleast_1d(_flx)[0])))
        for _nm in _diff_names:
            _spec = _xm.get(_nm)
            if getattr(sc, "sources", None) and _spec is not None:
                _fs = np.asarray(_spec.freqs, float)
                _diff_pinc[_nm] = np.asarray(
                    normalize.plane_flux_divisor(sc, _spec, _fs), float)
            else:
                # With no plane wave source (a point dipole, a mode source and so on), tidy3d's amps
                # are absolute quantities and "the sum of the flux over the orders" must not be the
                # divisor: that would force sum |amps|^2 to 1.
                # The divisor is UM^2 * |field_norm|^2. The first factor converts the area element
                # from m^2 to um^2 (our power is 1e12 smaller than tidy3d's convention, the same
                # root cause as the FluxTimeMonitor bug), and the second is the source-spectrum
                # normalization: _diffraction_data only divided out the phase of field_norm, whose
                # modulus is 1, and the magnitude used to be absorbed by the efficiency
                # normalization where it could not be seen.
                from openem.units import UM as _UM
                _mo = next((m for m in sim.monitors if m.name == _nm), None)
                _fs2 = np.atleast_1d(np.asarray(
                    getattr(_mo, "freqs", [0.0]), float))
                _fn = np.abs(np.asarray(normalize.field_norm(sc, _fs2)))
                _diff_pinc[_nm] = np.asarray(
                    float(_UM * _UM) * (_fn ** 2), float)

    # A 3D FluxMonitor is computed as six faces, and the readout side adds them back under their
    # outward normals
    _boxflux = {b.name: b for b in (getattr(sc, "box_flux_monitors", None) or [])}
    _faces = {m.name: m for m in (getattr(sc, "flux_monitors", None) or [])}

    for mon in sim.monitors:
        nm = mon.name
        tname = type(mon).__name__
        if tname == "DiffractionMonitor":
            dd = _diffraction_data(sim, sc, mon, ours,
                                   _diff_pinc.get(mon.name, 0.0))
            if dd is not None:
                data.append(dd)
        elif nm in _boxflux and f"mon:{_boxflux[nm].face_names[0]}" in ours:
            data.append(_box_flux_data(sim, sc, mon, ours, _boxflux[nm], _faces))
        elif f"mon:{nm}" in ours and nm in xmon:
            data.append(_plane_flux_data(sim, sc, mon, ours, xmon[nm]))
        elif f"tdf:{nm}|Ex" in ours and nm in fmon:
            # Already in Tidy3D convention, so copy directly; the coordinates come from tdc: and
            # match the reference point for point
            fs = np.asarray(fmon[nm].freqs, float)
            data.append(_field_data_td(sim, mon, ours, nm, fs))
        elif f"fld:{nm}" in ours and nm in fmon:
            m = fmon[nm]
            fs = np.asarray(m.freqs, float)
            _n = normalize.field_norm(sc, fs)
            _s = _denorm_spectrum(sim, fs)
            if _s is not None:
                _n = _n / _s          # dividing it out of the divisor multiplies it into the data
            data.append(_field_data(sim, sc, mon, ours[f"fld:{nm}"], fs, _n))
        elif tname == "ModeMonitor":
            md = _mode_data(sim, sc, mon, ours)
            if md is not None:
                data.append(md)
        elif tname == "PermittivityMonitor":
            data.append(_perm_data(sim, sc, mon))
        elif tname == "ModeSolverMonitor":
            msd = _mode_solver_monitor_data(sim, mon)
            if msd is not None:
                data.append(msd)
        elif tname == "FieldTimeMonitor" and f"tim:{nm}" in ours:
            data.append(_fieldtime_data(sim, sc, mon, ours))
        elif tname in ("FieldProjectionKSpaceMonitor", "FieldProjectionAngleMonitor",
                       "FieldProjectionCartesianMonitor"):
            kd = _kspace_data(sim, sc, mon, ours)
            if kd is not None:
                data.append(kd)
        elif tname == "FluxTimeMonitor" and f"ftm:{nm}" in ours:
            ftd = _flux_time_data(sim, sc, mon, ours)
            if ftd is not None:
                data.append(ftd)
    # A minimal log: some notebooks parse sim_data.log (one reads "Field projection time (s):").
    # We have no corresponding real timing, so 0 is filled in to keep the parsing from breaking.
    _log = ("[openem] SimulationData assembled from local OpenEM solve\n"
            "Field projection time (s):    0.0\n"
            "Solver time (s):    0.0\n")
    sd = td.SimulationData(simulation=sim, data=tuple(data), log=_log)
    try:
        _check_divergence(sd)
    except DivergenceError as _e:
        raise _e.attach(sim, sc, ours)
    return sd


#: Upper bound on the magnitude of normalized monitor data. In Tidy3D's convention the field is in
#: V/um and a mode amplitude is in sqrt(W), and a healthy example peaks around 1e3. One diverging
#: scene reached 1e13, finite but enormous, which the NaN guard cannot catch.
DIVERGENCE_ABS_MAX = 1.0e6


class DivergenceError(RuntimeError):
    """The solve output looks divergent: NaN or Inf, finite but enormous, or the energy rebounded.

    The caller may retry with staircasing or a lower Courant number. When ``build`` raises, it
    attaches the assembly ingredients ``(sim, sc, ours)`` to the exception (:meth:`attach`), and once
    the fallbacks are exhausted :func:`failed_data` uses them to zero the monitor data and assemble
    as usual.
    """

    sim = None
    sc = None
    ours = None

    def attach(self, sim, sc, ours):
        if self.ours is None:
            self.sim, self.sc, self.ours = sim, sc, ours
        return self


#: The two thresholds of the energy-rebound guard: once the trough has fallen below GROWTH_FLOOR of
#: the peak, a rise back to GROWTH_REBOUND times that trough counts as divergence. In a passive
#: medium the total energy after the source stops only decays monotonically (residual PML reflection
#: is of order 1e-4), and dropping three orders of magnitude then climbing three back has no
#: physical explanation. EnergyDecay is the total of E and H, which does not oscillate at 2*omega,
#: so the sampling phase cannot manufacture a false rebound.
GROWTH_FLOOR = 1e-3
GROWTH_REBOUND = 1e3
#: Physical ceiling on the mode amplitude in a pure mode-source scene: a mode source injects 1 W, so
#: in a passive structure no mode can have |amp|^2 above 1. A factor of 10 margin is allowed, and N
#: coherently superposed mode sources relax it by N^2. One diverging member measured |amp|^2 = 1.5e7.
MODE_AMP2_MAX = 10.0


def _refuse_energy_growth(ours) -> None:
    """The second divergence guard at the assembly entry point: look at the time trace of the total
    energy and catch a divergence that is still finite but already exploding.

    The trace comes from ``ours['__meta']``, where the solver records the total energy every
    interval steps along with the step at which the source stops. Without a trace, i.e. with a
    pinned step count, early termination off, or an older result file, this is skipped. Two rules:

    A. **The total energy must not rise after the source stops**: exceeding GROWTH_REBOUND times the
       pre-stop peak counts as divergence. In a passive medium it only decays monotonically once the
       source stops (residual PML reflection is of order 1e-4), so a rise of three orders of
       magnitude can only be an instability in the time stepping.
    B. **Decay to a trough and then rebound**, the backstop for a source that never stops or an
       unknown source_end: the trough falls below GROWTH_FLOOR of the peak and then climbs back to
       GROWTH_REBOUND times the trough. EnergyDecay is the total of E and H, which does not
       oscillate at 2*omega, so the sampling phase cannot manufacture a false rebound.

    Measured on one optimization: two members grew by a factor of e every 50 steps or so, taking over
    before the source had even stopped and climbing monotonically to 1e26 (rule A); another was
    slower, ran its full 9,205 steps without terminating early with |amp| only reaching a finite
    3.9e3, invisible to both the NaN guard and the 1e6 ceiling, yet had its objective written as
    -69 dB and handed to the optimizer; a third was slower still, ending with a total energy eight
    orders of magnitude above the source peak while its mode amplitudes still looked normal, which a
    few thousand more steps would have turned to garbage (rules A and B both catch it).
    OPENEM_NB_ALLOW_NONFINITE=1 degrades this to a warning.
    """
    keys = getattr(ours, "files", None)
    if keys is None:
        keys = list(ours.keys())
    if "__meta" not in keys:
        return
    try:
        import json as _json
        meta = _json.loads(np.asarray(ours["__meta"]).item())
    except Exception:
        return
    tr = meta.get("energy_trace") or []
    if len(tr) < 3:
        return
    steps = np.array([t[0] for t in tr], dtype=np.int64)
    e = np.array([t[1] for t in tr], dtype=np.float64)
    if not np.all(np.isfinite(e)) or float(e.max()) <= 0.0:
        return                                   # NaN and Inf belong to _refuse_nonfinite_ours
    why = None
    src_end = int(meta.get("source_end", -1))
    if src_end >= 0:
        before = steps <= src_end
        # At least three samples before the source stops are needed to speak of a peak: with only the
        # one or two points at step 0, the energy is essentially zero (measured 3.8e-33 on one scene)
        # and any later energy would read as an explosion. That case is left to rule B.
        if int(before.sum()) >= 3 and (~before).any():
            e_src = float(e[before].max())
            hit = np.flatnonzero(~before & (e >= GROWTH_REBOUND * e_src))
            if e_src > 0 and hit.size:
                j = int(hit[0])
                why = (f"the source finished injecting at step {src_end} with a prior total-energy "
                       f"peak of {e_src:.2e}; by step {steps[j]} it had reached "
                       f"{e[j]:.2e} ({e[j] / e_src:.1e} times), ending at {e[-1]:.2e} at step {steps[-1]}")
    if why is None:
        r = e / np.maximum.accumulate(e)          # relative to the running peak
        prev_min = np.concatenate(([r[0]], np.minimum.accumulate(r)[:-1]))   # trough up to the previous sample
        rebound = r / np.maximum(prev_min, 1e-300)
        hit = np.flatnonzero((prev_min <= GROWTH_FLOOR) & (rebound >= GROWTH_REBOUND))
        if hit.size:
            j = int(hit[0])
            why = (f"at step {steps[j]} energy over running peak = {r[j]:.2e}, with a prior trough of "
                   f"{prev_min[j]:.2e} (a rebound of {rebound[j]:.1e} times), ending at "
                   f"{r[-1]:.2e} at step {steps[-1]}")
    if why is None:
        return
    msg = ("the solve output looks divergent: the total energy exploded after the source stopped. " +
           why +
           f"; early stop {meta.get('stop') or 'none'}.\n"
           "In a passive medium the total energy only decays monotonically once the source stops, so "
           "a rebound can only be an instability in the time stepping; on one scene it was a "
           "cell-scale mode on the top face of the output waveguide, in the outermost cell of the "
           "+x PML.\n"
           "To assemble with it deliberately, set OPENEM_NB_ALLOW_NONFINITE=1.")
    if os.environ.get("OPENEM_NB_ALLOW_NONFINITE"):
        print("[openem] WARNING (allowed by OPENEM_NB_ALLOW_NONFINITE=1): " + msg, flush=True)
        return
    raise DivergenceError(msg)


def failed_data(exc, tag=""):
    """What happens after all three fallbacks (as given, staircased, staircased plus Courant 0.6)
    have diverged: stop raising and **mark the solve failed**.

    Every monitor array is zeroed and the SimulationData is assembled as usual, so mode amplitudes
    are 0, transmission is 0 and the objective becomes +inf or the worst possible value. The
    optimizer then discards that design on its own rather than taking an exploded value seriously;
    one -69 dB value had stood as the "best" for nine generations.

    Set OPENEM_NB_DIVERGED_RAISE=1 to restore raising. The in-job assembly path
    (simdata.DIVERGENCE) has no ours in hand and still raises.
    """
    ours = getattr(exc, "ours", None)
    if os.environ.get("OPENEM_NB_DIVERGED_RAISE") == "1" or ours is None:
        raise exc
    keys = getattr(ours, "files", None)
    if keys is None:
        keys = list(ours.keys())
    zeros = {}
    for k in keys:
        v = ours[k]
        zeros[k] = (np.zeros_like(np.asarray(v))
                    if str(k).startswith(("mon:", "fld:", "tim:", "ftm:", "tdf:")) else v)
    print(f"[openem] {tag}: all three fallbacks are exhausted (as given, staircased and Courant 0.6 "
          "all diverged). Marking it failed: the monitor data is zeroed and assembled as usual, so "
          "transmission and flux come out 0 and the optimizer discards this design. "
          "OPENEM_NB_DIVERGED_RAISE=1 restores raising.\n"
          f"    the last verdict: {str(exc).splitlines()[0][:160]}", flush=True)
    _prev = os.environ.get("OPENEM_NB_ALLOW_NONFINITE")
    os.environ["OPENEM_NB_ALLOW_NONFINITE"] = "1"
    try:
        return build(exc.sim, exc.sc, zeros)
    finally:
        if _prev is None:
            os.environ.pop("OPENEM_NB_ALLOW_NONFINITE", None)
        else:
            os.environ["OPENEM_NB_ALLOW_NONFINITE"] = _prev


def _check_divergence(sd) -> None:
    if os.environ.get("OPENEM_NB_ALLOW_NONFINITE"):
        return
    bad = []
    # In a current-source scene (PointDipole, UniformCurrentSource, CustomCurrentSource), the
    # normalized fields and mode amplitudes are of order 1e7 to 1e9 in Tidy3D's convention to begin
    # with; two such cases measured 2.07e7 and 2.53e7, identical across three reruns, and were not
    # diverging at all. Such scenes are checked for NaN and Inf only, with no 1e6 ceiling.
    _srcs = getattr(getattr(sd, "simulation", None), "sources", None) or []
    current_like = any(type(s_).__name__ in ("PointDipole", "UniformCurrentSource", "CustomCurrentSource") for s_ in _srcs)
    limit = np.inf if current_like else DIVERGENCE_ABS_MAX
    # A pure mode-source scene injects 1 W, so in a passive structure |amp|^2 of any mode is at most
    # 1, with a factor of 10 margin and N sources relaxing it by N^2. One diverging member had
    # |amp|^2 = 1.5e7, finite and below the 1e6 amplitude ceiling, yet wrote its objective as -69 dB.
    _mode_only = bool(_srcs) and all(type(s_).__name__ == "ModeSource" for s_ in _srcs)
    amp2_limit = MODE_AMP2_MAX * len(_srcs) ** 2 if _mode_only else np.inf
    for d in sd.data:
        if type(d).__name__ in ("FieldTimeData", "FluxTimeData"):
            # Time-domain data is a raw quantity, not normalized by the source spectrum, and in
            # Tidy3D's convention the E(t) of a point dipole is already of order 1e9 (one example
            # peaked at 2.5e9 and was once misjudged as divergent), so the 1e6 ceiling does not apply
            continue
        for comp in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "amps", "flux"):
            a = getattr(d, comp, None) if comp != "flux" or type(d).__name__.startswith("Flux") else None
            if a is None:
                continue
            try:
                v = np.asarray(a.values)
            except Exception:
                continue
            if v.size and not np.all(np.isfinite(v)):
                bad.append(f"{d.monitor.name}.{comp} (contains NaN/Inf)")
            elif comp == "amps" and v.size and float(np.nanmax(np.abs(v))) ** 2 > amp2_limit:
                bad.append(f"{d.monitor.name}.amps (|amp|^2 max {float(np.nanmax(np.abs(v))) ** 2:.2e} "
                           f"> {amp2_limit:g}, impossible in a passive structure fed 1 W by a mode source)")
            elif v.size and np.nanmax(np.abs(v)) > limit:
                bad.append(f"{d.monitor.name}.{comp} (max {np.nanmax(np.abs(v)):.2e})")
    if bad:
        raise DivergenceError(
            "the solve output looks divergent: after normalization |value| > 1e6, or a mode "
            "amplitude exceeds its physical ceiling: " + ", ".join(bad[:4]) + ".\n"
            "Measured on one scene: with subpixel on, the time stepping diverges on the mixed "
            "interface cells of a dispersive medium, while subpixel=False or OPENEM_DISP_STAIRCASE=1 "
            "(staircasing the dispersive interface cells) matches the reference.")


def _mode_data(sim, sc, mon, ours):
    """ModeData: decompose the mode-plane field phasors into amps and n_complex."""
    try:
        from openem.scene import modes
    except Exception:
        return None
    mspec = {m.name: m for m in sc.mode_monitors}.get(mon.name)
    if mspec is None or f"fld:{mspec.plane_name}" not in ours:
        return None
    fs = np.atleast_1d(np.asarray(mon.freqs, float))
    # The mode-plane field must use **the same** normalization as FieldData and flux, i.e. the
    # source-type field_norm with the adjoint coefficient or source spectrum multiplied back, not
    # tidy3d's normalization of 1/spectrum. The two are equal only in a mode-source scene, where
    # field_norm = D. In a dipole scene field_norm = D * 1e6 (the field is in V/um), and with
    # 1/spectrum the mode-plane field would be in V/m, making |amps|^2 about 1e12 larger than the
    # flux; one coupling efficiency then came out as 1.4e11 against a reference value of 0.97.
    # The k-space projection path (_kspace_data) has passed norm this way for some time.
    _n = np.asarray(normalize.field_norm(sc, fs), dtype=np.complex128)
    _s = _denorm_spectrum(sim, fs)
    if _s is not None:
        _n = _n / _s
    # A mode plane builds its FieldData through **the same** path as an ordinary FieldMonitor.
    # It used to go through projection.field_data, whose colocation points are "primal boundaries
    # within the monitor bounds, inset by COORD_ATOL at each end, with the normal snapped to the
    # floor boundary". Those are not the same points as the ModeSolver plane grid (the colocation
    # boundaries of grid_expanded), so outer_dot interpolated a second time inside the colocation
    # path: H, at cell centers, took one more half-cell average along the normal than E, on cell
    # lines, and the net effect was an H/E ratio short by cos^2(beta*delta/2). A purely forward mode
    # therefore leaked into a_minus:
    #     |a_minus/a_plus| = (1 - cos^2)/(1 + cos^2) = (beta*delta)^2 / 4
    # For a silicon waveguide at 25 cells per wavelength that is exactly -48 dB, the floor seen in
    # one comparison that depended on neither the PML layer count nor run_time (measured
    # c_eff = 0.992354 against cos^2 = 0.992290).
    # _field_data's target coordinates are tidy3d's own colocation boundaries, the same grid as the
    # mode basis, so outer_dot no longer interpolates twice. Measured on a straight-waveguide probe
    # at 10/20/25/30/40 cells per wavelength, the backward leakage went from
    # -32.4/-54.1/-48.3/-61.3/-56.8 to -72.3/-70.1/-87.0/-81.3/-82.9 dB, and |a_plus| went from
    # 1.0119/1.0041/1.0016/1.0017/1.0003 to 0.9996 through 0.9998 across the board.
    fd = _field_data(sim, sc, modes.plane_monitor(mon),
                     ours[f"fld:{mspec.plane_name}"], fs, _n)
    amp = modes.amplitudes(sim, mon, fd)   # a dict of '+'/'-', each (nmode, nf)
    # ``amp['+']`` is (nmode, nf). coords are given as (direction, mode_index, f), and
    # ``ModeAmpsDataArray`` converts them into its own canonical order (direction, f, mode_index).
    # Remember that when reading its values: slicing directly by (mode, f) slices the frequency axis.
    nmode = amp["+"].shape[0]
    arr = np.stack([amp["+"], amp["-"]], axis=0)
    amps = ModeAmpsDataArray(
        arr, coords=dict(direction=["+", "-"],
                         mode_index=np.arange(nmode), f=fs))
    # **The effective index has to be the real one.** This used to be hardcoded to 1.0, which the
    # validation set never exposed because the scoring only looked at amps, but tidy3d's autograd
    # uses n_complex to build the **mode adjoint source**, and a hardcoded value makes the
    # inverse-design gradient wrong. The true value is right there in the mode solve result.
    md_full = modes.mode_solver_data(sim, mon)
    n_cx = np.asarray(md_full.n_complex.transpose("mode_index", "f").data,
                      dtype=np.complex128)
    if n_cx.shape != (nmode, len(fs)):
        # **Do not use np.resize**: it tiles the data into the new shape and silently scrambles it,
        # which is how modes 1 and 2 used to receive wrong values. A shape mismatch means the
        # convention changed, so report it plainly.
        raise RuntimeError(
            f"{mon.name}: n_complex has shape {n_cx.shape}, which does not match "
            f"({nmode} modes, {len(fs)} frequency points)")
    n_complex = ModeIndexDataArray(
        n_cx, coords=dict(mode_index=np.arange(nmode), f=fs))
    return td.ModeData(monitor=mon, amps=amps, n_complex=n_complex)


def _eps_is_freq_independent(sim) -> bool:
    """See :func:`openem.scene.media.is_freq_independent`, where the implementation now lives, shared
    with nb/shapegrad.

    This name is kept for tests/test_fasteps.py and for older notebooks.
    """
    from openem.scene import media as _media
    return _media.is_freq_independent(sim)


def _mode_solver_monitor_data(sim, mon):
    """ModeSolverMonitor to ModeSolverData: solve the modes on the monitor plane with tidy3d's own
    local mode solver.

    OpenEM does not do mode solving; this only connects a notebook's incidental mode-solve need to
    Tidy3D's solver.
    """
    try:
        from tidy3d.plugins.mode import ModeSolver
        ms = ModeSolver(simulation=sim, plane=td.Box(center=mon.center, size=mon.size), mode_spec=mon.mode_spec,
                        freqs=list(np.atleast_1d(np.asarray(mon.freqs, dtype=float))), direction=mon.direction,
                        colocate=bool(getattr(mon, "colocate", True)))
        from openem.scene.modes import _memo_eps as _me
        with _me(ms.simulation):
            msd = ms.solve()
        try:
            return msd.updated_copy(monitor=mon)
        except Exception:
            return msd
    except Exception as e:
        print(f"[openem] the local mode solve for ModeSolverMonitor {mon.name} failed: {str(e)[:100]}")
        return None


def _perm_data(sim, sc, mon):
    """PermittivityData: **the complex eps(omega) at the monitor frequencies**, sampled onto the
    monitor grid with Yee staggering (nearest neighbour, since permittivity is piecewise constant).
    eps_xx uses the Ex staggering, eps_yy the Ey one and eps_zz the Ez one.

    **sc.eps_ex/ey/ez must not be used**: those are the eps_inf of the FDTD update, with the
    dispersive part carried by the polarization currents. The forward physics is still right, but
    Tidy3D's PermittivityMonitor stores the complex eps at each frequency, and autograd's
    **geometry gradient is proportional to delta_eps = eps_inside - eps_outside**.

    Measured on a boundary-gradient case with a Sellmeier silicon nitride whose eps_inf is 1:
    submitting eps_inf gave 1.0000 to 1.0266 against a true value near 4, so delta_eps collapsed
    from about 3 to about 0.027 and the gradient collapsed with it, a ratio of 2.5e5 against the
    reference, while the correlation coefficient was still 0.98, i.e. the shape was essentially
    right. A scale difference would not matter on its own, since Adam is invariant to per-coordinate
    scale, but Adam's eps defaults to 1e-8 and our gradient was only of order 1e-7, which breaks
    that invariance.
    """
    from tidy3d.components.data.data_array import ScalarFieldDataArray
    from openem.scene import media as _media
    grid = sim.discretize_monitor(mon)
    fs = np.atleast_1d(np.asarray(mon.freqs, float))
    comp_of = {"eps_xx": ("Ex", sc.eps_ex, sc.pec_ex),
               "eps_yy": ("Ey", sc.eps_ey, sc.pec_ey),
               "eps_zz": ("Ez", sc.eps_ez, sc.pec_ez)}
    _flat_eps = _eps_is_freq_independent(sim)
    out = {}
    for key, (comp, epsarr, pecarr) in comp_of.items():
        xc, yc, zc = _comp_coords(grid, comp)          # µm
        # The scene's Yee coordinates on each axis, in metres
        sax = []
        for ax in range(3):
            a = sc.grid.axes[ax]
            base = a.edges if YEE_ON_EDGE[comp][ax] else a.centers
            sax.append(np.asarray(base, float))
        # A genuine nearest-neighbour index. This used to call searchsorted directly, which returns
        # an **insertion position** (i+1 for a point between sax[i] and sax[i+1]) and is
        # systematically off by half a cell. eps is piecewise constant, so a half-cell offset moves
        # the whole interface by one cell, and the geometry gradient is very sensitive to that.
        idx = []
        for ax in range(3):
            tgt = np.atleast_1d((xc, yc, zc)[ax]) * 1e-6
            j = np.clip(np.searchsorted(sax[ax], tgt), 1, len(sax[ax]) - 1)
            near = np.where(np.abs(tgt - sax[ax][j - 1])
                            <= np.abs(sax[ax][j] - tgt), j - 1, j)
            idx.append(np.clip(near, 0, epsarr.shape[ax] - 1))
        # Take eps(omega) with subpixel averaging and the disk cache, the same source as Tidy3D's own
        # path. When the medium is frequency-independent it is computed once and broadcast;
        # computing it per frequency drags a multi-frequency example to a halt (one went from 41
        # minutes to over 175 without finishing). The shape must match the eps_inf array; a mismatch
        # means the grid convention changed, so report it plainly.
        _fs_calc = fs[:1] if _flat_eps else fs
        per_f = []
        for _f in _fs_calc:
            e_w = _media.epsilon_complex(sim, comp, float(_f))
            if e_w.shape != epsarr.shape:
                raise RuntimeError(
                    f"{comp}: eps(omega) has shape {e_w.shape}, which does not match eps_inf {epsarr.shape}")
            per_f.append(e_w[np.ix_(idx[0], idx[1], idx[2])])
        sub = np.stack(per_f, axis=-1)                # (nx, ny, nz, nf or 1)
        if _flat_eps and len(fs) > 1:
            sub = np.repeat(sub, len(fs), axis=-1)
        if pecarr is not None:
            # PEC mask cells: the eps array is filled with 1.0 (the kernel zeroes them through the
            # mask), and the readout has to restore tidy3d's convention, constants.pec_val = -1e8,
            # which the client always uses. pec_ex/ey/ez are **flat cell index lists** (from the eps
            # decision in media.py), so they are scattered back to 3D first.
            mask3 = np.zeros(epsarr.shape, bool)
            mask3.ravel()[np.asarray(pecarr, dtype=np.int64)] = True
            sub = np.where(mask3[np.ix_(idx[0], idx[1], idx[2])][..., None],
                           -1.0e8, sub)
        out[key] = ScalarFieldDataArray(
            sub.astype(complex), coords=dict(x=xc, y=yc, z=zc, f=fs))
    return td.PermittivityData(monitor=mon, grid_expanded=grid, **out)


def _diffraction_data(sim, sc, mon, ours, pinc):
    """DiffractionMonitor: FFT the colocated tangential flux phasors into the s and p amplitudes of
    each diffraction order.

    The (4, nf, n1, n2) phasors on the near-field plane (E_t1, E_t2, H_t1, H_t2, already colocated
    along the normal) are decomposed into plane wave orders by a 2D FFT along the periodic
    directions. Each order recovers its normal component from transversality, and car_2_sph_field
    converts to spherical Etheta and Ephi. The **phase** first has the phase of ``field_norm``
    divided out, putting it back into Tidy3D's convention; see the comments in the body, since the
    adjoint source amplitude is set directly by the complex value of amps.
    Normalization: |amps|^2 is that order's power divided by the incident power (``pinc``, a scalar
    or a per-frequency array), so each order's power becomes an efficiency. **pinc must be the real
    incident power**; using "the sum of the flux over all diffraction surfaces" is only correct with
    monitors on both sides (see the comment in build).
    Parseval and energy conservation have both been verified.
    """
    key = f"mon:{mon.name}"
    pinc = np.atleast_1d(np.asarray(pinc, dtype=np.float64))
    if key not in ours or float(pinc.max()) <= 0:
        return None
    try:
        from tidy3d.components.data.data_array import DiffractionDataArray
        ph = np.asarray(ours[key])                      # (4, nf, n1, n2) raw phasors
        _dbg(f"{mon.name} diffraction phasor shape {ph.shape} (4, nf, n_t1, n_t2)")
        # **Phase convention**: the raw phasors carry one extra arg(field_norm) relative to Tidy3D's
        # fields, namely the phase of the injection-table DFT. On the magnitude side it is later
        # rescaled into an efficiency by s and cannot be seen, but on the phase side it reaches amps,
        # and tidy3d's diffraction adjoint source amplitude is exactly
        # `1j * grad_const * amps` (source_factory._diffraction_plane_wave), so the whole gradient is
        # rotated by it. Measured, arg(D) = +pi/2, which together with that 1j makes -1, turning the
        # gradient into Re[i*Z] = -Im Z instead of Re Z.
        # Only the phase is divided out, whose modulus is 1: Szm has the form E * conj(H), so a
        # common unit phase cancels exactly and power, efficiency, T and R do not change by a bit.
        _fnorm = np.asarray(normalize.field_norm(
            sc, np.atleast_1d(np.asarray(mon.freqs, float))))
        _uph = _fnorm / np.abs(_fnorm)                  # (nf,)
        ph = ph / _uph[None, :, None, None]
        ax = int(np.argmin([abs(s) for s in mon.size])) # the normal axis, where size=0
        _xm_spec = next((m for m in (getattr(sc, "flux_monitors", None) or []) if m.name == mon.name), None)
        t1, t2 = grid_mod.transverse_axes(ax)
        UM_ = 1e-6
        # Refractive index of the medium the monitor sits in (in Tidy3D, DiffractionData.medium is
        # the medium at the monitor): the square root of the median eps on the monitor plane, the
        # same approach as normalize.background_index.
        n_mon = None
        try:
            _mm = sim.monitor_medium(mon)              # Tidy3D's own semantics: the medium at the monitor center
            _f0 = float(np.atleast_1d(np.asarray(mon.freqs, float))[0])
            n_mon = float(np.sqrt(complex(_mm.eps_model(_f0))).real)
        except Exception as _me:
            _dbg(f"{mon.name}: monitor_medium is unavailable ({type(_me).__name__}); falling back to "
                 "the median eps on the monitor plane")
        if os.environ.get("OPENEM_DIFF_MEDIUM") == "0":     # diagnostic knob: treat it as vacuum (the old behaviour)
            n_mon = 1.0
        if n_mon is None or not np.isfinite(n_mon) or n_mon <= 0:
            from openem.scene import sources as _srcmod
            _e = np.real(np.asarray(sc.eps_ex))
            _ki = int(_srcmod.nearest_edge_index(sc.grid.axes[ax].edges, float(mon.center[ax]) * UM_))
            _sl = [slice(None)] * 3
            _sl[ax] = max(0, min(_e.shape[ax] - 1, _ki))
            n_mon = float(np.sqrt(max(float(np.median(_e[tuple(_sl)])), 1e-12)))
        _dbg(f"{mon.name}: refractive index at the monitor n={n_mon:.4f} (k, allowed orders and "
             "angles are computed in the medium)")
        # The transverse extent follows **the same area rule as flux** (flux.axis_dl: a genuinely 2D
        # flat axis counts as FLAT_AXIS_M = 1 um, everything else its actual extent). The
        # normalization denominator pinc = plane_flux_divisor uses that same rule. This used to take
        # the actual grid extent, and on a 2D example the flat axis is only one cell of about
        # lambda/60, which squeezed every order's power down by dl/1um (a probe measured a 2D vacuum
        # plane wave with T flux 0.982 against a diffraction power sum of 0.025).
        from openem import flux as _flux_mod2
        Lt1 = float(np.sum(_flux_mod2.axis_dl(sc.grid.axes[t1])[0]))
        Lt2 = float(np.sum(_flux_mod2.axis_dl(sc.grid.axes[t2])[0]))
        if os.environ.get("OPENEM_DIFF_AREA") == "edges":    # diagnostic knob: the old behaviour, the actual grid extent
            _e1 = np.asarray(sc.grid.axes[t1].edges); _e2 = np.asarray(sc.grid.axes[t2].edges)
            Lt1 = float(_e1[-1] - _e1[0]); Lt2 = float(_e2[-1] - _e2[0])
        bl = sim.boundary_spec
        def _bv(a):
            b = getattr(bl, "xyz"[a]).minus
            return float(getattr(b, "bloch_vec", 0.0)) if type(b).__name__ == "BlochBoundary" else 0.0
        # A zero-size transverse dimension (a genuinely 2D sim) counts as a single order: dk=0 and
        # that axis has only order 0
        dk1 = 2*np.pi/Lt1 if Lt1 > 0 else 0.0
        dk2 = 2*np.pi/Lt2 if Lt2 > 0 else 0.0
        k1_0 = _bv(t1)*dk1; k2_0 = _bv(t2)*dk2
        fs = np.atleast_1d(np.asarray(mon.freqs, float))
        C0 = 2.99792458e8
        # Collect the allowed orders across frequencies, with the first frequency fixing the set;
        # every diffraction example in the validation set is single-frequency or shares one set
        Et1, Et2, Ht1, Ht2 = (ph[i] for i in range(4))   # (nf,n1,n2)
        n1, n2 = ph.shape[2], ph.shape[3]
        # The transverse grid may be **non-uniform**, since an auto grid refines inside a structure
        # footprint. Each order's amplitude must be a DFT over the samples' actual coordinates and
        # area elements, not an equispaced np.fft: order 0 is an area-weighted average, while an
        # equispaced FFT gives the arithmetic mean of the samples, so the refined region (which in
        # one case sat right around the dipole and dielectric block, where the field is strongest)
        # was counted about twice over: order 0 came out 6.07 W against a reference around 2.64. On a
        # uniform grid the expression below is bitwise identical to the old FFT.
        # The area elements follow the same rule as flux.area_weights: the E_t1 and H_t2 terms use
        # dl1 * dl2_dual, the E_t2 and H_t1 terms dl1_dual * dl2.
        def _pos_dl(axis_obj, n):
            dl, dld = _flux_mod2.axis_dl(axis_obj)
            if axis_obj.is_flat:
                return np.zeros(1), dl, dld
            e = np.asarray(axis_obj.edges, dtype=np.float64)
            return (e[:-1] - e[0])[:n], dl[:n], dld[:n]
        p1, dl1, dld1 = _pos_dl(sc.grid.axes[t1], n1)
        p2, dl2, dld2 = _pos_dl(sc.grid.axes[t2], n2)
        wA = np.outer(dl1, dld2) / (Lt1*Lt2)      # area weights for the E_t1 and H_t2 terms
        wB = np.outer(dld1, dl2) / (Lt1*Lt2)      # area weights for the E_t2 and H_t1 terms
        # The orders must be exact integers: for odd n, fftfreq(n)*n gives 0.99999..., which a later
        # cast to int truncates to three zeros. In one case an orders_y of 49 made sel(orders_y=0)
        # return a (3, nf) array and the notebook indexed out of range. rint first, then cast.
        m1 = np.rint(np.fft.fftshift(np.fft.fftfreq(n1, d=1.0))*n1).astype(int)
        m2 = np.rint(np.fft.fftshift(np.fft.fftfreq(n2, d=1.0))*n2).astype(int)
        k1 = k1_0 + m1*dk1; k2 = k2_0 + m2*dk2
        ndir = 1.0 if str(mon.normal_dir) == "+" else -1.0
        # Assemble frequency by frequency
        nf = len(fs)
        # The center frequency fixes the allowed orders first
        allowed_x = None
        acc = {}
        for fi, f0 in enumerate(fs):
            k0 = n_mon * 2*np.pi*f0/C0          # wavenumber in the medium
            # The order sequence matches the previous fftshift(fft2): order m corresponds to
            # k1[m] = k1_0 + m1[m] * dk1
            K1e = np.exp(-1j*np.outer(k1, p1))     # (n1, n1)
            K2e = np.exp(-1j*np.outer(k2, p2))     # (n2, n2)
            def fo(a, w):
                return K1e @ (a[fi]*w) @ K2e.T
            E1m, E2m, H1m, H2m = fo(Et1, wA), fo(Et2, wB), fo(Ht1, wB), fo(Ht2, wA)
            K1, K2 = np.meshgrid(k1, k2, indexing="ij")
            kz2 = k0**2 - K1**2 - K2**2
            prop = kz2 > 0
            kz = ndir*np.sqrt(np.abs(kz2))
            if os.environ.get("OPENEM_DIFF_SPLIT", "1") != "0":
                # Keep only the branch propagating along normal_dir: in a total-field region where
                # incident and reflected waves superpose, the net Poynting vector gives 1-R, not R
                E1m, E2m, H1m, H2m = _split_forward(E1m, E2m, H1m, H2m, K1, K2, kz2, f0, n_mon,
                                                    ndir, handed=(-1.0 if ax == 1 else 1.0))
            # Move the phase reference to the monitor center. The phasors are recorded on the snapped
            # grid line x_e, while tidy3d's amps take their phase from the monitor center x_c that
            # the user gave (the adjoint PlaneWave it builds for a diffraction objective is also
            # referenced to x_c, and our beam carrier is evaluated from the analytic profile at x_e
            # but likewise referenced to x_c). Without the move the whole adjoint gradient rotates by
            # k * (x_c - x_e): one case with delta = -0.030 um gave -7.1 degrees, and a
            # complex-amplitude probe measured the adjoint against the difference at -7.9 degrees.
            # That objective is the small real part of a large complex number, and 7 degrees flips
            # its sign. Each propagating order is translated by x_c - x_e under its own kz, signed by
            # normal_dir; evanescent orders are left alone. Power and efficiency are unaffected.
            _x_e = float(sc.grid.axes[ax].edges[int(_xm_spec.plane_index)]) if _xm_spec is not None else None
            if _x_e is not None and os.environ.get("OPENEM_DIFF_CENTER_PHASE", "1") != "0":
                _shift = np.where(prop, np.exp(1j * kz * (float(mon.center[ax]) * UM_ - _x_e)), 1.0)
                E1m, E2m, H1m, H2m = E1m * _shift, E2m * _shift, H1m * _shift, H2m * _shift
            with np.errstate(divide="ignore", invalid="ignore"):
                Ezm = -(K1*E1m + K2*E2m)/kz
                Hzm = -(K1*H1m + K2*H2m)/kz
            # Power of each order along the normal (the normal Poynting component), times the whole face
            Szm = 0.5*np.real(E1m*np.conj(H2m) - E2m*np.conj(H1m)) * (Lt1*Lt2)
            if allowed_x is None:
                xm = prop.any(1); allowed_x = np.where(xm)[0]
            acc[fi] = (E1m, E2m, Ezm, H1m, H2m, Hzm, Szm, prop, K1, K2, k0)
        ox = m1[allowed_x].astype(int); oy = m2.astype(int)
        nox, noy = len(ox), len(oy)
        Et = np.zeros((nox, noy, nf), complex); Ep = np.zeros_like(Et)
        Ht = np.zeros_like(Et); Hp = np.zeros_like(Et)
        for fi in range(nf):
            E1m, E2m, Ezm, H1m, H2m, Hzm, Szm, prop, K1, K2, k0 = acc[fi]
            for ii, ix in enumerate(allowed_x):
                for jj in range(noy):
                    if not prop[ix, jj]:
                        continue
                    th, phi = mon.kspace_2_sph(K1[ix,jj]/k0, K2[ix,jj]/k0, axis=2)
                    # (theta, phi) are the **local** angles of an axis=2 frame, which is how tidy3d's
                    # DiffractionData.compute_angles treats any normal: theta from the normal, phi
                    # from the first tangential axis t1 towards t2. So the components are fed to
                    # car_2_sph_field in the local frame (t1, t2, normal) too. They used to be laid
                    # out in the global (x, y, z) frame, which happens to agree for a z normal but
                    # gets the s and p decomposition entirely wrong for an x or y normal: at normal
                    # incidence p should correspond to E_t1 but received the normal component, which
                    # is about 0. The per-order power still came out right through the eff rescaling
                    # below, but the s and p shares and phases were wrong, and since tidy3d builds
                    # the diffraction adjoint source from amps, the gradient came out near zero.
                    # Per tidy3d's documentation, at normal incidence an x normal has P and S
                    # corresponding to Ey and Ez, a y normal to Ex and Ez, and a z normal to Ex and
                    # Ey, i.e. P to E_t1 and S to E_t2. In the local frame at theta=0, Etheta = f_x =
                    # E_t1 and Ephi = f_y = E_t2, exactly as required.
                    _, et, ep = mon.car_2_sph_field(E1m[ix,jj], E2m[ix,jj], Ezm[ix,jj], th, phi)
                    _, ht, hp = mon.car_2_sph_field(H1m[ix,jj], H2m[ix,jj], Hzm[ix,jj], th, phi)
                    # tidy3d's DiffractionData.amps uses the wave impedance of the medium at the
                    # monitor, eta = eta0/n (monitor_data.py: eta = ETA_0/sqrt(eps)), with
                    # power = |amps|^2. This used to hardcode vacuum eta0, so in a medium the power
                    # read out as the true value times n; in one case with n=3.5 the diffraction R
                    # came out 3.5 times too large.
                    eta = 376.730313 / n_mon
                    norm = 1.0/np.sqrt(2.0*eta)/np.sqrt(max(np.cos(th), 1e-9))
                    mag2 = abs(et)**2 + abs(ep)**2
                    _pi = float(pinc[fi] if pinc.size > 1 else pinc[0])
                    eff = abs(Szm[ix,jj])/_pi
                    # scale makes |amps|^2 equal eff, preserving the s and p shares and the phase
                    s = np.sqrt(eff/(mag2*norm**2)) if mag2 > 0 else 0.0
                    Et[ii,jj,fi] = et*s; Ep[ii,jj,fi] = ep*s
                    Ht[ii,jj,fi] = ht*s; Hp[ii,jj,fi] = hp*s
        coords = dict(orders_x=ox, orders_y=oy, f=fs)
        da = lambda v: DiffractionDataArray(v, coords=coords)
        z0 = np.zeros((nox, noy, nf), complex)
        return td.DiffractionData(
            monitor=mon, Etheta=da(Et), Ephi=da(Ep), Er=da(z0),
            Htheta=da(Ht), Hphi=da(Hp), Hr=da(z0),
            sim_size=(Lt1/UM_, Lt2/UM_), bloch_vecs=(_bv(t1), _bv(t2)),
            medium=td.Medium(permittivity=n_mon ** 2))
    except Exception as e:
        print(f"[openem] assembling the diffraction data of {mon.name} failed: {type(e).__name__}: {str(e)[:80]}")
        return None


def _split_forward(E1, E2, H1, H2, K1, K2, kz2, f0, n_mon, ndir, handed=1.0):
    """Split the tangential phasors (E_t1, E_t2, H_t1, H_t2) of each diffraction order into the two
    branches propagating along +- the normal, returning only the ``ndir`` one.

    Tidy3D's convention is e^{i(k.r - omega t)} with k = k_t p_hat + s kz n_hat (s = +-1), where
    p_hat = k_t/|k_t| and s_hat = n_hat x p_hat:
    - s polarization: E = E_s s_hat, so H = (k x E)/(omega mu0) and the tangential
      H_t = -s (kz/omega mu0) E_s p_hat
    - p polarization: H = H_s s_hat, so E = -(k x H)/(omega eps) and the tangential
      E_t = +s (kz/omega eps) H_s p_hat
    Superposing the two branches:
      E.s_hat = E_s+ + E_s-,  H.p_hat = -(kz/omega mu0)(E_s+ - E_s-)
      H.s_hat = H_s+ + H_s-,  E.p_hat = (kz/omega eps)(H_s+ - H_s-)
    which inverts to
      E_s^+- = (E.s_hat -+ (omega mu0/kz) H.p_hat) / 2
      H_s^+- = (H.s_hat +- (omega eps/kz) E.p_hat) / 2
    ``handed`` is +1 when (t1, t2, n_hat) is right-handed; for ax=1, (x, z, y) is left-handed and
    s_hat flips sign. An evanescent order, where kz^2 <= 0, is returned unchanged.
    """
    kt = np.hypot(K1, K2)
    with np.errstate(divide="ignore", invalid="ignore"):
        p1 = np.where(kt > 0, K1 / np.where(kt > 0, kt, 1.0), 1.0)
        p2 = np.where(kt > 0, K2 / np.where(kt > 0, kt, 1.0), 0.0)
    s1, s2 = -handed * p2, handed * p1                 # ŝ = n̂ × p̂
    ok = kz2 > 0
    kz = np.sqrt(np.where(ok, kz2, 1.0))
    w = 2 * np.pi * f0
    mu0 = 4e-7 * np.pi
    eps = 8.8541878128e-12 * n_mon ** 2
    Zs = w * mu0 / kz                                  # ωμ₀/kz
    Yp = w * eps / kz                                  # ωε/kz
    Es = E1 * s1 + E2 * s2; Ep = E1 * p1 + E2 * p2
    Hs = H1 * s1 + H2 * s2; Hp = H1 * p1 + H2 * p2
    sg = 1.0 if ndir > 0 else -1.0
    Es_f = 0.5 * (Es - sg * Zs * Hp)
    Hs_f = 0.5 * (Hs + sg * Yp * Ep)
    Ep_f = sg * Hs_f / Yp                              # E_p = s (kz/ωε) H_s
    Hp_f = -sg * Es_f / Zs                             # H_p = −s (kz/ωμ₀) E_s
    E1f = Es_f * s1 + Ep_f * p1; E2f = Es_f * s2 + Ep_f * p2
    H1f = Hs_f * s1 + Hp_f * p1; H2f = Hs_f * s2 + Hp_f * p2
    return (np.where(ok, E1f, E1), np.where(ok, E2f, E2),
            np.where(ok, H1f, H1), np.where(ok, H2f, H2))


def _dbg(msg: str) -> None:
    """Print, and also append to openem_debug.log in the current directory.

    When a whole-notebook conversion fails, a cell's stdout is never written back into the notebook,
    and this file is the only diagnostic left to look at afterwards.
    """
    print(f"[openem] {msg}")
    try:
        with open(os.environ.get("OPENEM_DEBUG_LOG", "openem_debug.log"), "a") as f:
            import time as _t
            f.write(f"{_t.strftime('%H:%M:%S')} {msg}\n")
    except Exception:
        pass


def _kspace_data(sim, sc, mon, ours):
    """FieldProjectionKSpaceData: near-field surfaces to a near-to-far transform.

    ``surface_monitors`` enumerates the projection surfaces, ``field_data`` colocates the native Yee
    surface phasors into a td.FieldData, and ``project`` runs the FieldProjector. The surface field
    keys are ``fld:<monitor>__n2f_<axis><index>``, the solver's naming template.
    The absolute-magnitude convention carries the same reservation it always did, and is good enough
    for the figures a notebook produces.
    """
    try:
        from openem.scene import projection
    except Exception:
        return None
    specs = {m.name: m for m in sc.field_monitors}
    faces = projection.surface_monitors(mon, sim)
    datas = []
    for fm, _nd, _ax in faces:
        key = f"fld:{fm.name}"
        if key not in ours or fm.name not in specs:
            _dbg(f"{mon.name} projection: missing surface field {key} (in ours={key in ours}, "
                 f"in specs={fm.name in specs})")
            return None                      # a single missing surface makes the projection impossible
        try:
            # The projection input uses the same normalization as FieldData: raw / field_norm, which
            # depends on the source type (a dipole gives D * 1e6, a plane wave D * 1e6 *
            # sqrt(n*A/2*eta0), TFSF D/sqrt(2*eta0) times a phase, a mode source D). It used to use
            # normalization * FAR_FIELD_SCALE = (1/D) * 1e-6, which happens to be equivalent only for
            # a dipole; a plane wave scene was short by 1e6/sqrt(n*A/2*eta0) (a probe measured all
            # four kinds of projection monitor at 1.8e-7, and another case at 1/1.36e6), while the
            # near field was right all along.
            _fs = np.atleast_1d(np.asarray(mon.freqs, float))
            _inv = 1.0 / np.asarray(normalize.field_norm(sc, _fs), dtype=np.complex128)
            datas.append(projection.field_data(
                sim, fm, specs[fm.name], _read_face_with_retry(ours, key, mon.name, fm.name), norm=_inv))
        except Exception as e:
            _dbg(f"{mon.name}: colocating surface field {fm.name} failed: {type(e).__name__}: {str(e)[:300]}")
            return None
    try:
        return projection.project(sim, mon, faces, datas)
    except Exception as e:
        _dbg(f"{mon.name}: the projection failed: {type(e).__name__}: {str(e)[:300]}")
        return None


def _read_face_with_retry(ours, key, mon_name, face_name, tries: int = 6):
    """Fetch one surface-field array from an npz, retrying with backoff when the file is not fully
    written yet.

    These arrays are npz files a solve job wrote to shared storage, and a network filesystem can
    expose the directory entry before the contents land, so a read halfway through raises
    ``BadZipFile: Truncated file header`` or ``EOFError``. Those are **transient** errors, entirely
    unlike a genuine colocation failure (mismatched geometry, a missing component), and must not
    make a whole projection monitor disappear. Two full notebook reruns lost four projection
    monitors exactly that way and ended in a KeyError inside the notebook.
    """
    import time as _t
    import zipfile as _zf
    last = None
    for i in range(tries):
        try:
            return np.asarray(ours[key])
        except (_zf.BadZipFile, EOFError, OSError, ValueError) as e:
            last = e
            _dbg(f"{mon_name}: reading surface field {face_name} failed on attempt {i + 1} "
                 f"({type(e).__name__}); waiting {5 + 5 * i}s before retrying")
            _t.sleep(5 + 5 * i)
    raise RuntimeError(f"reading surface field {face_name} failed on all {tries} attempts: {type(last).__name__}: {last}")

def _downsample_inds(n: int, k: int) -> np.ndarray:
    """Point selection for Tidy3D's ``interval_space``: take every k-th point, **always keeping the
    last one**.

    The rule was reverse-engineered from reference data (a time movie monitor with
    interval_space=(3,3,1) and colocate=True): an x axis of 467 points gives [0,3,...,465] plus
    [466], 157 points in all; a y axis of 241 points gives [0,3,...,240], 81 points, where the last
    point is already in the stepped sequence and is not repeated.
    """
    n, k = int(n), max(int(k), 1)
    inds = np.arange(0, n, k)
    if n > 0 and inds[-1] != n - 1:
        inds = np.append(inds, n - 1)
    return inds


def _apply_interval_space(mon, xc, yc, zc):
    """Thin the coordinates of all three axes by the monitor's ``interval_space``; the buffer is at
    full resolution and the writeback takes a subset.
    """
    step = tuple(int(v) for v in getattr(mon, "interval_space", (1, 1, 1)))
    if step == (1, 1, 1):
        return xc, yc, zc
    return tuple(np.asarray(c)[_downsample_inds(len(c), k)]
                 for c, k in zip((xc, yc, zc), step))


def _fieldtime_data(sim, sc, mon, ours):
    """FieldTimeData: ``tim:<name>`` of shape (nc, nslots, ni, nj, nk) turned into a time-domain
    array per component.
    """
    from tidy3d.components.data.data_array import ScalarFieldTimeDataArray
    spec = {m.name: m for m in sc.field_time_monitors}.get(mon.name)
    if spec is None:
        return None
    arr = ours[f"tim:{mon.name}"]              # (nc, nslots, ni, nj, nk)
    nslots = arr.shape[1]
    lam = _time_scale(sim, sc)
    ts = (spec.step_begin + np.arange(nslots) * spec.interval) * sc.dt
    grid = sim.discretize_monitor(mon)
    from openem import colocate
    comps = {}
    names = ["Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]
    _colo = bool(getattr(mon, "colocate", False))
    for ci, cidx in enumerate(spec.comps):
        comp = names[cidx]
        xc, yc, zc = _apply_interval_space(
            mon, *_comp_coords(grid, comp, colocate=_colo))
        v = arr[ci]                            # (nslots, ni, nj, nk)
        # The time-domain field is the raw value divided by lambda (normalize.time_scale), and lambda
        # depends on the source type: 1e6 for a current source (one case had Re{Ey} at 4e9 against a
        # reference 4e3), 1 for a mode source, and 1e6 * sqrt(n*A/2*eta0) * sqrt(cos) * (table
        # amplitude) for a plane wave. This used to be a hardcoded /1e6.
        v = v / lam
        # Native Yee positions to the target coordinates. With colocate=True all six components share
        # one grid; with False the targets are a subset of the sources, where linear interpolation is
        # bitwise the same as taking the original value. The targets are clipped into the source
        # range, taking the nearest value outside the boundary.
        src = colocate.sample_coords(sc.grid, spec.origin, spec.box, comp)
        tgt = []
        for s_um, sa in zip((xc, yc, zc), src):
            s_m = np.asarray(s_um, float) * 1e-6
            tgt.append(np.clip(s_m, float(sa.min()), float(sa.max())))
        v = colocate.interp_to(v, src, tgt)               # (nslots, x, y, z)
        v = np.moveaxis(v, 0, -1)                          # (x, y, z, t)
        comps[comp] = ScalarFieldTimeDataArray(
            v.astype(float),
            coords=dict(x=xc, y=yc, z=zc, t=ts))
    if not comps:
        return None
    return td.FieldTimeData(monitor=mon, grid_expanded=grid, **comps)



# Helper: take a monitor's origin and box from sc (the field names may differ after folding)
def mon_origin(sc, mon):
    for m in sc.field_monitors:
        if m.name == mon.name:
            return m.origin
    return (0, 0, 0)


def mon_box(sc, mon):
    for m in sc.field_monitors:
        if m.name == mon.name:
            return m.box
    return (1, 1, 1)
