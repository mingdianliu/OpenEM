# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""How the autograd hook makes a traced Simulation static.

Autograd31 (GratingCouplerWithBeamOptimization) goes through the high-level wrapper in
``tidy3d.plugins.autograd.optimize``, so the ``run(sim)`` inside the objective receives a
Simulation **that still carries tracers**: GaussianBeam.angle_theta, for instance, is still an
ArrayBox. On the normal autograd path tidy3d calls ``to_static`` before ``_run_tidy3d``, but the
high-level wrapper and a direct ``web.run`` slip past it. The hook's ``_solve`` has to strip it to
a numeric sim right at the entry, otherwise ``float(ArrayBox)`` inside ``from_simulation`` raises
TypeError immediately.
"""

from __future__ import annotations

import types

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")
ag = pytest.importorskip("autograd")

# This module imports openem.solver at the top level, which imports cupy at its own top level: on a
# CPU-only machine that would crash during collection.
pytest.importorskip("cupy", reason="needs CUDA and cupy")

from openem import autograd_hook as H
from openem import scene as scene_mod

F0 = 2.0e14


def _traced_sim(p):
    """A GaussianBeam sim whose angle_theta carries a tracer."""
    src = td.GaussianBeam(
        center=(0, 0, 1.0), size=(td.inf, td.inf, 0), direction="-",
        source_time=td.GaussianPulse(freq0=F0, fwidth=F0 / 10),
        angle_theta=p * 0.1, pol_angle=np.pi / 2, waist_radius=1.5)
    box = td.Structure(geometry=td.Box(center=(0, 0, 0), size=(1.5, td.inf, 0.22)),
                       medium=td.Medium(permittivity=12.1))
    mon = td.ModeMonitor(center=(-1.5, 0, 0), size=(0, td.inf, 2), freqs=[F0],
                         mode_spec=td.ModeSpec(num_modes=1), name="m")
    return td.Simulation(
        size=(4, 0, 3),
        grid_spec=td.GridSpec.auto(wavelength=td.C_0 / F0, min_steps_per_wvl=8),
        structures=[box], sources=[src], monitors=[mon], run_time=1e-13,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.pml(), y=td.Boundary.periodic(), z=td.Boundary.pml()))


def test_solve_statically_fies_traced_sim(monkeypatch):
    """The hook's _solve entry must strip the traced sim clean, so from_simulation sees no
    tracers."""
    seen = {}

    def fake_from_sim(sim):
        seen["has_tracers"] = sim._has_tracers
        return types.SimpleNamespace(shape=(4, 4, 1))

    monkeypatch.setattr(scene_mod, "from_simulation", fake_from_sim)
    monkeypatch.setattr(H.scene_mod, "from_simulation", fake_from_sim)
    monkeypatch.setattr(H.solver, "run", lambda sc, **kw: types.SimpleNamespace(
        steps_run=1, phasors={}, field_phasors={}, time_samples={}, flux_time={}))
    from openem.nb import backend as nb
    monkeypatch.setattr(nb, "build", lambda sim, sc, arrs: "SD")
    monkeypatch.setitem(H._STATE, "kernels", object())
    monkeypatch.setitem(H._STATE, "max_cells", None)

    def f(p):
        H._solve(_traced_sim(p), "t")            # go through the real hook entry
        return p

    ag.grad(f)(1.0)                              # must not raise an ArrayBox TypeError
    assert seen.get("has_tracers") is False


def test_to_static_is_noop_for_plain_sim():
    """A sim with no tracers: to_static should return it unchanged, so the hook has zero effect on
    the normal path."""
    sim = _traced_sim(1.0)                        # called directly, so p is a plain float
    assert sim.to_static() is sim or sim.to_static()._has_tracers is False


def test_traced_web_run_routes_through_autograd_primitive(monkeypatch):
    """A single traced web.run has to be handed back to tidy3d's own autograd run (the primitive);
    otherwise the output carries no tracer and the objective breaks the gradient as soon as it
    reads from sim_data (the root cause in Autograd12/20)."""
    routed = {"called": False}

    def fake_autograd_run(sim, **kw):
        routed["called"] = True
        routed["sim"] = sim
        return "SD_traced"

    monkeypatch.setitem(H._STATE, "autograd_run", fake_autograd_run)
    sim = _traced_sim(1.0)  # built from a plain float, so it carries no tracer

    # No tracer: must not be routed, goes to the local _solve
    plain_solved = {"n": 0}
    monkeypatch.setattr(H, "_solve", lambda s, tag, **kw: plain_solved.__setitem__("n", plain_solved["n"] + 1) or "SD_plain")
    assert H._run_plain(sim, "t") == "SD_plain"
    assert routed["called"] is False
    assert plain_solved["n"] == 1

    # With a tracer: must be routed to the native autograd run
    def f(p):
        out = H._run_plain(_traced_sim(p), "t")
        assert out == "SD_traced"
        return p

    ag.grad(f)(1.0)
    assert routed["called"] is True
    assert routed["sim"]._has_tracers is True


def test_run_plain_batch_and_list_still_local(monkeypatch):
    """The dict and list forms (the multi-objective batch of Autograd4) still solve locally one by
    one, with no routing."""
    routed = {"called": False}
    monkeypatch.setitem(H._STATE, "autograd_run",
                        lambda *a, **k: routed.__setitem__("called", True))
    # The whole container (dict/list) is handed to _solve to dispatch (since 2026-09-05 it can be
    # dispatched in parallel in remote mode); this stand-in returns a result of the matching shape
    def _fake_solve(s, tag, **kw):
        if isinstance(s, dict):
            return {k: f"SD:{tag}" for k in s}
        if isinstance(s, (list, tuple)):
            return [f"SD:{tag}" for _ in s]
        return f"SD:{tag}"
    monkeypatch.setattr(H, "_solve", _fake_solve)
    sim = _traced_sim(1.0)

    out_dict = H._run_plain({"a": sim}, "t")
    assert set(out_dict.keys()) == {"a"}
    out_list = H._run_plain([sim, sim], "t")
    assert len(out_list) == 2
    assert routed["called"] is False


def test_openem_run_entry_routes_traced_sim(monkeypatch):
    """Once port_cat has rewritten the notebook's ``web.run`` into ``openem_run``, a traced sim must
    be handed from this entry back to ``autograd_hook._run_plain`` (and on to tidy3d's own autograd
    run); it must not go straight into ``_openem_run_inner`` (the root cause of TidyFab0GC's
    grad_norm being zero for all 40 steps). A sim with no tracers still goes through
    ``_openem_run_inner``, so the original path is untouched."""
    from openem.nb import backend as nb
    calls = []

    def _traced(sim):
        return any(s._has_tracers for s in sim.values()) if isinstance(sim, dict) else sim._has_tracers

    def fake_inner(sim, **kw):
        calls.append(("inner", _traced(sim)))
        return "SD_inner"

    def fake_plain(sim, **kw):
        calls.append(("plain", _traced(sim)))
        return "SD_plain"

    monkeypatch.setattr(nb, "_openem_run_inner", fake_inner)
    monkeypatch.setattr(H, "_run_plain", fake_plain)
    monkeypatch.setitem(H._STATE, "autograd_run", object())   # marks the hook already installed

    # No tracer: the original path is unchanged (including the form that passes path=)
    assert nb.openem_run(_traced_sim(1.0), task_name="t", verbose=False) == "SD_inner"
    assert nb.openem_run({"a": _traced_sim(1.0)}, path="x", verbose=False) == "SD_inner"

    # With a tracer: both a single sim and a dict batch are handed back to _run_plain
    def f(p):
        assert nb.openem_run(_traced_sim(p), task_name="t", verbose=False) == "SD_plain"
        assert nb.openem_run({"a": _traced_sim(p)}, verbose=False) == "SD_plain"
        return p

    ag.grad(f)(1.0)
    assert calls == [("inner", False), ("inner", False), ("plain", True), ("plain", True)]
