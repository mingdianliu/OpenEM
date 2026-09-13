# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""The openem/nb/ package: two constraints that hold now that autograd_hook and shapegrad live there.

1. The old paths ``openem.autograd_hook`` and ``openem.shapegrad`` are forwarding shims and must be
   the **same object** as the new modules. External launch hooks may do ``autograd_hook._solve = ...``
   and then ``install``; external probes assign to ``_ah._solve`` and read ``_ah.CALLS``. A copy
   made by ``from ... import *`` would land those assignments on a different module and fail
   silently.
2. ``openem/nb/__init__.py`` imports no submodules: ``import openem.nb.shapegrad`` must not drag in
   solver, cupy or tidy3d. autograd_hook pulls solver at top level, and a CPU-only scene job also
   has to import nb.backend.

Both run in a subprocess, since other tests have long since imported these modules into this one.
"""
import importlib.util
import subprocess
import sys

import pytest

#: Symbols referenced by name from launch hooks, external probes, the tests and tools/golden.py
HOOK_SYMBOLS = ["install", "_solve", "_run_plain", "_run_one", "_STATE", "CALLS", "_FWD_STEPS",
                "_adjoint_steps", "_BatchShim", "_patch_costs", "scene_mod", "solver"]
SHAPEGRAD_SYMBOLS = ["apply", "revert", "_STATE", "_numerical_geometry_vjp"]


def _run(code: str) -> None:
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0, f"stdout={r.stdout}\nstderr={r.stderr}"


@pytest.mark.skipif(importlib.util.find_spec("cupy") is None,
                    reason="comparing module objects requires actually importing openem.solver, "
                           "which imports cupy at top level")
def test_old_paths_are_the_same_module_objects():
    _run(
        "import openem.autograd_hook as old_h, openem.nb.autograd_hook as new_h\n"
        "import openem.shapegrad as old_s, openem.nb.shapegrad as new_s\n"
        "from openem import autograd_hook as from_h, shapegrad as from_s\n"
        "assert old_h is new_h is from_h, (old_h, new_h, from_h)\n"
        "assert old_s is new_s is from_s, (old_s, new_s, from_s)\n"
        "old_h._solve = 'patched'\n"
        "assert new_h._solve == 'patched'\n"
        f"for n in {HOOK_SYMBOLS!r}: getattr(new_h, n)\n"
        f"for n in {SHAPEGRAD_SYMBOLS!r}: getattr(new_s, n)\n"
    )


def test_nb_package_init_stays_light():
    _run(
        "import openem.nb, openem.nb.shapegrad, sys\n"
        "heavy = sorted(m for m in sys.modules if m.split('.')[0] in ('cupy', 'tidy3d')"
        " or m in ('openem.solver', 'openem.nb.autograd_hook'))\n"
        "assert not heavy, heavy\n"
    )
