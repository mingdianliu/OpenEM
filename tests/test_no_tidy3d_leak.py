# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Lock the module boundary: the solver chain must not import tidy3d, not even transitively.

Why this test has to exist: a GPU job container may not be able to install tidy3d (the preinstalled
packages conflict with its dependencies), so the solving path has to stay clean. It went wrong once
in practice: ``serialize.py`` imported a dataclass from ``scene.py``, ``scene.py`` had
``import tidy3d`` at the top level, so ``import openem.serialize`` dragged tidy3d in and the job
died with ``ModuleNotFoundError``.

Only ``scene/`` (which builds the scene) and ``nb/`` (the notebook-side hooks autograd_hook and
shapegrad) are allowed to import tidy3d; the old paths ``openem.autograd_hook`` and
``openem.shapegrad`` are forwarding shims (tests/test_nb_shim.py).

Run in a subprocess, because in this process tidy3d was long since imported by another test, which
makes checking sys.modules meaningless.
"""

import subprocess
import sys

import pytest

#: The modules a job container uses; not one of them may touch tidy3d
SOLVER_CHAIN = ["openem.model", "openem.grid", "openem.cpml", "openem.waveform",
                "openem.serialize", "openem.solver", "openem.flux",
                "openem.shutoff", "openem.device"]


def test_solver_chain_is_tidy3d_free():
    """This runs on machines without cupy too: it checks tidy3d, which has nothing to do with
    having a GPU.

    Import them one at a time and skip whichever stops at ``import cupy``: ``openem.solver`` needs
    cupy at the top level, so on a CPU-only machine (public CI included) it cannot be imported at
    all. On CPU this therefore covers what each module does before it reaches cupy, and only a
    machine with cupy gets full coverage. But when the boundary is broken it is almost always an
    extra ``import tidy3d`` at some module's top level, and either kind of machine catches that on
    the spot.
    """
    code = (
        "import importlib, sys\n"
        f"for name in {SOLVER_CHAIN!r}:\n"
        "    try:\n"
        "        importlib.import_module(name)\n"
        "    except ModuleNotFoundError as e:\n"
        "        if e.name != 'cupy':\n"
        "            raise\n"
        "leaked = sorted(m for m in sys.modules if m.split('.')[0] == 'tidy3d')\n"
        "assert not leaked, 'these modules dragged tidy3d in: ' + repr(leaked)\n"
        "print('OK')\n"
    )
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0, f"stdout={r.stdout}\nstderr={r.stderr}"


def test_scene_builder_does_import_tidy3d():
    """The other direction: confirm that scene.py really **is** that boundary, since it is supposed
    to depend on tidy3d.

    This can only run where tidy3d is installed. A job container has no tidy3d, so it simply skips;
    otherwise it would fail because tidy3d is absent, and that is exactly the case this test cannot
    tell apart.
    """
    pytest.importorskip("tidy3d", reason="tidy3d is needed to show that scene.py depends on it")
    code = (
        "import openem.scene, sys\n"
        "assert 'tidy3d' in sys.modules, 'scene.py should be the only module importing tidy3d'\n"
    )
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0, f"stdout={r.stdout}\nstderr={r.stderr}"
