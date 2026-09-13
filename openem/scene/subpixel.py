# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""tidy3d's global local-subpixel switch: setting it explicitly, marking the scope in which we
stand in for the hosted solver, and the notebook hooks.

This is switch management at notebook runtime rather than "eps and conductivity", which is why it
lives here instead of in media.py; media and ``scene/__init__`` re-export the same names.
"""

from __future__ import annotations

import contextlib as _contextlib
import contextvars as _contextvars
import functools as _functools

from tidy3d import config


def set_local_subpixel(enable: bool) -> list[str] | None:
    """Set tidy3d's local-subpixel switch **explicitly** to ``enable``.

    The default of ``use_local_subpixel=None``, meaning "use it if present", must not be relied on:
    it makes the eps convention depend on whether extras happens to be installed, which fails
    silently in both directions.

    Args:
        enable: should equal ``bool(sim.subpixel)``.

    Returns:
        The list of licensed features when ``enable=True``, otherwise None.
    """
    if not enable:
        config.simulation.use_local_subpixel = False
        return None
    return require_local_subpixel()


@_contextlib.contextmanager
def local_subpixel(enable: bool):
    """Within the scope, set ``use_local_subpixel`` explicitly to ``enable`` (see
    :func:`set_local_subpixel`), and restore the previous value on exit.

    This is the shared "save, change, restore in try/finally" spelling.
    """
    prev = config.simulation.use_local_subpixel
    set_local_subpixel(enable)
    try:
        yield
    finally:
        config.simulation.use_local_subpixel = prev


#: Marks the scope in which OpenEM stands in for the hosted solver: the web.run replacement, scene
#: construction, the mode basis, and web.run(ModeSolver).
_CLOUD_EMU = _contextvars.ContextVar("openem_cloud_emulation", default=False)


@_contextlib.contextmanager
def cloud_emulation():
    """Wraps everything OpenEM does in place of the hosted solver.

    Inside this scope the hooks installed by ``install_notebook_staircase`` do not take effect, so
    eps and the mode solve continue to follow ``sim.subpixel`` through extras, as the hosted solver
    does. On exit ``use_local_subpixel`` is restored, so the value that ``set_local_subpixel`` set
    during scene construction no longer leaks back into the notebook.

    Also usable as a decorator: ``@cloud_emulation``.
    """
    prev = config.simulation.use_local_subpixel
    tok = _CLOUD_EMU.set(True)
    try:
        yield
    finally:
        _CLOUD_EMU.reset(tok)
        config.simulation.use_local_subpixel = prev


def install_notebook_staircase() -> list[str]:
    """Staircase every local mode solve and ``sim.epsilon`` the notebook calls **directly**.

    That covers ``ModeSolver.solve`` and ``.data``, the waveguide plugin,
    ``ModeSimulation.run_local`` and ``Simulation.epsilon``. The environment the reference
    notebooks ran in had no tidy3d-extras, so those numbers were staircased there to begin with;
    one of the reference outputs still carries the warning "Use the remote mode solver with
    subpixel averaging".

    Nothing inside the ``cloud_emulation`` scope, i.e. OpenEM itself, is affected. Install this only
    from the notebook startup hook; installing it twice is a no-op.

    Returns:
        The names of the methods actually wrapped.
    """
    from tidy3d.components.mode.mode_solver import ModeSolver
    from tidy3d.components.mode.simulation import ModeSimulation
    from tidy3d.components.simulation import Simulation

    done = []
    for cls, name in ((ModeSolver, "_solve_all_freqs"), (ModeSolver, "_solve_all_freqs_relative"),
                      (ModeSimulation, "run_local"), (Simulation, "epsilon")):
        owner = next(k for k in cls.__mro__ if name in k.__dict__)   # epsilon is defined on AbstractYeeGridSimulation
        orig = owner.__dict__[name]
        if getattr(orig, "_openem_staircase", False):
            continue

        def wrapped(self, *a, _orig=orig, **k):
            if _CLOUD_EMU.get():
                return _orig(self, *a, **k)
            with local_subpixel(False):
                return _orig(self, *a, **k)

        wrapped = _functools.wraps(orig)(wrapped)
        wrapped._openem_staircase = True
        setattr(owner, name, wrapped)
        done.append(f"{owner.__name__}.{name}")
    return done


def require_local_subpixel() -> list[str]:
    """Confirm that ``local_subpixel`` from ``tidy3d-extras`` is available, turn it on, or raise.

    Returns:
        The list of features available under the licence.

    Raises:
        NotImplementedError: extras is not installed, or the licence does not include
            ``local_subpixel``.
    """
    try:
        import tidy3d_extras
    except ImportError as exc:
        raise NotImplementedError(
            "subpixel=True needs tidy3d-extras to obtain the subpixel-averaged eps tensor, and it "
            "is not installed here. Two ways forward: install tidy3d-extras (Flexcompute's add-on "
            "package), or build the simulation with subpixel=False, which staircases and matches "
            "the convention of a reference environment without extras."
        ) from exc
    feats = sorted(tidy3d_extras.extension._features())
    if "local_subpixel" not in feats:
        raise NotImplementedError(
            f"tidy3d-extras is installed, but the licence does not include 'local_subpixel'; it has {feats}"
        )
    config.simulation.use_local_subpixel = True
    return feats
