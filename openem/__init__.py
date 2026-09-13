# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""OpenEM: a CUDA FDTD solver that takes its input from Tidy3D.

The module boundaries are deliberately narrow:

- ``scene/`` is the package that imports tidy3d (``td.Simulation`` -> Scene) and the one place
  where lengths change from um to metres. ``td_readout`` has to run inside a job container that
  has no ``scene/``, so it carries its own copy of the same ``_UM`` constant.
  ``nb/`` holds the notebook-side hooks (autograd_hook, shapegrad) and is the other place allowed
  to import tidy3d. The old paths ``openem.autograd_hook`` and ``openem.shapegrad`` are forwarding
  shims pointing at the very same module objects under ``nb/`` (tests/test_nb_shim.py).
  No other module may import tidy3d **at top level**; tests/test_no_tidy3d_leak.py locks the solver
  chain. The one deferred import inside a function lives in ``normalize.py`` and is reached only by
  the handful of cases where the reference implementation did not normalize the source.
- The solver chain, ``solver`` / ``device_tables`` / ``sources_setup`` / ``monitors_setup`` /
  ``dispersion_setup`` / ``pitch`` / ``coeffs`` / ``fusion`` / ``readout`` / ``device``, imports
  cupy at top level (``setup_tables`` is only a re-export shim over those). ``tfsf1d`` imports it
  locally, inside the optional GPU path.
- Everything else (``model grid cpml waveform serialize flux fold normalize colocate td_readout
  shutoff apodization tfsf_oblique results units knobs``) sees numpy only and has to run on a
  machine without a GPU; tests/test_no_cupy_leak.py locks that.
- Every ``OPENEM_*`` environment knob in the package registers its default and its purpose in
  ``knobs.py`` and is read live through ``knobs.env``, on the solver side as well as in ``scene/``
  and ``nb/``. The notebook integration layer ``nb/backend.py`` reads os.environ directly.
"""

__all__ = ["cpml", "grid", "scene", "waveform", "install", "run"]


def install(steps: int | None = None, max_cells: int | None = 80_000_000,
            staircase_local_modesolver: bool = True) -> None:
    """Install OpenEM as tidy3d's solver backend.

    After this call, ``tidy3d.web.run``, ``web.Batch`` and autograd's ``grad`` all land on the
    local GPU, with no change to the notebook.

    Args:
        steps: pin the number of time steps. The default, None, means the scene's nominal step
            count plus early termination. The hook keeps the adjoint solve on the same count as
            the forward one.
        max_cells: per-simulation cell ceiling. Anything above it is refused outright (fail
            closed), which stops an oversized scene from being submitted by accident.
        staircase_local_modesolver: staircase the local ``ModeSolver`` and ``sim.epsilon`` that
            the notebook calls **directly**, matching an environment without tidy3d-extras, which
            is what the reference notebooks ran against. Work OpenEM does internally on the
            reference solver's behalf (mode source basis, ``web.run(ModeSolver)``) is unaffected.
    """
    from openem.nb import autograd_hook
    autograd_hook.install(steps=steps, max_cells=max_cells)
    if staircase_local_modesolver:
        from openem.scene import media
        media.install_notebook_staircase()


def run(sim, task_name: str = "openem", verbose: bool = True, **kw):
    """Solve a ``td.Simulation`` on the local GPU directly, without installing the hook.

    Also accepts a list of simulations, or a ``{name: sim}`` dict. Returns a
    ``td.SimulationData``, or the matching list or dict.
    """
    from openem.nb import backend
    return backend.run(sim, task_name=task_name, verbose=verbose, **kw)
