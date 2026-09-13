# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Notebook-side hooks: wire tidy3d's web.run and autograd through to OpenEM.

- :mod:`~openem.nb.autograd_hook`  takes over tidy3d's forward and adjoint entry points along with
  the Job and Batch shells; once installed, web.run and ag.grad land on OpenEM
- :mod:`~openem.nb.shapegrad`      patches for adjoint gradients with respect to geometry
  parameters (zeroing the mirror side, perturbing the eps map)

Like ``scene/``, this package is allowed to import tidy3d; tests/test_no_tidy3d_leak.py locks only
the solver chain. **This __init__ imports no submodule**: ``autograd_hook`` pulls in solver and
cupy at top level, while a CPU-only scene job also has to ``import openem.nb.backend`` (which
imports this package only on demand). tests/test_nb_shim.py locks that. The old paths
``openem.autograd_hook`` and ``openem.shapegrad`` are forwarding shims pointing at the very same
module objects here, so that an external startup hook assigning to ``autograd_hook._solve`` takes
effect.
"""
