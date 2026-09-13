# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Forwarding shim: the implementation moved to :mod:`openem.nb.autograd_hook`, since only
``scene/`` and ``nb/`` are allowed to import tidy3d.

This aliases through ``sys.modules`` rather than ``from ... import *`` because launch hooks and
external probe scripts **assign to module attributes** (``autograd_hook._solve = _solve_remote``,
``_ah.CALLS``). Both names have to refer to the same module object for those assignments to be
visible to ``_run_one``; tests/test_nb_shim.py locks that.
"""
import sys

from openem.nb import autograd_hook as _impl

sys.modules[__name__] = _impl
