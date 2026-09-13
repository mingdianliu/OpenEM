# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Forwarding shim: the implementation moved to :mod:`openem.nb.shapegrad`.

Aliased through ``sys.modules`` so both names refer to the same module object, which keeps
``_STATE`` single; tests/test_nb_shim.py locks that.
"""
import sys

from openem.nb import shapegrad as _impl

sys.modules[__name__] = _impl
