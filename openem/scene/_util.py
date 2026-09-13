# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Small helpers shared by the scene-side modules: Yee coordinates, component keys, signs, plane
indices.

Indices, strings and coordinate lookups only, no physics. The functions that take a ``sim`` (a
tidy3d object) stay inside ``scene/``, which is part of the package boundary
(tests/test_no_tidy3d_leak.py).
"""

from __future__ import annotations

import numpy as np

from openem.grid import Grid, nearest_edge_index
from openem.units import UM

#: Keys of the three Yee components in the array dicts (``eps["ex"]`` and so on).
E_KEYS = ("ex", "ey", "ez")

#: Component keys for tidy3d's ``sim.epsilon(coord_key=...)``.
COORD_KEYS = ("Ex", "Ey", "Ez")


def yee_coords(sim, comp: str, magnetic: bool = False, meters: bool = True) -> list[np.ndarray]:
    """Coordinates ``[x, y, z]`` of ``sim.grid.yee.E/H.<comp>`` along all three axes, as float64.

    With ``meters=True`` they are multiplied by :data:`UM` to become metres; the same
    multiplication in the same order, so this is bitwise stable. With ``False`` they stay in
    tidy3d's um. ``comp`` is ``"x"``, ``"y"`` or ``"z"``.
    """
    yee = getattr(sim.grid.yee.H if magnetic else sim.grid.yee.E, comp)
    if meters:
        return [np.asarray(getattr(yee, a), dtype=np.float64) * UM for a in "xyz"]
    return [np.asarray(getattr(yee, a), dtype=np.float64) for a in "xyz"]


def sign_of(direction) -> int:
    """tidy3d's ``"+"`` or ``"-"`` mapped to ``+1`` or ``-1``."""
    return +1 if str(direction) == "+" else -1


def plane_index(grid: Grid, ax: int, center_um: float) -> int:
    """Snap a source or monitor plane, given by its center in um, to the nearest primal edge index
    on axis ``ax``.
    """
    return nearest_edge_index(grid.axes[ax].edges, center_um * UM)


def upstream_h_index(ks: int, direction: int) -> int:
    """Cell index of the H just upstream of source plane ``ks``: ``ks-1`` for the + direction,
    ``ks`` for the - direction.
    """
    return ks - 1 if direction > 0 else ks
