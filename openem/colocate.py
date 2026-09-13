# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Colocation from native Yee positions onto monitor coordinates. **Imports neither tidy3d nor cupy.**

Kernels store values only at their native Yee positions; colocation happens here. It is a linear
operation and commutes with the DFT, so the order does not matter, and keeping it in Python makes
changing the convention cheap.

Where the convention comes from: `tidy3d/components/data/dataset.py:319` is
``field_data.interp(...)``, xarray's **linear** interpolation. The target coordinates are
``colocation_boundaries`` from ``monitor_data.py:634``, i.e. the primal cell boundaries. A
zero-size dimension is not colocated but snapped to the monitor center.

Array convention: **the three spatial axes are the last three dimensions**; whatever precedes them
may be frequency or time.
"""

from __future__ import annotations

import numpy as np

from openem.grid import YEE_ON_EDGE, Grid

#: Whether each component sits at the cell center along a given axis, as opposed to on the cell
#: boundary. This is :data:`openem.grid.YEE_ON_EDGE` inverted; that table is the single source.
AT_CENTER: dict[str, tuple[bool, bool, bool]] = {
    k: tuple(not v for v in vs) for k, vs in YEE_ON_EDGE.items()   # type: ignore[misc]
}


def sample_coords(
    grid: Grid, origin: tuple[int, int, int], box: tuple[int, int, int], comp: str
) -> list[np.ndarray]:
    """Sample coordinates of this component along all three axes within the storage box, in metres."""
    out = []
    for ax in range(3):
        a = grid.axes[ax]
        base = a.centers if AT_CENTER[comp][ax] else a.edges
        out.append(base[origin[ax]: origin[ax] + box[ax]])
    return out


def _lerp_axis(vals: np.ndarray, src: np.ndarray, dst: np.ndarray, axis: int) -> np.ndarray:
    """Linearly interpolate along ``axis`` from the ``src`` coordinates to ``dst``.

    When src holds a single point, the value is broadcast as a constant.

    Going more than **one sample spacing** beyond ``src`` raises; within that margin the nearest
    endpoint is used.

    The split exists because one kind of overshoot is unavoidable: a ``colocate=True`` monitor asks
    for primal boundaries, and a component that sits at a **cell center** along that axis (Ex along
    x) has its outermost center half a cell inside the outermost boundary. There is no sample to
    interpolate across that half cell, so the nearest value is all there is. "The storage box was
    opened too small" is a different problem entirely: one scene was short by two full cells, got
    silently clamped, and reported a relative difference of 3.2 that looked like an error in
    epsilon.
    """
    if src.size == 1:
        return vals if dst.size == 1 else np.repeat(vals, dst.size, axis=axis)
    span = max(abs(float(src[1] - src[0])), abs(float(src[-1] - src[-2])))
    if dst[0] < src[0] - span or dst[-1] > src[-1] + span:
        raise ValueError(
            f"target coordinates [{dst[0]!r}, {dst[-1]!r}] lie outside the sampled range "
            f"[{src[0]!r}, {src[-1]!r}] by more than one sample spacing ({span!r}): the "
            "storage box is too small, and no extrapolation is done"
        )
    j = np.clip(np.searchsorted(src, dst) - 1, 0, src.size - 2)
    w = (dst - src[j]) / (src[j + 1] - src[j])
    shape = [1] * vals.ndim
    shape[axis] = dst.size
    w = w.reshape(shape)
    return np.take(vals, j, axis=axis) * (1.0 - w) + np.take(vals, j + 1, axis=axis) * w


def interp_to(
    vals: np.ndarray, coords: list[np.ndarray], dst: list[np.ndarray]
) -> np.ndarray:
    """Linearly interpolate the last three dimensions from ``coords`` to ``dst``.

    Both are lists of three coordinate arrays.
    """
    first = vals.ndim - 3
    out = vals
    for ax in range(3):
        out = _lerp_axis(out, coords[ax], dst[ax], first + ax)
    return out


def snap_floor(edges, center: float, tol: float) -> float:
    """tidy3d's snapping rule for a zero-thickness axis.

    Take the nearest primal boundary that is **not greater than the center** (``ind_min`` in
    ``discretize_inds`` is the last boundary that is ``<=``), falling back to the nearest one when
    the center is below every boundary. An earlier version used "nearest", which picked a different
    layer than tidy3d whenever the offset exceeded half a cell.

    ``tol`` is the tolerance added to the center when testing ``<=``. It is **passed by the caller
    and deliberately not unified**: scene/projection uses 1e-12 (um) while nb.backend and
    td_readout use 1e-15 (m). All three used to carry their own copy, and forcing one shared
    tolerance would switch to a different layer in edge cases.

    ``edges`` is ascending. Verified case by case against all three older spellings.
    """
    e = np.asarray(edges, dtype=np.float64)
    t0 = float(center) + tol
    below = e[e <= t0]
    return float(below[-1]) if below.size else float(e[int(np.argmin(np.abs(e - t0)))])


def at_point(
    vals: np.ndarray,
    coords: list[np.ndarray],
    target: tuple[float, float, float],
    mode: str = "interp",
) -> np.ndarray:
    """Take a single spatial point, collapsing the last three dimensions.

    Only tests/test_timemonitor.py uses this at present.

    Args:
        mode: ``"interp"`` trilinearly interpolates to the exact point; ``"snap"`` takes the
            nearest Yee point. **Measurement settled on interp**: on a resonance-finding case the
            interpolated ratios all sat on 1 with 6 of 7 points within 5e-3, while snapping
            scattered between 0.935 and 1.150. ``snap`` is kept for reproducing that comparison.
    """
    if mode not in ("interp", "snap"):
        raise ValueError(f"mode must be interp or snap, got {mode!r}")
    first = vals.ndim - 3
    out = vals
    for ax in range(3):
        c, x = coords[ax], float(target[ax])
        a = first + ax
        if mode == "snap" or c.size == 1:
            out = np.take(out, [int(np.argmin(np.abs(c - x)))], axis=a)
        else:
            out = _lerp_axis(out, c, np.array([x]), a)
    return out.reshape(out.shape[:first])
