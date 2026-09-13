# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Convert field phasors into Tidy3D's output convention. **Imports neither tidy3d nor cupy.**

The fields Tidy3D hands the user (``colocate=True``, the tutorial default) do not sit at their
native Yee positions: it linearly interpolates all six components onto one common set of cell
boundaries (``field_data.interp`` in its source), takes a single layer along any zero-thickness
dimension while labelling it with the monitor center, and stores only the reduced domain for
symmetric cases. All of these **output definitions** are applied once, here. The resulting
``tdf:<monitor>`` entries are already the final numbers, normalization included, so everything
downstream just copies them and never interpolates again. The target coordinates come from a
pre-generated ``td_coords.npz`` in the case directory, point for point.

Internal planes carrying the ``__mode_plane`` suffix are not converted: mode decomposition and
near-to-far projection both need the native Yee positions.
"""
from __future__ import annotations

import numpy as np

from openem import colocate, normalize
from openem import knobs
from openem.model import COMPONENTS, MODE_PLANE_SUFFIX

#: Same source as PLANE_SUFFIX in openem/scene/modes.py (model.MODE_PLANE_SUFFIX). We do not import
#: the scene package: its __init__ drags in tidy3d, and this module has to run inside a job
#: container where tidy3d cannot be installed.
_PLANE_SUFFIX = MODE_PLANE_SUFFIX

#: um -> m, the same value as openem.units.UM, copied here to avoid implying a dependency
_UM = 1e-6


def _snap_thin_axis_target(t_m, srcax, edges_m):
    """On a zero-thickness axis, where the target is a single point, move the interpolation target
    to the nearest primal boundary, i.e. the E layer.

    This matches tidy3d's "take the layer, label it with the center" behaviour; the caller keeps
    the label. Controlled by OPENEM_FIELDMON_SNAP, default 0 meaning no snapping; see the inline
    comment below.
    """
    # Its own knob, OPENEM_FIELDMON_SNAP, default 0: across the real comparison cases it improved 5
    # field rows and worsened 13, so it was rolled back. Snapping on mode planes goes through a
    # separate knob, OPENEM_COLOC_SNAP, and is unaffected.
    if t_m.size != 1 or knobs.env("FIELDMON_SNAP") != "1":
        return t_m
    e = np.asarray(edges_m, dtype=np.float64)
    if e.size < 2:
        return t_m
    # tidy3d's rule: take the boundary that is not greater than the center (ind_min in
    # discretize_inds is the last boundary <= center), falling back to the nearest one when the
    # center is below every boundary. An earlier version used "nearest", which picked a different
    # layer than tidy3d whenever the offset exceeded half a cell.
    snapped = colocate.snap_floor(e, float(t_m[0]), 1e-15)      # tolerance 1e-15 m (scene/projection uses 1e-12 um)
    lo, hi = float(np.min(srcax)), float(np.max(srcax))
    return np.array([min(max(snapped, lo), hi)])


def finalize(sc, field_phasors: dict, plans: dict) -> dict:
    """Convert each monitor's phasors into Tidy3D's convention, using pre-generated target coordinates.

    Args:
        sc: the Scene, for the grid, the monitor declarations and the normalization.
        field_phasors: ``Result.field_phasors``, ``{name: (6, nf, ni, nj, nk)}``.
        plans: the mapping from ``np.load("td_coords.npz")``, keyed ``"<monitor>|<component>:<axis>"``.

    Returns:
        ``{"tdf:<name>|<component>": complex array (x,y,z,f),
        "tdc:<name>|<component>:<axis>": coordinates in um}``. A monitor with no plan does not
        appear, and everything downstream falls back to the older path on its own.
    """
    out: dict[str, np.ndarray] = {}
    fmons = {m.name: m for m in sc.field_monitors}
    for name, raw in field_phasors.items():
        if name.endswith(_PLANE_SUFFIX) or name not in fmons:
            continue
        m = fmons[name]
        if not any(f"{name}|{c}:x" in plans for c in COMPONENTS):
            continue
        fs = np.asarray(m.freqs, float)
        Dn = normalize.field_norm(sc, fs)
        # With colocate=False the six components sit at different coordinates along a
        # zero-thickness axis (E boundaries vs H cell centers), so no snapping. With
        # colocate=True all six share one label coordinate and must snap to the E layer.
        _thin_same = True
        for ax in "xyz":
            cs = [np.asarray(plans[f"{name}|{c}:{ax}"], np.float64) for c in COMPONENTS if f"{name}|{c}:{ax}" in plans]
            if cs and cs[0].size == 1 and any(c.size != 1 or abs(float(c[0]) - float(cs[0][0])) > 1e-9 for c in cs):
                _thin_same = False
        for ci, comp in enumerate(COMPONENTS):
            key = f"{name}|{comp}"
            if f"{key}:x" not in plans:
                continue
            src = colocate.sample_coords(sc.grid, m.origin, m.box, comp)
            tgt = []
            for ai, (ax, sa) in enumerate(zip("xyz", src)):
                t = np.asarray(plans[f"{key}:{ax}"], np.float64) * _UM
                lo, hi = float(sa.min()), float(sa.max())
                t = np.clip(t, lo, hi)
                # Zero-thickness axis: snap the interpolation target to the nearest E layer
                # (tidy3d takes the layer and labels it with the center)
                tgt.append(_snap_thin_axis_target(t, sa, sc.grid.axes[ai].edges) if _thin_same else t)
            oc = colocate.interp_to(raw[ci], src, tgt)      # shape (nf, x, y, z)
            oc = oc / Dn[:, None, None, None]
            out[f"tdf:{key}"] = np.moveaxis(oc, 0, -1)      # shape (x, y, z, f)
            for ax in "xyz":
                out[f"tdc:{key}:{ax}"] = np.asarray(plans[f"{key}:{ax}"],
                                                    np.float64)
    return out
