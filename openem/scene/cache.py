# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Shared pieces of the disk cache: hashing a model's spatial arrays, the cache root, atomic writes.

Shared by the eps cache (media.epsilon_complex) and the mode basis cache (modes._cached_solve).
"""

from __future__ import annotations

import os
import pathlib

import numpy as np


def cache_root() -> pathlib.Path:
    """The ``cache/`` directory beside the repository: eps under ``cache/epsilon``, mode bases
    under ``cache/modes``.
    """
    return pathlib.Path(__file__).resolve().parents[3] / "cache"


def atomic_write(path: pathlib.Path, writer) -> None:
    """Write a temporary file and then rename it, so an interruption cannot leave half a file to be
    mistaken for a cache hit.

    The temporary name carries the pid: two processes computing the same data write the same key,
    and sharing one ``.part`` name would let the first rename win while the second replace raised
    FileNotFoundError. A rename within a directory is atomic, so whoever arrives later overwrites,
    with identical content; if the rename fails the temporary file is removed. ``writer(tmp_path)``
    writes the content.
    """
    tmp = path.with_name(f"{path.stem}.{os.getpid()}.part{path.suffix}")
    writer(tmp)
    try:
        tmp.replace(path)
    except OSError:
        tmp.unlink(missing_ok=True)


def hash_data_arrays(obj, h, depth: int = 0) -> None:
    """Feed a model's spatial arrays into a hash.

    **Why this is necessary**: tidy3d's ``json`` **silently drops** the spatial data of things like
    CustomMedium and TriangleMesh (the warning in `scene/build.py::from_file` is the same issue).
    Using only ``structure.json`` as the cache key makes different designs collide on **the same
    key** and receive the same eps. Measured on an inverse-design case, four different parameter
    sets, including all zeros and all ones, produced element-wise identical eps grids, so the
    optimization responded to its parameters not at all.
    """
    if depth > 6 or obj is None:
        return
    if isinstance(obj, np.ndarray):
        h.update(np.ascontiguousarray(obj).tobytes())
        return
    if isinstance(obj, dict):
        for v in obj.values():
            hash_data_arrays(v, h, depth + 1)
        return
    if isinstance(obj, (list, tuple)):
        for v in obj:
            hash_data_arrays(v, h, depth + 1)
        return
    vals = getattr(obj, "values", None)          # an xarray DataArray
    if isinstance(vals, np.ndarray):
        h.update(np.ascontiguousarray(vals).tobytes())
        return
    flds = getattr(obj, "model_fields", None) or getattr(obj, "__fields__", None)
    if flds:
        for name in flds:
            hash_data_arrays(getattr(obj, name, None), h, depth + 1)
