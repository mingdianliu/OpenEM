# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Array pitch, padding and upload helpers. The z-dimension pitch is defined here and nowhere else.

Nothing here depends on the scene, only on the scalars nz, nzp and ny, so every table-building
module and the solver itself take the pitch from this module.
"""

from __future__ import annotations

from typing import NamedTuple

import cupy as cp
import numpy as np

from openem import knobs
from openem.grid import ravel_cells, unravel_cells


def _pitch(nz: int) -> int:
    """z-dimension pitch: 32 floats, i.e. 128-byte alignment.

    OPENEM_PAD_TEST forces non-zero padding; it is a debug knob for isolating the pitch mechanism.
    Every construction site that addresses by flat index must use the same pitch, which is why it
    is defined only here.
    """
    nzp = ((nz + 31) // 32) * 32
    if knobs.env("PAD_TEST"):
        nzp += 32
        return nzp
    # Memory guard. On a thin z domain (nz around 3, as in a 2D-material scene) padding up to 32
    # inflates the volume arrays by roughly 10x, which was measured to OOM an 80 GB H800. The
    # threshold is 2x: a case going 35 -> 64 (1.83x) was measured 16% *slower* when dropped back to
    # unaligned, since row alignment is worth about 1.4x and the main kernel is 91% of the time, so
    # that one has to be kept. The 10.7x case is still caught.
    if nzp > 2 * nz:
        return nz
    return nzp


class Dims(NamedTuple):
    """Domain dimensions.

    ``nz`` is the true number of z cells; ``nzp`` is the allocated length after the P11 pitch
    (``nzp >= nz``). Where they differ, ``[nz, nzp)`` is padding that no physical cell ever reads.

    Almost every function in the table-building phase needs two to four of these. Passing them
    flat means ``nx, ny, nz, nzp`` all over the place, so they travel as one group and each caller
    unpacks what it needs.
    """

    nx: int
    ny: int
    nz: int
    nzp: int

    @classmethod
    def of(cls, sc) -> "Dims":
        """Derive the dimensions from the Scene's logical cell counts and :func:`_pitch`."""
        nx, ny, nz = sc.shape
        return cls(nx, ny, nz, _pitch(nz))

    @property
    def shape(self) -> tuple[int, int, int]:
        """Allocated 3D shape after pitching, ``(nx, ny, nzp)``."""
        return (self.nx, self.ny, self.nzp)

    @property
    def n32(self) -> tuple:
        """The trailing three arguments of a volume kernel: nx, ny, nzp as int32."""
        return (np.int32(self.nx), np.int32(self.ny), np.int32(self.nzp))

    @property
    def ncell(self) -> int:
        """Allocated cell count after pitching, ``nx * ny * nzp``."""
        return self.nx * self.ny * self.nzp


def _f32(x) -> cp.ndarray:
    return cp.asarray(np.ascontiguousarray(x, dtype=np.float32))


def _f64(x) -> cp.ndarray:
    return cp.asarray(np.ascontiguousarray(x, dtype=np.float64))


def _pad_z_to(a: np.ndarray, nz: int, nzp: int) -> np.ndarray:
    """Pad a volume array along z up to nzp, filling with zeros. The padding is never read."""
    if nzp == nz:
        return np.ascontiguousarray(a)
    out = np.zeros(a.shape[:-1] + (nzp,), dtype=a.dtype)
    out[..., :nz] = a
    return out


def _repitch_cells(cell, ny: int, nz: int, nzp: int):
    """Convert a logical flat index (i*ny+j)*nz+k into its pitched form (i*ny+j)*nzp+k.

    Past a hundred million elements this runs on the GPU; on the host it was measured at 1.6 to
    1.9 s for the larger scenes. It returns a device array, which costs the caller nothing: the
    caller was going to ``cp.asarray`` it anyway, and that is a no-op on a device array. Integer
    arithmetic throughout, so the result is bitwise identical either way.
    """
    big = (getattr(cell, "size", 0) >= 1_000_000
           and knobs.env("INIT_GPU") != "0")
    xp = cp if big else np
    ii, jj, kk = unravel_cells(xp.asarray(cell).astype(xp.int64), (ny, nz))
    return ravel_cells(ii, jj, kk, (ny, nzp)).astype(xp.int32)


def _padz1_to(a, fill, nz: int, nzp: int):
    """Pad a 1D z table up to nzp.

    Rounding the launch up to whole blocks lets padding threads in [nz, ceil) pass the guard, since
    the guard uses the pitch, and those threads do read these tables. Safe pad values: index tables
    pad with nz-1 (inside the domain), masks pad with 0, spacings pad with 1, PML pads with the
    identity. A padding cell has ca=cb=0, so its E stays 0: the junk lands only in the padding and
    no physical cell ever reads it.
    """
    a = np.asarray(a)
    if nzp == nz:
        return np.ascontiguousarray(a)
    out = np.full(nzp, fill, dtype=a.dtype)
    out[:nz] = a
    return out


def _pad32_dev_to(v_dev, nz: int, nzp: int):
    """Pad a device array along z up to nzp and downcast it to float32."""
    p = cp.zeros(v_dev.shape[:-1] + (nzp,), cp.float32)
    p[..., :nz] = v_dev.astype(cp.float32)
    return p


def _pad32_host_to(a, nz: int, nzp: int):
    """Upload a host array in chunks and downcast it, so peak memory is one chunk.

    Used by the very large scenes.
    """
    out = cp.zeros((a.shape[0], a.shape[1], nzp), cp.float32)
    step = max(1, min(a.shape[0], int(256e6 // max(a.shape[1] * nz * 8, 1))))
    for i0 in range(0, a.shape[0], step):
        i1 = min(a.shape[0], i0 + step)
        blk = cp.asarray(a[i0:i1])
        out[i0:i1, :, :nz] = blk.astype(cp.float32)
        del blk
    return out


def _d2h_pinned(dev_arr, cache: dict, key: str = "_pin_buf"):
    """Device-to-host copy reusing a **pinned** buffer, one per monitor, rebuilt on demand.

    A non-pinned destination has to go through CUDA's staging buffer, measured at only 1.55 GB/s
    for 500 MB; pinned memory saturates PCIe. The bytes copied are the same either way, so the
    result is bitwise identical to ``cp.asnumpy``.
    """
    buf = cache.get(key)
    if buf is None or buf.shape != dev_arr.shape or buf.dtype != dev_arr.dtype:
        mem = cp.cuda.alloc_pinned_memory(dev_arr.nbytes)
        buf = np.frombuffer(mem, dev_arr.dtype, dev_arr.size).reshape(dev_arr.shape)
        cache[key] = buf
    dev_arr.get(out=buf)
    return buf


def _plain_neighbors(nxt, prv, dims) -> bool:
    """Whether the neighbour tables are the plain "one cell forward, one cell back" pattern.

    The endpoints are handled by a guard inside the kernel.

    The P41 fusion avoids a race only because the E value read across a block boundary is the one
    cell in the nxt direction, and that cell lies on the ring which is skipped and left to the edge
    kernel. **Symmetry folding points nxt at a mirrored position instead**, and that cell may be
    under write by another block, so a folded scene is never fused.
    """
    for a in range(3):
        n = dims[a]
        nx_ = cp.asnumpy(nxt[a]).astype(np.int64)
        pv = cp.asnumpy(prv[a]).astype(np.int64)
        if nx_.size < n or pv.size < n:
            return False
        if not (nx_[:n - 1] == np.arange(1, n)).all():
            return False
        if not (pv[1:n] == np.arange(n - 1)).all():
            return False
    return True


def _padz_axis(v, axis, fill, nz, nzp):
    """Pad a 1D table up to the pitch on the z axis (``axis==2``); return the other two unchanged.

    Rounding the launch up to whole blocks spills into padding threads in [nz, nzp), since the
    guard uses the pitch, and those threads read these tables. Safe pad values: index tables pad
    with nz-1 (inside the domain), masks pad with 0, spacings pad with 1, PML pads with the
    identity, pec pads with 1. A padding cell has ca=cb=0, so its E stays 0: the junk lands only in
    the padding and no physical cell ever reads it.
    """
    return v if axis != 2 else _padz1_to(v, fill, nz, nzp)
