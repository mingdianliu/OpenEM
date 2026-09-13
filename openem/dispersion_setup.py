# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Reordering of dispersion entries and the absorber tables: (am1, b) deduplication, pole
bucketing, deferred inlining, the interior / boundary-shell split, absorber slabs and the
in-layer P decay factors. The thresholds live in fusion.py.
(Moved verbatim out of setup_tables.py, 2026-09.)
"""

from __future__ import annotations

from typing import NamedTuple

import cupy as cp
import numpy as np

from openem import knobs
from openem.device import grid_1d
from openem.grid import unravel_cells
from openem.model import Scene
from openem.pitch import Dims


class LeanTables(NamedTuple):
    """The per-axis tables that P25b works backwards from when marking out the interior. The E
    side looks at ``(pml_e, pec, prv, mpv)`` and the H side at ``(pml_h, nxt, mnx)``; an
    "identity update" simply means those coefficients take their identity values."""

    pml_e: list
    pec: list
    prv: list
    mpv: list
    pml_h: list
    nxt: list
    mnx: list


def _coef_lut(am1: np.ndarray, b: np.ndarray) -> tuple[dict, np.ndarray]:
    """Lossless deduplication table for the complex pairs (am1, b). Returns
    ``({"am1": table, "b": table}, per-pole index)``.

    Deduplication is by the quadruple (am1.re, am1.im, b.re, b.im); the values go into the table
    **unchanged** and a lookup returns the same float, so it is bitwise equivalent (unlike a
    quantizing lookup table). The cast to float32 happens at upload time and is exactly the same
    conversion as the old path's ``_f32(d.am1.real)``.
    """
    # With hundreds of millions of rows, numpy's axis=0 deduplication has to sort everything
    # (measured argsort: 28.8 s for hBN, 191 s for Ring, 208 s for GroupDelay, with the H800 idle
    # throughout). The GPU path uses a hash plus a bit-pattern check.
    # The memory layout of complex128 is exactly [re, im, re, im, ...], so it is viewed as two
    # columns directly, saving the np.stack on the host (measured 3.9 s for GroupDelay).
    if (am1.shape[0] >= 1_000_000
            and knobs.env("INIT_GPU") != "0"):
        try:
            _a = cp.asarray(am1).view(cp.float64).reshape(-1, 2)
            _b = cp.asarray(b).view(cp.float64).reshape(-1, 2)
            qd = cp.concatenate([_a, _b], axis=1)
            del _a, _b
            u64 = qd.view(cp.uint64).reshape(qd.shape)
            h = cp.zeros(qd.shape[0], cp.uint64)
            for _c in range(4):          # splitmix64-style mixing, integer-only, deterministic
                x = u64[:, _c] ^ h
                x = (x ^ (x >> cp.uint64(30))) * cp.uint64(0xBF58476D1CE4E5B9)
                x = (x ^ (x >> cp.uint64(27))) * cp.uint64(0x94D049BB133111EB)
                h = x ^ (x >> cp.uint64(31))
            _uh, first, inv_d = cp.unique(h, return_index=True,
                                          return_inverse=True)
            rep_u64 = u64[first]
            ok = bool((rep_u64[inv_d] == u64).all())   # byte-exact bit-pattern check: no clash
            if ok:
                uq = cp.asnumpy(qd[first])
                # inv is handed back while it is still in device memory (the caller uploads it
                # as cidx right afterwards), which saves two big transfers, an int64 download
                # and an int32 upload.
                inv = inv_d.astype(cp.int32)
                del qd, u64, h, rep_u64, inv_d, _uh, first
                cp.get_default_memory_pool().free_all_blocks()
                return ({"am1": uq[:, 0] + 1j * uq[:, 1],
                         "b": uq[:, 2] + 1j * uq[:, 3]}, inv)
            del qd, u64, h, rep_u64, inv_d, _uh, first
            cp.get_default_memory_pool().free_all_blocks()
        except Exception:
            cp.get_default_memory_pool().free_all_blocks()
    quad = np.stack([am1.real, am1.imag, b.real, b.imag], axis=1)
    uniq, inv = np.unique(quad, axis=0, return_inverse=True)
    return ({"am1": uniq[:, 0] + 1j * uniq[:, 1], "b": uniq[:, 2] + 1j * uniq[:, 3]},
            inv)


def _absorber_slabs(sc: Scene) -> list[dict]:
    """Device-side data and launch configuration of the absorber slabs. An empty list means the
    scene has no Absorber.

    The kernel is only launched when there are slabs, so a scene without an Absorber is bitwise
    unchanged; that is what keeps it a single path (for the mechanism, see kernels/absorber.cu).
    """
    nx, ny, nz = sc.shape
    out = []
    for sl in sc.absorbers:
        nlay = sl.decay_int.size
        # Bounding box of the slab: the absorbing axis has only nlay layers, the other two span
        # everything.
        # (j,k) is flattened into one dimension (k innermost) and i goes to blockIdx.y, so memory
        # coalescing no longer depends on the absorbing axis, and lane utilization is not hurt by
        # nlay being only a few dozen layers (see absorber.cu)
        box = [nx, ny, nz]
        box[sl.axis] = nlay
        org = [0, 0, 0]
        org[sl.axis] = int(sl.g0)
        nthr = 256
        grid = ((box[1] * box[2] + nthr - 1) // nthr, box[0], 1)
        block = (nthr, 1, 1)
        out.append({
            "grid": grid, "block": block,
            "decay_int": cp.asarray(sl.decay_int, dtype=cp.float32),
            "decay_half": cp.asarray(sl.decay_half, dtype=cp.float32),
            "axis": np.int32(sl.axis),
            "org": tuple(np.int32(v) for v in org),
            "box": tuple(np.int32(v) for v in box),
        })
    return out


def _absorber_disp_decay(sc: Scene) -> dict | None:
    """Indices and per-step P decay factors of the dispersion terms inside an absorber layer.
    ``None`` = there are no dispersion terms inside the layer.

    P hangs off an E component, so the factor takes the **E-family** profile value at that
    component's staggered position (half a cell along the normal, a whole cell tangentially, the
    same as ``is_e=1`` in ``absorb_slab``; in a corner several slabs multiply together). The
    kernel (``absorb_disp_hist``) is only launched when this returns non-None, so a scene with no
    absorber, with no dispersion, or with all its dispersion outside the layer is bitwise
    unchanged (a single path, the same approach as the absorber).
    """
    if not sc.absorbers or sc.dispersion is None or not sc.dispersion.n_entry:
        return None
    d = sc.dispersion
    _, ny, nz = sc.shape
    # The whole block runs on the GPU (Ring has 100 million entries, measured at 5.96 s on the
    # host). The formula is unchanged entry by entry, only the place it runs has changed.
    # Same size threshold as P53: on small entry counts, moving it to the GPU does not pay
    _min = int(knobs.env("ABSDISP_GPU_MIN"))
    xp = (cp if (knobs.env("ABSDISP_GPU") != "0"
                 and int(d.n_entry) >= _min) else np)
    _dev = xp is cp
    cell = cp.asarray(d.cell) if _dev else d.cell
    comp = cp.asarray(d.comp) if _dev else d.comp
    coords: dict[int, np.ndarray] = {}

    def _coord(ax: int):
        # cell = (i*ny + j)*nz + k; take the global index along the absorbing axis, shared by
        # both ends of that axis
        if ax not in coords:
            coords[ax] = (cell // (ny * nz) if ax == 0
                          else (cell // nz) % ny if ax == 1
                          else cell % nz)
        return coords[ax]

    factor = xp.ones(d.n_entry)
    inside = xp.zeros(d.n_entry, dtype=bool)
    # The fused single-pass kernel has to replay the E-family decay in the order the slabs are
    # applied. An entry cell lands in at most 3 slabs (the two ends of one axis do not overlap),
    # so record the factors in that order (1.0f by default)
    fs = xp.ones((3, d.n_entry), dtype=xp.float32)
    cnt = xp.zeros(d.n_entry, dtype=xp.int64)
    for sl in sc.absorbers:
        lay = _coord(sl.axis) - sl.g0
        m = (lay >= 0) & (lay < sl.decay_int.size)
        if not bool(m.any()):
            continue
        lm = lay[m]
        d_int = xp.asarray(sl.decay_int) if _dev else sl.decay_int
        d_half = xp.asarray(sl.decay_half) if _dev else sl.decay_half
        vals = xp.where(comp[m] == sl.axis, d_half[lm], d_int[lm])
        factor[m] *= vals
        idxs = xp.flatnonzero(m)
        assert int(cnt[idxs].max()) <= 2, "an entry cell landed in more than 3 slabs"
        fs[cnt[idxs], idxs] = vals.astype(xp.float32)
        cnt[idxs] += 1
        inside |= m
    n_in = int(inside.sum())
    if not n_in:
        return None
    # Returned dense (identically 1.0 outside the layer): dispersion_post indexes straight by
    # entry index, with no second indirection entry[i] -> ofs[e] (see kernels/dispersion.cu)
    return {
        "dense": cp.asarray(factor, dtype=cp.float32),
        "f1": cp.asarray(fs[0]), "f2": cp.asarray(fs[1]),
        "f3": cp.asarray(fs[2]),
        "n_in": np.int32(n_in),
        # (n_entry,) bool: entries inside the layer cannot be deferred into update_e; the ade
        # path uses it on the host
        "inside": cp.asnumpy(inside) if _dev else inside,
    }


def _shell_segments(span, d: Dims, block=(128, 2, 1)):
    """Three-segment launch of the boundary shell (so no threads are started for the interior).
    Returns [(grid, block, (z0, z1))].

    ``span`` is the interior's ``(lo, hi)``. seg0/seg1 are the two segments outside the
    interior's z range (whole xy cross-sections); seg2 lies inside the interior's z range and
    covers only the xy border (the kernel's interior skip filters out the middle by itself). The
    union of the three segments is the boundary shell.
    """
    lo, hi = span
    nx, ny, nzp = d.shape
    segs = []
    for z0, z1 in ((0, lo[2]), (hi[2], nzp), (lo[2], hi[2])):
        if z1 <= z0:
            continue
        nzs = z1 - z0
        segs.append((((nzs + block[0] - 1) // block[0],
                      (ny + block[1] - 1) // block[1], nx),
                     block, (np.int32(z0), np.int32(z1))))
    return segs


def _pole_bucket_setup(disp, pabs):
    """Bucket the dispersion entries by pole count in a **stable ascending** order; the 2-pole
    bucket takes the float2 fast path.

    This is a pure permutation (the pole order inside an entry does not move) and it runs during
    initialization while P is all zeros, so it is bitwise unchanged (the same argument as for
    _ade_defer_setup). The slots of the 2-pole bucket are naturally contiguous (slot0 + 2j); an
    **orphan slot** is padded in front of the bucket to align slot0 to an even index (float2
    needs 8-byte alignment), and that orphan slot belongs to no entry range and is never read or
    written. cellc = comp<<30 | cell (packed on encoding, matched on decoding).

    Returns:
        ``{"lo","n2","n_pre","n_post","slot0","g2","g_pre","g_post"}``, or ``None`` (no 2-pole
        entry, or cell >= 2^30 so it does not fit into cellc; then no array is touched).
    """
    # The whole permutation is done on the GPU. The original moved hundreds of millions of
    # entries back to the host and expanded them with a numpy sort, measured at 8-9 s of the
    # preparation phase for Ring/GroupDelay.
    # Size threshold: GPU bucketing has a fixed overhead (allocation plus a dozen or so
    # synchronizations to fetch scalars). Measured: 2.2 s slower for Topo at 11.1 million
    # entries, 2.1 s faster for hBN at 20.5 million, 8 s faster for Ring/GroupDelay at 100
    # million, so the crossover is around 15 million.
    _n_ent = int(disp["comp"].size)
    _min = int(knobs.env("DISP_P2_GPU_MIN"))
    xp = (cp if (knobs.env("DISP_P2_GPU") != "0"
                 and _n_ent >= _min) else np)
    _dev = xp is cp

    ofs64 = (disp["ofs"] if _dev else cp.asnumpy(disp["ofs"])).astype(xp.int64)
    counts = xp.diff(ofs64)
    n2 = int((counts == 2).sum())
    if n2 == 0:
        return None
    ncnt = int(counts.size)
    comp = (disp["comp"] if _dev else cp.asnumpy(disp["comp"])).astype(xp.int64)
    cell = (disp["cell"] if _dev else cp.asnumpy(disp["cell"])).astype(xp.int64)
    if int(cell.max()) >= 1 << 30:  # cellc has only 30 bits for cell; above that, generic kernel
        return None
    cidx = disp["cidx"] if _dev else cp.asnumpy(disp["cidx"])
    # Stable ascending order: key = counts * ncnt + position. No two keys are equal, so the sort
    # result is unique and in the same order as numpy's kind="stable"; cupy.argsort itself gives
    # no stability guarantee.
    order = xp.argsort(counts * ncnt + xp.arange(ncnt, dtype=xp.int64))
    cs = counts[order]
    # cs is already ascending, so these two counts are exactly searchsorted's left and right
    # insertion points (cupy.searchsorted does not accept a Python scalar, so an equivalent form
    # is used that works on both sides)
    lo = int((cs < 2).sum())
    hi = int((cs <= 2).sum())

    new_ofs = xp.zeros(ncnt + 1, xp.int64)
    new_ofs[1:] = xp.cumsum(cs)
    pad = int(new_ofs[lo]) & 1
    if pad:
        new_ofs[lo:] += 1              # one orphan slot before the bucket, even alignment
    total = int(new_ofs[-1])

    new_cidx = xp.zeros(total, xp.int32)   # orphan slot has cidx=0 and is never read
    src_base = ofs64[order]
    for a, b in ((0, lo), (lo, hi), (hi, ncnt)):
        if a == b:
            continue
        seg_cs = cs[a:b]
        seg_total = int(seg_cs.sum())
        tgt0 = int(new_ofs[a])
        # Segment expansion: cupy.repeat does not accept an array of repeat counts, so a
        # cumulative sum plus searchsorted works out which segment each output position falls in
        # and then reads the value of that segment, which is equivalent to np.repeat(x, counts)
        cum = xp.cumsum(seg_cs)
        pos = xp.arange(seg_total, dtype=xp.int64)
        seg = xp.searchsorted(cum, pos, side="right")
        rel = pos - (new_ofs[a:b][seg] - tgt0)
        src = src_base[a:b][seg] + rel
        new_cidx[tgt0:tgt0 + seg_total] = cidx[src]

    comp_o, cell_o = comp[order], cell[order]
    disp["comp"] = cp.asarray(comp_o.astype(xp.int32))
    disp["cell"] = cp.asarray(cell_o.astype(xp.int32))
    disp["ofs"] = cp.asarray(new_ofs.astype(xp.int32))
    disp["cidx"] = cp.asarray(new_cidx)
    disp["cellc"] = cp.asarray(
        ((comp_o.astype(xp.uint32) << 30) | cell_o.astype(xp.uint32)
         ).astype(xp.int32))
    # P is all zeros at this point: just reallocate at the new slot count, nothing to move
    disp["p_re"] = cp.zeros(total, cp.float32)
    disp["p_im"] = cp.zeros(total, cp.float32)
    if pabs is not None:
        for key in ("dense", "f1", "f2", "f3"):
            pabs[key] = cp.asarray(
                (pabs[key] if _dev else cp.asnumpy(pabs[key]))[order])
        # inside is a host-side bool array: on the GPU path it goes up, is permuted, comes back
        pabs["inside"] = (cp.asnumpy(cp.asarray(pabs["inside"])[order])
                          if _dev else pabs["inside"][order])
    return {
        "lo": lo, "n2": n2, "slot0": np.int32(int(new_ofs[lo])),
        "g2": grid_1d(n2),
        "g_pre": grid_1d(max(lo, 1)),
        "g_post": grid_1d(max(ncnt - hi, 1)),
        "n_pre": lo, "n_post": ncnt - hi,
    }


def _ade_defer_setup(disp, pabs, d: Dims):
    """Move the dispersion entries that can be deferred into update_e to the front, and build the
    per-cell pole slice tables.

    **The deferral identity**: the E that ``dispersion_step`` reads at the end of step n is the
    same piece of memory that update_e reads at step n+1, with only the absorber-layer decay in
    between. So entries **outside** the layer can be moved wholesale into update_e and inlined
    (saving one write and one read of hist, plus the three streams comp/cell/ofs), while entries
    inside the layer still go through the old ``dispersion_step_dec``.

    The reordering is a **pure permutation** (the pole order inside an entry does not move), and
    P is still all zeros at initialization, so nothing has to be moved and the result is bitwise
    unchanged. ``None`` means the conditions are not met and everything falls back.

    Returns:
        ``{"estart", "ecount", "eb", "n_def"}``, or ``None``.
    """
    comp = cp.asnumpy(disp["comp"]).astype(np.int64)
    cell = cp.asnumpy(disp["cell"]).astype(np.int64)
    ofs = cp.asnumpy(disp["ofs"]).astype(np.int64)
    cidx = cp.asnumpy(disp["cidx"])
    n_entry = comp.size
    counts = np.diff(ofs)
    if counts.max() > 255:
        return None                      # ecount is uint8

    defer = np.ones(n_entry, dtype=bool) if pabs is None else ~pabs["inside"]
    n_def = int(defer.sum())
    if n_def == 0:
        return None

    order = np.lexsort((cell, comp))     # canonical (comp, cell) order
    d_ord = defer[order]
    perm = np.concatenate([order[d_ord], order[~d_ord]])

    new_counts = counts[perm]
    new_ofs = np.zeros(n_entry + 1, dtype=np.int64)
    np.cumsum(new_counts, out=new_ofs[1:])
    pole_src = _gather_poles(ofs, perm, new_counts, new_ofs)

    tables = _entry_tables(comp[perm[:n_def]], cell[perm[:n_def]],
                           new_ofs[:n_def], new_counts[:n_def], d)
    if tables is None:
        return None          # several entries share a (comp,cell): per-cell lookup breaks

    disp["comp"] = cp.asarray(np.ascontiguousarray(comp[perm], dtype=np.int32))
    disp["cell"] = cp.asarray(np.ascontiguousarray(cell[perm], dtype=np.int32))
    disp["ofs"] = cp.asarray(np.ascontiguousarray(new_ofs, dtype=np.int32))
    disp["cidx"] = cp.asarray(np.ascontiguousarray(cidx[pole_src], dtype=np.int32))
    if pabs is not None:
        for key in ("dense", "f1", "f2", "f3"):
            pabs[key] = cp.asarray(cp.asnumpy(pabs[key])[perm])
        pabs["inside"] = pabs["inside"][perm]

    estart, ecount, eb = tables
    return {"estart": estart, "ecount": ecount, "eb": eb, "n_def": n_def}


def _gather_poles(ofs, perm, new_counts, new_ofs) -> np.ndarray:
    """After the entries have been reordered by ``perm``, where new pole slot ``j`` sits in the
    old pole array (a pole-level gather index; the pole order inside an entry does not move)."""
    return (np.repeat(ofs[perm], new_counts)
            + np.arange(int(new_ofs[-1]), dtype=np.int64)
            - np.repeat(new_ofs[:-1], new_counts))


def _entry_tables(pc, pcell, ofs_n, cnt_n, d: Dims):
    """Per-cell pole slice tables ``(estart, ecount, eb)``: the poles of cell ``cell`` of
    component c start at ``estart[c][cell]`` and there are ``ecount[c][cell]`` of them (uint8);
    ``eb`` is the six-tuple bounding box of those cells. Returns None when several entries share
    one (comp, cell), because then direct per-cell lookup does not hold.

    Duplicate detection is an O(n) scatter (np.unique used to hash-deduplicate 81 million cells,
    measured at 50 s for GroupDelay). A duplicate means fewer landing spots than assignments.
    """
    ncell, ny, nzp = d.ncell, d.ny, d.nzp
    est = [np.zeros(ncell, np.uint32) for _ in range(3)]
    ecn = [np.zeros(ncell, np.uint8) for _ in range(3)]
    for c in range(3):
        sel = pc == c
        cc = pcell[sel]
        seen = np.zeros(ncell, dtype=bool)
        seen[cc] = True
        if int(seen.sum()) != cc.size:
            del seen
            return None
        del seen
        est[c][cc] = ofs_n[sel].astype(np.uint32)
        ecn[c][cc] = cnt_n[sel].astype(np.uint8)
    ii, jj, kk = unravel_cells(pcell, (ny, nzp))
    eb = tuple(np.int32(int(v)) for v in
               (ii.min(), ii.max() + 1, jj.min(), jj.max() + 1,
                kk.min(), kk.max() + 1))
    return (tuple(cp.asarray(x) for x in est),
            tuple(cp.asarray(x) for x in ecn), eb)


def _identity_span(mask: np.ndarray) -> tuple[int, int, bool]:
    """Search for the identity region: find the contiguous interior ``[lo, hi)`` in the 1D mask
    ``mask`` (True = that position gets an identity update). n//2 is the dividing line: it starts
    after the last non-identity position in the lower half and ends at the first non-identity
    position in the upper half. Returns ``(lo, hi, ok)``, with ``ok=False`` when the interval is
    empty or still contains a non-identity position."""
    n = int(mask.size)
    bad = np.where(~mask)[0]
    lo_i = int(bad[bad < n // 2].max()) + 1 if (bad < n // 2).any() else 0
    hi_i = int(bad[bad >= n // 2].min()) if (bad >= n // 2).any() else n
    ok = not (lo_i >= hi_i or not mask[lo_i:hi_i].all())
    return lo_i, hi_i, ok


def _shell_boxes(span, d: Dims, what: str) -> list:
    """The six boxes of the boundary shell, ``[(o, e), ...]``: they do not overlap and their
    union is the whole domain minus the interior ``span=(lo, hi)``; volume conservation is pinned
    down by an assertion."""
    lo, hi = span
    nx, ny, nzp = d.shape
    boxes = []

    def _add(o, e):
        if all(ee > oo for oo, ee in zip(o, e)):
            boxes.append((tuple(o), tuple(e)))

    _add((0, 0, 0), (lo[0], ny, nzp))
    _add((hi[0], 0, 0), (nx, ny, nzp))
    _add((lo[0], 0, 0), (hi[0], lo[1], nzp))
    _add((lo[0], hi[1], 0), (hi[0], ny, nzp))
    _add((lo[0], lo[1], 0), (hi[0], hi[1], lo[2]))
    _add((lo[0], lo[1], hi[2]), (hi[0], hi[1], nzp))
    n_int = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2])
    vol = n_int + sum((e[0] - o[0]) * (e[1] - o[1]) * (e[2] - o[2])
                      for o, e in boxes)
    assert vol == nx * ny * nzp, f"{what}: six boxes + interior do not conserve volume"
    return boxes


def _lean_side_out(span, boxes, d: Dims) -> dict:
    """Launch geometry of the interior ``span=(lo, hi)`` on one side (E or H): origin and end,
    the row-aligned z launch start (rounded down to a multiple of 32, with the kernel guarding
    k>=oz_real so that no layer is lost), the interior grid/block, the launch tuple of each of
    the six boxes, and the three-segment shell launch."""
    lo, hi = span
    oz = (lo[2] // 32) * 32
    span = [hi[i] - lo[i] for i in range(3)]
    return {
        "o": tuple(np.int32(v) for v in lo),
        "e": tuple(np.int32(v) for v in hi),
        "oz_launch": np.int32(oz), "oz_real": np.int32(lo[2]),
        "grid": ((hi[2] - oz + 127) // 128, span[1], span[0]),
        "block": (128, 1, 1),
        "boxes": [((((e[2] - o[2]) + 127) // 128, ((e[1] - o[1]) + 1) // 2, e[0] - o[0]),
                   (128, 2, 1),
                   tuple(np.int32(v) for v in (*o, *e)))
                  for o, e in boxes],
        "segs": _shell_segments(span, d),
    }


def _lean_setup(tabs: LeanTables, abs_slabs, disp, pabs, d: Dims):
    """Initialization of the interior / boundary-shell split.

    The interior is marked out from the data: for each axis the identity region is worked back
    from the coefficient arrays (a=0, b=1, invk=1, pec=1, mpv=1, prv=i-1), and the dispersive
    case then shrinks it by the face depth of the absorber slabs; the z start is aligned to 32
    (row alignment).
    Reordering of the dispersion entries: the interior prefix goes into canonical (comp,cell)
    order (a pure permutation, so P is bitwise unchanged) while the boundary-shell suffix keeps
    its original relative order. If any threshold or assertion fails, it returns None and
    everything falls back.

    Three parts: :func:`_lean_geometry` (the E-side interior plus the thresholds), the H-side
    interior (marked out independently, optional), and :func:`_lean_disp_partition` (grouping of
    the dispersion entries).
    """
    def _rej(why):
        # The reason for the fallback stays in the f-string at the call site, where it is visible
        # when reading the code; it is no longer printed
        return None

    nx, ny, nzp = d.shape
    dims = d.shape
    ncell = d.ncell
    geo = _lean_geometry(tabs, abs_slabs, disp, d, _rej)
    if geo is None:
        return None
    lo, hi, boxes, n_int_cells = geo
    # ---- H-side interior (marked out independently: identity = a=0,b=1,invk=1,mnx=1,nxt=i+1) --
    h_lo, h_hi = [0, 0, 0], [nx, ny, nzp]
    h_ok = True
    for ax in range(3):
        a = cp.asnumpy(tabs.pml_h[ax][0]).astype(np.float64)
        b = cp.asnumpy(tabs.pml_h[ax][1]).astype(np.float64)
        kk = cp.asnumpy(tabs.pml_h[ax][2]).astype(np.float64)
        nxv = cp.asnumpy(tabs.nxt[ax]).astype(np.int64)
        mv = cp.asnumpy(tabs.mnx[ax]).astype(np.float64)
        n = dims[ax]
        m = ((a == 0.0) & (b == 1.0) & (kk == 1.0) & (mv == 1.0)
             & (nxv == np.arange(n) + 1))
        lo_i, hi_i, ok = _identity_span(m)
        if not ok:
            h_ok = False
            break
        h_lo[ax], h_hi[ax] = lo_i, hi_i
    h_out = None
    if h_ok:
        h_span = [h_hi[i] - h_lo[i] for i in range(3)]
        h_cells = h_span[0] * h_span[1] * h_span[2]
        if min(h_span) > 0 and h_span[2] >= 32 and h_cells / ncell >= 0.50:
            h_boxes = _shell_boxes((h_lo, h_hi), d, "P25c: H")
            h_out = _lean_side_out((h_lo, h_hi), h_boxes, d)
    out = {"h": h_out, **_lean_side_out((lo, hi), boxes, d),
           "n_int": 0, "n_skin": 0, "frac": n_int_cells / ncell}
    if disp is None:
        return out
    return _lean_disp_partition(out, disp, pabs, (lo, hi), d)


def _lean_geometry(tabs: LeanTables, abs_slabs, disp, d: Dims, _rej):
    """E-side interior: the identity region per axis -> shrunk by the slabs in the dispersive
    case -> the span, fraction and domain-size thresholds -> the six boundary-shell boxes.
    Returns ``(lo, hi, boxes, n_int_cells)``, or None when a threshold is not met."""
    dims = d.shape
    nx, ny, nzp = dims
    lo, hi = [0, 0, 0], [nx, ny, nzp]
    for ax in range(3):
        a = cp.asnumpy(tabs.pml_e[ax][0]).astype(np.float64)
        b = cp.asnumpy(tabs.pml_e[ax][1]).astype(np.float64)
        kk = cp.asnumpy(tabs.pml_e[ax][2]).astype(np.float64)
        pe = cp.asnumpy(tabs.pec[ax]).astype(np.float64)
        pv = cp.asnumpy(tabs.prv[ax]).astype(np.int64)
        mv = cp.asnumpy(tabs.mpv[ax]).astype(np.float64)
        n = dims[ax]
        m = ((a == 0.0) & (b == 1.0) & (kk == 1.0) & (pe == 1.0)
             & (mv == 1.0) & (pv == np.arange(n) - 1))
        lo_i, hi_i, ok = _identity_span(m)
        if not ok:
            return _rej(f"axis {ax}: identity region not contiguous l={lo_i} h={hi_i}")
        lo[ax], hi[ax] = lo_i, hi_i
    if disp is not None:
        # Dispersive case: the slabs (the absorbing decay regions) all go into the boundary
        # shell, since entries inside the shell cannot be deferred
        for sl in abs_slabs:
            ax = int(sl["axis"])
            g0 = int(sl["org"][ax])
            nl = int(sl["box"][ax])
            if g0 <= lo[ax]:
                lo[ax] = max(lo[ax], g0 + nl)
            if g0 + nl >= hi[ax]:
                hi[ax] = min(hi[ax], g0)
            if g0 < hi[ax] and g0 + nl > lo[ax]:
                return _rej(f"slab in mid-domain ax={ax} g0={g0} nl={nl} lo={lo} hi={hi}")
    span = [hi[i] - lo[i] for i in range(3)]
    ncell = d.ncell
    _min_z = int(knobs.env("LEAN_MINZ"))
    if min(span) <= 0 or span[2] < _min_z:
        return _rej(f"span too small span={span} lo={lo} hi={hi} (z needs >={_min_z})")
    n_int_cells = span[0] * span[1] * span[2]
    # Thresholds (calibrated by measurement): there is a net gain only when the interior is at
    # least 70% of the domain and the domain is at least 8M cells. The single boundary-shell
    # kernel still has to launch over the whole domain once, and when the interior fraction is
    # not high enough the machine saved does not pay for that pass; on a small domain (after
    # folding) a step is already down at a few hundred microseconds, where the extra launches of
    # the split are a noticeable share.
    _min_frac = float(knobs.env("LEAN_MINFRAC"))
    if n_int_cells / ncell < _min_frac:
        return _rej(f"interior fraction {n_int_cells/ncell:.0%} < {_min_frac:.0%} "
                    f"lo={lo} hi={hi}")
    if ncell < 8_000_000:
        return _rej(f"domain too small {ncell:,} < 8M (split overhead too large a share)")
    boxes = _shell_boxes((lo, hi), d, "P25b")
    return lo, hi, boxes, n_int_cells


def _lean_disp_partition(out, disp, pabs, span, d: Dims):
    """Grouping of the dispersion entries: the interior prefix goes into canonical (comp, cell)
    order while the boundary-shell suffix keeps its original relative order; then the pole-level
    gather, the matching reordering of the pabs factors, and the per-cell table of the interior.
    Returns None if any threshold is not met. ``span`` is the interior's ``(lo, hi)``."""
    lo, hi = span
    ncell, ny, nzp = d.ncell, d.ny, d.nzp
    comp = cp.asnumpy(disp["comp"])
    cell = cp.asnumpy(disp["cell"])
    ofs = cp.asnumpy(disp["ofs"]).astype(np.int64)
    cidx = cp.asnumpy(disp["cidx"])
    ii, jj, kk2 = unravel_cells(cell, (ny, nzp))
    is_int = ((ii >= lo[0]) & (ii < hi[0]) & (jj >= lo[1]) & (jj < hi[1])
              & (kk2 >= lo[2]) & (kk2 < hi[2]))
    idx_int = np.where(is_int)[0]
    idx_skin = np.where(~is_int)[0]
    # A single-key GPU argsort replaces the two-key lexsort (measured 8.3 s for GroupDelay).
    # comp is in {0,1,2} and cell < ncell, so key = comp*ncell + cell is unique and keeps the
    # order.
    if idx_int.size >= 1_000_000:
        _k = cp.asarray(comp[idx_int].astype(np.int64) * ncell
                        + cell[idx_int].astype(np.int64))
        order_int = idx_int[cp.asnumpy(cp.argsort(_k, kind="stable"))]
        del _k
    else:
        order_int = idx_int[np.lexsort((cell[idx_int], comp[idx_int]))]
    perm = np.concatenate([order_int, idx_skin])
    n_int, n_skin = idx_int.size, idx_skin.size
    counts = ofs[1:] - ofs[:-1]
    new_counts = counts[perm]
    new_ofs = np.zeros(perm.size + 1, dtype=np.int64)
    np.cumsum(new_counts, out=new_ofs[1:])
    # Pole-level gather (cidx; p starts from zero, so nothing has to be migrated)
    pole_src = _gather_poles(ofs, perm, new_counts, new_ofs)
    # estart/ecount per-cell tables (the interior prefix; (comp,cell) has to be unique). The two
    # fallback checks sit **before any reordering** (the same order as in _ade_defer_setup):
    # previously the four pabs tables were permuted by perm and written back before the fallback,
    # after which the entry order of pabs and disp no longer matched, with no error raised
    # (2026-09-10 moved the statements only; the out[...] assignment stays where it was).
    tables = None
    if n_int:
        cnt_int = new_counts[:n_int]
        if cnt_int.max() > 255:
            return None
        tables = _entry_tables(comp[perm[:n_int]], cell[perm[:n_int]],
                               new_ofs[:n_int], cnt_int, d)
        if tables is None:
            return None      # several entries share a (comp,cell): layout assumption breaks
    # Interior entries must have identity decay (dec=1, f=1), otherwise fall back
    if pabs is not None and n_int:
        dec = cp.asnumpy(pabs["dense"])
        f1 = cp.asnumpy(pabs["f1"])
        f2 = cp.asnumpy(pabs["f2"])
        f3 = cp.asnumpy(pabs["f3"])
        pi = perm[:n_int]
        if not ((dec[pi] == 1.0).all() and (f1[pi] == 1.0).all()
                and (f2[pi] == 1.0).all() and (f3[pi] == 1.0).all()):
            return None
        pabs["dense"] = cp.asarray(dec[perm])
        pabs["f1"] = cp.asarray(f1[perm])
        pabs["f2"] = cp.asarray(f2[perm])
        pabs["f3"] = cp.asarray(f3[perm])
    elif pabs is not None:
        for key in ("dense", "f1", "f2", "f3"):
            pabs[key] = cp.asarray(cp.asnumpy(pabs[key])[perm])
    if n_int:
        out["estart"], out["ecount"], out["eb"] = tables
    # Entry-level arrays are replaced in place (the kernel of the old full path handles entries
    # independently, so a permutation does not change the arithmetic)
    disp["comp"] = cp.asarray(np.ascontiguousarray(comp[perm], dtype=np.int32))
    disp["cell"] = cp.asarray(np.ascontiguousarray(cell[perm]))
    disp["ofs"] = cp.asarray(np.ascontiguousarray(new_ofs, dtype=np.int32))
    disp["cidx"] = cp.asarray(np.ascontiguousarray(cidx[pole_src], dtype=np.int32))
    out["n_int"], out["n_skin"] = int(n_int), int(n_skin)
    return out
