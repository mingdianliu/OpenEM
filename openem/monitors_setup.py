# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Monitor tables: accumulation buffers for the four monitor types, the DFT phase table, launch
geometry, shared snapshots and the batched sampling table.
(Moved verbatim out of setup_tables.py and readout.py, 2026-09; the readout itself is still in
readout.py.)
"""

from __future__ import annotations

import cupy as cp
import numpy as np

from openem import apodization as apod_mod
from openem import flux as flux_mod
from openem import knobs
from openem.device import grid_1d, grid_2d, grid_3d
from openem.grid import transverse_axes
from openem.model import Scene
from openem.pitch import Dims, _f32, _f64



_DFT_TMOD = 80   # rows of the phase ring = upper bound on the batch size; 80 divides the
                 # shutoff interval of 400, and the larger the batch the thinner the fp64
                 # accumulator traffic spreads (on Multipole with T=16 it is still the big item)


#: The name template of an n2f surface monitor is ``{name}__n2f_{axis}{sign}``
#: (scene/projection.py). Near-to-far only uses the tangential E/H (the equivalent surface
#: currents J=n x H and M=-n x E of Tidy3D's FieldProjector) and never reads the normal
#: components, so only four are accumulated, saving a third of the device memory and a third of
#: the DFT work.
_N2F_TAN = {0: (1, 2, 4, 5), 1: (0, 2, 3, 5), 2: (0, 1, 3, 4)}


def _flux_time_buffers(sc: Scene) -> list:
    """Buffers and geometry of the FluxTimeMonitor.

    Every step reduces the instantaneous Poynting vector over the colocated plane to a single
    scalar on the GPU, so only a ``(num_slots,)`` sequence is stored. The geometry
    (interpolation weights, signed trapezoidal area elements) all lives in
    ``flux.flux_time_geometry``; sampling and the tests share ``flux.colocated_component``
    (same expression under numpy and cupy). The time average of H keeps the two-shot structure
    of tmons, and the second shot has to be paired with "the colocated E plane of the previous
    step", so two small planes are kept.
    """
    ftmons = []
    for m in sc.flux_time_monitors:
        g = flux_mod.flux_time_geometry(sc.grid, m)
        for key in ("wlo", "whi"):
            g[key] = [cp.asarray(w) for w in g[key]]
        for _p in (0, 1):                       # parity weights of the folded axis
            if g["fold_lo"][_p]:
                g["idx_lo"][_p] = cp.asarray(g["idx_lo"][_p])
                g["wlo_c"][_p] = {k: cp.asarray(v)
                                  for k, v in g["wlo_c"][_p].items()}
        g["dS"] = cp.asarray(g["dS"])
        ftmons.append({
            "name": m.name,
            "geom": g,
            "beg": m.step_begin,
            "end": m.step_end,
            "iv": m.interval,
            "slots": m.num_slots,
            "buf": cp.zeros(m.num_slots, cp.float64),
            "e1": None,
            "e2": None,
            "filled": 0,
        })
    return ftmons


def _monitor_buffers(sc: Scene, d: Dims, n_total: int,
                     batch_cap: int = 0) -> tuple[list, list, list, list]:
    """Accumulation buffers of the four monitor types, ``(flux, frequency-domain field,
    time-domain field, time-domain flux)``.

    The shapes and index boxes are pinned down here; the main loop only writes into them.

    The apodization window is budgeted per monitor into two per-step scalar tables: E at
    ``t=(n+1)·dt`` and H at ``(n+0.5)·dt``, each windowed at its own real sampling instants
    (the same convention as the DFT phase). Without apodization it is identically 1.0, and
    multiplying by 1.0 in the kernel is exactly the identity, so there is one code path.
    """
    t_eh = ((np.arange(n_total, dtype=np.float64) + 1.0) * sc.dt,
            (np.arange(n_total, dtype=np.float64) + 0.5) * sc.dt)
    dftab = _dft_table(sc, n_total)

    # ---- monitors ----
    mons, fmons = [], []
    _off = 0
    for m in sc.flux_monitors:
        mons.append(_flux_mon_entry(m, d, t_eh, _off, batch_cap))
        _off += m.freqs.size
    for m in sc.field_monitors:
        fmons.append(_field_mon_entry(m, t_eh, _off, batch_cap))
        _off += m.freqs.size
    tmons = [_tmon_entry(m) for m in sc.field_time_monitors]

    return mons, fmons, tmons, _flux_time_buffers(sc), dftab


def _dft_table(sc: Scene, n_total: int) -> dict:
    """Phase table (one table shared by every DFT monitor).

    cos/sin(2πft) only depends on (f, t), not on the monitor. The frequency points of all DFT
    monitors are concatenated into one list so the table kernel is launched once per step, and
    each monitor takes a view of its own segment (the ``f0`` offset; the order is flux monitors
    first, then field monitors). The "one launch per monitor" version measured 1.76 times
    slower on scenes with many monitors and few frequency points.
    """
    _dft_mons = list(sc.flux_monitors) + list(sc.field_monitors)
    _freq_all = (np.concatenate([np.asarray(m.freqs, np.float64) for m in _dft_mons])
                 if _dft_mons else np.zeros(0, np.float64))
    _n_dft = int(_freq_all.size)
    return {
        "n": np.int32(_n_dft),
        "freqs": _f64(_freq_all),
        # ring table (_DFT_TMOD, nf, 4), flattened; row = step % _DFT_TMOD
        "tab": cp.zeros(_DFT_TMOD * max(_n_dft, 1) * 4, cp.float64),
        # fp32 copy for the batched flush kernel; the fp64 table stays for the legacy per-step path
        "tab32": cp.zeros(_DFT_TMOD * max(_n_dft, 1) * 4, cp.float32),
        "tmod": np.int32(_DFT_TMOD),
        "launch": grid_1d(max(_n_dft, 1)),
    }


def _pick_T(points: int, nf: int, batch_cap: int) -> int:
    """Pick the batch depth. On DFT-heavy cases the 128 B the fp64 accumulator reads and writes
    per (point, frequency) per step is the big item (90.6% of the step time on Multipole), so
    the raw f32 fields are snapshotted for T steps instead and settled every T steps in step
    order (the sequence of additions is unchanged, so the result is bitwise identical). On
    small monitors the launch overhead outweighs the gain and the original path is taken."""
    _batch_on = knobs.env("DFT_BATCH") != "0"
    _batch_min = int(knobs.env("DFT_BATCH_MIN"))
    _ring_mb = float(knobs.env("DFT_RING_MB"))
    if not _batch_on or batch_cap < 2 or points * nf < _batch_min:
        return 0
    T = min(batch_cap, 32)          # smem limit (tile>=64 needs T<=32)
    while T >= 2 and T * 6 * 4 * points > _ring_mb * 1e6:
        T //= 2
    return T if T >= 2 else 0


def _flush_launch(points: int, nf: int):
    """P20c flush launch: grid=(frequency group, point tile_y, point tile_z), 256-thread blocks."""
    tile = 256
    nb = (points + tile - 1) // tile
    return (((nf + 7) // 8, min(nb, 65535), (nb + 65534) // 65535),
            (tile, 1, 1))


def _use_sm(points: int, nf: int) -> bool:
    """The frequency-group version is only used when both the point and the frequency count are
    large enough (at small sizes the lost parallelism costs +11%). The thresholds can be pushed
    to 0 through the environment, which is how the guard tests force this path."""
    return (nf >= int(knobs.env("DFT_SM_NF"))
            and points >= int(knobs.env("DFT_SM_PTS")))


def _flux_mon_entry(m, d: Dims, t_eh, f0: int, batch_cap: int) -> dict:
    """Buffer entry of one FluxMonitor (planar DFT, ``(4, nf, n1, n2)``).
    ``t_eh`` is ``(t_e_all, t_h_all)``, the two sampling-instant tables with E at ``(n+1)·dt``
    and H at ``(n+0.5)·dt``, each apodized at its own instants."""
    nx, ny, nz = d.nx, d.ny, d.nz
    t_e_all, t_h_all = t_eh
    nf = m.freqs.size
    # the two transverse axes in ascending order, matching (t1, t2) in the kernel
    t1, t2 = transverse_axes(m.axis)
    n1, n2 = (nx, ny, nz)[t1], (nx, ny, nz)[t2]
    # frequency points sit on grid.z (see the notes on accumulate_dft), capped by CUDA at 65535
    if nf > 65535:
        raise NotImplementedError(
            f"monitor {m.name} has {nf} frequency points, over the grid.z limit of 65535"
        )
    g2, b2 = grid_2d(n1, n2)
    # the window table is computed once and reused everywhere (same pure function, same input,
    # so the values are bitwise identical)
    we = apod_mod.window(t_e_all, m.apodization)
    wh = apod_mod.window(t_h_all, m.apodization)
    T = _pick_T(n1 * n2, nf, batch_cap)
    return {
        "name": m.name,
        "km": m.plane_index,
        "axis": m.axis,
        "nf": nf,
        "launch": ((g2[0], g2[1], nf), (b2[0], b2[1], 1)),
        "launch_fl": _flush_launch(n1 * n2, nf) if T else None,
        "sm": _use_sm(n1 * n2, nf),
        "launch_pf": ((g2[0], g2[1], nf), (b2[0], b2[1], 1)),
        "launch_snap": ((g2[0], g2[1], 1), (b2[0], b2[1], 1)),
        "f0": np.int32(f0),
        "n1": np.int32(n1), "n2": np.int32(n2),
        "T": np.int32(T),
        "snap": cp.zeros(T * 6 * n1 * n2, cp.float32) if T else None,
        "freqs": _f64(m.freqs),
        "win_e": we,
        "win_h": wh,
        # device copy of the window table; win[*step] in the kernel is bitwise the same as
        # taking the scalar on the host step by step, as it used to
        "win_e_dev": _f64(we),
        "win_h_dev": _f64(wh),
        "win_e_dev32": _f32(we),
        "win_h_dev32": _f32(wh),
        "re": cp.zeros((4, nf, n1, n2), cp.float64),
        "im": cp.zeros((4, nf, n1, n2), cp.float64),
    }


def _field_mon_entry(m, t_eh, f0: int, batch_cap: int) -> dict:
    """Buffer entry of one FieldMonitor (box DFT, ``(ncomp, nf, ni, nj, nk)``; an n2f surface
    accumulates only the four tangential components, see :func:`_n2f_cmap`). ``t_eh`` is as in
    :func:`_flux_mon_entry`."""
    t_e_all, t_h_all = t_eh
    nf = m.freqs.size
    ni, nj, nk = m.box
    # (j,k) flattened onto grid.x, frequencies onto grid.y, i onto grid.z (accumulate_dft_box)
    if nf > 65535 or int(ni) > 65535:
        raise NotImplementedError(
            f"monitor {m.name}: nf={nf}, ni={int(ni)}, over the CUDA grid dimension limit of 65535"
        )
    _nthr = 256
    _gjk = (int(nj) * int(nk) + _nthr - 1) // _nthr
    cells = int(ni) * int(nj) * int(nk)
    T = _pick_T(cells, nf, batch_cap)
    cmap = _n2f_cmap(m.name, bool(T))
    we = apod_mod.window(t_e_all, m.apodization)
    wh = apod_mod.window(t_h_all, m.apodization)
    return {
        "name": m.name,
        "origin": tuple(np.int32(v) for v in m.origin),
        "box": (np.int32(ni), np.int32(nj), np.int32(nk)),
        "nf": np.int32(nf),
        "launch": ((_gjk, nf, int(ni)), (_nthr, 1, 1)),
        "launch_snap": ((_gjk, 1, int(ni)), (_nthr, 1, 1)),
        # flush threads = (grid point, 8 frequency groups), occupancy first
        "launch_flush": _flush_launch(cells, nf) if T else None,
        "sm": _use_sm(cells, nf) and len(cmap) == 6,
        "launch_pf": (((cells + _nthr - 1) // _nthr, nf, 1), (_nthr, 1, 1)),
        "f0": np.int32(f0),
        "cells": np.int32(cells),
        "T": np.int32(T),
        "snap": cp.zeros(T * len(cmap) * cells, cp.float32) if T else None,
        "freqs": _f64(m.freqs),
        "win_e": we,
        "win_h": wh,
        "win_e_dev": _f64(we),
        "win_h_dev": _f64(wh),
        "win_e_dev32": _f32(we),
        "win_h_dev32": _f32(wh),
        "cmap": cmap,
        "re": cp.zeros((len(cmap), nf, ni, nj, nk), cp.float64),
        "im": cp.zeros((len(cmap), nf, ni, nj, nk), cp.float64),
    }


def _tmon_entry(m, complex_buf: bool = False) -> dict:
    """One buffer entry of a FieldTimeMonitor: the real path gets a single ``buf``, the complex
    Bloch path gets ``buf_re`` / ``buf_im`` (solver._run_bloch), and every other field is the
    same."""
    ni, nj, nk = m.box
    out = {
        "name": m.name,
        "beg": m.step_begin,
        "end": m.step_end,
        "iv": m.interval,
        "slots": m.num_slots,
        "nc": np.int32(len(m.comps)),
        "comps": cp.asarray(np.asarray(m.comps, dtype=np.int32)),
        # E takes full weight on the whole step; H is split into two shots of half each, see
        # the notes in kernels
        "w_a": _f32([1.0 if c < 3 else 0.5 for c in m.comps]),
        "w_b": _f32([0.0 if c < 3 else 0.5 for c in m.comps]),
        "origin": tuple(np.int32(v) for v in m.origin),
        "box": (np.int32(ni), np.int32(nj), np.int32(nk)),
    }
    shape = (len(m.comps), m.num_slots, ni, nj, nk)
    if complex_buf:
        out["buf_re"] = cp.zeros(shape, cp.float32)
        out["buf_im"] = cp.zeros(shape, cp.float32)
    else:
        out["buf"] = cp.zeros(shape, cp.float32)
    out["launch"] = grid_3d(ni, nj, nk)
    out["filled"] = 0
    return out


def _n2f_cmap(name: str, batched: bool = True) -> tuple:
    """The physical components this field monitor has to accumulate (0..2=E, 3..5=H). Anything
    that is not an n2f surface => all six.

    With batched=False (_pick_T returned 0, so accumulate_dft_box writes directly) six must be
    returned: that kernel writes re/im in a fixed 6-component layout, and 4 components would
    run out of bounds. Small monitors take little device memory anyway, so skipping the
    optimization there costs nothing.
    """
    if not batched or knobs.env("N2F_TANGENTIAL") == "0":
        return (0, 1, 2, 3, 4, 5)
    i = name.find("__n2f_")
    if i < 0:
        return (0, 1, 2, 3, 4, 5)
    ax = "xyz".find(name[i + 6:i + 7])
    return _N2F_TAN[ax] if ax >= 0 else (0, 1, 2, 3, 4, 5)


def _share_snapshots(fmons, verbose) -> None:
    """Field monitors on the same box share one snapshot buffer.

    On Near2Far the near_field_* and far_field__n2f_* monitors form six pairs on the same box,
    and each of them used to snapshot separately every step. The first monitor on a box becomes
    the group leader, the rest point their ``snap`` at the leader's buffer and set
    ``skip_snap``, so a box is snapshotted only once per step.
    """
    owner = {}
    for m in fmons:
        if m.get("snap") is None:
            continue
        # since P33 the snap layout is (T, ncomp, cells); different component tables cannot share
        key = (m["origin"], m["box"], int(m["T"]),
               tuple(m.get("cmap", (0, 1, 2, 3, 4, 5))))
        own = owner.get(key)
        if own is None:
            owner[key] = m
            m["skip_snap"] = False
        else:
            m["snap"] = own["snap"]        # share the same buffer
            m["skip_snap"] = True          # this monitor no longer snapshots
    n = sum(1 for m in fmons if m.get("skip_snap"))
    if n and verbose:
        print(f"  [P26] field monitors share snapshots per box: of {len(fmons)}, {n} reuse "
              f"the group leader's buffer")


def _fmon_launch_tables(fmons, verbose) -> None:
    """Fill in the launch tables of every field monitor (modifies the dicts inside ``fmons`` in
    place).

    The launch geometry and shared-memory byte count of the smem tiled version (``_pick_T``
    guarantees T<=32, so smem is <=12 KB).

    A tangential-only surface monitor accumulates just four components, so it carries a
    component mapping table and launches the ``_sub`` variant of the kernel.
    """
    for _m in fmons:
        if _m.get("T"):
            _m["launch_sm2"] = ((( int(_m["cells"]) + 15) // 16, 1, 1),
                                (16, 16, 1))
            _m["smem2"] = int(_m["T"]) * 6 * 16 * 4   # T<=32 (_pick_T) -> <=12KB
    for _m in fmons:
        _cm = _m.get("cmap", (0, 1, 2, 3, 4, 5))
        _m["ncomp"] = np.int32(len(_cm))
        _m["cmap_dev"] = cp.asarray(np.asarray(_cm, dtype=np.int32))
    _n_tan = sum(1 for _m in fmons if int(_m["ncomp"]) < 6)
    if _n_tan and verbose:
        print(f"  [P33] n2f surface monitors accumulate tangential only: {_n_tan} of them, "
              f"4/6 components")


def _tmon_batch(tmons, verbose):
    """Launch table for batched sampling of the time monitors; returns None with fewer than two
    monitors.

    13 monitors x 2 phases = 26 launches per step, and of the 530 µs measured, 444 µs was pure
    launch overhead (the point monitors are 8-12 cells each). Packed into one launch.
    """
    if not (tmons and len(tmons) >= 2
            and knobs.env("TMON_BATCH") != "0"):
        return None

    def _ptr(arrs):
        return cp.asarray(np.array([int(a.data.ptr) for a in arrs],
                                   dtype=np.uint64))

    def _i32(vals):
        return cp.asarray(np.array(vals, dtype=np.int32))

    cells = [int(tm["box"][0]) * int(tm["box"][1]) * int(tm["box"][2])
             for tm in tmons]
    nthr = 128
    if verbose:
        print(f"  [P27] batched time monitor sampling: {len(tmons)} monitors x 2 phases "
              f"-> 2 launches (largest box {max(cells):,} cells)")
    return {
        "outs": _ptr([tm["buf"] for tm in tmons]),
        "comps": _ptr([tm["comps"] for tm in tmons]),
        "w_a": _ptr([tm["w_a"] for tm in tmons]),
        "w_b": _ptr([tm["w_b"] for tm in tmons]),
        "nc": _i32([int(tm["nc"]) for tm in tmons]),
        "beg": _i32([int(tm["beg"]) for tm in tmons]),
        "iv": _i32([int(tm["iv"]) for tm in tmons]),
        "end": _i32([int(tm["end"]) for tm in tmons]),
        "slots": _i32([int(tm["slots"]) for tm in tmons]),
        "i0": _i32([int(tm["origin"][0]) for tm in tmons]),
        "j0": _i32([int(tm["origin"][1]) for tm in tmons]),
        "k0": _i32([int(tm["origin"][2]) for tm in tmons]),
        "ni": _i32([int(tm["box"][0]) for tm in tmons]),
        "nj": _i32([int(tm["box"][1]) for tm in tmons]),
        "nk": _i32([int(tm["box"][2]) for tm in tmons]),
        "cells": _i32(cells),
        "n": np.int32(len(tmons)),
        "launch": (((max(cells) + nthr - 1) // nthr, len(tmons), 1),
                   (nthr, 1, 1)),
    }


def _phase_table(te: np.ndarray, th: np.ndarray, freqs, apod) -> np.ndarray:
    """Phase table of the per-step DFT on the Bloch path, ``(n_total, nf, 4)`` =
    ``[cosωt_e·win, sinωt_e·win, cosωt_h·win, sinωt_h·win]``, float64.
    Order: ``np.outer`` first, then cos/sin, then multiply by the apodization window (the same
    as the two inlined versions this replaced)."""
    w_e = apod_mod.window(te, apod)
    w_h = apod_mod.window(th, apod)
    om = 2.0 * np.pi * np.asarray(freqs, dtype=np.float64)   # (nf,)
    ph = np.empty((te.size, om.size, 4), dtype=np.float64)
    ang_e = np.outer(te, om)   # (n_total, nf)
    ang_h = np.outer(th, om)
    ph[:, :, 0] = np.cos(ang_e) * w_e[:, None]
    ph[:, :, 1] = np.sin(ang_e) * w_e[:, None]
    ph[:, :, 2] = np.cos(ang_h) * w_h[:, None]
    ph[:, :, 3] = np.sin(ang_h) * w_h[:, None]
    return ph


def _bloch_field_monitors(sc: Scene, t_eh) -> list[dict]:
    """Frequency-domain accumulation of FieldMonitors (complex field phasors: xy/xz/yz_freq and
    the like at oblique incidence). On the Bloch path every monitor gets one per-step phase
    table (:func:`_phase_table`) plus re/im float64 accumulators of shape ``(6, nf, box)``.
    ``t_eh`` is ``(te, th)``, the sampling-instant tables of E and H."""
    te, th = t_eh
    fmons_c = []
    for m in sc.field_monitors:
        ni, nj, nk = (int(v) for v in m.box)
        nf = int(m.freqs.size)
        ph = _phase_table(te, th, m.freqs, m.apodization)
        _gjk = (nj * nk + 255) // 256
        fmons_c.append({
            "name": m.name,
            "origin": tuple(np.int32(v) for v in m.origin),
            "box": (np.int32(ni), np.int32(nj), np.int32(nk)),
            "nf": np.int32(nf),
            "ph": cp.asarray(np.ascontiguousarray(ph)),
            "re": cp.zeros((6, nf, ni, nj, nk), cp.float64),
            "im": cp.zeros((6, nf, ni, nj, nk), cp.float64),
            "launch": ((_gjk, nf, ni), (256, 1, 1)),
        })
    return fmons_c


def _bloch_flux_monitors(sc: Scene, t_eh) -> list[dict]:
    """Complex colocated flux DFT of FluxMonitors (T/R at oblique incidence, the layer under
    DiffractionMonitor). On the Bloch path every monitor gets re/im float64 accumulators of
    shape ``(4, nf, n1, n2)``. ``t_eh`` is as in :func:`_bloch_field_monitors`."""
    te, th = t_eh
    nx, ny, nz = sc.shape
    mons_c = []
    for m in sc.flux_monitors:
        nf = int(m.freqs.size)
        t1a, t2a = transverse_axes(m.axis)
        n1 = (nx, ny, nz)[t1a]
        n2 = (nx, ny, nz)[t2a]
        ph = _phase_table(te, th, m.freqs, m.apodization)
        g2, b2 = grid_2d(n1, n2)
        mons_c.append({
            "name": m.name, "km": np.int32(m.plane_index), "axis": np.int32(m.axis),
            "nf": np.int32(nf), "n1": n1, "n2": n2, "freqs": _f64(m.freqs),
            "ph": cp.asarray(np.ascontiguousarray(ph)),
            "re": cp.zeros((4, nf, n1, n2), cp.float64),
            "im": cp.zeros((4, nf, n1, n2), cp.float64),
            "launch": ((g2[0], g2[1], nf), (b2[0], b2[1], 1)),
        })
    return mons_c
