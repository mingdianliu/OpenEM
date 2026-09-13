# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Per-step monitor sampling, batched flushing and final readout.

Table construction lives in monitors_setup.py.
"""

from __future__ import annotations

from typing import NamedTuple

import cupy as cp
import numpy as np

from openem import flux as flux_mod
from openem.model import COMPONENTS
from openem.results import Phasors


def _sample_flux_time(m_step: int, ftmons, e_arr, h_arr, nz) -> None:
    """Time-domain flux sampling; the two-shot structure matches :func:`_sample_time`.

    Instantaneous flux is linear in H, so the time average splits into two terms::

        s(m) = E^m x H_bar^m = 0.5*E^m x H^{m-1/2} + 0.5*E^m x H^{m+1/2}

    The first shot, at step m where E^m and H^{m-1/2} are in hand, computes the first term and
    stores the colocated E planes into e1/e2 while it is there. The second shot, at step m+1 where
    H^{m+1/2} is in hand, uses those stored E planes for the second term. Both shots fire under
    exactly the same sampling condition, so e1/e2 can never be read stale. But the **second shot
    has to go first**, because it reads the E planes stored on the previous step: with interval=1
    both shots happen in the same call, and the other order would have the first shot overwrite
    e1/e2 before the second could read them.

    ``filled`` only advances on the second shot, so if early termination lands between the two, the
    last slot is dropped rather than kept half-formed.

    The spatial convention -- colocated integration over each primal plane, flux composed linearly
    by frac plus trapezoidal area elements -- lives entirely in ``flux.colocated_component`` and
    ``flux_time_geometry``, the same implementation the tests use as their host reference
    (``flux.flux_time_sample``). Composition, colocation and time averaging are all linear, so
    their order can be exchanged.
    """
    for fm in ftmons:
        g = fm["geom"]
        t1, t2 = g["tax"]
        for phase, m_ref in ((1, m_step - 1), (0, m_step)):
            off = m_ref - fm["beg"]
            if off < 0 or m_ref >= fm["end"] or off % fm["iv"]:
                continue
            slot = off // fm["iv"]
            if slot >= fm["slots"]:
                continue
            if phase == 0:
                fm["e1"] = [flux_mod.colocated_component(
                                e_arr[t1][..., :nz], g, "Eu", p)
                            for p in range(len(g["planes"]))]
                fm["e2"] = [flux_mod.colocated_component(
                                e_arr[t2][..., :nz], g, "Ev", p)
                            for p in range(len(g["planes"]))]
            s = 0.0
            for p, pe in enumerate(g["planes"]):
                hu = flux_mod.colocated_component(h_arr[t1][..., :nz], g, "Hu", p)
                hv = flux_mod.colocated_component(h_arr[t2][..., :nz], g, "Hv", p)
                # The 0.5 is H's half share of the time average; dS already carries normal_dir
                # and the flux sign convention
                s = s + pe["w"] * 0.5 * cp.sum(
                    (fm["e1"][p] * hv - fm["e2"][p] * hu) * g["dS"])
            if phase == 0:
                fm["buf"][slot] = s
            else:
                fm["buf"][slot] += s
                fm["filled"] = slot + 1


def _flush_args(m, dftab, stepdev, rem, tail):
    """Common arguments for the tail flush kernels.

    ``tail`` carries the size terms specific to each kind: ``(cells, nf)`` for a box monitor,
    ``(n1, n2, nf)`` for a plane monitor.
    """
    return (m["re"], m["im"], m["snap"], dftab["tab32"],
            m["win_e_dev32"], m["win_h_dev32"], stepdev,
            m["T"], np.int32(rem), m["f0"], dftab["n"],
            dftab["tmod"], *tail)


class FlushKernels(NamedTuple):
    """Handles for the six DFT flush kernels.

    Batch boundaries (every T steps) and the tail batch (the one after the run ends) go through the
    same dispatch, so the whole set travels together. ``box_sm2`` is None when the shared-memory
    tiled variant is switched off, in which case the general shared-memory version is used.
    """

    box: object
    box_sm2: object
    box_pf: object
    box_pf_sub: object
    plane: object
    plane_pf: object


def _flush_box(m, kf: FlushKernels, argsf) -> None:
    """Launch dispatch for flushing one box monitor.

    The variants are the P49 shared-memory tile, the general shared-memory version, the P33
    four-component ``_sub`` version, and the per-point one. A batch boundary and the tail batch
    differ only in the remainder inside ``argsf``; the dispatch itself is the same.
    """
    if m["sm"]:
        if kf.box_sm2 is not None:
            kf.box_sm2(m["launch_sm2"][0], m["launch_sm2"][1],
                       argsf, shared_mem=m["smem2"])
        else:
            kf.box(m["launch_flush"][0], m["launch_flush"][1], argsf)
    elif int(m["ncomp"]) < 6:
        # A four-component surface monitor takes the _sub version; its snap layout is
        # (T, ncomp, cells)
        kf.box_pf_sub(m["launch_pf"][0], m["launch_pf"][1],
                      (*argsf, m["cmap_dev"], m["ncomp"]))
    else:
        kf.box_pf(m["launch_pf"][0], m["launch_pf"][1], argsf)


def _flush_plane(m, kf: FlushKernels, argsf) -> None:
    """Launch dispatch for flushing one plane monitor, shared-memory or per-point.

    Same convention as :func:`_flush_box`.
    """
    if m["sm"]:
        kf.plane(m["launch_fl"][0], m["launch_fl"][1], argsf)
    else:
        kf.plane_pf(m["launch_pf"][0], m["launch_pf"][1], argsf)


def _flush_tail(step, fmons, mons, dftab, stepdev, kf: FlushKernels) -> None:
    """Flush the tail batch: the steps still unaccounted for when a pinned step count is not a
    multiple of T.

    Early termination only lands on multiples of interval, and T divides interval, so an early stop
    falls exactly on a batch boundary with remainder 0 and nothing is launched here. By this point
    ``stepdev`` holds the number of steps actually run (incremented at the end of each step), and
    the kernel uses ``*stp-1`` as the end of the batch.
    """
    for m in fmons:
        _T = int(m["T"])
        if not _T or step % _T == 0:
            continue
        _flush_box(m, kf, _flush_args(m, dftab, stepdev, step % _T,
                                      (m["cells"], m["nf"])))
    for m in mons:
        _T = int(m["T"])
        if not _T or step % _T == 0:
            continue
        _flush_plane(m, kf, _flush_args(m, dftab, stepdev, step % _T,
                                        (m["n1"], m["n2"], m["nf"])))


def _phas6(m):
    """Expand (ncomp, nf, ...) to (6, nf, ...), zero-filling the normal components that P33 skipped.

    Those components are never read anyway. Each component is written **directly** into the real or
    imaginary view of the destination array through a pinned staging buffer. The earlier approach
    assembled the whole complex array on the host and copied it in one go, which moved an extra
    57 GB on the largest scene, measured at 70.6 s.
    """
    cm = m.get("cmap", (0, 1, 2, 3, 4, 5))
    re_d, im_d = m["re"], m["im"]
    out = np.zeros((6,) + re_d.shape[1:], dtype=np.complex128)
    buf = None
    for _s, _c in enumerate(cm):
        for _part, _src in ((out[_c].real, re_d[_s]),
                            (out[_c].imag, im_d[_s])):
            if buf is None or buf.shape != _src.shape:
                _mem = cp.cuda.alloc_pinned_memory(_src.nbytes)
                buf = np.frombuffer(_mem, np.float64,
                                    _src.size).reshape(_src.shape)
            _src.get(out=buf)
            _part[...] = buf
    return out


def _readout(mons, fmons, tmons, ftmons, steps_run, e, h, nz,
             return_fields):
    """Move every monitor buffer back to the host once time stepping has finished.

    Returns ``(phasors, field phasors, time samples, time slot counts, time flux, full fields)``.

    ``filled`` for ``tmons`` is recomputed from the formula: the second shot last happens at
    ``m_ref = steps_run - 1``, and only slots inside [beg, end) and aligned to interval are
    counted. This matches the older host version's step-by-step advance exactly.
    """
    for tm in tmons:
        m_last = min(steps_run - 1, tm["end"] - 1)
        off = m_last - tm["beg"]
        tm["filled"] = 0 if off < 0 else min(tm["slots"], off // tm["iv"] + 1)

    out = {}
    for m in mons:
        out[m["name"]] = Phasors(
            name=m["name"],
            data=cp.asnumpy(m["re"]) + 1j * cp.asnumpy(m["im"]),
            freqs=cp.asnumpy(m["freqs"]),
        )
    fphas = {m["name"]: _phas6(m) for m in fmons}
    tsamp = {tm["name"]: cp.asnumpy(tm["buf"][:, : tm["filled"]])
             for tm in tmons}
    tfill = {tm["name"]: tm["filled"] for tm in tmons}
    ftime = {fm["name"]: cp.asnumpy(fm["buf"][: fm["filled"]])
             for fm in ftmons}
    fields = {}
    if return_fields:
        for name, arr in zip(COMPONENTS, (*e, *h)):
            fields[name] = cp.asnumpy(arr[..., :nz])
    return out, fphas, tsamp, tfill, ftime, fields
