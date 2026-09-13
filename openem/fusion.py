# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Thresholds and launch parameters of the individual fusions and fast paths.

There is plenty to save per step: absorber decay merged into update_h, the dispersion pole
update merged into update_e, two-pole entries on the float2 fast path, the H step and the E step
merged into one kernel. But each one has its own preconditions, and missing a single one breaks
bitwise equivalence. What lives here is the on-or-off decision, plus the static parameters that
have to be prepared once it is on.

The decisions are pure functions and can be tested on their own; the parameter helpers read the
device arrays in place and return static tuples that a Graph can capture.
"""

from __future__ import annotations

from typing import NamedTuple

import cupy as cp
import numpy as np

from openem import knobs

from openem.dispersion_setup import _pole_bucket_setup
from openem.pitch import Dims


def _no_h_source(plane_src, tfsf, mode_src, dip_h) -> bool:
    """No source injection at all on the H step.

    This is the shared precondition of the optimizations that fuse a source injection and some
    other H-side operation into one kernel: only when nothing sits between the source and the
    step being fused in can the fusion avoid missing the source's share.
    """
    return (plane_src is None and tfsf is None and mode_src is None
            and dip_h is None)


def _plain_step(absb, abs_slabs, disp, mix, mod, ten) -> bool:
    """The step does nothing besides update_h / update_e.

    That is, no absorber, dispersion, mixed medium, time-varying medium or tensor medium: only
    when the gap between the H step and the E step is empty can the two be fused into one kernel.
    """
    return (absb is None and not abs_slabs and disp is None
            and mix is None and mod is None and ten is None)


def _graph_allowed(ten, mix, ftmons, n_total) -> bool:
    """Whether the step body can be captured in a CUDA Graph.

    Tensor and mixed media need a host-side decision every step, and time-domain flux monitors
    need a readback every step; neither fits into a static graph. Below three steps the capture
    does not pay for itself anyway.
    """
    return (knobs.env("NO_GRAPH") != "1"
            and ten is None and mix is None and not ftmons
            and n_total > 2)


def _ade_fuse_allowed(disp, has_lut, pabs, ten, mod, mix) -> bool:
    """Threshold for P43 full ADE fusion (pole update inlined into update_e, hist gone entirely).

    update_e_ade has no hist argument => it is only usable when every entry can be deferred,
    that is a dispersive scene with no absorber (in-layer entries have to stay in
    dispersion_step_dec).

    **Off** by default: measured at 0.82-0.90x. The fusion forces the pole update from per-entry
    parallelism over to per-cell parallelism, and warp divergence plus the strided access into
    the p arrays costs 27%, while the hist round trip and index traffic it saves come to only
    18%; the ceiling is below the cost.
    """
    return (disp is not None and has_lut and pabs is None
            and knobs.env("ADE_FUSE") == "1"
            and ten is None and mod is None and mix is None)


def _lean_allowed(has_lut, ade, mod, ten, disp) -> bool:
    """Outer threshold of P25b lean mode (interior / boundary shell split).

    If P43 has already folded the dispersion into update_e (``ade`` is not None), lean is not
    used.

    Tensor scenes are let through only when there is **no dispersion**: the two tensor passes
    (ktdot / ktscat) only touch the cells their entries sit in and do not change E anywhere
    else, so their one conflict with lean is the deferred dispersion (the second pass would
    change E), which does not exist without dispersion. Dispersion plus tensor is still
    excluded.
    """
    return (knobs.env("LEAN") != "0" and has_lut
            and ade is None and mod is None
            and (ten is None or disp is None))


def _lean_report(lean, verbose) -> None:
    """Print the lean verdict (when verbose)."""
    if lean is not None and verbose:
        print(f"  [P25b] interior/boundary-shell split: interior is {lean['frac']:.0%}, "
              f"dispersion entries interior {lean['n_int']:,} / shell {lean['n_skin']:,}")


def _h_fuse_allowed(absb, plane_src, tfsf, mode_src, dip_h, verbose) -> bool:
    """Whether update_h and the H-family absorber decay can be fused into one kernel (and print
    one line about it).

    The gate: there is an absorber, and no source injection on the H step; only with nothing
    between the two can the fusion avoid missing the source's share of the decay. Measured at
    1.35-1.53x, bitwise identical.
    """
    ok = (absb is not None
          and knobs.env("H_FUSE") != "0"
          and _no_h_source(plane_src, tfsf, mode_src, dip_h))
    if ok and verbose:
        print("  [P37] update_h fused with the H-family absorber decay (no H-side source)")
    return ok


def _abs_e_fold(disp_fused_dec, absb, disp, d: Dims, verbose) -> bool:
    """With full coverage by the dispersion entries, the E-side absorber pass folds into
    dispersion_step_dec.

    "Full coverage" means exactly one entry per (comp, cell), 3 * cell count in total. The en2
    that kdstepd computes and writes back into E is then equivalent to absorb's whole E decay
    pass, and that pass can be skipped on every step but the last.
    """
    nx, ny, nz, nzp = d
    if not (disp_fused_dec and absb is not None
            and knobs.env("ABS_E_FOLD") != "0"
            # entries count with the real nz, not the pitch
            and int(disp["n"]) == 3 * nx * ny * nz):
        return False
    _c = cp.asnumpy(disp["comp"]).astype(np.int64)
    _id = cp.asnumpy(disp["cell"]).astype(np.int64)
    _seen = np.zeros(3 * nx * ny * nzp, dtype=bool)
    _seen[_c * (nx * ny * nzp) + _id] = True
    ok = bool(int(_seen.sum()) == int(disp["n"]))
    if ok and verbose:
        print("  [P44] full dispersion coverage: E-side absorber decay folded into "
              "dispersion_step_dec")
    return ok


def _absorber_1d_profiles(abs_slabs, d: Dims):
    """  Six 1D decay profiles [xi, xh, yi, yh, zi, zh]; returns None when same-axis slabs
    overlap (or the knob is off)."""
    if knobs.env("ABS_1D") == "0":
        return None
    dims = d.shape
    for a in range(3):
        ivs = sorted((int(s["org"][a]), int(s["org"][a]) + int(s["box"][a]))
                     for s in abs_slabs if int(s["axis"]) == a)
        if any(y[0] < x[1] for x, y in zip(ivs, ivs[1:])):
            return None
    g = [np.ones(dims[a], np.float32) for a in range(3) for _ in (0, 1)]
    for s in abs_slabs:
        a = int(s["axis"])
        b0 = int(s["org"][a])
        n = int(s["box"][a])
        g[2 * a][b0:b0 + n] = cp.asnumpy(s["decay_int"])[:n]
        g[2 * a + 1][b0:b0 + n] = cp.asnumpy(s["decay_half"])[:n]
    return tuple(cp.asarray(x) for x in g)


def _absorber_batch(abs_slabs, d: Dims, verbose):
    """Merge the decay of the individual absorber layers into a single launch and return the
    launch parameters; returns None when merging is not worth it.

    Launching per slab (10 times per step on Ring, 8 on hBN) measures only 27% of roofline; the
    merged single launch is 1.45-1.56x (the multiplication order is unchanged => bitwise
    identical).

    When no two same-axis slabs intersect, six 1D profiles are premultiplied on top of that
    (3 axes x integer/half cell). Each cell gets at most one factor per axis => the value is the
    original single factor and the order follows the slab order, so that is bitwise identical
    too. If non-intersection does not hold, or ``OPENEM_ABS_1D=0``, ``g6`` is left None and the
    original path is taken.
    """
    if len(abs_slabs) < 2 or knobs.env("ABS_BATCH") == "0":
        return None
    nx, ny, nzp = d.shape

    def _ptrs(arrs):
        return cp.asarray(np.array([int(a.data.ptr) for a in arrs],
                                   dtype=np.uint64))

    def _i32(v):
        return cp.asarray(np.array(v, dtype=np.int32))

    _o = [min(int(s["org"][a]) for s in abs_slabs) for a in range(3)]
    _e = [max(int(s["org"][a]) + int(s["box"][a]) for s in abs_slabs)
          for a in range(3)]
    _nt = 128
    absb = {
        "di": _ptrs([s["decay_int"] for s in abs_slabs]),
        "dh": _ptrs([s["decay_half"] for s in abs_slabs]),
        "ax": _i32([int(s["axis"]) for s in abs_slabs]),
        "i0": _i32([int(s["org"][0]) for s in abs_slabs]),
        "j0": _i32([int(s["org"][1]) for s in abs_slabs]),
        "k0": _i32([int(s["org"][2]) for s in abs_slabs]),
        "bi": _i32([int(s["box"][0]) for s in abs_slabs]),
        "bj": _i32([int(s["box"][1]) for s in abs_slabs]),
        "bk": _i32([int(s["box"][2]) for s in abs_slabs]),
        "ns": np.int32(len(abs_slabs)),
        "o": tuple(np.int32(v) for v in _o),
        "e": tuple(np.int32(v) for v in _e),
        "grid": (((_e[2] - _o[2]) + _nt - 1) // _nt,
                 ((_e[1] - _o[1]) + 1) // 2, _e[0] - _o[0]),
        "block": (_nt, 2, 1),
    }
    absb["g6"] = _absorber_1d_profiles(abs_slabs, d)
    if absb["g6"] is not None and verbose:
        print("  [P50] absorb 1D premultiplied profiles: the slab loop disappears (bitwise)")
    return absb


def _pole_bucket(disp, ade, lean, ten, pabs, verbose):
    """float2 fast path for the 2-pole bucket (bitwise identical; 1.35x on the benchmark).

    ade / lean permute the entries again, which is mutually exclusive with bucketing, so this is
    not enabled while either of those paths is active.
    """
    if not (disp is not None and ade is None and lean is None and ten is None
            and knobs.env("DISP_P2") != "0"):
        return None
    p2 = _pole_bucket_setup(disp, pabs)
    if p2 is not None and verbose:
        print(f"  [P48] 2-pole bucket: {p2['n2']:,}/{int(disp['n']):,} entries on the float2 "
              f"fast path ({p2['n_pre']:,} before / {p2['n_post']:,} after on the general path)")
    return p2


class DispPlan(NamedTuple):
    """The device pointers and thresholds a dispersion kernel launch needs:
    :func:`_dispersion_args` and :func:`_pole_kernels` read the same bundle, which spelled out
    flat would be nine parameters listed on each side.

    ``e`` is the three E components ``(ex, ey, ez)``; ``dec`` / ``has_dec`` are the in-layer P
    decay factors (a placeholder array plus 0 when there is no in-layer dispersion);
    ``fused`` / ``fused_dec`` are the thresholds of the two paths that do post+pre in one pass.
    """

    disp: dict | None
    hist: tuple
    e: tuple
    dec: object
    has_dec: object
    pabs: dict | None
    fused: bool
    fused_dec: bool


def _pole_kernels(p2, k, dp: DispPlan):
    """Static tuples for the segmented launch of the 2-pole bucket (capturable by a Graph).

    The general kernel runs a subrange off **sliced pointers**: per-entry arrays are sliced
    ``[a:b]`` and ``ofs`` is sliced ``[a:b+1]`` (its values are still absolute slot numbers, so
    one extra boundary is needed at the end), and ``n`` is passed the length of the range.

    The order of the returned nine-tuple matches the unpacking at the call site; when ``p2`` is
    None all nine are None.
    """
    disp, hist, pabs = dp.disp, dp.hist, dp.pabs
    ex, ey, ez = dp.e
    _dec = dp.dec
    kdstep_p2 = kdstepd_p2 = kdstepd_we_p2 = None
    args_dstep_p2 = args_dstepd_p2 = None
    args_dstep_pre = args_dstep_post = None
    args_dstepd_pre = args_dstepd_post = None
    if p2 is not None:
        _lo, _hi = p2["lo"], p2["lo"] + p2["n2"]
        _lut4 = (disp["lut_am1_re"], disp["lut_am1_im"],
                 disp["lut_b_re"], disp["lut_b_im"])
        if dp.fused:
            kdstep_p2 = k["dispersion_step_p2"]
            args_dstep_p2 = (disp["p_re"], disp["p_im"], *hist, *_lut4,
                             disp["cidx"], disp["cellc"][_lo:_hi],
                             ex, ey, ez, p2["slot0"], np.int32(p2["n2"]))
            args_dstep_pre = (disp["p_re"], disp["p_im"], *hist, *_lut4,
                              disp["cidx"], disp["comp"][:_lo],
                              disp["cell"][:_lo], disp["ofs"][:_lo + 1],
                              ex, ey, ez, np.int32(p2["n_pre"]))
            args_dstep_post = (disp["p_re"], disp["p_im"], *hist, *_lut4,
                               disp["cidx"], disp["comp"][_hi:],
                               disp["cell"][_hi:], disp["ofs"][_hi:],
                               ex, ey, ez, np.int32(p2["n_post"]))
        if dp.fused_dec:
            kdstepd_p2 = k["dispersion_step_dec_p2"]
            kdstepd_we_p2 = k["dispersion_step_dec_we_p2"]
            args_dstepd_p2 = (disp["p_re"], disp["p_im"], *hist, *_lut4,
                              disp["cidx"], disp["cellc"][_lo:_hi],
                              ex, ey, ez, _dec[_lo:_hi],
                              pabs["f1"][_lo:_hi], pabs["f2"][_lo:_hi],
                              pabs["f3"][_lo:_hi], p2["slot0"],
                              np.int32(p2["n2"]))
            args_dstepd_pre = (disp["p_re"], disp["p_im"], *hist, *_lut4,
                               disp["cidx"], disp["comp"][:_lo],
                               disp["cell"][:_lo], disp["ofs"][:_lo + 1],
                               ex, ey, ez, _dec[:_lo], pabs["f1"][:_lo],
                               pabs["f2"][:_lo], pabs["f3"][:_lo],
                               np.int32(p2["n_pre"]))
            args_dstepd_post = (disp["p_re"], disp["p_im"], *hist, *_lut4,
                                disp["cidx"], disp["comp"][_hi:],
                                disp["cell"][_hi:], disp["ofs"][_hi:],
                                ex, ey, ez, _dec[_hi:], pabs["f1"][_hi:],
                                pabs["f2"][_hi:], pabs["f3"][_hi:],
                                np.int32(p2["n_post"]))
    return (kdstep_p2, kdstepd_p2, kdstepd_we_p2,
            args_dstep_p2, args_dstepd_p2,
            args_dstep_pre, args_dstep_post,
            args_dstepd_pre, args_dstepd_post)


def _dispersion_args(dp: DispPlan):
    """Actual arguments of the four dispersion kernel groups: pre / step / step_dec / post.

    ``pre`` and ``step`` take the same argument list (the same pointers in the same order) and
    differ only in the threshold: ``pre`` is needed every time, ``step`` is only launched on the
    fused post+pre path. So they simply share one tuple object here.

    Returns ``(args_dpre, args_dstep, args_dstepd, args_dpost)``; the group whose threshold is
    not met is None and its kernel is not launched. When ``disp`` is None all four are None.
    """
    disp, pabs, _dec = dp.disp, dp.pabs, dp.dec
    if disp is None:
        return None, None, None, None
    common = (disp["p_re"], disp["p_im"], *dp.hist,
              disp["lut_am1_re"], disp["lut_am1_im"],
              disp["lut_b_re"], disp["lut_b_im"], disp["cidx"],
              disp["comp"], disp["cell"], disp["ofs"], *dp.e)
    args_dpre = (*common, disp["n"])
    args_dstepd = ((*common, _dec, pabs["f1"], pabs["f2"], pabs["f3"],
                    disp["n"])
                   if dp.fused_dec else None)
    args_dpost = (disp["p_re"], disp["p_im"], disp["lut_b_re"],
                  disp["lut_b_im"], disp["cidx"], disp["comp"], disp["cell"],
                  disp["ofs"], *dp.e, _dec, dp.has_dec, disp["n"])
    return (args_dpre, args_dpre if dp.fused else None,
            args_dstepd, args_dpost)
