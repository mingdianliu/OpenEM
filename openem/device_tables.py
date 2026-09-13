# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Device memory allocation and upload: the host tables the scene side computed are padded to
the pitch and moved onto the device here, and the device-side state buffers the main loop needs
(fields, psi, polarization, step counter) are allocated here too.

It only allocates, pads and uploads, and does no numerics (padding moves addresses, it does not
change values); how things are launched is up to the solver. The tables shared by the real path
(:func:`solver.run`) and the complex Bloch path (:func:`solver._run_bloch`) live in one
function that differs by a single switch argument.
"""

from __future__ import annotations

import cupy as cp
import numpy as np

from openem import cpml
from openem.coeffs import (
    _apply_pec,
    _curl_passthrough,
    _e_coeffs,
    _psi_alloc,
    _psi_flags,
)
from openem.device import grid_1d
from openem.dispersion_setup import _coef_lut
from openem.grid import unravel_cells
from openem.model import Scene
from openem.pitch import (
    Dims,
    _f32,
    _f64,
    _pad32_dev_to,
    _pad32_host_to,
    _pad_z_to,
    _padz_axis,
    _repitch_cells,
)


def _alloc_fields(sc: Scene) -> dict:
    """Fields and auxiliary fields of the real path: the six E/H components and two sets of
    psi, all allocated on the pitch.

    There is one psi per field component per derivative direction, twelve in all. The ones on
    periodic axes are identically 0 yet still take memory and take part in the arithmetic,
    which is known debt.
    The z dimension is pitched to 32 floats (row starts aligned to 128 B, measured at
    1.44-1.48x on the update kernels). The launch only covers the logical nz and the padding
    columns are never written or read; not one value changes, only addresses move (bitwise
    identical).

    Returns ``nx ny nz nzp shape n32``, the six fields, ``e_arr/h_arr`` (triples), ``e3/h3``
    (the per-axis form used by mode source injection), ``psi_h/psi_e/on32`` and the PML
    coefficients ``pe_np/ph_np`` (for :func:`_axis_device_tables`).
    """
    d = Dims.of(sc)
    nx, ny, nz, nzp = d
    shape = d.shape
    n32 = d.n32                                         # trailing three args of the bulk kernel

    ex, ey, ez = (cp.zeros(shape, cp.float32) for _ in range(3))
    hx, hy, hz = (cp.zeros(shape, cp.float32) for _ in range(3))
    # Whether each derivative axis really has a PML: **look at the coefficient arrays
    # directly**, not at sc.pml, so that a case like "there is a PML face but its coefficients
    # happen to be the identity" is not misjudged.
    # On an axis that is off, psi is identically 0 (see stretch in kernels/yee.cu), so only one
    # placeholder element is allocated and the kernel skips the reads and writes. Ring has
    # absorbers on all six faces -> all three axes off, which saves 96 B per cell per step and
    # twelve nx*ny*nz float32 arrays (1.6 GB).
    _pe_np, _ph_np, psi_on, on32 = _psi_flags(sc)

    psi_h = _psi_alloc(shape, psi_on)   # hx_y hx_z hy_z hy_x hz_x hz_y
    psi_e = _psi_alloc(shape, psi_on)   # ex_y ex_z ey_z ey_x ez_x ez_y
    return {
        "dims": d,
        "nx": nx, "ny": ny, "nz": nz, "nzp": nzp, "shape": shape, "n32": n32,
        "ex": ex, "ey": ey, "ez": ez, "hx": hx, "hy": hy, "hz": hz,
        "e_arr": (ex, ey, ez), "h_arr": (hx, hy, hz),
        "e3": ((ex,), (ey,), (ez,)), "h3": ((hx,), (hy,), (hz,)),
        "pe_np": _pe_np, "ph_np": _ph_np, "on32": on32,
        "psi_h": psi_h, "psi_e": psi_e,
    }


def _alloc_fields_c(sc: Scene) -> dict:
    """Fields of the complex Bloch path: two float32 sets (re/im) per component (twice the
    device memory), and two groups of 12 psi.

    Returns ``nx ny nz nzp shape n32`` and several groupings: ``e6/h6`` are the twelve fields
    as they appear in the kernel arguments (exr, exi, eyr, eyi, ezr, ezi / hxr, ...),
    ``e3c/h3c`` are the per-axis ``(re, im)`` pairs (used by mode source injection, the early
    shutoff test and the final readout), and ``er3/ei3/hr3/hi3`` are the real and the imaginary
    triples (the absorber decay multiplies each group once).
    """
    d = Dims.of(sc)           # same as the main path
    nx, ny, nz, nzp = d
    shape = d.shape

    exr, exi, eyr, eyi, ezr, ezi = (cp.zeros(shape, cp.float32) for _ in range(6))
    hxr, hxi, hyr, hyi, hzr, hzi = (cp.zeros(shape, cp.float32) for _ in range(6))
    psi_h = tuple(cp.zeros(shape, cp.float32) for _ in range(12))  # (hx_y hx_z hy_z hy_x hz_x hz_y) x (re im)
    psi_e = tuple(cp.zeros(shape, cp.float32) for _ in range(12))
    return {
        "dims": d,
        "nx": nx, "ny": ny, "nz": nz, "nzp": nzp, "shape": shape,
        "n32": d.n32,
        "e6": (exr, exi, eyr, eyi, ezr, ezi), "h6": (hxr, hxi, hyr, hyi, hzr, hzi),
        "e3c": ((exr, exi), (eyr, eyi), (ezr, ezi)),   # per-axis (re, im) for mode injection
        "h3c": ((hxr, hxi), (hyr, hyi), (hzr, hzi)),
        "er3": (exr, eyr, ezr), "ei3": (exi, eyi, ezi),
        "hr3": (hxr, hyr, hzr), "hi3": (hxi, hyi, hzi),
        "psi_h": psi_h, "psi_e": psi_e,
    }


def _axis_device_tables(sc: Scene, tabs, pml_np, d: Dims,
                        complex_masks: bool = False) -> dict:
    """Upload the 1D tables of the three axes to device memory (shared by the real path and the
    Bloch path).

    The 1D z tables are padded to nzp as well: rounding the launch up to whole blocks spills
    into padding threads in [nz, nzp) (the guard uses the pitch), and those threads read these
    tables. Safe padding values: index tables pad with nz-1 (inside the domain), masks with 0,
    spacings with 1, PML with the identity, pec with 1. A padding cell has ca=cb=0 so its E is
    identically 0; the junk only lands in the padding and is never read by a physical cell.

    ``pml_np`` is ``(pe_np, ph_np)``, the per-axis ``(a, b, invk)`` of the E side and the H
    side (the output of _pml_axis_arrays), and ``tabs`` is the per-axis ``Axis.index_tables``.
    With ``complex_masks=True`` (Bloch), ``mnx`` / ``mpv`` each split into ``(re, im)``; on an
    axis with k=0 the imaginary part is all 0 (np.imag returns zeros for a real array).
    """
    pe_np, ph_np = pml_np
    nz, nzp = d.nz, d.nzp
    _pml_fill = (0.0, 1.0, 1.0)
    out = {
        "idl_p": [_f32(_padz_axis(v, a, 1.0, nz, nzp))
                  for a, v in enumerate([1.0 / x.dl for x in sc.grid.axes])],
        "idl_d": [_f32(_padz_axis(v, a, 1.0, nz, nzp))
                  for a, v in enumerate([1.0 / x.dl_dual for x in sc.grid.axes])],
        "pml_e": [[_f32(_padz_axis(v, a, _pml_fill[i], nz, nzp))
                   for i, v in enumerate(x)] for a, x in enumerate(pe_np)],
        "pml_h": [[_f32(_padz_axis(v, a, _pml_fill[i], nz, nzp))
                   for i, v in enumerate(x)] for a, x in enumerate(ph_np)],
        "nxt": [cp.asarray(_padz_axis(t["nxt"], a, nz - 1, nz, nzp))
                for a, t in enumerate(tabs)],
        "prv": [cp.asarray(_padz_axis(t["prv"], a, nz - 1, nz, nzp))
                for a, t in enumerate(tabs)],
    }
    if complex_masks:
        out["mnx"] = [(_f32(_padz_axis(np.real(t["mnx"]), a, 0.0, nz, nzp)),
                       _f32(_padz_axis(np.imag(t["mnx"]), a, 0.0, nz, nzp)))
                      for a, t in enumerate(tabs)]
        out["mpv"] = [(_f32(_padz_axis(np.real(t["mpv"]), a, 0.0, nz, nzp)),
                       _f32(_padz_axis(np.imag(t["mpv"]), a, 0.0, nz, nzp)))
                      for a, t in enumerate(tabs)]
    else:
        out["mnx"] = [_f32(_padz_axis(t["mnx"], a, 0.0, nz, nzp))
                      for a, t in enumerate(tabs)]
        out["mpv"] = [_f32(_padz_axis(t["mpv"], a, 0.0, nz, nzp))
                      for a, t in enumerate(tabs)]
    out["pec"] = [_f32(_padz_axis(t["pec"], a, 1.0, nz, nzp))
                  for a, t in enumerate(tabs)]
    return out


def _e_coeff_tables(sc: Scene, tabs, d: Dims) -> dict:
    """E update coefficients ``ca`` / ``cb`` (one array per component, already padded to the
    pitch) and the tensor cell table ``ten``; the passthrough cells (mix / tensor) and the PEC
    cells rewrite ca/cb here. Every coefficient is prepared here so that the kernel only does
    multiply-adds. Returns ``{"ch", "ca", "cb", "ten"}``."""
    nx, ny, nz, nzp = d
    ch = np.float32(sc.dt / cpml.MU_0)

    g_dense = (sc.dispersion.g_dense(sc.shape) if sc.any_dispersion
               else (np.zeros((nx, ny, nz)),) * 3)
    ca, cb = [], []
    for eps, sig, gg in ((sc.eps_ex, sc.sigma_ex, g_dense[0]),
                         (sc.eps_ey, sc.sigma_ey, g_dense[1]),
                         (sc.eps_ez, sc.sigma_ez, g_dense[2])):
        a_, b_ = _e_coeffs(eps, sig, gg, sc, d)
        ca.append(a_)
        cb.append(b_)
    # Cells that take the D form: ``ca=0, cb=1`` makes update_e write the **stretched curl**
    # into the E array unchanged, and update_e_from_d turns it into the real E afterwards. No
    # branch is needed inside update_e.
    # A side benefit: source injection goes through cb as well, so ``-cb·coef·amp`` becomes
    # exactly ``-coef·amp`` here and the D update multiplies by dt/ε₀ later, which makes
    # current sources come out right on both paths by themselves.
    if sc.any_dispersion_mix:
        sc.dispersion_mix.validate(sc.shape, sc.dispersion)
        _curl_passthrough(sc.dispersion_mix, ca, cb, d)

    # Tensor cells ride the same passthrough track: update_e writes "stretched curl - J" into
    # the E array unchanged and tensor_dot/tensor_scatter turns it into the real E afterwards
    # (the derivation is at the top of kernels/tensor.cu). Coexistence with dispersion or mix
    # already fails closed in build, so the three kinds of passthrough cell cannot overwrite
    # each other.
    ten = None
    if sc.any_tensor:
        tt = sc.tensor
        tt.validate(sc.shape, tabs)
        _curl_passthrough(tt, ca, cb, d)
        m1_t, m2_t = tt.rows(sc.dt)
        ten = {
            "n": np.int32(tt.n_entry),
            "comp": cp.asarray(np.ascontiguousarray(tt.comp, dtype=np.int32)),
            "cell": cp.asarray(_repitch_cells(tt.cell, ny, nz, nzp)),
            "m1": _f32(m1_t), "m2": _f32(m2_t),
            "rhs": cp.zeros(tt.n_entry, cp.float32),
            "enew": cp.zeros(tt.n_entry, cp.float32),
            "launch": grid_1d(tt.n_entry),
        }

    # PEC cells: E is identically 0 (staircasing). ``ca=cb=0`` makes update_e write E back to
    # 0 on every step;
    _apply_pec(sc, ca, cb, d)
    return {"ch": ch, "ca": ca, "cb": cb, "ten": ten}


def _bloch_e_coeffs(sc: Scene, d: Dims) -> dict:
    """``ca ≡ 1`` and ``cb = dt/(ε₀ε)`` for the Bloch subset (lossless and dispersionless:
    s = G = 0). Returns ``{"ch", "ca", "cb"}``."""
    nz, nzp = d.nz, d.nzp
    ch = np.float32(sc.dt / cpml.MU_0)
    ca, cb = [], []
    for eps in (sc.eps_ex, sc.eps_ey, sc.eps_ez):
        den = eps  # the subset is lossless and dispersionless: s = G = 0
        ca.append(_f32(_pad_z_to(np.ones_like(eps), nz, nzp)))
        cb.append(_f32(_pad_z_to((sc.dt / cpml.EPSILON_0) / den, nz, nzp)))
    return {"ch": ch, "ca": ca, "cb": cb}


def _dispersion_tables(sc: Scene, d: Dims) -> dict:
    """Dispersion: the polarization auxiliary fields plus the history term.

    Without dispersion the three hist entries point at **one and the same** all-zero array: the
    kernel does its multiply-adds as usual, the result is unchanged, and no branch is needed in
    the kernel. The cost is one nx*ny*nz float32 array (known debt).

    Returns ``{"hist", "disp", "hist_box", "cw"}``; ``disp`` is None without dispersion.
    """
    shape = d.shape
    ny, nz, nzp = d.ny, d.nz, d.nzp
    zeros_hist = cp.zeros(shape, cp.float32)
    hist = (zeros_hist, zeros_hist, zeros_hist)
    disp = None
    if sc.any_dispersion:
        d = sc.dispersion
        d.validate(sc.shape)
        hist = tuple(cp.zeros(shape, cp.float32) for _ in range(3))
        # Lossless deduplication of the (am1, b) pairs: the table holds the original values
        # (unique does not touch a bit), so a lookup is bitwise identical. The index is int32;
        # the table peaks at 8.4 MB (Topo) and lives in L2.
        uniq, inv = _coef_lut(d.am1, d.b)
        disp = {
            "n": np.int32(d.n_entry),
            "p_re": cp.zeros(d.n_pole, cp.float32),
            "p_im": cp.zeros(d.n_pole, cp.float32),
            "lut_am1_re": _f32(uniq["am1"].real),
            "lut_am1_im": _f32(uniq["am1"].imag),
            "lut_b_re": _f32(uniq["b"].real),
            "lut_b_im": _f32(uniq["b"].imag),
            # inv may already be on the device: astype works on both sides and asarray is a
            # no-op for a device array
            "cidx": cp.asarray(inv.astype(np.int32)),
            "comp": cp.asarray(np.ascontiguousarray(d.comp, dtype=np.int32)),
            "cell": cp.asarray(_repitch_cells(d.cell, ny, nz, nzp)),
            "ofs": cp.asarray(np.ascontiguousarray(d.pole_ofs, dtype=np.int32)),
            "launch": grid_1d(d.n_entry),
        }
    # Bounding box of the dispersion hist (dispersion only writes the cells it has entries
    # for). An empty box means no dispersion, and update_e skips the hist read everywhere.
    # cell has already been re-encoded on nzp, so it is decoded with the same stride.
    hist_box = (np.int32(0),) * 6
    if disp is not None:
        _cell = disp["cell"]
        _ii, _jj, _kk = unravel_cells(_cell, (ny, nzp))
        hist_box = tuple(np.int32(int(v)) for v in
                         (_ii.min(), _ii.max() + 1, _jj.min(), _jj.max() + 1,
                          _kk.min(), _kk.max() + 1))
        del _cell, _ii, _jj, _kk
    cw = np.float32(cpml.EPSILON_0 / sc.dt)
    return {"hist": hist, "disp": disp, "hist_box": hist_box, "cw": cw}


def _mix_table(sc: Scene, d: Dims) -> dict | None:
    """Entry table of the dispersion mix plus the state buffers of its two polarizations;
    returns None when there is no mix."""
    if not sc.any_dispersion_mix:
        return None
    ny, nz, nzp = d.ny, d.nz, d.nzp
    mm = sc.dispersion_mix
    k1, k2, omb, invd = mm.coefficients()
    return {
        "n": np.int32(mm.n_entry),
        "d_tot": cp.zeros(mm.n_entry, cp.float32),
        "d_h": cp.zeros(mm.n_entry, cp.float32),
        "sp": cp.zeros(mm.n_entry, cp.float32),
        "sq": cp.zeros(mm.n_entry, cp.float32),
        "p_re": cp.zeros(mm.pa.size, cp.float32),
        "p_im": cp.zeros(mm.pa.size, cp.float32),
        "q_re": cp.zeros(mm.qa.size, cp.float32),
        "q_im": cp.zeros(mm.qa.size, cp.float32),
        "pa_re": _f32(mm.pa.real), "pa_im": _f32(mm.pa.imag),
        "pb_re": _f32(mm.pb.real), "pb_im": _f32(mm.pb.imag),
        "qa_re": _f32(mm.qa.real), "qa_im": _f32(mm.qa.imag),
        "qb_re": _f32(mm.qb.real), "qb_im": _f32(mm.qb.imag),
        "p_ofs": cp.asarray(np.ascontiguousarray(mm.p_ofs, dtype=np.int32)),
        "q_ofs": cp.asarray(np.ascontiguousarray(mm.q_ofs, dtype=np.int32)),
        "comp": cp.asarray(np.ascontiguousarray(mm.comp, dtype=np.int32)),
        "cell": cp.asarray(_repitch_cells(mm.cell, ny, nz, nzp)),
        "k1": _f32(k1), "k2": _f32(k2), "omb": _f32(omb), "invd": _f32(invd),
        "dte": np.float32(sc.dt / cpml.EPSILON_0),
        "launch": grid_1d(mm.n_entry),
    }


def _dual_volumes_eps(sc: Scene, d: Dims, gpu: bool) -> dict:
    """Per-component Yee dual volumes and the three ε arrays, used for the volume weighting of
    the shutoff (grid.yee_dual_volumes is where they come from). Returns
    ``{"vols": (volx, voly, volz), "eps_dev"}``.

    With ``gpu=True`` (OPENEM_INIT_GPU, the default on the real path) the outer product of the
    volumes and the padding and narrowing of eps move onto the GPU. Elementwise multiplication
    and f64->f32 rounding are both exactly specified by IEEE754 => bitwise identical to the CPU
    version. The Bloch path always takes the CPU version.
    """
    nz, nzp = d.nz, d.nzp
    if gpu:
        _xp, _yp, _zp = (cp.asarray(a.dl) for a in sc.grid.axes)
        _xd, _yd, _zd = (cp.asarray(a.dl_dual) for a in sc.grid.axes)
        volx = _pad32_dev_to(_xp[:, None, None] * _yd[None, :, None] * _zd[None, None, :], nz, nzp)
        voly = _pad32_dev_to(_xd[:, None, None] * _yp[None, :, None] * _zd[None, None, :], nz, nzp)
        volz = _pad32_dev_to(_xd[:, None, None] * _yd[None, :, None] * _zp[None, None, :], nz, nzp)
        del _xp, _yp, _zp, _xd, _yd, _zd
        eps_dev = tuple(_pad32_host_to(e, nz, nzp)
                        for e in (sc.eps_ex, sc.eps_ey, sc.eps_ez))
    else:
        volx, voly, volz = (_f32(_pad_z_to(v, nz, nzp)) for v in sc.grid.yee_dual_volumes())
        eps_dev = tuple(_f32(_pad_z_to(e, nz, nzp)) for e in (sc.eps_ex, sc.eps_ey, sc.eps_ez))
    return {"vols": (volx, voly, volz), "eps_dev": eps_dev}


def _plane_wave_tables(plane_src, ex, ey, hx, hy, cb) -> list:
    """Plane wave sources: the main source and the multi-polarization "extra" group
    (sources_setup._plane_source) go into one list with the main source first, so the step body
    has a single for loop and the launch order is the same as when it was written as two
    separate blocks. Each group uploads its own two incident tables to device memory; the
    target arrays of the injection follow the polarization: Ex polarization -> (Ex, Hy), Ey
    polarization -> (Ey, Hx), and the sign difference is already folded into coef_h (see
    _plane_source).

    Each item is ``(src, ex_tab, hy_tab, H target, E target, cb)``; the list is empty when
    there is no plane wave.
    """
    pw_all = []
    if plane_src is not None:
        for _px in (plane_src, *plane_src.get("extra", [])):
            pw_all.append((_px, cp.asarray(_px["ex_inc"]), cp.asarray(_px["hy_inc"]),
                           hy if _px["pol"] == 0 else hx,
                           ex if _px["pol"] == 0 else ey, cb[_px["pol"]]))
    return pw_all


def _step_clock(n_total: int, dt: float) -> dict:
    """Device step counter and the per-step instant tables.

    ``stepdev`` starts at -1, so the _sample_time before the loop corresponds to
    m_step = *step+1 = 0 (the same meaning as _sample_time(0) in the old host version), and it
    is zeroed before entering the loop.
    Instant tables: the same formula as t_e_all/t_h_all in _monitor_buffers (bitwise
    identical).
    Returns ``{"stepdev", "t_e_dev", "t_h_dev"}``.
    """
    stepdev = cp.full(1, -1, dtype=cp.int32)
    t_e_dev = _f64((np.arange(n_total, dtype=np.float64) + 1.0) * dt)
    t_h_dev = _f64((np.arange(n_total, dtype=np.float64) + 0.5) * dt)
    return {"stepdev": stepdev, "t_e_dev": t_e_dev, "t_h_dev": t_h_dev}
