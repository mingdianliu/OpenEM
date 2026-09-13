# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""E/H update coefficients, PML axis tables, PEC masks, the cacb lookup table, psi knobs.

All of this is used by run's setup phase but has nothing to do with its local scope: it was
lifted out of run in place, and only later moved into this file.

"""

from __future__ import annotations

import cupy as cp
import numpy as np

from openem import cpml
from openem import knobs
from openem.model import Scene
from openem.pitch import Dims, _f32, _pad_z_to, _repitch_cells


def _pml_axis_arrays(sc: Scene, axis: int, at: str) -> tuple[np.ndarray, ...]:
    """Spread the PML coefficients into 1D arrays of length n, **filling the identity value
    outside the PML**.

    This is what keeps the kernel branch-free: ``a=0, b=1, invk=1`` makes psi identically 0 and
    lets the derivative through unchanged, so interior cells and PML cells take the same code
    path.
    """
    n = sc.grid.axes[axis].n
    a = np.zeros(n)
    b = np.ones(n)
    invk = np.ones(n)
    for side in ("lo", "hi"):
        c = sc.pml.get((axis, side))
        if c is None:
            continue
        m = c.num_layers
        sl = slice(0, m) if side == "lo" else slice(n - m, n)
        a[sl] = getattr(c, f"a_{at}")
        b[sl] = getattr(c, f"b_{at}")
        invk[sl] = getattr(c, f"inv_kappa_{at}")
    return a, b, invk


#: derivative axis of each of psi's six slots (hx_y hx_z hy_z hy_x hz_x hz_y, same order on E)
_PSI_AX = (1, 2, 2, 0, 0, 1)


def _idfrac(arrs) -> float:
    """Fraction of identity elements (a=0, b=1) in the PML coefficient pair ``(a, b)``, used as
    the threshold of the P38 per-axis fast path."""
    m = (np.asarray(arrs[0]) == 0.0) & (np.asarray(arrs[1]) == 1.0)
    return float(m.mean())


def _psi_alloc(shape, psi_on) -> tuple:
    """The six psi slots; an axis with no PML gets a single element."""
    return tuple(cp.zeros(shape if psi_on[ax] else 1, cp.float32)
                 for ax in _PSI_AX)


def _psi_flags(sc):
    """The psi knobs of the three axes and the per-axis fast-path level of P38, returning
    ``(pe_np, ph_np, psi_on, on32)``.

    If the PML coefficients on an axis are constantly (a=0, b=1), that axis has no PML and the
    kernel can skip reading and writing psi altogether: on Ring, with absorbers on all six faces,
    all three axes are off, which saves 96 B per cell per step and twelve nx*ny*nz float32 arrays
    (1.6 GB).

    Even with a PML, psi is identically 0 over the identity region (a==0 and b==1), so reading
    and writing it there is wasted as well; but when the hit rate is low, that branch gets in the
    way of issuing the psi loads early and measures 6% slower instead. So ``on32`` is decided per
    axis from the hit rate:

      0 = no PML on this axis (psi gets a single placeholder element)
      1 = PML present, read and write psi as usual
      2 = PML present and the identity region is large enough, take the fast path
    """
    pe_np = [list(_pml_axis_arrays(sc, a, "E")) for a in range(3)]
    ph_np = [list(_pml_axis_arrays(sc, a, "H")) for a in range(3)]
    psi_on = tuple(
        bool(np.any(pe_np[a][0] != 0.0) or np.any(pe_np[a][1] != 1.0)
             or np.any(ph_np[a][0] != 0.0) or np.any(ph_np[a][1] != 1.0))
        for a in range(3)
    )
    skip_min = float(knobs.env("PSI_SKIP_FRAC"))
    on = []
    for a in range(3):
        if not psi_on[a]:
            on.append(np.int32(0))
            continue
        frac = min(_idfrac(pe_np[a]), _idfrac(ph_np[a]))
        on.append(np.int32(2 if frac >= skip_min else 1))
    return pe_np, ph_np, psi_on, tuple(on)


def _e_coeffs(eps: np.ndarray, sigma: np.ndarray | None, gg: np.ndarray,
              sc, d: Dims):
    """``(ca, cb)``, see the formula at the top of kernels/yee.cu.

        ca = (eps_inf - s)/(eps_inf + s + G),  cb = (dt/eps0)/(eps_inf + s + G)

    ``s = σ dt/(2 ε₀)``, and ``G`` is the implicit coupling of the dispersion. With all three at
    zero this degenerates into the explicit update.
    """
    nx, ny, nz, nzp = d
    if knobs.env("INIT_GPU") != "0":
        # Whole-domain elementwise math moved onto the GPU. Same expression, same IEEE754
        # double-precision order => bitwise identical to the CPU version (+,-,*,/ and the
        # f64->f32 rounding are all exactly specified).
        # The **blocking** cuts along axis 0: whole-domain float64 temporaries cost nothing on
        # the CPU, but stuffing them into 80 GB of device memory runs out of memory (1.94 GB
        # each on the PEC case, ~58 GB stacked up). Every block does the same elementwise math
        # independently, bitwise unchanged, and the peak is down to one block.
        _ca = cp.zeros((nx, ny, nzp), cp.float32)
        _cb = cp.zeros((nx, ny, nzp), cp.float32)
        _bytes_plane = ny * nz * 8
        _step = max(1, min(nx, int(256e6 // max(_bytes_plane, 1))))
        for _i0 in range(0, nx, _step):
            _i1 = min(nx, _i0 + _step)
            _e = cp.asarray(eps[_i0:_i1])
            _sv = (cp.zeros_like(_e) if sigma is None
                   else cp.asarray(sigma[_i0:_i1]) * sc.dt
                   / (2.0 * cpml.EPSILON_0))
            _d = _e + _sv + cp.asarray(gg[_i0:_i1])
            _ca[_i0:_i1, :, :nz] = ((_e - _sv) / _d).astype(cp.float32)
            _cb[_i0:_i1, :, :nz] = ((sc.dt / cpml.EPSILON_0)
                                    / _d).astype(cp.float32)
            del _e, _sv, _d
        return (_ca, _cb)
    s_ = np.zeros_like(eps) if sigma is None else sigma * sc.dt / (2.0 * cpml.EPSILON_0)
    den = eps + s_ + gg
    # Pad before uploading (values unchanged, the padding is 0 and is never read)
    return (_f32(_pad_z_to((eps - s_) / den, nz, nzp)),
            _f32(_pad_z_to((sc.dt / cpml.EPSILON_0) / den, nz, nzp)))


def _apply_pec(sc, ca, cb, d: Dims) -> None:
    """Zero the E field update coefficients of PEC cells (modifies ``ca`` / ``cb`` in place).

    Every source injection (plane wave, TFSF, mode, dipole) is multiplied by cb, so they all go
    to zero with it: no new kernel is needed and there is no second code path (with no PEC this
    block never runs at all). The mask comes from exactly one place, the ε test in
    scene/media.py; PEC combined with dispersion already fails closed in build, so nothing here
    can collide with the cells of _curl_passthrough.
    """
    nx, ny, nz, nzp = d
    for c, pidx in enumerate((sc.pec_ex, sc.pec_ey, sc.pec_ez)):
        if pidx is None or pidx.size == 0:
            continue
        pidx = np.ascontiguousarray(pidx, dtype=np.int64)
        if pidx.min() < 0 or pidx.max() >= nx * ny * nz:
            raise ValueError(
                f"pec_e{'xyz'[c]} flat index out of range ({pidx.min()}..{pidx.max()}, "
                f"cell count {nx * ny * nz})")
        dev = cp.asarray(_repitch_cells(pidx, ny, nz, nzp).astype(np.int64))
        ca[c].ravel()[dev] = np.float32(0.0)
        cb[c].ravel()[dev] = np.float32(0.0)


def _curl_passthrough(tab, ca, cb, d: Dims) -> None:
    """Set ``ca/cb`` of these cells to 0/1, so update_e passes the stretched curl straight out."""
    ny, nz, nzp = d.ny, d.nz, d.nzp
    for c in range(3):
        m = tab.comp == c
        if not m.any():
            continue
        idx = cp.asarray(_repitch_cells(tab.cell[m], ny, nz, nzp))   # flat index, pitch layout
        ca[c].ravel()[idx] = np.float32(0.0)
        cb[c].ravel()[idx] = np.float32(1.0)


def _cacb_lut(ca, cb, mod):
    """Lossless deduplication of update_e's (ca, cb), returning ``(idx, lut, has_lut)``.

    The pair is bit-packed into a uint64 and deduplicated; with at most 65535 unique pairs,
    update_e reads a uint16 index and looks the pair up (a gather returns the same bit pattern
    as a dense read => bitwise identical). The dense arrays are kept for the injection kernels.
    Time-varying media (``mod``) rewrite ca/cb every step and are excluded.
    """
    idx = [cp.zeros(1, cp.uint16)] * 3
    lut = [cp.zeros(1, cp.float32)] * 6
    if mod is not None or knobs.env("CACB_LUT") == "0":
        return idx, lut, np.int32(0)
    _idx3, _lut6 = [], []
    for _c in range(3):
        _a32 = ca[_c].ravel().view(cp.uint32).astype(cp.uint64)
        _b32 = cb[_c].ravel().view(cp.uint32).astype(cp.uint64)
        try:
            _uniq, _inv = cp.unique((_a32 << 32) | _b32, return_inverse=True)
        except Exception:
            # On a very large domain the temporary buffer of cub's radix sort may fail to
            # allocate (cudaErrorInvalidValue measured on the 243-million-cell PEC case);
            # fall back to dense ca/cb.
            del _a32, _b32
            cp.get_default_memory_pool().free_all_blocks()
            return idx, lut, np.int32(0)
        del _a32, _b32
        if _uniq.size > 65535:
            return idx, lut, np.int32(0)
        _idx3.append(cp.ascontiguousarray(_inv.astype(cp.uint16)))
        _lut6.append(cp.ascontiguousarray(
            (_uniq >> 32).astype(cp.uint32).view(cp.float32)))
        _lut6.append(cp.ascontiguousarray(
            (_uniq & cp.uint64(0xFFFFFFFF)).astype(cp.uint32)
            .view(cp.float32)))
    return _idx3, _lut6, np.int32(1)
