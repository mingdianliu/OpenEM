# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""1D auxiliary grid of the TFSF box: the incident plane wave propagates along a single 1D Yee
line. **Does not import tidy3d.**

Why a 1D auxiliary grid rather than an analytic incident field:

- The 1D line uses **exactly the same discrete operator as the 3D grid** (same dt, same dl along
  z, same ε profile), so the incident field it produces satisfies the discrete Maxwell equations
  of the 3D grid bitwise: the box-face correction and the propagation inside agree by
  construction, and the only leakage left is float32 rounding. The analytic route would have to
  solve the discrete dispersion for k frequency by frequency, and a broadband pulse has a
  different k at every frequency.
- The MultipoleExpansion box **cuts through an infinite substrate**: the incident field is the
  plane wave plus the reflection and transmission at the interface. Give the 1D line an ε(z)
  profile and it gets that right by itself; the analytic route would need hand-written layered
  reflection and transmission, with the discrete dispersion on top.

Termination: a **long enough runway**, no absorbing boundary. The runway is sized as vacuum
speed of light * total steps * 1.05 + 64 cells, so a physical signal cannot reach the end within
n_steps. Outside the light cone the discrete scheme still has a **super-exponentially decaying
precursor** (the domain of influence grows by one cell per step; measured around 1e-271 at the
far end of the runway), so the end-point check is not "exactly 0" but "< 1e-30 * injection
peak": 250 orders of magnitude below the smallest float32 subnormal, so anything that reflects
back is strictly unable to add into the float32 fields in 3D. If the check fails it raises: no
tuning parameters, and falsifiable. A 1D cell is 8 bytes and a step is a few vector
multiply-adds, so even 100k cells x 90k steps takes only minutes.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

from openem import cpml
from openem import knobs

#: Signs of the two TF/SF correction terms, by propagation direction. The 3D single-face
#: injection (sources_setup._plane_source_one) and this 1D line read the **same table** (they
#: used to keep a copy each, and an environment override reached only one of them, leaving the
#: two inconsistent). The ``+1`` entry is the one that was verified term by term; the ``-1``
#: entry was measured on the 3D single-face injection against "directionality in a uniform
#: lossless medium should be ~1e-8", and the 1D line uses the same difference expressions, so it
#: reuses them directly.
SIGNS: dict[int, tuple[float, float]] = {+1: (+1.0, +1.0), -1: (-1.0, -1.0)}

#: Safety factor and constant padding (cells) of the runway length. The numerical group velocity
#: of the 1D Yee scheme does not exceed the vacuum speed of light by much, but that is not argued
#: in detail here: the margin is generous and the zero check at the far end backs it up.
RUNWAY_MARGIN = 1.05
RUNWAY_PAD = 64

#: How many cells at each end are checked when the run finishes, and the relative threshold for
#: "counts as zero" (see the module docstring).
UNTOUCHED_CELLS = 8
TAIL_MAX_REL = 1e-30


def matched_dl(k_hat, dl_axes, dt: float, freq: float,
               eps_r: float = 1.0) -> float:
    """Equivalent cell size of the 1D auxiliary line: it makes the numerical phase velocity along
    k̂ agree with the 3D grid.

    At oblique incidence the 1D line runs along k̂, but its own discrete dispersion differs from
    the 3D one: without matching, the incident field injected on the box faces travels at a
    different speed from the propagation inside the box, the TF/SF cancellation is not clean and
    the field leaks outside the box. The fix is to make the two give the **same k** at a given
    frequency (Taflove's standard treatment).

    Discrete dispersion of the lossless Yee scheme (dl may differ per axis)::

        [sin(ω dt/2)/(c̃ dt)]² = Σ_i [sin(k k̂_i dl_i/2)/dl_i]²      (3D)
        [sin(ω dt/2)/(c̃ dt)]² = [sin(k dl₁/2)/dl₁]²                 (1D)

    ``c̃ = c/√ε_r``. First solve the 3D expression for k (the left side is known, the right side
    increases monotonically in k up to the first peak), then invert the 1D expression for dl₁
    (its right side decreases monotonically in dl₁). Both steps use bisection; the signs at the
    interval end points are checked explicitly below, and if they do not hold it raises instead
    of guessing.

    Args:
        k_hat: unit propagation direction (xyz, norm 1).
        dl_axes: primal spacing of each of the three axes [m] (has to be uniform over the box,
            which the caller guarantees).
        dt: time step [s].
        freq: matching frequency [Hz], normally ``source_time.freq0``.
        eps_r: background relative permittivity.

    Returns:
        dl₁ [m]. When k̂ lies along a coordinate axis this is **exactly equal** to the dl of that
        axis (the two expressions then have the same form).

    Raises:
        ValueError: the bisection interval holds no root (dt too large, frequency too high and so
            on), or the result violates the 1D Courant condition ``c̃·dt ≤ dl₁``.
    """
    k_hat = np.asarray(k_hat, dtype=np.float64)
    dl_axes = np.asarray(dl_axes, dtype=np.float64)
    c = cpml.C_0 / np.sqrt(float(eps_r))
    omega = 2.0 * np.pi * float(freq)
    lhs = (np.sin(omega * dt / 2.0) / (c * dt)) ** 2

    def rhs3(k):
        return float(np.sum((np.sin(k * k_hat * dl_axes / 2.0) / dl_axes) ** 2))

    # Upper bound on k: rhs3 increases monotonically until the sin argument of some axis reaches
    # π/2; take the smallest of those bounds
    nz = np.abs(k_hat) > 1e-15
    k_hi = float(np.min(np.pi / (np.abs(k_hat[nz]) * dl_axes[nz])))
    if rhs3(k_hi) < lhs:
        raise ValueError(
            f"on k∈(0, {k_hi:.4e}] the 3D discrete dispersion never reaches {lhs:.6e} "
            f"(at most {rhs3(k_hi):.6e}): dt or the frequency is beyond what this grid can hold")
    lo, hi = 0.0, k_hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if rhs3(mid) < lhs:
            lo = mid
        else:
            hi = mid
    k = 0.5 * (lo + hi)

    # 1D: sin(k dl₁/2)/dl₁ decreases monotonically in dl₁ over (0, π/k), from k/2 down to 0
    target = np.sqrt(lhs)
    if k / 2.0 < target:
        raise ValueError(
            f"the 1D expression has no solution: k/2={k / 2:.6e} is below the target {target:.6e}")
    lo, hi = 1e-12, float(np.pi / k)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if np.sin(k * mid / 2.0) / mid > target:
            lo = mid
        else:
            hi = mid
    dl1 = 0.5 * (lo + hi)
    if c * dt > dl1:
        raise ValueError(
            f"matched dl₁={dl1:.6e} m violates the 1D Courant limit (c·dt={c * dt:.6e} m)")
    return float(dl1)


def runway_cells(n_steps: int, dt: float, dl: float) -> int:
    """Cells needed for one side of the runway: the signal travels at most c·dt·n_steps (the
    vacuum speed of light is the upper bound)."""
    return int(np.ceil(RUNWAY_MARGIN * n_steps * cpml.C_0 * dt / dl)) + RUNWAY_PAD


#:  CUDA source. The two kernels are the H phase and the E phase of the leapfrog, and each one
#: does "sweep + point injection + store the table row" in one go; see `_gpu_tables` below.
_CU_SRC = r"""
extern "C" __global__ void tfsf1d_h(
    double* __restrict__ hy1, const double* __restrict__ ex1,
    const double* __restrict__ ch, double* __restrict__ hy_tab,
    const double inj_val, const int kh,
    const int a, const int bh, const int col0, const int ncol_h,
    const long row)
{
    const int j = a + blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= bh) return;
    double v = hy1[j] - ch[j] * (ex1[j + 1] - ex1[j]);
    if (j == kh) v += inj_val;              // only this thread writes that cell, no race
    hy1[j] = v;
    const int c = j - col0;                 // store: write this cell's own new value
    if (c >= 0 && c < ncol_h) hy_tab[row * ncol_h + c] = v;
}

extern "C" __global__ void tfsf1d_e(
    double* __restrict__ ex1, const double* __restrict__ hy1,
    const double* __restrict__ ca, const double* __restrict__ cb,
    double* __restrict__ ex_tab,
    const double inj_val, const int inj,
    const int ae, const int b, const int col0, const int ncol_e,
    const long row)
{
    const int j = ae + blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= b) return;
    double v = ca[j] * ex1[j] - cb[j] * (hy1[j] - hy1[j - 1]);
    if (j == inj) v += inj_val;
    ex1[j] = v;
    const int c = j - col0;
    if (c >= 0 && c < ncol_e) ex_tab[row * ncol_e + c] = v;
}
"""

_MOD = None


def _gpu_mod():
    """Compile and cache the RawModule. -fmad=false: numpy does the multiply and the subtract
    separately, so contraction has to be off for the result to stay bitwise identical."""
    global _MOD
    if _MOD is None:
        import cupy as cp
        _MOD = cp.RawModule(code=_CU_SRC, options=("-std=c++17", "--fmad=false"))
    return _MOD


class _Line1D(NamedTuple):
    """Everything the main loop of the 1D auxiliary line needs; the numpy path and the GPU path
    read the same bundle.

    ``ch`` is the per-cell H update coefficient, ``ca`` / ``cb`` are the per-edge E update
    coefficients, ``coef_h`` / ``coef_e`` are the scalars at the two injection points, and
    ``win`` turns on the P51 causal window.
    """

    n1: int
    n_steps: int
    inj: int
    kh: int
    col0: int
    n_col_e: int
    ch: np.ndarray
    ca: np.ndarray
    cb: np.ndarray
    coef_h: float
    coef_e: float
    ex_inc: np.ndarray
    hy_inc: np.ndarray
    win: bool

    def window(self, n: int) -> tuple[int, int, int, int]:
        """Ranges to update on step n, ``(a, bh, ae, b)``: H updates ``[a, bh)`` and E updates
        ``[ae, b)`` (``bh`` / ``ae`` trim the edges exactly as ``hy1[:-1]`` / ``ex1[1:]`` do).

        P51 causal window: at the start of step n the non-zero region reaches at most n cells
        either side of the injection point (a discrete Yee step advances by one cell at most),
        and outside the window everything is identically zero, so skipping ``x -= c*(0-0)``
        leaves the result bitwise unchanged. The half width keeps another RUNWAY_PAD cells of
        margin. The runway is n1 ~ steps * Courant number * 2, so over the first half of the run
        the window is far smaller than the full grid (measured 10.5x on Near2Far). If a sampled
        column is still outside the window, what is read is the initial value 0, the same value
        the full-grid version gives (the signal really has not arrived). With ``win=False`` this
        is just the full grid.
        """
        if self.win:
            w = n + 2 + RUNWAY_PAD
            a = max(0, self.inj - w)
            b = min(self.n1, self.inj + w + 1)
        else:
            a, b = 0, self.n1
        return a, min(b, self.n1 - 1), max(a, 1), b


def _dual_dl(dl1: np.ndarray) -> np.ndarray:
    """Dual spacing, with the same absorbing-end rule as :attr:`openem.grid.Axis.dl_dual`."""
    dld1 = np.empty_like(dl1)
    dld1[1:] = 0.5 * (dl1[1:] + dl1[:-1])
    dld1[0] = dl1[0]
    return dld1


def _check_runway(ex1, hy1, ex_inc, n_steps: int) -> None:
    """Whether the runway is long enough: the :data:`UNTOUCHED_CELLS` cells at each end must
    still "count as zero" (threshold and reasoning in the module docstring). If they have been
    touched, the incident table is already contaminated by a reflection, so it raises."""
    m = UNTOUCHED_CELLS
    tol = TAIL_MAX_REL * float(np.abs(ex_inc[: n_steps + 1]).max())
    for name, a in (("ex1", ex1), ("hy1", hy1)):
        for side, seg in (("low end", a[:m]), ("high end", a[-m:])):
            if float(np.abs(seg).max()) > tol:
                raise ValueError(
                    f"1D auxiliary grid {name}, {side}, {m} cells carry {np.abs(seg).max():.3e} "
                    f"(> {tol:.1e}): the runway is too short, the incident table is unreliable")


def _gpu_tables(ln: _Line1D):
    """GPU version of the main loop. Returns (ex_tab, hy_tab, ex1, hy1), all back on the host."""
    import cupy as cp
    mod = _gpu_mod()
    kh_k = mod.get_function("tfsf1d_h")
    ke_k = mod.get_function("tfsf1d_e")
    n1, n_steps, n_col_e = ln.n1, ln.n_steps, ln.n_col_e
    d_ch, d_ca, d_cb = cp.asarray(ln.ch), cp.asarray(ln.ca), cp.asarray(ln.cb)
    d_ex1, d_hy1 = cp.zeros(n1, cp.float64), cp.zeros(n1, cp.float64)
    d_et = cp.zeros((n_steps + 1) * n_col_e, cp.float64)
    d_ht = cp.zeros(n_steps * (n_col_e + 1), cp.float64)
    nch, nce = np.int32(n_col_e + 1), np.int32(n_col_e)
    ch0, ce0 = np.int32(ln.col0 - 1), np.int32(ln.col0)
    BLK = 256
    for n in range(n_steps):
        a, bh, ae, b = ln.window(n)
        kh_k((max((bh - a + BLK - 1) // BLK, 1), 1, 1), (BLK, 1, 1),
             (d_hy1, d_ex1, d_ch, d_ht,
              ln.coef_h * float(ln.ex_inc[n]), np.int32(ln.kh),
              np.int32(a), np.int32(bh), ch0, nch, np.int64(n)))
        ke_k((max((b - ae + BLK - 1) // BLK, 1), 1, 1), (BLK, 1, 1),
             (d_ex1, d_hy1, d_ca, d_cb, d_et,
              ln.coef_e * float(ln.hy_inc[n]), np.int32(ln.inj),
              np.int32(ae), np.int32(b), ce0, nce, np.int64(n + 1)))
    return (cp.asnumpy(d_et).reshape(n_steps + 1, n_col_e),
            cp.asnumpy(d_ht).reshape(n_steps, n_col_e + 1),
            cp.asnumpy(d_ex1), cp.asnumpy(d_hy1))


def _numpy_tables(ln: _Line1D):
    """numpy version of the main loop, returning the same (ex_tab, hy_tab, ex1, hy1) as
    :func:`_gpu_tables`.

    Each step does the H half step first (inject ``ex_inc[n]`` at cell ``kh``, store the H table
    row), then the full E step (inject ``hy_inc[n]`` at edge ``inj``, store the E table row),
    term for term the same as the two functions in the kernel.
    """
    n1, n_steps, n_col_e = ln.n1, ln.n_steps, ln.n_col_e
    ch, ca, cb, inj, kh = ln.ch, ln.ca, ln.cb, ln.inj, ln.kh
    ex1 = np.zeros(n1, dtype=np.float64)
    hy1 = np.zeros(n1, dtype=np.float64)
    ex_tab = np.zeros((n_steps + 1, n_col_e), dtype=np.float64)
    hy_tab = np.zeros((n_steps, n_col_e + 1), dtype=np.float64)
    se_ = slice(ln.col0, ln.col0 + n_col_e)           # E columns: 3D edges klo..khi
    sh_ = slice(ln.col0 - 1, ln.col0 + n_col_e)       # H columns: 3D cells klo-1..khi
    for n in range(n_steps):
        a, bh, ae, b = ln.window(n)
        # H half step: Hy -= ch * dEx/dz (the sz term of yee.cu, transverse terms are 0)
        hy1[a:bh] -= ch[a:bh] * (ex1[a + 1:bh + 1] - ex1[a:bh])
        hy1[kh] += ln.coef_h * ln.ex_inc[n]
        hy_tab[n] = hy1[sh_]
        # E full step: Ex -= cb * dHy/dz (backward difference / dual spacing)
        ex1[ae:b] = (ca[ae:b] * ex1[ae:b]
                     - cb[ae:b] * (hy1[ae:b] - hy1[ae - 1:b - 1]))
        ex1[inj] += ln.coef_e * ln.hy_inc[n]
        ex_tab[n + 1] = ex1[se_]
    return ex_tab, hy_tab, ex1, hy1


def tables(
    dl1: np.ndarray, eps1: np.ndarray, dt: float,
    inj: int, direction: int,
    ex_inc: np.ndarray, hy_inc: np.ndarray,
    n_steps: int, col0: int, n_col_e: int,
    sigma1: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Run the 1D auxiliary grid and return the two incident field tables on the box columns.

    The 1D layout is isomorphic to the z axis of the 3D grid: ``ex1[q]`` sits on "edge q" and
    ``hy1[q]`` on "cell q", and the update expressions match the z differences of Hy/Ex in
    kernels/yee.cu term for term (the transverse derivatives are 0).

    The order within a step matches the leapfrog of solver.run: H half step first (injecting
    ``ex_inc[n]``), then the full E step (injecting ``hy_inc[n]``). The recorded instants match
    as well:

    - ``ex_tab[n]``  = E^n (the 3D H correction is applied after update_h and consumes the
      incident field at instant E^n)
    - ``hy_tab[n]``  = H^{n+1/2} (the 3D E correction consumes the incident field at the half
      step)

    Args:
        dl1: ``(n1,)`` 1D primal spacing [m].
        eps1: ``(n1,)`` relative permittivity at the edge positions (the 3D eps_ex column).
        inj: 1D injection plane (E edge index).
        direction: +1 / -1.
        ex_inc / hy_inc: ``(n_steps+1,)`` injection tables (scene.sources._incident_tables).
        col0: 1D edge index that the lower box edge klo maps to.
        n_col_e: length of the E column = khi - klo + 1.

    Returns:
        ``(ex_tab, hy_tab)``: ``ex_tab`` has shape ``(n_steps+1, n_col_e)`` and row n holds the
        incident value of E^n on the 3D edges ``klo..khi``; ``hy_tab`` has shape
        ``(n_steps, n_col_e+1)`` and row n holds the incident value of H^{n+1/2} on the 3D cells
        ``klo-1..khi``.

    Raises:
        ValueError: the runway is too short (its ends were touched), or the injection tables are
            not long enough.
    """
    n1 = int(dl1.size)
    if eps1.shape != (n1,):
        raise ValueError(f"eps1 has shape {eps1.shape}, expected ({n1},)")
    if ex_inc.size < n_steps + 1 or hy_inc.size < n_steps + 1:
        raise ValueError(
            f"the injection table holds only {ex_inc.size} points, running {n_steps} steps needs "
            f"at least {n_steps + 1} of them")
    if not (UNTOUCHED_CELLS < inj < n1 - UNTOUCHED_CELLS):
        raise ValueError(f"injection plane {inj} is too close to the 1D end points (n1={n1})")
    if not (col0 > 0 and col0 + n_col_e < n1):
        raise ValueError(f"box columns [{col0}, {col0 + n_col_e}) fall outside the 1D grid (n1={n1})")

    dl1 = np.ascontiguousarray(dl1, dtype=np.float64)
    dld1 = _dual_dl(dl1)
    ch = dt / cpml.MU_0 / dl1               # H update coefficient (per cell)
    # Lossy: s = σ dt/(2ε₀), ca = (ε-s)/(ε+s), cb = (dt/ε₀)/((ε+s)·dl_dual).
    # With σ=0, ca = ε/ε is **exactly 1.0**, and multiplying by 1.0 is exact in IEEE754, so the
    # lossless case stays bitwise identical to what it was before.
    _s = (np.zeros_like(eps1) if sigma1 is None
          else np.ascontiguousarray(sigma1, dtype=np.float64)
          * dt / (2.0 * cpml.EPSILON_0))
    _den = eps1 + _s
    ca = (eps1 - _s) / _den                 # decay factor of the E update (per edge)
    cb = dt / cpml.EPSILON_0 / _den / dld1  # E update coefficient (per edge)
    sh, se = SIGNS[direction]
    kh = inj - 1 if direction > 0 else inj  # H injection cell: as sources_setup._plane_source_one
    coef_h = sh * dt / cpml.MU_0 / dl1[kh]
    coef_e = se * dt / cpml.EPSILON_0 / _den[inj] / dld1[inj]

    ln = _Line1D(n1, n_steps, inj, kh, col0, n_col_e, ch, ca, cb,
                 coef_h, coef_e, ex_inc, hy_inc,
                 win=knobs.env("TFSF_WINDOW") != "0")
    # On small cases the launch overhead outweighs the gain (0.7x measured on PECSphere's 7,494
    # steps), so the GPU is only used for long runs
    _gpu_min = int(knobs.env("TFSF_GPU_MIN"))
    out = None
    if knobs.env("TFSF_GPU") != "0" and n_steps >= _gpu_min:
        try:
            out = _gpu_tables(ln)
        except Exception:
            out = None             # no GPU, or the compile failed -> just go through numpy
    ex_tab, hy_tab, ex1, hy1 = out if out is not None else _numpy_tables(ln)

    _check_runway(ex1, hy1, ex_inc, n_steps)
    return ex_tab, hy_tab
