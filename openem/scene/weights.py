# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Per-cell inversion of the geometric weights. **Geometry only, no poles involved**.

Subpixel averaging mixes two materials into one cell and we need the mixing weights. The
``PolarizedAveraging`` form we measured:

    ε_eff = (1−β)·Σfᵢεᵢ + β/(Σfᵢ/εᵢ),      β = n_i² (normal projection on this component, squared)

β=0 is a parallel interface (pure arithmetic mean, linear in ε), β=1 is a perpendicular interface
(pure harmonic mean, linear in 1/ε).

How: sample ``sim.epsilon`` at nfreq frequencies and solve a least-squares problem per cell. The
design matrix is **the same for every cell**, so one pseudo-inverse plus one matrix multiply solves
the entire field, and the residual doubles as the verdict.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
import tidy3d as td

from openem import knobs
from openem.scene.media import epsilon_complex, media_structures, monitor_freqs
from openem.scene.poles import ANISOTROPIC


#: Relative ceiling on the least-squares residual. ``sim.epsilon`` returns complex64, whose
#: single-precision machine epsilon is 1.2e-07; times the frequency count, plus a decade of margin.
FIT_RTOL = 1e-4

#: Weights below this count as 0. A weight is a geometric fill fraction, and a fill of 1e-6
#: contributes to ε far below the noise.
WEIGHT_FLOOR = 1e-6

#: Tolerance for negative weights. The least-squares noise we observe sits in −1e-6..−1.1e-4
#: (measured across the whole RingResonator domain) and is clamped to zero. The filter/project
#: greyscale of a CustomMedium overshoots into negative fills of order −0.7% on interface cells
#: (measured on Autograd15Antenna; still present after detach, so it is interpolation residual and
#: not a model error), so the bound is set at −2e-2. Anything below that means the fit genuinely
#: needs a negative weight (a gain medium), which raises.
NEG_TOL = 2e-2

#: Separate tolerance for the "per-cell columns" path (Custom dispersive media, _decompose_percell).
#: Negative weights on that path come from interpolation overshoot of the design region's
#: filter/project greyscale on interface cells, which is not a model error, and there is no mixing
#: model to fall back on: Autograd15Antenna still reaches −2.52e-2 after the reference cell was
#: reverted to the conservative choice (an earlier measurement gave −2.8e-2), so exceeding NEG_TOL
#: would kill the whole case. The bound goes to −5e-2 with clamping to zero; once NEG_TOL is
#: exceeded the number of affected cells is printed, and only below −5e-2 does this raise.
PERCELL_NEG_TOL = 5e-2

#: Residual ceiling for the mixing model. Looser than :data:`FIT_RTOL` because (f, β) comes from a
#: one-dimensional scan plus an analytic β, so the resolution of the f grid alone contributes of
#: order 1e-5 (median residual 7.9e-06 measured on a sphere).
MIX_RTOL = 3e-4

#: **Lenient ceiling** on the mixing-model residual. An interface cell occasionally holds three
#: materials (three-phase points on edges and corners), which the weighted sum of the arithmetic
#: and harmonic means of two materials cannot describe, so the residual cannot be pushed below
#: :data:`MIX_RTOL`. As long as only a small fraction of cells is over tolerance and the residual
#: stays within this ceiling, take the (f, β) with the smallest residual, carry on, and print the
#: cell count; above it, or with too many cells over tolerance, still fail closed.
#: MetalOxideSunscreen lands here with 33174/12441600 cells (0.27%), worst residual 8.19e-4.
MIX_RTOL_SOFT = 5e-3

#: Largest fraction of over-tolerance cells still allowed onto the lenient path.
MIX_SOFT_FRAC = 0.02

#: Number of points on the 1-D grid used to invert f. β is analytic for every f (the model is
#: linear in β), so only f is scanned.
MIX_GRID = 2001

#: Ceiling on the condition number of the design matrix. Going over it means two materials have
#: ε(ω) too close together at the sample frequencies, so any split of the weight fits: their
#: **sum** is still right, but which of them gets which pole becomes arbitrary.
COND_MAX = 1e8

#: Element budget for blocking the per-cell least squares (_decompose_percell): each block keeps
#: ``step·nmat·nf`` under it.
_PERCELL_BLOCK_ELEMS = 4_000_000

#: Relative tolerance for merging columns: two columns that differ by a real factor with a relative
#: residual ≤ this value are merged into one. Columns that really are proportional differ by
#: ~1e-15, while physically distinct materials differ by at least the dielectric contrast, which
#: leaves 6 orders of magnitude of margin on either side.
MERGE_RTOL = 1e-9

def _eps_col(m, freqs: np.ndarray) -> np.ndarray:
    """The ``ε(f)`` column of one scalar medium at the sample frequencies (complex128, 1-D)."""
    return np.asarray(m.eps_model(freqs), dtype=np.complex128).ravel()


def column_groups(mediums: list, freqs: np.ndarray) -> list[list[int]]:
    """Group the ε(f) columns by "real multiples of each other"; the first member of each group is
    the representative column.

    Why "proportional" and not merely "identical": the columns of any two **non-dispersive**
    materials are constant vectors and are therefore necessarily real multiples of each other, so
    in such a scene the design matrix is always rank deficient, whatever the sample frequencies
    (this is what ailed PlasmonicNanorodArray). After merging, the whole group's weight is recorded
    on the representative column: the reconstructed ε(f) is **exactly the same** as under any other
    split within the group (the weights were never unique anyway), and downstream the ε_∞/σ/pole
    tables only ever consume linear combinations of the weights.
    """
    cols = [_eps_col(m, freqs) for m in mediums]
    groups: list[list[int]] = []
    for i, e in enumerate(cols):
        for g in groups:
            r = cols[g[0]]
            c = float(np.real(np.vdot(r, e)) / np.real(np.vdot(r, r)))
            if np.max(np.abs(e - c * r)) <= MERGE_RTOL * np.max(np.abs(e)):
                g.append(i)
                break
        else:
            groups.append([i])
    return groups

#: Custom dispersive medium with per-cell coefficients: the explanation printed when the mixing
#: model does not apply (shared by the two guards)
_PERCELL_MSG = (
    "{coord_key}: the scene contains a Custom dispersive medium with **per-cell coefficients**, "
    "and the mixing decomposition model does not apply.\n"
    "The model assumes that every cell is a set of non-negative fill fractions over a few "
    "**discrete** media, whereas a medium like this is **its own medium in every cell** "
    "(measured on an 8³ case: sim.epsilon() had 393 distinct values and the design matrix only "
    "2 columns).\n"
    "Expanding the medium over the corners of its coefficient box does not work either: ε is "
    "affine in eps_inf and in the residues, so the corner columns are **exactly linearly "
    "dependent** on the background column (v₂ = 1.3·v₁ − 0.6·v_bg), the least-squares solution is "
    "not unique, and the one the pseudo-inverse picks carries a negative weight (= a gain "
    "medium).\n"
    "The right approach is to **read the medium's own coefficient arrays directly** (the per-cell "
    "eps_inf and residues are already there) and keep the inversion only for subpixel interface "
    "cells."
)


def media_columns(sim: td.Simulation, axis: int) -> list:
    """All media after de-duplication, background first. All of them are **scalar** media, so
    ``eps_model`` can be called on them directly.

    For an ``AnisotropicMedium`` we take the diagonal component on ``axis`` (0/1/2 → x/y/z):
    ``sim.epsilon(coord_key='Ex')`` returns its xx component for it, so when fitting weights it is
    simply "three ordinary media, one per component". That is also why the pole table is organized
    by (medium, axis).

    De-duplication goes by serialized content, not by object identity: when one material appears on
    several structures it has to collapse into a single column, otherwise the design matrix ends up
    with two identical columns and the pseudo-inverse degenerates.

    A ``Medium2D`` enters the columns only after
    :func:`openem.scene.media.media_structures` has replaced it with its equivalent volume, which
    is the same conversion ``sim.epsilon`` does internally; that is what lets the weights fit.
    """
    out, seen = [], set()
    for med in [sim.medium] + [s.medium for s in media_structures(sim)]:
        if type(med) in ANISOTROPIC:
            med = (med.xx, med.yy, med.zz)[axis]
        key = med.json() if hasattr(med, "json") else repr(med)
        if key in seen:
            continue
        seen.add(key)
        out.append(med)
    return out

#: Ceiling on the relative drift of the per-cell residue strength s_c with frequency. By definition
#: s_c does not depend on frequency (ε is affine in the residues); more drift than this means **the
#: pole positions vary per cell as well**, which is not affine, so reading the coefficients
#: directly no longer holds.
PERCELL_S_RTOL = 1e-6


def percell_indices(mediums: list) -> list[int]:
    """Indices, among the media columns, of Custom dispersive media with per-cell coefficients."""
    out = []
    for i, m in enumerate(mediums):
        pr = getattr(m, "pole_residue", None)
        if pr is None:
            continue
        if hasattr(getattr(pr, "eps_inf", None), "values"):
            out.append(i)
    return out


def _has_percell_dispersive(mediums: list) -> bool:
    """Whether the media columns contain a Custom dispersive medium with per-cell coefficients."""
    return bool(percell_indices(mediums))


def percell_columns(sim: td.Simulation, coord_key: str, freqs: np.ndarray,
                    med) -> np.ndarray:
    """``(nf, ncell)``: **read** this medium's own ε(f) **directly**, interpolated onto the Yee
    coordinates of this component.

    Uses tidy3d's own ``eps_comp_on_grid``, the same source ``sim.epsilon`` uses internally for
    this medium, so on pure cells the two are bitwise identical and the geometric weight comes out
    exactly 1 (residual measured at 1.8e-15).

    This is the "read directly" path: the per-cell eps_inf and residues are already in the medium's
    coefficient arrays, so there is no need to invert them out of ε(f); inversion is only needed on
    subpixel interface cells.
    """
    ax = "xyz".index(coord_key[1])
    coords = sim.grid[coord_key]
    return np.stack([
        np.asarray(med.eps_comp_on_grid(ax, ax, float(f), coords),
                   dtype=np.complex128).ravel() for f in freqs])


def percell_strength(col: np.ndarray, eps_inf: float,
                     ref: np.ndarray) -> np.ndarray:
    """Per-cell residue strength ``s_c``, plus a check that the pole positions are uniform across
    cells.

    ``ε_c(f) = eps_inf + s_c·L(f)``, where ``L(f) = ref(f) − eps_inf`` is the pole shape of the
    reference cell. ε is affine in the **residues**, so ``s_c = (ε_c − eps_inf)/L`` **does not
    depend on frequency**: if it drifts with frequency, the pole **positions** vary per cell too, ε
    is not affine in those, and reading the coefficients directly does not hold.

    Args:
        col: ``(nf, ncell)`` per-cell ε(f) of this medium.
        eps_inf: the uniform ε_∞.
        ref: ``(nf,)`` ε(f) of the reference cell.

    Returns:
        ``(ncell,)`` real strengths.

    Raises:
        NotImplementedError: the pole positions vary per cell.
    """
    L = ref - eps_inf
    keep = np.abs(L) > 0.0
    if not keep.any():
        # The reference cell contributes no dispersion. If the **whole column** is also ≈ eps_inf
        # (as in the "empty antenna" reference run of an inverse design: density all 0, every
        # residue 0), the strengths are legitimately all 0; only when the column still carries a
        # dispersive difference is the shape genuinely unobtainable.
        dev = float(np.max(np.abs(col - eps_inf))) if col.size else 0.0
        if dev <= 1e-9 * max(abs(float(eps_inf)), 1.0):
            return np.zeros(col.shape[1], dtype=np.float64)
        raise NotImplementedError(
            "the reference cell contributes no dispersion, cannot extract a pole shape")
    s = (col[keep] - eps_inf) / L[keep, None]          # (nkeep, ncell)
    s_mean = s.mean(axis=0)
    scale = np.maximum(np.abs(s_mean), 1e-12)
    drift = float(np.max(np.abs(s - s_mean) / scale))
    if drift > PERCELL_S_RTOL:
        raise NotImplementedError(
            f"per-cell residue strength drifts with frequency by {drift:.3e} > "
            f"{PERCELL_S_RTOL:g}: the pole **positions** vary per cell as well. ε is not affine "
            "in the positions, so reading the coefficients directly does not hold; that would "
            "need a separate per-cell ADE table.")
    im = float(np.max(np.abs(s_mean.imag)))
    if im > PERCELL_S_RTOL * max(float(np.max(np.abs(s_mean.real))), 1e-12):
        raise NotImplementedError(
            f"per-cell residue strength has an imaginary part ({im:.3e}), the model does not apply")
    return s_mean.real


def scalar_pole_residue(pr):
    """A ``PoleResidue`` with per-cell coefficients → the scalar ``(eps_inf, poles)`` of the
    **reference cell**.

    The reference cell is the one with the largest first-pole residue ``|c|``. As long as
    `_material_tables` and :func:`percell_strength` use **the same** reference, the strength
    ``s_c`` and the pole table ``b`` follow the same convention (``s_c`` is a ratio relative to the
    reference cell).

    ``eps_inf`` has to be constant across cells: ε is affine in ``eps_inf`` too, but that term is
    carried by the **geometric weight** (``w·eps_inf``), so letting it vary per cell couples it to
    the residue term, leaving two unknowns in one equation.
    """
    ei = np.asarray(getattr(pr.eps_inf, "values", pr.eps_inf), dtype=np.float64)
    if ei.size > 1 and float(np.ptp(ei)) > 1e-12 * max(float(np.max(np.abs(ei))), 1.0):
        raise NotImplementedError(
            "the eps_inf of a Custom dispersive medium with per-cell coefficients varies per cell "
            f"as well ({ei.min():.6g} ~ {ei.max():.6g}): it couples to the residues, which leaves "
            "two unknowns in one equation.")
    eps_inf = float(ei.ravel()[0])
    a0, c0 = pr.poles[0]
    c0v = np.asarray(getattr(c0, "values", c0)).ravel()
    # Pick the reference cell by |c|, not by |Re c|: a lossless Lorentz medium (γ=0) turns into a
    # **purely imaginary** residue once converted to pole-residue form, |Re c| is 0 everywhere and
    # argmax falls on cell 0. CustomMediumTutorial's graded-index lens puts the whole index
    # difference into the oscillator strength and keeps eps_inf uniform, and cell 0 happens to have
    # zero strength, so we got "the reference cell contributes no dispersion" (2026-09-05).
    # To stay conservative, switch to |c| only when |Re c| is ≈0 everywhere (lossless Lorentz);
    # lossy media still pick the reference cell by |Re c| as before. Changing the reference cell
    # for the per-cell CustomMedium of Autograd15Antenna made the fill fraction fit at −2.8e-2
    # (threshold −2e-2) and building the scene failed (2026-09-05)
    re_abs = np.abs(np.real(c0v))
    idx = int(np.argmax(re_abs)) if float(re_abs.max()) > 1e-12 * max(float(np.abs(c0v).max()), 1e-300) else int(np.argmax(np.abs(c0v)))
    poles = []
    for a, c in pr.poles:
        av = np.asarray(getattr(a, "values", a)).ravel()
        cv = np.asarray(getattr(c, "values", c)).ravel()
        poles.append((complex(av[idx if av.size > 1 else 0]),
                      complex(cv[idx if cv.size > 1 else 0])))
    return eps_inf, poles


def _percell_eps_model(eps_inf: float, poles: list,
                       freqs: np.ndarray) -> np.ndarray:
    """``ε(f)`` of the reference cell, paired with the ``(eps_inf, poles)`` from
    :func:`scalar_pole_residue`.

    Evaluated with Tidy3D's own ``PoleResidue`` so that rewriting the formula here cannot drift.
    """
    return np.asarray(
        td.PoleResidue(eps_inf=eps_inf, poles=poles).eps_model(freqs),
        dtype=np.complex128).ravel()


def fit_frequencies(sim: td.Simulation, nmat: int) -> np.ndarray:
    """Sample frequencies: they cover the **monitor** band, 3 more of them than the design matrix
    has columns.

    Pass ``nmat`` as the column count **after collinear columns have been merged** (the number of
    groups from :func:`column_groups`); the sample count exists only to make the residual
    meaningful, so it follows the actual degrees of freedom.

    Taking a few extra is what makes the residual meaningful: at exactly ``nfreq = nmat`` least
    squares can always reach zero residual, the verdict becomes useless, and interface cells
    silently pass as pure arithmetic means.

    **The monitor band must be used, not ``sim.frequency_range``.** The latter is the full width of
    the source spectrum and often spans two or three orders of magnitude (MultipoleExpansion covers
    1.9–1080 THz while the monitors only need 214–545 THz). Those far-out frequencies are well
    beyond the valid range of the material fit, ``eps_model`` cannot be trusted there, and the fit
    on interface cells is thrown off completely: the measured mixing-model residual degraded from
    **2.3e-05** to **4.3e-01**, and the count of cells "needing mixing" was overstated by a factor
    of 4 (11,540 vs 2,970).

    With no frequency-carrying monitor, fall back to the source spectrum and narrow the endpoints
    the way it did before.
    """
    fs = monitor_freqs(sim)
    if fs:
        allf = np.unique(np.concatenate(fs))
        lo, hi = float(allf.min()), float(allf.max())
        if hi > lo:
            return np.linspace(lo, hi, nmat + 3)
        # Single-frequency monitor: open a narrow ±10% window around it, or the design matrix
        # degenerates
        return np.linspace(0.9 * lo, 1.1 * lo, nmat + 3)
    lo, hi = (float(v) for v in sim.frequency_range)
    if not (np.isfinite(lo) and np.isfinite(hi)) or hi <= lo:
        raise ValueError(f"no usable frequency band: {sim.frequency_range}")
    pad = 0.05 * (hi - lo)
    return np.linspace(lo + pad, hi - pad, nmat + 3)

def widen_if_ill_conditioned(freqs: np.ndarray, rep_meds: list, coord_key: str = "") -> np.ndarray:
    """When the monitor band is too narrow, widen the sample window around its center until the
    condition number of the design matrix is back under ``COND_MAX``.

    An adjoint simulation (the one tidy3d autograd builds) only has monitors at the few frequency
    points of the objective function. TidyFab0GC's forward monitors cover 1.5–1.6 µm (condition
    number 4.8e5) while the adjoint only has 1.545–1.555 µm; in a window that narrow the ε(f)
    columns of the same set of media (air, cSi, two SiO2 PoleResidue) are nearly collinear, the
    condition number 1.9e8 goes past the ceiling, and the gradients for the whole notebook cannot
    be computed. This is the same degeneracy the "open a ±10% window for a single-frequency
    monitor" rule in :func:`fit_frequencies` deals with, and it is handled the same way: widen the
    window around its center, trying ±1% up to ±10% in steps and stopping as soon as it is enough
    (the closer the window stays to the monitor band, the smaller the difference between the
    forward and adjoint weights).

    **Only act when the condition number is already over the ceiling**: a forward simulation with
    enough bandwidth is returned unchanged, so existing cases keep the same numbers. The static
    rank deficiency :func:`_static_null_vector` lets through (constant columns that are multiples
    of each other, where widening the window would not help anyway) is also returned unchanged and
    left to :func:`decompose` to handle as before.
    """
    if len(freqs) < 2:
        return freqs
    mat = design_matrix(rep_meds, freqs)
    if np.linalg.cond(mat) <= COND_MAX or _static_null_vector(mat, rep_meds) is not None:
        return freqs
    center = 0.5 * (float(freqs[0]) + float(freqs[-1]))
    for half in (0.01, 0.02, 0.03, 0.05, 0.10):
        wide = np.linspace(center * (1 - half), center * (1 + half), len(freqs))
        cond = float(np.linalg.cond(design_matrix(rep_meds, wide)))
        if cond <= COND_MAX:
            print(f"[openem] {coord_key} monitor band too narrow ({td.C_0 / freqs[-1]:.4f}–"
                  f"{td.C_0 / freqs[0]:.4f} µm), medium design matrix ill-conditioned; widening "
                  f"the sample window to ±{half * 100:.0f}% (condition number {cond:.1e})")
            return wide
    return freqs


def _eps_batch(sim: td.Simulation, coord_key: str, freqs: np.ndarray):
    """Sample ε frequency by frequency, **in parallel batches** (2026-09-05).

    Most of ``epsilon_complex`` is the gencoeffs **subprocess** from tidy3d_extras (a fresh TempDir
    each time, with the frequency in the cache key), so a thread pool is enough to run them in
    parallel and they share no files. The number of lanes is bounded by a memory budget: each lane
    holds one complex128 full-grid array (16 B/cell) at a time, plus the same again as margin; the
    default budget is OPENEM_EPS_MEM_GB=48 and the cap is OPENEM_EPS_WORKERS=6. Run serially,
    CMOSRGBSensor (19 dispersive media → 22 sample frequencies × 3 components) and
    PlasmonicWaveguideCO2Sensor (108 million cells, 8 frequencies × 3) used to spend 40 min without
    finishing the scene. Yielding per batch keeps at most `workers` arrays in memory at any moment,
    only a constant factor away from the original discipline of "sample one, write it into rhs,
    release it immediately".
    """
    from concurrent.futures import ThreadPoolExecutor
    freqs = [float(f) for f in np.atleast_1d(freqs)]
    try:
        ncell = int(np.prod([int(n) for n in sim.grid.num_cells]))
    except Exception:
        ncell = 0
    budget = float(knobs.env("EPS_MEM_GB")) * 1e9
    cap = max(1, int(knobs.env("EPS_WORKERS")))
    workers = max(1, min(cap, len(freqs), int(budget // max(1, ncell * 32))))
    if workers <= 1:
        for f in freqs:
            yield np.asarray(epsilon_complex(sim, coord_key, f))
        return
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for i in range(0, len(freqs), workers):
            futs = [ex.submit(epsilon_complex, sim, coord_key, f) for f in freqs[i:i + workers]]
            for fu in futs:
                yield np.asarray(fu.result())


def design_matrix(mediums: list, freqs: np.ndarray,
                  invert: bool = False) -> np.ndarray:
    """``(2·nfreq, nmat)`` real design matrix: the complex equation is split into a block of real
    rows and a block of imaginary rows.

    Splitting instead of running a complex least squares is what **forces the weights to be real**:
    a fill fraction has to be real, and allowing it to go complex would disguise "does not fit" as
    "fits".
    """
    cols, mags = [], []
    for med in mediums:
        e = _eps_col(med, freqs)
        mags.append(np.abs(e))          # the |ε| the row weights need, same array row_scale builds
        if invert:
            e = 1.0 / e
        cols.append(np.concatenate([e.real, e.imag]))
    mat = np.stack(cols, axis=1)
    return mat * row_scale(mediums, freqs, invert=invert, mags=np.stack(mags))[:, None]

def row_scale(mediums: list, freqs: np.ndarray,
              invert: bool = False, mags: np.ndarray | None = None) -> np.ndarray:
    """``(2·nfreq,)`` row weights = the reciprocal of the largest ``|εᵢ|`` at that frequency.

    ``sim.frequency_range`` often spans two orders of magnitude (the full width of the source
    spectrum), while the ε of a Drude medium grows as 1/ω²: without normalization the low-frequency
    rows would completely dominate the least squares, and the residual verdict would fail along
    with it (the tolerance gets inflated by the largest |ε|). After normalization the fit is done
    in terms of **relative error**.

    ``mags``: the ``(nmat, nf)`` array of ``|ε(f)|`` the caller has already computed
    (:func:`design_matrix` passes it in to save one ``eps_model`` call); computed here when not
    given.
    """
    if mags is None:
        mags = np.stack([np.abs(_eps_col(m, freqs)) for m in mediums])
    if invert:
        mags = 1.0 / mags
    s = 1.0 / np.maximum(mags.max(axis=0), 1e-30)
    return np.concatenate([s, s])

def _decompose_percell(sim: td.Simulation, coord_key: str, freqs: np.ndarray,
                       mediums: list, pc: list[int]) -> tuple:
    """Decomposition when a Custom dispersive medium has per-cell coefficients.

    The only difference from the normal path: the columns of those media **differ per cell**, so we
    take their own ε(f) directly (``percell_columns``) instead of inverting against one
    representative value. The design matrix then varies per cell and a small least-squares problem
    has to be solved for every cell; the weights go back to their original meaning of a geometric
    fill fraction, and stay non-negative.

    The per-cell residue strength ``s_c`` goes to the pole table as the fourth return value
    ``pole_scale`` (``b`` has to be multiplied by it), while the ``eps_inf`` term is carried by the
    geometric weight itself.
    """
    nf = len(freqs)
    rows, shape = [], None
    for o in _eps_batch(sim, coord_key, freqs):
        o = np.asarray(o, dtype=np.complex128)
        if shape is None:
            shape = o.shape          # for the final reshape; re-querying ε reads npz even on a hit
        rows.append(o.ravel())
    obs = np.stack(rows)
    ncell = obs.shape[1]
    nmat = len(mediums)

    cols = np.empty((nmat, nf, ncell), dtype=np.complex128)
    scale = np.zeros((nmat, ncell), dtype=np.float64)
    for i, m in enumerate(mediums):
        if i in pc:
            ci = percell_columns(sim, coord_key, freqs, m)
            if ci.shape != (nf, ncell):
                raise RuntimeError(
                    f"{coord_key} per-cell column shape {ci.shape} does not match the "
                    f"observation {(nf, ncell)}")
            cols[i] = ci
            eps_inf, poles = scalar_pole_residue(m.pole_residue)
            ref = _percell_eps_model(eps_inf, poles, freqs)
            scale[i] = percell_strength(ci, eps_inf, ref)
        else:
            cols[i] = _eps_col(m, freqs)[:, None]
            scale[i] = 1.0

    # Per-cell least squares: A_c (2nf, nmat), b_c (2nf,). Blocked to keep memory in check.
    w = np.empty((nmat, ncell), dtype=np.float64)
    resid = np.empty(ncell, dtype=np.float64)
    step = max(1, int(_PERCELL_BLOCK_ELEMS / max(nmat * nf, 1)))
    for s0 in range(0, ncell, step):
        s1 = min(s0 + step, ncell)
        A = np.concatenate([cols[:, :, s0:s1].real, cols[:, :, s0:s1].imag],
                           axis=1)                       # (nmat, 2nf, n)
        A = np.transpose(A, (2, 1, 0))                   # (n, 2nf, nmat)
        b = np.concatenate([obs[:, s0:s1].real, obs[:, s0:s1].imag], axis=0)
        b = b.T[:, :, None]                              # (n, 2nf, 1)
        G = np.matmul(np.swapaxes(A, 1, 2), A)           # (n, nmat, nmat)
        r = np.matmul(np.swapaxes(A, 1, 2), b)           # (n, nmat, 1)
        G += np.eye(nmat) * (1e-14 * np.trace(G, axis1=1, axis2=2)[:, None, None]
                             / max(nmat, 1))
        wc = np.linalg.solve(G, r)                       # (n, nmat, 1)
        w[:, s0:s1] = wc[:, :, 0].T
        resid[s0:s1] = np.abs(np.matmul(A, wc) - b).max(axis=(1, 2))

    denom = np.maximum(np.abs(obs).max(axis=0), 1.0)
    bad = resid > FIT_RTOL * denom
    if bad.any():
        raise NotImplementedError(
            f"{coord_key}: {int(bad.sum())}/{ncell} cells do not fit even with per-cell columns "
            f"(worst relative residual {float((resid / denom).max()):.3e} > {FIT_RTOL:g}). On "
            "those subpixel interface cells this medium really is mixed with another one, and the "
            "mixing model has not been wired up to per-cell columns yet.")
    if float(w.min()) < -PERCELL_NEG_TOL:
        raise ValueError(
            f"{coord_key}: per-cell columns invert to a geometric fill fraction below "
            f"-{PERCELL_NEG_TOL:g} (smallest {w.min():.3e}); the model or the material columns "
            "are wrong")
    if float(w.min()) < -NEG_TOL:
        print(f"[openem] {coord_key}: per-cell columns invert "
              f"{int((w.min(axis=0) < -NEG_TOL).sum())}/{ncell} cells to a fill fraction below "
              f"-{NEG_TOL:g} (smallest {w.min():.3e}, greyscale interpolation overshoot in the "
              "design region); clamping to zero and carrying on", flush=True)
    np.clip(w, 0.0, None, out=w)
    return (w.reshape((nmat,) + shape), None,
            np.zeros(shape, dtype=bool), scale)


def _static_null_vector(mat: np.ndarray, rep_meds: list) -> np.ndarray | None:
    """Return the null-space vector when the rank deficiency is 1 and the null space only involves
    **non-dispersive** ``td.Medium`` (σ allowed).

    Columns like these all lie in span{[1], [1/f]} (ε = ε_∞ + iσ/(ωε₀)), so three or more of them
    are necessarily linearly dependent (Dispersion, BroadbandPlaneWave and similar cases: ε=1,
    ε=4+σ, ε=8.99+σ). The split of the weight is not unique, but what downstream really consumes,
    ε_∞ = Σwᵢεᵢ and σ = Σwᵢσᵢ, are exactly the coordinates along [1] and [1/f], and those are
    **uniquely determined by the data**, independently of the split. So this kind of rank
    deficiency can be let through: take the minimum-norm solution, then project it back to
    non-negative along the null space (:func:`_project_nonneg_along`).
    Return None as soon as the null space touches any medium with poles: how the poles are split
    changes ε(f) between the sample frequencies, and that is the case COND_MAX is there to catch.
    """
    s = np.linalg.svd(mat, compute_uv=False)
    small = s < s[0] / COND_MAX
    if int(small.sum()) != 1:
        return None
    _, _, vt = np.linalg.svd(mat)
    v = vt[-1]
    support = np.flatnonzero(np.abs(v) > 1e-6 * np.abs(v).max())
    if all(type(rep_meds[i]) is td.Medium for i in support):
        return v / np.linalg.norm(v)
    return None


def _project_nonneg_along(w: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Shift every column along the null-space direction ``v`` so the weights become non-negative
    with the smallest possible shift; cells with no feasible shift are left alone (NEG_TOL or the
    residual then marks them as interface cells and sends them to the mixing model). ``w`` has
    shape (nrep, ncell)."""
    pos, neg = v > 1e-12, v < -1e-12
    lo = (np.max(-w[pos] / v[pos][:, None], axis=0) if pos.any()
          else np.full(w.shape[1], -np.inf))
    hi = (np.min(-w[neg] / v[neg][:, None], axis=0) if neg.any()
          else np.full(w.shape[1], np.inf))
    t = np.where(lo > 0, lo, np.where(hi < 0, hi, 0.0))
    t = np.where(lo <= hi + 1e-12, t, 0.0)
    return w + t[None, :] * v[:, None]


def decompose(
    sim: td.Simulation, coord_key: str, freqs: np.ndarray, mediums: list,
    groups: list | None = None,
) -> tuple:
    """Solve the geometric weight of every medium per cell, ``(nmat, nx, ny, nz)``: the three-tuple
    version of :func:`decompose_full` without its fourth item ``pole_scale`` (older callers and
    tests use this form)."""
    w, mix, harm, _ = decompose_full(sim, coord_key, freqs, mediums, groups)
    return w, mix, harm


class _ArithFit(NamedTuple):
    """Result of the per-cell least squares for the arithmetic mean.

    ``w`` is the weight of every medium (after clamping to zero, ``(nmat, ncell)``), ``wr`` is the
    raw solution on the representative columns (not clamped; the mixing fit picks its material
    pairs from it), ``bad`` marks the cells whose residual is over tolerance or that need a
    negative fill fraction, ``neg_cell`` and ``wmin`` are the diagnostics for negative weights, and
    ``reps`` are the representative columns left after merging collinear ones.
    """

    w: np.ndarray
    wr: np.ndarray
    bad: np.ndarray
    neg_cell: np.ndarray
    wmin: float
    shape: tuple
    ncell: int
    reps: list


def _design_or_null(coord_key: str, mediums: list, rep_meds: list,
                    freqs: np.ndarray):
    """The design matrix ``mat``, plus the null-space vector to project back to non-negative along
    when the rank deficiency is 1 (``null``, None when the condition number is fine).

    Raises:
        NotImplementedError: a Custom dispersive medium with per-cell coefficients (it has its own
            precise diagnosis).
        ValueError: the condition number is still too large after collinear columns were merged, so
            the inverted weights cannot be trusted.
    """
    mat = design_matrix(rep_meds, freqs)
    cond = float(np.linalg.cond(mat))
    if cond <= COND_MAX:
        return mat, None
    if _has_percell_dispersive(mediums):
        raise NotImplementedError(_PERCELL_MSG.format(coord_key=coord_key))
    null = _static_null_vector(mat, rep_meds)
    if null is None:
        raise ValueError(
            f"{coord_key}: after merging collinear columns ({len(mediums)}→{len(rep_meds)}), the "
            f"medium design matrix still has condition number {cond:.2e}, which is too large. "
            "These materials are nearly linearly dependent at the sample frequencies, so the "
            "inverted weights cannot be trusted.")
    print(f"[openem] {coord_key} design matrix has rank deficiency 1 ({len(rep_meds)} columns, "
          f"condition number {cond:.1e}); the null space holds only non-dispersive media (ε_∞ and "
          "σ are still unique): take the minimum-norm solution and project it back to non-negative "
          "along the null space")
    return mat, null


def _sample_rhs(sim: td.Simulation, coord_key: str, freqs: np.ndarray,
                rep_meds: list):
    """Right-hand side of the least squares, ``(rhs, shape, ncell)``: real parts in the first nf
    rows, imaginary parts in the last nf, all multiplied by the same row weights as the design
    matrix.

    Incremental: sample ε one frequency at a time, write it into rhs and release it immediately
    (with a CustomMedium on a large grid, keeping nf full-grid ε arrays alive at once runs out of
    memory). The obs the mixing step needs is resampled later, on the bad cells only.
    """
    nf = len(freqs)
    shape = None
    rhs = None
    for j, o in enumerate(_eps_batch(sim, coord_key, freqs)):
        if shape is None:
            shape = o.shape
            ncell = int(np.prod(shape))
            rhs = np.empty((2 * nf, ncell), dtype=np.float64)
        rhs[j] = o.real.ravel()
        rhs[nf + j] = o.imag.ravel()
        del o
    rhs *= row_scale(rep_meds, freqs)[:, None]     # same row weights as the design matrix
    return rhs, shape, ncell


def _arithmetic_fit(sim: td.Simulation, coord_key: str, freqs: np.ndarray,
                    mediums: list, groups: list) -> _ArithFit:
    """Per-cell least squares for the arithmetic mean (linear in ε); see :class:`_ArithFit`."""
    reps = [g[0] for g in groups]
    rep_meds = [mediums[i] for i in reps]
    mat, null = _design_or_null(coord_key, mediums, rep_meds, freqs)
    rhs, shape, ncell = _sample_rhs(sim, coord_key, freqs, rep_meds)

    wr = np.linalg.pinv(mat, rcond=1.0 / COND_MAX) @ rhs  # (nrep, ncell); min-norm if rank defic
    if null is not None:
        wr = _project_nonneg_along(wr, null)
    resid = np.abs(mat @ wr - rhs).max(axis=0)
    scale = np.maximum(np.abs(rhs).max(axis=0), 1.0)
    bad = resid > FIT_RTOL * scale
    w = np.zeros((len(mediums), ncell), dtype=np.float64)
    w[reps] = wr

    # A weight is a geometric fill fraction and **has to be non-negative**. Unconstrained least
    # squares leaves ±1e-6..1e-4 of noise weight on materials that are not present, and the
    # negative half of that is a reversed oscillator, i.e. a narrowband gain medium: the SiO2
    # background cells of RingResonator carried a cSi pole of order −3e-6, which over the narrow k
    # band where that pole crosses the grid dispersion branch gave ~6e-4 of gain per step and
    # overflowed to NaN after 159,200 steps.
    # More negative than NEG_TOL is no longer noise: the fit genuinely needs a negative weight, so
    # we have to stop and raise. Clamping to zero silently would disguise "the model is wrong" as
    # "it fits".
    # **A fit that needs a negative fill fraction has not fitted.** Those cells are the same kind
    # as the cells with an over-tolerance residual (interface cells), and both go to the mixing
    # model below.
    #
    # The residual alone does not catch them: once there are more material columns, unconstrained
    # least squares can fit an interface cell very well with one small negative weight. Measured on
    # Autograd19ApodizedCoupler (3 material columns, 6 frequency points): 177 cells had residuals
    # **all** below FIT_RTOL, not one of them in `bad`, yet they wanted a weight of −3.99e-02.
    # Under the mixing model the median residual is 3.07e-07 (threshold 3e-04, three orders of
    # magnitude lower), with f ∈ [0.014, 0.877] and β ∈ [0.128, 1.000], all legitimate geometric
    # quantities.
    #
    # This used to raise right here, which killed the whole case before those cells could be
    # handled correctly (it blocked four of them: Autograd19ApodizedCoupler, BlueMicroLED,
    # AnisotropicMetamaterialBroadbandPBS and CMOSRGBSensor).
    # The real backstop is on the mixing-model side: when that does not fit, it raises on the
    # MIX_RTOL residual.
    wmin = float(w.min())
    neg_cell = w.min(axis=0) < -NEG_TOL
    if neg_cell.any() and len(reps) < 2:
        if _has_percell_dispersive(mediums):
            raise NotImplementedError(_PERCELL_MSG.format(coord_key=coord_key))
        raise ValueError(
            f"{coord_key}: {int(neg_cell.sum())} cells invert to a weight below -{NEG_TOL:g} "
            f"(smallest {w.min():.3e}), but only {len(reps)} columns are left after merging "
            "collinear ones, not enough to form the material pair the mixing model needs")
    bad = bad | neg_cell
    # Downstream does not use the arithmetic weights of the mixed cells (dispersion.py zeroes their
    # ca/cb and switches them to the mix path), so clamping here only affects the cells that really
    # use the arithmetic model. After the check above those cells are non-negative anyway, and all
    # that gets clamped is the −1e-6..−1e-4 of least-squares noise.
    np.clip(w, 0.0, None, out=w)
    return _ArithFit(w, wr, bad, neg_cell, wmin, shape, ncell, reps)


def _mix_fit(sim: td.Simulation, coord_key: str, freqs: np.ndarray, mediums: list,
             groups: list, wr: np.ndarray, bad: np.ndarray, nf: int) -> dict:
    """The **interface cells** that do not fit: Tidy3D uses a weighted sum of the arithmetic and
    harmonic means there,
        ε_eff = (1−β)·Σfᵢεᵢ + β/(Σfᵢ/εᵢ),   β = n_i²
    so (f, β) is inverted per cell. The material pair is the two largest weights from the
    arithmetic fit: an interface cell only involves two materials to begin with, and although the
    arithmetic fit has a large residual it still tells the major material from the minor one.
    Returns ``{"cells", "mi", "mj", "f", "beta", "resid"}``; the residual check is left to
    :func:`_report_mix_residual`.
    """
    cols = np.flatnonzero(bad)
    obs_c = np.empty((nf, cols.size), dtype=np.complex128)   # resample on the bad cells only
    for j, o in enumerate(_eps_batch(sim, coord_key, freqs)):
        obs_c[j] = np.asarray(o, dtype=np.complex128).ravel()[cols]
        del o
    order = np.argsort(-np.abs(wr[:, cols]), axis=0)
    gi, gj = order[0], order[1]                           # group indices
    eps_cols = np.stack([_eps_col(m, freqs) for m in mediums])   # (nmat, nf)
    best = _MixBest.empty(cols.size)
    key = gi.astype(np.int64) * len(groups) + gj
    for kk in np.unique(key):
        sel = np.flatnonzero(key == kk)
        ga, gb = int(gi[sel[0]]), int(gj[sel[0]])
        # Try the candidate material pairs in **medium space**: merging erased which material of
        # the group sits in this cell, and the mixing model's f+(1−f)=1 constraint is not
        # equivalent for different members of a group. The residual is the only physical verdict,
        # so take the candidate pair with the smallest residual per cell, with MIX_RTOL as backstop.
        for i in groups[ga]:
            for j in groups[gb]:
                best.take_better(sel, _fit_mix(obs_c[:, sel], eps_cols[i],
                                               eps_cols[j]), i, j)
    worst = float(best.resid.max())
    if worst > MIX_RTOL:
        # The within-group candidates have already been tried; this catches the cells where the two
        # largest weights picked **the wrong group**. Only the cells that still do not fit go
        # through an exhaustive sweep over all material pairs; the residual threshold is not
        # loosened one bit, the candidate list is merely completed.
        retry = np.flatnonzero(best.resid > MIX_RTOL)
        for i in range(len(mediums)):
            for j in range(i + 1, len(mediums)):    # (j,i) differs from (i,j) by f↔1−f
                best.take_better(retry, _fit_mix(obs_c[:, retry], eps_cols[i],
                                                 eps_cols[j]), i, j)
    return {"cells": cols, "mi": best.mi, "mj": best.mj, "f": best.f,
            "beta": best.beta, "resid": best.resid}


def _report_mix_residual(coord_key: str, r_out: np.ndarray, ncell: int,
                         neg_cell: np.ndarray, wmin: float, mediums: list) -> None:
    """Three-way triage of the mixing-model residual: the precise diagnosis for per-cell Custom
    dispersive media, the lenient-band message, or fail closed."""
    worst = float(r_out.max())
    if worst > MIX_RTOL:
        n_bad = int(np.sum(r_out > MIX_RTOL))
        # Triage first: a Custom dispersive medium with per-cell coefficients has its own precise
        # diagnosis, so do not let it degrade into the generic "does not fit" (that message does
        # not tell anyone what to do). The other two guards (condition number, only one column
        # left) follow the same pattern.
        if _has_percell_dispersive(mediums):
            raise NotImplementedError(_PERCELL_MSG.format(coord_key=coord_key))
        if worst <= MIX_RTOL_SOFT and n_bad <= MIX_SOFT_FRAC * ncell:
            print(f"[openem] {coord_key}: {n_bad}/{ncell} cells "
                  f"({n_bad / ncell:.2%}) have a mixing-model residual above {MIX_RTOL:g}, "
                  f"worst {worst:.3e}, still within the lenient ceiling {MIX_RTOL_SOFT:g}; "
                  "taking the (f, β) with the smallest residual and carrying on. Most likely "
                  "three materials squeezed into one cell.")
        else:
            raise NotImplementedError(
                f"{coord_key}: {n_bad}/{ncell} cells do not fit even a weighted sum of the "
                f"arithmetic and harmonic means (worst relative residual {worst:.3e} > "
                f"{MIX_RTOL:g}). The cells that fail are exactly the interface cells subpixel "
                "averaging produced, most likely with more than three materials in one cell. Of "
                f"the cells sent in, {int(neg_cell.sum())} are ones where **the arithmetic fit "
                f"needs a negative fill fraction** (smallest weight {wmin:.3e}).")


def decompose_full(
    sim: td.Simulation, coord_key: str, freqs: np.ndarray, mediums: list,
    groups: list | None = None,
) -> tuple:
    """Solve the geometric weight of every medium per cell, ``(nmat, nx, ny, nz)``.

    Collinear ε(f) columns are first merged into representative columns by :func:`column_groups`
    and only then fitted; the whole group's weight is recorded on the representative column, so the
    reconstructed ε is unchanged and non-representative columns always carry weight 0.

    Returns:
        ``(arithmetic weights, harmonic weights or None, harmonic-cell mask, pole_scale)``. The
        first three are ``(nmat, nx, ny, nz)`` or ``(nx, ny, nz)``; the harmonic weights only mean
        anything on the cells where the mask is true. ``pole_scale`` is non-None only on the
        per-cell path (:func:`_decompose_percell`).

    Raises:
        NotImplementedError: some cells fit neither average.
        ValueError: the design matrix is still ill-conditioned after collinear columns were merged,
            so the weights cannot be trusted.
    """
    pc = percell_indices(mediums)
    if pc:
        # Custom dispersive medium with per-cell coefficients: **read** its own ε(f) **directly**
        # as the column, so the weights go back to meaning a geometric fill fraction.
        return _decompose_percell(sim, coord_key, freqs, mediums, pc)
    if groups is None:
        groups = column_groups(mediums, freqs)
    fit = _arithmetic_fit(sim, coord_key, freqs, mediums, groups)
    # np.clip already ran at the end of _arithmetic_fit (after the neg_cell check, before mix);
    # that order must not be changed
    mix = None
    if fit.bad.any():
        if len(fit.reps) < 2:
            raise NotImplementedError(
                f"{coord_key} has cells that do not fit, but only {len(fit.reps)} columns are "
                "left after merging collinear ones, not enough to form the material pair the "
                "mixing model needs.")
        mix = _mix_fit(sim, coord_key, freqs, mediums, groups, fit.wr, fit.bad,
                       len(freqs))
        _report_mix_residual(coord_key, mix["resid"], fit.ncell, fit.neg_cell,
                             fit.wmin, mediums)
    return (fit.w.reshape((len(mediums),) + fit.shape), mix,
            fit.bad.reshape(fit.shape), None)

class _MixBest(NamedTuple):
    """The per-cell "best so far" mixing-model solution: fill fraction ``f``, mixing coefficient
    ``beta``, relative residual ``resid``, and the material pair ``(mi, mj)`` that achieved it. All
    five arrays share length and indexing and are updated together."""

    f: np.ndarray
    beta: np.ndarray
    resid: np.ndarray
    mi: np.ndarray
    mj: np.ndarray

    @classmethod
    def empty(cls, n: int) -> "_MixBest":
        """The residual starts at ``inf``, so the first candidate is always taken."""
        return cls(np.zeros(n), np.zeros(n), np.full(n, np.inf),
                   np.zeros(n, dtype=np.int64), np.zeros(n, dtype=np.int64))

    def take_better(self, sel, fit, i: int, j: int) -> None:
        """``fit=(f, β, residual)`` is the result for the candidate material pair ``(i, j)`` on the
        batch of cells ``sel``; the cells with a smaller residual are overwritten in place."""
        f_, b_, r_ = fit
        take = r_ < self.resid[sel]
        idx = sel[take]
        self.f[idx], self.beta[idx], self.resid[idx] = f_[take], b_[take], r_[take]
        self.mi[idx], self.mj[idx] = i, j


def _fit_mix(obs: np.ndarray, e1: np.ndarray, e2: np.ndarray):
    """Invert ``(f, β, residual)`` per cell. ``obs`` has shape ``(nfreq, ncell)``.

    The model ``ε = (1−β)·A(f) + β·H(f)`` is **linear in β**, so for every candidate f the β is
    solved analytically (a real least squares) and only f needs a one-dimensional scan. That is far
    more accurate than a two-dimensional grid search: refined down to 1e-3, the 2-D search stalled
    at a residual of 6.6e-04, while this approach reaches 7.9e-06.

    For memory, the loop goes over f one value at a time and keeps only the current best per cell,
    so the peak is ``O(nfreq × ncell)``.
    """
    fs = np.linspace(0.0, 1.0, MIX_GRID)
    best_e = np.full(obs.shape[1], np.inf)
    best_f = np.zeros(obs.shape[1])
    best_b = np.zeros(obs.shape[1])
    inv_obs = 1.0 / np.abs(obs)
    for f in fs:
        ari = f * e1 + (1.0 - f) * e2
        har = 1.0 / (f / e1 + (1.0 - f) / e2)
        diff = har - ari
        den = float(np.sum(np.abs(diff) ** 2))
        if den <= 0.0:
            continue
        res = obs - ari[:, None]
        beta = np.real(np.sum(np.conj(diff)[:, None] * res, axis=0)) / den
        err = np.max(np.abs(res - beta[None, :] * diff[:, None]) * inv_obs, axis=0)
        take = err < best_e
        best_e[take] = err[take]
        best_f[take] = f
        best_b[take] = beta[take]
    return best_f, np.clip(best_b, 0.0, 1.0), best_e
