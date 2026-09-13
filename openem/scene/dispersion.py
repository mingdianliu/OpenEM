# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""**Assemble** the material pole tables and the per-cell geometric weights into the CSR tables the
solver wants.

Division of labour: :mod:`openem.scene.poles` deals only with materials,
:mod:`openem.scene.weights` only with geometry, and this module only assembles.

It produces two tables: :class:`~openem.model.Dispersion` (the poles of eps, for pure material
cells and cells on an interface parallel to an axis) and :class:`~openem.model.DispersionMix` (for
slanted interface cells, the weighted sum of an arithmetic and a harmonic branch).
"""

from __future__ import annotations

import numpy as np

from openem import cpml, knobs
import tidy3d as td

from openem.grid import Grid
from openem.model import Dispersion, DispersionMix
from openem.scene import media as _media
from openem.scene._util import COORD_KEYS, E_KEYS
from openem.scene.poles import (
    DISPERSIVE, pole_residue, sigma_si, split_dc_pole,
    stabilize_lossless, trap_coeffs, trap_from_poles, zeta_poles, check_split,
)
from openem.scene.weights import (
    WEIGHT_FLOOR, column_groups, decompose_full, fit_frequencies, media_columns,
    scalar_pole_residue, widen_if_ill_conditioned)


#: dtype of one component's CSR fragment, ``(comp, cell, pole_ofs count, am1, b, g)``.
_ENTRY_DTYPES = (np.int32, np.int32, np.int64, np.complex128, np.complex128, np.float64)


def _empty_entries() -> tuple[np.ndarray, ...]:
    """A component with no dispersive cells: six empty arrays."""
    return tuple(np.zeros(0, t) for t in _ENTRY_DTYPES)


def _entries_for_component(
    comp: int, w: np.ndarray, disp: list[tuple[int, np.ndarray, np.ndarray]],
    pole_scale: np.ndarray | None = None,
) -> tuple[np.ndarray, ...]:
    """Every dispersive cell of one component, turned into a CSR fragment.

    **Grouped by which set of dispersive materials is active**, then vectorized: within a group
    every cell has the same number of poles in the same arrangement, so tile and repeat do it in
    one go with no per-cell Python loop. The largest example in the validation set has several
    million dispersive cells, where a loop would run all night.
    """
    ncell = w.shape[1] * w.shape[2] * w.shape[3]
    wf = w.reshape(w.shape[0], ncell)
    # For a Custom dispersive medium whose coefficients vary cell by cell, the pole table belongs to
    # a **reference cell**, so the residues are additionally multiplied by that cell's strength s_c.
    # The eps_inf term is not: the geometric weight already carries it.
    sf = (np.ones_like(wf) if pole_scale is None
          else np.asarray(pole_scale, dtype=np.float64).reshape(wf.shape))
    # A dispersive material counts as "active" in a cell when its **actual residue strength** is
    # non-zero, i.e. geometric weight times per-cell strength s_c exceeds the floor, not merely the
    # geometric weight. For a Custom dispersive medium with per-cell coefficients, a cell of zero
    # density has s_c=0, so the residue b = base*w*s_c is identically zero and contributes nothing
    # to eps(omega). Such a cell must not enter the dispersive set, or a plane wave source plane
    # landing on it would be wrongly refused as "inside a dispersive medium"; that is exactly what
    # happened to the normalization run of an all-zero-density design. A non-per-cell medium has
    # sf=1, which degenerates back to the original "geometric weight > floor".
    sig = np.zeros(ncell, dtype=np.int64)
    for bit, (i, _, _) in enumerate(disp):
        sig |= (np.abs(wf[i] * sf[i]) > WEIGHT_FLOOR).astype(np.int64) << bit

    out_cell, out_ofs, out_am1, out_b, out_g = [], [], [], [], []
    for s in np.unique(sig):
        if s == 0:
            continue
        members = [d for bit, d in enumerate(disp) if s >> bit & 1]
        cells = np.flatnonzero(sig == s)
        am1_1 = np.concatenate([d[1] for d in members])            # (npole,)
        b_base = np.concatenate([d[2] for d in members])
        counts = [d[1].size for d in members]
        npole = int(am1_1.size)
        # The material weight for each pole slot of each cell
        wrows = np.stack([wf[d[0]][cells] for d in members])       # (nmat_s, ncells)
        srows = np.stack([sf[d[0]][cells] for d in members])
        wpole = np.repeat(wrows * srows, counts, axis=0).T.ravel()  # (ncells*npole,)
        out_cell.append(cells)
        out_ofs.append(np.full(cells.size, npole, dtype=np.int64))
        out_am1.append(np.tile(am1_1, cells.size))
        out_b.append(np.tile(b_base, cells.size) * wpole)
        # g = 2 sum Re(B), which is linear in the weights once aggregated by material
        per_mat = np.array([2.0 * float(np.sum(d[2].real)) for d in members])
        out_g.append(per_mat @ (wrows * srows))
    if not out_cell:
        return _empty_entries()
    cell = np.concatenate(out_cell).astype(np.int32)
    return (np.full(cell.size, comp, dtype=np.int32), cell,
            np.concatenate(out_ofs), np.concatenate(out_am1),
            np.concatenate(out_b), np.concatenate(out_g))

def _arith_pole_terms(sigma: float, poles: list, wt: float, dt: float):
    """The ``(a sequence, b sequence)`` of the arithmetic branch (the poles of **eps**) at fill
    fraction ``wt``.

    Both sequences have the same length and run in the order "sigma pole first, material poles
    after". The concatenation order *is* the CSR order and must not be changed.

    Conductivity has to become a **pole** again here: the solver fills these cells' ``ca/cb`` with
    0/1, which disables the semi-implicit sigma path. sigma corresponds to a pole at ``s=0`` with
    residue ``sigma/(2*eps0)`` (``sigma = 2*Re(c)*eps0``, see :func:`split_dc_pole`), a pure
    integrator with ``A=1``. The kernel composes the arithmetic branch as ``SP = sum 2*Re(P)``, and
    that single slot gives exactly sigma/eps0. The harmonic branch needs nothing: zeta_poles has
    always folded sigma into the poles of 1/eps. At beta=0 this differs in origin from the E-form
    semi-implicit sigma (an integrator versus a semi-implicit update), and their agreement is held
    down by tests/test_dispersion.py::test_mix_sigma_beta_zero_matches_semi_implicit.
    """
    a_l, b_l = [], []
    if sigma:
        a_dc, b_dc = trap_coeffs([(0.0, sigma / (2.0 * cpml.EPSILON_0))], dt)
        a_l.append(a_dc)
        b_l.append(b_dc * wt)
    if poles:
        a_, b_ = trap_coeffs(poles, dt)
        a_l.append(a_)
        b_l.append(b_ * wt)
    return a_l, b_l


def _harm_pole_terms(zeta, wt: float, dt: float):
    """The ``(a sequence, b sequence)`` of the harmonic branch (the poles of **1/eps**) at fill
    fraction ``wt``.

    Marginally stable zeta poles, the +-i*omega_L of a lossless material (see poles.zeta_poles),
    used to be rejected here on the grounds that long-term boundedness in fp32 had not been
    verified. **It was verified, they are bounded, and they are now allowed through.**

    Why it was a concern: the trapezoidal rule maps the imaginary axis exactly onto the unit circle
    (|A| = 1), and the continuous system is marginally stable to begin with. Let fp32 rounding push
    |A| above 1 and after N steps the growth is |A|^N.

    Why it turns out to be fine: ``trap_coeffs`` stores **A-1** rather than A (see
    kernels/dispersion.cu). At small omega*dt, A is about 1 + i*omega*dt, and storing A-1 puts the
    rounding error at the scale of omega*dt instead of 1. Measured in fp32, |A|-1 comes out between
    1e-13 and 6e-9, far below the naive estimate of 6e-8.

    How it was measured end to end: an undriven free oscillation, stepped the way the kernel does it
    (``P += (A-1)*P``) entirely in complex64 for **250,000 steps**, which is more than the longest
    run in the validation set at 248,863, measuring |P_N/P_0|::

        omega*dt   0.005    0.05     0.2      0.5      1.0 (extreme)
        ratio      0.99999  0.99996  1.00011  1.00018  1.00161

    Worst case 1.0016. **No clamping is applied**: the smallest safe clamp would introduce an
    artificial decay of about 1e-7 per step, losing 2.5% over 250k steps, which is worse than this
    0.16% growth. tests/test_zeta_stability.py pins the property down.
    """
    _, pp, rr = zeta
    if not pp.size:
        return [], []
    a2, b2 = trap_from_poles(pp, rr * wt, dt)
    return [a2], [b2]


def _mix_entries(comp: int, mixp: dict, kept: list, sigma_of: np.ndarray,
                 eps_inf_of: np.ndarray, zeta_of: list, dt: float) -> dict:
    """The slanted-interface cells of one component, turned into the raw material for a
    :class:`DispersionMix`.

    The arithmetic branch uses the poles of **eps** (:func:`_arith_pole_terms`) and the harmonic
    branch the poles of **1/eps** (:func:`_harm_pole_terms`), with the residues of both scaled by
    the fill fraction.

    The four material tables are passed flat rather than as one group, because
    tests/test_dispersion.py calls this positionally.
    """
    cells, f = mixp["cells"], mixp["f"]
    beta = np.array(mixp["beta"], dtype=np.float64, copy=True)
    mi, mj = mixp["mi"], mixp["mj"]
    n = cells.size
    p_cnt, q_cnt = np.zeros(n, np.int64), np.zeros(n, np.int64)
    pa_l, pb_l, qa_l, qb_l = [], [], [], []
    eps_inf = np.zeros(n)
    zeta_inf = np.zeros(n)
    for t in range(n):
        i, j, ft = int(mi[t]), int(mj[t]), float(f[t])
        if zeta_of[i] is None or zeta_of[j] is None:
            beta[t] = 0.0            # harmonic branch unusable here: fall back to pure arithmetic
        eps_inf[t] = ft * eps_inf_of[i] + (1 - ft) * eps_inf_of[j]
        zeta_inf[t] = ft / eps_inf_of[i] + (1 - ft) / eps_inf_of[j]
        for m, wt in ((i, ft), (j, 1.0 - ft)):
            if abs(wt) <= WEIGHT_FLOOR:
                continue
            pa, pb = _arith_pole_terms(sigma_of[m], kept[m], wt, dt)
            pa_l.extend(pa)
            pb_l.extend(pb)
            p_cnt[t] += sum(a.size for a in pa)
            if zeta_of[m] is None:
                continue                 # beta already 0, so skip the harmonic branch
            qa, qb = _harm_pole_terms(zeta_of[m], wt, dt)
            qa_l.extend(qa)
            qb_l.extend(qb)
            q_cnt[t] += sum(a.size for a in qa)
    return {
        "comp": np.full(n, comp, np.int32), "cell": cells.astype(np.int32),
        "p_cnt": p_cnt, "q_cnt": q_cnt,
        "pa": np.concatenate(pa_l) if pa_l else np.zeros(0, np.complex128),
        "pb": np.concatenate(pb_l) if pb_l else np.zeros(0, np.complex128),
        "qa": np.concatenate(qa_l) if qa_l else np.zeros(0, np.complex128),
        "qb": np.concatenate(qb_l) if qb_l else np.zeros(0, np.complex128),
        "beta": beta, "eps_inf": eps_inf, "zeta_inf": zeta_inf,
    }

def _staircase_weights(w: np.ndarray, harm: np.ndarray):
    """Experimental knob ``OPENEM_DISP_STAIRCASE=1``: staircase the interface cells of a dispersive
    medium, keeping only the highest-weight material in each cell.

    That is equivalent to the material layout at subpixel=False. Returns ``(w, harm, mixp)``, with
    ``w`` one-hot, ``harm`` all False, and ``mixp`` None, so every mixed cell falls back to the
    arithmetic path.

    One scene diverged with subpixel on and matched the reference with it off (amplitude ratio
    0.997), and this knob was used to localize and work around that instability in dispersive mixed
    cells. A weaker knob once existed that only disabled the harmonic branch on slanted interface
    cells, to isolate the suspicion that the harmonic branch was the unstable one; this knob
    subsumes it and it has been removed.
    """
    w_flat = w.reshape(w.shape[0], -1)
    onehot = np.zeros_like(w_flat)
    onehot[np.argmax(w_flat, axis=0), np.arange(w_flat.shape[1])] = 1.0
    return onehot.reshape(w.shape), np.zeros_like(harm, dtype=bool), None


def _stack_mix(parts: list) -> DispersionMix:
    """Join the fragments of all three components into one :class:`DispersionMix`."""
    def cat(key):
        return np.concatenate([p[key] for p in parts])
    return DispersionMix(
        comp=cat("comp"), cell=cat("cell"),
        p_ofs=np.concatenate([[0], np.cumsum(cat("p_cnt"))]).astype(np.int32),
        pa=cat("pa"), pb=cat("pb"),
        q_ofs=np.concatenate([[0], np.cumsum(cat("q_cnt"))]).astype(np.int32),
        qa=cat("qa"), qb=cat("qb"),
        beta=cat("beta"), eps_inf=cat("eps_inf"), zeta_inf=cat("zeta_inf"))

def _zeta_or_none(eps_inf, kept, sigma):
    """The poles of zeta = 1/eps; returns ``None`` when the fitted eps has a zero in the right half
    plane.

    For such a material, typically a broadband pole fit to tabulated data, as with the metals in
    some scenes, 1/eps is not a stable system and the harmonic branch cannot be used. Returning None
    makes :func:`_mix_entries` fall back to pure arithmetic averaging (beta=0) on **the slanted
    interface cells involving that material**; every other cell still takes the eps branch and is
    unaffected. Only the accuracy of the normal component on slanted interface cells is affected.
    """
    try:
        return zeta_poles(eps_inf, kept, sigma)
    except Exception as e:
        if "right half plane" not in str(e):
            raise
        print(f"[openem] the poles of 1/eps fall in the right half plane; slanted interface cells "
              f"involving this material fall back to arithmetic averaging ({str(e)[:70]})")
        return None


def _material_tables(sim: td.Simulation, mediums: list, freqs: np.ndarray) -> dict:
    """One column of scalar media, turned into every pole table for that column.

    ``mediums`` is selected per axis by :func:`media_columns`, so these tables are naturally a
    mapping from (medium, axis) to pole table: an anisotropic medium gets one per axis, while an
    isotropic one gets the same input on all three axes and produces bitwise identical results.
    """
    prs = [pole_residue(m) for m in mediums]
    # For a Custom dispersive medium with per-cell coefficients, the pole table is the scalar
    # version of a **reference cell**, and the per-cell variation is carried by the pole_scale that
    # weights computes. Both sides pick the reference cell the same way
    # (weights.scalar_pole_residue), which is what makes the convention self-consistent.
    scal = []
    for m, pr in zip(mediums, prs):
        if pr is not None and hasattr(getattr(pr, "eps_inf", None), "values"):
            ei, pl = scalar_pole_residue(pr)
            pr = td.PoleResidue(eps_inf=ei, poles=pl)
            m = pr
        scal.append((m, pr))
    mediums = [x[0] for x in scal]
    prs = [x[1] for x in scal]
    eps_inf_of = np.array([
        float(pr.eps_inf) if pr is not None else float(m.permittivity)
        for m, pr in zip(mediums, prs)])
    # A pole at a=0 inside a dispersive medium is first converted to sigma (see :func:`split_dc_pole`)
    kept, sigma_dc = [], []
    for pr in prs:
        if pr is None:
            kept.append([])
            sigma_dc.append(0.0)
        else:
            k_, s_ = split_dc_pole(pr, float(sim.dt))
            kept.append(k_)
            sigma_dc.append(s_)
    sigma_of = np.array([sigma_si(m, pr, float(freqs[len(freqs) // 2])) + sd
                         for m, pr, sd in zip(mediums, prs, sigma_dc)])
    # The self-check uses the original poles **before any loss is injected**: the correctness of the
    # split, at a 1e-10 tolerance, has nothing to do with stabilization and cannot absorb a
    # perturbation of order gamma_rel ~ 1e-6.
    check_split(mediums, prs, kept, sigma_dc, eps_inf_of, sigma_of, freqs)

    # Inject a tiny artificial loss into lossless poles (poles.stabilize_lossless; bitwise unchanged
    # when gamma_rel=0). A marginally stable pole with |A| identically 1 diverges over a long run on
    # the GPU solver; one case went NaN after 62k steps. The magnitude argument is in
    # poles.LOSSLESS_DAMPING_REL. Everything downstream starts from this stabilized table: the eps
    # branch (trap_coeffs), the 1/eps branch (zeta_poles) and the slanted interface cells (the
    # t["kept"] that _mix_entries receives), so both branches share one convention, and the zeta
    # poles of a lossless material land strictly in the left half plane.
    kept = [stabilize_lossless(k_, dt=float(sim.dt)) for k_ in kept]

    disp = [(i, *trap_coeffs(k_, float(sim.dt)))
            for i, k_ in enumerate(kept) if k_]
    if len(disp) > 63:
        raise NotImplementedError(f"{len(disp)} dispersive materials exceeds the 63-bit limit of the "
                                  "grouping signature")

    zeta_of = [_zeta_or_none(eps_inf_of[i], kept[i], sigma_of[i]) if prs[i] is not None
               or sigma_of[i] else (1.0 / eps_inf_of[i], np.zeros(0, complex),
                                    np.zeros(0, complex))
               for i in range(len(mediums))]
    return {"eps_inf_of": eps_inf_of, "kept": kept, "sigma_of": sigma_of,
            "disp": disp, "zeta_of": zeta_of}


def build(
    sim: td.Simulation, grid: Grid
) -> tuple:
    """``(eps_inf, sigma, Dispersion, DispersionMix | None)`` at the three Yee component positions.

    ``eps_inf`` and ``sigma`` are likewise linear combinations of the weights, so arithmetic
    averaging holds for them too.

    The pole tables are **built once per axis** (:func:`media_columns` takes the diagonal component
    per axis for an anisotropic medium). In an isotropic scene all three axes get the same input and
    produce bitwise identical tables, which is no different from sharing one.

    An **axis with no dispersive medium**, which only arises when an anisotropic medium is
    dispersive on some axes only, does not go through the least-squares inversion: a column of
    entirely static media has a necessarily collinear design matrix, since the columns of lossless
    constant media are multiples of one another. That case is exactly what the already-verified
    static extraction handles, so it calls :func:`openem.scene.media.eps_and_sigma` directly.
    """
    axes: list[dict | None] = []
    for c in range(3):
        mediums = media_columns(sim, axis=c)
        if not any(isinstance(m, DISPERSIVE) for m in mediums):
            axes.append(None)              # static axis: use the static extraction
            continue
        # Group the collinear columns first: the number of sample frequencies must follow the
        # **merged** column count. Which frequencies the probe itself uses does not affect the
        # result, since proportionality is independent of the frequency grid (see
        # weights.column_groups).
        groups = column_groups(mediums, fit_frequencies(sim, len(mediums)))
        freqs = fit_frequencies(sim, len(groups))
        freqs = widen_if_ill_conditioned(freqs, [mediums[g[0]] for g in groups], COORD_KEYS[c])
        axes.append({"mediums": mediums, "freqs": freqs, "groups": groups,
                     **_material_tables(sim, mediums, freqs)})
    if not any(t["disp"] for t in axes if t is not None):
        raise ValueError("dispersion.build was called but the scene has no non-DC poles")
    if not all(axes):
        freq_ref, freq2_ref = _media.reference_freq_pair(sim)

    eps: dict[str, np.ndarray] = {}
    sigma: dict[str, np.ndarray | None] = {}
    parts: list[tuple[np.ndarray, ...]] = []
    mparts: list[dict] = []
    for c, key in enumerate(E_KEYS):
        t = axes[c]
        if t is None:
            eps[key], sigma[key] = _media.eps_and_sigma(
                sim, COORD_KEYS[c], freq_ref, freq2_ref)
            parts.append(_empty_entries())
            continue
        w, mixp, harm, pole_scale = decompose_full(
            sim, COORD_KEYS[c], t["freqs"], t["mediums"], t["groups"])
        if knobs.env("DISP_STAIRCASE") == "1":
            # Compute one version and then override it: the mixed fit for mixp still runs, so its
            # residual check and error behaviour are unchanged, and the result is then overridden by
            # the staircased one. See _staircase_weights.
            w, harm, mixp = _staircase_weights(w, harm)
        eps[key] = np.tensordot(t["eps_inf_of"], w, axes=(0, 0))
        sig = np.tensordot(t["sigma_of"], w, axes=(0, 0))
        sigma[key] = sig if float(np.max(np.abs(sig))) > 0.0 else None
        # Harmonic cells take no part in the E form: the solver fills their ca/cb with 0/1
        w_arith = w.copy()
        w_arith[:, harm] = 0.0
        parts.append(_entries_for_component(c, w_arith, t["disp"], pole_scale))
        if mixp is not None:
            mparts.append(_mix_entries(c, mixp, t["kept"], t["sigma_of"],
                                       t["eps_inf_of"], t["zeta_of"],
                                       float(sim.dt)))

    comp = np.concatenate([p[0] for p in parts])
    cell = np.concatenate([p[1] for p in parts])
    per = np.concatenate([p[2] for p in parts])
    d = Dispersion(
        comp=comp, cell=cell,
        pole_ofs=np.concatenate([[0], np.cumsum(per)]).astype(np.int32),
        am1=np.concatenate([p[3] for p in parts]),
        b=np.concatenate([p[4] for p in parts]),
        g=np.concatenate([p[5] for p in parts]),
    )
    d.validate(grid.shape)
    mm = None
    if mparts:
        mm = _stack_mix(mparts)
        mm.validate(grid.shape, d)
    return eps, sigma, d, mm
