# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Main loop. **Method B**: the kernels live in CUDA, and Python drives every step.

Timing of one step (standard leapfrog; it matches the launch order in ``_step_body``)::

    update_h     advance E^n to H^{n+1/2} (or update_h_absorb / update_he_fused / lean variant)
    inject_h     magnetic dipoles, plane wave Ex_inc(n·dt), TFSF H correction, mode source H side
    absorb(H)    absorber decay of the H family (launched on its own when not fused into
                 update_h)
    disp_pre     two-pass path only: advance the polarization to P_half with E^n, write the
                 history term
    update_e     after mix_pre / modulate_coeffs / tensor_dot, advance H^{n+1/2} to E^{n+1}
    inject_e     plane wave Hy_inc((n+½)·dt), TFSF E correction, mode source E side, electric
                 dipoles
    disp_post    after mix_post: complete P^{n+1} with E^{n+1}; the fused path does post(n) and
                 pre(n+1) in one pass (dispersion_step*), and the last step falls back to a
                 plain post
    absorb(E)    absorber decay of the E family, placed after every E-side update
    sample/dft   time-domain sampling; E^{n+1} @ (n+1)dt and H^{n+1/2} @ (n+½)dt accumulate into
                 the DFT in batches; finally the device step counter is incremented
"""

from __future__ import annotations

import time
from types import SimpleNamespace

import cupy as cp
import numpy as np

from openem import cpml, shutoff as shutoff_mod
from openem import flux as flux_mod
from openem import knobs
from openem.coeffs import _cacb_lut, _pml_axis_arrays
from openem.device import Kernels, grid_1d, grid_2d, grid_3d
from openem.device_tables import (
    _alloc_fields,
    _alloc_fields_c,
    _axis_device_tables,
    _bloch_e_coeffs,
    _dispersion_tables,
    _dual_volumes_eps,
    _e_coeff_tables,
    _mix_table,
    _plane_wave_tables,
    _step_clock,
)
from openem.grid import cyclic_axes
from openem.fusion import (
    DispPlan,
    _abs_e_fold,
    _absorber_batch,
    _ade_fuse_allowed,
    _dispersion_args,
    _graph_allowed,
    _h_fuse_allowed,
    _lean_allowed,
    _lean_report,
    _no_h_source,
    _plain_step,
    _pole_bucket,
    _pole_kernels,
)
from openem.model import COMPONENTS, Scene
from openem.readout import (
    FlushKernels,
    _flush_args,
    _flush_box,
    _flush_plane,
    _flush_tail,
    _readout,
    _sample_flux_time,
)
from openem.results import Phasors, Result
from openem.dispersion_setup import (
    LeanTables,
    _absorber_disp_decay,
    _absorber_slabs,
    _ade_defer_setup,
    _lean_setup,
)
from openem.monitors_setup import (
    _DFT_TMOD,
    _bloch_field_monitors,
    _bloch_flux_monitors,
    _fmon_launch_tables,
    _monitor_buffers,
    _share_snapshots,
    _tmon_batch,
    _tmon_entry,
)
from openem.pitch import _d2h_pinned, _plain_neighbors
from openem.sources_setup import (
    _dipole_tables,
    _mode_launch_tables,
    _modulation_phase_table,
    _modulation_setup,
    _plane_source,
    _tfsf_box,
)

#: Shutoff check interval. Tidy3D checks every 400 steps (docs/architecture.md §4.1).
SHUTOFF_INTERVAL = 400


def _sum2(w, a) -> float:
    """Volume-weighted |u|² reduction. ``a`` is a field component: an array on the real path, a
    ``(real part, imaginary part)`` pair on the complex Bloch path.

    The reduction **must be done in float64**: cupy sums a float32 array into a float32 by
    default, and multiplying that by ε₀=8.85e-12 pushes it below float32's smallest normal
    number, 1.2e-38, where it underflows to 0. The ratio is then 0 < threshold, so "the field
    has decayed" is declared by mistake. Measured consequence: DielectricMetasurfaceAbsorber
    stopped as early as step 5200, while the peak of its source envelope is at step 83,400, so
    the source amplitude was then only 3.1e-05 of the peak.
    """
    if isinstance(a, tuple):
        re, im = a
        return float(cp.sum(w * (re * re + im * im), dtype=cp.float64))
    return float(cp.sum(w * a * a, dtype=cp.float64))


def _shutoff_setup(P, policy, mons, vols, eps_dev) -> SimpleNamespace:
    """Early-termination context, built once during setup: the bundle of criteria ``pol`` (by
    default :func:`shutoff.default_policy`), the source end step, the time-monitor guard, and
    everything that does not change from check to check (dual volume elements, ε, plane
    monitors, cell count, step limit). ``trace`` collects ``(step, field_decay.ratio)`` at every
    check and goes into ``Result.decay_trace``. The stepping loop only passes ``(step, e, h)``
    (see :func:`_early_stop`).

    ``mons`` / ``vols`` / ``eps_dev`` are passed separately: on the Bloch path the flux criterion
    gets an empty monitor list, and the volume elements and ε are computed on the spot (they are
    not part of :class:`StepPlan`); everything else is taken from ``P``.
    """
    sc, verbose = P.sc, P.verbose
    return SimpleNamespace(
        sc=sc,
        pol=policy or shutoff_mod.default_policy(sc.shutoff, interval=SHUTOFF_INTERVAL),
        source_end=shutoff_mod.source_end_step(sc),
        time_mon_guard=_time_mon_guard(sc),
        vols=vols, eps_dev=eps_dev, mons=mons, nz=P.nz, nzp=P.nzp,
        n_total=P.n_total, verbose=verbose, trace=[],
    )


def _early_stop(S: SimpleNamespace, step: int, e, h) -> bool:
    """On every step that is a multiple of ``pol.interval``, evaluate the three criteria
    (:func:`_shutoff_check`) and print the reason when the decision is to stop; every other step
    returns False immediately."""
    if step % S.pol.interval:
        return False
    stop = _shutoff_check(S, step, e, h)
    if stop and S.verbose:
        print(f"  early termination at step {step}: {S.pol.reason}"
              f" (stopping is only allowed after the source ends at step {S.source_end})")
    return stop


def _shutoff_check(S: SimpleNamespace, step: int, e, h) -> bool:
    """Evaluate the three early-termination criteria for this step, ask ``S.pol`` whether to
    stop, and return ``stop``.

    ``e`` / ``h`` are each three field components (see the two forms in :func:`_sum2`), and
    ``S.vols`` holds the matching dual volume elements. The three criteria: (1) the
    volume-weighted instantaneous Σ|E|² over the whole domain, in Tidy3D's convention; (2) the
    total energy including H (whose envelope does not oscillate at 2ω); (3) convergence of the
    observable itself, by turning the current phasors into a flux (when ``S.mons`` is empty this
    one is an empty array, which is what the Bloch path does).
    """
    pol, vols, mons = S.pol, S.vols, S.mons
    m_e2 = sum(_sum2(v, a) for v, a in zip(vols, e))
    if knobs.env("SHUTOFF_DEBUG"):
        _shutoff_debug(step, pol, S.sc, e, m_e2, S.nz, S.nzp)
    w_e = sum(_sum2(v * ep, a) for v, ep, a in zip(vols, S.eps_dev, e))
    w_h = sum(_sum2(v, a) for v, a in zip(vols, h))
    m_u = 0.5 * (cpml.EPSILON_0 * w_e + cpml.MU_0 * w_h)
    # **axis must be passed**: every buffer is allocated in the shape of its own normal
    # ((n1, n2) are the two transverse axes), and the default axis=2 would multiply the
    # (ny, nz) data of an x normal by an (nx, ny) area element. The box-face monitors of
    # PlasmonicYagiUda blew up right here, and only with use_shutoff=True, so a diagnostic
    # script with early termination switched off could not reproduce it. The convergence
    # criterion does not need transverse: it only asks whether this number is still changing,
    # and the integral over the whole face converges just as well.
    # The complex value is assembled on the GPU (re + 1j·im is exactly complex(re, im), bitwise
    # identical) so there is only one D2H. The D2H lands in a preallocated **pinned** buffer;
    # a non-pinned one has to go through a staging buffer, measured at only 1.55 GB/s. Not a
    # single byte changes ⇒ bitwise safe.
    fl = (np.concatenate([
        flux_mod.plane_flux(_d2h_pinned(m["re"] + 1j * m["im"], m), S.sc.grid,
                            axis=m["axis"])
        for m in mons]) if mons else np.zeros(0))
    # **The peak has to be tracked starting from step 0**, so the check runs as usual; only the
    # decision to stop is held back until source injection is over. Skipping the whole check
    # would miss the ramp-up part of the peak and inflate the ratio (measured on
    # test_all_pml_scene_actually_absorbs: 6e-3 turned into 1.4e-2).
    stop = (pol.check(step, m_e2, m_u, fl) and step > S.source_end
            and step >= S.time_mon_guard)
    S.trace.append((step,
                    pol.field_decay.ratio if pol.field_decay else float("nan")))
    if S.verbose and len(S.trace) % 25 == 0:
        print(f"  step {step:7d} / {S.n_total}  {pol.summary()}", flush=True)
    return stop


def _shutoff_debug(step, pol, sc, e, m_e2, nz, nzp):
    """Adjudication experiment: the readings of three candidate weightings (dual volume
    weighted / unweighted / primal volume weighted) at the same instant. Only runs when
    ``OPENEM_SHUTOFF_DEBUG`` is on."""
    if not hasattr(pol, "_dbg_peaks"):
        pol._dbg_peaks = [0.0, 0.0, 0.0]
    _u = sum(_sum2(1.0, a) for a in e)
    # Third candidate (primal volume weighted): grid never had a yee_primal_volume method, so
    # that branch was always false and the reading was always 0. Written as a plain 0 here; the
    # output format is unchanged.
    _pw = 0.0
    _cands = [m_e2, _u, _pw]
    for _ci in range(3):
        pol._dbg_peaks[_ci] = max(pol._dbg_peaks[_ci], _cands[_ci])
    print("SHUTDBG|%d|%.6e|%.6e|%.6e" % (
        step,
        _cands[0] / pol._dbg_peaks[0] if pol._dbg_peaks[0] else 0,
        _cands[1] / pol._dbg_peaks[1] if pol._dbg_peaks[1] else 0,
        _cands[2] / pol._dbg_peaks[2] if pol._dbg_peaks[2] else 0),
        flush=True)




def _launch_mode_inject(kern, fields, ms_launch, mode_srcs, *, side, step, n32,
                        cb=None) -> None:
    """Shared double loop for mode source injection, used in all four places: H side / E side,
    real / complex Bloch.

    ``side="h"``: the H step injects the **E profile**, so it pairs with the whole-step amp_e.
    Wired the other way round (the E profile with amp_h, the H profile with amp_e) the two sides
    are each half a step off, in opposite directions, for a net phase error of ω·dt. That does
    not diverge, it only raises the backward leakage from 1e-8 to 3e-2.
    ``side="e"``: the E step injects the **H profile**, pairs with the half-step amp_h
    (tabulated: the kernel reads row *step), and passes two extra cb tables.

    ``fields`` is a tuple of field arrays per axis: ``((hx,), (hy,), (hz,))`` on the real path,
    ``((hxr, hxi), (hyr, hyi), (hzr, hzi))`` on the Bloch path. Both kernels take their
    arguments in the same order: field of t1, field of t2, [cb_t1, cb_t2], the four profiles,
    the amplitude table, the step, the two coefficients, the plane index, dims, and the three
    trailing arguments.
    One launch per injection term (multiple faces and multiple sources superpose by
    themselves); within a term, the base term plus the broadband corrections. The injection
    kernel accumulates, so calling it once per term sums them naturally; with a single source
    and a single profile both loops run exactly once, bitwise identical to before the change.
    """
    amp_i = 0 if side == "h" else 1
    prof = (("ey_re", "ey_im", "ez_re", "ez_im") if side == "h"
            else ("hy_re", "hy_im", "hz_re", "hz_im"))
    coef_key, plane_key = ("coef_h", "kh") if side == "h" else ("coef_e", "ks")
    for ms, lc in zip(mode_srcs, ms_launch):
        t1, t2 = cyclic_axes(int(ms["axis"]))
        head = (*fields[t1], *fields[t2])
        if cb is not None:
            head = (*head, cb[t1], cb[t2])
        for bp, ba in ((ms, (lc.amp_e, lc.amp_h)[amp_i]),
                       *zip(ms.get("bb", ()), [b[amp_i] for b in lc.bb])):
            kern(lc.grid, lc.block, (*head,
                                     bp[prof[0]], bp[prof[1]],
                                     bp[prof[2]], bp[prof[3]],
                                     ba, step,
                                     ms[coef_key][0], ms[coef_key][1],
                                     ms[plane_key], *lc.dims, *n32))


def _launch_absorb(kabs3, absb, abs_slabs, f3, is_e: int, n32) -> None:
    """The three absorber-decay branches (1D premultiplied profile / single batched launch /
    slab by slab) written in one place, called once each for the H side, the E side, and the re
    and im parts on the Bloch path.

    ``kabs3 = (kabsb1d, kabsb, kabs)``; which path is available depends on whether the handle is
    None (only the batched one is wired up for Bloch). ``f3`` is the three field components, and
    ``is_e`` picks the E-family or H-family profile.
    """
    kabsb1d, kabsb, kabs = kabs3
    if absb is not None:
        if kabsb1d is not None:
            kabsb1d(absb["grid"], absb["block"],
                    (*f3, *absb["g6"], np.int32(is_e),
                     *absb["o"], *absb["e"],
                     *n32))
        else:
            kabsb(absb["grid"], absb["block"],
                  (*f3, absb["di"], absb["dh"], absb["ax"],
                   absb["i0"], absb["j0"], absb["k0"],
                   absb["bi"], absb["bj"], absb["bk"],
                   np.int32(is_e), absb["ns"], *absb["o"], *absb["e"],
                   *n32))
    else:
        for sl in abs_slabs:
            kabs(sl["grid"], sl["block"],
                 (*f3, sl["decay_int"], sl["decay_half"], np.int32(is_e),
                  sl["axis"], *sl["org"], *sl["box"],
                  *n32))


def _launch_tfsf_corr(kern, corr, fields, tab, step, cb=None) -> None:
    """Shared loop for the TFSF box-face correction: the H step (``cb=None``) consumes the
    incident field at time E^n (row *step of ``ex_tab``), while the E step additionally passes
    the per-cell cb of that component and consumes the incident field at time H^{n+1/2} (row
    *step of ``hy_tab``). For the descriptor fields see sources_setup._tfsf_box."""
    ncol = np.int32(tab.shape[1])
    for d in corr:
        head = ((fields[d["comp"]],) if cb is None
                else (fields[d["comp"]], cb[d["comp"]]))
        kern(d["launch"][0], d["launch"][1],
             (*head, d["coef"], tab, step, ncol,
              d["c0"], d["cp"], d["cq"],
              d["base"], d["sp"], d["sq"], d["np"], d["nq"]))


def _launch_p2(kgen, kp2, p2, args_pre, args_mid, args_post) -> None:
    """The 2-pole bucket takes the float2 fast path ``kp2``; the odd buckets at either end go
    through the generic kernel ``kgen`` (sliced range pointers, see fusion._pole_kernels). The
    three segments are launched in the order pre -> p2 -> post."""
    if p2["n_pre"]:
        kgen(p2["g_pre"][0], p2["g_pre"][1], args_pre)
    kp2(p2["g2"][0], p2["g2"][1], args_mid)
    if p2["n_post"]:
        kgen(p2["g_post"][0], p2["g_post"][1], args_post)



def _time_mon_guard(sc: Scene) -> int:
    """Early-termination guard: a time-domain monitor has to record **at least one** sample,
    otherwise monitors that take a snapshot at run_time/2 and at run_time end up with an empty
    time axis (NanobeamCavity's fieldProfileMon only starts at step 269,133; AndersonLocalization
    broke the same way).

    It only guards against an empty axis, not against truncation: truncation is the existing
    semantics of Tidy3D's solver (its early termination only looks at field decay and keeps
    whatever was recorded by then), and demanding a complete recording would make
    BoundaryConditions run 34,000 steps for nothing. OPENEM_TMON_GUARD=full restores "must be
    recorded in full".
    """
    _tm = list(sc.field_time_monitors) + list(sc.flux_time_monitors)
    if knobs.env("TMON_GUARD") == "full":
        return max([int(m.step_end) for m in _tm] + [0])
    return max([int(m.step_begin) + 1 for m in _tm] + [0])


def _warn_nonfinite(phasors, field_phasors) -> None:
    """Divergence sentinel at the exit of the solver data: non-finite values only **warn**, they
    are not blocked and the data is not modified.

    Once a NaN/Inf leaks out silently, it spreads through the notebook's post-processing and
    only blows up far away from where it came from (Autograd20: mode amps of NaN -> np.angle ->
    np.interp -> l_slots=NaN -> td.Box(size=nan) reporting "1 validation error for Box").
    Naming it here at the exit keeps the origin visible at a glance even when the job log has
    been cut down to error lines only. For healthy data this costs one extra sum (a NaN
    propagates into the sum, and so does an Inf; a float64 accumulation cannot overflow
    spuriously), and behavior is bitwise unchanged.
    """
    for tag, d in (("mon", phasors), ("fld", field_phasors)):
        for nm, obj in (d or {}).items():
            a = np.asarray(obj.data if hasattr(obj, "data") else obj)
            if not a.size:
                continue
            s = complex(np.sum(a, dtype=np.complex128))
            if not (np.isfinite(s.real) and np.isfinite(s.imag)):
                bad = int(np.count_nonzero(~np.isfinite(a)))
                print(f"[openem] warning: solver data of monitor {tag}:{nm} contains "
                      f"{bad}/{a.size} non-finite values (NaN/Inf); the field has most likely "
                      "diverged, and downstream post-processing or optimization will fail far "
                      "away from the origin", flush=True)


def _result(S: SimpleNamespace, steps_run: int, triggered: bool, wall: float,
            fields: dict, **monitor_data) -> Result:
    """:class:`Result` assembly shared by both paths. The monitor data (``phasors`` /
    ``field_phasors`` / ``time_samples`` / ``time_slots`` / ``flux_time``) is passed by keyword;
    whatever is not passed takes Result's default (the Bloch path has no flux_time)."""
    return Result(
        steps_run=steps_run,
        shutoff_triggered=triggered,
        decay_trace=S.trace,
        energy_trace=(S.pol.energy_decay.abs_history if S.pol.energy_decay is not None else []),
        source_end=int(S.source_end),
        wall_seconds=wall,
        fields=fields,
        stop_reason=S.pol.reason,
        shutoff_summary=S.pol.summary(),
        **monitor_data,
    )


class StepPlan(SimpleNamespace):
    """Everything a run prepares statically: the device tables of the fields / coefficients /
    sources / monitors, the kernel handles, the argument tuples and the launch geometry. The
    setup functions fill it in (``P.xxx = ...`` or ``P.update(table)``), and the stepping body
    (:func:`_step_body` and the phase functions) only reads it. These are exactly the local
    variables ``run`` used to hold, name for name."""

    def update(self, table: dict) -> None:
        self.__dict__.update(table)


def run(
    sc: Scene,
    num_steps: int | None = None,
    use_shutoff: bool = True,
    kernels: Kernels | None = None,
    verbose: bool = True,
    return_fields: bool = False,
    policy: shutoff_mod.Policy | None = None,
    force_complex: bool = False,
) -> Result:
    """Run one FDTD simulation.

    Args:
        sc: The scene.
        num_steps: Upper bound on the number of steps. Defaults to ``sc.num_time_steps`` (the
            nominal step count).
        use_shutoff: Whether early termination is enabled.
        kernels: Reuse already compiled kernels (saves the NVRTC compile when running several
            times in a row).
        return_fields: Whether to bring the final E/H back as well (for diagnostics).
        policy: Bundle of early-termination criteria, :func:`shutoff.default_policy` by default.
        force_complex: Take the complex-field path even at k=0 (only for testing the criteria:
            at k=0 the complex path must agree with the real one).

    When ``sc.any_bloch`` (the Bloch wave vector is nonzero on some axis), the complex-field
    path of :func:`_run_bloch` is taken automatically; every other scene is unaffected (the
    complex kernels are never launched).
    """
    _ph_t0 = time.time()          # start of the setup phase
    if sc.any_bloch or force_complex:
        return _run_bloch(sc, num_steps=num_steps, use_shutoff=use_shutoff,
                          kernels=kernels, verbose=verbose,
                          return_fields=return_fields, policy=policy)
    k = kernels or Kernels()
    n_total = int(num_steps if num_steps is not None else sc.num_time_steps)
    P = StepPlan(sc=sc, k=k, n_total=n_total, verbose=verbose)

    # ---- Fields and auxiliary fields; the 1D tables of the three axes ----
    P.update(_alloc_fields(sc))
    P.tabs = [a.index_tables() for a in sc.grid.axes]
    P.update(_axis_device_tables(sc, P.tabs, (P.pe_np, P.ph_np), P.dims))
    # ---- Coefficients (all precomputed here; the kernels only multiply and add) ----
    P.update(_e_coeff_tables(sc, P.tabs, P.dims))                       # ch ca cb ten
    P.mod = _modulation_setup(sc, P.ny, P.nz, P.nzp)
    P.update(_dispersion_tables(sc, P.dims))                            # hist disp hist_box cw
    # ---- Sources and monitors: pure setup data, lifted out of run (see each function) ----
    P.plane_src = _plane_source(sc)
    P.tfsf = _tfsf_box(sc, n_total)
    P.update(_dipole_tables(sc, n_total, P.dims))                       # dip dip_h g_dip g_dip_h
    P.mons, P.fmons, P.tmons, P.ftmons, P.dftab = _monitor_buffers(
        sc, P.dims, n_total, batch_cap=_dft_batch_cap(use_shutoff))
    P.update(_dual_volumes_eps(sc, P.dims, gpu=knobs.env("INIT_GPU") != "0"))   # vols eps_dev
    P.update(_mode_launch_tables(sc, P.dims))                           # mode_srcs ms_launch
    P.abs_slabs = _absorber_slabs(sc)
    P.pabs = _absorber_disp_decay(sc)
    P.absb = _absorber_batch(P.abs_slabs, P.dims, verbose)
    S = _shutoff_setup(P, policy, P.mons, P.vols, P.eps_dev)

    # ---- Kernel handles, monitor launch tables, the remaining device tables ----
    P.g3, P.b3 = grid_3d(P.nx, P.ny, P.nz)
    P.g2, P.b2 = grid_2d(P.nx, P.ny)
    _bind_kernels(P)
    # Field monitors on the same box share a snapshot: it depends only on (origin, box), not on
    # the frequency points or the window.
    _share_snapshots(P.fmons, verbose)
    P.hfuse = _h_fuse_allowed(P.absb, P.plane_src, P.tfsf, P.mode_srcs, P.dip_h, verbose)
    P.khf = k["update_h_absorb"] if P.hfuse else None
    _fmon_launch_tables(P.fmons, verbose)
    P.tmon_batch = _tmon_batch(P.tmons, verbose)
    P.ktbatch = k["sample_time_batch"] if P.tmon_batch else None
    P.ten_tables = (*P.nxt, *P.prv, *P.mnx, *P.mpv) if P.ten is not None else None
    P.n32_yz = (np.int32(P.ny), np.int32(P.nzp))   # tensor kernels decode/encode with the pitch
    P.mix = _mix_table(sc, P.dims)
    P.update(_step_clock(n_total, sc.dt))                               # stepdev t_e_dev t_h_dev
    P.pw_all = _plane_wave_tables(P.plane_src, P.ex, P.ey, P.hx, P.hy, P.cb)
    P.mod_cs = _modulation_phase_table(P.mod, n_total, sc.dt)

    # ---- Thresholds, argument tuples and launch geometry of the fused paths. The order is
    # fixed: p2 reorders the disp entries, so it has to run after the other arguments and
    # before the dispersion arguments (see _plan_dispersion) ----
    _plan_update_args(P)
    _plan_ade(P)
    _plan_lean(P)
    _plan_dispersion(P)
    _plan_hefuse(P)
    _plan_lean_args(P)
    _plan_mix_args(P)

    t0 = time.time()
    # The t=0 sample is the zero field before the loop: inside the loop we always hold E^{n+1},
    # so m=0 is out of reach
    _sample_time(P)
    _sample_flux_time(0, P.ftmons, P.e_arr, P.h_arr, P.nz)
    _ph_t1 = time.time()          # end of the setup phase
    P.stepdev[:] = 0          # entering the loop: *step == n
    steps_run, triggered = _time_loop(P, S, use_shutoff)
    wall = time.time() - t0
    _ph_t2 = time.time()          # end of the stepping phase
    out, fphas, tsamp, tfill, ftime, fields = _readout(
        P.mons, P.fmons, P.tmons, P.ftmons, steps_run,
        P.e_arr, P.h_arr, P.nz, return_fields)

    if knobs.env("PHASES") == "1":
        print(f"PHASES|{_ph_t1 - _ph_t0:.3f}|{wall:.3f}"
              f"|{time.time() - _ph_t2:.3f}|{steps_run}", flush=True)
    _warn_nonfinite(out, fphas)
    return _result(S, steps_run, triggered, wall, fields,
                   phasors=out, field_phasors=fphas, time_samples=tsamp,
                   time_slots=tfill, flux_time=ftime)


# ---------------------------------------------------------------------------
# Setup phase: kernel handles and the argument tuples of the fused paths
# (fills StepPlan only, launches nothing)
# ---------------------------------------------------------------------------


def _dft_batch_cap(use_shutoff: bool) -> int:
    """Upper bound on the batch size: an early-termination check step must land on a batch
    boundary (the criteria read the accumulator and need it up to date), so with early
    termination on, ``_DFT_TMOD`` is halved until it divides ``SHUTOFF_INTERVAL``."""
    _cap = _DFT_TMOD
    while _cap > 1 and use_shutoff and SHUTOFF_INTERVAL % _cap:
        _cap //= 2      # a check step must land on a batch boundary (the criteria need it fresh)
    return _cap


def _bind_kernels(P: StepPlan) -> None:
    """Fetch the kernel handles once and hang them on P (``Kernels.__getitem__`` is a dict
    lookup). Only the ones that do not depend on a threshold go here; those that switch with a
    threshold (khf / ktbatch / ke_lut / kdstep* / klean*) are attached by run or by the
    relevant ``_plan_*``."""
    k = P.k
    P.kh = k["update_h"]
    P.ke = k["update_e"]
    # H and E in a single pass (operator fusion plus blocking along i). The H neighbors the E
    # part needs come from shared memory or from the previous iteration's registers, which
    # saves update_e a whole re-read of H; the ring of cells whose neighbors are unavailable is
    # finished off by update_e_lut_edge. Updating E in place is safe: the E read across a block
    # boundary only comes from the three nxt directions, which land exactly on that skipped
    # ring, and nothing in this launch writes there.
    P.khe_f = k["update_he_fused"]
    P.ke_edge = k["update_e_lut_edge"]
    P.kdip = k["inject_dipole"]
    P.kdip_h = k["inject_dipole_h"]
    P.kbox = k["accumulate_dft_box"]
    P.ktab = k["dft_phase_table"]
    P.kabsb = k["absorb_batch"] if P.absb else None
    P.kabsb1d = (k["absorb_batch_1d"]
                 if P.absb is not None and P.absb.get("g6") is not None else None)
    P.ksnapb, P.kflushb = k["dft_snap_box"], k["dft_flush_box"]
    # Shared-memory tiled version (bitwise identical, 2.17x on the benchmark).
    # OPENEM_DFT_SM2=0 falls back.
    P.kflushb_sm2 = (k["dft_flush_box_sm2"]
                     if knobs.env("DFT_SM2") != "0" else None)
    # Plane monitors that need tangential components only use the _sub kernel (with a
    # component map)
    P.ksnapb_sub = k["dft_snap_box_sub"]
    P.kflushb_pf_sub = k["dft_flush_box_pf_sub"]
    P.ksnapp, P.kflushp = k["dft_snap_plane"], k["dft_flush_plane"]
    P.kflushb_pf, P.kflushp_pf = k["dft_flush_box_pf"], k["dft_flush_plane_pf"]
    P.kflush = FlushKernels(P.kflushb, P.kflushb_sm2, P.kflushb_pf,
                            P.kflushb_pf_sub, P.kflushp, P.kflushp_pf)
    P.ktime = k["sample_time_box"]
    P.kih, P.kie = k["inject_h"], k["inject_e"]
    P.ktfe, P.ktfh = k["tfsf_corr_e"], k["tfsf_corr_h"]
    P.kmh, P.kme = k["inject_mode_h"], k["inject_mode_e"]
    P.kabs = k["absorb_slab"]
    P.kabs3 = (P.kabsb1d, P.kabsb, P.kabs)
    P.kdpre, P.kdpost = k["dispersion_pre"], k["dispersion_post"]
    P.kmpre, P.kmpost = k["mix_pre"], k["mix_post"]
    P.kmod = k["modulate_coeffs"]
    P.ktdot, P.ktscat = k["tensor_dot"], k["tensor_scatter"]
    P.kdft = k["accumulate_dft"]
    P.kadv = k["step_advance"]
    P.ke_ade = k["update_e_ade"]


def _plan_update_args(P: StepPlan) -> None:
    """The static arguments ``args_h`` / ``args_e`` of update_h / update_e; whether ca/cb can be
    turned into a LUT (``has_lut``) decides whether the E update uses ``ke_lut`` or ``ke``."""
    P.args_h = (
        P.hx, P.hy, P.hz, P.ex, P.ey, P.ez, *P.psi_h, *P.on32,
        *P.idl_p, *P.pml_h[0], *P.pml_h[1], *P.pml_h[2],
        *P.nxt, *P.mnx, P.ch,
        *P.n32,
    )
    P.cacb_idx, P.cacb_lut, P.has_lut = _cacb_lut(P.ca, P.cb, P.mod)

    # The LUT gets its own kernel (two paths in one kernel make nvcc reorder the FMAs of the
    # dense expression, and the sentinel caught a drift of 1e-7 per step). The arguments take
    # one of the two in the ca/cb slots.
    P.ke_lut = P.k["update_e_lut"] if P.has_lut else None
    P.args_e = (
        P.ex, P.ey, P.ez, P.hx, P.hy, P.hz, *P.psi_e, *P.on32,
        *((*P.cacb_idx, *P.cacb_lut) if P.has_lut else (*P.ca, *P.cb)),
        *P.idl_d, *P.pml_e[0], *P.pml_e[1], *P.pml_e[2],
        *P.prv, *P.mpv, *P.pec,
        *P.hist, P.cw, *P.hist_box,
        *P.n32,
    )


def _plan_ade(P: StepPlan) -> None:
    """Full ADE fusion (the pole update is inlined into update_e and hist disappears entirely).

    update_e_ade has no hist argument ⇒ it is only usable when **every** entry can be deferred,
    that is, a dispersive scene with no absorber layers (entries inside a layer have to stay in
    dispersion_step_dec).
    """
    ade = None
    if _ade_fuse_allowed(P.disp, P.has_lut, P.pabs, P.ten, P.mod, P.mix):
        ade = _ade_defer_setup(P.disp, P.pabs, P.dims)
        if ade is not None and ade["n_def"] != int(P.disp["n"]):
            ade = None            # an entry cannot be deferred ⇒ its hist is unread, fall back
    P.ade = ade
    P.args_ade = None
    if ade is not None:
        disp = P.disp
        P.args_ade = P.args_e[:-13] + (
            disp["p_re"], disp["p_im"],
            disp["lut_am1_re"], disp["lut_am1_im"],
            disp["lut_b_re"], disp["lut_b_im"], disp["cidx"],
            *ade["estart"], *ade["ecount"], *ade["eb"], P.cw,
            *P.n32)
        if P.verbose:
            print(f"  [P43] dispersion folded into update_e: all {ade['n_def']:,} entries inlined")


def _plan_lean(P: StepPlan) -> None:
    """Interior / boundary-shell split (measured). If the threshold is not met, lean=None and
    everything falls back."""
    P.lean = None
    if _lean_allowed(P.has_lut, P.ade, P.mod, P.ten, P.disp):
        P.lean = _lean_setup(
            LeanTables(pml_e=P.pml_e, pec=P.pec, prv=P.prv, mpv=P.mpv,
                       pml_h=P.pml_h, nxt=P.nxt, mnx=P.mnx),
            P.abs_slabs, P.disp, P.pabs, P.dims)
    _lean_report(P.lean, P.verbose)


def _plan_dispersion(P: StepPlan) -> None:
    """Thresholds and arguments of the dispersion kernels: pole bucketing, the in-layer decay
    factor, the two post+pre fused paths (``disp_fused`` / ``disp_fused_dec``) and their static
    argument tuples."""
    # Must happen **before** any args tuple is built: those tuples capture array pointers.
    P.p2 = _pole_bucket(P.disp, P.ade, P.lean, P.ten, P.pabs, P.verbose)

    # The in-layer decay of P is folded into dispersion_post: with no in-layer dispersion entry
    # a placeholder array plus has_dec=0 is passed and the kernel does not even read it
    # (results bitwise identical to before)
    P.dec = P.pabs["dense"] if P.pabs is not None else cp.zeros(1, cp.float32)
    P.has_dec = np.int32(1 if P.pabs is not None else 0)
    # Conditions for the post+pre fusion (for the proof of bitwise equivalence see
    # dispersion_step in kernels/dispersion.cu): no dispersion inside the layers, no tensor.
    # When it applies, one full read/write pass over P is saved per step.
    P.disp_fused = P.disp is not None and P.pabs is None and P.ten is None
    # In-layer dispersion (pabs) also takes the fused single pass. Tensor scenes stay excluded
    # (pass2 modifies E).
    P.disp_fused_dec = P.disp is not None and P.pabs is not None and P.ten is None
    k = P.k
    P.kdstep = k["dispersion_step"] if P.disp_fused else None
    P.kdstepd = k["dispersion_step_dec"] if P.disp_fused_dec else None
    P.kdstepd_we = k["dispersion_step_dec_we"] if P.disp_fused_dec else None
    P.absorb_e_fold = _abs_e_fold(P.disp_fused_dec, P.absb, P.disp,
                                  P.dims, P.verbose)
    P.dplan = DispPlan(P.disp, P.hist, P.e_arr, P.dec, P.has_dec, P.pabs,
                       P.disp_fused, P.disp_fused_dec)
    P.args_dpre, P.args_dstep, P.args_dstepd, P.args_dpost = _dispersion_args(P.dplan)
    (P.kdstep_p2, P.kdstepd_p2, P.kdstepd_we_p2, P.args_dstep_p2, P.args_dstepd_p2,
     P.args_dstep_pre, P.args_dstep_post, P.args_dstepd_pre,
     P.args_dstepd_post) = _pole_kernels(P.p2, k, P.dplan)


def _plan_hefuse(P: StepPlan) -> None:
    """Thresholds and parameters of the H/E fusion.

    Fusion is only used when nothing at all happens between update_h and update_e.
    **Off** by default: measured about 20% slower than launching the two separately. It is kept
    because it is bitwise correct and guarded by tests/test_he_fuse.py, so on new hardware or
    with a different implementation idea it can be measured again with OPENEM_HE_FUSE=1,
    without a rewrite.
    """
    nx, ny, nz, nzp = P.nx, P.ny, P.nz, P.nzp
    nxt, prv = P.nxt, P.prv
    _FBJ, _FBK = 8, 32                # must match BJ/BK in yee_fused.cu
    _FNI = int(knobs.env("HE_NI"))
    P.hefuse = (P.ke_lut is not None and P.lean is None
                and knobs.env("HE_FUSE") == "1"
                and _no_h_source(P.plane_src, P.tfsf, P.mode_srcs, P.dip_h)
                and _plain_step(P.absb, P.abs_slabs, P.disp, P.mix, P.mod, P.ten)
                and nx >= 2 * _FNI and _plain_neighbors(nxt, prv, (nx, ny, nz)))
    P.args_hef = P.args_edge = None
    if P.hefuse:
        P.g_hef = ((nzp + _FBK - 1) // _FBK, (ny + _FBJ - 1) // _FBJ,
                   (nx + _FNI - 1) // _FNI)
        P.b_hef = (_FBK, _FBJ, 1)
        P.args_hef = ((P.hx, P.hy, P.hz, P.ex, P.ey, P.ez) + P.args_h[6:-3] + P.args_e[6:-3]
                      + (*P.n32, np.int32(_FNI)))
        # Which cells still need to be filled in: the strict complement of the write-E
        # condition in yee_fused.cu
        _ii = cp.arange(nx, dtype=cp.int64)[:, None, None]
        _jj = cp.arange(ny, dtype=cp.int64)[None, :, None]
        _kk = cp.arange(nzp, dtype=cp.int64)[None, None, :]
        _pl = ((nxt[0][_ii] == _ii + 1) & (prv[0][_ii] == _ii - 1)
               & (nxt[1][_jj] == _jj + 1) & (prv[1][_jj] == _jj - 1)
               & (nxt[2][_kk] == _kk + 1) & (prv[2][_kk] == _kk - 1))
        _inner = ((_jj % _FBJ > 0) & (_kk % _FBK > 0) & (_ii % _FNI > 0) & _pl)
        _eidx = cp.flatnonzero(~cp.broadcast_to(_inner, (nx, ny, nzp)).ravel()
                               ).astype(cp.int32)
        del _ii, _jj, _kk, _pl, _inner
        P.args_edge = (*P.args_e, _eidx, np.int32(_eidx.size))
        P.g_edge = grid_1d(int(_eidx.size))
        if P.verbose:
            print(f"  [P41] H/E fusion: blocks {_FBJ}x{_FBK}, {_FNI} layers along i per block; "
                  f"edge pass covers {_eidx.size / (nx * ny * nzp):.0%} of the cells")


def _plan_lean_args(P: StepPlan) -> None:
    """Parameters of lean mode (static tuples, capturable by a Graph): the lean interior kernel,
    the lean H side and the boundary shell (a single kernel), then on to
    :func:`_plan_lean_disp_args`. With ``lean is None`` all of them stay None."""
    P.klean = P.args_lean = None
    P.klean_h = P.args_lean_h = None
    P.keshell = P.args_eshell = P.khshell = P.args_hshell = None
    P.lean_dskin = P.lean_dskin_dec = P.lean_dpost_skin = P.lean_dpost_int = None
    lean = P.lean
    if lean is None:
        return
    k, nzp = P.k, P.nzp
    P.klean = k["update_e_ade_lean"]
    _t1 = cp.zeros(1, cp.float32)
    _d = P.disp if P.disp is not None else {
        "p_re": _t1, "p_im": _t1, "lut_am1_re": _t1, "lut_am1_im": _t1,
        "lut_b_re": _t1, "lut_b_im": _t1, "cidx": cp.zeros(1, cp.int32)}
    _est = lean.get("estart", (cp.zeros(1, cp.uint32),) * 3)
    _ecn = lean.get("ecount", (cp.zeros(1, cp.uint8),) * 3)
    _eb = lean.get("eb", (np.int32(0),) * 6)
    P.args_lean = (
        P.ex, P.ey, P.ez, P.hx, P.hy, P.hz,
        *P.cacb_idx, *P.cacb_lut, *P.idl_d,
        P.pml_e[0][2], P.pml_e[1][2], P.pml_e[2][2],
        _d["p_re"], _d["p_im"],
        _d["lut_am1_re"], _d["lut_am1_im"],
        _d["lut_b_re"], _d["lut_b_im"], _d["cidx"],
        *_est, *_ecn, *_eb, P.cw,
        lean["o"][0], lean["o"][1], lean["oz_launch"], lean["oz_real"],
        lean["e"][0], lean["e"][1], lean["e"][2],
        *P.n32,
    )
    if lean.get("h"):
        _lh = lean["h"]
        P.klean_h = k["update_h_lean"]
        P.args_lean_h = (
            P.hx, P.hy, P.hz, P.ex, P.ey, P.ez, *P.idl_p,
            P.pml_h[0][2], P.pml_h[1][2], P.pml_h[2][2], P.ch,
            _lh["o"][0], _lh["o"][1], _lh["oz_launch"], _lh["oz_real"],
            _lh["e"][0], _lh["e"][1], _lh["e"][2],
            *P.n32,
        )
        P.khshell = k["update_h_shell"]
        _hbox = (_lh["o"][0], _lh["o"][1], _lh["o"][2],
                 _lh["e"][0], _lh["e"][1], _lh["e"][2])
        P.args_hshell = [(P.g3, P.b3, P.args_h[:-3] + _hbox
                          + (np.int32(0), np.int32(nzp)) + P.args_h[-3:])]
    # One kernel for the boundary shell (launch over the whole domain, skip the interior),
    # replacing six fragmented box launches
    P.keshell = k["update_e_lut_shell"]
    _ebox = (lean["o"][0], lean["o"][1], lean["o"][2],
             lean["e"][0], lean["e"][1], lean["e"][2])
    # Measured: a three-segment launch (73% fewer threads) made it worse, 1.23x -> 1.02x. The
    # xy frame of seg2 has a bad shape (most of a block still skips) plus two extra launches.
    # Reverted to a single launch that passes the full-domain z range.
    P.args_eshell = [(P.g3, P.b3, P.args_e[:-3] + _ebox
                      + (np.int32(0), np.int32(nzp)) + P.args_e[-3:])]
    _plan_lean_disp_args(P)


def _plan_lean_disp_args(P: StepPlan) -> None:
    """Launch parameters for the boundary-shell dispersion entries in lean mode:
    ``lean_dskin`` / ``lean_dskin_dec`` (the post+pre fusion on every step but the last) and
    ``lean_dpost_skin`` / ``lean_dpost_int`` (the last step falls back to a plain post and needs
    both). Without dispersion all four stay None."""
    lean, disp, pabs = P.lean, P.disp, P.pabs
    if disp is None:
        return
    ni, ns = lean["n_int"], lean["n_skin"]
    if ns:
        _sl = slice(ni, None)
        _launch_s = grid_1d(ns)
        _skin_common = (disp["p_re"], disp["p_im"], *P.hist,
                        disp["lut_am1_re"], disp["lut_am1_im"],
                        disp["lut_b_re"], disp["lut_b_im"], disp["cidx"],
                        disp["comp"][_sl], disp["cell"][_sl],
                        disp["ofs"][ni:], P.ex, P.ey, P.ez)
        P.lean_dskin = (_launch_s, (*_skin_common, np.int32(ns)))
        if pabs is not None:
            P.lean_dskin_dec = (_launch_s, (*_skin_common,
                                            P.dec[_sl], pabs["f1"][_sl],
                                            pabs["f2"][_sl], pabs["f3"][_sl],
                                            np.int32(ns)))
        P.lean_dpost_skin = (_launch_s, (
            disp["p_re"], disp["p_im"], disp["lut_b_re"], disp["lut_b_im"],
            disp["cidx"], disp["comp"][_sl], disp["cell"][_sl],
            disp["ofs"][ni:], P.ex, P.ey, P.ez,
            P.dec[_sl] if pabs is not None else P.dec, P.has_dec,
            np.int32(ns)))
    if ni:
        P.lean_dpost_int = (grid_1d(ni), (
            disp["p_re"], disp["p_im"], disp["lut_b_re"], disp["lut_b_im"],
            disp["cidx"], disp["comp"], disp["cell"],
            disp["ofs"], P.ex, P.ey, P.ez, P.dec, P.has_dec, np.int32(ni)))


def _plan_mix_args(P: StepPlan) -> None:
    """Static arguments of mix_pre / mix_post; None when there is no mix."""
    mix, ex, ey, ez = P.mix, P.ex, P.ey, P.ez
    P.args_mpre = (
        (mix["p_re"], mix["p_im"], mix["q_re"], mix["q_im"], mix["sp"], mix["sq"],
         mix["pa_re"], mix["pa_im"], mix["pb_re"], mix["pb_im"],
         mix["qa_re"], mix["qa_im"], mix["qb_re"], mix["qb_im"],
         mix["p_ofs"], mix["q_ofs"], mix["comp"], mix["cell"], mix["d_h"],
         ex, ey, ez, mix["n"]) if mix else None
    )
    P.args_mpost = (
        (ex, ey, ez, mix["d_tot"], mix["d_h"], mix["p_re"], mix["p_im"],
         mix["q_re"], mix["q_im"], mix["pb_re"], mix["pb_im"],
         mix["qb_re"], mix["qb_im"], mix["p_ofs"], mix["q_ofs"],
         mix["comp"], mix["cell"], mix["sp"], mix["sq"],
         mix["k1"], mix["k2"], mix["omb"], mix["invd"], mix["dte"], mix["n"])
        if mix else None
    )


# ---------------------------------------------------------------------------
# Stepping body: the phases of one time step
# (reads StepPlan only, launched in the order of the module docstring)
# ---------------------------------------------------------------------------


def _sample_time(P: StepPlan) -> None:
    """Sampling. The slot and phase logic is derived inside the kernel from the device step
    counter (``m_step = *step + 1``), so both launches go out unconditionally here and the
    kernel guards itself: the same sequence every step, which a graph can capture. ``filled`` is
    instead recomputed at the end of run from the same formula.
    """
    fields = (*P.e_arr, *P.h_arr)
    if P.tmon_batch is not None:
        _tb = P.tmon_batch
        for phase in (0, 1):
            P.ktbatch(_tb["launch"][0], _tb["launch"][1],
                      (_tb["outs"], *fields,
                       _tb["comps"], _tb["w_b"] if phase else _tb["w_a"],
                       _tb["nc"], P.stepdev,
                       _tb["beg"], _tb["iv"], _tb["end"], _tb["slots"],
                       np.int32(phase),
                       _tb["i0"], _tb["j0"], _tb["k0"],
                       _tb["ni"], _tb["nj"], _tb["nk"], _tb["cells"],
                       _tb["n"], *P.n32))
    else:
        for tm in P.tmons:
            for phase in (0, 1):
                g, b = tm["launch"]
                i0, j0, k0 = tm["origin"]
                P.ktime(g, b, (tm["buf"], *fields, tm["comps"],
                               tm["w_b"] if phase else tm["w_a"], tm["nc"],
                               P.stepdev, np.int32(tm["beg"]), np.int32(tm["iv"]),
                               np.int32(tm["end"]), np.int32(tm["slots"]),
                               np.int32(phase),
                               i0, j0, k0,
                               *tm["box"], *P.n32))


def _launch_h_phase(P: StepPlan) -> None:
    """H half step: update_h (one of four paths) -> H-side source injection (magnetic dipoles,
    plane wave, TFSF, mode sources) -> absorber decay of the H family (already done inside
    update_h_absorb when fused)."""
    if P.hefuse:
        # H and E in one pass, then fill in the ring of cells whose neighbors are unavailable
        P.khe_f(P.g_hef, P.b_hef, P.args_hef)
        P.ke_edge(P.g_edge[0], P.g_edge[1], P.args_edge)
    elif P.hfuse:
        # curl update and absorber decay in one pass (saves a whole read/write pass over H)
        absb = P.absb
        P.khf(P.g3, P.b3, P.args_h[:-3]
              + (absb["di"], absb["dh"], absb["ax"],
                 absb["i0"], absb["j0"], absb["k0"],
                 absb["bi"], absb["bj"], absb["bk"], absb["ns"])
              + P.args_h[-3:])
    elif P.args_lean_h is not None:
        P.klean_h(P.lean["h"]["grid"], P.lean["h"]["block"], P.args_lean_h)
        for g_, b_, a_ in P.args_hshell:
            P.khshell(g_, b_, a_)
    else:
        P.kh(P.g3, P.b3, P.args_h)
    if P.dip_h is not None:
        # Magnetic dipoles are injected with the H step, amplitude from the whole-step table
        # (for the timing see the comments in inject_dipole_h)
        dh = P.dip_h
        P.kdip_h(P.g_dip_h[0], P.g_dip_h[1],
                 (P.hx, P.hy, P.hz, P.ch, dh["flat"], dh["comp"], dh["coef"],
                  dh["src_of"], dh["amp"], P.stepdev,
                  dh["nsteps1"], dh["n"]))
    for _px, _ext, _hyt, _ph, _pe, _pcb in P.pw_all:
        P.kih(P.g2, P.b2, (_ph, _px["coef_h"], _ext, P.stepdev,
                           _px["kh"], *P.n32))
    if P.tfsf is not None:
        # The H-step correction consumes the incident field at time E^n, same timing as
        # single-plane injection
        _launch_tfsf_corr(P.ktfh, P.tfsf["h_corr"], P.h_arr, P.tfsf["ex_tab"], P.stepdev)
    if P.mode_srcs is not None:
        _launch_mode_inject(P.kmh, P.h3, P.ms_launch, P.mode_srcs,
                            side="h", step=P.stepdev, n32=P.n32)
    if not P.hfuse:     # when fused, already done inside update_h_absorb
        _launch_absorb(P.kabs3, P.absb, P.abs_slabs, P.h_arr, 0, P.n32)


def _launch_e_phase(P: StepPlan) -> None:
    """E half step: disp_pre of the two-pass path, mix_pre, time-varying ca/cb, the first tensor
    pass -> update_e (one of four paths) -> E-side source injection (plane wave, TFSF, mode
    sources, electric dipoles)."""
    disp, mix, mod, ten = P.disp, P.mix, P.mod, P.ten
    if disp is not None and P.ade is None and not (P.disp_fused or P.disp_fused_dec):
        P.kdpre(disp["launch"][0], disp["launch"][1], P.args_dpre)
    if mix is not None:
        P.kmpre(mix["launch"][0], mix["launch"][1], P.args_mpre)
    if mod is not None:
        # The cos/sin table is precomputed in float64 outside the loop and then cast down to
        # f32 (the same as converting step by step, bitwise identical)
        P.kmod(mod["launch"][0], mod["launch"][1],
               (*P.ca, *P.cb, mod["amp"], mod["cph"], mod["sph"], mod["eps"],
                mod["comp"], mod["cell"],
                P.mod_cs, mod["sd"], mod["gg"], P.stepdev,
                mod["dte"], mod["n"]))
    if ten is not None:
        # First pass: while the E array still holds E^n, store M1·⟨E^n⟩
        P.ktdot(ten["launch"][0], ten["launch"][1],
                (ten["rhs"], ten["rhs"], np.float32(0.0), P.ex, P.ey, P.ez,
                 ten["m1"], ten["comp"], ten["cell"], *P.ten_tables,
                 *P.n32_yz, ten["n"]))
    if P.hefuse:
        pass          # already done in update_he_fused + edge
    elif P.lean is not None:
        # Lean interior kernel (ADE inlined, deferral is bitwise identical) + the six boxes of
        # the boundary shell (a clone of the fat kernel)
        P.klean(P.lean["grid"], P.lean["block"], P.args_lean)
        for g_, b_, a_ in P.args_eshell:
            P.keshell(g_, b_, a_)
    elif P.ade is not None:
        P.ke_ade(P.g3, P.b3, P.args_ade)        # pole update already inlined
    else:
        (P.ke_lut or P.ke)(P.g3, P.b3, P.args_e)
    for _px, _ext, _hyt, _ph, _pe, _pcb in P.pw_all:
        P.kie(P.g2, P.b2, (_pe, _pcb, _px["coef_e"], _hyt, P.stepdev,
                           _px["ks"], *P.n32))
    if P.tfsf is not None:
        # The E-step correction consumes the incident field at time H^{n+1/2}
        _launch_tfsf_corr(P.ktfe, P.tfsf["e_corr"], P.e_arr, P.tfsf["hy_tab"], P.stepdev,
                          cb=P.cb)
    if P.mode_srcs is not None:
        _launch_mode_inject(P.kme, P.e3, P.ms_launch, P.mode_srcs,
                            side="e", step=P.stepdev, n32=P.n32, cb=P.cb)
    if P.dip is not None:
        d = P.dip
        P.kdip(P.g_dip[0], P.g_dip[1],
               (P.ex, P.ey, P.ez, P.cb[0], P.cb[1], P.cb[2], d["flat"], d["comp"], d["coef"],
                d["src_of"], d["amp"], P.stepdev, d["nsteps1"], d["n"]))


def _launch_disp_post(P: StepPlan, last: bool) -> None:
    """The post half step of the dispersion: the ade / lean / fused single pass / two-pass paths
    each launch their own version; on the last step (``last``) all of them fall back to a plain
    post, so P at the end of the loop is bitwise aligned with the two-pass path."""
    disp, lean = P.disp, P.lean
    if P.ade is not None:
        # The pole update is already inlined in update_e_ade (deferred). The last step adds one
        # plain post so that P at the end of the loop is bitwise aligned with the two-pass path.
        if last:
            P.kdpost(disp["launch"][0], disp["launch"][1], P.args_dpost)
    elif lean is not None:
        # The interior entries are already inlined in klean (deferred), so only the boundary
        # shell is handled here; on the last step both groups fall back to a plain post, leaving
        # the final P bitwise aligned with the old path
        if not last:
            if P.lean_dskin_dec is not None:
                P.kdstepd(P.lean_dskin_dec[0][0], P.lean_dskin_dec[0][1],
                          P.lean_dskin_dec[1])
            elif P.lean_dskin is not None:
                P.kdstep(P.lean_dskin[0][0], P.lean_dskin[0][1], P.lean_dskin[1])
        else:
            if P.lean_dpost_skin is not None:
                P.kdpost(P.lean_dpost_skin[0][0], P.lean_dpost_skin[0][1],
                         P.lean_dpost_skin[1])
            if P.lean_dpost_int is not None:
                P.kdpost(P.lean_dpost_int[0][0], P.lean_dpost_int[0][1],
                         P.lean_dpost_int[1])
    elif P.disp_fused and not last:
        # post(n) + pre(n+1) in one pass; the last step falls back to a plain post, so P at the
        # end of the loop is bitwise identical to the two-pass path
        if P.kdstep_p2 is not None:
            # Fast path for the 2-pole bucket + the generic kernel for the odd buckets at either
            # end (sliced range pointers)
            _launch_p2(P.kdstep, P.kdstep_p2, P.p2,
                       P.args_dstep_pre, P.args_dstep_p2, P.args_dstep_post)
        else:
            P.kdstep(disp["launch"][0], disp["launch"][1], P.args_dstep)
    elif P.disp_fused_dec and not last:
        # As above; in a pabs scene it carries the P decay plus a replay of the E-family decay
        # factors. With full coverage the write-back variant is used and the E-side absorb pass
        # below is skipped accordingly
        if P.kdstepd_p2 is not None and not P.absorb_e_fold:
            _launch_p2(P.kdstepd, P.kdstepd_p2, P.p2,
                       P.args_dstepd_pre, P.args_dstepd_p2, P.args_dstepd_post)
        elif P.kdstepd_we_p2 is not None and P.absorb_e_fold:
            _launch_p2(P.kdstepd_we, P.kdstepd_we_p2, P.p2,
                       P.args_dstepd_pre, P.args_dstepd_p2, P.args_dstepd_post)
        else:
            (P.kdstepd_we if P.absorb_e_fold else P.kdstepd)(
                disp["launch"][0], disp["launch"][1], P.args_dstepd)
    else:
        P.kdpost(disp["launch"][0], disp["launch"][1], P.args_dpost)


def _launch_post_e(P: StepPlan, last: bool) -> None:
    """E-side wrap-up: the D->E conversion of mix, the dispersion post, the second tensor pass,
    the absorber decay of the E family, all placed after every E-side injection (for the reasons
    see the comments on each part)."""
    mix, ten = P.mix, P.ten
    # D form: in those cells the E array currently holds "stretched curl + source term"; convert
    # it into the real E. Placed **after** source injection so that current sources enter the D
    # update automatically (see the note on ca/cb).
    if mix is not None:
        P.kmpost(mix["launch"][0], mix["launch"][1], P.args_mpost)
    # **Must come after source injection**: what the source adds is part of E^{n+1} too, and the
    # polarization has to respond to it. Placed before, the polarization in the source cells
    # would miss the source's own contribution: two equivalent scenes were measured 62% apart,
    # against 1e-5 in the correct order.
    if P.disp is not None:
        _launch_disp_post(P, last)
    if ten is not None:
        # The second pass comes **after every E-side injection**: a pass-through cell currently
        # holds "stretched curl − J", so the source current enters the D update automatically
        # (the same mechanism as the ca/cb note for mix). The gather only writes the private
        # enew and only the scatter writes back into E: the neighbor average reads cells that
        # other entries are about to write, so writing in place would create a read/write race.
        P.ktdot(ten["launch"][0], ten["launch"][1],
                (ten["enew"], ten["rhs"], np.float32(1.0), P.ex, P.ey, P.ez,
                 ten["m2"], ten["comp"], ten["cell"], *P.ten_tables,
                 *P.n32_yz, ten["n"]))
        P.ktscat(ten["launch"][0], ten["launch"][1],
                 (P.ex, P.ey, P.ez, ten["enew"], ten["comp"], ten["cell"],
                  *P.pec, *P.n32_yz, ten["n"]))
    # The E-family absorber decay goes **after every E-side update** (source, mix, ADE post), so
    # that the in-layer ADE sees the decayed E. What is approximated is the passive dissipative
    # system "a dispersive medium inside a lossy host", not the mismatched stretched system of
    # the old ψ approach.
    if not (P.absorb_e_fold and not last):
        # with absorb_e_fold the E decay was already written back by dispersion_step_dec_we
        # (the last step behaves as before)
        _launch_absorb(P.kabs3, P.absb, P.abs_slabs, P.e_arr, 1, P.n32)
    # The P of the in-layer ADE decays in lockstep with E (after dispersion_post): this keeps
    # the decaying frame consistent and amounts to damping the in-layer oscillators (passive).
    # The factors are precomputed in _absorber_disp_decay and multiplied in once per entry here.


def _launch_dft(P: StepPlan) -> None:
    """Frequency-domain accumulation: the phase table (shared by every DFT monitor, once per
    step) -> box monitors -> plane monitors."""
    dftab = P.dftab
    # Phase table: shared by every DFT monitor, once per step (time table + device step counter)
    if dftab["n"]:
        P.ktab(dftab["launch"][0], dftab["launch"][1],
               (dftab["tab"], dftab["freqs"],
                P.t_e_dev, P.t_h_dev, P.stepdev, dftab["n"], dftab["tmod"],
                dftab["tab32"]))
    _launch_dft_boxes(P)
    _launch_dft_planes(P)


def _launch_dft_boxes(P: StepPlan) -> None:
    """Field monitors (box): snapshot plus a flush at the batch boundary, or direct
    accumulation."""
    dftab, fields = P.dftab, (*P.e_arr, *P.h_arr)
    for m in P.fmons:
        i0, j0, k0 = m["origin"]
        ni, nj, nk = m["box"]
        if m["T"]:
            # Snapshot plus flush at the batch boundary (on steps inside a batch the flush
            # kernel guards itself and returns immediately)
            if m.get("skip_snap"):
                pass          # the leader on the same box snapshotted; we share its buffer
            elif int(m["ncomp"]) < 6:
                P.ksnapb_sub(m["launch_snap"][0], m["launch_snap"][1],
                             (m["snap"], *fields, P.stepdev,
                              m["T"], m["cmap_dev"], m["ncomp"],
                              i0, j0, k0, ni, nj, nk,
                              *P.n32))
            else:
                P.ksnapb(m["launch_snap"][0], m["launch_snap"][1],
                         (m["snap"], *fields, P.stepdev, m["T"],
                          i0, j0, k0, ni, nj, nk,
                          *P.n32))
            # Same argument list as the tail flush (readout._flush_args) and the same dispatch
            # (readout._flush_box), only with rem=0
            _flush_box(m, P.kflush,
                       _flush_args(m, dftab, P.stepdev, 0, (m["cells"], m["nf"])))
        else:
            gb, bb = m["launch"]
            P.kbox(gb, bb, (m["re"], m["im"], *fields,
                            dftab["tab"],
                            m["win_e_dev"], m["win_h_dev"], P.stepdev,
                            i0, j0, k0,
                            ni, nj, nk, m["nf"],
                            *P.n32,
                            m["f0"], dftab["n"], dftab["tmod"]))


def _launch_dft_planes(P: StepPlan) -> None:
    """Flux monitors (plane): snapshot plus a flush at the batch boundary, or direct
    accumulation."""
    dftab, fields, nz = P.dftab, (*P.e_arr, *P.h_arr), P.nz
    for m in P.mons:
        if m["T"]:
            P.ksnapp(m["launch_snap"][0], m["launch_snap"][1],
                     (m["snap"], *fields, P.stepdev, m["T"],
                      np.int32(m["km"]), np.int32(m["axis"]), np.int32(nz),
                      *P.n32))
            _flush_plane(m, P.kflush,
                         _flush_args(m, dftab, P.stepdev, 0,
                                     (m["n1"], m["n2"], m["nf"])))
        else:
            P.kdft(
                m["launch"][0], m["launch"][1],
                (m["re"], m["im"], *fields, dftab["tab"],
                 m["win_e_dev"], m["win_h_dev"], P.stepdev,
                 np.int32(m["km"]), np.int32(m["axis"]), np.int32(m["nf"]),
                 np.int32(nz),
                 *P.n32,
                 m["f0"], dftab["n"], dftab["tmod"]),
            )


def _step_body(P: StepPlan, last: bool) -> None:
    """Every kernel launch of one time step: the unit a graph captures.

    The body does not depend on the iteration number: everything that varies per step goes
    through the device step counter, and the degeneration of kdpost on the last step is carried
    by ``last``. For the launch order of the phases see the module docstring.
    """
    _launch_h_phase(P)
    _launch_e_phase(P)
    _launch_post_e(P, last)
    # At this point we hold E^{n+1} (t=(n+1)dt) and H^{n+1/2} (t=(n+0.5)dt)
    _sample_time(P)
    _launch_dft(P)
    # Increment the device step counter: this must come after **every** kernel in the body that
    # uses *step. (The first version put it after sampling but before the DFT, so ktab/kbox/kdft
    # read the window of n+1; being one step off is ≈1e-2, and the bitwise gate of the three
    # sentinels caught it on the spot.)
    P.kadv((1,), (1,), (P.stepdev,))


def _time_loop(P: StepPlan, S: SimpleNamespace, use_shutoff: bool) -> tuple[int, bool]:
    """Main time loop: graph capture, per-step launching, time-domain flux sampling,
    early-termination checks, tail-batch flush. Returns ``(steps_run, triggered)``."""
    # ---- CUDA Graph ----
    # Captured once, then one graph.launch per step (10-25 launches collapse into 1). Capturing
    # does not execute, so *step stays 0. The graph is captured on a non-default stream and
    # launched on the current one, which orders it naturally against the default-stream work of
    # shutoff and readout. Numerically bitwise identical to launching kernel by kernel (locked
    # by the criteria).
    n_total = P.n_total
    use_graph = _graph_allowed(P.ten, P.mix, P.ftmons, n_total)
    graph = None
    if use_graph:
        _gs = cp.cuda.Stream(non_blocking=True)
        with _gs:
            _gs.begin_capture()
            _step_body(P, False)
            graph = _gs.end_capture()

    # pre(0) of the fused dispersion: it used to sit inside the body of iteration 0 (E≡0, so it
    # only writes zeros); hoisting it here keeps the body independent of n, bitwise harmless
    if P.disp is not None and (P.disp_fused or P.disp_fused_dec):
        P.kdpre(P.disp["launch"][0], P.disp["launch"][1], P.args_dpre)

    _st_timer = knobs.env("STEP_TIMER") == "1"
    if _st_timer:
        cp.cuda.runtime.deviceSynchronize()
        _st_ev0, _st_ev1 = cp.cuda.Event(), cp.cuda.Event()
        _st_ev0.record()
    triggered = False
    for n in range(n_total):
        if graph is not None and n + 1 < n_total:
            graph.launch()
        else:
            _step_body(P, n + 1 >= n_total)
        _sample_flux_time(n + 1, P.ftmons, P.e_arr, P.h_arr, P.nz)
        if use_shutoff and _early_stop(S, n + 1, P.e_arr, P.h_arr):
            triggered = True
            break

    _flush_tail(n + 1, P.fmons, P.mons, P.dftab, P.stepdev, P.kflush)

    cp.cuda.Stream.null.synchronize()
    if _st_timer:
        _st_ev1.record()
        _st_ev1.synchronize()
        _st_ms = cp.cuda.get_elapsed_time(_st_ev0, _st_ev1) / max(n_total, 1)
        print(f"STEPTIME|{_st_ms:.6f}|{n_total}", flush=True)
    return n + 1, triggered


# ---------------------------------------------------------------------------
# Bloch (k≠0): the complex-field path.
# ---------------------------------------------------------------------------


def _refuse_unsupported_bloch(sc: Scene) -> None:
    """The complex-field path only covers the Bandstructure-shaped subset; anything else fails
    closed.

    Supported: electric and magnetic point dipoles, FieldTimeMonitor, CPML (a z-PML, say),
    periodic / Bloch axes, symmetry planes (the solver is unaware of them; the mirror images are
    already folded into the Scene's dipole table). The rest of the physics has no complex
    variant, and running it silently would drop that physics without raising, so it must be
    refused.
    """
    bad = []
    if sc.sources:
        bad.append("plane wave sources")
    if sc.tfsf_sources:
        bad.append("TFSF sources")
    # Mode sources (oblique plane wave / beam carriers) are supported
    if sc.any_dispersion or sc.any_dispersion_mix:
        bad.append("dispersive media")
    if sc.any_modulation:
        bad.append("time-varying media")
    if sc.any_tensor:
        bad.append("tensor media")
    # Absorber is supported (solver._run_bloch runs absorb_batch once for re and once for im)
    if sc.any_loss:
        bad.append("conductivity")
    if sc.any_pec:
        bad.append("PEC cells")
    # FluxMonitor with the complex colocated flux DFT is supported
    # FieldMonitor frequency-domain accumulation is supported
    if sc.flux_time_monitors:
        bad.append("FluxTimeMonitor")
    if bad:
        raise NotImplementedError(
            "the Bloch complex-field path does not support: " + ", ".join(bad)
            + ". The supported subset is point dipoles (electric/magnetic) + FieldTimeMonitor "
            "+ CPML / periodic boundaries")
    if not sc.dipoles and not sc.mode_sources:
        raise NotImplementedError(
            "the Bloch complex-field path needs a point dipole source or an (oblique) mode "
            "source, and this scene has neither")


def _run_bloch(
    sc: Scene,
    num_steps: int | None = None,
    use_shutoff: bool = True,
    kernels: Kernels | None = None,
    verbose: bool = True,
    return_fields: bool = False,
    policy: shutoff_mod.Policy | None = None,
) -> Result:
    """Complex-field FDTD: Bloch boundaries with k≠0 (kernels/bloch.cu).

    It corresponds line by line to the real path of :func:`run` (same leapfrog timing, same
    two-launch time average of H, same shutoff convention), with only three differences:

    - the fields and psi are split into re/im pairs of float32 (twice the device memory);
    - the ``mnx``/``mpv`` masks turn complex through ``Axis.index_tables(bloch_k)``: the
      wrap-around phase exp(±i·2π·k) shows up here and nowhere else;
    - the dipoles inject the full complex amplitude, and FieldTime samples are complex64.

    At k=0 every phase factor is 1+0i and the injected imaginary part corresponds to an
    independent real solution, so "complex solution = real solution of the Re injection +
    i·(real solution of the Im injection)" holds component by component (tests/test_bloch.py
    locks this down).
    """
    _refuse_unsupported_bloch(sc)
    k = kernels or Kernels()
    n_total = int(num_steps if num_steps is not None else sc.num_time_steps)
    P = StepPlan(sc=sc, k=k, n_total=n_total, verbose=verbose)

    P.update(_alloc_fields_c(sc))
    P.tabs = [a.index_tables(bloch_k=kv)
              for a, kv in zip(sc.grid.axes, sc.bloch_k)]
    # The complex masks are split into re/im (complex_masks); every other table shares the code
    # of the real path
    P.update(_axis_device_tables(
        sc, P.tabs,
        ([_pml_axis_arrays(sc, a, "E") for a in range(3)],
         [_pml_axis_arrays(sc, a, "H") for a in range(3)]),
        P.dims, complex_masks=True))
    P.update(_bloch_e_coeffs(sc, P.dims))                                 # ch ca cb
    P.update(_dipole_tables(sc, n_total, P.dims, complex_amp=True))        # dip dip_h g_dip g_dip_h

    # ---- FieldTimeMonitor: entries as in _monitor_buffers, with two buffers, re and im ----
    P.tmons = [_tmon_entry(m, complex_buf=True) for m in sc.field_time_monitors]
    # ---- FieldMonitor frequency-domain accumulation (complex-field phasors, e.g. the
    # xy/xz/yz_freq monitors at oblique incidence) ----
    _te = (np.arange(n_total, dtype=np.float64) + 1.0) * sc.dt
    _th = (np.arange(n_total, dtype=np.float64) + 0.5) * sc.dt
    P.fmons_c = _bloch_field_monitors(sc, (_te, _th))
    # Absorber (MaxwellStressTensor: Bloch + Absorber): the decay is a real per-cell
    # multiplication of the field, which for a complex field just means multiplying re and im
    # separately, so the absorb_batch tables and kernel of the real path are reused as they are.
    _abs_slabs = _absorber_slabs(sc)
    P.absb = _absorber_batch(_abs_slabs, P.dims, verbose) if _abs_slabs else None
    _bind_kernels_c(P)
    # ---- FluxMonitor complex colocated flux DFT (oblique T/R, the basis of
    # DiffractionMonitor) ----
    P.mons_c = _bloch_flux_monitors(sc, (_te, _th))

    _v = _dual_volumes_eps(sc, P.dims, gpu=False)
    # This path has no plane monitors, so the flux criterion is given an empty monitor list
    # (same value as the real path when mons is empty)
    S = _shutoff_setup(P, policy, (), _v["vols"], _v["eps_dev"])
    del _v

    # ---- Single-plane TF/SF mode source (oblique plane wave / beam, the ModeSource of Bloch
    # scenes) ----
    # Reuses the _mode_source tables of the real path; the injection kernel is swapped for the
    # complex version (inject_mode_*_c), which writes the full complex product profile·amp into
    # the re and im halves. The padding is the same as on the main path.
    P.update(_mode_launch_tables(sc, P.dims, broadband_ok=False))       # mode_srcs ms_launch

    P.g3, P.b3 = grid_3d(P.nx, P.ny, P.nz)
    _plan_update_args_c(P)

    triggered = False
    t0 = time.time()

    _sample_time_c(P, 0)

    for n in range(n_total):
        _step_body_c(P, n)
        # Bloch fields are complex, so every component is passed as (real, imaginary)
        if use_shutoff and _early_stop(S, n + 1, P.e3c, P.h3c):
            triggered = True
            break

    cp.cuda.Stream.null.synchronize()
    wall = time.time() - t0
    steps_run = n + 1

    tsamp, tfill, fields, fphas_c, ph_c = _readout_c(P, return_fields)
    _warn_nonfinite(ph_c, fphas_c)
    return _result(S, steps_run, triggered, wall, fields,
                   phasors=ph_c, field_phasors=fphas_c, time_samples=tsamp,
                   time_slots=tfill)


def _bind_kernels_c(P: StepPlan) -> None:
    """Kernel handles of the complex Bloch path (kernels/bloch.cu); only the batched absorber is
    wired up."""
    k = P.k
    P.kh = k["update_h_c"]
    P.ke = k["update_e_c"]
    P.kabsb = k["absorb_batch"] if P.absb is not None else None
    P.kdip = k["inject_dipole_c"]
    P.kdip_h = k["inject_dipole_h_c"]
    P.kmh_c = k["inject_mode_h_c"]
    P.kme_c = k["inject_mode_e_c"]
    P.kdft_c = k["accumulate_dft_box_c"]
    P.kflux_c = k["accumulate_dft_c"]
    P.ktime = k["sample_time_box_c"]


def _plan_update_args_c(P: StepPlan) -> None:
    """Static arguments of update_h_c / update_e_c: the twelve fields (``h6`` / ``e6``), psi,
    the axis tables, and the ``mnx`` / ``mpv`` masks expanded per axis as (re, im)."""
    mnx_ri, mpv_ri = P.mnx, P.mpv
    P.args_h = (
        *P.h6,
        *P.e6,
        *P.psi_h,
        *P.idl_p, *P.pml_h[0], *P.pml_h[1], *P.pml_h[2],
        *P.nxt,
        mnx_ri[0][0], mnx_ri[0][1], mnx_ri[1][0], mnx_ri[1][1],
        mnx_ri[2][0], mnx_ri[2][1],
        P.ch, *P.n32,
    )
    P.args_e = (
        *P.e6,
        *P.h6,
        *P.psi_e,
        *P.ca, *P.cb,
        *P.idl_d, *P.pml_e[0], *P.pml_e[1], *P.pml_e[2],
        *P.prv,
        mpv_ri[0][0], mpv_ri[0][1], mpv_ri[1][0], mpv_ri[1][1],
        mpv_ri[2][0], mpv_ri[2][1],
        *P.pec, *P.n32,
    )


def _launch_absorb_c(P: StepPlan, fr, fi, which: int) -> None:
    """Absorber decay for Bloch: multiply re and im separately. Only the batched path is wired
    up (the 1D premultiplied profile and slab-by-slab handles are passed as None); when absb is
    None (fewer than two slabs) nothing decays, exactly as before the change."""
    for f3 in (fr, fi):
        _launch_absorb((None, P.kabsb, None), P.absb, [], f3, which, P.n32)


def _sample_time_c(P: StepPlan, m_step: int) -> None:
    """Two launches for the time average of H, same semantics as :func:`_sample_time` on the
    real path; this path uses no graph, the slot is computed on the host, and ``filled`` advances
    with the second launch."""
    for tm in P.tmons:
        for phase, m_ref in ((0, m_step), (1, m_step - 1)):
            off = m_ref - tm["beg"]
            if off < 0 or m_ref >= tm["end"] or off % tm["iv"]:
                continue
            slot = off // tm["iv"]
            if slot >= tm["slots"]:
                continue
            g, b = tm["launch"]
            i0, j0, k0 = tm["origin"]
            P.ktime(g, b, (tm["buf_re"], tm["buf_im"],
                           *P.e6,
                           *P.h6,
                           tm["comps"], tm["w_b"] if phase else tm["w_a"],
                           tm["nc"], np.int32(slot), np.int32(tm["slots"]),
                           np.int32(phase), i0, j0, k0,
                           *tm["box"], *P.n32))
            if phase:
                tm["filled"] = slot + 1


def _step_body_c(P: StepPlan, n: int) -> None:
    """Every kernel launch of one time step on the complex Bloch path (no graph; the step number
    ``n`` is passed directly as an argument): update_h_c -> magnetic dipoles / mode source H
    side -> H-family decay -> update_e_c -> electric dipoles / mode source E side -> E-family
    decay -> time-domain sampling -> complex DFT of the field and flux monitors. The leapfrog
    timing corresponds line by line to :func:`_step_body`."""
    P.kh(P.g3, P.b3, P.args_h)
    if P.dip_h is not None:
        dh = P.dip_h
        P.kdip_h(P.g_dip_h[0], P.g_dip_h[1],
                 (*P.h6, P.ch,
                  dh["flat"], dh["comp"], dh["coef"], dh["src_of"],
                  dh["amp_re"], dh["amp_im"], np.int32(n),
                  dh["nsteps1"], dh["n"]))
    if P.mode_srcs is not None:
        _launch_mode_inject(P.kmh_c, P.h3c, P.ms_launch, P.mode_srcs,
                            side="h", step=np.int32(n), n32=P.n32)
    if P.kabsb is not None:
        _launch_absorb_c(P, P.hr3, P.hi3, 0)
    P.ke(P.g3, P.b3, P.args_e)
    if P.dip is not None:
        d = P.dip
        P.kdip(P.g_dip[0], P.g_dip[1],
               (*P.e6, P.cb[0], P.cb[1], P.cb[2],
                d["flat"], d["comp"], d["coef"], d["src_of"],
                d["amp_re"], d["amp_im"], np.int32(n),
                d["nsteps1"], d["n"]))
    if P.mode_srcs is not None:
        _launch_mode_inject(P.kme_c, P.e3c, P.ms_launch, P.mode_srcs,
                            side="e", step=np.int32(n), n32=P.n32, cb=P.cb)

    if P.kabsb is not None:
        _launch_absorb_c(P, P.er3, P.ei3, 1)
    _sample_time_c(P, n + 1)

    for fm in P.fmons_c:
        P.kdft_c(fm["launch"][0], fm["launch"][1],
                 (fm["re"], fm["im"], *P.e6,
                  *P.h6, fm["ph"][n], np.int32(n),
                  *fm["origin"], *fm["box"], fm["nf"],
                  *P.n32))
    for fx in P.mons_c:
        P.kflux_c(fx["launch"][0], fx["launch"][1],
                  (fx["re"], fx["im"], *P.e6,
                   *P.h6, fx["ph"][n], np.int32(n),
                   fx["km"], fx["axis"], fx["nf"], np.int32(P.nz),
                   *P.n32))


def _readout_c(P: StepPlan, return_fields: bool) -> tuple:
    """Readout of the Bloch path: the time-domain samples are assembled into complex64 (only the
    slots that were filled), and the final fields and the phasors of both monitor kinds are
    assembled into complex arrays. Returns ``(tsamp, tfill, fields, fphas_c, ph_c)``."""
    tsamp = {}
    for tm in P.tmons:
        m_ = tm["filled"]
        tsamp[tm["name"]] = (cp.asnumpy(tm["buf_re"][:, :m_])
                             + 1j * cp.asnumpy(tm["buf_im"][:, :m_])
                             ).astype(np.complex64)
    tfill = {tm["name"]: tm["filled"] for tm in P.tmons}
    fields = {}
    if return_fields:
        for name, (ar, ai) in zip(COMPONENTS, (*P.e3c, *P.h3c)):
            fields[name] = cp.asnumpy(ar[..., :P.nz]) + 1j * cp.asnumpy(ai[..., :P.nz])

    fphas_c = {fm["name"]: (cp.asnumpy(fm["re"]) + 1j * cp.asnumpy(fm["im"]))
               for fm in P.fmons_c}
    ph_c = {fx["name"]: Phasors(name=fx["name"],
                                data=cp.asnumpy(fx["re"]) + 1j * cp.asnumpy(fx["im"]),
                                freqs=cp.asnumpy(fx["freqs"])) for fx in P.mons_c}
    return tsamp, tfill, fields, fphas_c, ph_c
