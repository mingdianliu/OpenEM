# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Tidy3D's native normalization: convert the solver's raw phasors into Tidy3D's convention,
with **no fitting**.

The rule: divide by the naive DFT of "the time table the normalization source actually injects"
(the same convention as the monitor DFT), then divide by the unit difference of 1e6
(V/m -> V/µm). Which injection table is used depends on the source type:

    PlaneWave / TFSF -> ``ex_inc`` (the incident field table; per-case constants such as
                        amplitude, impedance and background index are absorbed
                        **automatically**, with no per-source-type derivation)
    ModeSource       -> ``amp_e``
    PointDipole      -> ``waveform.amp_int``

Evidence (measured): the fitted constant k was worked back out class by class. Dipole cases give
k=1.000 (the analytic constant is the DFT of the injection table); TFSF cases give
k=1e6*√(2η₀) (the incident table carries that constant itself); plane wave cases give a
different k per case (the incident table absorbs it case by case). The solver itself is
untouched, this module is only used in the readout, scoring and export layers.
"""
from __future__ import annotations

import numpy as np

from openem import flux as _flux
from openem.grid import transverse_axes

#: V/m → V/µm
UNIT_E = 1.0e6
ETA0 = 376.730313668
C_0 = 299792458.0


def _src_half_cell_cos(sc, freqs) -> np.ndarray:
    """The plane wave source's "discrete injected flux / continuous power" = cos(k·dz_h/2), per
    frequency (2026-09-05).

    What a Yee grid conserves exactly is the discrete Poynting flux ½Re(E_k·H*_{k+½}); for a
    travelling wave that is smaller than the continuous flux |E|²/(2η) by cos(k·dz/2) (E and H
    are half a cell apart). The two-point colocated flux monitor ½Re[E_k·(H_{k-½}+H_{k+½})*] is
    exactly the average of two neighboring discrete fluxes, so it measures that same conserved
    quantity, whether or not the grid is uniform. Measured on a 1D probe: vacuum transmission at
    10/20/40/60 cells/λ = 0.9511/0.9877/0.9969/0.9986, and over five frequency points, 2D and
    3D, uniform and automatic grids, all of it agrees with cos(k·dz/2) digit for digit.

    Tidy3D normalizes the source to 1 W of **discrete** flux (measured against its solver: the
    flux stays at 1.000 and |E| is larger than the analytic value by 1/√cos, 22.99 at the center
    frequency with 10 cells/λ versus 22.41 analytic and 22.98 predicted), so its T is a ratio of
    two conserved quantities and does not drift with the grid. Our 1D auxiliary grid injects by
    continuous power (E analytic), so the normalization divisor becomes the discrete injected
    flux on the source plane instead: the flux divisor gets *cos and the field divisor *√cos
    (the field is larger than the analytic value by 1/√cos, as on the reference side).
    dz is taken from the cell the H correction face sits in (ks-1 for the + direction, ks for
    the -, the same dl_h as sources.plane_wave); k is the continuous wavenumber at the
    background index of the source plane (the probe shows cos agrees with the continuous k to
    four digits).

    The reverse approach was tried as well, four-point Lagrange colocation in the monitor with
    the source kept on continuous power: a vacuum plane wave does reach 0.9965, but a mode
    source on a straight Si waveguide at 9 cells/λ reads 1.042 (1.000 on the reference).
    Four-point colocation measures the continuous flux, which is **not conserved** on a Yee
    grid, so the injection error shows up in full; under the discrete flux convention the same
    scene gives 1.007. Dropped.
    """
    src = sc.sources[0]
    dl = np.asarray(sc.grid.axes[int(src.axis)].dl, dtype=np.float64)
    kh = int(src.plane_index) - 1 if int(src.direction) > 0 else int(src.plane_index)
    dz = float(dl[min(max(kh, 0), dl.size - 1)])
    k = (2.0 * np.pi * np.atleast_1d(np.asarray(freqs, dtype=np.float64)) / C_0
         * background_index(sc))
    return np.cos(k * dz / 2.0)


def _src_snap_phase(sc, freqs) -> np.ndarray:
    """Snapping phase of the 1D plane wave source plane (2026-09-05): Tidy3D refers the plane
    wave phase to **the source center the user gave**, z_c, while we inject on the snapped grid
    line z_e=edges[plane_index], so the field carries an extra global phase
    e^{ik·dir·(z_c-z_e)} (measured: z_c-z_e=+0.0161 µm predicts +3.75°, and the reference solver
    gives +3.79° point by point; mode sources have the same kind of problem, while beam sources
    evaluate the analytic field profile on the real injection plane and never have it).
    Multiplying the divisor by this factor divides it back out. An old scene.npz has no center
    stored (nan) -> factor 1."""
    src = sc.sources[0]
    zc = float(getattr(src, "center", float("nan")))
    freqs = np.atleast_1d(np.asarray(freqs, dtype=np.float64))
    if not np.isfinite(zc):
        return np.ones(freqs.shape, dtype=np.complex128)
    ze = float(sc.grid.axes[src.axis].edges[src.plane_index])
    k = 2.0 * np.pi * freqs / C_0 * background_index(sc)
    return np.exp(1j * float(src.direction) * k * (zc - ze))


def source_table(sc) -> np.ndarray:
    """Injection time table of the scene's normalization source (Tidy3D's normalize_index=0
    convention).

    TFSF takes priority over sources (in a TFSF scene both can be non-empty).
    """
    if getattr(sc, "tfsf_sources", None):
        return np.asarray(sc.tfsf_sources[0].ex_inc, dtype=np.float64)
    if getattr(sc, "sources", None):
        return np.asarray(sc.sources[0].ex_inc, dtype=np.float64)
    if getattr(sc, "mode_sources", None):
        # The ½ convention of mode injection: the amp_e table is twice the source time series
        # (calibration measured exactly 2.000000)
        return np.asarray(sc.mode_sources[0].amp_e, dtype=np.complex128) / 2.0
    if getattr(sc, "dipoles", None):
        return np.asarray(sc.dipoles[0].waveform.amp_int, dtype=np.float64)
    raise NotImplementedError("the scene has no usable normalization source")


def norm_spectrum(sc, freqs) -> np.ndarray:
    """Naive DFT of the injection table (same convention as the monitor DFT: E is taken at
    t=(n+1)·dt, without a dt factor)."""
    freqs = np.atleast_1d(np.asarray(freqs, dtype=np.float64))
    s = source_table(sc)
    n = np.arange(s.size, dtype=np.float64)
    t = n * sc.dt          # calibration: δ ≡ 0 for every source type
    return np.array([np.sum(s * np.exp(2j * np.pi * f * t)) for f in freqs])


#: √(2η₀): the TFSF incident amplitude convention difference (fitted 27.43-27.45 vs √753.5=27.45)
SQRT_2ETA0 = float(np.sqrt(2.0 * 376.730313668))


def field_norm(sc, freqs) -> np.ndarray:
    """Divisor of the field phasors (the rule per source type is in the module docstring)."""
    D = norm_spectrum(sc, freqs)
    if getattr(sc, "tfsf_sources", None):
        # divisor = D·1e6/k, fitted k=27.43e6=1e6·√(2η₀) => divisor = D/√(2η₀)
        # Phase: the TFSF injection table is shifted by 3.5 steps against the monitor DFT
        # (diagnosed as δs=-3.497 on Near2Far)
        ph = np.exp(2j * np.pi * np.atleast_1d(np.asarray(freqs, float))
                    * 3.5 * sc.dt)
        return D / SQRT_2ETA0 * ph
    if getattr(sc, "mode_sources", None) and not getattr(sc, "sources", None):
        # Unit-power normalization => no 1e6.
        #
        # **The phase closes** (an older comment claimed it did not and applied no correction;
        # that was out of date).
        # The test: ``field_norm * projection.normalization`` is identically 1∠0, and the
        # second factor is Tidy3D's own ``(dt/√2π)/source_time.spectrum``, which says the DFT
        # of our injection table and the reference source spectrum are bit for bit the same
        # thing. Measured on 2026-08-31 over 8 mode source cases: 7 give 1.00000000∠+0.000000
        # and only 90OpticalHybrid gives 1.02922∠0 (2.9% off in amplitude, phase still 0; it
        # is also the one case whose mode amps cannot be computed).
        #
        # This one matters: the adjoint source amplitude of a mode objective is
        # ``1j·(k0/4/η₀)·amps``, the phase of amps goes into the adjoint source unchanged, and
        # the gradient rotates by however much the phase is off (compare the diffraction case,
        # which was off by exactly π/2 and flipped the sign of the whole gradient).
        return D
    if getattr(sc, "sources", None):
        # Plane wave field divisor = D·1e6·√(n_bg·A/(2η₀)).
        #
        # That √A is **not a leftover of fitting, do not drop it**: Tidy3D's PlaneWave is
        # normalized to "unit power through the source cross-section" (lengths in µm), so the
        # field it gives is |E| = √(2η₀/(n·A_µm²)), which is meant to shrink as the simulation
        # gets wider.
        #
        # Three independent pieces of evidence (2026-08-30):
        # 1. The uniform reference solution of MIM has raw/D = 0.999878; we inject at **unit E
        #    amplitude**, the current divisor reads the same incident field as 63.31, and
        #    √(2η₀/(1.45·0.1296)) = 63.32.
        # 2. Two transverse widths in vacuum (1 µm / 2 µm) halve |E| exactly, = 1/√A.
        #    At A = 1 µm² it degenerates to √(2η₀) = 27.45, matching the TFSF branch (unit
        #    intensity).
        # 3. Point-by-point ratio against the reference output of MIM: the median |ratio| is
        #    1.0015 and the **median phase is -0.019 rad**, so amplitude and phase both close
        #    already and no time shift correction like the TFSF one is needed.
        src = sc.sources[0]
        t1, t2 = transverse_axes(src.axis)
        A = (float(np.sum(_flux.axis_dl(sc.grid.axes[t1])[0]))
             * float(np.sum(_flux.axis_dl(sc.grid.axes[t2])[0])))
        # *√cos: the discrete injected flux convention (_src_half_cell_cos); the field is larger
        # than the analytic value by 1/√cos, as on the reference side
        # *snapping phase: Tidy3D refers the phase to the user's source center, we refer it to
        # the snapped grid line (_src_snap_phase)
        return (D * UNIT_E * np.sqrt(background_index(sc) * A / (2.0 * ETA0))
                * np.sqrt(_src_half_cell_cos(sc, freqs)) * _src_snap_phase(sc, freqs))
    return D * UNIT_E


def time_scale(sc, source_time, dt, num_steps) -> float:
    """Conversion factor λ for time-domain quantities: the solver's raw time-domain field is
    E_ours(t) = λ·E_td3d(t) (2026-09-05).

    Derivation: once normalized, the two sides agree in the frequency domain,
    ``E_ours(f)/field_norm(f) = E_td3d(f)/S_real(f)``, where S_real is the naive DFT of the
    **real signal** Re(amp_time) of the reference normalization source (the same DFT convention
    as :func:`norm_spectrum`, t = n·dt). Both sides inject the same waveform, so the
    time-domain ratio is a constant ``λ = field_norm(f0)/S_real(f0)`` (frequency points differ
    only by slow variations under 0.3%, of the √cos kind, so f0 is used).
    Checked against small reference runs of all three source types: plane wave λ=0.0484, mode
    source λ=1, point dipole λ=1e6, measurement within 0.5% of the prediction and no time lag.
    Readout: FieldTime = raw/λ; FluxTime = raw/(λ²·UM²) (E and H each pick up 1/λ and the area
    element goes m²->µm²).
    The FieldTime/1e6 that used to be hard-coded was only right for current sources (the plane
    wave time-domain field of AndersonLocalization came out 2e7 times too small, the
    AnimationTutorial mode source 1e6 times too small), and FluxTime/UM² only for plane wave
    and mode sources (the NanobeamCavity dipole time-domain flux came out 1e12 times too big).
    """
    f0 = float(source_time.freq0)
    t = np.arange(int(num_steps), dtype=np.float64) * float(dt)
    s_real = np.sum(np.real(np.asarray(source_time.amp_time(t))) * np.exp(2j * np.pi * f0 * t))
    if abs(s_real) == 0.0:
        raise ValueError("the real-signal spectrum of the normalization source is zero at f0, "
                         "so the time-domain conversion factor cannot be computed")
    return float(abs(field_norm(sc, [f0])[0] / s_real))


def flux_norm(sc, freqs) -> np.ndarray:
    """Flux divisor: |field_norm|²/1e12 (the two 1e6 factors cancel the 1e-12 of the µm² area)."""
    return np.abs(field_norm(sc, freqs)) ** 2 / 1.0e12


def td_source_spectrum(source_time, freqs, dt, num_steps) -> np.ndarray:
    """The reference solver's own source spectrum ``S(f)``, which has to be **multiplied back
    in** when it ran with ``normalize_index=None``.

    The divisors in this module line up with the reference default of ``normalize_index=0``
    (where the data has already been divided by the source spectrum). A few scenes had that
    normalization switched off (3 cases: CavityFOM, ResonanceFinder, TimeModulationTutorial) and
    store the **raw DFT**; there the readout side has to multiply the source spectrum back in,
    fields by ``S(f)`` and flux by ``|S(f)|²``.

    The factor is derived, not fitted: the calibration says the normalized quantities on both
    sides are equal, ``E_ours/divisor == E_td_raw/S``, hence
    ``E_td_raw = E_ours/divisor · S``. ``S`` calls the reference ``SourceTime.spectrum``
    directly, the very function it normalizes with.

    Measured (no fitting): after conversion, the flux of TimeModulation over 181 significant
    frequency points is within 1.2% of the reference output; the field amplitude ratio of
    CavityFOM is 4.5387e13 against 4.5504e13 for ``1/|S|``, 0.26% apart.

    A step or two of difference in ``num_steps`` does not matter, the tail of the pulse has
    already decayed to zero. ``post_norm`` is 1.0 in every reference case, so there is no
    instance to check against and it is not handled.
    """
    t = np.arange(int(num_steps), dtype=np.float64) * float(dt)
    return np.asarray(source_time.spectrum(
        t, np.atleast_1d(np.asarray(freqs, dtype=np.float64)),
        float(dt))).ravel()


def td_denorm_from_archive(sim_json, dt, freqs):
    """:func:`td_source_spectrum` in terms of the stored JSON; returns ``None`` when the
    reference run did normalize.

    It only consumes the ``simulation`` section of ``JSON_STRING`` inside ``data.hdf5``: neither
    the scorer nor the injector builds a ``td.Simulation``, so assembling the source-time object
    from the JSON on the spot is enough. The step count is worked back out from ``run_time/dt``
    (a step or two does not matter, see :func:`td_source_spectrum`).
    """
    if (sim_json.get("normalize_index") is not None
            or not sim_json.get("sources")):
        return None
    import tidy3d as td                  # only those 3 cases get here, so import it lazily
    d = sim_json["sources"][0]["source_time"]
    st = getattr(td, d["type"]).parse_obj(d)
    n = int(round(float(sim_json["run_time"]) / float(dt)))
    return td_source_spectrum(st, freqs, dt, n)


def background_index(sc) -> float:
    """Background index: the square root of the median eps **on the upstream side of the source
    plane** (3 cells back against the propagation direction), so the sample is taken in the
    medium the source sits in and stays clear of the structure layer (the AllDielectric
    lesson)."""
    src = sc.sources[0]
    e = sc.eps_ex
    k = int(src.plane_index) - 3 * int(src.direction)
    k = max(0, min(e.shape[src.axis] - 1, k))
    sl = [slice(None)] * 3
    sl[src.axis] = k
    return float(np.sqrt(np.median(e[tuple(sl)])))


def monitor_area(sc, m) -> float:
    """Transverse area of a flux monitor (m², including the transverse crop).

    The flat axis in 2D follows the **same** rule as :func:`openem.flux.area_weights`
    (``flux.FLAT_AXIS_M``); otherwise the flux would change while the divisor did not, and the
    convention of the ratio would go wrong with it.
    """
    t1, t2 = transverse_axes(m.axis)
    d1 = _flux.axis_dl(sc.grid.axes[t1])[0]
    d2 = _flux.axis_dl(sc.grid.axes[t2])[0]
    if getattr(m, "transverse", None) is not None:
        (a1, b1), (a2, b2) = m.transverse
        d1 = d1[a1:b1]
        d2 = d2[a2:b2]
    return float(d1.sum() * d2.sum())


def plane_flux_divisor(sc, m, freqs) -> np.ndarray:
    """Flux divisor for plane wave scenes: |D|²·n_bg·A/(2η₀)·cos(k·dz_h/2) (fully analytic, no
    reference run needed; for cos see _src_half_cell_cos)."""
    D = norm_spectrum(sc, freqs)
    return (np.abs(D) ** 2 * background_index(sc)
            * monitor_area(sc, m) / (2.0 * ETA0)
            * _src_half_cell_cos(sc, freqs))       # discrete injected flux convention
