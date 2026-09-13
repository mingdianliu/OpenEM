# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Saving and loading a ``Scene`` (npz). **Does not import tidy3d.**

A GPU job container cannot always have tidy3d installed, so the cut is made along the module
boundary that already exists:

- **where tidy3d is available**: ``td.Simulation`` -> ``Scene`` -> ``scene.npz``
- **in the job container** (numpy and cupy only): ``scene.npz`` -> ``Scene`` -> solve

Saving and loading come in pairs, one pair per entity: ``_save_xxx(d, sc)`` writes keys into
the dict ``d``, and ``_load_xxx(z, ...)`` reads the same keys back out of the npz ``z``.
Writing them as pairs makes a mismatched key name obvious at a glance.
"""

from __future__ import annotations

import os
import pathlib

import numpy as np

from openem import cpml, waveform
from openem.grid import Axis, Grid
from openem.model import (
    AbsorberSlab, BoxFluxMonitor, Dispersion, DispersionMix, FieldMonitor,
    FieldTimeMonitor, FluxMonitor, FluxTimeMonitor, ModeMonitorSpec, ModeSource,
    PermittivityMonitor, PlaneWaveSource, PointDipole, ProjectionMonitor, Scene,
    TensorEps, TFSFSource, TimeModulation)

#: From 2 on, conductivity; from 3 on, the dispersion CSR, absorber slabs and mode sources;
#: from 4 on, TFSF sources, the transverse range of a FluxMonitor, and FieldMonitor /
#: PermittivityMonitor (the monitors of the TFSF cases are all of those two kinds); from 5 on,
#: FieldTimeMonitor and FluxTimeMonitor; from 6 on, the PEC masks and the projection monitors
#: (the solver does not consume projection monitors, but dropping them would silently leave a
#: hole in the write-back and the post-processing); from 7 on, time-varying media
#: (TimeModulation); from 8 on, full tensor medium cells (TensorEps) and box flux monitors;
#: from 9 on, bloch_k and the imaginary part of a dipole's complex amplitude table (the
#: complex field path has to inject the whole exp(iφ)·g(t)); from 10 on, the oblique-incidence
#: geometry of TFSF (k̂/ê/dl₁/proj0/number of table columns, not written at normal incidence).
#: Reading accepts 1..10, writing is always the newest version. The bloch_k of an old npz
#: reads as all zeros: those could only ever be k=0 scenes (k≠0 already failed closed on the
#: extraction side in v8 and earlier).
FORMAT_VERSION = 10

_PML_FIELDS = ("a_E", "b_E", "inv_kappa_E", "a_H", "b_H", "inv_kappa_H")
_DMIX_FIELDS = ("comp", "cell", "p_ofs", "pa", "pb", "q_ofs", "qa", "qb",
                "beta", "eps_inf", "zeta_inf")
_MSRC_BB_FIELDS = ("bb_ey", "bb_ez", "bb_hy", "bb_hz", "bb_amp_e", "bb_amp_h")


# ---------------------------------------------------------------------------
# Small helpers on the read side
# ---------------------------------------------------------------------------

def _npz_opt(z, key):
    """Optional arrays in scene.npz: a missing key, or one stored as an object, both count as
    None (savez stores None as a 0-d object array, which cannot be read back under
    allow_pickle=False; the Medium2D graphene of GrapheneMetamaterial only has in-plane σ, so
    its sigma_ez is None)."""
    if key not in z.files:
        return None
    try:
        a = z[key]
    except ValueError:
        return None
    return None if a.dtype == object or a.shape == () and a.dtype.kind not in "fc" else a


def _n(z, key: str) -> int:
    """A count key (``n_absorbers`` and the like); 0 when an old npz does not have it."""
    return int(z[key]) if key in z else 0


def _names(z, key: str):
    """A name table (``tmon_names`` and the like); an empty list when an old npz lacks it."""
    return z[key] if key in z else []


def _load_waveform(z, prefix: str, dt: float) -> waveform.Waveform:
    """``<prefix>_amp_int`` / ``<prefix>_amp_half`` -> a real waveform (plane wave, TFSF); the
    complex table is rebuilt as purely real (both source types already fail closed on the
    complex path)."""
    amp_int = z[f"{prefix}_amp_int"]
    return waveform.Waveform(
        amp_int=amp_int,
        amp_half=z[f"{prefix}_amp_half"],
        amp_int_complex=amp_int.astype(np.complex128),
        dt=dt,
    )


def _tfsf_obl(z, i: int) -> dict:
    """The five oblique-incidence fields (from v10 on). A normal-incidence npz has no such
    entry, and an empty dict is returned."""
    key = f"tfsf{i}_obl"
    if key not in z:
        return {}
    o = np.asarray(z[key], dtype=np.float64)
    return {"k_hat": tuple(float(x) for x in o[0:3]),
            "e_hat": tuple(float(x) for x in o[3:6]),
            "dl1_step": float(o[6]), "proj0": float(o[7]),
            "n_col_e": int(round(float(o[8]))),
            "amp_scale": float(o[9]) if o.size > 9 else 1.0}


def _apod_pack(apod) -> np.ndarray:
    """``(start, end, width)`` -> three floats, with None encoded as nan."""
    return np.array([np.nan if v is None else float(v) for v in apod])


def _apod_unpack(arr) -> tuple:
    return tuple(None if np.isnan(v) else float(v) for v in arr)


# ---------------------------------------------------------------------------
# Grid / medium arrays / PML
# ---------------------------------------------------------------------------

def _save_grid(d: dict, sc: Scene) -> None:
    for name, ax in zip("xyz", sc.grid.axes):
        d[f"edges_{name}"] = ax.edges
        if getattr(ax, "flat", None) is not None:      # flat-axis flag in 2D (2026-09-05)
            d[f"flat_{name}"] = np.array(bool(ax.flat))
        d[f"bnd_{name}"] = np.array([ax.boundary_lo, ax.boundary_hi])


def _load_grid(z) -> Grid:
    axes = []
    for name in "xyz":
        lo, hi = (str(v) for v in z[f"bnd_{name}"])
        axes.append(Axis(z[f"edges_{name}"], lo, hi,
                         flat=(bool(z[f"flat_{name}"]) if f"flat_{name}" in z else None)))
    return Grid(*axes)


def _save_arrays(d: dict, sc: Scene) -> None:
    """The three groups of volume arrays: ε, σ and the PEC masks."""
    d["eps_ex"] = sc.eps_ex
    d["eps_ey"] = sc.eps_ey
    d["eps_ez"] = sc.eps_ez
    if sc.any_loss:
        d["sigma_ex"] = sc.sigma_ex
        d["sigma_ey"] = sc.sigma_ey
        d["sigma_ez"] = sc.sigma_ez
    if sc.any_pec:
        d["pec_ex"] = sc.pec_ex
        d["pec_ey"] = sc.pec_ey
        d["pec_ez"] = sc.pec_ez


def _load_arrays(z) -> dict:
    """The inverse of :func:`_save_arrays`, returning the matching ``Scene`` keywords.
    Each of the three pec keys is tested on its own (pec_ey/pec_ez used to be tested with
    ``"pec_ex" in z``, which only held up because all three are always written together)."""
    return {
        "eps_ex": z["eps_ex"],
        "eps_ey": z["eps_ey"],
        "eps_ez": z["eps_ez"],
        "sigma_ex": _npz_opt(z, "sigma_ex"),
        "sigma_ey": _npz_opt(z, "sigma_ey"),
        "sigma_ez": _npz_opt(z, "sigma_ez"),
        "pec_ex": z["pec_ex"] if "pec_ex" in z else None,
        "pec_ey": z["pec_ey"] if "pec_ey" in z else None,
        "pec_ez": z["pec_ez"] if "pec_ez" in z else None,
    }


def _save_pml(d: dict, sc: Scene) -> None:
    d["pml_keys"] = np.array([f"{a}:{s}" for (a, s) in sorted(sc.pml)])
    for (a, s), c in sc.pml.items():
        for f in _PML_FIELDS:
            d[f"pml_{a}_{s}_{f}"] = getattr(c, f)


def _load_pml(z) -> dict:
    pml = {}
    for key in z["pml_keys"]:
        a_s, side = str(key).split(":")
        a = int(a_s)
        pml[(a, side)] = cpml.PMLCoeffs(
            **{f: z[f"pml_{a}_{side}_{f}"] for f in _PML_FIELDS}
        )
    return pml


# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------

def _save_sources(d: dict, sc: Scene) -> None:
    d["n_sources"] = np.array(len(sc.sources))
    for i, src in enumerate(sc.sources):
        d[f"src{i}_meta"] = np.array(
            [src.axis, src.plane_index, src.direction, src.pol_axis, src.num_freqs]
        )
        d[f"src{i}_angle_theta"] = np.array(src.angle_theta)
        d[f"src{i}_ex_inc"] = src.ex_inc
        d[f"src{i}_hy_inc"] = src.hy_inc
        d[f"src{i}_amp_int"] = src.waveform.amp_int
        d[f"src{i}_amp_half"] = src.waveform.amp_half
        d[f"src{i}_center"] = np.array(float(getattr(src, "center", float("nan"))))


def _load_sources(z, dt: float) -> list:
    sources = []
    for i in range(int(z["n_sources"])):
        axis, plane_index, direction, pol_axis, num_freqs = (
            int(v) for v in z[f"src{i}_meta"]
        )
        sources.append(
            PlaneWaveSource(
                axis=axis,
                plane_index=plane_index,
                direction=direction,
                pol_axis=pol_axis,
                waveform=_load_waveform(z, f"src{i}", dt),
                angle_theta=float(z[f"src{i}_angle_theta"]),
                num_freqs=num_freqs,
                ex_inc=z[f"src{i}_ex_inc"],
                hy_inc=z[f"src{i}_hy_inc"],
                center=float(z[f"src{i}_center"]) if f"src{i}_center" in z else float("nan"),
            )
        )
    return sources


def _save_mode_sources(d: dict, sc: Scene) -> None:
    d["n_mode_sources"] = np.array(len(sc.mode_sources))
    for i, m in enumerate(sc.mode_sources):
        d[f"msrc{i}_meta"] = np.array([m.plane_index, m.direction, m.axis])
        d[f"msrc{i}_n_eff"] = np.array([m.n_eff.real, m.n_eff.imag])
        for f in ("ey_inc", "ez_inc", "hy_inc", "hz_inc", "amp_e", "amp_h"):
            d[f"msrc{i}_{f}"] = getattr(m, f)
        # Broadband correction terms: if absent the key is not written and reads back as None
        # (the common K=1 case)
        for f in _MSRC_BB_FIELDS:
            v = getattr(m, f, None)
            if v is not None:
                d[f"msrc{i}_{f}"] = v


def _load_mode_sources(z) -> list:
    mode_sources = []
    for i in range(_n(z, "n_mode_sources")):
        meta = [int(v) for v in z[f"msrc{i}_meta"]]
        pi, dr = meta[0], meta[1]
        ne = z[f"msrc{i}_n_eff"]
        mode_sources.append(ModeSource(
            plane_index=pi, direction=dr, n_eff=complex(ne[0], ne[1]),
            axis=meta[2] if len(meta) > 2 else 0,
            ey_inc=z[f"msrc{i}_ey_inc"], ez_inc=z[f"msrc{i}_ez_inc"],
            hy_inc=z[f"msrc{i}_hy_inc"], hz_inc=z[f"msrc{i}_hz_inc"],
            amp_e=z[f"msrc{i}_amp_e"], amp_h=z[f"msrc{i}_amp_h"],
            **{f: z[f"msrc{i}_{f}"] for f in _MSRC_BB_FIELDS
               if f"msrc{i}_{f}" in z}))
    return mode_sources


def _save_tfsf(d: dict, sc: Scene) -> None:
    d["n_tfsf"] = np.array(len(sc.tfsf_sources))
    for i, t in enumerate(sc.tfsf_sources):
        d[f"tfsf{i}_meta"] = np.array(
            [t.direction, *t.box_lo, *t.box_hi, t.inj, t.col0, t.num_freqs, t.axis,
             sum(1 << q for q in t.open_axes)])
        d[f"tfsf{i}_pol"] = np.array([t.pol_u, t.pol_v], dtype=np.float64)
        for f in ("dl1", "eps1", "ex_inc", "hy_inc"):
            d[f"tfsf{i}_{f}"] = getattr(t, f)
        if getattr(t, "sigma1", None) is not None:
            d[f"tfsf{i}_sigma1"] = t.sigma1
        if getattr(t, "k_hat", None) is not None:      # only at oblique incidence
            d[f"tfsf{i}_obl"] = np.array(
                [*t.k_hat, *t.e_hat, t.dl1_step, t.proj0, t.n_col_e,
                 t.amp_scale],
                dtype=np.float64)
        d[f"tfsf{i}_amp_int"] = t.waveform.amp_int
        d[f"tfsf{i}_amp_half"] = t.waveform.amp_half


def _load_tfsf(z, dt: float) -> list:
    tfsf_sources = []
    for i in range(_n(z, "n_tfsf")):
        m = [int(v) for v in z[f"tfsf{i}_meta"]]
        pol = z[f"tfsf{i}_pol"] if f"tfsf{i}_pol" in z else np.array([1.0, 0.0])
        tfsf_sources.append(TFSFSource(
            direction=m[0], box_lo=tuple(m[1:4]), box_hi=tuple(m[4:7]),
            inj=m[7], col0=m[8], num_freqs=m[9],
            axis=m[10] if len(m) > 10 else 2,
            open_axes=tuple(q for q in range(3) if len(m) > 11 and (m[11] >> q) & 1),
            pol_u=float(pol[0]), pol_v=float(pol[1]),
            dl1=z[f"tfsf{i}_dl1"], eps1=z[f"tfsf{i}_eps1"],
            ex_inc=z[f"tfsf{i}_ex_inc"], hy_inc=z[f"tfsf{i}_hy_inc"],
            waveform=_load_waveform(z, f"tfsf{i}", dt),
            sigma1=(z[f"tfsf{i}_sigma1"] if f"tfsf{i}_sigma1" in z else None),
            **_tfsf_obl(z, i)))
    return tfsf_sources


def _save_dipoles(d: dict, sc: Scene) -> None:
    d["n_dipoles"] = np.array(len(sc.dipoles))
    for i, dd in enumerate(sc.dipoles):
        d[f"dip{i}_component"] = np.array(dd.component)
        d[f"dip{i}_indices"] = dd.indices
        d[f"dip{i}_coef"] = dd.coef
        d[f"dip{i}_amp_int"] = dd.waveform.amp_int
        d[f"dip{i}_amp_half"] = dd.waveform.amp_half
        # Only the imaginary part of the complex table is stored (the real part is
        # amp_int/amp_half, guaranteed by the way sample builds them)
        d[f"dip{i}_amp_int_im"] = dd.waveform.amp_int_complex.imag
        if dd.waveform.amp_half_complex is not None:
            d[f"dip{i}_amp_half_im"] = dd.waveform.amp_half_complex.imag
        d[f"dip{i}_magnetic"] = np.array(int(dd.magnetic))


def _load_dipoles(z, dt: float) -> list:
    dips = []
    for i in range(_n(z, "n_dipoles")):
        amp_int = z[f"dip{i}_amp_int"]
        amp_half = z[f"dip{i}_amp_half"]
        # From v9 on the complex table carries an imaginary part; an old npz has none and is
        # rebuilt as purely real (old files can only be k=0 scenes, and the real path consumes
        # the real part only).
        a_int_c = (amp_int + 1j * z[f"dip{i}_amp_int_im"]
                   if f"dip{i}_amp_int_im" in z else amp_int.astype(np.complex128))
        a_half_c = (amp_half + 1j * z[f"dip{i}_amp_half_im"]
                    if f"dip{i}_amp_half_im" in z else None)
        dips.append(PointDipole(
            component=int(z[f"dip{i}_component"]),
            indices=z[f"dip{i}_indices"],
            coef=z[f"dip{i}_coef"],
            waveform=waveform.Waveform(
                amp_int=amp_int, amp_half=amp_half,
                amp_int_complex=a_int_c, dt=dt,
                amp_half_complex=a_half_c),
            magnetic=bool(int(z[f"dip{i}_magnetic"]))
                if f"dip{i}_magnetic" in z else False))
    return dips


# ---------------------------------------------------------------------------
# Medium tables: dispersion / mix / time-varying / tensor / absorber
# ---------------------------------------------------------------------------

def _save_media(d: dict, sc: Scene) -> None:
    if sc.dispersion_mix is not None and sc.dispersion_mix.n_entry > 0:
        mm = sc.dispersion_mix
        for f in _DMIX_FIELDS:
            d[f"dmix_{f}"] = getattr(mm, f)
    if sc.dispersion is not None and sc.dispersion.n_entry > 0:
        dd = sc.dispersion
        d["disp_comp"] = dd.comp
        d["disp_cell"] = dd.cell
        d["disp_pole_ofs"] = dd.pole_ofs
        d["disp_am1"] = dd.am1
        d["disp_b"] = dd.b
        d["disp_g"] = dd.g

    if sc.any_modulation:
        mo = sc.modulation
        d["mod_freq"] = np.array(mo.freq)
        d["mod_comp"] = mo.comp
        d["mod_cell"] = mo.cell
        d["mod_amp"] = mo.amp
        d["mod_phase"] = mo.phase

    if sc.tensor is not None and sc.tensor.n_entry > 0:
        for f in ("comp", "cell", "eps", "sigma"):
            d[f"tensor_{f}"] = getattr(sc.tensor, f)


def _load_media(z) -> dict:
    """Return the ``dispersion / dispersion_mix / modulation / tensor`` keywords of ``Scene``."""
    dispersion = None
    if "disp_comp" in z:
        dispersion = Dispersion(
            comp=z["disp_comp"], cell=z["disp_cell"], pole_ofs=z["disp_pole_ofs"],
            am1=z["disp_am1"], b=z["disp_b"], g=z["disp_g"])

    modulation = None
    if "mod_comp" in z:
        modulation = TimeModulation(
            freq=float(z["mod_freq"]), comp=z["mod_comp"], cell=z["mod_cell"],
            amp=z["mod_amp"], phase=z["mod_phase"])

    tensor = None
    if "tensor_comp" in z:
        tensor = TensorEps(comp=z["tensor_comp"], cell=z["tensor_cell"],
                           eps=z["tensor_eps"], sigma=z["tensor_sigma"])

    dispersion_mix = None
    if "dmix_comp" in z:
        dispersion_mix = DispersionMix(**{f: z[f"dmix_{f}"] for f in _DMIX_FIELDS})
    return {"dispersion": dispersion, "dispersion_mix": dispersion_mix,
            "modulation": modulation, "tensor": tensor}


def _save_absorbers(d: dict, sc: Scene) -> None:
    d["n_absorbers"] = np.array(len(sc.absorbers))
    for i, sl in enumerate(sc.absorbers):
        d[f"abs{i}_meta"] = np.array([sl.axis, sl.g0])
        d[f"abs{i}_decay_int"] = sl.decay_int
        d[f"abs{i}_decay_half"] = sl.decay_half


def _load_absorbers(z) -> list:
    absorbers = []
    for i in range(_n(z, "n_absorbers")):
        ax_, g0_ = (int(v) for v in z[f"abs{i}_meta"])
        absorbers.append(AbsorberSlab(
            axis=ax_, g0=g0_,
            decay_int=z[f"abs{i}_decay_int"], decay_half=z[f"abs{i}_decay_half"]))
    return absorbers


# ---------------------------------------------------------------------------
# Monitors (eight kinds)
# ---------------------------------------------------------------------------

def _save_monitors(d: dict, sc: Scene) -> None:
    d["mon_names"] = np.array([m.name for m in sc.flux_monitors])
    for i, m in enumerate(sc.flux_monitors):
        d[f"mon{i}_meta"] = np.array([m.axis, m.plane_index, m.normal_dir])
        d[f"mon{i}_freqs"] = m.freqs
        d[f"mon{i}_spectrum"] = m.source_spectrum
        # Without the transverse range the whole plane would be integrated silently, so it has
        # to be stored along with the rest
        if m.transverse is not None:
            d[f"mon{i}_transverse"] = np.array(m.transverse)
        if m.apodization is not None:
            d[f"mon{i}_apod"] = _apod_pack(m.apodization)
        # Tangential physical bounds: needed by the colocated integration convention
        # (2026-09-06); an old npz has no such entry and falls back to the in-place Yee
        # convention
        if getattr(m, "t_bounds", None) is not None:
            d[f"mon{i}_tbounds"] = np.array(m.t_bounds, dtype=np.float64)

    d["fmon_names"] = np.array([m.name for m in sc.field_monitors])
    for i, m in enumerate(sc.field_monitors):
        d[f"fmon{i}_meta"] = np.array([*m.origin, *m.box, *m.interval_space])
        d[f"fmon{i}_freqs"] = m.freqs
        if m.apodization is not None:
            d[f"fmon{i}_apod"] = _apod_pack(m.apodization)

    d["modemon_names"] = np.array([m.name for m in sc.mode_monitors])
    d["modemon_planes"] = np.array([m.plane_name for m in sc.mode_monitors])

    d["pmon_names"] = np.array([m.name for m in sc.permittivity_monitors])
    for i, m in enumerate(sc.permittivity_monitors):
        d[f"pmon{i}_meta"] = np.array([*m.origin, *m.box])
        d[f"pmon{i}_freqs"] = m.freqs

    d["projmon_names"] = np.array([m.name for m in sc.projection_monitors])
    for i, m in enumerate(sc.projection_monitors):
        d[f"projmon{i}_surfaces"] = np.array(m.surface_names)
        d[f"projmon{i}_dirs"] = np.array(m.normal_dirs)

    d["bfmon_names"] = np.array([m.name for m in sc.box_flux_monitors])
    for i, m in enumerate(sc.box_flux_monitors):
        d[f"bfmon{i}_faces"] = np.array(m.face_names)

    d["tmon_names"] = np.array([m.name for m in sc.field_time_monitors])
    for i, m in enumerate(sc.field_time_monitors):
        d[f"tmon{i}_meta"] = np.array([m.step_begin, m.step_end, m.interval,
                                       m.num_slots])
        d[f"tmon{i}_comps"] = np.array(m.comps)
        d[f"tmon{i}_origin"] = np.array(m.origin)
        d[f"tmon{i}_box"] = np.array(m.box)

    d["ftmon_names"] = np.array([m.name for m in sc.flux_time_monitors])
    for i, m in enumerate(sc.flux_time_monitors):
        d[f"ftmon{i}_meta"] = np.array([m.axis, m.plane_index, m.normal_dir,
                                        m.step_begin, m.step_end, m.interval,
                                        m.num_slots])
        d[f"ftmon{i}_frac"] = np.array(m.frac)
        d[f"ftmon{i}_tbounds"] = np.array(m.t_bounds, dtype=np.float64)


def _load_flux_mons(z) -> list:
    """List of FluxMonitors."""
    mons = []
    for i, name in enumerate(z["mon_names"]):
        axis, plane_index, normal_dir = (int(v) for v in z[f"mon{i}_meta"])
        tkey = f"mon{i}_transverse"
        mons.append(
            FluxMonitor(
                name=str(name),
                axis=axis,
                plane_index=plane_index,
                normal_dir=normal_dir,
                freqs=z[f"mon{i}_freqs"],
                source_spectrum=z[f"mon{i}_spectrum"],
                transverse=tuple(map(tuple, z[tkey].tolist())) if tkey in z else None,
                apodization=(_apod_unpack(z[f"mon{i}_apod"])
                             if f"mon{i}_apod" in z else None),
                t_bounds=(tuple(map(tuple, z[f"mon{i}_tbounds"].tolist()))
                          if f"mon{i}_tbounds" in z else None),
            )
        )

    return mons


def _load_time_mons(z) -> list:
    """List of FieldTimeMonitors."""
    tmons = []
    for i, name in enumerate(_names(z, "tmon_names")):
        beg, end, iv, slots = (int(v) for v in z[f"tmon{i}_meta"])
        tmons.append(
            FieldTimeMonitor(
                name=str(name),
                comps=tuple(int(v) for v in z[f"tmon{i}_comps"]),
                origin=tuple(int(v) for v in z[f"tmon{i}_origin"]),
                box=tuple(int(v) for v in z[f"tmon{i}_box"]),
                step_begin=beg,
                step_end=end,
                interval=iv,
                num_slots=slots,
            )
        )

    return tmons


def _load_flux_time_mons(z) -> list:
    """List of FluxTimeMonitors."""
    ftmons = []
    for i, name in enumerate(_names(z, "ftmon_names")):
        axis, plane_index, normal_dir, beg, end, iv, slots = (
            int(v) for v in z[f"ftmon{i}_meta"])
        tb = z[f"ftmon{i}_tbounds"]
        ftmons.append(
            FluxTimeMonitor(
                name=str(name),
                axis=axis,
                plane_index=plane_index,
                normal_dir=normal_dir,
                step_begin=beg,
                step_end=end,
                interval=iv,
                num_slots=slots,
                t_bounds=tuple(tuple(float(v) for v in r) for r in tb),
                frac=float(z[f"ftmon{i}_frac"]),
            )
        )

    return ftmons


def _load_field_mons(z) -> list:
    """List of FieldMonitors."""
    fmons = []
    for i, name in enumerate(_names(z, "fmon_names")):
        m = [int(v) for v in z[f"fmon{i}_meta"]]
        fmons.append(FieldMonitor(
            name=str(name), freqs=z[f"fmon{i}_freqs"],
            origin=tuple(m[0:3]), box=tuple(m[3:6]),
            interval_space=tuple(m[6:9]) if len(m) > 6 else (1, 1, 1),
            apodization=(_apod_unpack(z[f"fmon{i}_apod"])
                         if f"fmon{i}_apod" in z else None)))

    return fmons


def _load_mode_mons(z) -> list:
    """List of ModeMonitorSpecs (an old npz has no such section)."""
    if "modemon_names" not in z:
        return []
    return [ModeMonitorSpec(name=str(n), plane_name=str(p))
            for n, p in zip(z["modemon_names"], z["modemon_planes"])]


def _load_perm_mons(z) -> list:
    """List of PermittivityMonitors."""
    pmons = []
    for i, name in enumerate(_names(z, "pmon_names")):
        m = [int(v) for v in z[f"pmon{i}_meta"]]
        pmons.append(PermittivityMonitor(
            name=str(name), freqs=z[f"pmon{i}_freqs"],
            origin=tuple(m[0:3]), box=tuple(m[3:6])))

    return pmons


def _load_proj_mons(z) -> list:
    """List of ProjectionMonitors."""
    projmons = []
    for i, name in enumerate(_names(z, "projmon_names")):
        projmons.append(ProjectionMonitor(
            name=str(name),
            surface_names=tuple(str(s) for s in z[f"projmon{i}_surfaces"]),
            normal_dirs=tuple(str(s) for s in z[f"projmon{i}_dirs"])))

    return projmons


def _load_box_flux_mons(z) -> list:
    """List of BoxFluxMonitors."""
    bfmons = []
    for i, name in enumerate(_names(z, "bfmon_names")):
        bfmons.append(BoxFluxMonitor(
            name=str(name),
            face_names=tuple(str(s) for s in z[f"bfmon{i}_faces"])))
    return bfmons


def _load_monitors(z) -> dict:
    """Return the eight monitor-list keywords of ``Scene`` (one ``_load_*`` per kind, matching
    the eight blocks of :func:`_save_monitors` one for one)."""
    return {
        "flux_monitors": _load_flux_mons(z),
        "field_monitors": _load_field_mons(z),
        "mode_monitors": _load_mode_mons(z),
        "field_time_monitors": _load_time_mons(z),
        "flux_time_monitors": _load_flux_time_mons(z),
        "permittivity_monitors": _load_perm_mons(z),
        "projection_monitors": _load_proj_mons(z),
        "box_flux_monitors": _load_box_flux_mons(z),
    }


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

def save(sc: Scene, path: str | pathlib.Path) -> pathlib.Path:
    """Write a Scene out as a single npz."""
    d: dict[str, np.ndarray] = {
        "format_version": np.array(FORMAT_VERSION),
        "dt": np.array(sc.dt),
        "num_time_steps": np.array(sc.num_time_steps),
        "shutoff": np.array(sc.shutoff),
        "symmetry": np.array(sc.symmetry),
        "bloch_k": np.array(sc.bloch_k, dtype=np.float64),
    }
    _save_arrays(d, sc)
    _save_grid(d, sc)
    _save_pml(d, sc)
    _save_sources(d, sc)
    _save_media(d, sc)
    _save_absorbers(d, sc)
    _save_mode_sources(d, sc)
    _save_tfsf(d, sc)
    _save_monitors(d, sc)
    _save_dipoles(d, sc)

    p = pathlib.Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Atomic write: write a .part file in the same directory and rename it, so a reader never
    # sees half a file
    tmp = p.with_name(p.name + f".{os.getpid()}.part")
    with open(tmp, "wb") as f:
        np.savez_compressed(f, **d)
    os.replace(tmp, p)
    return p


def load(path: str | pathlib.Path) -> Scene:
    """Read a Scene back from an npz. numpy only."""
    z = np.load(path, allow_pickle=False)
    ver = int(z["format_version"])
    if not 1 <= ver <= FORMAT_VERSION:
        raise ValueError(
            f"scene.npz format version {ver}, this code accepts 1..{FORMAT_VERSION}: export again")

    dt = float(z["dt"])
    return Scene(
        grid=_load_grid(z),
        **_load_media(z),
        absorbers=_load_absorbers(z),
        mode_sources=_load_mode_sources(z),
        dt=dt,
        num_time_steps=int(z["num_time_steps"]),
        shutoff=float(z["shutoff"]),
        symmetry=tuple(int(v) for v in z["symmetry"]) if "symmetry" in z else (0, 0, 0),
        bloch_k=(tuple(float(v) for v in z["bloch_k"])
                 if "bloch_k" in z else (0.0, 0.0, 0.0)),
        **_load_arrays(z),
        pml=_load_pml(z),
        sources=_load_sources(z, dt),
        dipoles=_load_dipoles(z, dt),
        tfsf_sources=_load_tfsf(z, dt),
        **_load_monitors(z),
    )
