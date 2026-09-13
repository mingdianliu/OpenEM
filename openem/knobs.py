# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Registry of every ``OPENEM_*`` environment knob (solver side plus the scene/ and nb/ side):
default value, category and purpose.

Three categories:

- ``perf``   performance knobs. The default path is untouched; turning one off only costs speed or
  device memory, the result stays bitwise identical (each one has its own guard test or a
  measurement on record).
- ``policy`` convention knobs. Changing one changes the result (flux integration convention, early
  shutoff guard, colocation snapping), so they are only used in controlled comparisons.
- ``debug``  diagnostics and timing printouts, result unchanged.

How to read it: ``knobs.env("HE_FUSE")`` returns a string and each call site keeps the comparison
it always used (``!= "0"`` means "on unless 0", ``== "1"`` means "on only if 1"). The value is
read **on every call**, not once at import time: tests/test_he_fuse.py, test_ade_fuse.py and
test_flux_*.py all setenv after the import and then run the solver, and reading once would make
them useless; in normal use the two are the same value.

Measurement back doors that have been removed (zero references from tests, notebooks or docs):
``OPENEM_MODE_SIGNS`` / ``OPENEM_TFSF_SIGNS`` (temporary override of the TF/SF sign tables) and
``OPENEM_LEAN_DEBUG`` (printout of the P25b fallback reason).
The knobs on the scene/ and nb/ side (DISP_STAIRCASE, MODE_CACHE, PML_*, SHAPEGRAD* and so on)
have also been registered here since 2026-09-10; those call sites likewise read through
``knobs.env`` every time, with their default values and type conversions unchanged one by one.
The three registered with a default of ``None`` (EPS_CACHE, MODE_CACHE_DIR, HALF_CELL_SIGN) are
the ones where the call site computes the default itself, or branches on set versus unset, so
``env`` returns ``None`` for them when they are unset (tests/test_knobs_scene.py). Knobs of the
notebook integration layer and of external probes are not in this table.
"""

from __future__ import annotations

import os

#: name -> (default, category, purpose)
KNOBS: dict[str, tuple[str | None, str, str]] = {
    # ---- absorber ----
    "ABS_1D": ("1", "perf", ("P50: when same-axis slabs do not overlap, premultiply six 1D decay "
                             "profiles and the slab loop disappears")),
    "ABS_BATCH": ("1", "perf", "P36: several absorber slabs merged into a single launch"),
    "ABS_E_FOLD": ("1", "perf", ("P44: with full dispersion coverage, fold the E-side absorber decay "
                                 "into dispersion_step_dec")),
    "ABSDISP_GPU": ("1", "perf", "P54: compute the decay factors of in-layer dispersion terms on the GPU"),
    "ABSDISP_GPU_MIN": ("15000000", "perf", "P54 size threshold (entry count); below it the host does it"),
    "H_FUSE": ("1", "perf", "P37: update_h and the H-family absorber decay fused into one kernel"),
    # ---- dispersion ----
    "ADE_FUSE": ("0", "perf", ("P43: pole update inlined into update_e (measured slower, off by "
                               "default, guarded by test_ade_fuse)")),
    "DISP_P2": ("1", "perf", "P48: the 2-pole bucket takes the float2 fast path"),
    "DISP_P2_GPU": ("1", "perf", "P53: the P48 bucket permutation is done on the GPU"),
    "DISP_P2_GPU_MIN": ("15000000", "perf", "P53 size threshold (entry count)"),
    "CACB_LUT": ("1", "perf", "P18: lossless deduplication of update_e's (ca, cb) into a lookup table"),
    # ---- interior / boundary shell split ----
    "LEAN": ("1", "perf", "P25b: lean interior kernel plus boundary shell"),
    "LEAN_MINZ": ("32", "perf", "P25b: lower bound on the z extent of the interior"),
    "LEAN_MINFRAC": ("0.70", "perf", ("P25b: lower bound on the interior fraction "
                                      "(calibrated by measurement)")),
    # ---- main update kernels ----
    "HE_FUSE": ("0", "perf", ("P41: H step and E step merged into one kernel (measured 20% slower, "
                              "off by default, guarded by test_he_fuse)")),
    "HE_NI": ("8", "perf", "P41: how many layers along i per block"),
    "BLOCK3": ("64x2", "perf", "thread block shape of the bulk kernel"),
    "PSI_SKIP_FRAC": ("0.5", "perf", "P38: identity-region fraction needed before the psi-skip fast path"),
    "NO_GRAPH": ("", "perf", ("set 1 to not capture the step body in a CUDA Graph (launch kernel by "
                              "kernel, for troubleshooting)")),
    "INIT_GPU": ("1", "perf", ("P28: elementwise math of the setup phase moved onto the GPU; set 0 "
                               "for the numpy fallback")),
    # ---- DFT monitors ----
    "DFT_BATCH": ("1", "perf", "P17: snapshot T steps, then flush them in one batch"),
    "DFT_BATCH_MIN": ("8192", "perf", "P17: fewer than this many points x frequencies, no batching"),
    "DFT_RING_MB": ("2048", "perf", "P17: cap on the snapshot ring buffer (MB)"),
    "DFT_SM2": ("1", "perf", "P49: shared-memory tiled flush kernel"),
    "DFT_SM_NF": ("16", "perf", ("P20b/c: lower bound on the frequency count for the "
                                 "frequency-group version")),
    "DFT_SM_PTS": ("65536", "perf", ("P20b/c: lower bound on the point count for the "
                                     "frequency-group version")),
    "N2F_TANGENTIAL": ("1", "perf", ("P33: n2f surface monitors accumulate only the four "
                                     "tangential components")),
    "TMON_BATCH": ("1", "perf", "P27: batched sampling of time-domain monitors"),
    # ---- TFSF 1D auxiliary grid ----
    "TFSF_GPU": ("1", "perf", "P52: the 1D auxiliary grid runs on the GPU"),
    "TFSF_GPU_MIN": ("20000", "perf", "P52 step-count threshold"),
    "TFSF_WINDOW": ("1", "perf", "P51: the 1D update only covers the causal window"),
    # ---- conventions ----
    "FLUX_COLOCATED": ("1", "policy", ("frequency-domain planar flux uses the colocated integration "
                                       "convention (0 falls back to in-place Yee integration, "
                                       "for comparison)")),
    "TMON_GUARD": ("", "policy", ("full = time-domain monitors must be fully recorded before early "
                                  "shutoff (the default only guards against an empty axis, P56)")),
    "FIELDMON_SNAP": ("0", "policy", ("the interpolation target on a field monitor's zero-thickness "
                                      "axis snaps to the E layer (reverted 2026-09-07, default 0)")),
    # ---- diagnostics ----
    "PAD_TEST": ("", "debug", "force non-zero z padding, to isolate the pitch machinery"),
    "SHUTOFF_DEBUG": ("", "debug", ("verdict experiment: print the readings of all three "
                                    "early-shutoff weightings")),
    "STEP_TIMER": ("", "debug", "set 1 to print STEPTIME|milliseconds per step"),
    "PHASES": ("", "debug", "set 1 to print PHASES|setup|stepping|readout, the three phase timings"),
    "PROFILE": ("0", "debug", "event budget of the kernel-level timer device.Profiler (0 = off)"),
    # ---- scene/: epsilon and media ----
    "DISP_STAIRCASE": ("0", "policy", ("weights of cells on a slanted dispersive interface are "
                                       "overridden with staircasing (scene/dispersion; first rung of "
                                       "the divergence fallback ladder, set and cleared by nb.backend)")),
    "EPS_CACHE": (None, "perf", ("disk cache directory for the subpixel ε (scene/media); unset = "
                                 "cache/epsilon next to the repository, empty string turns it off")),
    "EPS_MEM_GB": ("48", "perf", ("memory budget of the per-frequency ε sampling thread pool (GB): one "
                                  "complex128 full-grid array per worker plus as much again in "
                                  "headroom (scene/weights)")),
    "EPS_WORKERS": ("6", "perf", ("concurrency cap of the per-frequency ε sampling thread pool (lower "
                                  "it when the gencoeffs subprocess gets SIGKILLed)")),
    # ---- scene/: PML and poles ----
    "PML_DISP_ALPHA": ("0.05", "policy", ("CFS α injected on faces that have pole cells inside the "
                                          "ψ-PML band while the user set α_max=0 "
                                          "(= boundaries.PML_DISP_ALPHA), 0 = old behavior")),
    "PML_POLE_DAMP": ("0.3", "policy", ("depth-graded damping of poles inside the ψ-PML band "
                                        "(= boundaries.PML_POLE_DAMP)")),
    "LOSSLESS_DAMPING_REL": ("1e-6", "policy", ("relative artificial loss injected into lossless poles "
                                                "(= poles.LOSSLESS_DAMPING_REL), 0 = old behavior")),
    "LOSSLESS_LOWFREQ_DAMPING_REL": ("0.0", "policy", ("extra damping for low-frequency lossless poles "
                                                       "(= poles.LOSSLESS_LOWFREQ_DAMPING_REL), "
                                                       "default 0 = off")),
    # ---- scene/: mode sources and mode planes ----
    "MODE_CACHE": ("1", "perf", "disk cache for the mode basis (modes.disk_cached), 0 turns it off"),
    "MODE_CACHE_DIR": (None, "perf", ("cache directory for the mode basis; unset = the default "
                                      "directory the caller gives, failing that cache/modes next to "
                                      "the repository")),
    "HALF_CELL_SIGN": (None, "policy", ("measurement back door: override the sign of the mode source "
                                        "half-cell phase (modes.HALF_CELL_SIGN = -1); if set it prints "
                                        "a warning, the physics is not trustworthy")),
    "MODESRC_DISC": ("1", "policy", ("the mode source injects the fully discrete (Yee plus leapfrog) "
                                     "eigenmode; 0 falls back to the continuous Maxwell mode with "
                                     "analytic β")),
    "MODESRC_PAD": ("0", "policy", ("µm by which the mode solve plane is grown outward (this much on "
                                    "each of the two tangential axes, clipped to the simulation "
                                    "domain), 0 = old behavior")),
    "COLOC_SNAP": ("1", "policy", ("the normal coordinate of a mode plane's colocation points snaps to "
                                   "the largest primal boundary not above the center "
                                   "(projection.colocation_points)")),
    # ---- nb/: shapegrad geometry-gradient patch ----
    "SHAPEGRAD": ("mirror", "policy", ("geometry gradient patch mode: mirror (gradient of a "
                                       "mirrored-side geometry is recorded as 0) / eps (geometry "
                                       "gradients go through an ε-map perturbation instead) / 0 off")),
    "SHAPEGRAD_H": ("", "policy", ("ε-map perturbation step (µm); empty = smallest cell size * "
                                   "shapegrad.H_FRACTION")),
    "SHAPEGRAD_DEBUG": ("", "debug", ("set 1 to print the per-component contribution of "
                                      "every shapegrad task")),
}


def env(name: str) -> str | None:
    """Current value of ``OPENEM_<name>``; if unset, the default from the table (the few registered
    as ``None`` return ``None`` when unset, and their default is computed at the call site)."""
    return os.environ.get("OPENEM_" + name, KNOBS[name][0])
