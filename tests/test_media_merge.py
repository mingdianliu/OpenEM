# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Merging collinear medium columns (``weights.column_groups``).

Two criteria:

1. **Parameter-free**: define the same material twice under different names (``media_columns``
   deduplicates by JSON content, so it will not merge them) and the resulting Scene must be
   **bitwise identical** to a scene where they were merged into one material by hand.
2. **The real cause**: the columns of two non-dispersive materials are real multiples of each other
   (the rank deficiency of PlasmonicNanorodArray); after merging, the scene must build, and for a
   slanted interface the material pair must be recovered candidate by candidate from the merged
   group.
"""

from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import scene as scene_mod
from openem.scene import media
from openem.scene.weights import (
    FIT_RTOL, column_groups, decompose, fit_frequencies, media_columns)

#: The monitor band has to be **wide**: 500-1000 nm, over which cSi's ε moves from ~13 to ~18. In a
#: narrow band the harmonic mean is nearly affine in ε(f), least squares fits it anyway, and the
#: mixed cells are never exercised at all.
FREQ0 = 4.5e14
FREQS = (3.0e14, 4.5e14, 6.0e14)
DL = 0.05

#: The interface sits at 0.3 of a cell, which avoids both the primal grid face (0) and Ez's dual
#: face (0.5), so mixed cells are guaranteed and the mix branch is really taken
_ZA = -0.4 + 0.3 * DL               # bottom face of the cSi slab (interface with the material below)
_ZB = 0.4 + 0.3 * DL                # top face of the cSi slab (interface with vacuum)


def _csi():
    return td.material_library["cSi"]["Green2008"]


def _sim(structures: list) -> td.Simulation:
    return td.Simulation(
        size=(4 * DL, 4 * DL, 3.0),
        grid_spec=td.GridSpec.uniform(dl=DL),
        structures=structures,
        sources=[td.PlaneWave(
            center=(0, 0, 0.6), size=(td.inf, td.inf, 0), direction="-",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.1 * FREQ0))],
        monitors=[td.FluxMonitor(center=(0, 0, 0.9), size=(td.inf, td.inf, 0),
                                 freqs=list(FREQS), name="R")],
        run_time=1e-13,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.periodic(), y=td.Boundary.periodic(),
            z=td.Boundary.pml()))


def _structures(glass_low, glass_mid):
    """A glass lower half space, a cSi slab in the middle, and a thin glass sheet floating above.

    The bottom face of the cSi touches glass and the top face touches vacuum, so both kinds of
    "metal / constant dielectric" mixed cell appear and both legs of the candidate-by-candidate
    material-pair selection are exercised.
    """
    return [
        td.Structure(
            geometry=td.Box(center=(0, 0, _ZA - 5.0), size=(td.inf, td.inf, 10.0)),
            medium=glass_low),
        td.Structure(
            geometry=td.Box(center=(0, 0, 0.5 * (_ZA + _ZB)),
                            size=(td.inf, td.inf, _ZB - _ZA)),
            medium=_csi()),
        td.Structure(
            geometry=td.Box(center=(0, 0, 0.6 + 0.5 * DL + 0.1),
                            size=(td.inf, td.inf, 0.2)),
            medium=glass_mid),
    ]


def _flat(sc) -> dict[str, np.ndarray | None]:
    """Flatten the array fields of the Scene that take part in the criteria, for bitwise
    comparison."""
    out = {}
    for k in ("eps_ex", "eps_ey", "eps_ez", "sigma_ex", "sigma_ey", "sigma_ez"):
        out[k] = getattr(sc, k)
    d = sc.dispersion
    for k in ("comp", "cell", "pole_ofs", "am1", "b", "g"):
        out[f"disp.{k}"] = getattr(d, k)
    m = sc.dispersion_mix
    if m is not None:
        for k in ("comp", "cell", "p_ofs", "pa", "pb", "q_ofs", "qa", "qb",
                  "beta", "eps_inf", "zeta_inf"):
            out[f"mix.{k}"] = getattr(m, k)
    else:
        out["mix"] = None
    return out


def test_duplicate_medium_bitwise():
    """The same material under two names vs merged by hand into one: the Scenes are bitwise
    identical."""
    glass_a = td.Medium(permittivity=2.25, name="glass_a")
    glass_b = td.Medium(permittivity=2.25, name="glass_b")
    sc_dup = scene_mod.from_simulation(_sim(_structures(glass_a, glass_b)))
    sc_one = scene_mod.from_simulation(_sim(_structures(glass_a, glass_a)))

    a, b = _flat(sc_dup), _flat(sc_one)
    assert a.keys() == b.keys()
    for k in a:
        if a[k] is None or b[k] is None:
            assert a[k] is b[k], k
            continue
        assert np.array_equal(np.asarray(a[k]), np.asarray(b[k])), \
            f"{k} is not bitwise identical"


def test_proportional_constant_columns():
    """The vacuum and glass constant columns are collinear: after merging the scene builds, and the
    material pair is recovered as glass candidate by candidate."""
    glass = td.Medium(permittivity=2.5921)
    sim = _sim(_structures(glass, glass))
    media.set_local_subpixel(bool(sim.subpixel))

    mediums = media_columns(sim, axis=2)          # everything below is done on the Ez component
    assert len(mediums) == 3                      # vacuum, glass, cSi
    groups = column_groups(mediums, fit_frequencies(sim, len(mediums)))
    assert sorted(len(g) for g in groups) == [1, 2]   # {vacuum, glass} merge, cSi stays alone
    freqs = fit_frequencies(sim, len(groups))

    # Ez: both flat faces of the cSi slab carry the normal component, so mixed cells must appear
    w, mix, harm = decompose(sim, "Ez", freqs, mediums, groups)

    # The glass column is not the representative, so its weight is identically 0 and the ε of a
    # pure glass cell is recorded on the vacuum representative column
    assert np.all(w[1] == 0.0)
    ix, iy, iz = w.shape[1] // 2, w.shape[2] // 2, int(0.2 * w.shape[3])
    pure_glass = w[:, ix, iy, iz]                 # deep in the glass half space, far from any face
    assert abs(pure_glass[0] - 2.5921) < 1e-3 and abs(pure_glass[2]) < 1e-3

    # For non-mixed cells the reconstructed ε passes the residual criterion (a physical criterion;
    # it does not matter that the weights are not unique).
    # Tolerance 2xFIT_RTOL: decompose's criterion is row-weighted (1/max|ε|), while this is the bare
    # per-frequency relative difference. The metrics differ, so a boundary cell can be off by less
    # than a factor of two.
    eps_cols = np.stack([np.asarray(m.eps_model(freqs), dtype=np.complex128).ravel()
                         for m in mediums])
    for k, f in enumerate(freqs):
        obs = np.asarray(media.epsilon_complex(sim, "Ez", float(f)))
        rec = np.tensordot(eps_cols[:, k], w, axes=(0, 0))
        rel = np.abs(rec - obs) / np.maximum(np.abs(obs), 1.0)
        assert float(rel[~harm].max()) < 2 * FIT_RTOL

    # Mixed cells exist (the normal glass-cSi interface; on the vacuum-cSi face the harmonic mean is
    # nearly constant and least squares fits it away outright). The key assertion: the paired
    # material is **glass**. The group's representative column is vacuum, and without expanding
    # candidate by candidate these cells would be assigned to the vacuum-cSi pair and the residual
    # would blow up.
    assert mix is not None
    partners = set(int(v) for v in mix["mi"]) | set(int(v) for v in mix["mj"])
    assert partners == {1, 2}                     # glass + cSi, not the representative vacuum

    # End to end: the whole from_simulation is no longer stopped by the condition-number guard
    scene_mod.from_simulation(sim)


def test_widen_only_when_ill_conditioned():
    """A narrow band (an adjoint simulation only has monitors at the target frequency point) widens
    the sampling window around its centre until the condition number is back within COND_MAX; a
    window with enough bandwidth is returned unchanged, so existing examples keep their numbers.
    This is TidyFab0GC's combination of media."""
    from openem.scene.weights import COND_MAX, design_matrix, widen_if_ill_conditioned
    meds = [td.Medium(permittivity=1.0), td.material_library["cSi"]["Li1993_293K"],
            td.material_library["SiO2"]["Palik_LowLoss"]]

    wide = np.linspace(td.C_0 / 1.6, td.C_0 / 1.5, 6)        # forward: 1.5-1.6 µm
    assert widen_if_ill_conditioned(wide, meds, "Ex") is wide

    narrow = np.linspace(td.C_0 / 1.551, td.C_0 / 1.549, 6)  # adjoint: a 2 nm window
    assert np.linalg.cond(design_matrix(meds, narrow)) > COND_MAX
    out = widen_if_ill_conditioned(narrow, meds, "Ex")
    assert np.linalg.cond(design_matrix(meds, out)) <= COND_MAX
    assert len(out) == len(narrow)
    center = 0.5 * (narrow[0] + narrow[-1])
    assert abs(0.5 * (out[0] + out[-1]) - center) < 1e-6 * center   # still centred on the original
    assert out[-1] - out[0] <= 0.2 * center + 1e-6 * center           # at most ±10% wide
