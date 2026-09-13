# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Apodization: windowing the DFT.

Both main criteria are parameter-free:

* **Window equivalence**: put a FieldTimeMonitor (which stores the raw time series) and an
  apodized FieldMonitor (which accumulates the windowed sum on the device) at the same point, run
  the same windowed DFT over the time series with numpy, and the two must agree. Only E is
  compared: H in a time monitor follows the time-averaging convention, while H in the DFT uses the
  raw half-step value, so the two are not supposed to be equal.
* **An all-None spec is bitwise identical**: ``ApodizationSpec`` with everything None and no spec
  at all go down the same code path (the window is identically 1.0, and multiplying by 1.0 is exact
  under IEEE), so the results must be bitwise identical.
"""

from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import apodization as apod_mod
from openem import scene as scene_mod
from openem.model import COMPONENTS

FREQ0 = 3.0e14
#: The probe point is deliberately off the dipole, so the field is nonzero and not singular
POINT = (0.15, 0.0, 0.0)


def _sim(monitors, run_time=8e-14):
    return td.Simulation(
        size=(1.2, 1.2, 1.2), grid_spec=td.GridSpec.uniform(dl=0.04),
        sources=[td.PointDipole(
            center=(0, 0, 0), polarization="Ez",
            source_time=td.GaussianPulse(freq0=FREQ0, fwidth=0.2 * FREQ0))],
        monitors=monitors,
        run_time=run_time,
        boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=10)))


# ------------------------------------------------------------ the window function itself

def test_window_flat_top_and_gaussian_ramp():
    start, end, width = 2e-14, 5e-14, 1e-14
    t = np.array([start - width, start, 3.5e-14, end, end + 2 * width])
    w = apod_mod.window(t, (start, end, width))
    np.testing.assert_allclose(w[1:4], 1.0, rtol=0)          # flat top, endpoints included
    np.testing.assert_allclose(w[0], np.exp(-apod_mod.APOD_K), rtol=1e-12)
    np.testing.assert_allclose(w[4], np.exp(-4 * apod_mod.APOD_K), rtol=1e-12)


def test_window_none_is_exact_ones():
    t = np.linspace(0, 1e-13, 7)
    assert (apod_mod.window(t, None) == 1.0).all()


# ------------------------------------------------------------ refusals and extraction

def test_width_only_normalizes_to_none():
    """A width alone, with no ramp, is mathematically no apodization at all: normalize to None."""
    mon = td.FieldMonitor(center=POINT, size=(0, 0, 0), freqs=[FREQ0], name="f",
                          apodization=td.ApodizationSpec(width=1e-14))
    sc = scene_mod.from_simulation(_sim([mon]))
    assert sc.field_monitors[0].apodization is None


def test_missing_width_fails_closed():
    """A ramp with no width leaves the window undefined. tidy3d may refuse this itself, but failing
    closed must not depend on upstream."""
    from types import SimpleNamespace
    from openem.scene import monitors as mons_mod
    fake = SimpleNamespace(
        name="m", apodization=SimpleNamespace(start=1e-14, end=None, width=None))
    with pytest.raises(ValueError, match="width"):
        mons_mod._apodization(fake)


def test_cavityfom_spec_survives_to_scene():
    """CavityFOM's field_vol is the case this unlocks: the spec must reach the Scene untouched."""
    start, width = 5.52382141648136e-13, 5.5238214164813604e-14
    mon = td.FieldMonitor(
        center=POINT, size=(0, 0, 0), freqs=[FREQ0], name="f",
        apodization=td.ApodizationSpec(start=start, width=width))
    sc = scene_mod.from_simulation(_sim([mon]))
    assert sc.field_monitors[0].apodization == (start, None, width)


def test_projection_monitor_apodization_reaches_all_faces():
    """Apodization on a projection monitor must reach all six face monitors unchanged (this used to
    be a loud refusal)."""
    mon = td.FieldProjectionAngleMonitor(
        center=(0, 0, 0), size=(0.6, 0.6, 0.6), freqs=[FREQ0], name="n2f",
        theta=[0.0], phi=[0.0], proj_distance=1e6,
        apodization=td.ApodizationSpec(start=1e-14, width=1e-14))
    sc = scene_mod.from_simulation(_sim([mon]))
    assert len(sc.field_monitors) == 6
    for fm in sc.field_monitors:
        assert fm.apodization == (1e-14, None, 1e-14)


# ------------------------------------------------------------ serialization (for solver jobs)

def test_serialize_roundtrip_keeps_apodization(tmp_path):
    """A solver job rebuilds the Scene from an npz: apodization, dipoles and FieldMonitors all have
    to come back unchanged."""
    from openem import serialize
    spec = td.ApodizationSpec(start=8e-14, end=1.4e-13, width=2e-14)
    sc = scene_mod.from_simulation(_sim([
        td.FieldMonitor(center=POINT, size=(0, 0, 0), freqs=[FREQ0], name="f",
                        apodization=spec)]))
    back = serialize.load(serialize.save(sc, tmp_path / "s.npz"))
    bm, m = back.field_monitors[0], sc.field_monitors[0]
    assert (bm.name, bm.origin, bm.box, bm.apodization) == \
           (m.name, m.origin, m.box, m.apodization)
    np.testing.assert_array_equal(bm.freqs, m.freqs)
    assert back.dipoles[0].component == sc.dipoles[0].component
    np.testing.assert_array_equal(back.dipoles[0].indices, sc.dipoles[0].indices)
    np.testing.assert_array_equal(back.dipoles[0].coef, sc.dipoles[0].coef)
    np.testing.assert_array_equal(back.dipoles[0].waveform.amp_half,
                                  sc.dipoles[0].waveform.amp_half)


def test_serialize_roundtrips_mode_monitors(tmp_path):
    """serialize now supports every Scene field, so this test changed from "must refuse" to "if it
    can be written it must read back": mode_monitors was the last unsupported field (the old
    version asserted a raise here), and now the assertion is that ModeMonitorSpec survives the
    round trip field by field. Silently dropping a monitor is still the most dangerous failure
    mode; only the guard changed, from a refusal to a round-trip identity.
    """
    from openem import serialize
    sc = scene_mod.from_simulation(_sim([
        td.ModeMonitor(center=POINT, size=(0.4, 0.4, 0), freqs=[FREQ0],
                       mode_spec=td.ModeSpec(num_modes=1), name="m")]))
    assert sc.mode_monitors, "precondition: the scene really does hold a ModeMonitor"
    back = serialize.load(serialize.save(sc, tmp_path / "s.npz"))
    assert [(m.name, m.plane_name) for m in back.mode_monitors] ==         [(m.name, m.plane_name) for m in sc.mode_monitors]


# ------------------------------------------------------------ end to end (GPU)

@pytest.fixture(scope="module")
def windowed_run():
    """One run covering it all: the time series, an apodized FieldMonitor, and the all-None /
    no-spec pair."""
    pytest.importorskip("cupy")
    from openem import solver

    freqs = [FREQ0, 1.1 * FREQ0]
    # Put the window ramp where the pulse energy is densest, so the window really bites into the
    # main lobe
    spec = td.ApodizationSpec(start=8e-14, end=1.4e-13, width=2e-14)
    mons = [
        td.FieldTimeMonitor(center=POINT, size=(0, 0, 0), interval=1,
                            fields=["Ex", "Ey", "Ez"], name="t"),
        td.FieldMonitor(center=POINT, size=(0, 0, 0), freqs=freqs, name="f_apod",
                        apodization=spec),
        td.FieldMonitor(center=POINT, size=(0, 0, 0), freqs=freqs, name="f_plain"),
        td.FieldMonitor(center=POINT, size=(0, 0, 0), freqs=freqs, name="f_allnone",
                        apodization=td.ApodizationSpec()),
    ]
    sc = scene_mod.from_simulation(_sim(mons, run_time=2e-13))
    # return_fields=True: the last step's E sample is missing from the time-domain buffer (its slot
    # lacks the second H half-kick and is cut by time_slots), so take it from the final fields
    res = solver.run(sc, use_shutoff=False, verbose=False, return_fields=True)
    return sc, res


def test_windowed_dft_matches_numpy(windowed_run):
    """The windowed DFT on the device == the windowed DFT numpy computes from the time series
    (E components, ~1e-5 relative)."""
    sc, res = windowed_run
    mon = next(m for m in sc.field_monitors if m.name == "f_apod")
    tmon = sc.field_time_monitors[0]
    assert mon.origin == tmon.origin and mon.box == tmon.box

    n_steps = res.steps_run
    buf = res.time_samples["t"]                      # (3, slots, ni, nj, nk)
    # The time slots only reach n_steps-1 (the last slot lacks H's second half-kick and is cut),
    # so E^{n_steps} comes from the final fields
    assert buf.shape[1] >= n_steps, "the time series must cover all accumulation steps but the last"

    dt = sc.dt
    m = np.arange(1, n_steps + 1)                    # E^{m} sits at t = m·dt
    w = apod_mod.window(m * dt, mon.apodization)
    assert w.min() < 0.9 and w.max() == 1.0, "the window must actually do something"

    i0, j0, k0 = mon.origin
    ni, nj, nk = mon.box
    phasor = res.field_phasors["f_apod"]             # (6, nf, ni, nj, nk)
    for ci, c in enumerate(tmon.comps):
        assert COMPONENTS[c] in ("Ex", "Ey", "Ez")
        last = res.fields[COMPONENTS[c]][i0:i0 + ni, j0:j0 + nj, k0:k0 + nk]
        series = np.concatenate(
            [buf[ci, 1: n_steps], last[None]], axis=0).astype(np.float64)
        for fi, f in enumerate(mon.freqs):
            phase = np.exp(1j * 2 * np.pi * f * m * dt)
            ref = np.tensordot(w * phase, series, axes=(0, 0))
            got = phasor[c, fi]
            scale = np.max(np.abs(ref))
            assert scale > 0
            err = np.max(np.abs(got - ref)) / scale
            assert err < 1e-5, f"{COMPONENTS[c]} f{fi} relative residual {err:.3e}"


def test_all_none_spec_bitwise_identical(windowed_run):
    """An all-None spec and no spec must be bitwise identical: direct evidence of a single code
    path."""
    _, res = windowed_run
    np.testing.assert_array_equal(res.field_phasors["f_allnone"],
                                  res.field_phasors["f_plain"])


def test_apodized_differs_from_plain(windowed_run):
    """Reverse sanity check: the window really does change the result, otherwise the equivalence
    test above would be vacuous."""
    _, res = windowed_run
    a = res.field_phasors["f_apod"]
    b = res.field_phasors["f_plain"]
    assert np.max(np.abs(a - b)) > 0
