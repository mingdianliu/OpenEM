# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Source waveform and spectrum: confirm that we delegate to Tidy3D rather than reimplement it."""

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem import waveform


@pytest.fixture(scope="module")
def source_time():
    """Source waveform parameters taken from a real example."""
    return td.GaussianPulse(
        freq0=361221792828282.9,
        fwidth=56238123434343.44,
        offset=5.0,
        remove_dc_component=True,
    )


DT = 4.766437173827506e-17


def test_sample_matches_tidy3d_exactly(source_time):
    """The sampled table must equal amp_time bitwise: we do not reimplement the waveform."""
    w = waveform.sample(source_time, DT, 200)
    t = np.arange(201) * DT
    expect = np.asarray(source_time.amp_time(t))
    np.testing.assert_array_equal(w.amp_int_complex, expect.astype(np.complex128))
    np.testing.assert_array_equal(w.amp_int, np.real(expect))


def test_half_step_offset(source_time):
    """The half-step table is sampled at (n+0.5)*dt."""
    w = waveform.sample(source_time, DT, 200)
    expect = np.real(np.asarray(source_time.amp_time((np.arange(201) + 0.5) * DT)))
    np.testing.assert_array_equal(w.amp_half, expect)


def test_sample_covers_nominal_steps(source_time):
    """Sample num_steps+1 points: when shutoff ends the run early, the waveform used for
    normalization must not be truncated with it.
    """
    w = waveform.sample(source_time, DT, 1000)
    assert w.amp_int.size == 1001
    assert w.num_steps == 1000


def test_remove_dc_component_kills_zero_frequency(source_time):
    """The physical consequence of remove_dc_component=True: the DC component is suppressed.

    This is exactly why GaussianPulse should not be reimplemented; conventions like this one are
    buried inside Tidy3D.
    """
    n = 20000
    spec = waveform.spectrum(source_time, DT, n, np.array([0.0, source_time.freq0]))
    dc, at_f0 = np.abs(spec)
    assert dc < 1e-3 * at_f0, f"DC {dc:.3e} should be suppressed relative to {at_f0:.3e} at f0"


def test_spectrum_shape_and_finite(source_time):
    freqs = np.array([3.25e14, 3.43e14, 3.61e14, 3.79e14, 3.97e14])
    spec = waveform.spectrum(source_time, DT, 5000, freqs)
    assert spec.shape == (5,)
    assert np.iscomplexobj(spec)
    assert np.isfinite(spec).all()
    assert np.all(np.abs(spec) > 0)


def test_peak_index_near_offset(source_time):
    """A GaussianPulse should peak near offset/fwidth, neither at t=0 nor at the end."""
    n = 20000
    w = waveform.sample(source_time, DT, n)
    # Tidy3D's definition: twidth = 1/(2*pi*fwidth), with the peak at t = offset*twidth
    twidth = 1.0 / (2 * np.pi * source_time.fwidth)
    expect_step = source_time.offset * twidth / DT
    assert w.peak_index == pytest.approx(expect_step, rel=0.05)
