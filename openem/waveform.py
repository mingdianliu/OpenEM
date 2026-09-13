# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Source time waveform and spectrum: **delegate to Tidy3D, do not reimplement**.

Conventions like ``GaussianPulse``'s ``offset`` and ``remove_dc_component`` are buried inside
Tidy3D, and reimplementing them is asking for bugs. Instead we call ``source_time.amp_time(t)``
and sample it step by step into a table. That is a few hundred thousand floats, and the
correctness comes for free.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Waveform:
    """A pre-sampled source waveform.

    Attributes:
        amp_int: real amplitude at whole steps ``n*dt``, shape ``(num_steps+1,)``.
        amp_half: real amplitude at half steps ``(n+0.5)*dt``.
        amp_int_complex: the raw complex values at whole steps. A magnetic dipole injects these on
            the complex-field path (Bloch); the real path takes only the real part.
        amp_half_complex: the raw complex values at half steps, used by the complex path for an
            electric dipole. Reading an older npz back may leave this None, in which case the
            complex path refuses (fail closed) and the fix is to re-export.
    """

    amp_int: np.ndarray
    amp_half: np.ndarray
    amp_int_complex: np.ndarray
    dt: float
    amp_half_complex: np.ndarray | None = None

    @property
    def num_steps(self) -> int:
        return self.amp_int.size - 1

    @property
    def peak_index(self) -> int:
        """Index of the whole step with the largest absolute amplitude. Only the tests use it."""
        return int(np.argmax(np.abs(self.amp_int)))


def sample(source_time, dt: float, num_steps: int) -> Waveform:
    """Sample ``source_time`` into a table of ``num_steps+1`` points.

    It samples the **nominal** step count, not the actual one: when shutoff ends the run early,
    the waveform used for normalization must not be truncated with it.
    """
    n = np.arange(num_steps + 1, dtype=np.float64)
    t_int = n * dt
    t_half = (n + 0.5) * dt
    a_int = np.asarray(source_time.amp_time(t_int))
    a_half = np.asarray(source_time.amp_time(t_half))
    return Waveform(
        amp_int=np.real(a_int).astype(np.float64),
        amp_half=np.real(a_half).astype(np.float64),
        amp_int_complex=a_int.astype(np.complex128),
        dt=dt,
        amp_half_complex=a_half.astype(np.complex128),
    )


def spectrum(source_time, dt: float, num_steps: int, freqs: np.ndarray) -> np.ndarray:
    """Source spectrum, complex, shape ``(len(freqs),)``, used by the second flux-normalization path.

    Tidy3D defaults to ``normalize_index=0``, meaning monitor results are normalized by source 0's
    spectrum.
    """
    times = np.arange(num_steps + 1, dtype=np.float64) * dt
    return np.atleast_1d(np.asarray(source_time.spectrum(times, np.asarray(freqs), dt)))
