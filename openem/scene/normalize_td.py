# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""The source-spectrum divisor that takes our raw DFT to tidy3d monitor data: the global
normalization convention.

Shared by :func:`projection.field_data`, the docstring of ``modes.amplitudes`` and
``tests/test_dipole_h.py``; projection re-exports it under the same name.
"""

from __future__ import annotations

import numpy as np
import tidy3d as td


def normalization(sim: td.Simulation, freqs: np.ndarray) -> np.ndarray:
    """Conversion factor from our raw DFT to Tidy3D monitor data.

    We accumulate the bare sum_n F_n e^{i omega t_n}. Tidy3D stores that multiplied by
    dt/sqrt(2*pi) and divided by the **source spectrum**, of whichever source normalize_index
    names (source 0 by default). In a flux comparison this factor cancels between the structure run
    and the reference run, which is why it was not needed before; a far-field projection is an
    absolute quantity, so it has to be applied.

    **Only the source-spectrum step happens here.** The far field additionally multiplies by
    :data:`FAR_FIELD_SCALE`; **flux must not**. Folding that factor in here made one flux come out
    1e-12 too small.

    That spectrum carries dt/sqrt(2*pi) is locked by
    test_dipole.py::test_spectrum_is_dt_dft_over_sqrt_two_pi.
    """
    idx = 0 if sim.normalize_index is None else int(sim.normalize_index)
    st = sim.sources[idx].source_time
    times = np.arange(sim.num_time_steps + 1, dtype=np.float64) * sim.dt
    spec = np.asarray(st.spectrum(times, np.asarray(freqs, dtype=np.float64), sim.dt))
    return (sim.dt / np.sqrt(2 * np.pi)) / spec
