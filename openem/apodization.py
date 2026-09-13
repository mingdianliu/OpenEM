# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Apodization: the time window applied to a frequency-domain monitor's DFT accumulation.

Tidy3D's ``ApodizationSpec(start, end, width)``, all in seconds, defines a flat-top window with
Gaussian ramps: the window is 1 for ``start <= t <= end`` and decays as a Gaussian on either side::

    t < start:  exp(-K * ((t - start) / width)^2)
    t > end:    exp(-K * ((t - end) / width)^2)

The window multiplies only the monitor's DFT accumulation. It does not touch the source, nor the
source spectrum used for normalization.
"""

from __future__ import annotations

import numpy as np

#: Exponent coefficient K of the Gaussian ramp. The Tidy3D client carries no window formula; it
#: only defines the spec, and the window is applied by the solver. The reference output cannot
#: settle K either: 0.5 and 1.0 give shape residuals that agree to four significant figures,
#: because the ramp covers only about 1% of the integration interval. We take 0.5, the standard
#: Gaussian convention.
APOD_K = 0.5


def window(times: np.ndarray, spec: tuple | None) -> np.ndarray:
    """Map an array of times in seconds to window values.

    ``spec`` is ``(start, end, width)``, or ``None``. ``None`` -- and equally a spec whose start
    and end are both None, see monitors._apodization -- returns all ones, which is bitwise
    equivalent to no apodization at all.
    """
    t = np.asarray(times, dtype=np.float64)
    w = np.ones_like(t)
    if spec is None:
        return w
    start, end, width = spec
    if start is not None:
        m = t < start
        w[m] = np.exp(-APOD_K * ((t[m] - start) / width) ** 2)
    if end is not None:
        m = t > end
        w[m] = np.exp(-APOD_K * ((t[m] - end) / width) ** 2)
    return w
