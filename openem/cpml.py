# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""CPML profiles and recursion coefficients.

Implemented to **Tidy3D's own definition** rather than to the literature. Source:
``tidy3d/components/mode/derivatives.py:264-297``::

    s(x)     = kappa(x) + 1j * sigma(x) / (omega * EPSILON_0)
    kappa(x) = kappa_min + (kappa_max - kappa_min) * step ** order
    sigma(x) = sigma_max * avg_speed / (ETA_0 * dl) * step ** order

``dl`` is the **grid step**, not the PML thickness, and ``avg_speed`` is the wave speed inside the
PML divided by c0.

The recursion coefficients use the Roden-Gedney form, algebraically identical to FDTDX
``perfectly_matched_layer.py:123-127``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

#: SI constants. These are exactly equal to tidy3d.constants once converted; tests/test_cpml.py
#: locks that. tidy3d works in um, so EPSILON_0[F/um] * 1e6 is the EPSILON_0[F/m] used here.
EPSILON_0 = 8.8541878128e-12       # F/m
MU_0 = 1.25663706212e-6            # H/m
ETA_0 = 376.7303136668535          # ohm, independent of the unit system
C_0 = 299792458.0                  # m/s

Side = Literal["lo", "hi"]


@dataclass(frozen=True)
class PMLParams:
    """Tidy3D's ``PMLParams``, submitted with simulation.json.

    The defaults here are ``td.PML``'s own presets.
    """

    num_layers: int
    sigma_min: float = 0.0
    sigma_max: float = 1.5
    sigma_order: int = 3
    kappa_min: float = 1.0
    kappa_max: float = 3.0
    kappa_order: int = 3
    alpha_min: float = 0.0
    alpha_max: float = 0.0
    alpha_order: int = 1


def steps(n_pml: int, side: Side, at: Literal["E", "H"]) -> np.ndarray:
    """Normalized depth ``step = x/d`` per layer, ordered by increasing global cell index.

    Reproduces the values from tidy3d's ``create_sfactor_f`` (H positions, forward difference) and
    ``create_sfactor_b`` (E positions, backward difference).

    The E side is **asymmetric** between min and max: at the min end the outermost layer has
    step=1.0, at the max end the innermost has step=0. That falls out of the backward difference
    being off by half a cell. Tidy3D writes it that way, so we copy it rather than correct it.

    Returns:
        ``(n_pml,)``. ``side="lo"`` maps to global indices 0..n-1, ``"hi"`` to N-n..N-1.
    """
    n = n_pml
    if n <= 0:
        return np.zeros(0)
    i = np.arange(n, dtype=np.float64)
    if side == "lo":
        # min end: depth grows from the inside out, hence (n - i)
        return (n - i) / n if at == "E" else (n - i - 0.5) / n
    # max end: j = i - (N - n) runs 0..n-1, and depth grows with j
    return i / n if at == "E" else (i + 0.5) / n


def profile(
    params: PMLParams,
    dl: float,
    side: Side,
    at: Literal["E", "H"],
    avg_speed: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute ``(sigma, kappa, alpha)`` for each PML layer.

    Args:
        dl: **the representative grid step at this end, in metres**. Tidy3D uses only ``dls[0]``
            at the min end and only ``dls[-1]`` at the max end; even on a non-uniform grid it is
            one value, not one per layer.
        at: ``"E"`` or ``"H"``; the two are offset by half a cell.
        avg_speed: wave speed inside the PML divided by c0, i.e. ``1/n``.

    Returns:
        Three ``(num_layers,)`` arrays: sigma [S/m], kappa, alpha [S/m].
    """
    s = steps(params.num_layers, side, at)
    sigma = (
        params.sigma_max * avg_speed / (ETA_0 * dl) * s**params.sigma_order
    )
    if params.sigma_min:
        # Tidy3D's s_value has no min term for sigma; it always starts from 0. There is no
        # evidence for how it would handle a non-zero sigma_min, so raise rather than guess.
        raise NotImplementedError(
            f"sigma_min={params.sigma_min} is not 0. In Tidy3D's s_value(), sigma always starts "
            "from 0, and there is no evidence for how a non-zero sigma_min should be graded. "
            "Confirm the intended behaviour before implementing it."
        )
    kappa = params.kappa_min + (params.kappa_max - params.kappa_min) * s**params.kappa_order
    alpha = params.alpha_min + (params.alpha_max - params.alpha_min) * s**params.alpha_order
    return sigma, kappa, alpha


def recursion_coeffs(
    sigma: np.ndarray, kappa: np.ndarray, alpha: np.ndarray, dt: float
) -> tuple[np.ndarray, np.ndarray]:
    """Recursion coefficients ``(a, b)`` for psi: ``psi_new = b*psi + a*d_field``.

    Layers with ``sigma=0`` get ``a=0``, not nan.
    """
    b = np.exp(-dt / EPSILON_0 * (sigma / kappa + alpha))
    denom = sigma + alpha * kappa
    with np.errstate(divide="ignore", invalid="ignore"):
        a = (b - 1.0) * sigma / denom / kappa
    return np.nan_to_num(a, nan=0.0, posinf=0.0, neginf=0.0), b


@dataclass(frozen=True)
class PMLCoeffs:
    """Every coefficient for one PML face, ready to hand to a kernel.

    All arrays have length ``num_layers``. ``inv_kappa`` stores the reciprocal so the kernel only
    ever multiplies.
    """

    a_E: np.ndarray
    b_E: np.ndarray
    inv_kappa_E: np.ndarray
    a_H: np.ndarray
    b_H: np.ndarray
    inv_kappa_H: np.ndarray

    @property
    def num_layers(self) -> int:
        return self.a_E.size


def build(
    params: PMLParams, dl: float, side: Side, dt: float, avg_speed: float = 1.0
) -> PMLCoeffs:
    """Do it all in one step: parameters -> profile -> recursion coefficients.

    Every coefficient is computed on the Python side, which is why the kernel contains no ``powf``
    and no boundary branches.
    """
    out = {}
    for at in ("E", "H"):
        sigma, kappa, alpha = profile(params, dl, side, at, avg_speed)  # type: ignore[arg-type]
        a, b = recursion_coeffs(sigma, kappa, alpha, dt)
        out[f"a_{at}"] = a
        out[f"b_{at}"] = b
        out[f"inv_kappa_{at}"] = 1.0 / kappa
    return PMLCoeffs(**out)  # type: ignore[arg-type]
