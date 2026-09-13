# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Materials to pole tables. **About materials only; the grid does not enter here.**

Tidy3D reduces every dispersive material to a PoleResidue (``_pole_residue_dict`` is an abstract
method on ``DispersiveMedium`` that every subclass must implement), so only one model is
recognized here. There are two pole tables:

* **poles of eps**, used by the arithmetic-averaging branch
* **poles of 1/eps**, used by the harmonic-averaging branch

The derivation, the measured convergence order, and why the real second-order form is not used,
are all recorded alongside the ADE kernel.
"""

from __future__ import annotations

import numpy as np
import tidy3d as td
from tidy3d.components.medium import AnisotropicMediumFromMedium2D

from openem import cpml, knobs


#: Diagonally anisotropic media admitted by exact type. ``AnisotropicMediumFromMedium2D`` is the
#: equivalent volumetric type (a subclass) that ``sim.volumetric_structures`` produces from a
#: ``Medium2D``. The whitelist matches on exact type, so it has to be listed explicitly.
ANISOTROPIC = (td.AnisotropicMedium, AnisotropicMediumFromMedium2D,
               *(getattr(td, n) for n in ("CustomAnisotropicMedium",) if hasattr(td, n)))


def medium_parts(m) -> list:
    """The scalar parts of a medium: three diagonal components for an ``AnisotropicMedium``,
    otherwise the medium itself.

    Support for diagonal anisotropy rests on this: for an ``AnisotropicMedium``,
    ``sim.epsilon(coord_key='Ex')`` returns the xx component sampled at the Ex position, and
    likewise for Ey and Ez; both subpixel conventions were confirmed by measurement. So component
    by component it behaves as three ordinary media.
    """
    if type(m) in ANISOTROPIC:
        return [m.xx, m.yy, m.zz]
    return [m]


#: Medium types treated as staircased PEC. ``LossyMetalMedium``, which Tidy3D models as a good
#: conductor through a Leontovich surface-impedance boundary (the copper in the microwave
#: examples), is **approximated as PEC** here: at GHz the skin depth of copper is of order um, far
#: below the grid step of more than 100 um, so the ohmic loss is a second-order effect on the
#: transmission and reflection spectra. A known accuracy concession, printed once when the scene is
#: built.
PEC_LIKE = (td.PECMedium,
            *(getattr(td, n) for n in ("LossyMetalMedium",) if hasattr(td, n)))


def is_pec_like(medium) -> bool:
    return type(medium) in PEC_LIKE


#: Medium types that reduce to PoleResidue. ``_pole_residue_dict`` is an abstract method on
#: ``DispersiveMedium`` that every subclass must implement, so this list amounts to "every
#: dispersive medium".
DISPERSIVE = (td.PoleResidue, td.Lorentz, td.Drude, td.Sellmeier, td.Debye)

#: A pole counts as DC when ``|a|`` is below this factor times the reciprocal of one time step. A
#: Drude DC pole is **exactly** 0, so this threshold only guards against floating-point noise and
#: its precise value does not affect the outcome.
DC_POLE_RTOL = 1e-9

def sigma_from_eps(eps: np.ndarray, freq: float) -> np.ndarray:
    """Solve back for the conductivity in S/m from a complex relative eps.

    Tidy3D's definition (``medium.py:786``) is ``eps = eps_real + 1j*sigma/(omega*eps0)``, so
    ``sigma = Im(eps_r) * omega * eps0``. This sigma acts on **E**, the textbook convention. Meep's
    sigma acts on D and differs by a factor of eps; the two must not be mixed.

    Substituting our own SI eps0 gives S/m directly, with no um-to-m conversion needed, because
    ``Im(eps_r)`` is dimensionless.
    """
    return np.imag(eps) * (2.0 * np.pi * freq) * cpml.EPSILON_0


#: Spatially varying Custom dispersive media: per-cell residues at fixed pole positions.
CUSTOM_DISPERSIVE = tuple(
    getattr(td, n) for n in ("CustomLorentz", "CustomDrude", "CustomSellmeier",
                             "CustomPoleResidue") if hasattr(td, n))


def _custom_scalar_pole(med):
    """A spatially varying Custom dispersive medium, reduced to a **scalar reference** PoleResidue
    with fixed poles and reference residues.

    In the CustomLorentz cases in the validation set the pole positions (f0 and gamma, or a) are
    uniform across cells and only the residue (delep) varies in space. A reference residue is used
    to build the scalar poles; the per-cell delep is then recovered by the least-squares fit in the
    dispersion pipeline as a material weight. With a material table of [background, this reference
    pole], the fit yields w = delep_cell/delep_ref, and the pipeline scales the pole by w so that
    B_cell is proportional to delep_cell, which is correct.
    """
    if type(med).__name__ == "CustomLorentz":
        eps_inf = float(np.asarray(med.eps_inf.values).ravel()[0]) \
            if hasattr(med.eps_inf, "values") else float(med.eps_inf)
        coeffs = []
        for de, f0, g in med.coeffs:
            de_v = np.asarray(de.values)
            dr = float(de_v[de_v != 0].max()) if np.any(de_v != 0) else 1.0
            f0s = float(np.asarray(f0.values).ravel()[0])
            gs = float(np.asarray(g.values).ravel()[0])
            coeffs.append((dr, f0s, gs))
        return td.Lorentz(eps_inf=eps_inf, coeffs=tuple(coeffs)).pole_residue
    if type(med).__name__ == "CustomPoleResidue":
        a = np.asarray(med.poles[0][0].values)
        c = np.asarray(med.poles[0][1].values)
        eps_inf = float(np.asarray(med.eps_inf.values).ravel()[0])
        cr = c[np.abs(c) == np.abs(c).max()].ravel()[0]
        ar = a.ravel()[np.argmax(np.abs(c))]
        return td.PoleResidue(eps_inf=eps_inf, poles=[(complex(ar), complex(cr))])
    # For CustomDrude and CustomSellmeier, how to pick the scalar reference (which cell represents
    # the medium) is undecided. The old ``med.pole_residue.to_medium.pole_residue`` route was never
    # verified and no example reaches it, so raise loudly rather than silently take an unverified
    # path.
    raise NotImplementedError(
        f"how to pick the scalar reference poles of a {type(med).__name__} is undecided; the "
        "validation set contains only CustomLorentz and CustomPoleResidue")


def pole_residue(med):
    """Any dispersive medium reduced to a ``PoleResidue``; returns ``None`` if non-dispersive."""
    if isinstance(med, td.PoleResidue):
        return med
    if type(med).__name__ == "CustomMedium":
        return None                    # spatially static: eps comes from sim.epsilon, no poles
    if CUSTOM_DISPERSIVE and isinstance(med, CUSTOM_DISPERSIVE):
        return _custom_scalar_pole(med)
    if isinstance(med, DISPERSIVE):
        return med.pole_residue
    return None

def any_dispersive(sim: td.Simulation) -> bool:
    """Whether the scene contains any dispersive medium, checking an ``AnisotropicMedium`` axis by
    axis.

    It looks at ``volumetric_structures`` rather than ``structures``, because the equivalent
    volumetric medium that a ``Medium2D`` reduces to (anisotropic plus PoleResidue) is what the
    solver actually runs. Without a Medium2D the two lists are the same.
    """
    for med in [sim.medium] + [s.medium for s in sim.volumetric_structures]:
        # isinstance is used deliberately here rather than the exact-type table of medium_parts:
        # CustomAnisotropicMediumInternal is a subclass of AnisotropicMedium but is not in
        # ANISOTROPIC, and the other spelling would miss its components.
        parts = ([med.xx, med.yy, med.zz]
                 if isinstance(med, td.AnisotropicMedium) else [med])
        if any(isinstance(p, DISPERSIVE) for p in parts):
            return True
    return False

def split_dc_pole(pr, dt: float) -> tuple[list, float]:
    """``(non-DC poles, equivalent conductivity in S/m)``.

    A pole at ``a = 0`` is **pure conductivity**, not an oscillator: that pair contributes

        c/(jω) + c*/(jω) = 2Re(c)/(jω)   →   Im ε = 2Re(c)/ω

    Matching that against Tidy3D's ``eps = eps' + i*sigma/(omega*eps0)`` gives
    ``sigma = 2*Re(c)*eps0``. A Drude medium converted to PoleResidue **always** has such a pole;
    measured, ``a`` is exactly 0j and the equivalent sigma is 31.91 S/m.

    It has to be moved out, because ``a = 0`` gives ``A = 1`` and the recursion degenerates to
    ``P += B(E^{n+1}+E^n)``, **an integrator that never decays**. Run that in float32 for a few
    thousand steps and the error only accumulates. Moving it onto the sigma path instead puts it on
    the semi-implicit scheme already verified against the reference outputs, which is
    unconditionally stable.
    """
    keep, sigma = [], 0.0
    scale = 1.0 / dt
    for a, c in pr.poles:
        q, r = complex(a), complex(c)
        if abs(q) <= DC_POLE_RTOL * scale:
            sigma += 2.0 * r.real * cpml.EPSILON_0
        else:
            keep.append((q, r))
    return keep, sigma

#: Test for a "lossless pole": ``|Re a| <= this factor * |a|``. Measured magnitudes on both sides:
#: nominally lossless fitted materials sit at floating-point noise, 2.7e-16 for one silicon fit and
#: 7.2e-10 / 1.4e-9 for a silica one, while genuinely lossy materials are at least about 1e-2.
#: 1e-6 is at least three orders of magnitude from either side. Note that ZETA_AXIS_RTOL=1e-9
#: cannot be reused here: the silica fit's second pole, at 1.448e-9, sits right on that line.
LOSSLESS_AXIS_RTOL = 1e-6

#: Artificial loss injected into a lossless pole, as a relative amount: ``Re a <- -gamma_rel*|a|``.
#: The magnitude argument:
#:
#: * **Why inject anything**: after trapezoidal discretization a lossless pole has |A| identically
#:   1, i.e. marginal stability; check_pole_gain measures all 15 digits as 1. The isolated
#:   recursion stays bounded in fp32 over a long run (tests/test_zeta_stability.py), but the full
#:   GPU solver carries no such guarantee. On a slot waveguide of lossless silicon and silica run
#:   for 62,028 steps, a GPU probe bisection measured divergence to NaN: replacing the medium with
#:   a non-dispersive one was stable, removing the slot was stable, and slot plus lossless
#:   dispersion diverged, independently of the y boundary type. A term-by-term CPU reproduction,
#:   CPML and fp32 included, did not reproduce it, and the amplifying step was never localized on
#:   the GPU side. Marginal stability is a necessary condition for it, so pushing |A| inside the
#:   unit circle suppresses any marginal gain under an exponential decay. That is a backstop
#:   independent of the mechanism.
#: * **Upper bound, the physical tolerance**: the perturbation gamma causes in eps(omega) is about
#:   |c|*gamma/|omega-omega_0|^2. The working band at 193 THz is far enough from the silicon pole at
#:   1.02e15 Hz that at gamma_rel=1e-6 the change is about 8e-6 absolute, 7e-7 relative, far below
#:   any observable tolerance; the corresponding artificial absorption is below 2e-5 over a 2.7 um
#:   propagation distance.
#: * **Lower bound, representable in fp32 and enough to suppress the divergence**: the coefficient
#:   table stores A-1 in fp32, so the injected gamma*dt = gamma_rel*(|a|*dt) must exceed the fp32
#:   rounding of |A-1| ~ |a|*dt, which is 1.2e-7 relative. gamma_rel=1e-7 would lose half of itself
#:   to rounding; 1e-6 clears it by a factor of 8. The cumulative decay over 62k steps is
#:   e^(-gamma_rel*|a|*dt*N): with gamma_rel=1e-6, |a|*dt=0.31 and N=62,028 that is 1.9%, and it
#:   acts only on freely ringing pole states, which is exactly the numerical noise being
#:   suppressed. The driven response depends only on eps(omega); see the upper bound above.
#: * **Consistent with Tidy3D**: the tidy3d client does not allow zero damping on an "ideal
#:   lossless inductor" either, taking damping = frequency/1e6 (LOSS_FACTOR_INDUCTOR = 1e6 in
#:   components/lumped_element.py), the same relative magnitude of 1e-6.
#:
#: The environment variable ``OPENEM_LOSSLESS_DAMPING_REL`` overrides it; setting 0 restores the old
#: behaviour bitwise. A lossy pole is always returned unchanged, bitwise.
LOSSLESS_DAMPING_REL = 1e-6

def _lossless_damping_rel() -> float:
    return float(knobs.env("LOSSLESS_DAMPING_REL"))      # default registered in knobs.KNOBS, equal to LOSSLESS_DAMPING_REL

#: Extra damping for a "low-frequency lossless pole". Default **0, i.e. off and bitwise unchanged**.
#:
#: **Motivation**: on one scene the low-frequency CPML-by-ADE instability mode (spectral peak at
#: omega*dt about 0.0084) resonates with the low-frequency pole of the lossless silica fit
#: (|a|*dt about 0.0083). A gamma sweep gave a dose slope d(growth rate)/d(gamma_rel) of 8.6e-3 per
#: step, about that pole's |a|*dt, i.e. a participation factor near 1, so the unstable mode rides
#: on that pole. A uniform relative damping of gamma_rel*|a| is inherently unfavourable for a
#: low-|a| pole: damping it to about 1e-3 per step would need gamma_rel about 0.1, and at that level
#: the dominant silicon pole, whose |c| is 300 times larger and which sets eps in the working band,
#: would cost about 50% of the field amplitude. Injecting damping relative to its own |a| into
#: **only** the low-frequency poles (|a|*dt < :data:`LOSSLESS_LOWFREQ_ADT_MAX`, far below the
#: working band and contributing almost nothing to eps(omega_0)) is two orders of magnitude
#: cheaper: at gamma_lf=1.0 the damping is 8.3e-3 per step, with |delta eps| at most 7.4e-3 in the
#: working band and a 2.5% field-amplitude change over 2.7 um of propagation (gamma_lf=0.3 gives
#: 2.5e-3 per step and 0.77%).
#:
#: The test needs the grid dt (``stabilize_lossless(poles, dt=...)``); without dt, or with the
#: environment variable unset, the whole mechanism is off. The variable is
#: ``OPENEM_LOSSLESS_LOWFREQ_DAMPING_REL``.
LOSSLESS_LOWFREQ_DAMPING_REL = 0.0

#: A lossless pole counts as "low frequency" when ``|a|*dt`` falls below this value. Measured on one
#: grid: the silica low-frequency pole is 0.0083, the dominant silicon pole 0.31 and the dominant
#: silica pole 0.77, so 0.05 is at least a factor of 6 from either side. It also covers the
#: pathological region of an alpha=0 CPML, omega*dt below about 0.02 (see
#: boundaries.PML_DISP_ALPHA).
LOSSLESS_LOWFREQ_ADT_MAX = 0.05

def _lossless_lowfreq_damping_rel() -> float:
    return float(knobs.env("LOSSLESS_LOWFREQ_DAMPING_REL"))    # default registered in knobs.KNOBS, equal to LOSSLESS_LOWFREQ_DAMPING_REL

def stabilize_lossless(poles: list, dt: float | None = None) -> list:
    """Inject a tiny artificial loss into lossless poles; lossy poles are returned bitwise unchanged.

    With ``gamma_rel = 0`` and the low-frequency branch off, the whole table is returned bitwise
    unchanged. For the tests and the magnitudes see :data:`LOSSLESS_AXIS_RTOL`,
    :data:`LOSSLESS_DAMPING_REL` and :data:`LOSSLESS_LOWFREQ_DAMPING_REL`. The low-frequency branch
    needs ``dt`` in order to test |a|*dt; without ``dt`` only the uniform relative damping applies,
    which is bitwise identical to the old signature.
    """
    rel = _lossless_damping_rel()
    lf = _lossless_lowfreq_damping_rel() if dt is not None else 0.0
    if rel == 0.0 and lf == 0.0:
        return poles
    out = []
    for q, r in poles:
        q = complex(q)
        if abs(q.real) <= LOSSLESS_AXIS_RTOL * abs(q):
            g = rel
            if lf > 0.0 and abs(q) * dt < LOSSLESS_LOWFREQ_ADT_MAX:
                g = max(g, lf)
            if g > 0.0:
                q = complex(-g * abs(q), q.imag)
        out.append((q, r))
    return out

def trap_coeffs(poles: list, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """``(A-1, B)`` for a batch of poles, before any weight is applied.

    Trapezoidal discretization of ``dP/dt = qP + rE`` under the mapping ``q = a, r = c``; measured
    point by point, eps agrees on both sides to 8.6e-15 relative. ``A-1`` is stored rather than
    ``A``; see ``kernels/dispersion.cu`` for why.
    """
    am1, b = [], []
    for q, r in poles:
        den = 1.0 - q * dt / 2.0
        am1.append((1.0 + q * dt / 2.0) / den - 1.0)
        b.append((r * dt / 2.0) / den)
    return np.array(am1, dtype=np.complex128), np.array(b, dtype=np.complex128)

def eps_rational(eps_inf: float, poles: list, sigma: float = 0.0):
    """``eps(s) = N(u)/D(u)`` with ``u = s/scale``; returns ``(N, D, scale)``, highest power first.

    ``eps(s) = eps_inf + sigma/(eps0*s) + sum_k [c_k/(s-a_k) + conj(c_k)/(s-conj(a_k))]``.

    ``s = -i*omega`` is **the variable in which the verified mapping ``q=a, r=c`` holds**, i.e.
    Tidy3D's formula with ``j = +i`` and ``s = -j*omega``. Choosing the wrong variable moves the
    poles into the right half plane and makes a stable scheme look unstable.
    """
    qs = []
    for a, c in poles:
        a, c = complex(a), complex(c)
        qs += [(a, c), (np.conj(a), np.conj(c))]
    if sigma:
        qs.append((0.0 + 0j, sigma / cpml.EPSILON_0))     # conductivity is a pole at s=0
    # **Non-dimensionalizing is mandatory.** An optical-frequency pole has |q| ~ 1e16, and seven
    # poles give a degree-14 polynomial whose coefficients span (1e16)^14 ~ 1e224. np.roots is
    # completely untrustworthy at that scale; one scene reported "max Re(p) = 1.000e+00", an integer
    # value that is pure numerical noise. Working in u = s/scale makes every coefficient O(1).
    scale = max((abs(q) for q, _ in qs if abs(q) > 0), default=1.0)
    qs = [(q / scale, r / scale) for q, r in qs]
    den = np.array([1.0 + 0j])
    for q, _ in qs:
        den = np.convolve(den, [1.0, -q])
    num = eps_inf * den.copy()
    for m, (_, r) in enumerate(qs):
        rest = np.array([1.0 + 0j])
        for n2, (q2, _) in enumerate(qs):
            if n2 != m:
                rest = np.convolve(rest, [1.0, -q2])
        num = num + r * np.concatenate([np.zeros(len(num) - len(rest)), rest])
    return num, den, scale

#: Threshold below which the real part of a zeta pole counts as numerical noise, in u space where
#: |u| is O(1). Below it, the pole is taken to lie exactly on the imaginary axis.
#: Basis: the backward error of np.roots on a polynomial with O(1) coefficients is at machine
#: precision; the two roots of one lossless fluoride fit measured Re(u) = 0.0 and -8.3e-17, while a
#: genuinely lossy material has |Re(u)| of order at least 1e-2 (|Re(p)| >= 4e13 against
#: scale ~ 1e16). 1e-9 is seven orders of magnitude from either side.
ZETA_AXIS_RTOL = 1e-9

def zeta_poles(eps_inf: float, poles: list, sigma: float = 0.0):
    """Pole expansion of ``1/eps``: ``(zeta_inf, p, R)`` with ``zeta = zeta_inf + sum R_m/(s-p_m)``.

    ``p_m`` are the roots of ``N(s)`` and ``R_m = D(p_m)/N'(p_m)``. Measured across gold, silver,
    silicon, silicon nitride and silica, every ``Re(p_m) < 0``, with a zeta reconstruction error of
    5e-13.

    **The zeta poles of a lossless material lie exactly on the imaginary axis.** With eps(omega)
    real and zero damping, the zeros of eps(s), which are the poles of 1/eps, are the longitudinal
    resonances +-i*omega_L. That is analytic, not a numerical coincidence: one lossless fluoride
    fit has a single Sellmeier-type pole at a = -i*2.5358e16 and zeta poles at +-i*3.4945e16, with
    np.roots returning real parts of 0.0 and -8.3e-17. Roots whose |Re| falls below the noise
    threshold are projected back onto the imaginary axis before being returned: **building the
    table no longer refuses a marginally stable pole**, and trap_coeffs on the eps branch never
    refused an undamped eps pole either, so the two branches share one convention. Whether a
    marginal pole is safe where it is actually used, the D recursion on a slanted interface cell,
    is the caller's responsibility; see dispersion._mix_entries.

    Raises:
        ValueError: a pole lies **strictly** in the right half plane, beyond the noise threshold.
            The D-form recursion would then diverge, so stop rather than continue.
    """
    num, den, scale = eps_rational(eps_inf, poles, sigma)
    pu = np.roots(num)
    on_axis = np.abs(pu.real) <= ZETA_AXIS_RTOL * np.abs(pu)
    pu = np.where(on_axis, 1j * pu.imag, pu)
    dn = np.polyder(num)
    # Converting a residue from u space back to s space multiplies by scale:
    # R/(u-u_m) = (R*scale)/(s-u_m*scale)
    ru = np.array([np.polyval(den, pk) / np.polyval(dn, pk) for pk in pu])
    p, r = pu * scale, ru * scale
    worst = float(np.max(p.real)) if p.size else -1.0
    if worst > 0.0:
        raise ValueError(
            f"1/eps has a pole in the right half plane (max Re = {worst:.3e}); the D-form recursion "
            "would diverge")
    return 1.0 / eps_inf, p, r

def trap_from_poles(p: np.ndarray, r: np.ndarray, dt: float):
    """``(A-1, B)`` from the trapezoidal discretization of ``dQ/dt = pQ + rD``."""
    den = 1.0 - p * dt / 2.0
    return (1.0 + p * dt / 2.0) / den - 1.0, (r * dt / 2.0) / den

def sigma_si(med, pr, freq: float) -> float:
    """Conductivity of a medium, **in SI S/m**.

    ``med.conductivity`` must not be read directly: Tidy3D works in um, so that field is in
    **S/um** and using it as is comes out a factor of 1e6 too small. Measured, a medium with
    eps_r=1.72 has a field value of 8.6119e-06 against a true 8.6119 S/m. This solves back from
    ``Im(eps)`` instead, on the ``media.sigma_from_eps`` conversion already verified against the
    reference outputs.

    A dispersive medium returns 0: all of its loss lives in the poles, and Tidy3D's
    ``PoleResidue`` has no ``conductivity`` field either (it reads back as ``None``). Adding a
    sigma on top would count the loss twice.
    """
    if pr is not None:
        return 0.0
    e = np.asarray(med.eps_model(np.array([freq])), dtype=np.complex128)
    return float(sigma_from_eps(e, freq)[0])

def check_split(mediums, prs, kept, sigma_dc, eps_inf_of, sigma_of, freqs) -> None:
    """Self-check: after removing the DC pole and folding it into sigma, eps(omega) must still equal
    ``eps_model``.

    This is the correctness test for the split: one term fewer in the coefficients, one more in
    sigma, and the two must cancel exactly. A mismatch means the ``sigma = 2*Re(c)*eps0``
    conversion is wrong, so stop here rather than continue.
    """
    w = 2 * np.pi * freqs
    for i, (med, pr) in enumerate(zip(mediums, prs)):
        if pr is None:
            continue
        if type(med).__name__.startswith("Custom"):
            continue  # spatial Custom: scalar reference poles against a spatial eps_model; the fit
                      # supplies the per-cell values
        got = np.full(freqs.size, eps_inf_of[i], dtype=np.complex128)
        for q, r in kept[i]:
            got -= r / (1j * w + q) + np.conj(r) / (1j * w + np.conj(q))
        got += 1j * sigma_of[i] / (w * cpml.EPSILON_0)
        ref = np.asarray(med.eps_model(freqs), dtype=np.complex128).ravel()
        rel = float(np.max(np.abs(got - ref) / np.maximum(np.abs(ref), 1.0)))
        if rel > 1e-10:
            raise ValueError(
                f"medium {i} ({type(med).__name__}): after removing the DC pole, eps no longer "
                f"matches eps_model; max relative difference {rel:.3e}. DC poles removed: "
                f"{len(pr.poles) - len(kept[i])}, folded sigma = {sigma_dc[i]:.6g} S/m")
