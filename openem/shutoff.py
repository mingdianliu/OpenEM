# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Early shutoff tests. Three quantities; the run stops as soon as any one of them is met.

- :class:`FieldDecay`        instantaneous Σ V·|E|², Tidy3D's convention, used to reproduce its
                             stopping step
- :class:`EnergyDecay`       total energy including H; its envelope does not oscillate at 2ω,
                             and the decay constant tau is fitted from it
- :class:`PhasorConvergence` extrapolates "how much is still going to change" from tau and
                             stops once that is below the target tolerance
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class FieldDecay:
    """Instantaneous ``Σ V_cell·|E|²`` over its **historical** maximum (docs/architecture.md
    §4.1).

    Four things that are easy to get wrong: use |E|², not |E|; sum over the whole domain
    instead of taking the largest cell; take the maximum over time in the denominator; include
    the PML. No ε factor, and H is not included.
    """

    threshold: float
    peak: float = 0.0
    ratio: float = float("nan")
    history: list[tuple[int, float]] = field(default_factory=list)

    def update(self, step: int, metric: float) -> bool:
        self.peak = max(self.peak, metric)
        self.ratio = metric / self.peak if self.peak > 0 else 0.0
        self.history.append((step, self.ratio))
        return self.peak > 0 and self.ratio < self.threshold

    @property
    def deepest(self) -> tuple[int, float]:
        """The deepest check, ``(step, ratio)``, used to say how far it was from firing."""
        if not self.history:
            return (0, float("nan"))
        i = int(np.argmin([r for _, r in self.history]))
        return self.history[i]


@dataclass
class EnergyDecay:
    """Total electromagnetic energy ``½(ε₀ε|E|² + μ₀|H|²)`` over its historical maximum.

    Because H is included, it does not oscillate at 2ω on a resonant mode: it is an envelope
    and does not depend on the sampling phase. The quantity comes from FDTDX's
    ``EnergyThresholdCondition``, and taking the ratio comes from Tidy3D.
    """

    threshold: float
    peak: float = 0.0
    ratio: float = float("nan")
    history: list[tuple[int, float]] = field(default_factory=list)
    #: Absolute total energy ``(step, energy)``, which the divergence guard on the assembly side
    #: needs (the ratio loses the fact that it grew back after the source stopped)
    abs_history: list[tuple[int, float]] = field(default_factory=list)

    def update(self, step: int, energy: float) -> bool:
        self.abs_history.append((step, float(energy)))
        self.peak = max(self.peak, energy)
        self.ratio = energy / self.peak if self.peak > 0 else 0.0
        self.history.append((step, self.ratio))
        return self.peak > 0 and self.ratio < self.threshold

    @property
    def tau(self) -> float:
        """Decay time constant (in steps): a linear fit of ``ln(ratio)`` over the second half.

        Returns ``inf`` when the decay is extremely slow or the energy is still growing, which
        means "nowhere near converged".
        """
        h = self.history
        if len(h) < 8:
            return float("inf")
        arr = np.asarray(h[len(h) // 2 :], dtype=np.float64)
        s, r = arr[:, 0], arr[:, 1]
        good = r > 0
        if good.sum() < 4:
            return float("inf")
        slope = np.polyfit(s[good], np.log(r[good]), 1)[0]
        if slope >= 0 or not np.isfinite(slope):
            return float("inf")
        return float(-1.0 / slope)


@dataclass
class PhasorConvergence:
    """Stop once the flux estimate stops moving: watch the observable directly, not the field
    decay.

    What is tested is the **extrapolated remaining change**, not the change per window::

        remaining ≈ relative change per window * tau / window length

    The per-window change cannot be used directly: it accumulates coherently, and a few hundred
    windows add up to far more than a single window's worth (measured). At high Q, tau is
    enormous -> the extrapolated value is large -> it correctly reports that the run has not
    converged.

    Attributes:
        tol: threshold on the remaining relative change; defaults to our target tolerance 1e-3.
        patience: how many times in a row it has to hold, to stop a false trigger at an
            inflection point.
        window_steps: the check interval, needed by the extrapolation.
    """

    tol: float = 1e-3
    patience: int = 3
    #: Length of the check window (steps); the extrapolation uses it to turn "change per window"
    #: into "total remaining change"
    window_steps: int = 400
    _prev: np.ndarray | None = None
    streak: int = 0
    #: The extrapolated **remaining** relative change; this is what the test looks at
    rel: float = float("inf")
    #: Relative change within the most recent window, for diagnostics
    per_window: float = float("inf")
    history: list[tuple[int, float]] = field(default_factory=list)

    def update(self, step: int, flux: np.ndarray, decay_tau: float | None = None) -> bool:
        """``decay_tau`` comes from :class:`EnergyDecay.tau`; given one, the test extrapolates,
        otherwise it falls back to the per-window change."""
        cur = np.asarray(flux, dtype=np.float64)
        if self._prev is None:
            self._prev = cur
            return False
        scale = float(np.max(np.abs(cur)))
        per_window = (
            float(np.max(np.abs(cur - self._prev)) / scale) if scale > 0 else float("inf")
        )
        self._prev = cur
        self.per_window = per_window
        # Extrapolation: remaining contribution ~ change per window * (tau / window length).
        # When tau is enormous (a high-Q resonance) the extrapolated value is enormous too,
        # which correctly says that this run has not actually converged.
        if decay_tau is not None and np.isfinite(decay_tau) and decay_tau > 0:
            self.rel = per_window * decay_tau / max(self.window_steps, 1)
        else:
            self.rel = per_window
        self.history.append((step, self.rel))
        self.streak = self.streak + 1 if self.rel < self.tol else 0
        return self.streak >= self.patience


@dataclass
class Policy:
    """The three tests combined: any one of them firing stops the run, and ``reason`` records
    which one did.

    Each can be turned off by passing ``None``. ``interval`` is the check interval, which is
    400 steps in Tidy3D.
    """

    field_decay: FieldDecay | None = None
    energy_decay: EnergyDecay | None = None
    phasor: PhasorConvergence | None = None
    interval: int = 400
    reason: str = ""

    def check(
        self, step: int, e2_weighted: float, energy: float, flux: np.ndarray
    ) -> bool:
        """Return whether the run should stop. ``reason`` records which test fired.

        With an empty ``flux`` array the phasor test is skipped: a scene that only has
        FieldMonitors (the PointDipole group, for instance) has no flux to watch and relies on
        the field and energy tests.
        """
        fired = []
        if self.field_decay is not None and self.field_decay.update(step, e2_weighted):
            fired.append(f"field_decay<{self.field_decay.threshold:g}")
        if self.energy_decay is not None and self.energy_decay.update(step, energy):
            fired.append(f"energy_decay<{self.energy_decay.threshold:g}")
        tau = self.energy_decay.tau if self.energy_decay is not None else None
        if (self.phasor is not None and np.size(flux) > 0
                and self.phasor.update(step, flux, decay_tau=tau)):
            fired.append(
                f"phasor_converged(remaining≈{self.phasor.rel:.2e}<{self.phasor.tol:g}"
                f" x{self.phasor.patience})"
            )
        if fired:
            self.reason = " + ".join(fired)
            return True
        return False

    def summary(self) -> str:
        parts = []
        if self.field_decay is not None:
            s, r = self.field_decay.deepest
            parts.append(
                f"field_decay now {self.field_decay.ratio:.3e}, deepest {r:.3e}@{s}"
            )
        if self.energy_decay is not None:
            parts.append(f"energy_decay now {self.energy_decay.ratio:.3e}")
        if self.phasor is not None:
            parts.append(
                f"phasor per window {self.phasor.per_window:.3e}, extrapolated remaining "
                f"{self.phasor.rel:.3e}"
            )
        return "; ".join(parts)


def default_policy(tidy3d_shutoff: float, interval: int = 400) -> Policy:
    """Production default: all three on, the thresholds taken from Tidy3D's shutoff and the
    phasor tolerance set to 1e-3."""
    return Policy(
        field_decay=FieldDecay(threshold=tidy3d_shutoff),
        energy_decay=EnergyDecay(threshold=tidy3d_shutoff),
        phasor=PhasorConvergence(tol=1e-3, patience=3, window_steps=interval),
        interval=interval,
    )


#: Once the source envelope drops below this fraction of its peak, injection counts as over.
SOURCE_TAIL = 1e-3


def source_end_step(sc) -> int:
    """The step by which every source in the scene has finished injecting. **Early shutoff must
    never fire before this.**

    The test "the field has decayed to 1e-6 of its peak" is **meaningless** while a source is
    still injecting: the source has a small but non-zero turn-on value at t=0, and that
    broadband transient shoots up a peak of its own and is then absorbed by the PML while the
    real pulse has not arrived yet, so the ratio slides downhill all through the ramp-up.
    Measured on DielectricMetasurfaceAbsorber, the ``Σ|E|²`` ratio was down to 0.31 at step
    5200 while the peak of the source envelope is at step 83,400.

    The test uses the **waveform tables that have already been sampled** and does not re-parse
    the source parameters: those tables are the very ones used for the injection, so they
    cannot disagree with what is actually injected.
    """
    end = 0
    waves = ([s.waveform for s in sc.sources] + [d.waveform for d in sc.dipoles]
             + [t.waveform for t in sc.tfsf_sources])
    for w in waves:
        a = np.abs(np.asarray(w.amp_int, dtype=np.float64))
        peak = float(a.max())
        if peak <= 0:
            continue
        live = np.flatnonzero(a > SOURCE_TAIL * peak)
        if live.size:
            end = max(end, int(live[-1]))
    return end
