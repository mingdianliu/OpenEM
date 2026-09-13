# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""The two dataclasses that carry a solve's results.

These used to live in solver.py. They were split out because readout.py needs to construct
``Phasors``, and leaving them in solver.py made the two modules import each other in a cycle.
solver.py still re-exports both names.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Phasors:
    """Phasors accumulated over one monitor plane.

    ``data`` is complex with shape ``(4, nf, n1, n2)``, where ``(n1, n2)`` counts cells along the
    two transverse axes, the ones other than the monitor's normal, in ascending order ``(t1, t2)``.
    Component order is E_t1, E_t2, H_t1, H_t2, with H already colocated along the normal onto the
    plane where E sits (``accumulate_dft`` in kernels/source_dft.cu).
    """

    name: str
    data: np.ndarray
    freqs: np.ndarray


def result_meta(res) -> dict:
    """Summary written into ``ours['__meta']``: step count, early termination, energy trace.

    Both result-packing paths use this (autograd_hook._result_to_ours and
    nb/solve_worker.save_result), and the divergence guard on the assembling side reads only this.
    """
    return {
        "steps_run": int(getattr(res, "steps_run", 0)),
        "triggered": bool(getattr(res, "shutoff_triggered", False)),
        "stop": getattr(res, "stop_reason", "") or "",
        "wall": float(getattr(res, "wall_seconds", 0.0)),
        "source_end": int(getattr(res, "source_end", -1)),
        "energy_trace": [[int(s), float(r)]
                         for s, r in (getattr(res, "energy_trace", None) or [])],
    }


@dataclass
class Result:
    phasors: dict[str, Phasors]
    steps_run: int
    shutoff_triggered: bool
    decay_trace: list[tuple[int, float]] = field(default_factory=list)
    #: Sampled total energy as ``(step, energy)`` (shutoff.EnergyDecay.abs_history, one point every
    #: interval steps; empty when the step count is pinned or early termination is off). The
    #: assembling side uses it to catch divergence that is finite but growing exponentially. One
    #: optimization iteration ran the full 9,205 steps with |amp| only reaching 3.9e3, invisible to
    #: both the NaN guard and the 1e6 ceiling, yet it had already driven the objective to -69 dB
    #: (nb.backend._refuse_energy_growth).
    energy_trace: list[tuple[int, float]] = field(default_factory=list)
    #: Step at which every source has finished injecting (shutoff.source_end_step); -1 if unknown
    source_end: int = -1
    wall_seconds: float = 0.0
    #: Which criterion triggered early termination; empty string if none did
    stop_reason: str = ""
    #: Final state of the early-termination criteria, for reporting
    shutoff_summary: str = ""
    #: Final-state fields ``{"Ex": (nx,ny,nz), ...}``, filled only when ``return_fields=True``.
    fields: dict[str, np.ndarray] = field(default_factory=dict)
    #: FieldMonitor phasors ``{name: (6, nf, ni, nj, nk)}``, **at their native Yee positions, not
    #: colocated**. Component order Ex, Ey, Ez, Hx, Hy, Hz; for the index box see
    #: ``Scene.field_monitors``.
    field_phasors: dict[str, np.ndarray] = field(default_factory=dict)
    #: FieldTimeMonitor time samples ``{name: (nc, nslots, ni, nj, nk)}``, at native Yee positions.
    #: The components are whichever ``FieldTimeMonitor.comps`` lists, not all six.
    time_samples: dict[str, np.ndarray] = field(default_factory=dict)
    #: How many slots each time monitor **actually** filled. Only slots where both H half-updates
    #: completed are counted: if early termination lands between the two, the last slot holds only
    #: half of H, and dropping it is the safer choice.
    time_slots: dict[str, int] = field(default_factory=dict)
    #: FluxTimeMonitor instantaneous flux series ``{name: (n_filled,)}`` as float64, one scalar per
    #: slot, at the same times ``m*dt`` as :attr:`time_samples` (H already time-averaged).
    flux_time: dict[str, np.ndarray] = field(default_factory=dict)
