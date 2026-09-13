# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""The one entry point to the CUDA side.

**This module is the only place that imports cupy.** Switching to ctypes one day, or moving the
main loop into C++, changes only this one file; ``kernels/*.cu`` is real CUDA C written as
``extern "C" __global__``, and the very same source can be compiled straight into a standalone
C++ main without rewriting.

Kernel-level timing (``Profiler``) lives here too, precisely so that invariant is not broken.
"""

from __future__ import annotations

import pathlib
from collections import Counter, defaultdict
from typing import Callable

import cupy as cp

from openem import knobs

KERNEL_DIR = pathlib.Path(__file__).parent / "kernels"

#: NVRTC compile options. sm_90 = H800.
NVRTC_OPTIONS = ("-std=c++17", "--use_fast_math")


class Profiler:
    """Kernel-level timing and launch counting (for profiling the GPU acceleration).

    With ``budget=0``, ``Kernels.__getitem__`` returns the ``RawKernel`` directly at zero cost.
    That is the default state, and the normal solve path runs not one extra line of code.

    After ``arm(n)``, a CUDA event is recorded before and after each of the first n launches
    and the results are merged by kernel name. The events only record, they never synchronize,
    so the pipeline is not interrupted; the gap between two consecutive launches is counted
    separately, and that gap is the launch overhead of Python plus the driver. The event pool
    is allocated up front so the cost of creating an event does not land inside that gap.

    ``calls`` counts the whole run (n does not limit it) and is what gives "launches per step".
    """

    def __init__(self, budget: int = 0) -> None:
        self.budget = budget
        self.steps = 0          #: only affects the normalization in the printout
        self.cells = 0          #: only affects the normalization in the printout
        self.calls: Counter[str] = Counter()
        self._recs: list[tuple[str, cp.cuda.Event, cp.cuda.Event]] = []
        self._pool: list[tuple[cp.cuda.Event, cp.cuda.Event]] = []
        self._i = 0
        self._reported = False

    @property
    def enabled(self) -> bool:
        return self.budget > 0

    def arm(self, budget: int) -> None:
        """Clear the old records and preallocate the event pool. ``budget=1`` only turns
        counting on."""
        self.budget = budget
        self.calls.clear()
        self._recs.clear()
        self._pool = [(cp.cuda.Event(), cp.cuda.Event()) for _ in range(max(budget, 0))]
        self._i = 0
        self._reported = False

    def wrap(self, name: str, fn) -> Callable:
        """Wrap a RawKernel into a callable that counts launches and times a limited number."""

        def launch(grid, block, args, **kw) -> None:
            self.calls[name] += 1
            if self._i < self.budget:
                beg, end = self._pool[self._i]
                self._i += 1
                beg.record()
                fn(grid, block, args, **kw)
                end.record()
                self._recs.append((name, beg, end))
            else:
                fn(grid, block, args, **kw)

        return launch

    def report(self) -> str:
        """Merge and format. **A truncated run warns loudly**, since span/steps is unreliable
        then."""
        self._reported = True
        if not self._recs:
            return "(the profiler recorded no launches)"
        cp.cuda.runtime.deviceSynchronize()
        ker: dict[str, float] = defaultdict(float)
        cnt: Counter[str] = Counter()
        gap = 0.0
        prev = None
        for name, beg, end in self._recs:
            ker[name] += cp.cuda.get_elapsed_time(beg, end)
            cnt[name] += 1
            if prev is not None:
                gap += cp.cuda.get_elapsed_time(prev, beg)
            prev = end
        span = cp.cuda.get_elapsed_time(self._recs[0][1], prev)

        total_calls = sum(self.calls.values())
        timed = len(self._recs)
        out = [
            "=" * 92,
            f"KERNEL PROFILE  timed {timed:,d} / {total_calls:,d} launches in total"
            + (f"  covering {self.steps:,d} steps" if self.steps else ""),
        ]
        if timed < total_calls:
            out.append(
                f"!! the event pool only covers {timed} launches, the run has {total_calls}: "
                "span/steps is unreliable, raise --budget and run it again"
            )
        out.append(
            f"{'kernel':<26}{'timed':>9}{'total ms':>11}{'share':>8}"
            f"{'mean us':>10}{'all calls':>11}{'/step':>8}"
        )
        for name in sorted(ker, key=lambda x: -ker[x]):
            tot = ker[name]
            per_step = f"{self.calls[name] / self.steps:>8.2f}" if self.steps else " " * 8
            out.append(
                f"{name:<26}{cnt[name]:>9d}{tot:>11.2f}{tot / span * 100:>7.1f}%"
                f"{tot / cnt[name] * 1e3:>10.1f}{self.calls[name]:>11d}{per_step}"
            )
        out.append(f"{'[launch gaps]':<26}{'':>9}{gap:>11.2f}{gap / span * 100:>7.1f}%")
        out.append(f"{'[span total]':<26}{'':>9}{span:>11.2f}{100.0:>7.1f}%")
        if self.steps:
            per = span / self.steps
            out.append(f"per step {per:.4f} ms (span / steps)")
            if self.cells:
                out.append(
                    f"throughput {self.cells / per / 1e6:.3f} Gcell/s"
                    f" ({self.cells:,d} cells)"
                )
            out.append(f"per step {total_calls / self.steps:.2f} launches")
        return "\n".join(out)


#: The global profiler. ``OPENEM_PROFILE=<n>`` turns it on straight from the environment.
PROF = Profiler(int(knobs.env("PROFILE")))


class Kernels:
    """Load and cache every .cu under kernels/.

    When ``PROF`` is on (``OPENEM_PROFILE=<n>``) this returns a callable wrapped with timing,
    otherwise the RawKernel itself, so with it off there is no overhead.
    """

    def __init__(self, options: tuple[str, ...] = NVRTC_OPTIONS) -> None:
        self._fn: dict[str, cp.RawKernel] = {}
        for path in sorted(KERNEL_DIR.glob("*.cu")):
            code = path.read_text(encoding="utf-8")
            mod = cp.RawModule(code=code, options=options, backend="nvrtc",
                               name_expressions=None)
            for name in _kernel_names(code):
                if name in self._fn:
                    raise ValueError(f"duplicate kernel name: {name} ({path.name})")
                self._fn[name] = mod.get_function(name)

    def __getitem__(self, name: str):
        if name not in self._fn:
            raise KeyError(f"no kernel {name!r}; the loaded ones are {sorted(self._fn)}")
        fn = self._fn[name]
        return PROF.wrap(name, fn) if PROF.enabled else fn

    @property
    def names(self) -> list[str]:
        return sorted(self._fn)


def _kernel_names(code: str) -> list[str]:
    """Pick ``extern "C" __global__ [attributes] void <name>(`` out of the source.

    Deliberately not a general C parser: we control the format of the kernel declarations
    ourselves, so a simple scan is enough, and if the format ever changes the scan stops
    matching and the problem shows up at once, which beats silently missing a kernel.
    """
    out = []
    for line in code.splitlines():
        s = line.strip()
        if s.startswith("extern \"C\" __global__ ") and "void " in s:
            # __launch_bounds__(...) may sit between __global__ and void
            name = s.split("void ", 1)[1].split("(", 1)[0].strip()
            if name:
                out.append(name)
    return out


def grid_1d(n: int, block: int = 256) -> tuple:
    """Launch configuration of the 1D table kernels (one thread per entry)."""
    return ((n + block - 1) // block, 1, 1), (block, 1, 1)


def grid_2d(nx: int, ny: int, block: tuple[int, int] = (32, 8)) -> tuple:
    """Launch configuration of the planar kernels. The first dimension maps j (the contiguous
    one), the second maps i."""
    return ((ny + block[0] - 1) // block[0], (nx + block[1] - 1) // block[1]), block


#: Default thread block shape of the bulk kernels. ``OPENEM_BLOCK3="32x8"`` overrides it: the
#: trade-off between occupancy and redundancy can only be measured (512 threads measured 1.7x
#: faster than 128), so there is a knob to sweep it without editing code. The default value is
#: unchanged.
BLOCK3 = tuple(int(v) for v in knobs.env("BLOCK3").split("x"))


def grid_3d(nx: int, ny: int, nz: int, block: tuple[int, int] | None = None) -> tuple:
    """Launch configuration of the bulk kernels.

    The thread x dimension maps **k** (the contiguous dimension in memory) so that accesses
    coalesce; blockIdx.z maps i.
    """
    bx, by = block or BLOCK3
    if block is None and nz < bx:
        # Thin z domains and thin boxes (MoS2 has nz=3 over the whole domain, a planar
        # monitor box has nk=1): the default block lays bx lanes along z, of which
        # (bx-nz)/bx would sit idle (61/64 on MoS2), so it switches to the smallest power of
        # two that still covers nz and gives the lanes to y. k is the contiguous dimension in
        # memory and neighboring (j,k) rows are contiguous anyway, so coalescing is
        # unaffected; the block shape is only a launch configuration and the result stays
        # bitwise the same.
        bx = max(2, 1 << (int(nz) - 1).bit_length())
        by = max((BLOCK3[0] * BLOCK3[1]) // bx, 1)
    return (
        ((nz + bx - 1) // bx, (ny + by - 1) // by, nx),
        (bx, by, 1),
    )
