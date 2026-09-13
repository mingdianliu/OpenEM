# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Solve a serialized scene (``scene.npz``) on the local GPU and store every monitor output as
``ours.npz``.

In-process: ``solve_dir(work, out)``. It also runs as a subprocess,
``python -m openem.nb.solve_worker <work> <out>``. Set ``OPENEM_SUBPROCESS=1`` when a long notebook
does dozens of solves and you want device memory isolated between them; the integration layer then
takes the subprocess route on its own.
"""
from __future__ import annotations

import json
import os
import pathlib
import subprocess
import sys

import numpy as np

from openem import serialize, solver
from openem.device import Kernels


def save_result(res, path: pathlib.Path, extra: dict | None = None) -> float:
    """``solver.Result`` to an npz.

    Keys: ``mon:<name>`` phasors, ``monf:<name>`` frequencies, and ``fld:``, ``tim:``, ``ftm:`` for
    the field and time-domain data.
    """
    arrs = dict(extra or {})
    for n, p in res.phasors.items():
        arrs[f"mon:{n}"] = p.data
        arrs[f"monf:{n}"] = np.asarray(p.freqs)
    for attr, tag in (("field_phasors", "fld"), ("time_samples", "tim"), ("flux_time", "ftm")):
        for n, a in (getattr(res, attr, {}) or {}).items():
            arrs[f"{tag}:{n}"] = a
    from openem.results import result_meta
    meta = result_meta(res)
    total = sum(a.nbytes for a in arrs.values())
    if total > float(os.environ.get("OPENEM_OUT_MAX_GB", "20")) * 1e9:
        path.with_suffix(".TOOBIG").write_text(json.dumps(meta))
        return -total / 1e6
    np.savez(path, __meta=json.dumps(meta), **arrs)
    return total / 1e6


def solve_dir(work: pathlib.Path, out: pathlib.Path, steps: int | None = None, subprocess_: bool = False) -> pathlib.Path:
    """Solve ``work/scene.npz``, write ``out/ours.npz``, and return that path."""
    work, out = pathlib.Path(work), pathlib.Path(out)
    if subprocess_:
        cmd = [sys.executable, "-u", "-m", "openem.nb.solve_worker", str(work), str(out)]
        if steps is not None:
            cmd += ["--steps", str(steps)]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=float(os.environ.get("OPENEM_SOLVE_TIMEOUT", "86400")))
        if r.returncode != 0 or not (out / "ours.npz").exists():
            raise RuntimeError(f"subprocess solve failed:\n{(r.stdout or '')[-1500:]}\n{(r.stderr or '')[-800:]}")
        return out / "ours.npz"
    out.mkdir(parents=True, exist_ok=True)
    sc = serialize.load(work / "scene.npz")
    res = solver.run(sc, num_steps=steps, use_shutoff=steps is None, kernels=Kernels(), verbose=True)
    extra = {}
    plan_f = work / "td_coords.npz"          # optional: pre-generated readout coordinates in Tidy3D convention
    if plan_f.exists() and res.field_phasors:
        from openem import td_readout
        extra = td_readout.finalize(sc, res.field_phasors, dict(np.load(plan_f)))
    mb = save_result(res, out / "ours.npz", extra=extra)
    print(f"SOLVE|{work.name}|{res.steps_run}|{int(res.shutoff_triggered)}|{res.stop_reason or '-'}|"
          f"{int(np.prod(sc.shape))}|{res.wall_seconds:.1f}|{mb:.1f}", flush=True)
    return out / "ours.npz"


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work"); ap.add_argument("out"); ap.add_argument("--steps", type=int, default=None)
    a = ap.parse_args(argv)
    solve_dir(a.work, a.out, steps=a.steps)
    return 0


if __name__ == "__main__":
    sys.exit(main())
