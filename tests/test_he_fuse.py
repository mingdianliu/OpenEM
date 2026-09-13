# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
# -*- coding: utf-8 -*-
"""P41 (fusing the H and E updates, blocked along i) must be bitwise identical to launching them
separately.

It is off by default, measured about 20% slower, but the code stays in the tree and this test keeps
it from rotting. On different hardware, or with a different fusion idea, it can be switched on and
measured directly without re-establishing correctness first.
"""
from __future__ import annotations

import os

import numpy as np
import pytest

cp = pytest.importorskip("cupy", reason="fusion is a GPU path")
from openem import serialize, solver           # noqa: E402
from openem.device import Kernels              # noqa: E402

CASES = ("StartHere", "BoundaryConditions_ft")
ROOT = os.environ.get("OPENEM_FIXTURES", "/nonexistent")  # external scene fixtures; skip if unset


def _outputs(case: str, fuse: str) -> dict[str, np.ndarray]:
    """Run once and pull out every monitor output."""
    os.environ["OPENEM_HE_FUSE"] = fuse
    try:
        sc = serialize.load(f"{ROOT}/{case}/scene.npz")
        res = solver.run(sc, num_steps=25, use_shutoff=False,
                         kernels=Kernels(), verbose=False)
    finally:
        os.environ.pop("OPENEM_HE_FUSE", None)
    out: dict[str, np.ndarray] = {}
    for nm, ph in res.field_phasors.items():
        out[f"fieldmon_{nm}"] = np.asarray(ph)
    for nm, pz in res.phasors.items():
        out[f"fluxmon_{nm}"] = np.asarray(pz.data)
    for nm, v in res.time_samples.items():
        out[f"timemon_{nm}"] = np.asarray(v)
    return out


@pytest.mark.parametrize("case", CASES)
def test_he_fuse_bitwise(case: str) -> None:
    if not os.path.isdir(f"{ROOT}/{case}"):
        pytest.skip(f"no exported case {case}")
    a = _outputs(case, "0")
    b = _outputs(case, "1")
    assert sorted(a) == sorted(b) and a, f"{case}: output names do not line up"
    for nm in sorted(a):
        x, y = a[nm], b[nm]
        assert x.shape == y.shape and x.dtype == y.dtype
        assert x.tobytes() == y.tobytes(), (
            f"{case}, array {nm}: fused and separate launches are not bitwise identical, "
            f"max|delta|={np.abs(x - y).max():.3e}")
