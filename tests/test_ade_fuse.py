# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
# -*- coding: utf-8 -*-
"""P43 (folding dispersion into update_e) must be bitwise identical to the two-pass form,
update_e followed by dispersion_step.

It is off by default, measured at 0.82 to 0.90 times the speed. The code stays in the tree and
this test guards it.
"""
from __future__ import annotations

import os

import numpy as np
import pytest

cp = pytest.importorskip("cupy", reason="fusion is a GPU path")
from openem import serialize, solver           # noqa: E402
from openem.device import Kernels              # noqa: E402

ROOT = os.environ.get("OPENEM_FIXTURES", "/nonexistent")  # external scene fixtures; skip if unset
CASES = ("PlasmonicYagiUda", "MIMResonator")


def _out(case: str, fuse: str) -> dict[str, np.ndarray]:
    os.environ["OPENEM_ADE_FUSE"] = fuse
    try:
        sc = serialize.load(f"{ROOT}/{case}/scene.npz")
        res = solver.run(sc, num_steps=20, use_shutoff=False,
                         kernels=Kernels(), verbose=False)
    finally:
        os.environ.pop("OPENEM_ADE_FUSE", None)
    out: dict[str, np.ndarray] = {}
    for nm, ph in res.field_phasors.items():
        out[f"field_{nm}"] = np.asarray(ph)
    for nm, pz in res.phasors.items():
        out[f"flux_{nm}"] = np.asarray(pz.data)
    return out


@pytest.mark.parametrize("case", CASES)
def test_ade_fuse_bitwise(case: str) -> None:
    if not os.path.isdir(f"{ROOT}/{case}"):
        pytest.skip(f"no exported case {case}")
    a, b = _out(case, "0"), _out(case, "1")
    assert sorted(a) == sorted(b) and a
    for nm in sorted(a):
        assert a[nm].tobytes() == b[nm].tobytes(), (
            f"{case}, array {nm}: folded and two-pass are not bitwise identical, "
            f"max|delta|={np.abs(a[nm] - b[nm]).max():.3e}")
