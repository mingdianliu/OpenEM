# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
# -*- coding: utf-8 -*-
"""A tiny artificial loss for lossless poles (poles.stabilize_lossless).

Background: the poles of a lossless dispersive material sit on the real frequency axis, and after
trapezoidal discretization |A| is identically 1, which is marginally stable. The recursion on its
own stays bounded over a long fp32 run (test_zeta_stability), but coupled to the CPML over a long
time it diverges (Autograd20: NaN at 62,028 steps; the same structure reproduced without a CPML
does not diverge). The fix: while building the tables, inject an artificial loss
γ = LOSSLESS_DAMPING_REL·|a| into the lossless poles, which pushes |A| inside the unit circle.

Four criteria:

1. |A| of a lossless pole is strictly < 1, and the shrinkage is of the same order as γ·dt;
2. lossy poles are bitwise unchanged (damping off vs damping on);
3. the perturbation of ε(ω) over the working band is <= 1e-5 relative, far below any observable
   tolerance;
4. an environment variable can turn it off (γ=0 gives bitwise the old behaviour).

Pure CPU, no GPU needed.
"""
from __future__ import annotations

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from openem.scene import poles as poles_mod
from openem.scene.poles import (
    LOSSLESS_AXIS_RTOL, LOSSLESS_DAMPING_REL, split_dc_pole, stabilize_lossless,
    trap_coeffs,
)

#: The time step measured on Autograd20 (the 62,028-step run); used for the order-of-magnitude
#: assertions.
DT = 4.8366e-17


def _kept(med, dt=DT):
    k, _ = split_dc_pole(med.pole_residue, dt)
    return k


@pytest.mark.parametrize("mat,name", [
    ("cSi", "Palik_Lossless"), ("SiO2", "Palik_Lossless")])
def test_lossless_poles_pushed_inside_unit_circle(mat, name):
    """After stabilization |A| is strictly < 1, with a shrinkage of order γ_rel·|a|·dt."""
    kept = stabilize_lossless(_kept(td.material_library[mat][name]))
    am1, _ = trap_coeffs(kept, DT)
    for (q, _r), a1 in zip(kept, am1):
        absA = abs(1.0 + complex(a1))
        gap_want = LOSSLESS_DAMPING_REL * abs(complex(q)) * DT
        assert absA < 1.0, f"{mat}: |A|={absA!r} did not get inside the unit circle"
        assert 0.2 * gap_want < 1.0 - absA < 5.0 * gap_want, (
            f"{mat}: 1-|A|={1.0 - absA:.3e}, expected ~{gap_want:.3e}")


def test_lossy_poles_bitwise_unchanged(monkeypatch):
    """The poles of a genuinely lossy material (gold) are bitwise unchanged."""
    kept = _kept(td.material_library["Au"]["JohnsonChristy1972"])
    assert all(abs(q.real) > LOSSLESS_AXIS_RTOL * abs(q) for q, _ in kept), (
        "precondition fails: every pole of gold should be lossy")
    stab = stabilize_lossless(kept)
    for (q0, r0), (q1, r1) in zip(kept, stab):
        assert complex(q0) == complex(q1) and complex(r0) == complex(r1)
    # With damping off the whole table comes back as it was (same object, zero-overhead path)
    monkeypatch.setenv("OPENEM_LOSSLESS_DAMPING_REL", "0")
    assert stabilize_lossless(kept) is kept


def test_eps_perturbation_within_tolerance():
    """The perturbation γ causes in ε(ω) is <= 1e-5 relative over the working band."""
    freqs = td.C_0 / np.linspace(1.5, 1.6, 31)
    w = 2 * np.pi * freqs
    for mat in ("cSi", "SiO2"):
        med = td.material_library[mat]["Palik_Lossless"]
        kept = _kept(med)
        stab = stabilize_lossless(kept)

        def eps(pl):
            e = np.zeros(w.size, dtype=np.complex128)
            for q, r in pl:
                e -= r / (1j * w + q) + np.conj(r) / (1j * w + np.conj(q))
            return e

        ref = np.asarray(med.eps_model(freqs), dtype=np.complex128).ravel()
        rel = np.max(np.abs(eps(stab) - eps(kept)) / np.abs(ref))
        assert rel < 1e-5, f"{mat}: Δε/ε = {rel:.3e}"


def test_env_override_restores_marginal(monkeypatch):
    """OPENEM_LOSSLESS_DAMPING_REL=0 reproduces the old, marginal |A|=1 behaviour bitwise."""
    kept = _kept(td.material_library["cSi"]["Palik_Lossless"])
    monkeypatch.setenv("OPENEM_LOSSLESS_DAMPING_REL", "0")
    am1_off, b_off = trap_coeffs(stabilize_lossless(kept), DT)
    am1_ref, b_ref = trap_coeffs(kept, DT)
    assert np.array_equal(am1_off, am1_ref) and np.array_equal(b_off, b_ref)
    assert abs(abs(1.0 + complex(am1_off[0])) - 1.0) < 1e-12
    # Override with 10x: the shrinkage scales proportionally
    monkeypatch.setenv("OPENEM_LOSSLESS_DAMPING_REL",
                       str(10 * LOSSLESS_DAMPING_REL))
    am1_10, _ = trap_coeffs(stabilize_lossless(kept), DT)
    monkeypatch.delenv("OPENEM_LOSSLESS_DAMPING_REL")
    am1_1, _ = trap_coeffs(stabilize_lossless(kept), DT)
    gap10 = 1.0 - abs(1.0 + complex(am1_10[0]))
    gap1 = 1.0 - abs(1.0 + complex(am1_1[0]))
    assert 8.0 < gap10 / gap1 < 12.0


def test_material_tables_zeta_strictly_left_half_plane():
    """Through the full _material_tables: the ζ poles of a lossless material land strictly in the
    left half plane."""
    from openem.scene import dispersion as disp_mod

    sim = td.Simulation(
        size=(0.5, 0.5, 0.5), grid_spec=td.GridSpec.uniform(dl=0.05),
        run_time=1e-13, medium=td.material_library["SiO2"]["Palik_Lossless"],
        boundary_spec=td.BoundarySpec.all_sides(td.Periodic()),
        sources=[td.PointDipole(
            center=(0, 0, 0), polarization="Ez",
            source_time=td.GaussianPulse(freq0=2e14, fwidth=2e13))])
    mediums = [sim.medium, td.material_library["cSi"]["Palik_Lossless"]]
    freqs = np.linspace(1.8e14, 2.1e14, 4)
    t = disp_mod._material_tables(sim, mediums, freqs)
    for zinf, p, r in t["zeta_of"]:
        if p.size:
            assert float(np.max(p.real)) < 0.0, "ζ poles must lie strictly in the left half plane"
    for _i, am1, _b in t["disp"]:
        for a1 in am1:
            assert abs(1.0 + complex(a1)) < 1.0, "|A| of the ε branch must be strictly < 1"
