# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""The knobs on the scene and nb sides are registered in knobs.KNOBS.

Their defaults must equal the constants that carry the physical explanation in each module. The
type conversion at each call site must match the original ``float(os.environ.get(...))`` and must
be read live on every call. A knob registered as ``None`` branches on whether it is set at all.
"""
import os
import openem.knobs as knobs
from openem.scene import boundaries, modes, poles


def test_defaults_are_str_or_none():
    for name, (default, kind, _) in knobs.KNOBS.items():
        assert default is None or isinstance(default, str), name
        assert kind in ("perf", "policy", "debug"), name


def test_scene_defaults_match_module_constants():
    assert float(knobs.KNOBS["PML_DISP_ALPHA"][0]) == boundaries.PML_DISP_ALPHA
    assert float(knobs.KNOBS["PML_POLE_DAMP"][0]) == boundaries.PML_POLE_DAMP
    assert float(knobs.KNOBS["LOSSLESS_DAMPING_REL"][0]) == poles.LOSSLESS_DAMPING_REL
    assert float(knobs.KNOBS["LOSSLESS_LOWFREQ_DAMPING_REL"][0]) == poles.LOSSLESS_LOWFREQ_DAMPING_REL


def test_scene_getters_read_env_each_call(monkeypatch):
    monkeypatch.delenv("OPENEM_PML_DISP_ALPHA", raising=False)
    assert boundaries._pml_disp_alpha() == boundaries.PML_DISP_ALPHA
    monkeypatch.setenv("OPENEM_PML_DISP_ALPHA", "0")
    assert boundaries._pml_disp_alpha() == 0.0
    monkeypatch.delenv("OPENEM_LOSSLESS_LOWFREQ_DAMPING_REL", raising=False)
    assert poles._lossless_lowfreq_damping_rel() == 0.0
    monkeypatch.setenv("OPENEM_LOSSLESS_LOWFREQ_DAMPING_REL", "0.3")
    assert poles._lossless_lowfreq_damping_rel() == 0.3
    monkeypatch.delenv("OPENEM_MODESRC_DISC", raising=False)
    assert modes._disc_on() is True
    monkeypatch.setenv("OPENEM_MODESRC_DISC", "0")
    assert modes._disc_on() is False


def test_half_cell_sign_branches_on_presence(monkeypatch):
    monkeypatch.delenv("OPENEM_HALF_CELL_SIGN", raising=False)
    assert knobs.env("HALF_CELL_SIGN") is None
    assert modes._half_cell_sign() == modes.HALF_CELL_SIGN
    monkeypatch.setenv("OPENEM_HALF_CELL_SIGN", "1")
    assert modes._half_cell_sign() == 1.0
