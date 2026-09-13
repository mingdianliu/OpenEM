# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""``serialize.save`` then ``serialize.load`` round-trip: every field of ``Scene`` compared array by
array with ``array_equal``.

What this guards against is a key name on the load side drifting out of step with the save side.
``pec_ey`` and ``pec_ez`` were once decided by testing ``"pec_ex" in z``, and only survived because
all three keys were always written together. The comparator walks dataclass fields recursively, so
a new field needs no change to the test.
"""

from __future__ import annotations

import dataclasses
import re

import numpy as np
import pytest

td = pytest.importorskip("tidy3d")

from dataclasses import replace

from openem import scene as scene_mod
from openem import serialize
from openem.model import BoxFluxMonitor, ModeMonitorSpec, ProjectionMonitor

F0 = 2.0e14
DL = 0.05


#: Fields deliberately kept out of the npz: the ``incident`` of a plane wave or TFSF source, which
#: exists only so model.lossless can recompute the incident table, and the complex tables of their
#: waveforms, since both source kinds fail closed on the complex-field path and the load side
#: rebuilds them as purely real.
_NOT_SERIALIZED = re.compile(
    r"^scene\.(sources|tfsf_sources)\[\d+\]\.(incident|waveform\.amp_(int|half)_complex)$")


def _assert_equal(a, b, where: str = "scene") -> None:
    """Compare recursively: a dataclass field by field, a sequence element by element, an array with
    ``array_equal`` and NaN counted as equal.
    """
    if _NOT_SERIALIZED.match(where):
        return
    if dataclasses.is_dataclass(a):
        assert type(a) is type(b), f"{where}: {type(a)} vs {type(b)}"
        for f in dataclasses.fields(a):
            _assert_equal(getattr(a, f.name), getattr(b, f.name), f"{where}.{f.name}")
    elif isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        assert a is not None and b is not None, f"{where}: {a!r} vs {b!r}"
        a, b = np.asarray(a), np.asarray(b)
        assert a.dtype == b.dtype, f"{where}: dtype {a.dtype} vs {b.dtype}"
        assert np.array_equal(a, b, equal_nan=True), f"{where}: arrays differ"
    elif isinstance(a, dict):
        assert set(a) == set(b), f"{where}: keys {set(a)} against {set(b)}"
        for k in a:
            _assert_equal(a[k], b[k], f"{where}[{k!r}]")
    elif isinstance(a, (list, tuple)):
        assert len(a) == len(b), f"{where}: lengths {len(a)} against {len(b)}"
        for i, (x, y) in enumerate(zip(a, b)):
            _assert_equal(x, y, f"{where}[{i}]")
    elif isinstance(a, float) and isinstance(b, float) and np.isnan(a) and np.isnan(b):
        pass
    else:
        assert a == b, f"{where}: {a!r} vs {b!r}"


def _plane_wave_sim() -> "td.Simulation":
    """Plane wave, a PEC slab, a lossy slab, and flux, frequency-domain field, permittivity,
    time-domain field and time-domain flux monitors, with a PML at one end of z and an Absorber at
    the other.
    """
    freqs = [0.8 * F0, F0, 1.2 * F0]
    per = td.Boundary.periodic()
    return td.Simulation(
        size=(8 * DL, 8 * DL, 4.0), run_time=1.0e-13,
        grid_spec=td.GridSpec.uniform(dl=DL),
        structures=[
            td.Structure(geometry=td.Box(center=(0, 0, 0.5), size=(td.inf, td.inf, 0.3)),
                         medium=td.PECMedium()),
            td.Structure(geometry=td.Box(center=(0, 0, -0.3), size=(td.inf, td.inf, 0.2)),
                         medium=td.Medium(permittivity=2.25, conductivity=5.0)),
        ],
        sources=[td.PlaneWave(
            center=(0, 0, -0.9), size=(td.inf, td.inf, 0), direction="+",
            pol_angle=0.0,
            source_time=td.GaussianPulse(freq0=F0, fwidth=0.25 * F0))],
        monitors=[
            td.FluxMonitor(center=(0, 0, -1.3), size=(td.inf, td.inf, 0),
                           freqs=freqs, name="back"),
            td.FluxMonitor(center=(0, 0, 1.2), size=(4 * DL, 4 * DL, 0),
                           freqs=freqs, name="trans",
                           apodization=td.ApodizationSpec(start=2e-14, width=5e-15)),
            td.FieldMonitor(center=(0, 0, 1.0), size=(td.inf, td.inf, 0),
                            freqs=freqs, name="fld"),
            td.PermittivityMonitor(center=(0, 0, 0.5), size=(td.inf, td.inf, 0.4),
                                   freqs=[F0], name="eps"),
            td.FieldTimeMonitor(center=(0, 0, 1.0), size=(2 * DL, 2 * DL, 0),
                                fields=["Ex", "Hy"], interval=7, name="tm"),
            td.FluxTimeMonitor(center=(0.01, -0.02, 1.1), size=(5 * DL, 6 * DL, 0),
                               interval=3, name="ftm"),
        ],
        boundary_spec=td.BoundarySpec(
            x=per, y=per,
            z=td.Boundary(minus=td.PML(), plus=td.Absorber(num_layers=20))),
        subpixel=False)


def _dipole_sim() -> "td.Simulation":
    """An electric and a magnetic point dipole, a dispersive medium box, PML on all sides, and a
    frequency-domain field monitor.
    """
    med = td.material_library["cSi"]["Green2008"]
    return td.Simulation(
        size=(1.0, 1.0, 1.0), run_time=5e-14,
        grid_spec=td.GridSpec.uniform(dl=DL),
        structures=[td.Structure(
            geometry=td.Box(center=(0.2, 0, 0), size=(0.3, 0.3, 0.3)), medium=med)],
        sources=[
            td.PointDipole(center=(-0.21, 0.02, 0.03), polarization="Ey",
                           source_time=td.GaussianPulse(freq0=F0, fwidth=0.3 * F0)),
            td.PointDipole(center=(-0.21, 0.02, -0.13), polarization="Hz",
                           source_time=td.GaussianPulse(freq0=F0, fwidth=0.3 * F0,
                                                        phase=0.7)),
        ],
        monitors=[td.FieldMonitor(center=(0, 0, 0), size=(0.4, 0.4, 0),
                                  freqs=[F0], name="xy")],
        boundary_spec=td.BoundarySpec.all_sides(boundary=td.PML()),
        subpixel=False)


@pytest.mark.parametrize("build", [_plane_wave_sim, _dipole_sim])
def test_roundtrip_bitwise(tmp_path, build):
    sc = scene_mod.from_simulation(build())
    # Also include the metadata the solver does not produce and that only passes through
    # serialization
    sc = replace(
        sc,
        projection_monitors=[ProjectionMonitor(
            name="ff", surface_names=("ff__n2f_x-", "ff__n2f_x+"), normal_dirs=("-", "+"))],
        box_flux_monitors=[BoxFluxMonitor(name="box", face_names=("box_x-", "box_x+"))],
        mode_monitors=[ModeMonitorSpec(name="mm", plane_name="mm__mode_plane")],
    )
    back = serialize.load(serialize.save(sc, tmp_path / "s.npz"))
    _assert_equal(sc, back)


def test_roundtrip_pec_masks_are_not_mixed_up(tmp_path):
    """The three PEC masks are all different, and reading them back must not cross the keys."""
    sc = scene_mod.from_simulation(_plane_wave_sim())
    assert sc.pec_ex is not None
    sc = replace(sc,
                 pec_ex=np.array([1, 2, 3], dtype=np.int64),
                 pec_ey=np.array([4, 5], dtype=np.int64),
                 pec_ez=np.array([6], dtype=np.int64))
    back = serialize.load(serialize.save(sc, tmp_path / "s.npz"))
    for f in ("pec_ex", "pec_ey", "pec_ez"):
        np.testing.assert_array_equal(getattr(back, f), getattr(sc, f))
