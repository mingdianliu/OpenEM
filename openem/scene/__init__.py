# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""``td.Simulation`` -> :class:`~openem.model.Scene`.

**This package and ``nb/`` (the notebook-side hooks) are the only two places that import tidy3d**;
``tests/test_no_tidy3d_leak.py`` locks that boundary. Below this layer there is only numpy, which
keeps the kernel-side interface narrow and means the tests do not need tidy3d.

Grid generation, geometry rasterization and material painting are not done here. The Tidy3D client
has already done them, consistently with its own solver, so the results are taken as they are.

Split by physics, with each piece responsible for failing closed on its own:

- :mod:`~openem.scene.media`        eps and conductivity, material types, PEC, tensors
- :mod:`~openem.scene.weights`      per-cell geometric weights, solved back out (arithmetic,
                                    harmonic and mixed models)
- :mod:`~openem.scene.poles`        extracting, splitting and artificially damping dispersion poles
- :mod:`~openem.scene.dispersion`   poles -> ADE entry tables, including the mix entries on a
                                    slanted interface
- :mod:`~openem.scene.modulation`   time-modulated media
- :mod:`~openem.scene.boundaries`   boundary types, CPML coefficients, topology, Absorber
- :mod:`~openem.scene.sources`      one-way plane wave injection, TFSF boxes, point dipoles,
                                    mirroring across a symmetry plane
- :mod:`~openem.scene.modes`        mode sources and mode monitors, the ModeSolver route
- :mod:`~openem.scene.monitors`     flux, frequency-domain field, time-domain field and
                                    projection monitors
- :mod:`~openem.scene.projection`   far-field projection: splitting the near-field surfaces and
                                    the colocated readout
- :mod:`~openem.scene.subpixel`     tidy3d's global local-subpixel switch and the notebook hooks

Plus three shared pieces that carry no physics: :mod:`~openem.scene._util` (Yee coordinates,
component keys, indices), :mod:`~openem.scene.cache` (hashing and atomic writes for the disk
cache), and :mod:`~openem.scene.normalize_td` (the normalization convention on the readout side).

:mod:`~openem.scene.build` does nothing but assemble them.

**Notebook-side interface**, i.e. what the nb.backend and startup hooks call directly, as opposed
to internal functions:

- ``subpixel.cloud_emulation`` / ``install_notebook_staircase`` (re-exported by media)
- ``modes.amplitudes`` / ``plane_monitor`` / ``mode_solver_data`` / ``memo_eps`` / ``PLANE_SUFFIX``
- ``projection.field_data`` / ``project`` / ``surface_monitors`` / ``normalization``
- ``sources.incident_tables`` (replaying the incident tables for ``model.lossless``)
"""

from openem.model import (
    FieldMonitor,
    FieldTimeMonitor,
    FluxMonitor,
    FluxTimeMonitor,
    PermittivityMonitor,
    PlaneWaveSource,
    PointDipole,
    Scene,
    homogeneous,
    lossless,
)
from openem.scene.build import from_file, from_simulation
from openem.scene.subpixel import require_local_subpixel, set_local_subpixel
from openem.units import UM

__all__ = [
    "Scene", "PlaneWaveSource", "PointDipole",
    "FluxMonitor", "FluxTimeMonitor", "FieldMonitor", "FieldTimeMonitor",
    "PermittivityMonitor",
    "homogeneous", "lossless", "from_simulation", "from_file",
    "set_local_subpixel", "require_local_subpixel", "UM",
]
