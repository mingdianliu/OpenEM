#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Golden regression set: proof that the solver's output is bitwise identical across a change.

Usage, from the repository root, on a machine with a GPU::

    pip install -e .            # or put the repo root on PYTHONPATH
    python3 tools/golden.py list                 # list the cases
    python3 tools/golden.py run base             # run them all -> golden/base/<case>.npz
    python3 tools/golden.py run base --only iso_pml,tfsf_box   # run only a few
    python3 tools/golden.py compare base base2   # compare case by case, array by array;
                                                 # the last line is GOLDEN_SAME or GOLDEN_DIFF n

How it runs: each case is a small ``td.Simulation`` taken all the way through
``openem.nb.backend.run``, i.e. build the scene, solve on the local GPU, assemble a
td.SimulationData. The script points ``openem.nb.backend.WORK_DIR`` at ``<output dir>/<tag>/work``
for the intermediate files (scene.npz, ours.npz). The output directory defaults to ``golden/``
under the current directory and can be moved with ``OPENEM_GOLDEN_DIR``. It also ``setdefault``s
``OPENEM_MODE_CACHE=0`` (recompute the mode basis every time, since the cache is itself something
to verify), ``OPENEM_NB_NOCACHE=1`` and ``OPENEM_DEBUG_LOG``.

What is stored: the ``.values`` of every DataArray of every monitor in the SimulationData (key
``<monitor>/<field>``) and its coordinates (key ``<monitor>/<field>@<dim>``), plus metadata under
``__meta`` (step count, early termination, cell count, nominal step count, wall time). ``wall`` is
excluded from the comparison; everything else is compared.

Comparison rule: bitwise identical by default (``np.array_equal``, with NaN counted as equal). A
case whose ``atol``/``rtol`` in CASES is not None is allowed a tolerance, for quantities that are
inherently non-deterministic; passing within tolerance counts as SAME but is marked ``~``.

Of the 22 cases, only grad_flux needs a tolerance; see CASES for why. Everything else is bitwise
identical between runs.

The cases cover: far-field projection, diagonal and full-tensor anisotropy, PEC walls with a PEC
structure, oblique TFSF, GaussianBeam, a two-mode source, apodization with a ModeSolverMonitor, and
three adjoint-gradient cases. A case constructor may return a ``td.Simulation``, or a
``{"arrays", "meta"}`` dict carrying its own way of running. The three adjoint cases take the
latter route through ``openem/nb/autograd_hook.py``: tidy3d's autograd runs both the forward and
the adjoint solve on the local GPU and the tidy3d client assembles the gradient. They store ``J``,
``dJ_dp`` and the forward monitor data, and each solve's ``(label, cells, steps)`` goes into
``__meta['calls']`` to be compared as well.

Note that the production solve path does **not** fold symmetric domains: the only remaining callers
of ``fold.py`` are the speed benchmark and ``tests/test_fold.py``. Symmetry instead goes through
the mirrored sources in scene/sources.py and the mirror-side zeroing in nb/shapegrad.py, and
grad_symm covers those two.

Two of the cases use ``subpixel=False`` (nosubpx_perm, nosubpx_mode). They lock the rule that the
``use_local_subpixel`` value set by from_simulation must be left in place for the readout code
inside the same cloud_emulation scope; see the comment above those two cases.

About early termination: ``sim.shutoff`` only controls the field and energy decay criteria. The
phasor convergence criterion (tol=1e-3) is always on and no environment variable turns it off. So
"run the full nominal step count" is arranged by keeping run_time short (disp_dipole,
absorber_mode, symm_mode, bloch_diff, tfsf_box, custom_field), while "early termination fires" is
covered by shutoff_fires (field and energy) and by iso_pml / timemod (phasor). Both the step count
and the stop reason go into ``__meta`` and are compared.
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import re
import sys
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]          # OpenEM_refactor/
GOLDEN = pathlib.Path(os.environ.get("OPENEM_GOLDEN_DIR", "golden"))

# ----------------------------------------------------------------------------- case definitions

C0 = 2.99792458e8


def _td():
    import tidy3d as td
    return td


def case_iso_pml():
    """Isotropic dielectric box, PML in z with periodic x and y, plane wave, and flux /
    frequency-domain field / time-domain field / time-domain flux / permittivity monitors. The
    phasor criterion terminates it early at about 3200 steps.
    """
    td = _td()
    f0 = C0 / 1.0e-6
    st = td.GaussianPulse(freq0=f0, fwidth=0.15 * f0)
    per = td.Boundary(plus=td.Periodic(), minus=td.Periodic())
    return td.Simulation(
        size=(0.6, 0.6, 2.0), grid_spec=td.GridSpec.uniform(dl=0.025),
        structures=[td.Structure(geometry=td.Box(center=(0, 0, 0), size=(0.35, 0.35, 0.3)),
                                 medium=td.Medium(permittivity=4.0))],
        sources=[td.PlaneWave(center=(0, 0, -0.7), size=(td.inf, td.inf, 0), direction="+", source_time=st)],
        monitors=[
            td.FluxMonitor(center=(0, 0, 0.7), size=(td.inf, td.inf, 0), freqs=[0.9 * f0, f0, 1.1 * f0], name="T"),
            td.FluxMonitor(center=(0, 0, -0.85), size=(td.inf, td.inf, 0), freqs=[f0], name="R"),
            td.FieldMonitor(center=(0, 0, 0), size=(td.inf, 0, td.inf), freqs=[f0], name="fxz"),
            td.FieldMonitor(center=(0, 0, 0.5), size=(td.inf, td.inf, 0), freqs=[f0], name="fxy", colocate=False),
            td.FieldTimeMonitor(center=(0.05, 0.05, 0.6), size=(0, 0, 0), name="pt", interval=5),
            td.FluxTimeMonitor(center=(0, 0, 0.7), size=(td.inf, td.inf, 0), name="Tt", interval=10),
            td.PermittivityMonitor(center=(0, 0, 0), size=(td.inf, 0, 0.6), freqs=[f0], name="eps"),
        ],
        run_time=2.0e-13, shutoff=0.0,
        boundary_spec=td.BoundarySpec(x=per, y=per, z=td.Boundary(plus=td.PML(), minus=td.PML())))


def case_shutoff_fires():
    """Lossy dielectric box, PML on all sides, plane wave with periodic transverse boundaries.
    shutoff=1e-5 with a long run_time, so early termination should fire.
    """
    td = _td()
    f0 = C0 / 1.0e-6
    st = td.GaussianPulse(freq0=f0, fwidth=0.2 * f0)
    per = td.Boundary(plus=td.Periodic(), minus=td.Periodic())
    return td.Simulation(
        size=(0.4, 0.4, 1.6), grid_spec=td.GridSpec.uniform(dl=0.025),
        structures=[td.Structure(geometry=td.Box(center=(0, 0, 0), size=(td.inf, td.inf, 0.3)),
                                 medium=td.Medium(permittivity=2.25, conductivity=0.02))],
        sources=[td.PlaneWave(center=(0, 0, -0.55), size=(td.inf, td.inf, 0), direction="+", source_time=st)],
        monitors=[
            td.FluxMonitor(center=(0, 0, 0.55), size=(td.inf, td.inf, 0), freqs=[f0], name="T"),
            td.FieldMonitor(center=(0, 0, 0), size=(td.inf, 0, td.inf), freqs=[f0], name="fxz"),
        ],
        run_time=3.0e-12, shutoff=1e-5,
        boundary_spec=td.BoundarySpec(x=per, y=per, z=td.Boundary(plus=td.PML(), minus=td.PML())))


def case_disp_dipole():
    """Dispersive medium (material_library cSi Green2008), point dipole, frequency-domain field
    and box flux monitors, PML on all sides. Runs the full nominal step count.
    """
    td = _td()
    f0 = C0 / 1.0e-6
    st = td.GaussianPulse(freq0=f0, fwidth=0.1 * f0)
    return td.Simulation(
        size=(0.6, 0.6, 0.6), grid_spec=td.GridSpec.uniform(dl=0.025),
        structures=[td.Structure(geometry=td.Box(center=(0, 0, 0.15), size=(0.3, 0.3, 0.15)),
                                 medium=td.material_library["cSi"]["Green2008"])],
        sources=[td.PointDipole(center=(0, 0, -0.1), polarization="Ex", source_time=st)],
        monitors=[
            td.FieldMonitor(center=(0, 0, 0), size=(td.inf, 0, td.inf), freqs=[0.95 * f0, f0], name="fxz"),
            td.FluxMonitor(center=(0, 0, 0), size=(0.4, 0.4, 0.4), freqs=[f0], name="box"),
        ],
        run_time=1.2e-13, shutoff=0.0,
        boundary_spec=td.BoundarySpec.all_sides(td.PML()))


def _wg_sim(boundary_spec, symmetry=(0, 0, 0), mode_index=0):
    td = _td()
    f0 = 1.934e14
    st = td.GaussianPulse(freq0=f0, fwidth=0.1 * f0)
    return td.Simulation(
        size=(2.0, 1.6, 1.4), grid_spec=td.GridSpec.uniform(dl=0.05), symmetry=symmetry,
        structures=[td.Structure(geometry=td.Box(center=(0, 0, 0), size=(td.inf, 0.5, 0.22)),
                                 medium=td.Medium(permittivity=12.1))],
        sources=[td.ModeSource(center=(-0.6, 0, 0), size=(0, 1.2, 1.0), source_time=st, direction="+",
                               mode_spec=td.ModeSpec(num_modes=1), mode_index=mode_index)],
        monitors=[
            td.ModeMonitor(center=(0.5, 0, 0), size=(0, 1.2, 1.0), freqs=[0.98 * f0, f0, 1.02 * f0],
                           mode_spec=td.ModeSpec(num_modes=2), name="m"),
            td.FluxMonitor(center=(0.5, 0, 0), size=(0, 1.2, 1.0), freqs=[f0], name="fl"),
            td.FieldMonitor(center=(0, 0, 0), size=(td.inf, td.inf, 0), freqs=[f0], name="fxy"),
        ],
        run_time=2.0e-13, shutoff=0.0, boundary_spec=boundary_spec)


def case_absorber_mode():
    """td.Absorber on all six faces, silicon strip waveguide, mode source, ModeMonitor plus flux
    and field monitors.
    """
    td = _td()
    return _wg_sim(td.BoundarySpec.all_sides(td.Absorber(num_layers=20)))


def case_symm_mode():
    """Symmetry (0,-1,1), waveguide, mode source, ModeMonitor, PML on all sides."""
    td = _td()
    return _wg_sim(td.BoundarySpec.all_sides(td.PML(num_layers=8)), symmetry=(0, -1, 1))


def case_bloch_diff():
    """Bloch periodic in x and periodic in y, obliquely incident plane wave (angle_theta>0),
    grating, DiffractionMonitor on both sides, plus flux.
    """
    td = _td()
    f0 = 4.0e14                       # lambda=0.75 um, period 1.0 um, so orders +-1 propagate
    st = td.GaussianPulse(freq0=f0, fwidth=0.1 * f0)
    size = (1.0, 0.1, 2.4)
    src = td.PlaneWave(center=(0, 0, -0.8), size=(td.inf, td.inf, 0), direction="+", source_time=st,
                       angle_theta=np.pi / 9, angle_phi=0.0, pol_angle=0.0)
    return td.Simulation(
        size=size, grid_spec=td.GridSpec.uniform(dl=0.025),
        structures=[
            td.Structure(geometry=td.Box(center=(0, 0, 0.6), size=(td.inf, td.inf, 1.2)),
                         medium=td.Medium(permittivity=2.25)),
            td.Structure(geometry=td.Box(center=(0, 0, 0.1), size=(0.5, td.inf, 0.2)),
                         medium=td.Medium(permittivity=2.25)),
        ],
        sources=[src],
        monitors=[
            td.DiffractionMonitor(center=(0, 0, -0.9), size=(td.inf, td.inf, 0), freqs=[f0], name="dR", normal_dir="-"),
            td.DiffractionMonitor(center=(0, 0, 0.9), size=(td.inf, td.inf, 0), freqs=[f0], name="dT", normal_dir="+"),
            td.FluxMonitor(center=(0, 0, 0.9), size=(td.inf, td.inf, 0), freqs=[f0], name="T"),
        ],
        run_time=3.0e-13, shutoff=0.0,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.bloch_from_source(source=src, domain_size=size[0], axis=0),
            y=td.Boundary.periodic(), z=td.Boundary.pml()))


def case_tfsf_box():
    """TFSF box at normal incidence, infinite substrate, PML on all sides, field and flux monitors."""
    td = _td()
    f0 = 6.0e14
    src = td.TFSF(center=(0, 0, 0), size=(0.4, 0.4, 0.4), source_time=td.GaussianPulse(freq0=f0, fwidth=f0 / 6),
                  injection_axis=2, direction="+", angle_theta=0.0, pol_angle=0.0)
    return td.Simulation(
        size=(0.8, 0.8, 0.8), grid_spec=td.GridSpec.uniform(dl=0.02),
        structures=[td.Structure(geometry=td.Box.from_bounds((-100, -100, -100), (100, 100, -0.1)),
                                 medium=td.Medium(permittivity=2.25))],
        sources=[src],
        monitors=[
            td.FieldMonitor(center=(0, 0, 0), size=(td.inf, 0, td.inf), freqs=[f0], name="fxz"),
            td.FluxMonitor(center=(0, 0, 0.3), size=(0.3, 0.3, 0), freqs=[f0], name="fl"),
        ],
        run_time=6.0e-14, shutoff=0.0,
        boundary_spec=td.BoundarySpec.all_sides(td.PML()))


def case_timemod():
    """Time-modulated medium (ContinuousWaveTimeModulation) slab, plane wave, multi-frequency flux,
    normalize_index=None.
    """
    td = _td()
    f0 = C0 / 1.0e-6
    fm = 0.1 * f0
    spec = td.ModulationSpec(permittivity=td.SpaceTimeModulation(
        time_modulation=td.ContinuousWaveTimeModulation(freq0=fm, amplitude=1.0, phase=0.0),
        space_modulation=td.SpaceModulation(amplitude=0.4, phase=0.0)))
    dl, lz = 0.05, 3.0
    return td.Simulation(
        size=(4 * dl, 4 * dl, lz), grid_spec=td.GridSpec.uniform(dl=dl),
        structures=[td.Structure(geometry=td.Box(center=(0, 0, 0), size=(td.inf, td.inf, 1.0)),
                                 medium=td.Medium(permittivity=2.0, modulation_spec=spec))],
        sources=[td.PlaneWave(source_time=td.GaussianPulse(freq0=f0, fwidth=f0 / 10),
                              size=(td.inf, td.inf, 0), center=(0, 0, -lz / 2 + 0.4), direction="+")],
        monitors=[td.FluxMonitor(center=(0, 0, lz / 2 - 0.4), size=(td.inf, td.inf, 0),
                                 freqs=[f0 - fm, f0, f0 + fm], name="flux")],
        run_time=8.0e-13, shutoff=0.0, normalize_index=None,
        boundary_spec=td.BoundarySpec.pml(x=False, y=False, z=True))


def case_custom_field():
    """CustomFieldSource built from an Ex and Hy plane dataset, PML on all sides, frequency-domain
    field monitor.
    """
    td = _td()
    f0 = C0 / 0.94e-6
    n = 6
    ax = np.linspace(-0.3, 0.3, n)
    coords = dict(x=ax, y=ax, z=np.array([0.0]), f=np.array([f0]))

    def sfa(v):
        return td.ScalarFieldDataArray(np.full((n, n, 1, 1), v, dtype=complex), coords=coords)

    src = td.CustomFieldSource(center=(0, 0, 0), size=(0.4, 0.4, 0),
                               source_time=td.GaussianPulse(freq0=f0, fwidth=f0 / 10),
                               field_dataset=td.FieldDataset(Ex=sfa(1.0), Hy=sfa(1.0 / 376.73)))
    return td.Simulation(
        size=(0.8, 0.8, 1.2), grid_spec=td.GridSpec.uniform(dl=0.05), sources=[src],
        monitors=[td.FieldMonitor(center=(0, 0, 0.4), size=(td.inf, td.inf, 0), freqs=[f0], name="f"),
                  td.FluxMonitor(center=(0, 0, 0.4), size=(td.inf, td.inf, 0), freqs=[f0], name="fl")],
        run_time=1.0e-13, shutoff=0.0,
        boundary_spec=td.BoundarySpec.all_sides(td.PML()))


def case_proj_far():
    """Point dipole in vacuum, FieldProjectionAngleMonitor over a closed box and
    FieldProjectionCartesianMonitor over a single face, PML on all sides. Exercises the six-face
    split in scene/projection.py and the client-side projection.
    """
    td = _td()
    f0 = 2.0e14
    return td.Simulation(
        size=(1.2, 1.2, 1.2), grid_spec=td.GridSpec.uniform(dl=0.04),
        sources=[td.PointDipole(center=(0, 0, 0), polarization="Ez",
                                source_time=td.GaussianPulse(freq0=f0, fwidth=0.2 * f0))],
        monitors=[
            td.FieldProjectionAngleMonitor(center=(0, 0, 0), size=(0.8, 0.8, 0.8), freqs=[f0], name="far",
                                           theta=np.linspace(0, np.pi, 7), phi=[0.0, np.pi / 2], proj_distance=1e6),
            td.FieldProjectionCartesianMonitor(center=(0, 0, 0.4), size=(0.8, 0.8, 0), freqs=[f0], name="cart",
                                               x=[-1.0, 0.0, 1.0], y=[-1.0, 0.0, 1.0], proj_axis=2,
                                               proj_distance=1e6, normal_dir="+"),
            td.FluxMonitor(center=(0, 0, 0), size=(0.8, 0.8, 0.8), freqs=[f0], name="box"),
        ],
        run_time=8.0e-14, shutoff=0.0,
        boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=10)))


def case_aniso_diag():
    """Diagonally anisotropic slab, td.AnisotropicMedium with distinct xx/yy/zz, plus an off-center
    Ey point dipole so all three components are exercised, PML on all sides.
    """
    td = _td()
    f0 = C0 / 1.0e-6
    med = td.AnisotropicMedium(xx=td.Medium(permittivity=2.0), yy=td.Medium(permittivity=3.0),
                               zz=td.Medium(permittivity=4.0))
    return td.Simulation(
        size=(0.8, 0.8, 1.6), grid_spec=td.GridSpec.uniform(dl=0.04),
        structures=[td.Structure(geometry=td.Box(center=(0, 0, 0.1), size=(td.inf, td.inf, 0.4)), medium=med)],
        sources=[td.PointDipole(center=(0.05, 0.03, -0.5), polarization="Ey",
                                source_time=td.GaussianPulse(freq0=f0, fwidth=0.15 * f0))],
        monitors=[
            td.FluxMonitor(center=(0, 0, 0.6), size=(td.inf, td.inf, 0), freqs=[0.9 * f0, f0], name="T"),
            td.FluxMonitor(center=(0, 0, -0.7), size=(td.inf, td.inf, 0), freqs=[f0], name="R"),
            td.FieldMonitor(center=(0, 0, 0), size=(td.inf, 0, td.inf), freqs=[f0], name="fxz"),
        ],
        run_time=1.5e-13, shutoff=0.0,
        boundary_spec=td.BoundarySpec.all_sides(td.PML()))


def case_tensor_full():
    """Full-tensor td.FullyAnisotropicMedium slab: a diagonal tensor rotated 30 degrees about z, so
    off-diagonal coupling and loss, with a plane wave. Goes through kernels/tensor.cu.
    subpixel=False, since tensors are only ever staircased.
    """
    td = _td()
    f0 = 2.0e14
    th = np.pi / 6
    R = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1.0]])
    eps = R @ np.diag([2.0, 3.0, 4.0]) @ R.T
    sig = R @ np.diag([1e-4, 2e-4, 3e-4]) @ R.T
    med = td.FullyAnisotropicMedium(permittivity=eps.tolist(), conductivity=sig.tolist())
    per = td.Boundary(plus=td.Periodic(), minus=td.Periodic())
    return td.Simulation(
        size=(0.4, 0.4, 2.4), grid_spec=td.GridSpec.uniform(dl=0.05), subpixel=False,
        structures=[td.Structure(geometry=td.Box(center=(0, 0, 0.3), size=(td.inf, td.inf, 0.6)), medium=med)],
        sources=[td.PlaneWave(center=(0, 0, -0.8), size=(td.inf, td.inf, 0), direction="+",
                              source_time=td.GaussianPulse(freq0=f0, fwidth=0.3 * f0))],
        monitors=[
            td.FieldMonitor(center=(0, 0, 0.3), size=(0.2, 0.2, 0.4), freqs=[f0], name="probe"),
            td.FluxMonitor(center=(0, 0, 0.9), size=(td.inf, td.inf, 0), freqs=[f0], name="T"),
        ],
        run_time=2.0e-13, shutoff=0.0,
        boundary_spec=td.BoundarySpec(x=per, y=per, z=td.Boundary(plus=td.PML(), minus=td.PML())))


def case_pec_walls():
    """PECBoundary on the four x and y walls, PML in z, a PEC structure box (td.PECMedium as a
    staircased mask), an off-center point dipole, and field / flux / time-domain point monitors.
    """
    td = _td()
    f0 = C0 / 1.0e-6
    pec = td.Boundary(plus=td.PECBoundary(), minus=td.PECBoundary())
    return td.Simulation(
        size=(0.8, 0.8, 1.2), grid_spec=td.GridSpec.uniform(dl=0.04),
        structures=[td.Structure(geometry=td.Box(center=(0.05, 0, 0.1), size=(0.3, 0.3, 0.2)), medium=td.PECMedium())],
        sources=[td.PointDipole(center=(0.1, 0.05, -0.3), polarization="Ez",
                                source_time=td.GaussianPulse(freq0=f0, fwidth=0.2 * f0))],
        monitors=[
            td.FieldMonitor(center=(0, 0, 0), size=(td.inf, 0, td.inf), freqs=[f0], name="fxz"),
            td.FluxMonitor(center=(0, 0, 0.4), size=(td.inf, td.inf, 0), freqs=[f0], name="T"),
            td.FieldTimeMonitor(center=(-0.2, 0.1, -0.1), size=(0, 0, 0), name="pt", interval=4),
        ],
        run_time=1.2e-13, shutoff=0.0,
        boundary_spec=td.BoundarySpec(x=pec, y=pec, z=td.Boundary(plus=td.PML(), minus=td.PML())))


def case_tfsf_oblique():
    """Oblique TFSF (angle_theta=pi/7, angle_phi=pi/6, pol_angle=0.3, through tfsf_oblique.py),
    with a small dielectric cube inside the box and PML on all sides.
    """
    td = _td()
    f0 = 6.0e14
    src = td.TFSF(center=(0, 0, 0), size=(0.4, 0.4, 0.4), source_time=td.GaussianPulse(freq0=f0, fwidth=f0 / 6),
                  injection_axis=2, direction="+", angle_theta=np.pi / 7, angle_phi=np.pi / 6, pol_angle=0.3)
    return td.Simulation(
        size=(0.8, 0.8, 0.8), grid_spec=td.GridSpec.uniform(dl=0.02),
        structures=[td.Structure(geometry=td.Box(center=(0.02, -0.02, 0), size=(0.16, 0.16, 0.16)),
                                 medium=td.Medium(permittivity=2.25))],
        sources=[src],
        monitors=[
            td.FieldMonitor(center=(0, 0, 0), size=(td.inf, 0, td.inf), freqs=[f0], name="fxz"),
            td.FluxMonitor(center=(0, 0, 0.3), size=(0.3, 0.3, 0), freqs=[f0], name="fl"),
        ],
        run_time=6.0e-14, shutoff=0.0,
        boundary_spec=td.BoundarySpec.all_sides(td.PML()))


def case_gauss_beam():
    """GaussianBeam at a small oblique angle with a 0.4 um waist, onto a dielectric slab, PML on
    all sides, flux and field monitors. Exercises the beam-to-ModeSource carrier in scene/modes.py.
    """
    td = _td()
    f0 = C0 / 1.0e-6
    src = td.GaussianBeam(center=(0, 0, -0.4), size=(td.inf, td.inf, 0), direction="+", waist_radius=0.4,
                          angle_theta=0.15, pol_angle=0.0,
                          source_time=td.GaussianPulse(freq0=f0, fwidth=0.1 * f0))
    return td.Simulation(
        size=(1.6, 1.6, 1.2), grid_spec=td.GridSpec.uniform(dl=0.04),
        structures=[td.Structure(geometry=td.Box(center=(0, 0, 0.1), size=(td.inf, td.inf, 0.3)),
                                 medium=td.Medium(permittivity=2.25))],
        sources=[src],
        monitors=[
            td.FluxMonitor(center=(0, 0, 0.45), size=(td.inf, td.inf, 0), freqs=[f0], name="T"),
            td.FieldMonitor(center=(0, 0, 0), size=(td.inf, 0, td.inf), freqs=[f0], name="fxz"),
        ],
        run_time=1.0e-13, shutoff=0.0,
        boundary_spec=td.BoundarySpec.all_sides(td.PML()))


def _wg_structure():
    td = _td()
    return td.Structure(geometry=td.Box(center=(0, 0, 0), size=(td.inf, 0.5, 0.22)),
                        medium=td.Medium(permittivity=12.1))


def case_two_mode_src():
    """Two ModeSources injecting at once, on opposite planes x=-0.6 and x=+0.6, facing each other,
    with different center frequencies and phases. ModeMonitor and field monitor in between, PML on
    all sides.
    """
    td = _td()
    f0 = 1.934e14

    def src(x, direction, freq0, phase):
        return td.ModeSource(center=(x, 0, 0), size=(0, 1.2, 1.0), direction=direction,
                             source_time=td.GaussianPulse(freq0=freq0, fwidth=0.1 * freq0, phase=phase),
                             mode_spec=td.ModeSpec(num_modes=1), mode_index=0)
    return td.Simulation(
        size=(2.0, 1.6, 1.4), grid_spec=td.GridSpec.uniform(dl=0.05),
        structures=[_wg_structure()],
        sources=[src(-0.6, "+", f0, 0.0), src(0.6, "-", 1.05 * f0, 0.7)],
        monitors=[
            td.ModeMonitor(center=(0, 0, 0), size=(0, 1.2, 1.0), freqs=[f0, 1.05 * f0],
                           mode_spec=td.ModeSpec(num_modes=1), name="m"),
            td.FieldMonitor(center=(0, 0, 0), size=(td.inf, td.inf, 0), freqs=[f0], name="fxy"),
        ],
        run_time=2.0e-13, shutoff=0.0, boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=8)))


def case_apod_msolver():
    """Monitor apodization: an ApodizationSpec with both start and end on the FluxMonitor and with
    start only on the FieldMonitor, plus a ModeSolverMonitor (local mode solve written back) and a
    waveguide mode source.
    """
    td = _td()
    f0 = 1.934e14
    return td.Simulation(
        size=(2.0, 1.6, 1.4), grid_spec=td.GridSpec.uniform(dl=0.05),
        structures=[_wg_structure()],
        sources=[td.ModeSource(center=(-0.6, 0, 0), size=(0, 1.2, 1.0), direction="+",
                               source_time=td.GaussianPulse(freq0=f0, fwidth=0.1 * f0),
                               mode_spec=td.ModeSpec(num_modes=1), mode_index=0)],
        monitors=[
            td.FluxMonitor(center=(0.5, 0, 0), size=(0, 1.2, 1.0), freqs=[f0], name="fl",
                           apodization=td.ApodizationSpec(start=4e-14, end=1.5e-13, width=2e-14)),
            td.FieldMonitor(center=(0, 0, 0), size=(td.inf, td.inf, 0), freqs=[f0], name="fxy",
                            apodization=td.ApodizationSpec(start=6e-14, width=3e-14)),
            td.ModeSolverMonitor(center=(0.3, 0, 0), size=(0, 1.2, 1.0), freqs=[f0, 1.02 * f0],
                                 mode_spec=td.ModeSpec(num_modes=2), name="ms"),
        ],
        run_time=2.0e-13, shutoff=0.0, boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=8)))


# --------------------------------------------------------------- readout convention with subpixel=False
#
# Once from_simulation sets the global use_local_subpixel switch to bool(sim.subpixel), it must
# **not** be restored: the readout code that runs after the scene is built, inside the same
# cloud_emulation scope, depends on it. That readout is _perm_data in openem.nb.backend calling
# sim.epsilon directly, and modes.amplitudes / mode_solver_data / _mode_solver_monitor_data going
# through the local ModeSolver. Restoring it to None makes tidy3d fall back to "use extras if
# present" and take the C++ engine, which changes the convention for PermittivityMonitor's
# epsilon(omega) and for the ModeMonitor / ModeSolverMonitor mode basis in a subpixel=False
# simulation. The first 20 cases that carry those monitors all use subpixel=True, where None and
# True produce the same numbers in the presence of extras, so they cannot catch it. The two cases
# below use subpixel=False precisely for that. Verified in reverse: restoring build.py to the
# ``with media.local_subpixel(...)`` spelling makes both of these cases DIFF.

def case_nosubpx_perm():
    """subpixel=False with a lossy dielectric cylinder (many mixed cells, epsilon(omega) varying
    with frequency), a small PEC box offset by half a cell, a three-frequency PermittivityMonitor
    on the xz plane, field and flux monitors, and a plane wave. Locks the epsilon convention in
    _perm_data.
    """
    td = _td()
    f0 = C0 / 1.0e-6
    per = td.Boundary(plus=td.Periodic(), minus=td.Periodic())
    return td.Simulation(
        size=(0.8, 0.8, 1.4), grid_spec=td.GridSpec.uniform(dl=0.025), subpixel=False,
        structures=[
            # Cylinder axis along y: a circle on the xz plane, mixed cells all around it. The
            # length of 0.33 puts the end faces inside a cell along y as well.
            td.Structure(geometry=td.Cylinder(center=(0, 0, 0), radius=0.21, length=0.33, axis=1),
                         medium=td.Medium(permittivity=4.0, conductivity=0.05)),
            # Small PEC box, with center and edge lengths deliberately off the grid
            # (0.21 / 0.01 / -0.31, edges 0.11 / 0.06). The monitor plane y=0 cuts through it.
            td.Structure(geometry=td.Box(center=(0.21, 0.01, -0.31), size=(0.11, 0.11, 0.06)), medium=td.PECMedium()),
        ],
        sources=[td.PlaneWave(center=(0, 0, -0.5), size=(td.inf, td.inf, 0), direction="+",
                              source_time=td.GaussianPulse(freq0=f0, fwidth=0.15 * f0))],
        monitors=[
            td.PermittivityMonitor(center=(0, 0, 0), size=(td.inf, 0, td.inf), freqs=[0.9 * f0, f0, 1.1 * f0], name="eps"),
            td.FieldMonitor(center=(0, 0, 0), size=(td.inf, 0, td.inf), freqs=[f0], name="fxz"),
            td.FluxMonitor(center=(0, 0, 0.55), size=(td.inf, td.inf, 0), freqs=[f0], name="T"),
        ],
        run_time=1.5e-13, shutoff=0.0,
        boundary_spec=td.BoundarySpec(x=per, y=per, z=td.Boundary(plus=td.PML(), minus=td.PML())))


def case_nosubpx_mode():
    """subpixel=False with a silicon strip waveguide (0.22 um thick at dl=0.05, so the interfaces
    fall inside cells), a mode source, a ModeMonitor (2 modes, 3 frequencies), a ModeSolverMonitor
    (2 modes, 2 frequencies) and flux, PML on all sides. Locks the local ModeSolver convention.
    """
    td = _td()
    f0 = 1.934e14
    return td.Simulation(
        size=(2.0, 1.6, 1.4), grid_spec=td.GridSpec.uniform(dl=0.05), subpixel=False,
        structures=[_wg_structure()],
        sources=[td.ModeSource(center=(-0.6, 0, 0), size=(0, 1.2, 1.0), direction="+",
                               source_time=td.GaussianPulse(freq0=f0, fwidth=0.1 * f0),
                               mode_spec=td.ModeSpec(num_modes=1), mode_index=0)],
        monitors=[
            td.ModeMonitor(center=(0.5, 0, 0), size=(0, 1.2, 1.0), freqs=[0.98 * f0, f0, 1.02 * f0],
                           mode_spec=td.ModeSpec(num_modes=2), name="m"),
            td.ModeSolverMonitor(center=(0.3, 0, 0), size=(0, 1.2, 1.0), freqs=[f0, 1.02 * f0],
                                 mode_spec=td.ModeSpec(num_modes=2), name="ms"),
            td.FluxMonitor(center=(0.5, 0, 0), size=(0, 1.2, 1.0), freqs=[f0], name="fl"),
        ],
        run_time=2.0e-13, shutoff=0.0, boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=8)))


# ----------------------------------------------------------------------------- adjoint gradient cases

def _unbox(x):
    """Strip autograd's ArrayBox. Under tracing, ``.values`` of a DataArray in the forward
    SimulationData is an ArrayBox.
    """
    while hasattr(x, "_value"):
        x = x._value
    return x


def _run_grad(name: str, make_sim, objective, p0) -> dict:
    """Run through ``autograd_hook``: tidy3d's autograd puts both the forward and the adjoint solve
    on the local GPU, and the tidy3d client assembles the gradient.

    ``make_sim(p)`` builds a ``td.Simulation`` from the parameter vector; when p carries a tracer,
    so does the structure. ``objective(sim_data)`` returns the scalar objective J. Stored: ``J``,
    ``dJ_dp``, and every monitor's data from the forward SimulationData.

    With ``install(steps=None)`` the forward solve stops on its own criterion and the adjoint
    copies its step count (``_FWD_STEPS``). Each solve's ``(label, cells, steps)`` goes into
    ``__meta['calls']`` and is compared.
    """
    import autograd as ag
    from openem.nb import backend as nb
    from openem.nb import autograd_hook as H
    from openem import scene as scene_mod

    if H._STATE.get("autograd_run") is None:
        H.install(steps=None)
    H.CALLS.clear()
    H._FWD_STEPS.clear()
    stash = {}

    def f(p):
        sd = nb.openem_run(make_sim(p), task_name=name, verbose=False)
        stash["sd"] = sd
        return objective(sd)

    p0 = np.asarray(p0, dtype=float)
    J, g = ag.value_and_grad(f)(p0)
    sd = stash["sd"]
    arrays = _collect_arrays(sd, sd.simulation)
    arrays["J"] = np.asarray(float(_unbox(J)))
    arrays["dJ_dp"] = np.asarray(_unbox(g), dtype=float)
    sc = scene_mod.from_simulation(make_sim(p0))
    calls = [[str(t), int(c), int(s)] for t, c, s, _ in H.CALLS]
    fwd = [c for c in H.CALLS if not str(c[0]).startswith("batch:")]
    meta = dict(steps_run=int(fwd[0][2]) if fwd else -1, triggered=None, stop="", source_end=-1,
                shape=[int(x) for x in sc.shape], cells=int(np.prod(sc.shape)),
                num_time_steps=int(sc.num_time_steps), calls=calls,
                wall=float(sum(c[3] for c in H.CALLS)))
    return dict(arrays=arrays, meta=meta)


def _grad_sim(p, symmetry, structures):
    """The waveguide simulation shared by the adjoint cases: silicon strip waveguide, mode source,
    ModeMonitor and FluxMonitor, PML on all sides.
    """
    td = _td()
    f0 = 1.934e14
    return td.Simulation(
        size=(2.0, 1.6, 1.4), grid_spec=td.GridSpec.uniform(dl=0.05), symmetry=symmetry,
        structures=[_wg_structure()] + structures,
        sources=[td.ModeSource(center=(-0.6, 0, 0), size=(0, 1.2, 1.0), direction="+",
                               source_time=td.GaussianPulse(freq0=f0, fwidth=0.1 * f0),
                               mode_spec=td.ModeSpec(num_modes=1), mode_index=0)],
        monitors=[
            td.ModeMonitor(center=(0.6, 0, 0), size=(0, 1.2, 1.0), freqs=[f0],
                           mode_spec=td.ModeSpec(num_modes=1), name="m"),
            # enable_adjoint: by default tidy3d 2.12 does not store the fields a FluxMonitor needs
            # for the adjoint, so it has to be turned on when the flux is the objective
            td.FluxMonitor(center=(0.6, 0, 0), size=(0, 1.2, 1.0), freqs=[f0], name="fl", enable_adjoint=True),
        ],
        run_time=2.0e-13, shutoff=0.0, boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=8)))


def case_grad_box():
    """Adjoint gradient with no symmetry: a td.Box above the waveguide, with center_x, size_x and
    Medium.permittivity all differentiable, against the ModeMonitor |amp|^2.
    """
    import autograd.numpy as anp
    td = _td()
    f0 = 1.934e14

    def make_sim(p):
        st = td.Structure(geometry=td.Box(center=(p[1], 0, 0.21), size=(p[2], 0.6, 0.2)),
                          medium=td.Medium(permittivity=p[0]))
        return _grad_sim(p, (0, 0, 0), [st])

    def objective(sd):
        amp = sd["m"].amps.sel(direction="+", f=f0, mode_index=0)
        return anp.sum(anp.abs(amp.values) ** 2)

    return _run_grad("grad_box", make_sim, objective, [4.0, 0.2, 0.4])


def case_grad_symm():
    """Adjoint gradient with symmetry=(0,-1,1): a mirrored pair of td.Box in a GeometryGroup on
    either side of the waveguide, with size and permittivity differentiable, against the
    ModeMonitor |amp|^2. The geometry gradient of the mirror-side member is zeroed by
    nb/shapegrad.py.
    """
    import autograd.numpy as anp
    td = _td()
    f0 = 1.934e14

    def make_sim(p):
        geo = td.GeometryGroup(geometries=[
            td.Box(center=(0.2, 0.45, 0), size=(p[1], p[2], 0.22)),
            td.Box(center=(0.2, -0.45, 0), size=(p[1], p[2], 0.22)),
        ])
        st = td.Structure(geometry=geo, medium=td.Medium(permittivity=p[0]))
        return _grad_sim(p, (0, -1, 1), [st])

    def objective(sd):
        amp = sd["m"].amps.sel(direction="+", f=f0, mode_index=0)
        return anp.sum(anp.abs(amp.values) ** 2)

    return _run_grad("grad_symm", make_sim, objective, [4.0, 0.4, 0.3])


def case_grad_flux():
    """Adjoint gradient with a FluxMonitor flux objective and no symmetry. The adjoint source is a
    td.CustomCurrentSource: point-by-point current expanded into a batch of PointDipole, through
    scene/sources.py::custom_current. See CASES for the rtol.
    """
    import autograd.numpy as anp
    td = _td()
    f0 = 1.934e14

    def make_sim(p):
        st = td.Structure(geometry=td.Box(center=(p[1], 0, 0.21), size=(p[2], 0.6, 0.2)),
                          medium=td.Medium(permittivity=p[0]))
        return _grad_sim(p, (0, 0, 0), [st])

    def objective(sd):
        return anp.sum(sd["fl"].flux.sel(f=f0).values)

    return _run_grad("grad_flux", make_sim, objective, [4.0, 0.2, 0.4])


#: (name, constructor, atol, rtol). Both atol and rtol None means bitwise identity is required.
#: The constructor returns either a ``td.Simulation``, which goes through openem_run, or an
#: ``{"arrays", "meta"}`` dict carrying its own way of running, as the adjoint cases do.
CASES = [
    ("iso_pml", case_iso_pml, None, None),
    ("shutoff_fires", case_shutoff_fires, None, None),
    ("disp_dipole", case_disp_dipole, None, None),
    ("absorber_mode", case_absorber_mode, None, None),
    ("symm_mode", case_symm_mode, None, None),
    ("bloch_diff", case_bloch_diff, None, None),
    ("tfsf_box", case_tfsf_box, None, None),
    ("timemod", case_timemod, None, None),
    ("custom_field", case_custom_field, None, None),
    ("proj_far", case_proj_far, None, None),
    ("aniso_diag", case_aniso_diag, None, None),
    ("tensor_full", case_tensor_full, None, None),
    ("pec_walls", case_pec_walls, None, None),
    ("tfsf_oblique", case_tfsf_oblique, None, None),
    ("gauss_beam", case_gauss_beam, None, None),
    ("two_mode_src", case_two_mode_src, None, None),
    ("apod_msolver", case_apod_msolver, None, None),
    ("nosubpx_perm", case_nosubpx_perm, None, None),
    ("nosubpx_mode", case_nosubpx_mode, None, None),
    ("grad_box", case_grad_box, None, None),
    ("grad_symm", case_grad_symm, None, None),
    # grad_flux's dJ_dp is not bitwise identical between runs; measured relative differences ran
    # from 4e-7 to 3.7e-6 over five runs. The adjoint source of a flux objective is a
    # CustomCurrentSource expanding into several thousand PointDipole, and the batched injection
    # kernel in kernels/source_dft.cu accumulates with atomicAdd because the stencils overlap, so
    # the float summation order varies. An in-process probe confirmed it: a mode objective
    # (ModeSource adjoint) is bitwise identical, while a flux objective is not, with or without
    # symmetry. The forward data and J stay bitwise identical; only dJ_dp takes this tolerance.
    # rtol is about 25 times the largest measured difference, so a real change in the numbers
    # (of order 1e-3 or more) is still caught.
    ("grad_flux", case_grad_flux, 1e-12, 1e-4),
]


# ----------------------------------------------------------------------------- run

def _label(task_name: str) -> str:
    """Matches the label rule in openem.nb.backend._openem_run_inner, with OPENEM_NB_TAG empty."""
    return "openem-nb-" + re.sub(r"[^a-z0-9-]", "-", task_name.lower())[:24]


def _collect_arrays(sim_data, sim) -> dict:
    """Flatten every monitor's DataArray in the SimulationData into {key: ndarray}."""
    out = {}
    for mon in sim.monitors:
        md = sim_data[mon.name]
        fields = getattr(md, "model_fields", None) or getattr(type(md), "__fields__", {})
        for fname in fields:
            if fname == "monitor":
                continue
            val = getattr(md, fname, None)
            if val is None:
                continue
            if hasattr(val, "values") and hasattr(val, "dims"):
                arr = np.asarray(_unbox(val.values))
                out[f"{mon.name}/{fname}"] = arr
                for d in val.dims:
                    out[f"{mon.name}/{fname}@{d}"] = np.asarray(val.coords[d].values)
            elif isinstance(val, (tuple, list)) and val and all(isinstance(x, (int, float, complex)) for x in val):
                out[f"{mon.name}/{fname}"] = np.asarray(val)
    return out


def run(tag: str, only: list[str] | None) -> int:
    outdir = GOLDEN / tag
    workdir = outdir / "work"
    outdir.mkdir(parents=True, exist_ok=True)
    workdir.mkdir(parents=True, exist_ok=True)

    # Solve locally, submit nothing. The mode solve skips the on-disk cache and is recomputed every
    # time, since the cache is itself something to verify.
    os.environ["OPENEM_INJOB"] = "1"
    os.environ.setdefault("OPENEM_MODE_CACHE", "0")
    os.environ.setdefault("OPENEM_NB_NOCACHE", "1")
    os.environ.setdefault("OPENEM_DEBUG_LOG", str(outdir / "openem_debug.log"))
    sys.path.insert(0, str(ROOT))
    from openem.nb import backend as nb
    from openem import serialize
    nb.WORK_DIR = workdir                       # put the work directory on local disk (see module docstring)

    sel = [c for c in CASES if not only or c[0] in only]
    if only:
        miss = set(only) - {c[0] for c in sel}
        if miss:
            print(f"unknown case: {sorted(miss)}; available: {[c[0] for c in CASES]}")
            return 2
    fails = 0
    summary = []
    for name, ctor, _a, _r in sel:
        t0 = time.time()
        try:
            obj = ctor()
            if isinstance(obj, dict):          # a case carrying its own run: {"arrays", "meta"}
                arrays = obj["arrays"]
                m = dict(case=name, **obj["meta"])
            else:
                sim = obj
                sd = nb.openem_run(sim, task_name=name, verbose=False)
                arrays = _collect_arrays(sd, sim)
                work = workdir / "" / _label(name)
                meta = json.loads(np.load(work / "ours.npz")["__meta"].item())
                sc = serialize.load(work / "scene.npz")
                m = dict(case=name, steps_run=int(meta["steps_run"]), triggered=bool(meta["triggered"]),
                         stop=str(meta.get("stop") or ""), source_end=int(meta.get("source_end", -1)),
                         shape=[int(x) for x in sc.shape], cells=int(np.prod(sc.shape)),
                         num_time_steps=int(sc.num_time_steps), wall=float(meta.get("wall", 0.0)))
            m.update(n_arrays=len(arrays), total_s=time.time() - t0)
            np.savez(outdir / f"{name}.npz", __meta=json.dumps(m), **arrays)
            summary.append(m)
            print(f"[golden] {name}: grid {m['shape']} = {m['cells']:,}, steps {m['steps_run']:,}/{m['num_time_steps']:,}"
                  f" (early stop: {m['stop'] or 'none'}), arrays {m['n_arrays']}, GPU {m['wall']:.1f}s, total {m['total_s']:.1f}s", flush=True)
        except Exception as e:
            fails += 1
            import traceback
            traceback.print_exc()
            print(f"[golden] {name}: failed {type(e).__name__}: {str(e)[:300]}", flush=True)
    (outdir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=1))
    print(f"[golden] finished {len(summary)}/{len(sel)} cases -> {outdir}; {fails} failed")
    return 1 if fails else 0


# ----------------------------------------------------------------------------- compare

def _cmp_array(a: np.ndarray, b: np.ndarray, atol, rtol):
    """Return (bitwise identical, within tolerance, max absolute difference, max relative difference)."""
    if a.shape != b.shape or a.dtype.kind != b.dtype.kind:
        return False, False, float("inf"), float("inf")
    if a.dtype.kind in "US":
        same = bool(np.array_equal(a, b))
        return same, same, 0.0 if same else float("inf"), 0.0 if same else float("inf")
    same = bool(np.array_equal(a, b, equal_nan=True))
    if same:
        return True, True, 0.0, 0.0
    af = np.asarray(a, dtype=np.complex128 if a.dtype.kind == "c" else np.float64)
    bf = np.asarray(b, dtype=np.complex128 if b.dtype.kind == "c" else np.float64)
    d = np.abs(af - bf)
    finite = np.isfinite(d)
    mad = float(d[finite].max()) if finite.any() else float("inf")
    scale = np.maximum(np.abs(af), np.abs(bf))
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(scale > 0, d / scale, 0.0)
    mrd = float(rel[finite].max()) if finite.any() else float("inf")
    close = bool(np.allclose(af, bf, atol=atol or 0.0, rtol=rtol or 0.0, equal_nan=True)) if (atol or rtol) else False
    return False, close, mad, mrd


def compare(tag_a: str, tag_b: str, verbose: bool) -> int:
    da, db = GOLDEN / tag_a, GOLDEN / tag_b
    tol = {c[0]: (c[2], c[3]) for c in CASES}
    names = sorted({p.stem for p in da.glob("*.npz")} | {p.stem for p in db.glob("*.npz")})
    if not names:
        print(f"GOLDEN_DIFF 0  (no case files under {da} or {db})")
        return 2
    n_diff = 0
    for name in names:
        fa, fb = da / f"{name}.npz", db / f"{name}.npz"
        if not (fa.exists() and fb.exists()):
            n_diff += 1
            print(f"== {name}: present only in {'A' if fa.exists() else 'B'}")
            continue
        za, zb = np.load(fa, allow_pickle=False), np.load(fb, allow_pickle=False)
        ma, mb = json.loads(za["__meta"].item()), json.loads(zb["__meta"].item())
        atol, rtol = tol.get(name, (None, None))
        keys = sorted((set(za.files) | set(zb.files)) - {"__meta"})
        bad = []
        lines = []
        for k in keys:
            if k not in za.files or k not in zb.files:
                bad.append(k)
                lines.append(f"   {k}: present only in {'A' if k in za.files else 'B'}")
                continue
            a, b = za[k], zb[k]
            same, close, mad, mrd = _cmp_array(a, b, atol, rtol)
            if same:
                mark = "="
            elif close:
                mark = "~"
            else:
                mark = "X"
                bad.append(k)
            if verbose or mark != "=":
                lines.append(f"   {mark} {k:40s} shape={a.shape} maxabs={mad:.3e} maxrel={mrd:.3e}")
        meta_keys = ("steps_run", "triggered", "stop", "source_end", "shape", "num_time_steps", "calls")
        meta_diff = [k for k in meta_keys if ma.get(k) != mb.get(k)]
        if meta_diff:
            bad.extend(f"meta:{k}" for k in meta_diff)
            lines.append("   meta differs: " + "; ".join(f"{k}: {ma.get(k)} vs {mb.get(k)}" for k in meta_diff))
        n_same = sum(1 for k in keys if k in za.files and k in zb.files and _cmp_array(za[k], zb[k], None, None)[0])
        status = "DIFF" if bad else "SAME"
        n_diff += bool(bad)
        print(f"== {name}: {status}  arrays {len(keys)}, bitwise identical {n_same}, steps {ma.get('steps_run')} vs {mb.get('steps_run')}"
              f"{'' if not (atol or rtol) else f'  (tolerance atol={atol} rtol={rtol})'}")
        for ln in lines:
            print(ln)
    print("GOLDEN_SAME" if n_diff == 0 else f"GOLDEN_DIFF {n_diff}")
    return 0 if n_diff == 0 else 1


# ----------------------------------------------------------------------------- cli

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="run every case, storing results under <output dir>/<tag>/")
    p.add_argument("tag")
    p.add_argument("--only", default="", help="comma-separated case names")
    p = sub.add_parser("compare", help="compare two tags")
    p.add_argument("tag_a")
    p.add_argument("tag_b")
    p.add_argument("-v", "--verbose", action="store_true", help="print one line per array")
    sub.add_parser("list", help="list the cases")
    args = ap.parse_args(argv)
    if args.cmd == "list":
        for name, ctor, atol, rtol in CASES:
            print(f"{name:16s} {ctor.__doc__.strip().splitlines()[0]}" + (f"  [atol={atol} rtol={rtol}]" if (atol or rtol) else ""))
        return 0
    if args.cmd == "run":
        return run(args.tag, [s for s in args.only.split(",") if s])
    return compare(args.tag_a, args.tag_b, args.verbose)


if __name__ == "__main__":
    sys.exit(main())
