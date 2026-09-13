# OpenEM

[English](README.md) · [简体中文](README.zh-CN.md)

A CUDA FDTD solver that takes Tidy3D `Simulation` objects as its input and solves them on one local
GPU — fp32, Yee grid, single card. After `openem.install()`, `tidy3d.web.run(sim)` and the `grad` of
tidy3d's autograd interface run on that GPU instead of the cloud; an existing notebook needs no edits.

> **The grid and the materials come from Tidy3D.** OpenEM reads a `td.Simulation` and takes the grid,
> the geometry rasterization and the material arrays **directly from what the Tidy3D client has already
> computed**. It only does the time stepping, the boundaries, the source injection and the monitor
> readout. This is deliberate: with both sides fed the same discretization, a difference in the results
> can only come from the solver itself — otherwise you cannot tell a modelling difference from a solver
> difference. A grid and material generator of our own is planned for later, and is not part of this
> release.

## Install

```bash
git clone <repository url>
cd OpenEM
pip install -e ".[gpu]"      # or: pip install -e .   if cupy is already installed
```

Requirements: Python ≥ 3.10, `tidy3d` ≥ 2.8, NumPy, SciPy, and an NVIDIA GPU with CUDA 12 or newer
plus `cupy-cuda12x` ≥ 13 (that is what the `gpu` extra pulls in). The CUDA kernels are compiled by
cupy on first use, so there is no build step and no `nvcc` needed at install time. On a machine
without a GPU the package still imports and the tests that do not need CUDA still run (that is the
configuration CI uses), but nothing can be solved.

## 30 seconds

A plane wave at normal incidence on a 0.3 µm slab of ε = 4, transmission read off a flux monitor:

```python
import tidy3d as td
import openem

openem.install()                       # from here on td.web.run solves locally

f0 = 2.998e14                          # 1 µm
periodic = td.Boundary(plus=td.Periodic(), minus=td.Periodic())
sim = td.Simulation(
    size=(0.6, 0.6, 3.0),
    grid_spec=td.GridSpec.uniform(dl=0.02),
    structures=[td.Structure(
        geometry=td.Box(center=(0, 0, 0), size=(td.inf, td.inf, 0.3)),
        medium=td.Medium(permittivity=4.0))],
    sources=[td.PlaneWave(
        center=(0, 0, -1.0), size=(td.inf, td.inf, 0), direction="+",
        source_time=td.GaussianPulse(freq0=f0, fwidth=0.2 * f0))],
    monitors=[td.FluxMonitor(center=(0, 0, 1.0), size=(td.inf, td.inf, 0),
                             freqs=[f0], name="T")],
    run_time=4e-13,
    boundary_spec=td.BoundarySpec(x=periodic, y=periodic, z=td.Boundary.pml()))

sim_data = td.web.run(sim, task_name="slab")
print("transmission =", float(sim_data["T"].flux.values[0]))
```

```
transmission = 0.8310044973694679
```

The analytic single-layer thin-film formula gives 0.8373 for this slab. `examples/01_slab_planewave.py`
is the same simulation over 11 frequencies and prints the reflection and the energy balance as well;
`docs/tutorial.md` walks through it, a mode source and a gradient.

`openem.install()` patches `web.Batch` and `web.Job` the same way, so a notebook that submits batches
keeps working, and the cost queries return a placeholder instead of talking to the cloud. If you would
rather not patch anything, `openem.run(sim)` solves a `Simulation` — or a list, or a `{name: sim}`
dict — directly and hands back a real `td.SimulationData`; that path keeps its intermediate files on
disk and reuses them, which the tutorial explains.

## What it covers

| | |
|---|---|
| **Media** | `Medium` (with `conductivity`), anything that reduces to `PoleResidue` (`Lorentz`, `Drude`, `Sellmeier`, `Debye`, and their `Custom*` spatially varying forms), `CustomMedium`, `AnisotropicMedium` (diagonal), `FullyAnisotropicMedium` (full tensor, staircased), `PECMedium` on structures, `Medium2D` through the client's volumetric equivalent, and time-modulated media (`modulation_spec`) |
| **Boundaries** | `PML` / `StablePML` (CPML), `Absorber`, `Periodic`, `BlochBoundary`, `PECBoundary`, `PMCBoundary`, and mirror `symmetry` planes |
| **Sources** | `PlaneWave` (normal and oblique), `ModeSource`, `PointDipole`, `UniformCurrentSource`, `GaussianBeam`, `TFSF`, `CustomFieldSource`, `CustomCurrentSource` |
| **Monitors** | `FluxMonitor` (plane or box), `FluxTimeMonitor`, `FieldMonitor`, `FieldTimeMonitor`, `ModeMonitor`, `ModeSolverMonitor`, `PermittivityMonitor`, `DiffractionMonitor`, and field projection monitors (angle domain and k-space) |
| **Gradients** | the discrete adjoint behind tidy3d's autograd interface: permittivity, `CustomMedium` pixels, and geometry parameters of `Box` / `Cylinder` / `PolySlab` / `GeometryGroup` |

Everything outside this list fails closed — the scene builder raises `NotImplementedError` rather than
quietly substituting something else.

## Validation

- **128 notebooks** from the Tidy3D example library, each run end to end and compared number by number
  and figure by figure with the reference output: **99 agree**, **11 agree within a wider band**, and
  **18 differ only in the trajectory** a non-convex optimizer walks (forward results and the first-step
  gradient still line up). The per-notebook table, with the actual numbers, is in
  [`docs/validation.md`](docs/validation.md); 11 further notebooks are held back and the reasons are
  listed there too.
- **22 golden simulations** (`tools/golden.py`) compared **bit for bit** against the baseline. Only the
  flux-objective gradient carries a tolerance, for the reason under Limitations.
- The test suite covers the solver piece by piece (`pytest tests`); the CUDA tests skip on a machine
  without a GPU.
- Timing against the reference implementation, 118 simulations, split into preparation / stepping /
  readout: [`docs/performance.md`](docs/performance.md).
- What "agree" means, and the band behind each verdict: the table at the top of
  [`docs/validation.md`](docs/validation.md).

## Limitations

- **One GPU.** No multi-GPU and no multi-node. A simulation has to fit on a single card;
  `install()` refuses anything above 80 M cells by default, and `install(max_cells=...)` moves that
  ceiling.
- **No grid and no material generation** — both come from the Tidy3D client, as described at the top.
  Subpixel averaging needs `tidy3d-extras`; where it is unavailable the scene builder falls back to
  staircasing and prints that it did.
- **Fail closed.** Unsupported physics raises `NotImplementedError`. There is no silent degradation,
  because a silent fallback makes every later discrepancy impossible to attribute.
- **One known source of non-determinism.** The adjoint source of a flux objective is a
  `CustomCurrentSource` that expands into several thousand point dipoles; the batched injection kernel
  accumulates with `atomicAdd`, and float addition is not associative, so the arrival order changes the
  last bits. Gradients of that kind can only be compared to a tolerance (measured 4e-7 to 3.7e-6
  relative over five runs). Mode-objective gradients, and all forward data, are bitwise reproducible.
- fp32 throughout.
- Code comments and docstrings are English; the pages under `docs/` are still Chinese.

## Documentation

| | |
|---|---|
| [`docs/tutorial.md`](docs/tutorial.md) | three worked steps: slab transmission against the analytic formula, a waveguide mode source and mode monitor, a gradient through autograd |
| [`docs/architecture.md`](docs/architecture.md) | directory tree, one line per module, and the import boundaries the tests lock down |
| [`docs/validation.md`](docs/validation.md) | the 128-notebook comparison, one row per notebook |
| [`docs/comparisons/`](docs/comparisons/) | the same comparison in full: every one of the 6,003 paired numbers and 1,125 figures, one page and one JSON file per notebook |
| [`docs/performance.md`](docs/performance.md) | per-example timing against the reference implementation |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | the rule for changing solver code: the golden set must stay bit for bit identical |

## License and citation

GPL-3.0-or-later — see [LICENSE](LICENSE). Derivative works must also be released under the GPL. If you use OpenEM in your work, please cite it with the metadata
in [CITATION.cff](CITATION.cff).
