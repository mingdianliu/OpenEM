# PhaseChangeAntennas

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/PhaseChangeAntennas/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 16 on the reference side, 16 on ours, 0 that did not line up |
| Numbers compared | 14 of 16 paired; the other 2 are near-zero, see below |
| Largest relative difference | 0.72% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 4 on the reference side, 4 on ours |
| Largest pixel difference | 25.01% |

Machine-readable form of everything below: [`data/PhaseChangeAntennas.json`](data/PhaseChangeAntennas.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 10 | ScatteringCylinder_0: -0.00 | 0 | 0 | 0 | near-zero |
| 10 | ScatteringCylinder_1: 38.03 | 38.03 | 38.01 | 0.05% |  |
| 10 | ScatteringCylinder_2: 107.05 | 107.05 | 107.1 | 0.05% |  |
| 10 | ScatteringCylinder_3: 162.07 | 162.07 | 162.09 | 0.01% |  |
| 10 | ScatteringCylinder_4: 180.00 | 180 | 180 | 0 |  |
| 10 | ScatteringCylinder_5: 218.03 | 218.03 | 218.01 | 9.17e-03% |  |
| 10 | ScatteringCylinder_6: 287.05 | 287.05 | 287.1 | 0.02% |  |
| 10 | ScatteringCylinder_7: 342.07 | 342.07 | 342.09 | 5.85e-03% |  |
| 15 | ScatteringCylinder_0: -0.00 | 0 | 0 | 0 | near-zero |
| 15 | ScatteringCylinder_1: 37.99 | 37.99 | 37.89 | 0.26% |  |
| 15 | ScatteringCylinder_2: 107.35 | 107.35 | 108.12 | 0.72% |  |
| 15 | ScatteringCylinder_3: 162.19 | 162.19 | 162.65 | 0.28% |  |
| 15 | ScatteringCylinder_4: 180.00 | 180 | 180 | 0 |  |
| 15 | ScatteringCylinder_5: 217.99 | 217.99 | 217.89 | 0.05% |  |
| 15 | ScatteringCylinder_6: 287.35 | 287.35 | 288.12 | 0.27% |  |
| 15 | ScatteringCylinder_7: 342.19 | 342.19 | 342.65 | 0.13% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 5 | 0 | 1.29% |
| 8 | 0 | 20.86% |
| 10 | 0 | 24.87% |
| 15 | 0 | 25.01% |
