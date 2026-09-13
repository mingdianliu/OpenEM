# BoundaryConditions

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/BoundaryConditions/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Tutorial |
| Verdict | **agree** |
| Reference output | official archived output |
| Cells | 17 on the reference side, 17 on ours, 0 that did not line up |
| Numbers compared | none printed |
| Largest relative difference | - |
| Over 5% / over 20% | 0 / 0 |
| Figures | 14 on the reference side, 14 on ours |
| Largest pixel difference | 37.87% |

Machine-readable form of everything below: [`data/BoundaryConditions.json`](data/BoundaryConditions.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 3 | 0 | 5.37% |
| 4 | 0 | 5.37% |
| 5 | 0 | 10.97% |
| 6 | 0 | 3.45% |
| 8 | 0 | 7.32% |
| 10 | 0 | 8.62% |
| 11 | 0 | 8.37% |
| 13 | 0 | 24.83% |
| 13 | 1 | 29.24% |
| 13 | 2 | 37.87% |
| 14 | 0 | 8.58% |
| 16 | 0 | 15.62% |
| 16 | 1 | 14.93% |
| 16 | 2 | 19.86% |
