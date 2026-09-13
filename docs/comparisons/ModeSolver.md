# ModeSolver

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/ModeSolver/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Tutorial |
| Verdict | **agree** |
| Reference output | official archived output |
| Cells | 34 on the reference side, 34 on ours, 1 that did not line up |
| Numbers compared | 12 of 13 paired; the other 1 are near-zero, see below |
| Largest relative difference | 1.13% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 19 on the reference side, 19 on ours |
| Largest pixel difference | 13.77% |

Machine-readable form of everything below: [`data/ModeSolver.json`](data/ModeSolver.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 9 | first mode effective index at freq0: n_eff = 1.77, k_eff = 0.00e+00 | 1.77 | 1.75 | 1.13% |  |
| 9 | first mode effective index at freq0: n_eff = 1.77, k_eff = 0.00e+00 | 0 | 0 | 0 | near-zero |
| 25 | Solving for resolution = 10 | 10 | 10 | 0 |  |
| 25 | Solving for resolution = 12 | 12 | 12 | 0 |  |
| 25 | Solving for resolution = 14 | 14 | 14 | 0 |  |
| 25 | Solving for resolution = 16 | 16 | 16 | 0 |  |
| 25 | Solving for resolution = 18 | 18 | 18 | 0 |  |
| 25 | Solving for resolution = 20 | 20 | 20 | 0 |  |
| 25 | Solving for resolution = 22 | 22 | 22 | 0 |  |
| 25 | Solving for resolution = 24 | 24 | 24 | 0 |  |
| 25 | Solving for resolution = 26 | 26 | 26 | 0 |  |
| 25 | Solving for resolution = 28 | 28 | 28 | 0 |  |
| 25 | Solving for resolution = 30 | 30 | 30 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 2 | 0 | 3.02% |
| 8 | 0 | 5.14% |
| 10 | 0 | 4.55% |
| 11 | 0 | 5.43% |
| 12 | 0 | 7.86% |
| 13 | 0 | 5.69% |
| 15 | 0 | 3.58% |
| 17 | 0 | 6.37% |
| 20 | 0 | 4.29% |
| 21 | 0 | 6.20% |
| 22 | 0 | 5.34% |
| 24 | 0 | 4.25% |
| 26 | 0 | 13.77% |
| 27 | 0 | 3.11% |
| 29 | 0 | 10.86% |
| 30 | 0 | 4.96% |
| 32 | 0 | 11.08% |
| 32 | 1 | 5.20% |
| 33 | 0 | 10.89% |
