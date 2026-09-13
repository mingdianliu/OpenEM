# WaveguidePluginDemonstration

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/WaveguidePluginDemonstration/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Tutorial |
| Verdict | **agree** |
| Reference output | official archived output |
| Cells | 24 on the reference side, 24 on ours, 4 that did not line up |
| Numbers compared | 32 |
| Largest relative difference | 3.92e-05% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 18 on the reference side, 18 on ours |
| Largest pixel difference | 7.20% |

Machine-readable form of everything below: [`data/WaveguidePluginDemonstration.json`](data/WaveguidePluginDemonstration.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 4 | Effective indices: [[2.53287852 1.89256577]] | 2.53288 | 2.53288 | 0 |  |
| 4 | Effective indices: [[2.53287852 1.89256577]] | 1.89257 | 1.89257 | 1.06e-06% |  |
| 14 | Solving for resolution = 15 | 15 | 15 | 0 |  |
| 14 | Solving for resolution = 18 | 18 | 18 | 0 |  |
| 14 | Solving for resolution = 21 | 21 | 21 | 0 |  |
| 14 | Solving for resolution = 24 | 24 | 24 | 0 |  |
| 14 | Solving for resolution = 27 | 27 | 27 | 0 |  |
| 14 | Solving for resolution = 30 | 30 | 30 | 0 |  |
| 14 | Solving for resolution = 33 | 33 | 33 | 0 |  |
| 14 | Solving for resolution = 36 | 36 | 36 | 0 |  |
| 14 | Solving for resolution = 39 | 39 | 39 | 0 |  |
| 14 | Solving for resolution = 42 | 42 | 42 | 0 |  |
| 14 | Solving for resolution = 45 | 45 | 45 | 0 |  |
| 17 | Curvature loss: 14.9 dB/cm | 14.9 | 14.9 | 0 |  |
| 18 | Solving for radius = 5.00 | 5 | 5 | 0 |  |
| 18 | Solving for radius = 6.46 | 6.46 | 6.46 | 0 |  |
| 18 | Solving for radius = 8.34 | 8.34 | 8.34 | 0 |  |
| 18 | Solving for radius = 10.77 | 10.77 | 10.77 | 0 |  |
| 18 | Solving for radius = 13.91 | 13.91 | 13.91 | 0 |  |
| 18 | Solving for radius = 17.97 | 17.97 | 17.97 | 0 |  |
| 18 | Solving for radius = 23.21 | 23.21 | 23.21 | 0 |  |
| 18 | Solving for radius = 29.97 | 29.97 | 29.97 | 0 |  |
| 18 | Solving for radius = 38.71 | 38.71 | 38.71 | 0 |  |
| 18 | Solving for radius = 50.00 | 50 | 50 | 0 |  |
| 21 | Effective indices: [[2.58790025 2.55348271 1.98156471 1.85005072]] | 2.5879 | 2.5879 | 0 |  |
| 21 | Effective indices: [[2.58790025 2.55348271 1.98156471 1.85005072]] | 2.55348 | 2.55348 | 0 |  |
| 21 | Effective indices: [[2.58790025 2.55348271 1.98156471 1.85005072]] | 1.98156 | 1.98156 | 5.05e-07% |  |
| 21 | Effective indices: [[2.58790025 2.55348271 1.98156471 1.85005072]] | 1.85005 | 1.85005 | 1.08e-06% |  |
| 21 | Effective mode areas (µm²): [[0.37130773 0.35534603 0.63694388 0.61225581]] | 0.371308 | 0.371308 | 0 |  |
| 21 | Effective mode areas (µm²): [[0.37130773 0.35534603 0.63694388 0.61225581]] | 0.355346 | 0.355346 | 0 |  |
| 21 | Effective mode areas (µm²): [[0.37130773 0.35534603 0.63694388 0.61225581]] | 0.636944 | 0.636944 | 2.67e-05% |  |
| 21 | Effective mode areas (µm²): [[0.37130773 0.35534603 0.63694388 0.61225581]] | 0.612256 | 0.612256 | 3.92e-05% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 2 | 0 | 3.16% |
| 3 | 0 | 5.72% |
| 4 | 0 | 5.40% |
| 5 | 0 | 3.44% |
| 6 | 0 | 3.82% |
| 7 | 0 | 3.66% |
| 8 | 0 | 4.77% |
| 10 | 0 | 2.79% |
| 11 | 0 | 4.19% |
| 12 | 0 | 2.80% |
| 13 | 0 | 6.29% |
| 15 | 0 | 7.20% |
| 16 | 0 | 4.13% |
| 19 | 0 | 6.71% |
| 20 | 0 | 2.73% |
| 21 | 0 | 5.42% |
| 22 | 0 | 4.43% |
| 23 | 0 | 4.56% |
