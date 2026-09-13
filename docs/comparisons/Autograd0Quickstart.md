# Autograd0Quickstart

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd0Quickstart/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **agree** |
| Reference output | official archived output |
| Cells | 5 on the reference side, 5 on ours, 0 that did not line up |
| Numbers compared | 21 |
| Largest relative difference | 1.59% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 1 on the reference side, 1 on ours |
| Largest pixel difference | 4.60% |

Machine-readable form of everything below: [`data/Autograd0Quickstart.json`](data/Autograd0Quickstart.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 3 | Iteration = 1 | 1 | 1 | 0 |  |
| 3 | size_box = 2.68 | 2.68 | 2.68 | 0 |  |
| 3 | intensity = 861 | 861 | 860 | 0.12% |  |
| 3 | Iteration = 2 | 2 | 2 | 0 |  |
| 3 | size_box = 2.89 | 2.89 | 2.9 | 0.35% |  |
| 3 | intensity = 1068 | 1068 | 1072 | 0.37% |  |
| 3 | Iteration = 3 | 3 | 3 | 0 |  |
| 3 | size_box = 3.07 | 3.07 | 3.07 | 0 |  |
| 3 | intensity = 1258 | 1258 | 1261 | 0.24% |  |
| 3 | Iteration = 4 | 4 | 4 | 0 |  |
| 3 | size_box = 3.25 | 3.25 | 3.26 | 0.31% |  |
| 3 | intensity = 1402 | 1402 | 1409 | 0.50% |  |
| 3 | Iteration = 5 | 5 | 5 | 0 |  |
| 3 | size_box = 3.72 | 3.72 | 3.72 | 0 |  |
| 3 | intensity = 1703 | 1703 | 1730 | 1.59% |  |
| 3 | Iteration = 6 | 6 | 6 | 0 |  |
| 3 | size_box = 3.90 | 3.9 | 3.89 | 0.26% |  |
| 3 | intensity = 2233 | 2233 | 2232 | 0.04% |  |
| 3 | Iteration = 7 | 7 | 7 | 0 |  |
| 3 | size_box = 4.08 | 4.08 | 4.09 | 0.25% |  |
| 3 | intensity = 2432 | 2432 | 2429 | 0.12% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 4 | 0 | 4.60% |
