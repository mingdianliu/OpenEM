# Autograd0QuickstartII

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd0QuickstartII/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **agree** |
| Reference output | official archived output |
| Cells | 6 on the reference side, 6 on ours, 0 that did not line up |
| Numbers compared | 14 |
| Largest relative difference | 0 |
| Over 5% / over 20% | 0 / 0 |
| Figures | 1 on the reference side, 1 on ours |
| Largest pixel difference | 3.52% |

Machine-readable form of everything below: [`data/Autograd0QuickstartII.json`](data/Autograd0QuickstartII.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 3 | step = 1 | 1 | 1 | 0 |  |
| 3 | step = 2 | 2 | 2 | 0 |  |
| 3 | step = 3 | 3 | 3 | 0 |  |
| 3 | step = 4 | 4 | 4 | 0 |  |
| 3 | step = 5 | 5 | 5 | 0 |  |
| 3 | step = 6 | 6 | 6 | 0 |  |
| 3 | step = 7 | 7 | 7 | 0 |  |
| 4 | step = 1 | 1 | 1 | 0 |  |
| 4 | step = 2 | 2 | 2 | 0 |  |
| 4 | step = 3 | 3 | 3 | 0 |  |
| 4 | step = 4 | 4 | 4 | 0 |  |
| 4 | step = 5 | 5 | 5 | 0 |  |
| 4 | step = 6 | 6 | 6 | 0 |  |
| 4 | step = 7 | 7 | 7 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 5 | 0 | 3.52% |
