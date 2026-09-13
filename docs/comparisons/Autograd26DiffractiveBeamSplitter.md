# Autograd26DiffractiveBeamSplitter

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd26DiffractiveBeamSplitter/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **trajectory divergence** |
| Reference output | official archived output |
| Cells | 21 on the reference side, 21 on ours, 0 that did not line up |
| Numbers compared | 11 |
| Largest relative difference | 5,179% |
| Over 5% / over 20% | 3 / 1 |
| Figures | 5 on the reference side, 6 on ours |
| Largest pixel difference | 29.56% |

Machine-readable form of everything below: [`data/Autograd26DiffractiveBeamSplitter.json`](data/Autograd26DiffractiveBeamSplitter.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 3 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 10 | Reference power (power0): 0.9675 | 0.9675 | 0.9676 | 0.01% |  |
| 13 | Loss value (maximize): 0.8790 | 0.879 | 0.8789 | 0.01% |  |
| 13 | Fab penalty: 9.9816e-01 | 0.99816 | 0.99816 | 0 |  |
| 13 | Efficiency: 0.9999 | 0.9999 | 0.9999 | 0 |  |
| 14 | Gradient norm: 1.249e-03 | 0.001249 | 0.001174 | 6.00% |  |
| 14 | Gradient range: [-1.059e-05, 5.702e-06] | -1.059e-05 | -1.015e-05 | 4.15% |  |
| 14 | Gradient range: [-1.059e-05, 5.702e-06] | 5.702e-06 | 5.862e-06 | 2.81% |  |
| 15 | Efficiency: 1.00 | 1 | 1 | 0 |  |
| 15 | RMSE: 13.68 | 13.68 | 13.68 | 0 |  |
| 20 | Efficiency: 0.89 | 0.89 | 0.94 | 5.62% |  |
| 20 | RMSE: 0.14 | 0.14 | 7.39 | 5,179% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 8 | 0 | 4.09% |
| 15 | 0 | 7.88% |
| 17 | 0 | 22.12% |
| 18 | 0 | 11.11% |
| 20 | 0 | 29.56% |
