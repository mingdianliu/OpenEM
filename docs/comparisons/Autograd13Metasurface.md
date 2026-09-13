# Autograd13Metasurface

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd13Metasurface/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **agree** |
| Reference output | official archived output |
| Cells | 25 on the reference side, 25 on ours, 1 that did not line up |
| Numbers compared | 10 |
| Largest relative difference | 374% |
| Over 5% / over 20% | 3 / 3 |
| Figures | 7 on the reference side, 7 on ours |
| Largest pixel difference | 22.10% |

Machine-readable form of everything below: [`data/Autograd13Metasurface.json`](data/Autograd13Metasurface.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 3 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 13 | Average intensity of '1.81' (a.u.) measured without any device. | 1.81 | 1.81 | 0 |  |
| 17 | initial loss value = 0.246 | 0.246 | 0.246 | 0 |  |
| 17 | gradient shape = (240, 240) | 240 | 240 | 0 |  |
| 17 | gradient shape = (240, 240) | 240 | 240 | 0 |  |
| 17 | norm of gradient = 1.327e-03 | 0.001327 | 0.001327 | 0 |  |
| 19 | step = (1 / 35) | 1 | 1 | 0 |  |
| 19 | step = (1 / 35) | 35 | 35 | 0 |  |
| 19 | loss = 1.185e-01 | 0.1185 | 0.2462 | 108% |  |
| 19 | beta = 18.00 | 18 | 1 | 94.44% |  |
| 19 | \|gradient\| = 2.799e-04 | 0.0002799 | 0.001327 | 374% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 8 | 0 | 4.66% |
| 11 | 0 | 3.47% |
| 18 | 0 | 22.10% |
| 20 | 0 | 4.42% |
| 21 | 0 | 17.27% |
| 23 | 0 | 11.23% |
| 24 | 0 | 11.21% |
