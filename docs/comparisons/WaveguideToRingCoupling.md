# WaveguideToRingCoupling

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/WaveguideToRingCoupling/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **close** |
| Reference output | reference rerun |
| Cells | 21 on the reference side, 21 on ours, 0 that did not line up |
| Numbers compared | 5 |
| Largest relative difference | 378% |
| Over 5% / over 20% | 5 / 4 |
| Figures | 6 on the reference side, 6 on ours |
| Largest pixel difference | 19.73% |

Machine-readable form of everything below: [`data/WaveguideToRingCoupling.json`](data/WaveguideToRingCoupling.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 5 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 17 | max \|through difference\|: 4.867e-05 | 4.867e-05 | 7.659e-05 | 57.37% |  |
| 17 | max \|drop difference\|: 3.060e-05 | 3.06e-05 | 9.043e-05 | 196% |  |
| 18 | max backward drop power, PML extrusion: 2.075e-05 | 2.075e-05 | 1.889e-05 | 8.96% |  |
| 18 | max backward drop power, absorber: 1.226e-05 | 1.226e-05 | 5.866e-05 | 378% |  |
| 18 | max \|backward drop power difference\|: 1.320e-05 | 1.32e-05 | 4.835e-05 | 266% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 8 | 0 | 4.83% |
| 11 | 0 | 3.04% |
| 12 | 0 | 19.73% |
| 17 | 0 | 7.76% |
| 18 | 0 | 5.43% |
| 19 | 0 | 1.96% |
