# VortexMetasurface

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/VortexMetasurface/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 14 on the reference side, 14 on ours, 0 that did not line up |
| Numbers compared | 1 |
| Largest relative difference | 99.96% |
| Over 5% / over 20% | 1 / 1 |
| Figures | 5 on the reference side, 5 on ours |
| Largest pixel difference | 9.83% |

Machine-readable form of everything below: [`data/VortexMetasurface.json`](data/VortexMetasurface.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 1 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 9 | 2.464539667330026 | 2.46454 | 0.001 | 99.96% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 4 | 0 | 1.53% |
| 6 | 0 | 9.83% |
| 8 | 0 | 9.66% |
| 12 | 0 | 5.24% |
| 13 | 0 | 0 |
