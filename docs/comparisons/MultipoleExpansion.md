# MultipoleExpansion

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/MultipoleExpansion/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 11 on the reference side, 11 on ours, 0 that did not line up |
| Numbers compared | 9 |
| Largest relative difference | 99.27% |
| Over 5% / over 20% | 1 / 1 |
| Figures | 5 on the reference side, 5 on ours |
| Largest pixel difference | 5.78% |

Machine-readable form of everything below: [`data/MultipoleExpansion.json`](data/MultipoleExpansion.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 1 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 4 | 0.13767145817929347 | 0.137671 | 0.001 | 99.27% |  |
| 7 | slice 1 done | 1 | 1 | 0 |  |
| 7 | slice 2 done | 2 | 2 | 0 |  |
| 7 | slice 3 done | 3 | 3 | 0 |  |
| 7 | slice 4 done | 4 | 4 | 0 |  |
| 7 | slice 5 done | 5 | 5 | 0 |  |
| 7 | slice 6 done | 6 | 6 | 0 |  |
| 7 | slice 7 done | 7 | 7 | 0 |  |
| 7 | slice 8 done | 8 | 8 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 8 | 0 | 5.78% |
| 10 | 0 | 4.17% |
| 10 | 1 | 3.61% |
| 10 | 2 | 3.01% |
| 10 | 3 | 4.07% |
