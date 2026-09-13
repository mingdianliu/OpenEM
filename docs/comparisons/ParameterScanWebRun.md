# ParameterScanWebRun

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/ParameterScanWebRun/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Tutorial |
| Verdict | **agree** |
| Reference output | official archived output |
| Cells | 23 on the reference side, 23 on ours, 4 that did not line up |
| Numbers compared | 30 of 32 paired; the other 2 are near-zero, see below |
| Largest relative difference | 680% |
| Over 5% / over 20% | 2 / 2 |
| Figures | 5 on the reference side, 5 on ours |
| Largest pixel difference | 9.32% |

Machine-readable form of everything below: [`data/ParameterScanWebRun.json`](data/ParameterScanWebRun.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 2 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 10 | amplitude^2 = 0.00 | 2 | 2 | 0 |  |
| 10 | amplitude^2 = 0.00 | 0 | 0 | 0 | near-zero |
| 10 | phase       = 0.05 (rad) | 0.05 | 0.39 | 680% |  |
| 10 | amplitude^2 = 0.00 | 2 | 2 | 0 |  |
| 10 | amplitude^2 = 0.00 | 0 | 0 | 0 | near-zero |
| 10 | phase       = -2.49 (rad) | -2.49 | 2.81 | 213% |  |
| 10 | amplitude^2 = 0.05 | 2 | 2 | 0 |  |
| 10 | amplitude^2 = 0.05 | 0.05 | 0.05 | 0 |  |
| 10 | phase       = -2.58 (rad) | -2.58 | -2.68 | 3.88% |  |
| 10 | amplitude^2 = 0.94 | 2 | 2 | 0 |  |
| 10 | amplitude^2 = 0.94 | 0.94 | 0.95 | 1.06% |  |
| 10 | phase       = 2.13 (rad) | 2.13 | 2.03 | 4.69% |  |
| 12 | (4, 11) | 4 | 4 | 0 |  |
| 12 | (4, 11) | 11 | 11 | 0 |  |
| 21 | [[{'wg_spacing_coup': 0.1, | 0.1 | 0.1 | 0 |  |
| 21 | 'coup_length': 6.0, | 6 | 6 | 0 |  |
| 21 | {'wg_spacing_coup': 0.1, | 0.1 | 0.1 | 0 |  |
| 21 | 'coup_length': 8.0, | 8 | 8 | 0 |  |
| 21 | {'wg_spacing_coup': 0.1, | 0.1 | 0.1 | 0 |  |
| 21 | 'coup_length': 10.0, | 10 | 10 | 0 |  |
| 21 | [{'wg_spacing_coup': 0.15, | 0.15 | 0.15 | 0 |  |
| 21 | 'coup_length': 6.5, | 6.5 | 6.5 | 0 |  |
| 21 | {'wg_spacing_coup': 0.15, | 0.15 | 0.15 | 0 |  |
| 21 | 'coup_length': 8.5, | 8.5 | 8.5 | 0 |  |
| 21 | {'wg_spacing_coup': 0.15, | 0.15 | 0.15 | 0 |  |
| 21 | 'coup_length': 10.5, | 10.5 | 10.5 | 0 |  |
| 21 | [{'wg_spacing_coup': 0.2, | 0.2 | 0.2 | 0 |  |
| 21 | 'coup_length': 7.0, | 7 | 7 | 0 |  |
| 21 | {'wg_spacing_coup': 0.2, | 0.2 | 0.2 | 0 |  |
| 21 | 'coup_length': 9.0, | 9 | 9 | 0 |  |
| 21 | {'wg_spacing_coup': 0.2, | 0.2 | 0.2 | 0 |  |
| 21 | 'coup_length': 11.0, | 11 | 11 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 5 | 0 | 4.83% |
| 11 | 0 | 6.71% |
| 14 | 0 | 5.13% |
| 18 | 0 | 9.32% |
| 22 | 0 | 5.65% |
