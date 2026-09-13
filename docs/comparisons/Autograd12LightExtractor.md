# Autograd12LightExtractor

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd12LightExtractor/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **trajectory divergence** |
| Reference output | official archived output |
| Cells | 30 on the reference side, 30 on ours, 0 that did not line up |
| Numbers compared | 5 of 6 paired; the other 1 are near-zero, see below |
| Largest relative difference | 94.27% |
| Over 5% / over 20% | 2 / 2 |
| Figures | 11 on the reference side, 11 on ours |
| Largest pixel difference | 53.01% |

Machine-readable form of everything below: [`data/Autograd12LightExtractor.json`](data/Autograd12LightExtractor.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 2 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 15 | Found 100 iterations previously completed out of 100 total. | 100 | 100 | 0 |  |
| 15 | Found 100 iterations previously completed out of 100 total. | 100 | 100 | 0 |  |
| 22 | Text(0.5, 0, 'iteration') | 0.5 | 0.5 | 0 |  |
| 22 | Text(0.5, 0, 'iteration') | 0 | 0 | 0 | near-zero |
| 28 | 0.9658721904728587 | 0.965872 | 0.0569368 | 94.11% |  |
| 28 | 0.8196551325491168 | 0.819655 | 0.0470058 | 94.27% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 8 | 0 | 6.41% |
| 10 | 0 | 4.86% |
| 12 | 0 | 2.95% |
| 19 | 0 | 53.01% |
| 20 | 0 | 45.99% |
| 21 | 0 | 16.22% |
| 22 | 0 | 11.01% |
| 23 | 0 | 38.97% |
| 24 | 0 | 13.38% |
| 26 | 0 | 19.65% |
| 27 | 0 | 4.65% |
