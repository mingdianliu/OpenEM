# GratingEfficiency

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/GratingEfficiency/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Tutorial |
| Verdict | **agree** |
| Reference output | official archived output |
| Cells | 17 on the reference side, 17 on ours, 2 that did not line up |
| Numbers compared | 19 of 25 paired; the other 6 are near-zero, see below |
| Largest relative difference | 0.13% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 7 on the reference side, 7 on ours |
| Largest pixel difference | 8.78% |

Machine-readable form of everything below: [`data/GratingEfficiency.json`](data/GratingEfficiency.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 6 | Total power: 0.9984229082058051 | 0.998423 | 0.999705 | 0.13% |  |
| 6 | Theta (degrees): 63.09, 41.98, 26.48, 12.88, 0.00, 12.88, 26.48, 41.98, 63.09 | 63.09 | 63.09 | 0 |  |
| 6 | Theta (degrees): 63.09, 41.98, 26.48, 12.88, 0.00, 12.88, 26.48, 41.98, 63.09 | 41.98 | 41.98 | 0 |  |
| 6 | Theta (degrees): 63.09, 41.98, 26.48, 12.88, 0.00, 12.88, 26.48, 41.98, 63.09 | 26.48 | 26.48 | 0 |  |
| 6 | Theta (degrees): 63.09, 41.98, 26.48, 12.88, 0.00, 12.88, 26.48, 41.98, 63.09 | 12.88 | 12.88 | 0 |  |
| 6 | Theta (degrees): 63.09, 41.98, 26.48, 12.88, 0.00, 12.88, 26.48, 41.98, 63.09 | 0 | 0 | 0 | near-zero |
| 6 | Theta (degrees): 63.09, 41.98, 26.48, 12.88, 0.00, 12.88, 26.48, 41.98, 63.09 | 12.88 | 12.88 | 0 |  |
| 6 | Theta (degrees): 63.09, 41.98, 26.48, 12.88, 0.00, 12.88, 26.48, 41.98, 63.09 | 26.48 | 26.48 | 0 |  |
| 6 | Theta (degrees): 63.09, 41.98, 26.48, 12.88, 0.00, 12.88, 26.48, 41.98, 63.09 | 41.98 | 41.98 | 0 |  |
| 6 | Theta (degrees): 63.09, 41.98, 26.48, 12.88, 0.00, 12.88, 26.48, 41.98, 63.09 | 63.09 | 63.09 | 0 |  |
| 6 | [[[ 0.00000000e+00+0.00000000e+00j, | 0 | 0 | 0 | near-zero |
| 6 | [[[ 0.00000000e+00+0.00000000e+00j, | 0 | 0 | 0 | near-zero |
| 6 | [[[ 0.00000000e+00+0.00000000e+00j, | 0 | 0 | 0 | near-zero |
| 6 | [[[ 0.00000000e+00+0.00000000e+00j, | 0 | 0 | 0 | near-zero |
| 6 | [[[ 0.00000000e+00+0.00000000e+00j, | 0 | 0 | 0 | near-zero |
| 13 | Total power: 1.0007346889466473 | 1.00073 | 1.00014 | 0.06% |  |
| 13 | Theta (degrees): 59.97, 41.10, 27.56, 18.00, 16.04, 23.30, 35.43, 51.66, 85.40 | 59.97 | 59.97 | 0 |  |
| 13 | Theta (degrees): 59.97, 41.10, 27.56, 18.00, 16.04, 23.30, 35.43, 51.66, 85.40 | 41.1 | 41.1 | 0 |  |
| 13 | Theta (degrees): 59.97, 41.10, 27.56, 18.00, 16.04, 23.30, 35.43, 51.66, 85.40 | 27.56 | 27.56 | 0 |  |
| 13 | Theta (degrees): 59.97, 41.10, 27.56, 18.00, 16.04, 23.30, 35.43, 51.66, 85.40 | 18 | 18 | 0 |  |
| 13 | Theta (degrees): 59.97, 41.10, 27.56, 18.00, 16.04, 23.30, 35.43, 51.66, 85.40 | 16.04 | 16.04 | 0 |  |
| 13 | Theta (degrees): 59.97, 41.10, 27.56, 18.00, 16.04, 23.30, 35.43, 51.66, 85.40 | 23.3 | 23.3 | 0 |  |
| 13 | Theta (degrees): 59.97, 41.10, 27.56, 18.00, 16.04, 23.30, 35.43, 51.66, 85.40 | 35.43 | 35.43 | 0 |  |
| 13 | Theta (degrees): 59.97, 41.10, 27.56, 18.00, 16.04, 23.30, 35.43, 51.66, 85.40 | 51.66 | 51.66 | 0 |  |
| 13 | Theta (degrees): 59.97, 41.10, 27.56, 18.00, 16.04, 23.30, 35.43, 51.66, 85.40 | 85.4 | 85.4 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 4 | 0 | 4.86% |
| 8 | 0 | 6.77% |
| 9 | 0 | 8.62% |
| 10 | 0 | 6.31% |
| 11 | 0 | 5.09% |
| 14 | 0 | 8.78% |
| 16 | 0 | 7.11% |
