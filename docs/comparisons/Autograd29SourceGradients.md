# Autograd29SourceGradients

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd29SourceGradients/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **trajectory divergence** |
| Reference output | official archived output |
| Cells | 17 on the reference side, 17 on ours, 0 that did not line up |
| Numbers compared | 91 |
| Largest relative difference | 53.68% |
| Over 5% / over 20% | 39 / 10 |
| Figures | 6 on the reference side, 6 on ours |
| Largest pixel difference | 8.54% |

Machine-readable form of everything below: [`data/Autograd29SourceGradients.json`](data/Autograd29SourceGradients.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 39 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 4 | Independent cylinders per lens: 69 | 69 | 69 | 0 |  |
| 4 | Physical cylinders per lens: 241 | 241 | 241 | 0 |  |
| 4 | focal_length = 6.00 um | 6 | 6 | 0 |  |
| 4 | bridge_distance = 2.89 um | 2.89 | 2.89 | 0 |  |
| 4 | focal_distance = 5.79 um | 5.79 | 5.79 | 0 |  |
| 11 | baseline focused power = 0.0413 | 0.0413 | 0.0396 | 4.12% |  |
| 11 | initial normalized objective = 0.637 | 0.637 | 0.601 | 5.65% |  |
| 11 | initial focused power = 0.0263 | 0.0263 | 0.0238 | 9.51% |  |
| 11 | initial gradient norm = 0.331 | 0.331 | 0.364 | 9.97% |  |
| 12 | step = 1 | 1 | 1 | 0 |  |
| 12 | normalized objective = 0.637 | 0.637 | 0.601 | 5.65% |  |
| 12 | focused power = 0.0263 | 0.0263 | 0.0238 | 9.51% |  |
| 12 | gradient norm = 0.331 | 0.331 | 0.364 | 9.97% |  |
| 12 | step = 2 | 2 | 2 | 0 |  |
| 12 | normalized objective = 1.493 | 1.493 | 1.456 | 2.48% |  |
| 12 | focused power = 0.0616 | 0.0616 | 0.0576 | 6.49% |  |
| 12 | gradient norm = 0.418 | 0.418 | 0.47 | 12.44% |  |
| 12 | step = 3 | 3 | 3 | 0 |  |
| 12 | normalized objective = 2.240 | 2.24 | 2.232 | 0.36% |  |
| 12 | focused power = 0.0924 | 0.0924 | 0.0884 | 4.33% |  |
| 12 | gradient norm = 0.643 | 0.643 | 0.655 | 1.87% |  |
| 12 | step = 4 | 4 | 4 | 0 |  |
| 12 | normalized objective = 2.765 | 2.765 | 2.687 | 2.82% |  |
| 12 | focused power = 0.1141 | 0.1141 | 0.1064 | 6.75% |  |
| 12 | gradient norm = 0.997 | 0.997 | 0.676 | 32.20% |  |
| 12 | step = 5 | 5 | 5 | 0 |  |
| 12 | normalized objective = 3.245 | 3.245 | 3.086 | 4.90% |  |
| 12 | focused power = 0.1339 | 0.1339 | 0.1222 | 8.74% |  |
| 12 | gradient norm = 0.544 | 0.544 | 0.836 | 53.68% |  |
| 12 | step = 6 | 6 | 6 | 0 |  |
| 12 | normalized objective = 3.563 | 3.563 | 3.392 | 4.80% |  |
| 12 | focused power = 0.1470 | 0.147 | 0.1343 | 8.64% |  |
| 12 | gradient norm = 0.571 | 0.571 | 0.584 | 2.28% |  |
| 12 | step = 7 | 7 | 7 | 0 |  |
| 12 | normalized objective = 3.919 | 3.919 | 3.681 | 6.07% |  |
| 12 | focused power = 0.1617 | 0.1617 | 0.1457 | 9.89% |  |
| 12 | gradient norm = 0.741 | 0.741 | 0.724 | 2.29% |  |
| 12 | step = 8 | 8 | 8 | 0 |  |
| 12 | normalized objective = 4.090 | 4.09 | 3.948 | 3.47% |  |
| 12 | focused power = 0.1688 | 0.1688 | 0.1563 | 7.41% |  |
| 12 | gradient norm = 1.379 | 1.379 | 0.705 | 48.88% |  |
| 12 | step = 9 | 9 | 9 | 0 |  |
| 12 | normalized objective = 4.407 | 4.407 | 4.175 | 5.26% |  |
| 12 | focused power = 0.1819 | 0.1819 | 0.1653 | 9.13% |  |
| 12 | gradient norm = 1.093 | 1.093 | 0.641 | 41.35% |  |
| 12 | step = 10 | 10 | 10 | 0 |  |
| 12 | normalized objective = 4.696 | 4.696 | 4.439 | 5.47% |  |
| 12 | focused power = 0.1938 | 0.1938 | 0.1757 | 9.34% |  |
| 12 | gradient norm = 0.867 | 0.867 | 0.617 | 28.84% |  |
| 12 | step = 11 | 11 | 11 | 0 |  |
| 12 | normalized objective = 4.907 | 4.907 | 4.729 | 3.63% |  |
| 12 | focused power = 0.2025 | 0.2025 | 0.1872 | 7.56% |  |
| 12 | gradient norm = 0.709 | 0.709 | 0.524 | 26.09% |  |
| 12 | step = 12 | 12 | 12 | 0 |  |
| 12 | normalized objective = 5.172 | 5.172 | 4.91 | 5.07% |  |
| 12 | focused power = 0.2134 | 0.2134 | 0.1944 | 8.90% |  |
| 12 | gradient norm = 0.708 | 0.708 | 0.447 | 36.86% |  |
| 12 | step = 13 | 13 | 13 | 0 |  |
| 12 | normalized objective = 5.301 | 5.301 | 5.063 | 4.49% |  |
| 12 | focused power = 0.2188 | 0.2188 | 0.2004 | 8.41% |  |
| 12 | gradient norm = 0.505 | 0.505 | 0.449 | 11.09% |  |
| 12 | step = 14 | 14 | 14 | 0 |  |
| 12 | normalized objective = 5.404 | 5.404 | 5.208 | 3.63% |  |
| 12 | focused power = 0.2230 | 0.223 | 0.2062 | 7.53% |  |
| 12 | gradient norm = 0.611 | 0.611 | 0.458 | 25.04% |  |
| 12 | step = 15 | 15 | 15 | 0 |  |
| 12 | normalized objective = 5.522 | 5.522 | 5.335 | 3.39% |  |
| 12 | focused power = 0.2279 | 0.2279 | 0.2112 | 7.33% |  |
| 12 | gradient norm = 0.483 | 0.483 | 0.501 | 3.73% |  |
| 12 | step = 16 | 16 | 16 | 0 |  |
| 12 | normalized objective = 5.620 | 5.62 | 5.479 | 2.51% |  |
| 12 | focused power = 0.2319 | 0.2319 | 0.2169 | 6.47% |  |
| 12 | gradient norm = 0.506 | 0.506 | 0.543 | 7.31% |  |
| 12 | step = 17 | 17 | 17 | 0 |  |
| 12 | normalized objective = 5.653 | 5.653 | 5.623 | 0.53% |  |
| 12 | focused power = 0.2333 | 0.2333 | 0.2226 | 4.59% |  |
| 12 | gradient norm = 0.638 | 0.638 | 0.569 | 10.82% |  |
| 12 | step = 18 | 18 | 18 | 0 |  |
| 12 | normalized objective = 5.705 | 5.705 | 5.73 | 0.44% |  |
| 12 | focused power = 0.2354 | 0.2354 | 0.2268 | 3.65% |  |
| 12 | gradient norm = 0.577 | 0.577 | 0.424 | 26.52% |  |
| 12 | step = 19 | 19 | 19 | 0 |  |
| 12 | normalized objective = 5.702 | 5.702 | 5.741 | 0.68% |  |
| 12 | focused power = 0.2353 | 0.2353 | 0.2273 | 3.40% |  |
| 12 | gradient norm = 0.630 | 0.63 | 0.473 | 24.92% |  |
| 12 | step = 20 | 20 | 20 | 0 |  |
| 12 | normalized objective = 5.757 | 5.757 | 5.749 | 0.14% |  |
| 12 | focused power = 0.2376 | 0.2376 | 0.2276 | 4.21% |  |
| 12 | gradient norm = 0.521 | 0.521 | 0.482 | 7.49% |  |
| 15 | final focused power = 0.2392 | 0.2392 | 0.2289 | 4.31% |  |
| 15 | final normalized objective = 5.796 | 5.796 | 5.782 | 0.24% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 3 | 0 | 4.01% |
| 5 | 0 | 3.09% |
| 9 | 0 | 5.27% |
| 13 | 0 | 3.43% |
| 14 | 0 | 8.54% |
| 16 | 0 | 5.52% |
