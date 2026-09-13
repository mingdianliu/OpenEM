# NanobeamCavity

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/NanobeamCavity/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **close** |
| Reference output | reference rerun |
| Cells | 26 on the reference side, 26 on ours, 0 that did not line up |
| Numbers compared | 17 of 29 paired; the other 12 are near-zero, see below |
| Largest relative difference | 730% |
| Over 5% / over 20% | 4 / 4 |
| Figures | 7 on the reference side, 7 on ours |
| Largest pixel difference | 10.76% |

Machine-readable form of everything below: [`data/NanobeamCavity.json`](data/NanobeamCavity.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 4 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 15 | 1.941258e+14  7.709018e+07  7.911051e+06  63674.827463  0.989458  0.00361 | 1.94126e+14 | 1.94124e+14 | 7.21e-04% |  |
| 15 | 1.941258e+14  7.709018e+07  7.911051e+06  63674.827463  0.989458  0.00361 | 7.70902e+07 | 8.67959e+07 | 5.00e-04% | near-zero |
| 15 | 1.941258e+14  7.709018e+07  7.911051e+06  63674.827463  0.989458  0.00361 | 7.91105e+06 | 7.02637e+06 | 4.56e-05% | near-zero |
| 15 | 1.941258e+14  7.709018e+07  7.911051e+06  63674.827463  0.989458  0.00361 | 63674.8 | 63469.2 | 1.06e-08% | near-zero |
| 15 | 1.941258e+14  7.709018e+07  7.911051e+06  63674.827463  0.989458  0.00361 | 0.989458 | 0.990619 | 5.98e-14% | near-zero |
| 15 | 1.941258e+14  7.709018e+07  7.911051e+06  63674.827463  0.989458  0.00361 | 0.00361 | 0.003783 | 8.91e-15% | near-zero |
| 15 | 1.941258e+14  1544.32022 | 1.94126e+14 | 1.94124e+14 | 7.21e-04% |  |
| 15 | 1.941258e+14  1544.32022 | 1544.32 | 1544.33 | 6.05e-13% | near-zero |
| 17 | Mode volume: 0.38 (lambda/n)^3 | 0.38 | 0.38 | 0 |  |
| 17 | Mode volume: 0.38 (lambda/n)^3 | 3 | 3 | 0 |  |
| 19 | Mode volume: 0.37 (lambda/n)^3 | 0.37 | 3.07 | 730% |  |
| 19 | Mode volume: 0.37 (lambda/n)^3 | 3 | 3 | 0 |  |
| 23 | Mode volume: 0.35 (lambda/n)^3 | 0.35 | 0.35 | 0 |  |
| 23 | Mode volume: 0.35 (lambda/n)^3 | 3 | 3 | 0 |  |
| 23 | 1.941612e+14  1.635266e+08  3.730129e+06  18817.985879 -2.816392  0.000098 | 1.94161e+14 | 1.94155e+14 | 3.09e-03% |  |
| 23 | 1.941612e+14  1.635266e+08  3.730129e+06  18817.985879 -2.816392  0.000098 | 1.63527e+08 | 1.44791e+08 | 9.65e-04% | near-zero |
| 23 | 1.941612e+14  1.635266e+08  3.730129e+06  18817.985879 -2.816392  0.000098 | 3.73013e+06 | 4.21266e+06 | 2.49e-05% | near-zero |
| 23 | 1.941612e+14  1.635266e+08  3.730129e+06  18817.985879 -2.816392  0.000098 | 18818 | 18736.7 | 4.19e-09% | near-zero |
| 23 | 1.941612e+14  1.635266e+08  3.730129e+06  18817.985879 -2.816392  0.000098 | -2.81639 | -2.81159 | 2.47e-13% | near-zero |
| 23 | 1.941612e+14  1.635266e+08  3.730129e+06  18817.985879 -2.816392  0.000098 | 9.8e-05 | 0.000251 | 7.88e-15% | near-zero |
| 23 | 1.941612e+14  1544.038947 | 1.94161e+14 | 1.94155e+14 | 3.09e-03% |  |
| 23 | 1.941612e+14  1544.038947 | 1544.04 | 1544.09 | 2.47e-12% | near-zero |
| 25 | Q x = 2.69 MM | 2.69 | 1.7 | 36.80% |  |
| 25 | Q -x = 81.60 MM | 81.6 | 59.4 | 27.21% |  |
| 25 | Q y = 16.45 MM | 16.45 | 16.46 | 0.06% |  |
| 25 | Q -y = 16.45 MM | 16.45 | 16.46 | 0.06% |  |
| 25 | Q z = 18.06 MM | 18.06 | 17.41 | 3.60% |  |
| 25 | Q -z = 18.06 MM | 18.06 | 17.41 | 3.60% |  |
| 25 | Q total Q = 1.62 MM | 1.62 | 1.19 | 26.54% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 6 | 0 | 7.78% |
| 6 | 1 | 10.76% |
| 6 | 2 | 8.35% |
| 10 | 0 | 5.72% |
| 13 | 0 | 3.25% |
| 20 | 0 | 5.11% |
| 22 | 0 | 3.26% |
