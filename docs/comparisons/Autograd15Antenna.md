# Autograd15Antenna

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd15Antenna/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **trajectory divergence** |
| Reference output | official archived output |
| Cells | 31 on the reference side, 31 on ours, 0 that did not line up |
| Numbers compared | 128 |
| Largest relative difference | 574% |
| Over 5% / over 20% | 36 / 24 |
| Figures | 28 on the reference side, 28 on ours |
| Largest pixel difference | 8.59% |

Machine-readable form of everything below: [`data/Autograd15Antenna.json`](data/Autograd15Antenna.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 36 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 21 | Intensity without structure = 2083.5806 (au) | 2083.58 | 2083.96 | 0.02% |  |
| 25 | enhancement = 1.20e+00 | 1.2 | 1.18 | 1.67% |  |
| 25 | penalty val = 9.98e-01 | 0.998 | 0.998 | 0 |  |
| 25 | objective = 2.17e-03 | 0.00217 | 0.00214 | 1.38% |  |
| 25 | starting objective function value = 0.0021721132071270124 | 0.00217211 | 0.00214152 | 1.41% |  |
| 26 | step = 1 | 1 | 1 | 0 |  |
| 26 | enhancement = 1.20e+00 | 1.2 | 1.18 | 1.67% |  |
| 26 | penalty val = 9.98e-01 | 0.998 | 0.998 | 0 |  |
| 26 | objective = 2.17e-03 | 0.00217 | 0.00214 | 1.38% |  |
| 26 | beta = 1.0000e+00 | 1 | 1 | 0 |  |
| 26 | grad_norm = 1.0915e-02 | 0.010915 | 0.010867 | 0.44% |  |
| 26 | step = 2 | 2 | 2 | 0 |  |
| 26 | enhancement = 2.33e+00 | 2.33 | 2.28 | 2.15% |  |
| 26 | penalty val = 9.98e-01 | 0.998 | 0.998 | 0 |  |
| 26 | objective = 4.51e-03 | 0.00451 | 0.00441 | 2.22% |  |
| 26 | beta = 3.7000e+00 | 3.7 | 3.7 | 0 |  |
| 26 | grad_norm = 2.4471e-02 | 0.024471 | 0.024179 | 1.19% |  |
| 26 | step = 3 | 3 | 3 | 0 |  |
| 26 | enhancement = 4.20e+00 | 4.2 | 4.1 | 2.38% |  |
| 26 | penalty val = 9.97e-01 | 0.997 | 0.997 | 0 |  |
| 26 | objective = 1.36e-02 | 0.0136 | 0.0133 | 2.21% |  |
| 26 | beta = 6.4000e+00 | 6.4 | 6.4 | 0 |  |
| 26 | grad_norm = 7.9013e-02 | 0.079013 | 0.078307 | 0.89% |  |
| 26 | step = 4 | 4 | 4 | 0 |  |
| 26 | enhancement = 7.29e+00 | 7.29 | 7.14 | 2.06% |  |
| 26 | penalty val = 9.90e-01 | 0.99 | 0.99 | 0 |  |
| 26 | objective = 7.31e-02 | 0.0731 | 0.0717 | 1.92% |  |
| 26 | beta = 9.1000e+00 | 9.1 | 9.1 | 0 |  |
| 26 | grad_norm = 4.4960e-01 | 0.4496 | 0.44662 | 0.66% |  |
| 26 | step = 5 | 5 | 5 | 0 |  |
| 26 | enhancement = 1.21e+01 | 12.1 | 11.8 | 2.48% |  |
| 26 | penalty val = 9.61e-01 | 0.961 | 0.961 | 0 |  |
| 26 | objective = 4.75e-01 | 0.475 | 0.464 | 2.32% |  |
| 26 | beta = 1.1800e+01 | 11.8 | 11.8 | 0 |  |
| 26 | grad_norm = 7.0505e+00 | 7.0505 | 7.0779 | 0.39% |  |
| 26 | step = 6 | 6 | 6 | 0 |  |
| 26 | enhancement = 2.01e+01 | 20.1 | 21.2 | 5.47% |  |
| 26 | penalty val = 8.40e-01 | 0.84 | 0.84 | 0 |  |
| 26 | objective = 3.22e+00 | 3.22 | 3.4 | 5.59% |  |
| 26 | beta = 1.4500e+01 | 14.5 | 14.5 | 0 |  |
| 26 | grad_norm = 2.8487e+01 | 28.487 | 29.13 | 2.26% |  |
| 26 | step = 7 | 7 | 7 | 0 |  |
| 26 | enhancement = 3.48e+01 | 34.8 | 35.9 | 3.16% |  |
| 26 | penalty val = 5.08e-01 | 0.508 | 0.509 | 0.20% |  |
| 26 | objective = 1.71e+01 | 17.1 | 17.6 | 2.92% |  |
| 26 | beta = 1.7200e+01 | 17.2 | 17.2 | 0 |  |
| 26 | grad_norm = 1.9941e+02 | 199.41 | 204.61 | 2.61% |  |
| 26 | step = 8 | 8 | 8 | 0 |  |
| 26 | enhancement = 6.64e+01 | 66.4 | 65.6 | 1.20% |  |
| 26 | penalty val = 2.80e-01 | 0.28 | 0.28 | 0 |  |
| 26 | objective = 4.79e+01 | 47.9 | 47.3 | 1.25% |  |
| 26 | beta = 1.9900e+01 | 19.9 | 19.9 | 0 |  |
| 26 | grad_norm = 4.8410e+02 | 484.1 | 488.84 | 0.98% |  |
| 26 | step = 9 | 9 | 9 | 0 |  |
| 26 | enhancement = 1.16e+02 | 116 | 118 | 1.72% |  |
| 26 | penalty val = 1.81e-01 | 0.181 | 0.181 | 0 |  |
| 26 | objective = 9.51e+01 | 95.1 | 96.5 | 1.47% |  |
| 26 | beta = 2.2600e+01 | 22.6 | 22.6 | 0 |  |
| 26 | grad_norm = 9.4315e+02 | 943.15 | 970.21 | 2.87% |  |
| 26 | step = 10 | 10 | 10 | 0 |  |
| 26 | enhancement = 2.89e+02 | 289 | 277 | 4.15% |  |
| 26 | penalty val = 1.09e-01 | 0.109 | 0.109 | 0 |  |
| 26 | objective = 2.58e+02 | 258 | 246 | 4.65% |  |
| 26 | beta = 2.5300e+01 | 25.3 | 25.3 | 0 |  |
| 26 | grad_norm = 1.6785e+03 | 1678.5 | 1660.6 | 1.07% |  |
| 26 | step = 11 | 11 | 11 | 0 |  |
| 26 | enhancement = 1.69e+02 | 169 | 178 | 5.33% |  |
| 26 | penalty val = 6.00e-02 | 0.06 | 0.0599 | 0.17% |  |
| 26 | objective = 1.58e+02 | 158 | 167 | 5.70% |  |
| 26 | beta = 2.8000e+01 | 28 | 28 | 0 |  |
| 26 | grad_norm = 2.4872e+03 | 2487.2 | 1856.2 | 25.37% |  |
| 26 | step = 12 | 12 | 12 | 0 |  |
| 26 | enhancement = 2.41e+02 | 241 | 239 | 0.83% |  |
| 26 | penalty val = 3.40e-02 | 0.034 | 0.034 | 0 |  |
| 26 | objective = 2.33e+02 | 233 | 231 | 0.86% |  |
| 26 | beta = 3.0700e+01 | 30.7 | 30.7 | 0 |  |
| 26 | grad_norm = 7.4299e+03 | 7429.9 | 2870.3 | 61.37% |  |
| 26 | step = 13 | 13 | 13 | 0 |  |
| 26 | enhancement = 4.95e+02 | 495 | 392 | 20.81% |  |
| 26 | penalty val = 2.07e-02 | 0.0207 | 0.02 | 3.38% |  |
| 26 | objective = 4.85e+02 | 485 | 384 | 20.82% |  |
| 26 | beta = 3.3400e+01 | 33.4 | 33.4 | 0 |  |
| 26 | grad_norm = 2.4180e+03 | 2418 | 8537 | 253% |  |
| 26 | step = 14 | 14 | 14 | 0 |  |
| 26 | enhancement = 5.05e+02 | 505 | 554 | 9.70% |  |
| 26 | penalty val = 1.48e-02 | 0.0148 | 0.0135 | 8.78% |  |
| 26 | objective = 4.97e+02 | 497 | 546 | 9.86% |  |
| 26 | beta = 3.6100e+01 | 36.1 | 36.1 | 0 |  |
| 26 | grad_norm = 2.9640e+03 | 2964 | 4020.6 | 35.65% |  |
| 26 | step = 15 | 15 | 15 | 0 |  |
| 26 | enhancement = 5.88e+02 | 588 | 610 | 3.74% |  |
| 26 | penalty val = 1.29e-02 | 0.0129 | 0.0115 | 10.85% |  |
| 26 | objective = 5.80e+02 | 580 | 603 | 3.97% |  |
| 26 | beta = 3.8800e+01 | 38.8 | 38.8 | 0 |  |
| 26 | grad_norm = 1.3945e+03 | 1394.5 | 3237.4 | 132% |  |
| 26 | step = 16 | 16 | 16 | 0 |  |
| 26 | enhancement = 6.60e+02 | 660 | 657 | 0.45% |  |
| 26 | penalty val = 1.35e-02 | 0.0135 | 0.012 | 11.11% |  |
| 26 | objective = 6.51e+02 | 651 | 649 | 0.31% |  |
| 26 | beta = 4.1500e+01 | 41.5 | 41.5 | 0 |  |
| 26 | grad_norm = 3.1391e+03 | 3139.1 | 4911.5 | 56.46% |  |
| 26 | step = 17 | 17 | 17 | 0 |  |
| 26 | enhancement = 6.67e+02 | 667 | 747 | 11.99% |  |
| 26 | penalty val = 1.64e-02 | 0.0164 | 0.0132 | 19.51% |  |
| 26 | objective = 6.56e+02 | 656 | 737 | 12.35% |  |
| 26 | beta = 4.4200e+01 | 44.2 | 44.2 | 0 |  |
| 26 | grad_norm = 2.5701e+03 | 2570.1 | 12680 | 393% |  |
| 26 | step = 18 | 18 | 18 | 0 |  |
| 26 | enhancement = 6.17e+02 | 617 | 950 | 53.97% |  |
| 26 | penalty val = 2.01e-02 | 0.0201 | 0.0131 | 34.83% |  |
| 26 | objective = 6.04e+02 | 604 | 937 | 55.13% |  |
| 26 | beta = 4.6900e+01 | 46.9 | 46.9 | 0 |  |
| 26 | grad_norm = 3.9547e+03 | 3954.7 | 26671 | 574% |  |
| 26 | step = 19 | 19 | 19 | 0 |  |
| 26 | enhancement = 6.23e+02 | 623 | 996 | 59.87% |  |
| 26 | penalty val = 2.20e-02 | 0.022 | 0.0131 | 40.45% |  |
| 26 | objective = 6.09e+02 | 609 | 983 | 61.41% |  |
| 26 | beta = 4.9600e+01 | 49.6 | 49.6 | 0 |  |
| 26 | grad_norm = 3.3044e+03 | 3304.4 | 18695 | 466% |  |
| 26 | step = 20 | 20 | 20 | 0 |  |
| 26 | enhancement = 6.63e+02 | 663 | 970 | 46.30% |  |
| 26 | penalty val = 5.72e-02 | 0.0572 | 0.0225 | 60.66% |  |
| 26 | objective = 6.25e+02 | 625 | 948 | 51.68% |  |
| 26 | beta = 5.2300e+01 | 52.3 | 52.3 | 0 |  |
| 26 | grad_norm = 4.3060e+03 | 4306 | 14586 | 239% |  |
| 27 | enhancement = 7.01e+02 | 701 | 1040 | 48.36% |  |
| 27 | penalty val = 1.31e-01 | 0.131 | 0.0402 | 69.31% |  |
| 27 | objective = 6.09e+02 | 609 | 993 | 63.05% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 5 | 0 | 6.01% |
| 6 | 0 | 5.80% |
| 11 | 0 | 5.44% |
| 19 | 0 | 4.86% |
| 20 | 0 | 4.74% |
| 25 | 0 | 4.68% |
| 26 | 0 | 0 |
| 26 | 1 | 0 |
| 26 | 2 | 0 |
| 26 | 3 | 0 |
| 26 | 4 | 0 |
| 26 | 5 | 0 |
| 26 | 6 | 0 |
| 26 | 7 | 0 |
| 26 | 8 | 0 |
| 26 | 9 | 0 |
| 26 | 10 | 0 |
| 26 | 11 | 0 |
| 26 | 12 | 0 |
| 26 | 13 | 0.07% |
| 26 | 14 | 0.12% |
| 26 | 15 | 0.13% |
| 26 | 16 | 0.16% |
| 26 | 17 | 0.51% |
| 26 | 18 | 1.14% |
| 26 | 19 | 1.58% |
| 28 | 0 | 4.68% |
| 30 | 0 | 8.59% |
