# Autograd21GaPLightExtractor

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd21GaPLightExtractor/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **agree** |
| Reference output | official archived output |
| Cells | 24 on the reference side, 24 on ours, 0 that did not line up |
| Numbers compared | 210 |
| Largest relative difference | 41.18% |
| Over 5% / over 20% | 13 / 1 |
| Figures | 73 on the reference side, 73 on ours |
| Largest pixel difference | 11.43% |

Machine-readable form of everything below: [`data/Autograd21GaPLightExtractor.json`](data/Autograd21GaPLightExtractor.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 13 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 12 | step 1: | 1 | 1 | 0 |  |
| 12 | beta = 1.000e+00 | 1 | 1 | 0 |  |
| 12 | objective = 9.460e+00 | 9.46 | 9.463 | 0.03% |  |
| 12 | step 2: | 2 | 2 | 0 |  |
| 12 | beta = 1.075e+00 | 1.075 | 1.075 | 0 |  |
| 12 | objective = 5.400e+01 | 54 | 53.88 | 0.22% |  |
| 12 | step 3: | 3 | 3 | 0 |  |
| 12 | beta = 1.156e+00 | 1.156 | 1.156 | 0 |  |
| 12 | objective = 1.387e+02 | 138.7 | 138.1 | 0.43% |  |
| 12 | step 4: | 4 | 4 | 0 |  |
| 12 | beta = 1.243e+00 | 1.243 | 1.243 | 0 |  |
| 12 | objective = 2.504e+02 | 250.4 | 249.6 | 0.32% |  |
| 12 | step 5: | 5 | 5 | 0 |  |
| 12 | beta = 1.337e+00 | 1.337 | 1.337 | 0 |  |
| 12 | objective = 4.400e+02 | 440 | 436.9 | 0.70% |  |
| 12 | step 6: | 6 | 6 | 0 |  |
| 12 | beta = 1.438e+00 | 1.438 | 1.438 | 0 |  |
| 12 | objective = 4.974e+02 | 497.4 | 493.4 | 0.80% |  |
| 12 | step 7: | 7 | 7 | 0 |  |
| 12 | beta = 1.546e+00 | 1.546 | 1.546 | 0 |  |
| 12 | objective = 5.442e+02 | 544.2 | 540.6 | 0.66% |  |
| 12 | step 8: | 8 | 8 | 0 |  |
| 12 | beta = 1.663e+00 | 1.663 | 1.663 | 0 |  |
| 12 | objective = 5.815e+02 | 581.5 | 577.4 | 0.71% |  |
| 12 | step 9: | 9 | 9 | 0 |  |
| 12 | beta = 1.788e+00 | 1.788 | 1.788 | 0 |  |
| 12 | objective = 6.259e+02 | 625.9 | 623.9 | 0.32% |  |
| 12 | step 10: | 10 | 10 | 0 |  |
| 12 | beta = 1.922e+00 | 1.922 | 1.922 | 0 |  |
| 12 | objective = 6.747e+02 | 674.7 | 673.2 | 0.22% |  |
| 12 | step 11: | 11 | 11 | 0 |  |
| 12 | beta = 2.067e+00 | 2.067 | 2.067 | 0 |  |
| 12 | objective = 4.429e+02 | 442.9 | 420.3 | 5.10% |  |
| 12 | step 12: | 12 | 12 | 0 |  |
| 12 | beta = 2.223e+00 | 2.223 | 2.223 | 0 |  |
| 12 | objective = 6.815e+02 | 681.5 | 653.1 | 4.17% |  |
| 12 | step 13: | 13 | 13 | 0 |  |
| 12 | beta = 2.390e+00 | 2.39 | 2.39 | 0 |  |
| 12 | objective = 6.765e+02 | 676.5 | 666.7 | 1.45% |  |
| 12 | step 14: | 14 | 14 | 0 |  |
| 12 | beta = 2.570e+00 | 2.57 | 2.57 | 0 |  |
| 12 | objective = 6.670e+02 | 667 | 658.6 | 1.26% |  |
| 12 | step 15: | 15 | 15 | 0 |  |
| 12 | beta = 2.764e+00 | 2.764 | 2.764 | 0 |  |
| 12 | objective = 6.584e+02 | 658.4 | 653.6 | 0.73% |  |
| 12 | step 16: | 16 | 16 | 0 |  |
| 12 | beta = 2.972e+00 | 2.972 | 2.972 | 0 |  |
| 12 | objective = 6.687e+02 | 668.7 | 663.9 | 0.72% |  |
| 12 | step 17: | 17 | 17 | 0 |  |
| 12 | beta = 3.196e+00 | 3.196 | 3.196 | 0 |  |
| 12 | objective = 6.984e+02 | 698.4 | 690.7 | 1.10% |  |
| 12 | step 18: | 18 | 18 | 0 |  |
| 12 | beta = 3.437e+00 | 3.437 | 3.437 | 0 |  |
| 12 | objective = 6.875e+02 | 687.5 | 683.6 | 0.57% |  |
| 12 | step 19: | 19 | 19 | 0 |  |
| 12 | beta = 3.696e+00 | 3.696 | 3.696 | 0 |  |
| 12 | objective = 6.762e+02 | 676.2 | 672.8 | 0.50% |  |
| 12 | step 20: | 20 | 20 | 0 |  |
| 12 | beta = 3.974e+00 | 3.974 | 3.974 | 0 |  |
| 12 | objective = 6.772e+02 | 677.2 | 674.2 | 0.44% |  |
| 12 | step 21: | 21 | 21 | 0 |  |
| 12 | beta = 4.273e+00 | 4.273 | 4.273 | 0 |  |
| 12 | objective = 6.920e+02 | 692 | 689.1 | 0.42% |  |
| 12 | step 22: | 22 | 22 | 0 |  |
| 12 | beta = 4.595e+00 | 4.595 | 4.595 | 0 |  |
| 12 | objective = 7.149e+02 | 714.9 | 711.3 | 0.50% |  |
| 12 | step 23: | 23 | 23 | 0 |  |
| 12 | beta = 4.941e+00 | 4.941 | 4.941 | 0 |  |
| 12 | objective = 7.323e+02 | 732.3 | 724.4 | 1.08% |  |
| 12 | step 24: | 24 | 24 | 0 |  |
| 12 | beta = 5.313e+00 | 5.313 | 5.313 | 0 |  |
| 12 | objective = 7.266e+02 | 726.6 | 723.5 | 0.43% |  |
| 12 | step 25: | 25 | 25 | 0 |  |
| 12 | beta = 5.713e+00 | 5.713 | 5.713 | 0 |  |
| 12 | objective = 7.301e+02 | 730.1 | 723.5 | 0.90% |  |
| 12 | step 26: | 26 | 26 | 0 |  |
| 12 | beta = 6.144e+00 | 6.144 | 6.144 | 0 |  |
| 12 | objective = 7.448e+02 | 744.8 | 739.7 | 0.68% |  |
| 12 | step 27: | 27 | 27 | 0 |  |
| 12 | beta = 6.607e+00 | 6.607 | 6.607 | 0 |  |
| 12 | objective = 7.289e+02 | 728.9 | 724.9 | 0.55% |  |
| 12 | step 28: | 28 | 28 | 0 |  |
| 12 | beta = 7.104e+00 | 7.104 | 7.104 | 0 |  |
| 12 | objective = 7.632e+02 | 763.2 | 757.6 | 0.73% |  |
| 12 | step 29: | 29 | 29 | 0 |  |
| 12 | beta = 7.639e+00 | 7.639 | 7.639 | 0 |  |
| 12 | objective = 7.029e+02 | 702.9 | 698.4 | 0.64% |  |
| 12 | step 30: | 30 | 30 | 0 |  |
| 12 | beta = 8.215e+00 | 8.215 | 8.215 | 0 |  |
| 12 | objective = 7.251e+02 | 725.1 | 717.6 | 1.03% |  |
| 12 | step 31: | 31 | 31 | 0 |  |
| 12 | beta = 8.833e+00 | 8.833 | 8.833 | 0 |  |
| 12 | objective = 7.248e+02 | 724.8 | 713.9 | 1.50% |  |
| 12 | step 32: | 32 | 32 | 0 |  |
| 12 | beta = 9.499e+00 | 9.499 | 9.499 | 0 |  |
| 12 | objective = 7.332e+02 | 733.2 | 713.2 | 2.73% |  |
| 12 | step 33: | 33 | 33 | 0 |  |
| 12 | beta = 1.021e+01 | 10.21 | 10.21 | 0 |  |
| 12 | objective = 7.866e+02 | 786.6 | 766.1 | 2.61% |  |
| 12 | step 34: | 34 | 34 | 0 |  |
| 12 | beta = 1.098e+01 | 10.98 | 10.98 | 0 |  |
| 12 | objective = 7.724e+02 | 772.4 | 774.3 | 0.25% |  |
| 12 | step 35: | 35 | 35 | 0 |  |
| 12 | beta = 1.181e+01 | 11.81 | 11.81 | 0 |  |
| 12 | objective = 6.893e+02 | 689.3 | 676.2 | 1.90% |  |
| 12 | step 36: | 36 | 36 | 0 |  |
| 12 | beta = 1.270e+01 | 12.7 | 12.7 | 0 |  |
| 12 | objective = 6.774e+02 | 677.4 | 689.2 | 1.74% |  |
| 12 | step 37: | 37 | 37 | 0 |  |
| 12 | beta = 1.366e+01 | 13.66 | 13.66 | 0 |  |
| 12 | objective = 7.939e+02 | 793.9 | 740.9 | 6.68% |  |
| 12 | step 38: | 38 | 38 | 0 |  |
| 12 | beta = 1.469e+01 | 14.69 | 14.69 | 0 |  |
| 12 | objective = 5.651e+02 | 565.1 | 797.8 | 41.18% |  |
| 12 | step 39: | 39 | 39 | 0 |  |
| 12 | beta = 1.579e+01 | 15.79 | 15.79 | 0 |  |
| 12 | objective = 6.405e+02 | 640.5 | 723.3 | 12.93% |  |
| 12 | step 40: | 40 | 40 | 0 |  |
| 12 | beta = 1.698e+01 | 16.98 | 16.98 | 0 |  |
| 12 | objective = 6.880e+02 | 688 | 735.6 | 6.92% |  |
| 12 | step 41: | 41 | 41 | 0 |  |
| 12 | beta = 1.826e+01 | 18.26 | 18.26 | 0 |  |
| 12 | objective = 6.188e+02 | 618.8 | 732.4 | 18.36% |  |
| 12 | step 42: | 42 | 42 | 0 |  |
| 12 | beta = 1.964e+01 | 19.64 | 19.64 | 0 |  |
| 12 | objective = 7.190e+02 | 719 | 767.3 | 6.72% |  |
| 12 | step 43: | 43 | 43 | 0 |  |
| 12 | beta = 2.111e+01 | 21.11 | 21.11 | 0 |  |
| 12 | objective = 7.708e+02 | 770.8 | 760.4 | 1.35% |  |
| 12 | step 44: | 44 | 44 | 0 |  |
| 12 | beta = 2.270e+01 | 22.7 | 22.7 | 0 |  |
| 12 | objective = 7.397e+02 | 739.7 | 736.4 | 0.45% |  |
| 12 | step 45: | 45 | 45 | 0 |  |
| 12 | beta = 2.441e+01 | 24.41 | 24.41 | 0 |  |
| 12 | objective = 7.303e+02 | 730.3 | 749.5 | 2.63% |  |
| 12 | step 46: | 46 | 46 | 0 |  |
| 12 | beta = 2.625e+01 | 26.25 | 26.25 | 0 |  |
| 12 | objective = 7.382e+02 | 738.2 | 773.1 | 4.73% |  |
| 12 | step 47: | 47 | 47 | 0 |  |
| 12 | beta = 2.823e+01 | 28.23 | 28.23 | 0 |  |
| 12 | objective = 8.222e+02 | 822.2 | 817.7 | 0.55% |  |
| 12 | step 48: | 48 | 48 | 0 |  |
| 12 | beta = 3.036e+01 | 30.36 | 30.36 | 0 |  |
| 12 | objective = 8.351e+02 | 835.1 | 775.3 | 7.16% |  |
| 12 | step 49: | 49 | 49 | 0 |  |
| 12 | beta = 3.264e+01 | 32.64 | 32.64 | 0 |  |
| 12 | objective = 7.341e+02 | 734.1 | 716.2 | 2.44% |  |
| 12 | step 50: | 50 | 50 | 0 |  |
| 12 | beta = 3.510e+01 | 35.1 | 35.1 | 0 |  |
| 12 | objective = 8.397e+02 | 839.7 | 761.3 | 9.34% |  |
| 12 | step 51: | 51 | 51 | 0 |  |
| 12 | beta = 3.775e+01 | 37.75 | 37.75 | 0 |  |
| 12 | objective = 8.323e+02 | 832.3 | 825.3 | 0.84% |  |
| 12 | step 52: | 52 | 52 | 0 |  |
| 12 | beta = 4.059e+01 | 40.59 | 40.59 | 0 |  |
| 12 | objective = 8.228e+02 | 822.8 | 793.7 | 3.54% |  |
| 12 | step 53: | 53 | 53 | 0 |  |
| 12 | beta = 4.365e+01 | 43.65 | 43.65 | 0 |  |
| 12 | objective = 8.273e+02 | 827.3 | 769 | 7.05% |  |
| 12 | step 54: | 54 | 54 | 0 |  |
| 12 | beta = 4.693e+01 | 46.93 | 46.93 | 0 |  |
| 12 | objective = 8.312e+02 | 831.2 | 769.3 | 7.45% |  |
| 12 | step 55: | 55 | 55 | 0 |  |
| 12 | beta = 5.047e+01 | 50.47 | 50.47 | 0 |  |
| 12 | objective = 8.298e+02 | 829.8 | 784.2 | 5.50% |  |
| 12 | step 56: | 56 | 56 | 0 |  |
| 12 | beta = 5.427e+01 | 54.27 | 54.27 | 0 |  |
| 12 | objective = 8.385e+02 | 838.5 | 812.3 | 3.12% |  |
| 12 | step 57: | 57 | 57 | 0 |  |
| 12 | beta = 5.836e+01 | 58.36 | 58.36 | 0 |  |
| 12 | objective = 8.426e+02 | 842.6 | 812.7 | 3.55% |  |
| 12 | step 58: | 58 | 58 | 0 |  |
| 12 | beta = 6.275e+01 | 62.75 | 62.75 | 0 |  |
| 12 | objective = 8.444e+02 | 844.4 | 830.8 | 1.61% |  |
| 12 | step 59: | 59 | 59 | 0 |  |
| 12 | beta = 6.748e+01 | 67.48 | 67.48 | 0 |  |
| 12 | objective = 8.481e+02 | 848.1 | 841 | 0.84% |  |
| 12 | step 60: | 60 | 60 | 0 |  |
| 12 | beta = 7.256e+01 | 72.56 | 72.56 | 0 |  |
| 12 | objective = 8.482e+02 | 848.2 | 817.2 | 3.65% |  |
| 12 | step 61: | 61 | 61 | 0 |  |
| 12 | beta = 7.803e+01 | 78.03 | 78.03 | 0 |  |
| 12 | objective = 8.469e+02 | 846.9 | 835.6 | 1.33% |  |
| 12 | step 62: | 62 | 62 | 0 |  |
| 12 | beta = 8.391e+01 | 83.91 | 83.91 | 0 |  |
| 12 | objective = 8.595e+02 | 859.5 | 844.1 | 1.79% |  |
| 12 | step 63: | 63 | 63 | 0 |  |
| 12 | beta = 9.023e+01 | 90.23 | 90.23 | 0 |  |
| 12 | objective = 8.660e+02 | 866 | 849.1 | 1.95% |  |
| 12 | step 64: | 64 | 64 | 0 |  |
| 12 | beta = 9.702e+01 | 97.02 | 97.02 | 0 |  |
| 12 | objective = 8.605e+02 | 860.5 | 825.7 | 4.04% |  |
| 12 | step 65: | 65 | 65 | 0 |  |
| 12 | beta = 1.043e+02 | 104.3 | 104.3 | 0 |  |
| 12 | objective = 8.615e+02 | 861.5 | 844 | 2.03% |  |
| 12 | step 66: | 66 | 66 | 0 |  |
| 12 | beta = 1.122e+02 | 112.2 | 112.2 | 0 |  |
| 12 | objective = 8.671e+02 | 867.1 | 838.6 | 3.29% |  |
| 12 | step 67: | 67 | 67 | 0 |  |
| 12 | beta = 1.206e+02 | 120.6 | 120.6 | 0 |  |
| 12 | objective = 8.722e+02 | 872.2 | 842.7 | 3.38% |  |
| 12 | step 68: | 68 | 68 | 0 |  |
| 12 | beta = 1.297e+02 | 129.7 | 129.7 | 0 |  |
| 12 | objective = 8.725e+02 | 872.5 | 843 | 3.38% |  |
| 12 | step 69: | 69 | 69 | 0 |  |
| 12 | beta = 1.395e+02 | 139.5 | 139.5 | 0 |  |
| 12 | objective = 8.722e+02 | 872.2 | 847.8 | 2.80% |  |
| 12 | step 70: | 70 | 70 | 0 |  |
| 12 | beta = 1.500e+02 | 150 | 150 | 0 |  |
| 12 | objective = 8.727e+02 | 872.7 | 826.9 | 5.25% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 7 | 0 | 3.38% |
| 12 | 0 | 0 |
| 12 | 1 | 0 |
| 12 | 2 | 0 |
| 12 | 3 | 0 |
| 12 | 4 | 0 |
| 12 | 5 | 0 |
| 12 | 6 | 0 |
| 12 | 7 | 0 |
| 12 | 8 | 0 |
| 12 | 9 | 0 |
| 12 | 10 | 0 |
| 12 | 11 | 0 |
| 12 | 12 | 0 |
| 12 | 13 | 0 |
| 12 | 14 | 0 |
| 12 | 15 | 0.06% |
| 12 | 16 | 0.11% |
| 12 | 17 | 0.17% |
| 12 | 18 | 0.27% |
| 12 | 19 | 0.39% |
| 12 | 20 | 0.57% |
| 12 | 21 | 0.68% |
| 12 | 22 | 0.83% |
| 12 | 23 | 0.98% |
| 12 | 24 | 1.37% |
| 12 | 25 | 1.91% |
| 12 | 26 | 2.33% |
| 12 | 27 | 2.79% |
| 12 | 28 | 3.46% |
| 12 | 29 | 4.11% |
| 12 | 30 | 4.86% |
| 12 | 31 | 5.41% |
| 12 | 32 | 5.78% |
| 12 | 33 | 6.17% |
| 12 | 34 | 6.49% |
| 12 | 35 | 6.68% |
| 12 | 36 | 6.90% |
| 12 | 37 | 8.91% |
| 12 | 38 | 10.04% |
| 12 | 39 | 9.77% |
| 12 | 40 | 10.08% |
| 12 | 41 | 9.95% |
| 12 | 42 | 9.81% |
| 12 | 43 | 9.75% |
| 12 | 44 | 9.42% |
| 12 | 45 | 9.22% |
| 12 | 46 | 9.56% |
| 12 | 47 | 10.14% |
| 12 | 48 | 10.21% |
| 12 | 49 | 9.49% |
| 12 | 50 | 9.60% |
| 12 | 51 | 9.80% |
| 12 | 52 | 9.66% |
| 12 | 53 | 9.93% |
| 12 | 54 | 9.88% |
| 12 | 55 | 10.05% |
| 12 | 56 | 10.31% |
| 12 | 57 | 10.41% |
| 12 | 58 | 10.39% |
| 12 | 59 | 10.42% |
| 12 | 60 | 10.38% |
| 12 | 61 | 10.53% |
| 12 | 62 | 10.70% |
| 12 | 63 | 10.72% |
| 12 | 64 | 10.78% |
| 12 | 65 | 11.01% |
| 12 | 66 | 11.20% |
| 12 | 67 | 11.43% |
| 12 | 68 | 11.42% |
| 12 | 69 | 11.26% |
| 13 | 0 | 4.85% |
| 20 | 0 | 10.27% |
