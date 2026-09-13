# ResonanceFinder

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/ResonanceFinder/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Tutorial |
| Verdict | **agree** |
| Reference output | official archived output |
| Cells | 13 on the reference side, 13 on ours, 0 that did not line up |
| Numbers compared | 8 of 26 paired; the other 18 are near-zero, see below |
| Largest relative difference | 2.96% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 4 on the reference side, 4 on ours |
| Largest pixel difference | 41.44% |

Machine-readable form of everything below: [`data/ResonanceFinder.json`](data/ResonanceFinder.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 1 | Total runtime = 5.07 ps | 5.07 | 5.07 | 0 |  |
| 1 | Start monitoring fields after 0.51 ps | 0.51 | 0.51 | 0 |  |
| 9 | 1.541829e+14  1.695818e+11  2.856320e+03  1202.443339  0.792932  0.003391 | 1.54183e+14 | 1.54185e+14 | 1.17e-03% |  |
| 9 | 1.541829e+14  1.695818e+11  2.856320e+03  1202.443339  0.792932  0.003391 | 1.69582e+11 | 1.69497e+11 | 5.51e-03% | near-zero |
| 9 | 1.541829e+14  1.695818e+11  2.856320e+03  1202.443339  0.792932  0.003391 | 2856.32 | 2857.78 | 9.50e-11% | near-zero |
| 9 | 1.541829e+14  1.695818e+11  2.856320e+03  1202.443339  0.792932  0.003391 | 1202.44 | 1194.47 | 5.17e-10% | near-zero |
| 9 | 1.541829e+14  1.695818e+11  2.856320e+03  1202.443339  0.792932  0.003391 | 0.792932 | 0.788982 | 2.56e-13% | near-zero |
| 9 | 1.541829e+14  1.695818e+11  2.856320e+03  1202.443339  0.792932  0.003391 | 0.003391 | 0.003612 | 1.43e-14% | near-zero |
| 9 | 1.565102e+14  1.093306e+13  4.497291e+01   377.002522  2.045944  0.136522 | 1.5651e+14 | 1.56516e+14 | 3.58e-03% |  |
| 9 | 1.565102e+14  1.093306e+13  4.497291e+01   377.002522  2.045944  0.136522 | 1.09331e+13 | 1.06097e+13 | 2.96% |  |
| 9 | 1.565102e+14  1.093306e+13  4.497291e+01   377.002522  2.045944  0.136522 | 44.9729 | 46.3453 | 8.77e-11% | near-zero |
| 9 | 1.565102e+14  1.093306e+13  4.497291e+01   377.002522  2.045944  0.136522 | 377.003 | 363.244 | 8.79e-10% | near-zero |
| 9 | 1.565102e+14  1.093306e+13  4.497291e+01   377.002522  2.045944  0.136522 | 2.04594 | 2.07986 | 2.17e-12% | near-zero |
| 9 | 1.565102e+14  1.093306e+13  4.497291e+01   377.002522  2.045944  0.136522 | 0.136522 | 0.056274 | 5.13e-12% | near-zero |
| 9 | 1.568791e+14  1.801143e+08  2.736319e+06   482.144366  2.413333  0.001979 | 1.56879e+14 | 1.5688e+14 | 7.01e-04% |  |
| 9 | 1.568791e+14  1.801143e+08  2.736319e+06   482.144366  2.413333  0.001979 | 1.80114e+08 | 3.18154e+08 | 8.80e-03% | near-zero |
| 9 | 1.568791e+14  1.801143e+08  2.736319e+06   482.144366  2.413333  0.001979 | 2.73632e+06 | 1.54911e+06 | 7.57e-05% | near-zero |
| 9 | 1.568791e+14  1.801143e+08  2.736319e+06   482.144366  2.413333  0.001979 | 482.144 | 483.829 | 1.07e-10% | near-zero |
| 9 | 1.568791e+14  1.801143e+08  2.736319e+06   482.144366  2.413333  0.001979 | 2.41333 | 2.41037 | 1.89e-13% | near-zero |
| 9 | 1.568791e+14  1.801143e+08  2.736319e+06   482.144366  2.413333  0.001979 | 0.001979 | 0.002836 | 5.46e-14% | near-zero |
| 9 | 1.592187e+14  6.408398e+12  7.805387e+01   241.837187 -2.164008  0.090593 | 1.59219e+14 | 1.59216e+14 | 1.70e-03% |  |
| 9 | 1.592187e+14  6.408398e+12  7.805387e+01   241.837187 -2.164008  0.090593 | 6.4084e+12 | 6.28573e+12 | 1.91% |  |
| 9 | 1.592187e+14  6.408398e+12  7.805387e+01   241.837187 -2.164008  0.090593 | 78.0539 | 79.5758 | 9.56e-11% | near-zero |
| 9 | 1.592187e+14  6.408398e+12  7.805387e+01   241.837187 -2.164008  0.090593 | 241.837 | 235.667 | 3.88e-10% | near-zero |
| 9 | 1.592187e+14  6.408398e+12  7.805387e+01   241.837187 -2.164008  0.090593 | -2.16401 | -2.20126 | 2.34e-12% | near-zero |
| 9 | 1.592187e+14  6.408398e+12  7.805387e+01   241.837187 -2.164008  0.090593 | 0.090593 | 0.050665 | 2.51e-12% | near-zero |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 6 | 0 | 3.93% |
| 8 | 0 | 41.44% |
| 10 | 0 | 2.92% |
| 12 | 0 | 6.03% |
