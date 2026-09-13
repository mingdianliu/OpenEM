# OptimizedL3

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/OptimizedL3/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 14 on the reference side, 14 on ours, 1 that did not line up |
| Numbers compared | 1 of 6 paired; the other 5 are near-zero, see below |
| Largest relative difference | 5.74e-04% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 5 on the reference side, 5 on ours |
| Largest pixel difference | 11.36% |

Machine-readable form of everything below: [`data/OptimizedL3.json`](data/OptimizedL3.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 10 | 1.917742e+14  6.795405e+08  886593.697701  12320.95201 -0.578757  0.000014 | 1.91774e+14 | 1.91775e+14 | 5.74e-04% |  |
| 10 | 1.917742e+14  6.795405e+08  886593.697701  12320.95201 -0.578757  0.000014 | 6.7954e+08 | 6.79291e+08 | 1.30e-05% | near-zero |
| 10 | 1.917742e+14  6.795405e+08  886593.697701  12320.95201 -0.578757  0.000014 | 886594 | 886925 | 1.73e-08% | near-zero |
| 10 | 1.917742e+14  6.795405e+08  886593.697701  12320.95201 -0.578757  0.000014 | 12321 | 12321.7 | 3.75e-11% | near-zero |
| 10 | 1.917742e+14  6.795405e+08  886593.697701  12320.95201 -0.578757  0.000014 | -0.578757 | -0.581122 | 1.23e-13% | near-zero |
| 10 | 1.917742e+14  6.795405e+08  886593.697701  12320.95201 -0.578757  0.000014 | 1.4e-05 | 3.3e-05 | 9.91e-16% | near-zero |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 4 | 0 | 0 |
| 7 | 0 | 0 |
| 9 | 0 | 0.13% |
| 11 | 0 | 1.17% |
| 12 | 0 | 11.36% |
