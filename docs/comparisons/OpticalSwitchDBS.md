# OpticalSwitchDBS

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/OpticalSwitchDBS/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 18 on the reference side, 18 on ours, 0 that did not line up |
| Numbers compared | 25 |
| Largest relative difference | 0.12% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 5 on the reference side, 5 on ours |
| Largest pixel difference | 17.74% |

Machine-readable form of everything below: [`data/OpticalSwitchDBS.json`](data/OpticalSwitchDBS.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 11 | Iteration 1 starts. | 1 | 1 | 0 |  |
| 11 | A best objective of 0.832 found after 1 evaluation(s). | 0.832 | 0.833 | 0.12% |  |
| 11 | A best objective of 0.832 found after 1 evaluation(s). | 1 | 1 | 0 |  |
| 11 | A best objective of 0.859 found after 2 evaluation(s). | 0.859 | 0.859 | 0 |  |
| 11 | A best objective of 0.859 found after 2 evaluation(s). | 2 | 2 | 0 |  |
| 11 | A best objective of 0.863 found after 6 evaluation(s). | 0.863 | 0.863 | 0 |  |
| 11 | A best objective of 0.863 found after 6 evaluation(s). | 6 | 6 | 0 |  |
| 11 | A best objective of 0.870 found after 7 evaluation(s). | 0.87 | 0.87 | 0 |  |
| 11 | A best objective of 0.870 found after 7 evaluation(s). | 7 | 7 | 0 |  |
| 11 | A best objective of 0.871 found after 8 evaluation(s). | 0.871 | 0.871 | 0 |  |
| 11 | A best objective of 0.871 found after 8 evaluation(s). | 8 | 8 | 0 |  |
| 11 | A best objective of 0.882 found after 9 evaluation(s). | 0.882 | 0.882 | 0 |  |
| 11 | A best objective of 0.882 found after 9 evaluation(s). | 9 | 9 | 0 |  |
| 11 | A best objective of 0.899 found after 10 evaluation(s). | 0.899 | 0.899 | 0 |  |
| 11 | A best objective of 0.899 found after 10 evaluation(s). | 10 | 10 | 0 |  |
| 11 | A best objective of 0.916 found after 11 evaluation(s). | 0.916 | 0.915 | 0.11% |  |
| 11 | A best objective of 0.916 found after 11 evaluation(s). | 11 | 11 | 0 |  |
| 11 | A best objective of 0.918 found after 13 evaluation(s). | 0.918 | 0.917 | 0.11% |  |
| 11 | A best objective of 0.918 found after 13 evaluation(s). | 13 | 13 | 0 |  |
| 11 | A best objective of 0.922 found after 18 evaluation(s). | 0.922 | 0.921 | 0.11% |  |
| 11 | A best objective of 0.922 found after 18 evaluation(s). | 18 | 18 | 0 |  |
| 11 | A best objective of 0.924 found after 24 evaluation(s). | 0.924 | 0.924 | 0 |  |
| 11 | A best objective of 0.924 found after 24 evaluation(s). | 24 | 24 | 0 |  |
| 11 | A best objective of 0.927 found after 26 evaluation(s). | 0.927 | 0.927 | 0 |  |
| 11 | A best objective of 0.927 found after 26 evaluation(s). | 26 | 26 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 5 | 0 | 17.74% |
| 12 | 0 | 3.82% |
| 14 | 0 | 4.04% |
| 15 | 0 | 4.01% |
| 16 | 0 | 5.99% |
