# VizData

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/VizData/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Tutorial |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 36 on the reference side, 36 on ours, 6 that did not line up |
| Numbers compared | 11 of 12 paired; the other 1 are near-zero, see below |
| Largest relative difference | 0.21% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 8 on the reference side, 8 on ours |
| Largest pixel difference | 1.11% |

Machine-readable form of everything below: [`data/VizData.json`](data/VizData.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 7 | (5.0, 5.0, 5.0) | 5 | 5 | 0 |  |
| 7 | (5.0, 5.0, 5.0) | 5 | 5 | 0 |  |
| 7 | (5.0, 5.0, 5.0) | 5 | 5 | 0 |  |
| 10 | shape of flux dataset = (11,) | 11 | 11 | 0 |  |
| 11 | at central frequency, flux is 9.77e-01 W | 0.977 | 0.978 | 0.10% |  |
| 12 | at 4th frequency, flux is 9.72e-01 W | 0.972 | 0.974 | 0.21% |  |
| 13 | at an intermediate (not stored) frequency, flux is 9.78e-01 W | 0.978 | 0.98 | 0.20% |  |
| 17 | (2, 11, 3) | 2 | 2 | 0 |  |
| 17 | (2, 11, 3) | 11 | 11 | 0 |  |
| 17 | (2, 11, 3) | 3 | 3 | 0 |  |
| 32 | [108324 rows x 5 columns] | 108324 | 108324 | 0 |  |
| 32 | [108324 rows x 5 columns] | 5 | 5 | 0 | near-zero |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 3 | 0 | 0 |
| 16 | 0 | 0.35% |
| 18 | 0 | 1.11% |
| 26 | 0 | 0.50% |
| 27 | 0 | 0.37% |
| 28 | 0 | 1.08% |
| 29 | 0 | 0.26% |
| 34 | 0 | 0.61% |
