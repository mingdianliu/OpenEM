# MaxwellStressTensor

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/MaxwellStressTensor/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 15 on the reference side, 15 on ours, 0 that did not line up |
| Numbers compared | 5 of 6 paired; the other 1 are near-zero, see below |
| Largest relative difference | 0 |
| Over 5% / over 20% | 0 / 0 |
| Figures | 4 on the reference side, 4 on ours |
| Largest pixel difference | 4.94% |

Machine-readable form of everything below: [`data/MaxwellStressTensor.json`](data/MaxwellStressTensor.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 6 | 0 | 0 | 0 | 0 | near-zero |
| 6 | 10 | 10 | 10 | 0 |  |
| 6 | 20 | 20 | 20 | 0 |  |
| 6 | 30 | 30 | 30 | 0 |  |
| 6 | 40 | 40 | 40 | 0 |  |
| 6 | 50 | 50 | 50 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 4 | 0 | 4.42% |
| 8 | 0 | 4.94% |
| 11 | 0 | 4.65% |
| 14 | 0 | 4.04% |
