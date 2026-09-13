# MMIPowerSplitter2x2

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/MMIPowerSplitter2x2/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 13 on the reference side, 13 on ours, 1 that did not line up |
| Numbers compared | 5 |
| Largest relative difference | 0 |
| Over 5% / over 20% | 0 / 0 |
| Figures | 4 on the reference side, 4 on ours |
| Largest pixel difference | 0.34% |

Machine-readable form of everything below: [`data/MMIPowerSplitter2x2.json`](data/MMIPowerSplitter2x2.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 11 | Average splitting ratio 1: 0.500 | 1 | 1 | 0 |  |
| 11 | Average splitting ratio 1: 0.500 | 0.5 | 0.5 | 0 |  |
| 11 | Average splitting ratio 2: 0.500 | 2 | 2 | 0 |  |
| 11 | Average splitting ratio 2: 0.500 | 0.5 | 0.5 | 0 |  |
| 11 | Average total power: 0.869 | 0.869 | 0.869 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 4 | 0 | 0 |
| 7 | 0 | 0 |
| 11 | 0 | 0.04% |
| 12 | 0 | 0.34% |
