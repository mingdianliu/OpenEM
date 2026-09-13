# AnisotropicMetamaterialBroadbandPBS

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/AnisotropicMetamaterialBroadbandPBS/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 20 on the reference side, 20 on ours, 1 that did not line up |
| Numbers compared | 2 |
| Largest relative difference | 0.07% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 8 on the reference side, 8 on ours |
| Largest pixel difference | 11.21% |

Machine-readable form of everything below: [`data/AnisotropicMetamaterialBroadbandPBS.json`](data/AnisotropicMetamaterialBroadbandPBS.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 7 | The coupling length for the TM mode is 6.95 μm. | 6.95 | 6.95 | 0 |  |
| 7 | The coupling length for the TE mode is 11296.20 μm. | 11296.2 | 11303.9 | 0.07% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 4 | 0 | 4.29% |
| 6 | 0 | 8.47% |
| 10 | 0 | 7.33% |
| 15 | 0 | 10.93% |
| 15 | 1 | 11.21% |
| 17 | 0 | 6.64% |
| 18 | 0 | 8.77% |
| 19 | 0 | 8.22% |
