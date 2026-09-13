# BroadbandDirectionalCoupler

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/BroadbandDirectionalCoupler/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 25 on the reference side, 25 on ours, 0 that did not line up |
| Numbers compared | 4 |
| Largest relative difference | 1.92e-09% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 13 on the reference side, 13 on ours |
| Largest pixel difference | 11.67% |

Machine-readable form of everything below: [`data/BroadbandDirectionalCoupler.json`](data/BroadbandDirectionalCoupler.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 5 | Effective indices: 2.459162628532339, 2.43821523251537 | 2.45916 | 2.45916 | 7.01e-11% |  |
| 5 | Effective indices: 2.459162628532339, 2.43821523251537 | 2.43822 | 2.43822 | 2.53e-11% |  |
| 7 | Effective indices: 2.5682454486522666, 2.2176841960818794 | 2.56825 | 2.56825 | 9.10e-12% |  |
| 7 | Effective indices: 2.5682454486522666, 2.2176841960818794 | 2.21768 | 2.21768 | 1.92e-09% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 4 | 0 | 3.04% |
| 5 | 0 | 4.70% |
| 6 | 0 | 3.12% |
| 7 | 0 | 5.18% |
| 9 | 0 | 4.74% |
| 11 | 0 | 3.46% |
| 12 | 0 | 4.73% |
| 13 | 0 | 2.91% |
| 15 | 0 | 11.67% |
| 17 | 0 | 4.17% |
| 20 | 0 | 3.47% |
| 23 | 0 | 3.93% |
| 24 | 0 | 6.39% |
