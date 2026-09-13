# SourceNormalization

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/SourceNormalization/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Tutorial |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 37 on the reference side, 37 on ours, 7 that did not line up |
| Numbers compared | 3 |
| Largest relative difference | 0.75% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 12 on the reference side, 12 on ours |
| Largest pixel difference | 6.73% |

Machine-readable form of everything below: [`data/SourceNormalization.json`](data/SourceNormalization.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 12 | Flux:			0.00134 [W] | 0.00134 | 0.00133 | 0.75% |  |
| 12 | Flux analytic:		0.00133 [W] | 0.00133 | 0.00133 | 0 |  |
| 12 | E0:			1.00918 [V/µm] | 1.00918 | 1.00616 | 0.30% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 4 | 0 | 0 |
| 7 | 0 | 0.04% |
| 8 | 0 | 0.39% |
| 15 | 0 | 0 |
| 19 | 0 | 2.17% |
| 20 | 0 | 6.73% |
| 21 | 0 | 6.19% |
| 25 | 0 | 3.86% |
| 27 | 0 | 0.02% |
| 30 | 0 | 0 |
| 33 | 0 | 1.81% |
| 36 | 0 | 1.40% |
