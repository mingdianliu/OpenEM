# DirectionalCoupler

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/DirectionalCoupler/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 18 on the reference side, 18 on ours, 1 that did not line up |
| Numbers compared | 2 |
| Largest relative difference | 0.11% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 5 on the reference side, 5 on ours |
| Largest pixel difference | 5.26% |

Machine-readable form of everything below: [`data/DirectionalCoupler.json`](data/DirectionalCoupler.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 10 | Crossing wavelength :1.53 | 1.53 | 1.53 | 0 |  |
| 10 | Total transmittance: 0.9910 | 0.991 | 0.9899 | 0.11% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 6 | 0 | 4.64% |
| 10 | 0 | 3.06% |
| 11 | 0 | 2.37% |
| 15 | 0 | 5.26% |
| 17 | 0 | 2.40% |
