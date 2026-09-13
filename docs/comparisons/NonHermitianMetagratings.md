# NonHermitianMetagratings

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/NonHermitianMetagratings/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 18 on the reference side, 18 on ours, 1 that did not line up |
| Numbers compared | 1 |
| Largest relative difference | 0.03% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 6 on the reference side, 6 on ours |
| Largest pixel difference | 14.22% |

Machine-readable form of everything below: [`data/NonHermitianMetagratings.json`](data/NonHermitianMetagratings.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 11 | C_exc = 0.9816 | 0.9816 | 0.9819 | 0.03% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 3 | 0 | 5.46% |
| 7 | 0 | 3.40% |
| 8 | 0 | 3.33% |
| 10 | 0 | 5.71% |
| 15 | 0 | 14.22% |
| 16 | 0 | 2.43% |
