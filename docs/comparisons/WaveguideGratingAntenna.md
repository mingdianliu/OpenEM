# WaveguideGratingAntenna

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/WaveguideGratingAntenna/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 15 on the reference side, 15 on ours, 0 that did not line up |
| Numbers compared | 1 |
| Largest relative difference | 0 |
| Over 5% / over 20% | 0 / 0 |
| Figures | 3 on the reference side, 3 on ours |
| Largest pixel difference | 4.76% |

Machine-readable form of everything below: [`data/WaveguideGratingAntenna.json`](data/WaveguideGratingAntenna.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 10 | The optimal offset is 170.0 nm | 170 | 170 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 10 | 0 | 4.76% |
| 12 | 0 | 3.80% |
| 13 | 0 | 3.19% |
