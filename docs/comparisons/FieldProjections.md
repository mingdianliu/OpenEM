# FieldProjections

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/FieldProjections/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Tutorial |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 38 on the reference side, 38 on ours, 0 that did not line up |
| Numbers compared | none printed |
| Largest relative difference | - |
| Over 5% / over 20% | 0 / 0 |
| Figures | 18 on the reference side, 18 on ours |
| Largest pixel difference | 18.02% |

Machine-readable form of everything below: [`data/FieldProjections.json`](data/FieldProjections.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 4 | 0 | 0 |
| 9 | 0 | 0.16% |
| 12 | 0 | 0.16% |
| 16 | 0 | 0.17% |
| 20 | 0 | 0.96% |
| 20 | 1 | 0.20% |
| 21 | 0 | 18.02% |
| 21 | 1 | 1.33% |
| 22 | 0 | 0.24% |
| 23 | 0 | 0 |
| 29 | 0 | 7.80% |
| 29 | 1 | 8.22% |
| 29 | 2 | 5.72% |
| 30 | 0 | 0 |
| 32 | 0 | 7.80% |
| 32 | 1 | 8.13% |
| 34 | 0 | 0 |
| 36 | 0 | 0.17% |
