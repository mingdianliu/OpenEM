# BilevelPSR

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/BilevelPSR/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 19 on the reference side, 19 on ours, 4 that did not line up |
| Numbers compared | none printed |
| Largest relative difference | - |
| Over 5% / over 20% | 0 / 0 |
| Figures | 10 on the reference side, 10 on ours |
| Largest pixel difference | 16.40% |

Machine-readable form of everything below: [`data/BilevelPSR.json`](data/BilevelPSR.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 6 | 0 | 3.65% |
| 7 | 0 | 8.68% |
| 8 | 0 | 16.40% |
| 9 | 0 | 15.58% |
| 10 | 0 | 15.68% |
| 11 | 0 | 14.16% |
| 15 | 0 | 8.18% |
| 16 | 0 | 4.55% |
| 17 | 0 | 9.58% |
| 18 | 0 | 4.50% |
