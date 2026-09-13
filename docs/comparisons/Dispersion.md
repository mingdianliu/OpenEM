# Dispersion

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Dispersion/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Tutorial |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 19 on the reference side, 19 on ours, 2 that did not line up |
| Numbers compared | none printed |
| Largest relative difference | - |
| Over 5% / over 20% | 0 / 0 |
| Figures | 7 on the reference side, 7 on ours |
| Largest pixel difference | 1.29% |

Machine-readable form of everything below: [`data/Dispersion.json`](data/Dispersion.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 8 | 0 | 0 |
| 9 | 0 | 0 |
| 10 | 0 | 0 |
| 10 | 1 | 0 |
| 12 | 0 | 1.16% |
| 14 | 0 | 1.29% |
| 17 | 0 | 1.26% |
