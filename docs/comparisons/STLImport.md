# STLImport

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/STLImport/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Tutorial |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 15 on the reference side, 15 on ours, 0 that did not line up |
| Numbers compared | none printed |
| Largest relative difference | - |
| Over 5% / over 20% | 0 / 0 |
| Figures | 11 on the reference side, 11 on ours |
| Largest pixel difference | 25.13% |

Machine-readable form of everything below: [`data/STLImport.json`](data/STLImport.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 3 | 0 | 0 |
| 4 | 0 | 0 |
| 7 | 0 | 0 |
| 8 | 0 | 25.12% |
| 8 | 1 | 25.13% |
| 9 | 0 | 0 |
| 9 | 1 | 0 |
| 10 | 0 | 0 |
| 12 | 0 | 0 |
| 13 | 0 | 0.17% |
| 13 | 1 | 0.17% |
