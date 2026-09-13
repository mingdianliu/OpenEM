# QMRSMetasurface

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/QMRSMetasurface/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 15 on the reference side, 15 on ours, 0 that did not line up |
| Numbers compared | none printed |
| Largest relative difference | - |
| Over 5% / over 20% | 0 / 0 |
| Figures | 7 on the reference side, 7 on ours |
| Largest pixel difference | 52.45% |

Machine-readable form of everything below: [`data/QMRSMetasurface.json`](data/QMRSMetasurface.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 3 | 0 | 52.45% |
| 5 | 0 | 9.95% |
| 7 | 0 | 10.99% |
| 9 | 0 | 11.67% |
| 11 | 0 | 11.89% |
| 13 | 0 | 13.38% |
| 14 | 0 | 34.23% |
