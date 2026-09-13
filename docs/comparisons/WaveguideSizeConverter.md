# WaveguideSizeConverter

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/WaveguideSizeConverter/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 26 on the reference side, 26 on ours, 0 that did not line up |
| Numbers compared | none printed |
| Largest relative difference | - |
| Over 5% / over 20% | 0 / 0 |
| Figures | 12 on the reference side, 12 on ours |
| Largest pixel difference | 7.86% |

Machine-readable form of everything below: [`data/WaveguideSizeConverter.json`](data/WaveguideSizeConverter.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 5 | 0 | 3.31% |
| 7 | 0 | 6.98% |
| 8 | 0 | 7.07% |
| 10 | 0 | 6.42% |
| 11 | 0 | 4.99% |
| 14 | 0 | 4.46% |
| 16 | 0 | 7.86% |
| 17 | 0 | 6.74% |
| 19 | 0 | 3.27% |
| 21 | 0 | 4.46% |
| 23 | 0 | 4.79% |
| 24 | 0 | 4.95% |
