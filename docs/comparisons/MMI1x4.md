# MMI1x4

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/MMI1x4/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 16 on the reference side, 16 on ours, 1 that did not line up |
| Numbers compared | none printed |
| Largest relative difference | - |
| Over 5% / over 20% | 0 / 0 |
| Figures | 9 on the reference side, 9 on ours |
| Largest pixel difference | 16.42% |

Machine-readable form of everything below: [`data/MMI1x4.json`](data/MMI1x4.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 3 | 0 | 4.78% |
| 5 | 0 | 6.10% |
| 7 | 0 | 9.33% |
| 8 | 0 | 6.56% |
| 9 | 0 | 5.10% |
| 12 | 0 | 16.42% |
| 13 | 0 | 3.67% |
| 14 | 0 | 5.77% |
| 15 | 0 | 5.63% |
