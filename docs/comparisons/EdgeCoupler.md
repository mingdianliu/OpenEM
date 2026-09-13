# EdgeCoupler

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/EdgeCoupler/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 22 on the reference side, 22 on ours, 1 that did not line up |
| Numbers compared | none printed |
| Largest relative difference | - |
| Over 5% / over 20% | 0 / 0 |
| Figures | 8 on the reference side, 8 on ours |
| Largest pixel difference | 7.76% |

Machine-readable form of everything below: [`data/EdgeCoupler.json`](data/EdgeCoupler.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 7 | 0 | 4.55% |
| 8 | 0 | 7.76% |
| 12 | 0 | 5.07% |
| 13 | 0 | 3.47% |
| 14 | 0 | 5.78% |
| 19 | 0 | 6.06% |
| 20 | 0 | 5.23% |
| 21 | 0 | 6.41% |
