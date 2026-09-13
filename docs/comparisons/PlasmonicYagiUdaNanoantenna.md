# PlasmonicYagiUdaNanoantenna

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/PlasmonicYagiUdaNanoantenna/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 17 on the reference side, 17 on ours, 4 that did not line up |
| Numbers compared | none printed |
| Largest relative difference | - |
| Over 5% / over 20% | 0 / 0 |
| Figures | 7 on the reference side, 7 on ours |
| Largest pixel difference | 0.94% |

Machine-readable form of everything below: [`data/PlasmonicYagiUdaNanoantenna.json`](data/PlasmonicYagiUdaNanoantenna.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 4 | 0 | 0 |
| 7 | 0 | 0.89% |
| 8 | 0 | 0.35% |
| 9 | 0 | 0 |
| 11 | 0 | 0.94% |
| 13 | 0 | 0 |
| 15 | 0 | 0.59% |
