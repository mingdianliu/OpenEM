# PlasmonicNanoparticle

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/PlasmonicNanoparticle/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Tutorial |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 13 on the reference side, 13 on ours, 1 that did not line up |
| Numbers compared | none printed |
| Largest relative difference | - |
| Over 5% / over 20% | 0 / 0 |
| Figures | 6 on the reference side, 6 on ours |
| Largest pixel difference | 17.65% |

Machine-readable form of everything below: [`data/PlasmonicNanoparticle.json`](data/PlasmonicNanoparticle.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 5 | 0 | 0 |
| 10 | 0 | 17.65% |
| 11 | 0 | 0.33% |
| 12 | 0 | 4.03% |
| 12 | 1 | 1.45% |
| 12 | 2 | 9.97% |
