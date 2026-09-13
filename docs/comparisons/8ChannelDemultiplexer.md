# 8ChannelDemultiplexer

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/8ChannelDemultiplexer/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 44 on the reference side, 44 on ours, 1 that did not line up |
| Numbers compared | none printed |
| Largest relative difference | - |
| Over 5% / over 20% | 0 / 0 |
| Figures | 17 on the reference side, 17 on ours |
| Largest pixel difference | 7.66% |

Machine-readable form of everything below: [`data/8ChannelDemultiplexer.json`](data/8ChannelDemultiplexer.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 7 | 0 | 7.66% |
| 9 | 0 | 6.71% |
| 14 | 0 | 3.02% |
| 16 | 0 | 7.66% |
| 19 | 0 | 7.08% |
| 21 | 0 | 7.29% |
| 23 | 0 | 7.10% |
| 25 | 0 | 7.28% |
| 27 | 0 | 7.19% |
| 29 | 0 | 7.03% |
| 33 | 0 | 7.28% |
| 35 | 0 | 4.30% |
| 38 | 0 | 4.11% |
| 39 | 0 | 3.82% |
| 40 | 0 | 4.27% |
| 42 | 0 | 4.02% |
| 43 | 0 | 4.46% |
