# Autograd1Intro

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd1Intro/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **agree** |
| Reference output | official archived output |
| Cells | 17 on the reference side, 17 on ours, 1 that did not line up |
| Numbers compared | 7 |
| Largest relative difference | 2.93% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 4 on the reference side, 4 on ours |
| Largest pixel difference | 3.98% |

Machine-readable form of everything below: [`data/Autograd1Intro.json`](data/Autograd1Intro.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 5 | 1.7015115293406988 | 1.70151 | 1.70151 | 0 |  |
| 6 | dgdx=1.0 | 1 | 1 | 0 |  |
| 6 | dgdx=1.0, dgdy=1.0, dgdz=2.0 | 1 | 1 | 0 |  |
| 6 | dgdx=1.0, dgdy=1.0, dgdz=2.0 | 1 | 1 | 0 |  |
| 6 | dgdx=1.0, dgdy=1.0, dgdz=2.0 | 2 | 2 | 0 |  |
| 16 | power = 0.552 | 0.552 | 0.547 | 0.91% |  |
| 16 | d_power/d_eps = -0.16138196406328512 | -0.161382 | -0.156655 | 2.93% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 2 | 0 | 2.18% |
| 3 | 0 | 3.07% |
| 4 | 0 | 3.06% |
| 9 | 0 | 3.98% |
