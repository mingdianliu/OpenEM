# CavityFOM

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/CavityFOM/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Tutorial |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 14 on the reference side, 14 on ours, 1 that did not line up |
| Numbers compared | 8 of 15 paired; the other 7 are near-zero, see below |
| Largest relative difference | 0.04% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 3 on the reference side, 3 on ours |
| Largest pixel difference | 0.18% |

Machine-readable form of everything below: [`data/CavityFOM.json`](data/CavityFOM.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 1 | Total runtime = 5.52 ps | 5.52 | 5.52 | 0 |  |
| 1 | Start monitoring fields after 0.55 ps | 0.55 | 0.55 | 0 |  |
| 10 | 3.297092e+14  1.091220e+10  94922.365484  144631.787403 -2.015378  0.000069 | 3.29709e+14 | 3.2971e+14 | 9.10e-05% |  |
| 10 | 3.297092e+14  1.091220e+10  94922.365484  144631.787403 -2.015378  0.000069 | 1.09122e+10 | 1.09175e+10 | 1.60e-04% | near-zero |
| 10 | 3.297092e+14  1.091220e+10  94922.365484  144631.787403 -2.015378  0.000069 | 94922.4 | 94876.5 | 1.39e-09% | near-zero |
| 10 | 3.297092e+14  1.091220e+10  94922.365484  144631.787403 -2.015378  0.000069 | 144632 | 144633 | 2.99e-11% | near-zero |
| 10 | 3.297092e+14  1.091220e+10  94922.365484  144631.787403 -2.015378  0.000069 | -2.01538 | -2.01593 | 1.69e-14% | near-zero |
| 10 | 3.297092e+14  1.091220e+10  94922.365484  144631.787403 -2.015378  0.000069 | 6.9e-05 | 7.8e-05 | 2.73e-16% | near-zero |
| 11 | V_eff = 0.014 um^3 | 0.014 | 0.014 | 0 | near-zero |
| 11 | V_eff = 0.014 um^3 | 3 | 3 | 0 |  |
| 11 | V_eff = 1.431e-20 m^3 | 1.431e-20 | 1.431e-20 | 0 | near-zero |
| 11 | V_eff = 1.431e-20 m^3 | 3 | 3 | 0 |  |
| 11 | V_eff = 0.79 (lambda/n)^3 | 0.79 | 0.79 | 0 |  |
| 11 | V_eff = 0.79 (lambda/n)^3 | 3 | 3 | 0 |  |
| 12 | F_p = 9173 | 9173 | 9169 | 0.04% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 7 | 0 | 0 |
| 9 | 0 | 0.18% |
| 13 | 0 | 0.02% |
