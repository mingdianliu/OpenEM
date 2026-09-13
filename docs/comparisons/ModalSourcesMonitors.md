# ModalSourcesMonitors

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/ModalSourcesMonitors/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Tutorial |
| Verdict | **agree** |
| Reference output | official archived output |
| Cells | 26 on the reference side, 26 on ours, 0 that did not line up |
| Numbers compared | 23 of 43 paired; the other 20 are near-zero, see below |
| Largest relative difference | 93.63% |
| Over 5% / over 20% | 3 / 1 |
| Figures | 11 on the reference side, 11 on ours |
| Largest pixel difference | 10.23% |

Machine-readable form of everything below: [`data/ModalSourcesMonitors.json`](data/ModalSourcesMonitors.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 3 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 4 | 1.934145e+14 0                 1.55  2.330600    0.0           0.0 | 1.93414e+14 | 1.93414e+14 | 0 |  |
| 4 | 1.934145e+14 0                 1.55  2.330600    0.0           0.0 | 0 | 0 | 0 | near-zero |
| 4 | 1.934145e+14 0                 1.55  2.330600    0.0           0.0 | 1.55 | 1.55 | 0 | near-zero |
| 4 | 1.934145e+14 0                 1.55  2.330600    0.0           0.0 | 2.3306 | 2.27697 | 2.77e-12% | near-zero |
| 4 | 1.934145e+14 0                 1.55  2.330600    0.0           0.0 | 0 | 0 | 0 | near-zero |
| 4 | 1.934145e+14 0                 1.55  2.330600    0.0           0.0 | 0 | 0 | 0 | near-zero |
| 4 | 1                 1.55  1.545621    0.0           0.0 | 1 | 1 | 0 |  |
| 4 | 1                 1.55  1.545621    0.0           0.0 | 1.55 | 1.55 | 0 |  |
| 4 | 1                 1.55  1.545621    0.0           0.0 | 1.54562 | 1.56277 | 1.11% |  |
| 4 | 1                 1.55  1.545621    0.0           0.0 | 0 | 0 | 0 | near-zero |
| 4 | 1                 1.55  1.545621    0.0           0.0 | 0 | 0 | 0 | near-zero |
| 4 | 2                 1.55  1.370211    0.0           0.0 | 2 | 2 | 0 |  |
| 4 | 2                 1.55  1.370211    0.0           0.0 | 1.55 | 1.55 | 0 |  |
| 4 | 2                 1.55  1.370211    0.0           0.0 | 1.37021 | 1.35585 | 1.05% |  |
| 4 | 2                 1.55  1.370211    0.0           0.0 | 0 | 0 | 0 | near-zero |
| 4 | 2                 1.55  1.370211    0.0           0.0 | 0 | 0 | 0 | near-zero |
| 4 | 1.934145e+14 0                   0.979003        0.708737        0.817743 | 1.93414e+14 | 1.93414e+14 | 0 |  |
| 4 | 1.934145e+14 0                   0.979003        0.708737        0.817743 | 0 | 0 | 0 | near-zero |
| 4 | 1.934145e+14 0                   0.979003        0.708737        0.817743 | 0.979003 | 0.969961 | 4.67e-13% | near-zero |
| 4 | 1.934145e+14 0                   0.979003        0.708737        0.817743 | 0.708737 | 0.68417 | 1.27e-12% | near-zero |
| 4 | 1.934145e+14 0                   0.979003        0.708737        0.817743 | 0.817743 | 0.817016 | 3.76e-14% | near-zero |
| 4 | 1                   0.068648        0.733946        0.882326 | 1 | 1 | 0 |  |
| 4 | 1                   0.068648        0.733946        0.882326 | 0.068648 | 0.07599 | 10.70% |  |
| 4 | 1                   0.068648        0.733946        0.882326 | 0.733946 | 0.737899 | 0.54% |  |
| 4 | 1                   0.068648        0.733946        0.882326 | 0.882326 | 0.878217 | 0.47% |  |
| 4 | 2                   0.012079        0.842603        0.936079 | 2 | 2 | 0 |  |
| 4 | 2                   0.012079        0.842603        0.936079 | 0.012079 | 0.012963 | 4.42% | near-zero |
| 4 | 2                   0.012079        0.842603        0.936079 | 0.842603 | 0.853585 | 1.30% |  |
| 4 | 2                   0.012079        0.842603        0.936079 | 0.936079 | 0.932849 | 0.35% |  |
| 4 | 1.934145e+14 0            0.186907 | 1.93414e+14 | 1.93414e+14 | 0 |  |
| 4 | 1.934145e+14 0            0.186907 | 0 | 0 | 0 | near-zero |
| 4 | 1.934145e+14 0            0.186907 | 0.186907 | 0.197228 | 5.34e-13% | near-zero |
| 4 | 1            0.407558 | 1 | 1 | 0 |  |
| 4 | 1            0.407558 | 0.407558 | 0.386959 | 5.05% |  |
| 4 | 2            1.784790 | 2 | 2 | 0 |  |
| 4 | 2            1.784790 | 1.78479 | 1.82588 | 2.30% |  |
| 10 | Flux at central frequency:  1.0000015 | 1 | 1.00001 | 9.64e-04% |  |
| 11 | positive dir.  [1.00000157e+00 6.10565813e-17 1.04270062e-17] | 1 | 1.00001 | 8.86e-04% |  |
| 11 | positive dir.  [1.00000157e+00 6.10565813e-17 1.04270062e-17] | 6.10566e-17 | 8.41633e-19 | 6.02e-13% | near-zero |
| 11 | positive dir.  [1.00000157e+00 6.10565813e-17 1.04270062e-17] | 1.0427e-17 | 2.39807e-20 | 1.04e-13% | near-zero |
| 11 | negative dir.  [1.49830280e-09 1.68449502e-18 3.02154446e-19] | 1.4983e-09 | 9.54956e-11 | 93.63% |  |
| 11 | negative dir.  [1.49830280e-09 1.68449502e-18 3.02154446e-19] | 1.6845e-18 | 3.03839e-19 | 9.21e-06% | near-zero |
| 11 | negative dir.  [1.49830280e-09 1.68449502e-18 3.02154446e-19] | 3.02154e-19 | 1.99143e-20 | 1.88e-06% | near-zero |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 3 | 0 | 4.23% |
| 4 | 0 | 5.01% |
| 8 | 0 | 5.31% |
| 9 | 0 | 5.01% |
| 11 | 0 | 6.93% |
| 12 | 0 | 5.42% |
| 15 | 0 | 5.43% |
| 17 | 0 | 5.25% |
| 18 | 0 | 8.22% |
| 22 | 0 | 10.23% |
| 24 | 0 | 8.94% |
