# Autograd22PhotonicCrystal

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd22PhotonicCrystal/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **trajectory divergence** |
| Reference output | official archived output |
| Cells | 27 on the reference side, 27 on ours, 0 that did not line up |
| Numbers compared | 31 |
| Largest relative difference | 68.12% |
| Over 5% / over 20% | 9 / 5 |
| Figures | 8 on the reference side, 8 on ours |
| Largest pixel difference | 12.38% |

Machine-readable form of everything below: [`data/Autograd22PhotonicCrystal.json`](data/Autograd22PhotonicCrystal.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 9 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 17 | 0.9999248769975477 | 0.999925 | 1.00204 | 0.21% |  |
| 21 | step = 1 | 1 | 1 | 0 |  |
| 21 | J = 5.0388e-01 | 0.50388 | 0.50031 | 0.71% |  |
| 21 | grad_norm = 4.5922e+00 | 4.5922 | 4.766 | 3.78% |  |
| 21 | step = 2 | 2 | 2 | 0 |  |
| 21 | J = 7.1618e-01 | 0.71618 | 0.71747 | 0.18% |  |
| 21 | grad_norm = 2.2663e+00 | 2.2663 | 2.3975 | 5.79% |  |
| 21 | step = 3 | 3 | 3 | 0 |  |
| 21 | J = 8.0503e-01 | 0.80503 | 0.7977 | 0.91% |  |
| 21 | grad_norm = 1.6398e+00 | 1.6398 | 1.8311 | 11.67% |  |
| 21 | step = 4 | 4 | 4 | 0 |  |
| 21 | J = 8.3990e-01 | 0.8399 | 0.84398 | 0.49% |  |
| 21 | grad_norm = 1.7218e+00 | 1.7218 | 1.8468 | 7.26% |  |
| 21 | step = 5 | 5 | 5 | 0 |  |
| 21 | J = 8.6579e-01 | 0.86579 | 0.89088 | 2.90% |  |
| 21 | grad_norm = 1.9119e+00 | 1.9119 | 1.4731 | 22.95% |  |
| 21 | step = 6 | 6 | 6 | 0 |  |
| 21 | J = 9.0569e-01 | 0.90569 | 0.92014 | 1.60% |  |
| 21 | grad_norm = 1.5462e+00 | 1.5462 | 0.68693 | 55.57% |  |
| 21 | step = 7 | 7 | 7 | 0 |  |
| 21 | J = 9.3005e-01 | 0.93005 | 0.91656 | 1.45% |  |
| 21 | grad_norm = 8.6348e-01 | 0.86348 | 1.394 | 61.44% |  |
| 21 | step = 8 | 8 | 8 | 0 |  |
| 21 | J = 9.3529e-01 | 0.93529 | 0.92282 | 1.33% |  |
| 21 | grad_norm = 9.1394e-01 | 0.91394 | 1.5365 | 68.12% |  |
| 21 | step = 9 | 9 | 9 | 0 |  |
| 21 | J = 9.3672e-01 | 0.93672 | 0.93914 | 0.26% |  |
| 21 | grad_norm = 1.3241e+00 | 1.3241 | 1.1713 | 11.54% |  |
| 21 | step = 10 | 10 | 10 | 0 |  |
| 21 | J = 9.4567e-01 | 0.94567 | 0.9568 | 1.18% |  |
| 21 | grad_norm = 1.3289e+00 | 1.3289 | 0.67052 | 49.54% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 12 | 0 | 12.38% |
| 14 | 0 | 4.06% |
| 15 | 0 | 5.17% |
| 16 | 0 | 5.17% |
| 22 | 0 | 3.60% |
| 23 | 0 | 5.90% |
| 25 | 0 | 4.50% |
| 26 | 0 | 8.95% |
