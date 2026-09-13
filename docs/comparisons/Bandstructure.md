# Bandstructure

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Bandstructure/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 16 on the reference side, 16 on ours, 1 that did not line up |
| Numbers compared | 2 of 7 paired; the other 5 are near-zero, see below |
| Largest relative difference | 318% |
| Over 5% / over 20% | 1 / 1 |
| Figures | 5 on the reference side, 5 on ours |
| Largest pixel difference | 1.15% |

Machine-readable form of everything below: [`data/Bandstructure.json`](data/Bandstructure.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 1 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 2 | Total runtime = 6.74 ps | 6.74 | 6.74 | 0 |  |
| 12 | 2.559011e+13  1.159709e+08  693222.812687   0.001382  1.241659  0.000016 | 2.55901e+13 | 1.06897e+14 | 318% |  |
| 12 | 2.559011e+13  1.159709e+08  693222.812687   0.001382  1.241659  0.000016 | 1.15971e+08 | 1.16898e+11 | 45.64% | near-zero |
| 12 | 2.559011e+13  1.159709e+08  693222.812687   0.001382  1.241659  0.000016 | 693223 | 2872.82 | 2.70e-04% | near-zero |
| 12 | 2.559011e+13  1.159709e+08  693222.812687   0.001382  1.241659  0.000016 | 0.001382 | 0.012521 | 4.35e-12% | near-zero |
| 12 | 2.559011e+13  1.159709e+08  693222.812687   0.001382  1.241659  0.000016 | 1.24166 | 2.52375 | 5.01e-10% | near-zero |
| 12 | 2.559011e+13  1.159709e+08  693222.812687   0.001382  1.241659  0.000016 | 1.6e-05 | 1.5e-05 | 3.91e-16% | near-zero |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 8 | 0 | 0 |
| 8 | 1 | 0 |
| 10 | 0 | 1.15% |
| 11 | 0 | 0.71% |
| 15 | 0 | 0.06% |
