# SbendCMAES

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/SbendCMAES/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 18 on the reference side, 18 on ours, 0 that did not line up |
| Numbers compared | 18 of 24 paired; the other 6 are near-zero, see below |
| Largest relative difference | 1.22% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 5 on the reference side, 5 on ours |
| Largest pixel difference | 8.69% |

Machine-readable form of everything below: [`data/SbendCMAES.json`](data/SbendCMAES.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 10 | Generation 0: Best function value = 0.271 | 0 | 0 | 0 | near-zero |
| 10 | Generation 0: Best function value = 0.271 | 0.271 | 0.271 | 0 |  |
| 10 | Generation 1: Best function value = 0.265 | 1 | 1 | 0 |  |
| 10 | Generation 1: Best function value = 0.265 | 0.265 | 0.265 | 0 |  |
| 10 | Generation 2: Best function value = 0.229 | 2 | 2 | 0 |  |
| 10 | Generation 2: Best function value = 0.229 | 0.229 | 0.229 | 0 |  |
| 10 | Generation 3: Best function value = 0.216 | 3 | 3 | 0 |  |
| 10 | Generation 3: Best function value = 0.216 | 0.216 | 0.217 | 0.46% |  |
| 10 | Generation 4: Best function value = 0.173 | 4 | 4 | 0 |  |
| 10 | Generation 4: Best function value = 0.173 | 0.173 | 0.174 | 0.58% |  |
| 10 | Generation 5: Best function value = 0.119 | 5 | 5 | 0 |  |
| 10 | Generation 5: Best function value = 0.119 | 0.119 | 0.12 | 0.84% |  |
| 10 | Generation 6: Best function value = 0.082 | 6 | 6 | 0 |  |
| 10 | Generation 6: Best function value = 0.082 | 0.082 | 0.083 | 1.22% |  |
| 10 | Generation 7: Best function value = 0.064 | 7 | 7 | 0 |  |
| 10 | Generation 7: Best function value = 0.064 | 0.064 | 0.065 | 1.43% | near-zero |
| 10 | Generation 8: Best function value = 0.064 | 8 | 8 | 0 |  |
| 10 | Generation 8: Best function value = 0.064 | 0.064 | 0.065 | 1.25% | near-zero |
| 10 | Generation 9: Best function value = 0.061 | 9 | 9 | 0 |  |
| 10 | Generation 9: Best function value = 0.061 | 0.061 | 0.061 | 0 | near-zero |
| 10 | Generation 10: Best function value = 0.044 | 10 | 10 | 0 |  |
| 10 | Generation 10: Best function value = 0.044 | 0.044 | 0.044 | 0 | near-zero |
| 10 | Generation 11: Best function value = 0.040 | 11 | 11 | 0 |  |
| 10 | Generation 11: Best function value = 0.040 | 0.04 | 0.041 | 0.91% | near-zero |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 4 | 0 | 2.50% |
| 5 | 0 | 8.69% |
| 11 | 0 | 3.32% |
| 14 | 0 | 2.59% |
| 15 | 0 | 4.44% |
