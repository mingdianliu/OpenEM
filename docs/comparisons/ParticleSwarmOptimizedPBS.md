# ParticleSwarmOptimizedPBS

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/ParticleSwarmOptimizedPBS/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **trajectory divergence** |
| Reference output | reference rerun |
| Cells | 19 on the reference side, 19 on ours, 1 that did not line up |
| Numbers compared | 19 of 20 paired; the other 1 are near-zero, see below |
| Largest relative difference | 33.13% |
| Over 5% / over 20% | 7 / 1 |
| Figures | 6 on the reference side, 6 on ours |
| Largest pixel difference | 14.89% |

Machine-readable form of everything below: [`data/ParticleSwarmOptimizedPBS.json`](data/ParticleSwarmOptimizedPBS.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 7 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 14 | Best Parameters: 0: 0.39604210811471385 1: 0.3419827901237928 2: | 0 | 0 | 0 | near-zero |
| 14 | Best Parameters: 0: 0.39604210811471385 1: 0.3419827901237928 2: | 0.396042 | 0.365695 | 7.66% |  |
| 14 | Best Parameters: 0: 0.39604210811471385 1: 0.3419827901237928 2: | 1 | 1 | 0 |  |
| 14 | Best Parameters: 0: 0.39604210811471385 1: 0.3419827901237928 2: | 0.341983 | 0.312229 | 8.70% |  |
| 14 | Best Parameters: 0: 0.39604210811471385 1: 0.3419827901237928 2: | 2 | 2 | 0 |  |
| 14 | 0.3429872664813658 3: 0.4484520145366082 4: 0.3681504370931801 5: | 0.342987 | 0.409119 | 19.28% |  |
| 14 | 0.3429872664813658 3: 0.4484520145366082 4: 0.3681504370931801 5: | 3 | 3 | 0 |  |
| 14 | 0.3429872664813658 3: 0.4484520145366082 4: 0.3681504370931801 5: | 0.448452 | 0.428283 | 4.50% |  |
| 14 | 0.3429872664813658 3: 0.4484520145366082 4: 0.3681504370931801 5: | 4 | 4 | 0 |  |
| 14 | 0.3429872664813658 3: 0.4484520145366082 4: 0.3681504370931801 5: | 0.36815 | 0.388036 | 5.40% |  |
| 14 | 0.3429872664813658 3: 0.4484520145366082 4: 0.3681504370931801 5: | 5 | 5 | 0 |  |
| 14 | 0.3335669914929738 6: 0.3689782147929135 7: 0.3230366586846238 8: | 0.333567 | 0.322397 | 3.35% |  |
| 14 | 0.3335669914929738 6: 0.3689782147929135 7: 0.3230366586846238 8: | 6 | 6 | 0 |  |
| 14 | 0.3335669914929738 6: 0.3689782147929135 7: 0.3230366586846238 8: | 0.368978 | 0.337337 | 8.58% |  |
| 14 | 0.3335669914929738 6: 0.3689782147929135 7: 0.3230366586846238 8: | 7 | 7 | 0 |  |
| 14 | 0.3335669914929738 6: 0.3689782147929135 7: 0.3230366586846238 8: | 0.323037 | 0.323365 | 0.10% |  |
| 14 | 0.3335669914929738 6: 0.3689782147929135 7: 0.3230366586846238 8: | 8 | 8 | 0 |  |
| 14 | 0.3040949366805732 9: 0.3379306749089194 | 0.304095 | 0.340036 | 11.82% |  |
| 14 | 0.3040949366805732 9: 0.3379306749089194 | 9 | 9 | 0 |  |
| 14 | 0.3040949366805732 9: 0.3379306749089194 | 0.337931 | 0.449888 | 33.13% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 5 | 0 | 5.53% |
| 7 | 0 | 8.40% |
| 10 | 0 | 4.84% |
| 15 | 0 | 4.81% |
| 17 | 0 | 14.89% |
| 18 | 0 | 4.74% |
