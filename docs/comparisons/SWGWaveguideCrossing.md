# SWGWaveguideCrossing

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/SWGWaveguideCrossing/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **close** |
| Reference output | reference rerun |
| Cells | 21 on the reference side, 21 on ours, 0 that did not line up |
| Numbers compared | 8 |
| Largest relative difference | 6.81% |
| Over 5% / over 20% | 2 / 0 |
| Figures | 3 on the reference side, 3 on ours |
| Largest pixel difference | 9.35% |

Machine-readable form of everything below: [`data/SWGWaveguideCrossing.json`](data/SWGWaveguideCrossing.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 2 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 1 | Effective index in regions B: 3.018, C: 3.131, D: 2.567. | 3.018 | 3.018 | 0 |  |
| 1 | Effective index in regions B: 3.018, C: 3.131, D: 2.567. | 3.131 | 3.131 | 0 |  |
| 20 | Average transmittance (dB) for through mode, TE0: -0.269, TE1: -0.343, TE2: -0.436 | -0.269 | -0.265 | 1.49% |  |
| 20 | Average transmittance (dB) for through mode, TE0: -0.269, TE1: -0.343, TE2: -0.436 | -0.343 | -0.331 | 3.50% |  |
| 20 | Average transmittance (dB) for through mode, TE0: -0.269, TE1: -0.343, TE2: -0.436 | -0.436 | -0.414 | 5.05% |  |
| 20 | Average transmittance (dB) into cross mode, TE0: -48.640, TE1: -44.789, TE2: -37.736 | -48.64 | -45.326 | 6.81% |  |
| 20 | Average transmittance (dB) into cross mode, TE0: -48.640, TE1: -44.789, TE2: -37.736 | -44.789 | -43.613 | 2.63% |  |
| 20 | Average transmittance (dB) into cross mode, TE0: -48.640, TE1: -44.789, TE2: -37.736 | -37.736 | -38.271 | 1.42% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 14 | 0 | 4.19% |
| 18 | 0 | 9.32% |
| 19 | 0 | 9.35% |
