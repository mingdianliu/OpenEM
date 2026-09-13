# EffectiveIndexApproximation

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/EffectiveIndexApproximation/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 23 on the reference side, 23 on ours, 1 that did not line up |
| Numbers compared | 12 |
| Largest relative difference | 3.32% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 10 on the reference side, 10 on ours |
| Largest pixel difference | 8.34% |

Machine-readable form of everything below: [`data/EffectiveIndexApproximation.json`](data/EffectiveIndexApproximation.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 17 | Mean FSR (nm): 21.083 | 21.083 | 21 | 0.39% |  |
| 17 | Resonance dip wavelengths (µm): [1.581   1.562   1.543   1.52467 1.50667] | 1.581 | 1.581 | 0 |  |
| 17 | Resonance dip wavelengths (µm): [1.581   1.562   1.543   1.52467 1.50667] | 1.562 | 1.562 | 0 |  |
| 17 | Resonance dip wavelengths (µm): [1.581   1.562   1.543   1.52467 1.50667] | 1.543 | 1.543 | 0 |  |
| 17 | Resonance dip wavelengths (µm): [1.581   1.562   1.543   1.52467 1.50667] | 1.52467 | 1.52467 | 0 |  |
| 17 | Resonance dip wavelengths (µm): [1.581   1.562   1.543   1.52467 1.50667] | 1.50667 | 1.50667 | 0 |  |
| 17 | FSR values (nm): [19.    19.    18.333 18.   ] | 19 | 19 | 0 |  |
| 17 | FSR values (nm): [19.    19.    18.333 18.   ] | 19 | 19 | 0 |  |
| 17 | FSR values (nm): [19.    19.    18.333 18.   ] | 18.333 | 18.333 | 0 |  |
| 17 | FSR values (nm): [19.    19.    18.333 18.   ] | 18 | 18 | 0 |  |
| 17 | Mean FSR (nm): 18.583 | 18.583 | 18.583 | 0 |  |
| 17 | Mean FSR difference: 2.500 nm | 2.5 | 2.417 | 3.32% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 5 | 0 | 4.39% |
| 6 | 0 | 6.18% |
| 10 | 0 | 7.47% |
| 11 | 0 | 5.20% |
| 12 | 0 | 8.18% |
| 14 | 0 | 5.48% |
| 16 | 0 | 8.34% |
| 18 | 0 | 7.98% |
| 19 | 0 | 5.50% |
| 21 | 0 | 5.26% |
