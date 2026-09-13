# Autograd30ParallelAdjoint

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd30ParallelAdjoint/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **agree** |
| Reference output | official archived output |
| Cells | 11 on the reference side, 11 on ours, 2 that did not line up |
| Numbers compared | 6 |
| Largest relative difference | 53.42% |
| Over 5% / over 20% | 1 / 1 |
| Figures | 1 on the reference side, 1 on ours |
| Largest pixel difference | 4.03% |

Machine-readable form of everything below: [`data/Autograd30ParallelAdjoint.json`](data/Autograd30ParallelAdjoint.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 1 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 7 | Targeting the mode that appeared as mode_index=2 in the initial solve. | 2 | 2 | 0 |  |
| 7 | Using ModeSpec(num_modes=1, target_neff=1.8282) for the differentiable run. | 1 | 1 | 0 |  |
| 7 | Using ModeSpec(num_modes=1, target_neff=1.8282) for the differentiable run. | 1.8282 | 1.8258 | 0.13% |  |
| 10 | parallel adjoint        3.8935e-01    3.5748e-01         35.21 | 0.38935 | 0.38929 | 0.02% |  |
| 10 | parallel adjoint        3.8935e-01    3.5748e-01         35.21 | 0.35748 | 0.35802 | 0.15% |  |
| 10 | parallel adjoint        3.8935e-01    3.5748e-01         35.21 | 35.21 | 16.4 | 53.42% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 5 | 0 | 4.03% |
