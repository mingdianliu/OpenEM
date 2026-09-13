# Autograd24DigitalSplitter

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd24DigitalSplitter/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **trajectory divergence** |
| Reference output | official archived output |
| Cells | 23 on the reference side, 23 on ours, 0 that did not line up |
| Numbers compared | 125 |
| Largest relative difference | 31.15% |
| Over 5% / over 20% | 74 / 49 |
| Figures | 19 on the reference side, 19 on ours |
| Largest pixel difference | 12.98% |

Machine-readable form of everything below: [`data/Autograd24DigitalSplitter.json`](data/Autograd24DigitalSplitter.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 74 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 10 | Iteration: 1 | 1 | 1 | 0 |  |
| 10 | FOM = 1.113e-01 | 0.1113 | 0.1107 | 0.54% |  |
| 10 | Iteration: 2 | 2 | 2 | 0 |  |
| 10 | FOM = 1.574e-01 | 0.1574 | 0.1475 | 6.29% |  |
| 10 | Iteration: 3 | 3 | 3 | 0 |  |
| 10 | FOM = 2.034e-01 | 0.2034 | 0.1846 | 9.24% |  |
| 10 | Iteration: 4 | 4 | 4 | 0 |  |
| 10 | FOM = 2.458e-01 | 0.2458 | 0.2194 | 10.74% |  |
| 10 | Iteration: 5 | 5 | 5 | 0 |  |
| 10 | FOM = 2.822e-01 | 0.2822 | 0.248 | 12.12% |  |
| 10 | Iteration: 6 | 6 | 6 | 0 |  |
| 10 | FOM = 3.130e-01 | 0.313 | 0.2661 | 14.98% |  |
| 10 | Iteration: 7 | 7 | 7 | 0 |  |
| 10 | FOM = 3.397e-01 | 0.3397 | 0.2746 | 19.16% |  |
| 10 | Iteration: 8 | 8 | 8 | 0 |  |
| 10 | FOM = 3.622e-01 | 0.3622 | 0.2805 | 22.56% |  |
| 10 | Iteration: 9 | 9 | 9 | 0 |  |
| 10 | FOM = 3.807e-01 | 0.3807 | 0.288 | 24.35% |  |
| 10 | Iteration: 10 | 10 | 10 | 0 |  |
| 10 | FOM = 3.961e-01 | 0.3961 | 0.2952 | 25.47% |  |
| 10 | Iteration: 11 | 11 | 11 | 0 |  |
| 10 | FOM = 4.086e-01 | 0.4086 | 0.3001 | 26.55% |  |
| 10 | Iteration: 12 | 12 | 12 | 0 |  |
| 10 | FOM = 4.176e-01 | 0.4176 | 0.3052 | 26.92% |  |
| 10 | Iteration: 13 | 13 | 13 | 0 |  |
| 10 | FOM = 4.253e-01 | 0.4253 | 0.31 | 27.11% |  |
| 10 | Iteration: 14 | 14 | 14 | 0 |  |
| 10 | FOM = 4.320e-01 | 0.432 | 0.3143 | 27.25% |  |
| 10 | Iteration: 15 | 15 | 15 | 0 |  |
| 10 | FOM = 4.380e-01 | 0.438 | 0.3171 | 27.60% |  |
| 10 | Iteration: 16 | 16 | 16 | 0 |  |
| 10 | FOM = 4.430e-01 | 0.443 | 0.3199 | 27.79% |  |
| 10 | Iteration: 17 | 17 | 17 | 0 |  |
| 10 | FOM = 4.466e-01 | 0.4466 | 0.3225 | 27.79% |  |
| 10 | Iteration: 18 | 18 | 18 | 0 |  |
| 10 | FOM = 4.497e-01 | 0.4497 | 0.3243 | 27.89% |  |
| 10 | Iteration: 19 | 19 | 19 | 0 |  |
| 10 | FOM = 4.524e-01 | 0.4524 | 0.3265 | 27.83% |  |
| 10 | Iteration: 20 | 20 | 20 | 0 |  |
| 10 | FOM = 4.548e-01 | 0.4548 | 0.3285 | 27.77% |  |
| 10 | Iteration: 21 | 21 | 21 | 0 |  |
| 10 | FOM = 4.570e-01 | 0.457 | 0.3307 | 27.64% |  |
| 10 | Iteration: 22 | 22 | 22 | 0 |  |
| 10 | FOM = 4.589e-01 | 0.4589 | 0.3326 | 27.52% |  |
| 10 | Iteration: 23 | 23 | 23 | 0 |  |
| 10 | FOM = 4.605e-01 | 0.4605 | 0.3343 | 27.40% |  |
| 10 | Iteration: 24 | 24 | 24 | 0 |  |
| 10 | FOM = 4.618e-01 | 0.4618 | 0.336 | 27.24% |  |
| 10 | Iteration: 25 | 25 | 25 | 0 |  |
| 10 | FOM = 4.628e-01 | 0.4628 | 0.3376 | 27.05% |  |
| 12 | Iteration: 1 | 1 | 1 | 0 |  |
| 12 | FOM = 4.638e-01 | 0.4638 | 0.3389 | 26.93% |  |
| 12 | MSE = 6.909e+00 | 6.909 | 6.113 | 11.52% |  |
| 12 | Iteration: 2 | 2 | 2 | 0 |  |
| 12 | FOM = 4.650e-01 | 0.465 | 0.3414 | 26.58% |  |
| 12 | MSE = 6.395e+00 | 6.395 | 5.622 | 12.09% |  |
| 12 | Iteration: 3 | 3 | 3 | 0 |  |
| 12 | FOM = 4.663e-01 | 0.4663 | 0.3436 | 26.31% |  |
| 12 | MSE = 6.016e+00 | 6.016 | 5.175 | 13.98% |  |
| 12 | Iteration: 4 | 4 | 4 | 0 |  |
| 12 | FOM = 4.672e-01 | 0.4672 | 0.346 | 25.94% |  |
| 12 | MSE = 5.628e+00 | 5.628 | 4.776 | 15.14% |  |
| 12 | Iteration: 5 | 5 | 5 | 0 |  |
| 12 | FOM = 4.681e-01 | 0.4681 | 0.3479 | 25.68% |  |
| 12 | MSE = 5.307e+00 | 5.307 | 4.424 | 16.64% |  |
| 12 | Iteration: 6 | 6 | 6 | 0 |  |
| 12 | FOM = 4.687e-01 | 0.4687 | 0.3492 | 25.50% |  |
| 12 | MSE = 5.003e+00 | 5.003 | 4.117 | 17.71% |  |
| 12 | Iteration: 7 | 7 | 7 | 0 |  |
| 12 | FOM = 4.693e-01 | 0.4693 | 0.3505 | 25.31% |  |
| 12 | MSE = 4.747e+00 | 4.747 | 3.85 | 18.90% |  |
| 12 | Iteration: 8 | 8 | 8 | 0 |  |
| 12 | FOM = 4.699e-01 | 0.4699 | 0.3518 | 25.13% |  |
| 12 | MSE = 4.450e+00 | 4.45 | 3.621 | 18.63% |  |
| 12 | Iteration: 9 | 9 | 9 | 0 |  |
| 12 | FOM = 4.705e-01 | 0.4705 | 0.353 | 24.97% |  |
| 12 | MSE = 4.189e+00 | 4.189 | 3.418 | 18.41% |  |
| 12 | Iteration: 10 | 10 | 10 | 0 |  |
| 12 | FOM = 4.709e-01 | 0.4709 | 0.3538 | 24.87% |  |
| 12 | MSE = 3.921e+00 | 3.921 | 3.236 | 17.47% |  |
| 12 | Iteration: 11 | 11 | 11 | 0 |  |
| 12 | FOM = 4.713e-01 | 0.4713 | 0.3545 | 24.78% |  |
| 12 | MSE = 3.720e+00 | 3.72 | 3.063 | 17.66% |  |
| 12 | Iteration: 12 | 12 | 12 | 0 |  |
| 12 | FOM = 4.716e-01 | 0.4716 | 0.3553 | 24.66% |  |
| 12 | MSE = 3.484e+00 | 3.484 | 2.904 | 16.65% |  |
| 12 | Iteration: 13 | 13 | 13 | 0 |  |
| 12 | FOM = 4.719e-01 | 0.4719 | 0.3559 | 24.58% |  |
| 12 | MSE = 3.326e+00 | 3.326 | 2.759 | 17.05% |  |
| 12 | Iteration: 14 | 14 | 14 | 0 |  |
| 12 | FOM = 4.719e-01 | 0.4719 | 0.3564 | 24.48% |  |
| 12 | MSE = 3.133e+00 | 3.133 | 2.627 | 16.15% |  |
| 12 | Iteration: 15 | 15 | 15 | 0 |  |
| 12 | FOM = 4.718e-01 | 0.4718 | 0.3569 | 24.35% |  |
| 12 | MSE = 3.032e+00 | 3.032 | 2.508 | 17.28% |  |
| 12 | Iteration: 16 | 16 | 16 | 0 |  |
| 12 | FOM = 4.716e-01 | 0.4716 | 0.357 | 24.30% |  |
| 12 | MSE = 2.855e+00 | 2.855 | 2.397 | 16.04% |  |
| 12 | Iteration: 17 | 17 | 17 | 0 |  |
| 12 | FOM = 4.708e-01 | 0.4708 | 0.3574 | 24.09% |  |
| 12 | MSE = 2.814e+00 | 2.814 | 2.289 | 18.66% |  |
| 12 | Iteration: 18 | 18 | 18 | 0 |  |
| 12 | FOM = 4.707e-01 | 0.4707 | 0.3575 | 24.05% |  |
| 12 | MSE = 2.617e+00 | 2.617 | 2.179 | 16.74% |  |
| 12 | Iteration: 19 | 19 | 19 | 0 |  |
| 12 | FOM = 4.714e-01 | 0.4714 | 0.3575 | 24.16% |  |
| 12 | MSE = 2.595e+00 | 2.595 | 2.072 | 20.15% |  |
| 12 | Iteration: 20 | 20 | 20 | 0 |  |
| 12 | FOM = 4.709e-01 | 0.4709 | 0.3575 | 24.08% |  |
| 12 | MSE = 2.447e+00 | 2.447 | 1.968 | 19.57% |  |
| 12 | Iteration: 21 | 21 | 21 | 0 |  |
| 12 | FOM = 4.717e-01 | 0.4717 | 0.3575 | 24.21% |  |
| 12 | MSE = 2.421e+00 | 2.421 | 1.867 | 22.88% |  |
| 12 | Iteration: 22 | 22 | 22 | 0 |  |
| 12 | FOM = 4.711e-01 | 0.4711 | 0.3574 | 24.14% |  |
| 12 | MSE = 2.313e+00 | 2.313 | 1.77 | 23.48% |  |
| 12 | Iteration: 23 | 23 | 23 | 0 |  |
| 12 | FOM = 4.720e-01 | 0.472 | 0.3574 | 24.28% |  |
| 12 | MSE = 2.288e+00 | 2.288 | 1.679 | 26.62% |  |
| 12 | Iteration: 24 | 24 | 24 | 0 |  |
| 12 | FOM = 4.714e-01 | 0.4714 | 0.3574 | 24.18% |  |
| 12 | MSE = 2.210e+00 | 2.21 | 1.594 | 27.87% |  |
| 12 | Iteration: 25 | 25 | 25 | 0 |  |
| 12 | FOM = 4.721e-01 | 0.4721 | 0.3575 | 24.27% |  |
| 12 | MSE = 2.199e+00 | 2.199 | 1.514 | 31.15% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 3 | 0 | 4.01% |
| 7 | 0 | 5.43% |
| 10 | 0 | 5.65% |
| 10 | 1 | 5.62% |
| 10 | 2 | 5.89% |
| 10 | 3 | 6.05% |
| 10 | 4 | 6.08% |
| 11 | 0 | 4.36% |
| 12 | 0 | 6.11% |
| 12 | 1 | 6.05% |
| 12 | 2 | 6.03% |
| 12 | 3 | 5.97% |
| 12 | 4 | 5.92% |
| 14 | 0 | 4.27% |
| 15 | 0 | 5.06% |
| 17 | 0 | 7.13% |
| 19 | 0 | 4.75% |
| 21 | 0 | 12.98% |
| 22 | 0 | 4.61% |
