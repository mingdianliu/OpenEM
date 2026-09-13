# Autograd25WaveguideCrossing

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd25WaveguideCrossing/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **close** |
| Reference output | official archived output |
| Cells | 21 on the reference side, 21 on ours, 1 that did not line up |
| Numbers compared | 79 |
| Largest relative difference | 77.05% |
| Over 5% / over 20% | 25 / 19 |
| Figures | 33 on the reference side, 33 on ours |
| Largest pixel difference | 14.59% |

Machine-readable form of everything below: [`data/Autograd25WaveguideCrossing.json`](data/Autograd25WaveguideCrossing.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 25 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 1 | floor_radius=15 | 15 | 15 | 0 |  |
| 1 | discrete_radius=15 | 15 | 15 | 0 |  |
| 1 | discrete_num_pixels=300 | 300 | 300 | 0 |  |
| 11 | Initial transmission at center wavelength: 0.7747789256088472 | 0.774779 | 0.774486 | 0.04% |  |
| 14 | step = 1 | 1 | 1 | 0 |  |
| 14 | J = 7.7478e-01 | 0.77478 | 0.77449 | 0.04% |  |
| 14 | grad_norm = 2.3355e-03 | 0.0023355 | 0.0027353 | 17.12% |  |
| 14 | step = 2 | 2 | 2 | 0 |  |
| 14 | J = 7.9722e-01 | 0.79722 | 0.79446 | 0.35% |  |
| 14 | grad_norm = 5.1836e-03 | 0.0051836 | 0.0048874 | 5.71% |  |
| 14 | step = 3 | 3 | 3 | 0 |  |
| 14 | J = 8.1713e-01 | 0.81713 | 0.81161 | 0.68% |  |
| 14 | grad_norm = 4.7390e-03 | 0.004739 | 0.003608 | 23.87% |  |
| 14 | step = 4 | 4 | 4 | 0 |  |
| 14 | J = 8.3984e-01 | 0.83984 | 0.82488 | 1.78% |  |
| 14 | grad_norm = 6.0274e-03 | 0.0060274 | 0.0050515 | 16.19% |  |
| 14 | step = 5 | 5 | 5 | 0 |  |
| 14 | J = 8.6152e-01 | 0.86152 | 0.84302 | 2.15% |  |
| 14 | grad_norm = 6.1428e-03 | 0.0061428 | 0.0036478 | 40.62% |  |
| 14 | step = 6 | 6 | 6 | 0 |  |
| 14 | J = 8.8062e-01 | 0.88062 | 0.86052 | 2.28% |  |
| 14 | grad_norm = 4.7030e-03 | 0.004703 | 0.0050133 | 6.60% |  |
| 14 | step = 7 | 7 | 7 | 0 |  |
| 14 | J = 8.7419e-01 | 0.87419 | 0.8628 | 1.30% |  |
| 14 | grad_norm = 3.4181e-02 | 0.034181 | 0.014076 | 58.82% |  |
| 14 | step = 8 | 8 | 8 | 0 |  |
| 14 | J = 8.9761e-01 | 0.89761 | 0.88238 | 1.70% |  |
| 14 | grad_norm = 3.7202e-03 | 0.0037202 | 0.0047497 | 27.67% |  |
| 14 | step = 9 | 9 | 9 | 0 |  |
| 14 | J = 8.9889e-01 | 0.89889 | 0.88888 | 1.11% |  |
| 14 | grad_norm = 4.7613e-03 | 0.0047613 | 0.0060517 | 27.10% |  |
| 14 | step = 10 | 10 | 10 | 0 |  |
| 14 | J = 9.0241e-01 | 0.90241 | 0.89761 | 0.53% |  |
| 14 | grad_norm = 4.4704e-03 | 0.0044704 | 0.006452 | 44.33% |  |
| 14 | step = 11 | 11 | 11 | 0 |  |
| 14 | J = 9.0740e-01 | 0.9074 | 0.90946 | 0.23% |  |
| 14 | grad_norm = 4.6739e-03 | 0.0046739 | 0.0076002 | 62.61% |  |
| 14 | step = 12 | 12 | 12 | 0 |  |
| 14 | J = 9.1445e-01 | 0.91445 | 0.92115 | 0.73% |  |
| 14 | grad_norm = 3.9454e-03 | 0.0039454 | 0.0027152 | 31.18% |  |
| 14 | step = 13 | 13 | 13 | 0 |  |
| 14 | J = 9.2114e-01 | 0.92114 | 0.92308 | 0.21% |  |
| 14 | grad_norm = 3.6759e-03 | 0.0036759 | 0.0032483 | 11.63% |  |
| 14 | step = 14 | 14 | 14 | 0 |  |
| 14 | J = 9.2644e-01 | 0.92644 | 0.92632 | 0.01% |  |
| 14 | grad_norm = 5.5271e-03 | 0.0055271 | 0.0026821 | 51.47% |  |
| 14 | step = 15 | 15 | 15 | 0 |  |
| 14 | J = 9.2947e-01 | 0.92947 | 0.93051 | 0.11% |  |
| 14 | grad_norm = 3.0017e-03 | 0.0030017 | 0.0021555 | 28.19% |  |
| 14 | step = 16 | 16 | 16 | 0 |  |
| 14 | J = 9.3158e-01 | 0.93158 | 0.93405 | 0.27% |  |
| 14 | grad_norm = 3.3266e-03 | 0.0033266 | 0.0018811 | 43.45% |  |
| 14 | step = 17 | 17 | 17 | 0 |  |
| 14 | J = 9.3527e-01 | 0.93527 | 0.93776 | 0.27% |  |
| 14 | grad_norm = 2.3057e-03 | 0.0023057 | 0.0019192 | 16.76% |  |
| 14 | step = 18 | 18 | 18 | 0 |  |
| 14 | J = 9.3906e-01 | 0.93906 | 0.94152 | 0.26% |  |
| 14 | grad_norm = 1.7436e-03 | 0.0017436 | 0.0013237 | 24.08% |  |
| 14 | step = 19 | 19 | 19 | 0 |  |
| 14 | J = 9.3952e-01 | 0.93952 | 0.94445 | 0.52% |  |
| 14 | grad_norm = 5.8448e-03 | 0.0058448 | 0.0013414 | 77.05% |  |
| 14 | step = 20 | 20 | 20 | 0 |  |
| 14 | J = 9.3894e-01 | 0.93894 | 0.94713 | 0.87% |  |
| 14 | grad_norm = 7.6984e-03 | 0.0076984 | 0.0018497 | 75.97% |  |
| 14 | step = 21 | 21 | 21 | 0 |  |
| 14 | J = 9.4332e-01 | 0.94332 | 0.95004 | 0.71% |  |
| 14 | grad_norm = 3.7946e-03 | 0.0037946 | 0.0012236 | 67.75% |  |
| 14 | step = 22 | 22 | 22 | 0 |  |
| 14 | J = 9.4567e-01 | 0.94567 | 0.95157 | 0.62% |  |
| 14 | grad_norm = 2.9461e-03 | 0.0029461 | 0.0015301 | 48.06% |  |
| 14 | step = 23 | 23 | 23 | 0 |  |
| 14 | J = 9.4783e-01 | 0.94783 | 0.95145 | 0.38% |  |
| 14 | grad_norm = 2.7046e-03 | 0.0027046 | 0.0018827 | 30.39% |  |
| 14 | step = 24 | 24 | 24 | 0 |  |
| 14 | J = 9.4879e-01 | 0.94879 | 0.95303 | 0.45% |  |
| 14 | grad_norm = 2.6637e-03 | 0.0026637 | 0.0016564 | 37.82% |  |
| 14 | step = 25 | 25 | 25 | 0 |  |
| 14 | J = 9.4947e-01 | 0.94947 | 0.95417 | 0.50% |  |
| 14 | grad_norm = 2.5249e-03 | 0.0025249 | 0.0017903 | 29.09% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 3 | 0 | 7.42% |
| 7 | 0 | 14.59% |
| 9 | 0 | 6.54% |
| 11 | 0 | 5.71% |
| 14 | 0 | 9.28% |
| 14 | 1 | 9.23% |
| 14 | 2 | 9.29% |
| 14 | 3 | 9.36% |
| 14 | 4 | 9.44% |
| 14 | 5 | 9.49% |
| 14 | 6 | 9.59% |
| 14 | 7 | 9.57% |
| 14 | 8 | 9.59% |
| 14 | 9 | 9.59% |
| 14 | 10 | 9.60% |
| 14 | 11 | 9.61% |
| 14 | 12 | 9.64% |
| 14 | 13 | 9.66% |
| 14 | 14 | 9.69% |
| 14 | 15 | 9.66% |
| 14 | 16 | 9.67% |
| 14 | 17 | 9.67% |
| 14 | 18 | 9.70% |
| 14 | 19 | 9.71% |
| 14 | 20 | 9.75% |
| 14 | 21 | 9.78% |
| 14 | 22 | 9.82% |
| 14 | 23 | 9.86% |
| 14 | 24 | 9.87% |
| 15 | 0 | 4.10% |
| 16 | 0 | 7.66% |
| 18 | 0 | 7.13% |
| 19 | 0 | 13.22% |
