# Autograd16BilayerCoupler

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd16BilayerCoupler/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **trajectory divergence** |
| Reference output | official archived output |
| Cells | 26 on the reference side, 26 on ours, 0 that did not line up |
| Numbers compared | 96 |
| Largest relative difference | 34.10% |
| Over 5% / over 20% | 29 / 8 |
| Figures | 37 on the reference side, 37 on ours |
| Largest pixel difference | 20.42% |

Machine-readable form of everything below: [`data/Autograd16BilayerCoupler.json`](data/Autograd16BilayerCoupler.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 29 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 8 | Effective index of computed modes:  [[2.70037873 2.18933548 1.81996644]] | 2.70038 | 2.70038 | 4.78e-05% |  |
| 8 | Effective index of computed modes:  [[2.70037873 2.18933548 1.81996644]] | 2.18934 | 2.18936 | 1.13e-03% |  |
| 8 | Effective index of computed modes:  [[2.70037873 2.18933548 1.81996644]] | 1.81997 | 1.81993 | 1.98e-03% |  |
| 18 | objective function value = -0.8573045319173542 | -0.857305 | -0.857352 | 5.50e-03% |  |
| 19 | step = 1 | 1 | 1 | 0 |  |
| 19 | objective = -8.5730e-01 | -0.8573 | -0.85735 | 5.83e-03% |  |
| 19 | grad_norm = 2.7646e-03 | 0.0027646 | 0.0026911 | 2.66% |  |
| 19 | step = 2 | 2 | 2 | 0 |  |
| 19 | objective = -7.0935e-02 | -0.070935 | -0.091272 | 28.67% |  |
| 19 | grad_norm = 1.7276e-02 | 0.017276 | 0.01607 | 6.98% |  |
| 19 | step = 3 | 3 | 3 | 0 |  |
| 19 | objective = 1.9784e-01 | 0.19784 | 0.17383 | 12.14% |  |
| 19 | grad_norm = 2.5687e-02 | 0.025687 | 0.022902 | 10.84% |  |
| 19 | step = 4 | 4 | 4 | 0 |  |
| 19 | objective = 4.1115e-01 | 0.41115 | 0.39455 | 4.04% |  |
| 19 | grad_norm = 2.7394e-02 | 0.027394 | 0.022883 | 16.47% |  |
| 19 | step = 5 | 5 | 5 | 0 |  |
| 19 | objective = 4.7146e-01 | 0.47146 | 0.44794 | 4.99% |  |
| 19 | grad_norm = 2.4444e-02 | 0.024444 | 0.025732 | 5.27% |  |
| 19 | step = 6 | 6 | 6 | 0 |  |
| 19 | objective = 5.4747e-01 | 0.54747 | 0.52778 | 3.60% |  |
| 19 | grad_norm = 1.8427e-02 | 0.018427 | 0.018475 | 0.26% |  |
| 19 | step = 7 | 7 | 7 | 0 |  |
| 19 | objective = 6.0636e-01 | 0.60636 | 0.59421 | 2.00% |  |
| 19 | grad_norm = 1.4281e-02 | 0.014281 | 0.0127 | 11.07% |  |
| 19 | step = 8 | 8 | 8 | 0 |  |
| 19 | objective = 6.3903e-01 | 0.63903 | 0.62843 | 1.66% |  |
| 19 | grad_norm = 1.6747e-02 | 0.016747 | 0.015642 | 6.60% |  |
| 19 | step = 9 | 9 | 9 | 0 |  |
| 19 | objective = 6.8609e-01 | 0.68609 | 0.67592 | 1.48% |  |
| 19 | grad_norm = 1.6021e-02 | 0.016021 | 0.011261 | 29.71% |  |
| 19 | step = 10 | 10 | 10 | 0 |  |
| 19 | objective = 7.2037e-01 | 0.72037 | 0.71183 | 1.19% |  |
| 19 | grad_norm = 1.0072e-02 | 0.010072 | 0.011527 | 14.45% |  |
| 19 | step = 11 | 11 | 11 | 0 |  |
| 19 | objective = 7.3898e-01 | 0.73898 | 0.7466 | 1.03% |  |
| 19 | grad_norm = 1.2367e-02 | 0.012367 | 0.010136 | 18.04% |  |
| 19 | step = 12 | 12 | 12 | 0 |  |
| 19 | objective = 7.6888e-01 | 0.76888 | 0.7734 | 0.59% |  |
| 19 | grad_norm = 7.8178e-03 | 0.0078178 | 0.0093459 | 19.55% |  |
| 19 | step = 13 | 13 | 13 | 0 |  |
| 19 | objective = 7.9041e-01 | 0.79041 | 0.79434 | 0.50% |  |
| 19 | grad_norm = 9.0470e-03 | 0.009047 | 0.010907 | 20.56% |  |
| 19 | step = 14 | 14 | 14 | 0 |  |
| 19 | objective = 8.0788e-01 | 0.80788 | 0.81749 | 1.19% |  |
| 19 | grad_norm = 8.4449e-03 | 0.0084449 | 0.0067576 | 19.98% |  |
| 19 | step = 15 | 15 | 15 | 0 |  |
| 19 | objective = 8.2130e-01 | 0.8213 | 0.8334 | 1.47% |  |
| 19 | grad_norm = 7.3545e-03 | 0.0073545 | 0.0082061 | 11.58% |  |
| 19 | step = 16 | 16 | 16 | 0 |  |
| 19 | objective = 8.3541e-01 | 0.83541 | 0.85145 | 1.92% |  |
| 19 | grad_norm = 8.2459e-03 | 0.0082459 | 0.0057766 | 29.95% |  |
| 19 | step = 17 | 17 | 17 | 0 |  |
| 19 | objective = 8.5142e-01 | 0.85142 | 0.8648 | 1.57% |  |
| 19 | grad_norm = 6.3082e-03 | 0.0063082 | 0.0059198 | 6.16% |  |
| 19 | step = 18 | 18 | 18 | 0 |  |
| 19 | objective = 8.6365e-01 | 0.86365 | 0.87817 | 1.68% |  |
| 19 | grad_norm = 5.7417e-03 | 0.0057417 | 0.0059407 | 3.47% |  |
| 19 | step = 19 | 19 | 19 | 0 |  |
| 19 | objective = 8.7403e-01 | 0.87403 | 0.89194 | 2.05% |  |
| 19 | grad_norm = 6.1316e-03 | 0.0061316 | 0.0041604 | 32.15% |  |
| 19 | step = 20 | 20 | 20 | 0 |  |
| 19 | objective = 8.8415e-01 | 0.88415 | 0.90101 | 1.91% |  |
| 19 | grad_norm = 4.9958e-03 | 0.0049958 | 0.0049727 | 0.46% |  |
| 19 | step = 21 | 21 | 21 | 0 |  |
| 19 | objective = 8.9300e-01 | 0.893 | 0.90938 | 1.83% |  |
| 19 | grad_norm = 4.6502e-03 | 0.0046502 | 0.0041589 | 10.57% |  |
| 19 | step = 22 | 22 | 22 | 0 |  |
| 19 | objective = 9.0172e-01 | 0.90172 | 0.91702 | 1.70% |  |
| 19 | grad_norm = 4.4607e-03 | 0.0044607 | 0.0047404 | 6.27% |  |
| 19 | step = 23 | 23 | 23 | 0 |  |
| 19 | objective = 9.0994e-01 | 0.90994 | 0.92432 | 1.58% |  |
| 19 | grad_norm = 3.9034e-03 | 0.0039034 | 0.0044219 | 13.28% |  |
| 19 | step = 24 | 24 | 24 | 0 |  |
| 19 | objective = 9.1736e-01 | 0.91736 | 0.93161 | 1.55% |  |
| 19 | grad_norm = 4.2583e-03 | 0.0042583 | 0.003835 | 9.94% |  |
| 19 | step = 25 | 25 | 25 | 0 |  |
| 19 | objective = 9.2464e-01 | 0.92464 | 0.93873 | 1.52% |  |
| 19 | grad_norm = 3.9619e-03 | 0.0039619 | 0.0034788 | 12.19% |  |
| 19 | step = 26 | 26 | 26 | 0 |  |
| 19 | objective = 9.3161e-01 | 0.93161 | 0.94581 | 1.52% |  |
| 19 | grad_norm = 4.1826e-03 | 0.0041826 | 0.0031499 | 24.69% |  |
| 19 | step = 27 | 27 | 27 | 0 |  |
| 19 | objective = 9.3748e-01 | 0.93748 | 0.95261 | 1.61% |  |
| 19 | grad_norm = 4.7450e-03 | 0.004745 | 0.0031269 | 34.10% |  |
| 19 | step = 28 | 28 | 28 | 0 |  |
| 19 | objective = 9.4364e-01 | 0.94364 | 0.95915 | 1.64% |  |
| 19 | grad_norm = 3.6218e-03 | 0.0036218 | 0.0031081 | 14.18% |  |
| 19 | step = 29 | 29 | 29 | 0 |  |
| 19 | objective = 9.4776e-01 | 0.94776 | 0.96497 | 1.82% |  |
| 19 | grad_norm = 3.9977e-03 | 0.0039977 | 0.0029993 | 24.97% |  |
| 19 | step = 30 | 30 | 30 | 0 |  |
| 19 | objective = 9.5310e-01 | 0.9531 | 0.96927 | 1.70% |  |
| 19 | grad_norm = 4.0515e-03 | 0.0040515 | 0.0032817 | 19.00% |  |
| 24 | final penalty (Si) = 0.110 | 0.11 | 0.101 | 8.18% |  |
| 24 | final penalty (SiN) = 0.022 | 0.022 | 0.022 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 6 | 0 | 5.58% |
| 8 | 0 | 15.15% |
| 11 | 0 | 7.59% |
| 13 | 0 | 6.81% |
| 18 | 0 | 3.86% |
| 19 | 0 | 3.08% |
| 19 | 1 | 7.44% |
| 19 | 2 | 9.26% |
| 19 | 3 | 10.82% |
| 19 | 4 | 11.88% |
| 19 | 5 | 12.45% |
| 19 | 6 | 12.82% |
| 19 | 7 | 13.14% |
| 19 | 8 | 13.45% |
| 19 | 9 | 13.72% |
| 19 | 10 | 14.16% |
| 19 | 11 | 14.47% |
| 19 | 12 | 14.79% |
| 19 | 13 | 15.16% |
| 19 | 14 | 15.41% |
| 19 | 15 | 15.77% |
| 19 | 16 | 16.11% |
| 19 | 17 | 16.53% |
| 19 | 18 | 17.10% |
| 19 | 19 | 17.57% |
| 19 | 20 | 17.89% |
| 19 | 21 | 18.19% |
| 19 | 22 | 18.50% |
| 19 | 23 | 18.87% |
| 19 | 24 | 19.14% |
| 19 | 25 | 19.50% |
| 19 | 26 | 19.65% |
| 19 | 27 | 19.95% |
| 19 | 28 | 20.21% |
| 19 | 29 | 20.42% |
| 21 | 0 | 3.71% |
| 23 | 0 | 14.30% |
