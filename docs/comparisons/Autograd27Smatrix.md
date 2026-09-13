# Autograd27Smatrix

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd27Smatrix/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **trajectory divergence** |
| Reference output | official archived output |
| Cells | 26 on the reference side, 26 on ours, 1 that did not line up |
| Numbers compared | 108 |
| Largest relative difference | 1,825% |
| Over 5% / over 20% | 47 / 29 |
| Figures | 33 on the reference side, 33 on ours |
| Largest pixel difference | 9.59% |

Machine-readable form of everything below: [`data/Autograd27Smatrix.json`](data/Autograd27Smatrix.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 47 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 4 | 451 451 | 451 | 451 | 0 |  |
| 4 | 451 451 | 451 | 451 | 0 |  |
| 17 | J = 0.2660 | 0.266 | 0.256 | 3.76% |  |
| 17 | gradient shape = (451, 451) | 451 | 451 | 0 |  |
| 17 | gradient shape = (451, 451) | 451 | 451 | 0 |  |
| 17 | gradient norm = 0.0286 | 0.0286 | 0.0287 | 0.35% |  |
| 18 | step = 1 | 1 | 1 | 0 |  |
| 18 | J = 2.8214e-01 | 0.28214 | 0.28363 | 0.53% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 1.8372e-02 | 0.018372 | 0.01952 | 6.25% |  |
| 18 | step = 2 | 2 | 2 | 0 |  |
| 18 | J = 2.3074e-01 | 0.23074 | 0.21601 | 6.38% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 9.9531e-03 | 0.0099531 | 0.0079353 | 20.27% |  |
| 18 | step = 3 | 3 | 3 | 0 |  |
| 18 | J = 2.0380e-01 | 0.2038 | 0.19003 | 6.76% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 6.5286e-03 | 0.0065286 | 0.0064624 | 1.01% |  |
| 18 | step = 4 | 4 | 4 | 0 |  |
| 18 | J = 1.8484e-01 | 0.18484 | 0.17531 | 5.16% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 5.9481e-03 | 0.0059481 | 0.0063969 | 7.55% |  |
| 18 | step = 5 | 5 | 5 | 0 |  |
| 18 | J = 1.7224e-01 | 0.17224 | 0.16073 | 6.68% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 5.9589e-03 | 0.0059589 | 0.0071999 | 20.83% |  |
| 18 | step = 6 | 6 | 6 | 0 |  |
| 18 | J = 1.5902e-01 | 0.15902 | 0.14566 | 8.40% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 6.5046e-03 | 0.0065046 | 0.0084379 | 29.72% |  |
| 18 | step = 7 | 7 | 7 | 0 |  |
| 18 | J = 1.4726e-01 | 0.14726 | 0.13026 | 11.54% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 7.5650e-03 | 0.007565 | 0.0099171 | 31.09% |  |
| 18 | step = 8 | 8 | 8 | 0 |  |
| 18 | J = 1.3329e-01 | 0.13329 | 0.11101 | 16.72% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 8.9116e-03 | 0.0089116 | 0.011 | 23.43% |  |
| 18 | step = 9 | 9 | 9 | 0 |  |
| 18 | J = 1.1570e-01 | 0.1157 | 0.090895 | 21.44% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 1.0077e-02 | 0.010077 | 0.010938 | 8.54% |  |
| 18 | step = 10 | 10 | 10 | 0 |  |
| 18 | J = 9.7908e-02 | 0.097908 | 0.071777 | 26.69% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 1.0771e-02 | 0.010771 | 0.0096569 | 10.34% |  |
| 18 | step = 11 | 11 | 11 | 0 |  |
| 18 | J = 7.7544e-02 | 0.077544 | 0.060675 | 21.75% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 1.0175e-02 | 0.010175 | 0.0085794 | 15.68% |  |
| 18 | step = 12 | 12 | 12 | 0 |  |
| 18 | J = 6.2473e-02 | 0.062473 | 0.051217 | 18.02% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 8.5627e-03 | 0.0085627 | 0.0091142 | 6.44% |  |
| 18 | step = 13 | 13 | 13 | 0 |  |
| 18 | J = 5.1949e-02 | 0.051949 | 0.038479 | 25.93% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 9.1905e-03 | 0.0091905 | 0.0058159 | 36.72% |  |
| 18 | step = 14 | 14 | 14 | 0 |  |
| 18 | J = 4.5226e-02 | 0.045226 | 0.038466 | 14.95% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 9.9427e-03 | 0.0099427 | 0.0058859 | 40.80% |  |
| 18 | step = 15 | 15 | 15 | 0 |  |
| 18 | J = 3.7023e-02 | 0.037023 | 0.043785 | 18.26% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 7.7900e-03 | 0.00779 | 0.01706 | 119% |  |
| 18 | step = 16 | 16 | 16 | 0 |  |
| 18 | J = 4.4086e-02 | 0.044086 | 0.036945 | 16.20% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 2.0016e-02 | 0.020016 | 0.01089 | 45.59% |  |
| 18 | step = 17 | 17 | 17 | 0 |  |
| 18 | J = 3.3325e-02 | 0.033325 | 0.053511 | 60.57% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 1.0578e-02 | 0.010578 | 0.041714 | 294% |  |
| 18 | step = 18 | 18 | 18 | 0 |  |
| 18 | J = 5.4691e-02 | 0.054691 | 0.028671 | 47.58% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 7.0848e-02 | 0.070848 | 0.010042 | 85.83% |  |
| 18 | step = 19 | 19 | 19 | 0 |  |
| 18 | J = 4.0703e-02 | 0.040703 | 0.043535 | 6.96% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 4.5776e-02 | 0.045776 | 0.04367 | 4.60% |  |
| 18 | step = 20 | 20 | 20 | 0 |  |
| 18 | J = 3.9999e-02 | 0.039999 | 0.018192 | 54.52% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 3.8556e-02 | 0.038556 | 0.0034499 | 91.05% |  |
| 18 | step = 21 | 21 | 21 | 0 |  |
| 18 | J = 1.9562e-02 | 0.019562 | 0.045633 | 133% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 3.4968e-03 | 0.0034968 | 0.067303 | 1,825% |  |
| 18 | step = 22 | 22 | 22 | 0 |  |
| 18 | J = 3.2912e-02 | 0.032912 | 0.016838 | 48.84% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 3.0316e-02 | 0.030316 | 0.0076124 | 74.89% |  |
| 18 | step = 23 | 23 | 23 | 0 |  |
| 18 | J = 2.6889e-02 | 0.026889 | 0.03645 | 35.56% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 2.1668e-02 | 0.021668 | 0.034029 | 57.05% |  |
| 18 | step = 24 | 24 | 24 | 0 |  |
| 18 | J = 2.1815e-02 | 0.021815 | 0.029924 | 37.17% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 3.9575e-03 | 0.0039575 | 0.023658 | 498% |  |
| 18 | step = 25 | 25 | 25 | 0 |  |
| 18 | J = 2.6270e-02 | 0.02627 | 0.018544 | 29.41% |  |
| 18 | beta = 30.00 | 30 | 30 | 0 |  |
| 18 | grad_norm = 1.5076e-02 | 0.015076 | 0.0039241 | 73.97% |  |
| 21 | Text(0.5, 1.0, 'final design') | 0.5 | 0.5 | 0 |  |
| 21 | Text(0.5, 1.0, 'final design') | 1 | 1 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 3 | 0 | 4.02% |
| 7 | 0 | 4.19% |
| 10 | 0 | 4.31% |
| 12 | 0 | 4.21% |
| 18 | 0 | 0 |
| 18 | 1 | 0 |
| 18 | 2 | 0 |
| 18 | 3 | 0.11% |
| 18 | 4 | 0.22% |
| 18 | 5 | 0.40% |
| 18 | 6 | 0.55% |
| 18 | 7 | 0.62% |
| 18 | 8 | 0.58% |
| 18 | 9 | 0.67% |
| 18 | 10 | 0.69% |
| 18 | 11 | 0.67% |
| 18 | 12 | 0.62% |
| 18 | 13 | 0.60% |
| 18 | 14 | 0.54% |
| 18 | 15 | 0.46% |
| 18 | 16 | 0.44% |
| 18 | 17 | 0.37% |
| 18 | 18 | 0.18% |
| 18 | 19 | 0.33% |
| 18 | 20 | 0.57% |
| 18 | 21 | 0.32% |
| 18 | 22 | 0.42% |
| 18 | 23 | 0.30% |
| 18 | 24 | 0.28% |
| 19 | 0 | 4.12% |
| 21 | 0 | 3.95% |
| 22 | 0 | 4.29% |
| 24 | 0 | 9.59% |
