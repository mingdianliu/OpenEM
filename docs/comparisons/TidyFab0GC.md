# TidyFab0GC

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/TidyFab0GC/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Tutorial |
| Verdict | **trajectory divergence** |
| Reference output | reference rerun |
| Cells | 33 on the reference side, 33 on ours, 0 that did not line up |
| Numbers compared | 157 of 184 paired; the other 27 are near-zero, see below |
| Largest relative difference | 1,091% |
| Over 5% / over 20% | 73 / 64 |
| Figures | 9 on the reference side, 9 on ours |
| Largest pixel difference | 12.73% |

Machine-readable form of everything below: [`data/TidyFab0GC.json`](data/TidyFab0GC.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 73 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 2 | WG_CLAD  (1, 0)  Waveguide clad                    #9da6a218     . | 1 | 1 | 0 |  |
| 2 | WG_CLAD  (1, 0)  Waveguide clad                    #9da6a218     . | 0 | 0 | 0 | near-zero |
| 2 | WG_CORE  (2, 0)  Waveguide core                    #6db5dd18     / | 2 | 2 | 0 |  |
| 2 | WG_CORE  (2, 0)  Waveguide core                    #6db5dd18     / | 0 | 0 | 0 | near-zero |
| 2 | SLAB     (3, 0)  Slab region                       #8851ad18     : | 3 | 3 | 0 |  |
| 2 | SLAB     (3, 0)  Slab region                       #8851ad18     : | 0 | 0 | 0 | near-zero |
| 2 | METAL    (5, 0)  Metal layer                       #b8a18b18     \ | 5 | 5 | 0 |  |
| 2 | METAL    (5, 0)  Metal layer                       #b8a18b18     \ | 0 | 0 | 0 | near-zero |
| 3 | 0  ()           -inf, 0         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0 | 0 | 0 | near-zero |
| 3 | 0  ()           -inf, 0         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0 | 0 | 0 | near-zero |
| 3 | 0  ()           -inf, 0         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0 | 0 | 0 | near-zero |
| 3 | 0  ()           -inf, 0         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 12.3 | 12.3 | 0 |  |
| 3 | 1  ()           -2, 1.72        0       SiO2_Palik_LowLoss     Medium(permittivity=4.2) | 1 | 1 | 0 |  |
| 3 | 1  ()           -2, 1.72        0       SiO2_Palik_LowLoss     Medium(permittivity=4.2) | -2 | -2 | 0 |  |
| 3 | 1  ()           -2, 1.72        0       SiO2_Palik_LowLoss     Medium(permittivity=4.2) | 1.72 | 1.72 | 0 |  |
| 3 | 1  ()           -2, 1.72        0       SiO2_Palik_LowLoss     Medium(permittivity=4.2) | 0 | 0 | 0 | near-zero |
| 3 | 1  ()           -2, 1.72        0       SiO2_Palik_LowLoss     Medium(permittivity=4.2) | 4.2 | 4.2 | 0 |  |
| 3 | 2  'WG_CORE'    0, 0.22         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 2 | 2 | 0 |  |
| 3 | 2  'WG_CORE'    0, 0.22         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0 | 0 | 0 | near-zero |
| 3 | 2  'WG_CORE'    0, 0.22         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0.22 | 0.22 | 0 |  |
| 3 | 2  'WG_CORE'    0, 0.22         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0 | 0 | 0 | near-zero |
| 3 | 2  'WG_CORE'    0, 0.22         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 12.3 | 12.3 | 0 |  |
| 3 | 3  'SLAB'       0, 0.15         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 3 | 3 | 0 |  |
| 3 | 3  'SLAB'       0, 0.15         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0 | 0 | 0 | near-zero |
| 3 | 3  'SLAB'       0, 0.15         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0.15 | 0.15 | 0 |  |
| 3 | 3  'SLAB'       0, 0.15         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0 | 0 | 0 | near-zero |
| 3 | 3  'SLAB'       0, 0.15         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 12.3 | 12.3 | 0 |  |
| 3 | 4  'METAL'     1.72, 2.22       0       Cu_JohnsonChristy1972  PEC | 4 | 4 | 0 |  |
| 3 | 4  'METAL'     1.72, 2.22       0       Cu_JohnsonChristy1972  PEC | 1.72 | 1.72 | 0 |  |
| 3 | 4  'METAL'     1.72, 2.22       0       Cu_JohnsonChristy1972  PEC | 2.22 | 2.22 | 0 |  |
| 3 | 4  'METAL'     1.72, 2.22       0       Cu_JohnsonChristy1972  PEC | 0 | 0 | 0 | near-zero |
| 3 | 5  'TRENCH'    -inf, inf        0       Medium()               Medium() | 5 | 5 | 0 |  |
| 3 | 5  'TRENCH'    -inf, inf        0       Medium()               Medium() | 0 | 0 | 0 | near-zero |
| 4 | 0  ()            0, inf         0       air                    air | 0 | 0 | 0 | near-zero |
| 4 | 0  ()            0, inf         0       air                    air | 0 | 0 | 0 | near-zero |
| 4 | 0  ()            0, inf         0       air                    air | 0 | 0 | 0 | near-zero |
| 4 | 1  ()           -inf, 0         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 1 | 1 | 0 |  |
| 4 | 1  ()           -inf, 0         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0 | 0 | 0 | near-zero |
| 4 | 1  ()           -inf, 0         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0 | 0 | 0 | near-zero |
| 4 | 1  ()           -inf, 0         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 12.3 | 12.3 | 0 |  |
| 4 | 2  ()           -2, 1.72        0       SiO2_Palik_LowLoss     Medium(permittivity=4.2) | 2 | 2 | 0 |  |
| 4 | 2  ()           -2, 1.72        0       SiO2_Palik_LowLoss     Medium(permittivity=4.2) | -2 | -2 | 0 |  |
| 4 | 2  ()           -2, 1.72        0       SiO2_Palik_LowLoss     Medium(permittivity=4.2) | 1.72 | 1.72 | 0 |  |
| 4 | 2  ()           -2, 1.72        0       SiO2_Palik_LowLoss     Medium(permittivity=4.2) | 0 | 0 | 0 | near-zero |
| 4 | 2  ()           -2, 1.72        0       SiO2_Palik_LowLoss     Medium(permittivity=4.2) | 4.2 | 4.2 | 0 |  |
| 4 | 3  'WG_CORE'    0, 0.22         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 3 | 3 | 0 |  |
| 4 | 3  'WG_CORE'    0, 0.22         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0 | 0 | 0 | near-zero |
| 4 | 3  'WG_CORE'    0, 0.22         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0.22 | 0.22 | 0 |  |
| 4 | 3  'WG_CORE'    0, 0.22         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0 | 0 | 0 | near-zero |
| 4 | 3  'WG_CORE'    0, 0.22         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 12.3 | 12.3 | 0 |  |
| 4 | 4  'SLAB'       0, 0.15         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 4 | 4 | 0 |  |
| 4 | 4  'SLAB'       0, 0.15         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0 | 0 | 0 | near-zero |
| 4 | 4  'SLAB'       0, 0.15         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0.15 | 0.15 | 0 |  |
| 4 | 4  'SLAB'       0, 0.15         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 0 | 0 | 0 | near-zero |
| 4 | 4  'SLAB'       0, 0.15         0       cSi_Li1993_293K        Medium(permittivity=12.3) | 12.3 | 12.3 | 0 |  |
| 4 | 5  'METAL'     1.72, 2.22       0       Cu_JohnsonChristy1972  PEC | 5 | 5 | 0 |  |
| 4 | 5  'METAL'     1.72, 2.22       0       Cu_JohnsonChristy1972  PEC | 1.72 | 1.72 | 0 |  |
| 4 | 5  'METAL'     1.72, 2.22       0       Cu_JohnsonChristy1972  PEC | 2.22 | 2.22 | 0 |  |
| 4 | 5  'METAL'     1.72, 2.22       0       Cu_JohnsonChristy1972  PEC | 0 | 0 | 0 | near-zero |
| 4 | 6  'TRENCH'    -inf, inf        0       Medium()               Medium() | 6 | 6 | 0 |  |
| 4 | 6  'TRENCH'    -inf, inf        0       Medium()               Medium() | 0 | 0 | 0 | near-zero |
| 14 | (33.565, 0.0, 3.72) | 33.565 | 33.565 | 0 |  |
| 14 | (33.565, 0.0, 3.72) | 0 | 0 | 0 | near-zero |
| 14 | (33.565, 0.0, 3.72) | 3.72 | 3.72 | 0 |  |
| 29 | step = 1 | 1 | 1 | 0 |  |
| 29 | J = 4.6606e-01 | 0.46606 | 0.45851 | 1.62% |  |
| 29 | grad_norm = 4.8987e+00 | 4.8987 | 7.1686 | 46.34% |  |
| 29 | step = 2 | 2 | 2 | 0 |  |
| 29 | J = 3.1737e-01 | 0.31737 | 0.31019 | 2.26% |  |
| 29 | grad_norm = 7.3309e+00 | 7.3309 | 5.9856 | 18.35% |  |
| 29 | step = 3 | 3 | 3 | 0 |  |
| 29 | J = 3.9765e-01 | 0.39765 | 0.30027 | 24.49% |  |
| 29 | grad_norm = 6.3420e+00 | 6.342 | 6.0438 | 4.70% |  |
| 29 | step = 4 | 4 | 4 | 0 |  |
| 29 | J = 5.0306e-01 | 0.50306 | 0.40242 | 20.01% |  |
| 29 | grad_norm = 5.7566e-01 | 0.57566 | 4.0778 | 608% |  |
| 29 | step = 5 | 5 | 5 | 0 |  |
| 29 | J = 4.2190e-01 | 0.4219 | 0.46403 | 9.99% |  |
| 29 | grad_norm = 6.5236e+00 | 6.5236 | 2.332 | 64.25% |  |
| 29 | step = 6 | 6 | 6 | 0 |  |
| 29 | J = 4.0963e-01 | 0.40963 | 0.41464 | 1.22% |  |
| 29 | grad_norm = 6.8633e+00 | 6.8633 | 5.8173 | 15.24% |  |
| 29 | step = 7 | 7 | 7 | 0 |  |
| 29 | J = 4.7728e-01 | 0.47728 | 0.40289 | 15.59% |  |
| 29 | grad_norm = 4.4487e+00 | 4.4487 | 5.4769 | 23.11% |  |
| 29 | step = 8 | 8 | 8 | 0 |  |
| 29 | J = 5.0180e-01 | 0.5018 | 0.41197 | 17.90% |  |
| 29 | grad_norm = 1.7149e+00 | 1.7149 | 2.9707 | 73.23% |  |
| 29 | step = 9 | 9 | 9 | 0 |  |
| 29 | J = 4.5748e-01 | 0.45748 | 0.38534 | 15.77% |  |
| 29 | grad_norm = 5.0710e+00 | 5.071 | 0.98505 | 80.57% |  |
| 29 | step = 10 | 10 | 10 | 0 |  |
| 29 | J = 4.4851e-01 | 0.44851 | 0.3388 | 24.46% |  |
| 29 | grad_norm = 5.3682e+00 | 5.3682 | 2.3768 | 55.72% |  |
| 29 | step = 11 | 11 | 11 | 0 |  |
| 29 | J = 4.8078e-01 | 0.48078 | 0.30794 | 35.95% |  |
| 29 | grad_norm = 3.9810e+00 | 3.981 | 2.8974 | 27.22% |  |
| 29 | step = 12 | 12 | 12 | 0 |  |
| 29 | J = 5.0654e-01 | 0.50654 | 0.30136 | 40.51% |  |
| 29 | grad_norm = 3.8822e-01 | 0.38822 | 2.7991 | 621% |  |
| 29 | step = 13 | 13 | 13 | 0 |  |
| 29 | J = 4.8664e-01 | 0.48664 | 0.31152 | 35.99% |  |
| 29 | grad_norm = 3.5928e+00 | 3.5928 | 1.9741 | 45.05% |  |
| 29 | step = 14 | 14 | 14 | 0 |  |
| 29 | J = 4.7143e-01 | 0.47143 | 0.32249 | 31.59% |  |
| 29 | grad_norm = 4.5558e+00 | 4.5558 | 1.3382 | 70.63% |  |
| 29 | step = 15 | 15 | 15 | 0 |  |
| 29 | J = 4.8606e-01 | 0.48606 | 0.32254 | 33.64% |  |
| 29 | grad_norm = 3.8446e+00 | 3.8446 | 1.7432 | 54.66% |  |
| 29 | step = 16 | 16 | 16 | 0 |  |
| 29 | J = 5.0806e-01 | 0.50806 | 0.31323 | 38.35% |  |
| 29 | grad_norm = 1.1430e+00 | 1.143 | 2.6477 | 132% |  |
| 29 | step = 17 | 17 | 17 | 0 |  |
| 29 | J = 5.0042e-01 | 0.50042 | 0.30882 | 38.29% |  |
| 29 | grad_norm = 2.4919e+00 | 2.4919 | 2.7655 | 10.98% |  |
| 29 | step = 18 | 18 | 18 | 0 |  |
| 29 | J = 4.8600e-01 | 0.486 | 0.30936 | 36.35% |  |
| 29 | grad_norm = 3.6897e+00 | 3.6897 | 2.1909 | 40.62% |  |
| 29 | step = 19 | 19 | 19 | 0 |  |
| 29 | J = 4.9044e-01 | 0.49044 | 0.30548 | 37.71% |  |
| 29 | grad_norm = 3.4153e+00 | 3.4153 | 1.7641 | 48.35% |  |
| 29 | step = 20 | 20 | 20 | 0 |  |
| 29 | J = 5.0664e-01 | 0.50664 | 0.29736 | 41.31% |  |
| 29 | grad_norm = 1.6898e+00 | 1.6898 | 1.6741 | 0.93% |  |
| 29 | step = 21 | 21 | 21 | 0 |  |
| 29 | J = 5.0872e-01 | 0.50872 | 0.286 | 43.78% |  |
| 29 | grad_norm = 1.5054e+00 | 1.5054 | 1.9999 | 32.85% |  |
| 29 | step = 22 | 22 | 22 | 0 |  |
| 29 | J = 4.9747e-01 | 0.49747 | 0.27511 | 44.70% |  |
| 29 | grad_norm = 3.2276e+00 | 3.2276 | 2.1721 | 32.70% |  |
| 29 | step = 23 | 23 | 23 | 0 |  |
| 29 | J = 4.9829e-01 | 0.49829 | 0.27032 | 45.75% |  |
| 29 | grad_norm = 3.1501e+00 | 3.1501 | 2.0376 | 35.32% |  |
| 29 | step = 24 | 24 | 24 | 0 |  |
| 29 | J = 5.0931e-01 | 0.50931 | 0.26979 | 47.03% |  |
| 29 | grad_norm = 1.4531e+00 | 1.4531 | 1.6791 | 15.55% |  |
| 29 | step = 25 | 25 | 25 | 0 |  |
| 29 | J = 5.0980e-01 | 0.5098 | 0.26951 | 47.13% |  |
| 29 | grad_norm = 1.2709e+00 | 1.2709 | 1.1905 | 6.33% |  |
| 29 | step = 26 | 26 | 26 | 0 |  |
| 29 | J = 5.0187e-01 | 0.50187 | 0.26487 | 47.22% |  |
| 29 | grad_norm = 2.6279e+00 | 2.6279 | 1.0449 | 60.24% |  |
| 29 | step = 27 | 27 | 27 | 0 |  |
| 29 | J = 5.0249e-01 | 0.50249 | 0.2561 | 49.03% |  |
| 29 | grad_norm = 2.5851e+00 | 2.5851 | 1.4088 | 45.50% |  |
| 29 | step = 28 | 28 | 28 | 0 |  |
| 29 | J = 5.1051e-01 | 0.51051 | 0.24891 | 51.24% |  |
| 29 | grad_norm = 1.2786e+00 | 1.2786 | 1.6643 | 30.17% |  |
| 29 | step = 29 | 29 | 29 | 0 |  |
| 29 | J = 5.1203e-01 | 0.51203 | 0.24611 | 51.93% |  |
| 29 | grad_norm = 1.0164e+00 | 1.0164 | 1.7184 | 69.07% |  |
| 29 | step = 30 | 30 | 30 | 0 |  |
| 29 | J = 5.0694e-01 | 0.50694 | 0.24632 | 51.41% |  |
| 29 | grad_norm = 2.3533e+00 | 2.3533 | 1.5704 | 33.27% |  |
| 29 | step = 31 | 31 | 31 | 0 |  |
| 29 | J = 5.0837e-01 | 0.50837 | 0.24639 | 51.53% |  |
| 29 | grad_norm = 2.1838e+00 | 2.1838 | 1.5161 | 30.58% |  |
| 29 | step = 32 | 32 | 32 | 0 |  |
| 29 | J = 5.1358e-01 | 0.51358 | 0.24686 | 51.93% |  |
| 29 | grad_norm = 6.5205e-01 | 0.65205 | 1.5609 | 139% |  |
| 29 | step = 33 | 33 | 33 | 0 |  |
| 29 | J = 5.1200e-01 | 0.512 | 0.24828 | 51.51% |  |
| 29 | grad_norm = 1.1880e+00 | 1.188 | 1.7251 | 45.21% |  |
| 29 | step = 34 | 34 | 34 | 0 |  |
| 29 | J = 5.0885e-01 | 0.50885 | 0.25095 | 50.68% |  |
| 29 | grad_norm = 1.9216e+00 | 1.9216 | 1.8424 | 4.12% |  |
| 29 | step = 35 | 35 | 35 | 0 |  |
| 29 | J = 5.1137e-01 | 0.51137 | 0.25575 | 49.99% |  |
| 29 | grad_norm = 1.5056e+00 | 1.5056 | 1.9026 | 26.37% |  |
| 29 | step = 36 | 36 | 36 | 0 |  |
| 29 | J = 5.1490e-01 | 0.5149 | 0.26327 | 48.87% |  |
| 29 | grad_norm = 1.6048e-01 | 0.16048 | 1.912 | 1,091% |  |
| 29 | step = 37 | 37 | 37 | 0 |  |
| 29 | J = 5.1329e-01 | 0.51329 | 0.2726 | 46.89% |  |
| 29 | grad_norm = 1.3641e+00 | 1.3641 | 1.8327 | 34.35% |  |
| 29 | step = 38 | 38 | 38 | 0 |  |
| 29 | J = 5.1245e-01 | 0.51245 | 0.28168 | 45.03% |  |
| 29 | grad_norm = 1.6722e+00 | 1.6722 | 1.6323 | 2.39% |  |
| 29 | step = 39 | 39 | 39 | 0 |  |
| 29 | J = 5.1499e-01 | 0.51499 | 0.28667 | 44.33% |  |
| 29 | grad_norm = 7.8625e-01 | 0.78625 | 1.3368 | 70.02% |  |
| 29 | step = 40 | 40 | 40 | 0 |  |
| 29 | J = 5.1518e-01 | 0.51518 | 0.28426 | 44.82% |  |
| 29 | grad_norm = 6.6040e-01 | 0.6604 | 1.1601 | 75.67% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 12 | 0 | 7.58% |
| 12 | 1 | 8.84% |
| 18 | 0 | 7.07% |
| 19 | 0 | 9.58% |
| 22 | 0 | 8.31% |
| 30 | 0 | 12.73% |
| 30 | 1 | 11.59% |
| 30 | 2 | 10.18% |
| 32 | 0 | 9.41% |
