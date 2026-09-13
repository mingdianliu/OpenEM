# Autograd31GratingCouplerWithBeamOptimization

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd31GratingCouplerWithBeamOptimization/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **trajectory divergence** |
| Reference output | official archived output |
| Cells | 15 on the reference side, 15 on ours, 0 that did not line up |
| Numbers compared | 126 of 177 paired; the other 51 are near-zero, see below |
| Largest relative difference | 138% |
| Over 5% / over 20% | 23 / 17 |
| Figures | 7 on the reference side, 7 on ours |
| Largest pixel difference | 8.92% |

Machine-readable form of everything below: [`data/Autograd31GratingCouplerWithBeamOptimization.json`](data/Autograd31GratingCouplerWithBeamOptimization.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 23 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 6 | Initial x_center = 6.000 um | 6 | 6 | 0 |  |
| 6 | Initial tilt = 8.000 deg | 8 | 8 | 0 |  |
| 6 | Initial waist = 2.250 um | 2.25 | 2.25 | 0 |  |
| 7 | epochs = 90, learning_rate = 0.0030 | 90 | 90 | 0 |  |
| 7 | epochs = 90, learning_rate = 0.0030 | 0.003 | 0.003 | 0 | near-zero |
| 7 | step 01/90 \| objective = 0.084270 \| grad norm = 1.229e+00 | 1 | 1 | 0 |  |
| 7 | step 01/90 \| objective = 0.084270 \| grad norm = 1.229e+00 | 90 | 90 | 0 |  |
| 7 | step 01/90 \| objective = 0.084270 \| grad norm = 1.229e+00 | 0.08427 | 0.08384 | 0.05% | near-zero |
| 7 | step 01/90 \| objective = 0.084270 \| grad norm = 1.229e+00 | 1.229 | 1.231 | 0.16% |  |
| 7 | step 05/90 \| objective = 0.140750 \| grad norm = 6.235e-01 | 5 | 5 | 0 |  |
| 7 | step 05/90 \| objective = 0.140750 \| grad norm = 6.235e-01 | 90 | 90 | 0 |  |
| 7 | step 05/90 \| objective = 0.140750 \| grad norm = 6.235e-01 | 0.14075 | 0.13932 | 0.16% | near-zero |
| 7 | step 05/90 \| objective = 0.140750 \| grad norm = 6.235e-01 | 0.6235 | 0.6533 | 3.31% | near-zero |
| 7 | step 10/90 \| objective = 0.145501 \| grad norm = 1.005e+00 | 10 | 10 | 0 |  |
| 7 | step 10/90 \| objective = 0.145501 \| grad norm = 1.005e+00 | 90 | 90 | 0 |  |
| 7 | step 10/90 \| objective = 0.145501 \| grad norm = 1.005e+00 | 0.145501 | 0.143031 | 0.27% | near-zero |
| 7 | step 10/90 \| objective = 0.145501 \| grad norm = 1.005e+00 | 1.005 | 1.131 | 12.54% |  |
| 7 | step 15/90 \| objective = 0.165367 \| grad norm = 8.510e-01 | 15 | 15 | 0 |  |
| 7 | step 15/90 \| objective = 0.165367 \| grad norm = 8.510e-01 | 90 | 90 | 0 |  |
| 7 | step 15/90 \| objective = 0.165367 \| grad norm = 8.510e-01 | 0.165367 | 0.163848 | 0.17% | near-zero |
| 7 | step 15/90 \| objective = 0.165367 \| grad norm = 8.510e-01 | 0.851 | 0.8593 | 0.92% | near-zero |
| 7 | step 20/90 \| objective = 0.176042 \| grad norm = 1.053e+00 | 20 | 20 | 0 |  |
| 7 | step 20/90 \| objective = 0.176042 \| grad norm = 1.053e+00 | 90 | 90 | 0 |  |
| 7 | step 20/90 \| objective = 0.176042 \| grad norm = 1.053e+00 | 0.176042 | 0.177889 | 0.21% | near-zero |
| 7 | step 20/90 \| objective = 0.176042 \| grad norm = 1.053e+00 | 1.053 | 0.8555 | 18.76% |  |
| 7 | step 25/90 \| objective = 0.199274 \| grad norm = 6.669e-01 | 25 | 25 | 0 |  |
| 7 | step 25/90 \| objective = 0.199274 \| grad norm = 6.669e-01 | 90 | 90 | 0 |  |
| 7 | step 25/90 \| objective = 0.199274 \| grad norm = 6.669e-01 | 0.199274 | 0.199794 | 0.06% | near-zero |
| 7 | step 25/90 \| objective = 0.199274 \| grad norm = 6.669e-01 | 0.6669 | 0.7418 | 8.32% | near-zero |
| 7 | step 30/90 \| objective = 0.214997 \| grad norm = 7.155e-01 | 30 | 30 | 0 |  |
| 7 | step 30/90 \| objective = 0.214997 \| grad norm = 7.155e-01 | 90 | 90 | 0 |  |
| 7 | step 30/90 \| objective = 0.214997 \| grad norm = 7.155e-01 | 0.214997 | 0.216544 | 0.17% | near-zero |
| 7 | step 30/90 \| objective = 0.214997 \| grad norm = 7.155e-01 | 0.7155 | 0.8516 | 15.12% | near-zero |
| 7 | step 35/90 \| objective = 0.234634 \| grad norm = 8.797e-01 | 35 | 35 | 0 |  |
| 7 | step 35/90 \| objective = 0.234634 \| grad norm = 8.797e-01 | 90 | 90 | 0 |  |
| 7 | step 35/90 \| objective = 0.234634 \| grad norm = 8.797e-01 | 0.234634 | 0.23384 | 0.09% | near-zero |
| 7 | step 35/90 \| objective = 0.234634 \| grad norm = 8.797e-01 | 0.8797 | 1.016 | 15.14% | near-zero |
| 7 | step 40/90 \| objective = 0.253888 \| grad norm = 7.545e-01 | 40 | 40 | 0 |  |
| 7 | step 40/90 \| objective = 0.253888 \| grad norm = 7.545e-01 | 90 | 90 | 0 |  |
| 7 | step 40/90 \| objective = 0.253888 \| grad norm = 7.545e-01 | 0.253888 | 0.25159 | 0.26% | near-zero |
| 7 | step 40/90 \| objective = 0.253888 \| grad norm = 7.545e-01 | 0.7545 | 0.7849 | 3.38% | near-zero |
| 7 | step 45/90 \| objective = 0.274034 \| grad norm = 9.434e-01 | 45 | 45 | 0 |  |
| 7 | step 45/90 \| objective = 0.274034 \| grad norm = 9.434e-01 | 90 | 90 | 0 |  |
| 7 | step 45/90 \| objective = 0.274034 \| grad norm = 9.434e-01 | 0.274034 | 0.269204 | 0.54% | near-zero |
| 7 | step 45/90 \| objective = 0.274034 \| grad norm = 9.434e-01 | 0.9434 | 1.005 | 6.53% |  |
| 7 | step 50/90 \| objective = 0.288923 \| grad norm = 1.119e+00 | 50 | 50 | 0 |  |
| 7 | step 50/90 \| objective = 0.288923 \| grad norm = 1.119e+00 | 90 | 90 | 0 |  |
| 7 | step 50/90 \| objective = 0.288923 \| grad norm = 1.119e+00 | 0.288923 | 0.286909 | 0.22% | near-zero |
| 7 | step 50/90 \| objective = 0.288923 \| grad norm = 1.119e+00 | 1.119 | 0.7145 | 36.15% |  |
| 7 | step 55/90 \| objective = 0.307696 \| grad norm = 8.266e-01 | 55 | 55 | 0 |  |
| 7 | step 55/90 \| objective = 0.307696 \| grad norm = 8.266e-01 | 90 | 90 | 0 |  |
| 7 | step 55/90 \| objective = 0.307696 \| grad norm = 8.266e-01 | 0.307696 | 0.304882 | 0.31% | near-zero |
| 7 | step 55/90 \| objective = 0.307696 \| grad norm = 8.266e-01 | 0.8266 | 0.7304 | 10.69% | near-zero |
| 7 | step 60/90 \| objective = 0.317485 \| grad norm = 1.209e+00 | 60 | 60 | 0 |  |
| 7 | step 60/90 \| objective = 0.317485 \| grad norm = 1.209e+00 | 90 | 90 | 0 |  |
| 7 | step 60/90 \| objective = 0.317485 \| grad norm = 1.209e+00 | 0.317485 | 0.316668 | 0.09% | near-zero |
| 7 | step 60/90 \| objective = 0.317485 \| grad norm = 1.209e+00 | 1.209 | 0.843 | 30.27% |  |
| 7 | step 65/90 \| objective = 0.329971 \| grad norm = 4.592e-01 | 65 | 65 | 0 |  |
| 7 | step 65/90 \| objective = 0.329971 \| grad norm = 4.592e-01 | 90 | 90 | 0 |  |
| 7 | step 65/90 \| objective = 0.329971 \| grad norm = 4.592e-01 | 0.329971 | 0.328109 | 0.21% | near-zero |
| 7 | step 65/90 \| objective = 0.329971 \| grad norm = 4.592e-01 | 0.4592 | 0.7169 | 28.63% | near-zero |
| 7 | step 70/90 \| objective = 0.334242 \| grad norm = 9.036e-01 | 70 | 70 | 0 |  |
| 7 | step 70/90 \| objective = 0.334242 \| grad norm = 9.036e-01 | 90 | 90 | 0 |  |
| 7 | step 70/90 \| objective = 0.334242 \| grad norm = 9.036e-01 | 0.334242 | 0.335319 | 0.12% | near-zero |
| 7 | step 70/90 \| objective = 0.334242 \| grad norm = 9.036e-01 | 0.9036 | 1.196 | 32.36% |  |
| 7 | step 75/90 \| objective = 0.342667 \| grad norm = 7.190e-01 | 75 | 75 | 0 |  |
| 7 | step 75/90 \| objective = 0.342667 \| grad norm = 7.190e-01 | 90 | 90 | 0 |  |
| 7 | step 75/90 \| objective = 0.342667 \| grad norm = 7.190e-01 | 0.342667 | 0.342074 | 0.07% | near-zero |
| 7 | step 75/90 \| objective = 0.342667 \| grad norm = 7.190e-01 | 0.719 | 0.4481 | 30.10% | near-zero |
| 7 | step 80/90 \| objective = 0.344083 \| grad norm = 9.081e-01 | 80 | 80 | 0 |  |
| 7 | step 80/90 \| objective = 0.344083 \| grad norm = 9.081e-01 | 90 | 90 | 0 |  |
| 7 | step 80/90 \| objective = 0.344083 \| grad norm = 9.081e-01 | 0.344083 | 0.348814 | 0.53% | near-zero |
| 7 | step 80/90 \| objective = 0.344083 \| grad norm = 9.081e-01 | 0.9081 | 0.3889 | 57.17% |  |
| 7 | step 85/90 \| objective = 0.352045 \| grad norm = 4.737e-01 | 85 | 85 | 0 |  |
| 7 | step 85/90 \| objective = 0.352045 \| grad norm = 4.737e-01 | 90 | 90 | 0 |  |
| 7 | step 85/90 \| objective = 0.352045 \| grad norm = 4.737e-01 | 0.352045 | 0.351247 | 0.09% | near-zero |
| 7 | step 85/90 \| objective = 0.352045 \| grad norm = 4.737e-01 | 0.4737 | 0.4611 | 1.40% | near-zero |
| 7 | step 90/90 \| objective = 0.354362 \| grad norm = 4.982e-01 | 90 | 90 | 0 |  |
| 7 | step 90/90 \| objective = 0.354362 \| grad norm = 4.982e-01 | 90 | 90 | 0 |  |
| 7 | step 90/90 \| objective = 0.354362 \| grad norm = 4.982e-01 | 0.354362 | 0.35344 | 0.10% | near-zero |
| 7 | step 90/90 \| objective = 0.354362 \| grad norm = 4.982e-01 | 0.4982 | 0.4014 | 10.76% | near-zero |
| 8 | Single-layer nominal CE = 0.352976 | 0.352976 | 0.353895 | 0.26% |  |
| 8 | Single-layer x_center = 3.569 um | 3.569 | 3.573 | 0.11% |  |
| 8 | Single-layer tilt = 9.690 deg | 9.69 | 9.846 | 1.61% |  |
| 8 | Single-layer waist = 2.379 um | 2.379 | 2.378 | 0.04% |  |
| 9 | epochs = 90, learning_rate = 0.0030 | 90 | 90 | 0 |  |
| 9 | epochs = 90, learning_rate = 0.0030 | 0.003 | 0.003 | 0 | near-zero |
| 9 | step 01/90 \| objective = 0.198748 \| grad norm = 2.288e+00 | 1 | 1 | 0 |  |
| 9 | step 01/90 \| objective = 0.198748 \| grad norm = 2.288e+00 | 90 | 90 | 0 |  |
| 9 | step 01/90 \| objective = 0.198748 \| grad norm = 2.288e+00 | 0.198748 | 0.197795 | 0.11% | near-zero |
| 9 | step 01/90 \| objective = 0.198748 \| grad norm = 2.288e+00 | 2.288 | 2.274 | 0.61% |  |
| 9 | step 05/90 \| objective = 0.280728 \| grad norm = 1.845e+00 | 5 | 5 | 0 |  |
| 9 | step 05/90 \| objective = 0.280728 \| grad norm = 1.845e+00 | 90 | 90 | 0 |  |
| 9 | step 05/90 \| objective = 0.280728 \| grad norm = 1.845e+00 | 0.280728 | 0.279992 | 0.08% | near-zero |
| 9 | step 05/90 \| objective = 0.280728 \| grad norm = 1.845e+00 | 1.845 | 1.761 | 4.55% |  |
| 9 | step 10/90 \| objective = 0.308899 \| grad norm = 1.517e+00 | 10 | 10 | 0 |  |
| 9 | step 10/90 \| objective = 0.308899 \| grad norm = 1.517e+00 | 90 | 90 | 0 |  |
| 9 | step 10/90 \| objective = 0.308899 \| grad norm = 1.517e+00 | 0.308899 | 0.306541 | 0.26% | near-zero |
| 9 | step 10/90 \| objective = 0.308899 \| grad norm = 1.517e+00 | 1.517 | 1.928 | 27.09% |  |
| 9 | step 15/90 \| objective = 0.345474 \| grad norm = 2.142e+00 | 15 | 15 | 0 |  |
| 9 | step 15/90 \| objective = 0.345474 \| grad norm = 2.142e+00 | 90 | 90 | 0 |  |
| 9 | step 15/90 \| objective = 0.345474 \| grad norm = 2.142e+00 | 0.345474 | 0.34023 | 0.58% | near-zero |
| 9 | step 15/90 \| objective = 0.345474 \| grad norm = 2.142e+00 | 2.142 | 2.602 | 21.48% |  |
| 9 | step 20/90 \| objective = 0.377360 \| grad norm = 2.077e+00 | 20 | 20 | 0 |  |
| 9 | step 20/90 \| objective = 0.377360 \| grad norm = 2.077e+00 | 90 | 90 | 0 |  |
| 9 | step 20/90 \| objective = 0.377360 \| grad norm = 2.077e+00 | 0.37736 | 0.377611 | 0.03% | near-zero |
| 9 | step 20/90 \| objective = 0.377360 \| grad norm = 2.077e+00 | 2.077 | 1.691 | 18.58% |  |
| 9 | step 25/90 \| objective = 0.418240 \| grad norm = 1.794e+00 | 25 | 25 | 0 |  |
| 9 | step 25/90 \| objective = 0.418240 \| grad norm = 1.794e+00 | 90 | 90 | 0 |  |
| 9 | step 25/90 \| objective = 0.418240 \| grad norm = 1.794e+00 | 0.41824 | 0.417448 | 0.09% | near-zero |
| 9 | step 25/90 \| objective = 0.418240 \| grad norm = 1.794e+00 | 1.794 | 2.016 | 12.37% |  |
| 9 | step 30/90 \| objective = 0.451548 \| grad norm = 1.873e+00 | 30 | 30 | 0 |  |
| 9 | step 30/90 \| objective = 0.451548 \| grad norm = 1.873e+00 | 90 | 90 | 0 |  |
| 9 | step 30/90 \| objective = 0.451548 \| grad norm = 1.873e+00 | 0.451548 | 0.44732 | 0.47% | near-zero |
| 9 | step 30/90 \| objective = 0.451548 \| grad norm = 1.873e+00 | 1.873 | 1.949 | 4.06% |  |
| 9 | step 35/90 \| objective = 0.475962 \| grad norm = 4.308e+00 | 35 | 35 | 0 |  |
| 9 | step 35/90 \| objective = 0.475962 \| grad norm = 4.308e+00 | 90 | 90 | 0 |  |
| 9 | step 35/90 \| objective = 0.475962 \| grad norm = 4.308e+00 | 0.475962 | 0.488038 | 1.34% | near-zero |
| 9 | step 35/90 \| objective = 0.475962 \| grad norm = 4.308e+00 | 4.308 | 2.272 | 47.26% |  |
| 9 | step 40/90 \| objective = 0.517014 \| grad norm = 2.262e+00 | 40 | 40 | 0 |  |
| 9 | step 40/90 \| objective = 0.517014 \| grad norm = 2.262e+00 | 90 | 90 | 0 |  |
| 9 | step 40/90 \| objective = 0.517014 \| grad norm = 2.262e+00 | 0.517014 | 0.498757 | 2.03% | near-zero |
| 9 | step 40/90 \| objective = 0.517014 \| grad norm = 2.262e+00 | 2.262 | 5.385 | 138% |  |
| 9 | step 45/90 \| objective = 0.549722 \| grad norm = 2.061e+00 | 45 | 45 | 0 |  |
| 9 | step 45/90 \| objective = 0.549722 \| grad norm = 2.061e+00 | 90 | 90 | 0 |  |
| 9 | step 45/90 \| objective = 0.549722 \| grad norm = 2.061e+00 | 0.549722 | 0.547652 | 0.23% | near-zero |
| 9 | step 45/90 \| objective = 0.549722 \| grad norm = 2.061e+00 | 2.061 | 2.737 | 32.80% |  |
| 9 | step 50/90 \| objective = 0.582705 \| grad norm = 1.334e+00 | 50 | 50 | 0 |  |
| 9 | step 50/90 \| objective = 0.582705 \| grad norm = 1.334e+00 | 90 | 90 | 0 |  |
| 9 | step 50/90 \| objective = 0.582705 \| grad norm = 1.334e+00 | 0.582705 | 0.577976 | 0.53% | near-zero |
| 9 | step 50/90 \| objective = 0.582705 \| grad norm = 1.334e+00 | 1.334 | 2.863 | 115% |  |
| 9 | step 55/90 \| objective = 0.609034 \| grad norm = 1.604e+00 | 55 | 55 | 0 |  |
| 9 | step 55/90 \| objective = 0.609034 \| grad norm = 1.604e+00 | 90 | 90 | 0 |  |
| 9 | step 55/90 \| objective = 0.609034 \| grad norm = 1.604e+00 | 0.609034 | 0.60511 | 0.44% | near-zero |
| 9 | step 55/90 \| objective = 0.609034 \| grad norm = 1.604e+00 | 1.604 | 1.277 | 20.39% |  |
| 9 | step 60/90 \| objective = 0.624934 \| grad norm = 1.807e+00 | 60 | 60 | 0 |  |
| 9 | step 60/90 \| objective = 0.624934 \| grad norm = 1.807e+00 | 90 | 90 | 0 |  |
| 9 | step 60/90 \| objective = 0.624934 \| grad norm = 1.807e+00 | 0.624934 | 0.622926 | 0.22% | near-zero |
| 9 | step 60/90 \| objective = 0.624934 \| grad norm = 1.807e+00 | 1.807 | 1.081 | 40.18% |  |
| 9 | step 65/90 \| objective = 0.630303 \| grad norm = 3.563e+00 | 65 | 65 | 0 |  |
| 9 | step 65/90 \| objective = 0.630303 \| grad norm = 3.563e+00 | 90 | 90 | 0 |  |
| 9 | step 65/90 \| objective = 0.630303 \| grad norm = 3.563e+00 | 0.630303 | 0.641619 | 1.26% | near-zero |
| 9 | step 65/90 \| objective = 0.630303 \| grad norm = 3.563e+00 | 3.563 | 1.119 | 68.59% |  |
| 9 | step 70/90 \| objective = 0.660135 \| grad norm = 1.067e+00 | 70 | 70 | 0 |  |
| 9 | step 70/90 \| objective = 0.660135 \| grad norm = 1.067e+00 | 90 | 90 | 0 |  |
| 9 | step 70/90 \| objective = 0.660135 \| grad norm = 1.067e+00 | 0.660135 | 0.648094 | 1.34% | near-zero |
| 9 | step 70/90 \| objective = 0.660135 \| grad norm = 1.067e+00 | 1.067 | 2.126 | 99.25% |  |
| 9 | step 75/90 \| objective = 0.668581 \| grad norm = 1.005e+00 | 75 | 75 | 0 |  |
| 9 | step 75/90 \| objective = 0.668581 \| grad norm = 1.005e+00 | 90 | 90 | 0 |  |
| 9 | step 75/90 \| objective = 0.668581 \| grad norm = 1.005e+00 | 0.668581 | 0.665753 | 0.31% | near-zero |
| 9 | step 75/90 \| objective = 0.668581 \| grad norm = 1.005e+00 | 1.005 | 0.6975 | 30.60% |  |
| 9 | step 80/90 \| objective = 0.671168 \| grad norm = 2.370e+00 | 80 | 80 | 0 |  |
| 9 | step 80/90 \| objective = 0.671168 \| grad norm = 2.370e+00 | 90 | 90 | 0 |  |
| 9 | step 80/90 \| objective = 0.671168 \| grad norm = 2.370e+00 | 0.671168 | 0.673279 | 0.23% | near-zero |
| 9 | step 80/90 \| objective = 0.671168 \| grad norm = 2.370e+00 | 2.37 | 1.309 | 44.77% |  |
| 9 | step 85/90 \| objective = 0.680527 \| grad norm = 1.543e+00 | 85 | 85 | 0 |  |
| 9 | step 85/90 \| objective = 0.680527 \| grad norm = 1.543e+00 | 90 | 90 | 0 |  |
| 9 | step 85/90 \| objective = 0.680527 \| grad norm = 1.543e+00 | 0.680527 | 0.6765 | 0.45% | near-zero |
| 9 | step 85/90 \| objective = 0.680527 \| grad norm = 1.543e+00 | 1.543 | 0.5397 | 65.02% |  |
| 9 | step 90/90 \| objective = 0.685291 \| grad norm = 1.049e+00 | 90 | 90 | 0 |  |
| 9 | step 90/90 \| objective = 0.685291 \| grad norm = 1.049e+00 | 90 | 90 | 0 |  |
| 9 | step 90/90 \| objective = 0.685291 \| grad norm = 1.049e+00 | 0.685291 | 0.681204 | 0.45% | near-zero |
| 9 | step 90/90 \| objective = 0.685291 \| grad norm = 1.049e+00 | 1.049 | 0.8816 | 15.96% |  |
| 10 | Reflector nominal CE = 0.685303 | 0.685303 | 0.682811 | 0.36% |  |
| 10 | Reflector x_center = 3.762 um | 3.762 | 3.779 | 0.45% |  |
| 10 | Reflector tilt = 9.858 deg | 9.858 | 9.795 | 0.64% |  |
| 10 | Reflector waist = 2.373 um | 2.373 | 2.373 | 0 |  |
| 10 | Reflector period = 0.339 um | 0.339 | 0.337 | 0.59% |  |
| 10 | Reflector duty cycle = 0.417 | 0.417 | 0.42 | 0.72% |  |
| 11 | Single layer \| nominal CE = 0.352976 \| x_center = 3.569 um \| tilt = 9.690 deg \| waist = 2. | 0.352976 | 0.353895 | 0.26% |  |
| 11 | Single layer \| nominal CE = 0.352976 \| x_center = 3.569 um \| tilt = 9.690 deg \| waist = 2. | 3.569 | 3.573 | 0.11% |  |
| 11 | Single layer \| nominal CE = 0.352976 \| x_center = 3.569 um \| tilt = 9.690 deg \| waist = 2. | 9.69 | 9.846 | 1.61% |  |
| 11 | Single layer \| nominal CE = 0.352976 \| x_center = 3.569 um \| tilt = 9.690 deg \| waist = 2. | 2.379 | 2.378 | 0.04% |  |
| 11 | With reflector \| nominal CE = 0.685303 \| x_center = 3.762 um \| tilt = 9.858 deg \| waist =  | 0.685303 | 0.682811 | 0.36% |  |
| 11 | With reflector \| nominal CE = 0.685303 \| x_center = 3.762 um \| tilt = 9.858 deg \| waist =  | 3.762 | 3.779 | 0.45% |  |
| 11 | With reflector \| nominal CE = 0.685303 \| x_center = 3.762 um \| tilt = 9.858 deg \| waist =  | 9.858 | 9.795 | 0.64% |  |
| 11 | With reflector \| nominal CE = 0.685303 \| x_center = 3.762 um \| tilt = 9.858 deg \| waist =  | 2.373 | 2.373 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 6 | 0 | 8.92% |
| 8 | 0 | 7.79% |
| 10 | 0 | 8.07% |
| 12 | 0 | 4.14% |
| 13 | 0 | 6.28% |
| 13 | 1 | 7.78% |
| 14 | 0 | 4.28% |
