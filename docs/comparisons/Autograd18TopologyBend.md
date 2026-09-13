# Autograd18TopologyBend

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd18TopologyBend/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **trajectory divergence** |
| Reference output | official archived output |
| Cells | 21 on the reference side, 21 on ours, 0 that did not line up |
| Numbers compared | 100 |
| Largest relative difference | 71.65% |
| Over 5% / over 20% | 16 / 4 |
| Figures | 32 on the reference side, 32 on ours |
| Largest pixel difference | 10.19% |

Machine-readable form of everything below: [`data/Autograd18TopologyBend.json`](data/Autograd18TopologyBend.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 16 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 13 | step = 1 | 1 | 1 | 0 |  |
| 13 | beta = 5.00 | 5 | 5 | 0 |  |
| 13 | J = -9.8748e-01 | -0.98748 | -0.98743 | 5.06e-03% |  |
| 13 | grad_norm = 1.1811e-02 | 0.011811 | 0.011862 | 0.43% |  |
| 13 | step = 2 | 2 | 2 | 0 |  |
| 13 | beta = 5.62 | 5.62 | 5.62 | 0 |  |
| 13 | J = -1.7849e-01 | -0.17849 | -0.18056 | 1.16% |  |
| 13 | grad_norm = 1.2433e-02 | 0.012433 | 0.012509 | 0.61% |  |
| 13 | step = 3 | 3 | 3 | 0 |  |
| 13 | beta = 6.25 | 6.25 | 6.25 | 0 |  |
| 13 | J = 1.0090e-02 | 0.01009 | 0.010592 | 4.98% |  |
| 13 | grad_norm = 1.3535e-02 | 0.013535 | 0.013518 | 0.13% |  |
| 13 | step = 4 | 4 | 4 | 0 |  |
| 13 | beta = 6.88 | 6.88 | 6.88 | 0 |  |
| 13 | J = 1.8153e-01 | 0.18153 | 0.19066 | 5.03% |  |
| 13 | grad_norm = 1.4126e-02 | 0.014126 | 0.014471 | 2.44% |  |
| 13 | step = 5 | 5 | 5 | 0 |  |
| 13 | beta = 7.50 | 7.5 | 7.5 | 0 |  |
| 13 | J = 3.0211e-01 | 0.30211 | 0.33717 | 11.61% |  |
| 13 | grad_norm = 1.7756e-02 | 0.017756 | 0.017768 | 0.07% |  |
| 13 | step = 6 | 6 | 6 | 0 |  |
| 13 | beta = 8.12 | 8.12 | 8.12 | 0 |  |
| 13 | J = 4.0307e-01 | 0.40307 | 0.41915 | 3.99% |  |
| 13 | grad_norm = 1.5350e-02 | 0.01535 | 0.015702 | 2.29% |  |
| 13 | step = 7 | 7 | 7 | 0 |  |
| 13 | beta = 8.75 | 8.75 | 8.75 | 0 |  |
| 13 | J = 4.7209e-01 | 0.47209 | 0.46814 | 0.84% |  |
| 13 | grad_norm = 1.4658e-02 | 0.014658 | 0.015178 | 3.55% |  |
| 13 | step = 8 | 8 | 8 | 0 |  |
| 13 | beta = 9.38 | 9.38 | 9.38 | 0 |  |
| 13 | J = 5.4106e-01 | 0.54106 | 0.53533 | 1.06% |  |
| 13 | grad_norm = 1.0828e-02 | 0.010828 | 0.0092249 | 14.81% |  |
| 13 | step = 9 | 9 | 9 | 0 |  |
| 13 | beta = 10.00 | 10 | 10 | 0 |  |
| 13 | J = 5.7739e-01 | 0.57739 | 0.57727 | 0.02% |  |
| 13 | grad_norm = 1.2721e-02 | 0.012721 | 0.009805 | 22.92% |  |
| 13 | step = 10 | 10 | 10 | 0 |  |
| 13 | beta = 10.62 | 10.62 | 10.62 | 0 |  |
| 13 | J = 6.1888e-01 | 0.61888 | 0.61476 | 0.67% |  |
| 13 | grad_norm = 9.5648e-03 | 0.0095648 | 0.008099 | 15.32% |  |
| 13 | step = 11 | 11 | 11 | 0 |  |
| 13 | beta = 11.25 | 11.25 | 11.25 | 0 |  |
| 13 | J = 6.4623e-01 | 0.64623 | 0.64206 | 0.65% |  |
| 13 | grad_norm = 6.8398e-03 | 0.0068398 | 0.0069365 | 1.41% |  |
| 13 | step = 12 | 12 | 12 | 0 |  |
| 13 | beta = 11.88 | 11.88 | 11.88 | 0 |  |
| 13 | J = 6.5977e-01 | 0.65977 | 0.66199 | 0.34% |  |
| 13 | grad_norm = 7.1572e-03 | 0.0071572 | 0.0075863 | 6.00% |  |
| 13 | step = 13 | 13 | 13 | 0 |  |
| 13 | beta = 12.50 | 12.5 | 12.5 | 0 |  |
| 13 | J = 6.7336e-01 | 0.67336 | 0.67664 | 0.49% |  |
| 13 | grad_norm = 4.9686e-03 | 0.0049686 | 0.004978 | 0.19% |  |
| 13 | step = 14 | 14 | 14 | 0 |  |
| 13 | beta = 13.12 | 13.12 | 13.12 | 0 |  |
| 13 | J = 6.8527e-01 | 0.68527 | 0.68719 | 0.28% |  |
| 13 | grad_norm = 4.6140e-03 | 0.004614 | 0.0041823 | 9.36% |  |
| 13 | step = 15 | 15 | 15 | 0 |  |
| 13 | beta = 13.75 | 13.75 | 13.75 | 0 |  |
| 13 | J = 6.9623e-01 | 0.69623 | 0.69804 | 0.26% |  |
| 13 | grad_norm = 5.1576e-03 | 0.0051576 | 0.0045454 | 11.87% |  |
| 13 | step = 16 | 16 | 16 | 0 |  |
| 13 | beta = 14.38 | 14.38 | 14.38 | 0 |  |
| 13 | J = 7.0596e-01 | 0.70596 | 0.7072 | 0.18% |  |
| 13 | grad_norm = 6.1804e-03 | 0.0061804 | 0.0057647 | 6.73% |  |
| 13 | step = 17 | 17 | 17 | 0 |  |
| 13 | beta = 15.00 | 15 | 15 | 0 |  |
| 13 | J = 7.1413e-01 | 0.71413 | 0.71724 | 0.44% |  |
| 13 | grad_norm = 5.7893e-03 | 0.0057893 | 0.0046452 | 19.76% |  |
| 13 | step = 18 | 18 | 18 | 0 |  |
| 13 | beta = 15.62 | 15.62 | 15.62 | 0 |  |
| 13 | J = 7.2290e-01 | 0.7229 | 0.72401 | 0.15% |  |
| 13 | grad_norm = 4.6860e-03 | 0.004686 | 0.00519 | 10.76% |  |
| 13 | step = 19 | 19 | 19 | 0 |  |
| 13 | beta = 16.25 | 16.25 | 16.25 | 0 |  |
| 13 | J = 7.2951e-01 | 0.72951 | 0.73101 | 0.21% |  |
| 13 | grad_norm = 4.7624e-03 | 0.0047624 | 0.0050916 | 6.91% |  |
| 13 | step = 20 | 20 | 20 | 0 |  |
| 13 | beta = 16.88 | 16.88 | 16.88 | 0 |  |
| 13 | J = 7.3623e-01 | 0.73623 | 0.73677 | 0.07% |  |
| 13 | grad_norm = 3.3963e-03 | 0.0033963 | 0.0047952 | 41.19% |  |
| 13 | step = 21 | 21 | 21 | 0 |  |
| 13 | beta = 17.50 | 17.5 | 17.5 | 0 |  |
| 13 | J = 7.4142e-01 | 0.74142 | 0.7433 | 0.25% |  |
| 13 | grad_norm = 2.8210e-03 | 0.002821 | 0.0033816 | 19.87% |  |
| 13 | step = 22 | 22 | 22 | 0 |  |
| 13 | beta = 18.12 | 18.12 | 18.12 | 0 |  |
| 13 | J = 7.4584e-01 | 0.74584 | 0.74784 | 0.27% |  |
| 13 | grad_norm = 3.0746e-03 | 0.0030746 | 0.0029994 | 2.45% |  |
| 13 | step = 23 | 23 | 23 | 0 |  |
| 13 | beta = 18.75 | 18.75 | 18.75 | 0 |  |
| 13 | J = 7.4950e-01 | 0.7495 | 0.75134 | 0.25% |  |
| 13 | grad_norm = 3.4878e-03 | 0.0034878 | 0.0035883 | 2.88% |  |
| 13 | step = 24 | 24 | 24 | 0 |  |
| 13 | beta = 19.38 | 19.38 | 19.38 | 0 |  |
| 13 | J = 7.5294e-01 | 0.75294 | 0.75468 | 0.23% |  |
| 13 | grad_norm = 2.5567e-03 | 0.0025567 | 0.0036158 | 41.42% |  |
| 13 | step = 25 | 25 | 25 | 0 |  |
| 13 | beta = 20.00 | 20 | 20 | 0 |  |
| 13 | J = 7.5589e-01 | 0.75589 | 0.75752 | 0.22% |  |
| 13 | grad_norm = 2.9656e-03 | 0.0029656 | 0.0050904 | 71.65% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 3 | 0 | 3.66% |
| 6 | 0 | 4.81% |
| 8 | 0 | 5.56% |
| 13 | 0 | 8.86% |
| 13 | 1 | 10.19% |
| 13 | 2 | 9.97% |
| 13 | 3 | 10.02% |
| 13 | 4 | 10.05% |
| 13 | 5 | 10.08% |
| 13 | 6 | 10.09% |
| 13 | 7 | 10.09% |
| 13 | 8 | 10.10% |
| 13 | 9 | 10.07% |
| 13 | 10 | 10.04% |
| 13 | 11 | 10.06% |
| 13 | 12 | 10.01% |
| 13 | 13 | 10.03% |
| 13 | 14 | 10.00% |
| 13 | 15 | 10.00% |
| 13 | 16 | 9.96% |
| 13 | 17 | 9.94% |
| 13 | 18 | 9.92% |
| 13 | 19 | 9.93% |
| 13 | 20 | 9.92% |
| 13 | 21 | 9.91% |
| 13 | 22 | 9.88% |
| 13 | 23 | 9.90% |
| 13 | 24 | 9.92% |
| 14 | 0 | 3.34% |
| 15 | 0 | 5.39% |
| 17 | 0 | 5.94% |
| 18 | 0 | 4.56% |
