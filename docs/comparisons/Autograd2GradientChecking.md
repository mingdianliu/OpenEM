# Autograd2GradientChecking

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd2GradientChecking/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Inverse Design |
| Verdict | **agree** |
| Reference output | official archived output |
| Cells | 18 on the reference side, 18 on ours, 2 that did not line up |
| Numbers compared | 32 |
| Largest relative difference | 0.14% |
| Over 5% / over 20% | 0 / 0 |
| Figures | 1 on the reference side, 1 on ours |
| Largest pixel difference | 5.00% |

Machine-readable form of everything below: [`data/Autograd2GradientChecking.json`](data/Autograd2GradientChecking.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 3 | T (tmm) = 0.786 | 0.786 | 0.786 | 0 |  |
| 5 | gradient w.r.t. eps (tmm)  = [-0.2766323   0.01377339 -0.2032054  -0.28999361] | -0.276632 | -0.276632 | 0 |  |
| 5 | gradient w.r.t. eps (tmm)  = [-0.2766323   0.01377339 -0.2032054  -0.28999361] | 0.0137734 | 0.0137734 | 0 |  |
| 5 | gradient w.r.t. eps (tmm)  = [-0.2766323   0.01377339 -0.2032054  -0.28999361] | -0.203205 | -0.203205 | 0 |  |
| 5 | gradient w.r.t. eps (tmm)  = [-0.2766323   0.01377339 -0.2032054  -0.28999361] | -0.289994 | -0.289994 | 0 |  |
| 5 | gradient w.r.t. ds  (tmm)  = [-1.75199732 -0.21552416  1.00729645 -2.08209951] | -1.752 | -1.752 | 0 |  |
| 5 | gradient w.r.t. ds  (tmm)  = [-1.75199732 -0.21552416  1.00729645 -2.08209951] | -0.215524 | -0.215524 | 0 |  |
| 5 | gradient w.r.t. ds  (tmm)  = [-1.75199732 -0.21552416  1.00729645 -2.08209951] | 1.0073 | 1.0073 | 0 |  |
| 5 | gradient w.r.t. ds  (tmm)  = [-1.75199732 -0.21552416  1.00729645 -2.08209951] | -2.0821 | -2.0821 | 0 |  |
| 10 | 0.7851633658072134 | 0.785163 | 0.785057 | 0.01% |  |
| 14 | T (tmm)  = 0.78581 | 0.78581 | 0.78581 | 0 |  |
| 14 | T (FDTD) = 0.78516 | 0.78516 | 0.78506 | 0.01% |  |
| 15 | grad_eps (tmm)  = [-0.2766323   0.01377339 -0.2032054  -0.28999361] | -0.276632 | -0.276632 | 0 |  |
| 15 | grad_eps (tmm)  = [-0.2766323   0.01377339 -0.2032054  -0.28999361] | 0.0137734 | 0.0137734 | 0 |  |
| 15 | grad_eps (tmm)  = [-0.2766323   0.01377339 -0.2032054  -0.28999361] | -0.203205 | -0.203205 | 0 |  |
| 15 | grad_eps (tmm)  = [-0.2766323   0.01377339 -0.2032054  -0.28999361] | -0.289994 | -0.289994 | 0 |  |
| 15 | grad_ds  (tmm)  = [-1.75199732 -0.21552416  1.00729645 -2.08209951] | -1.752 | -1.752 | 0 |  |
| 15 | grad_ds  (tmm)  = [-1.75199732 -0.21552416  1.00729645 -2.08209951] | -0.215524 | -0.215524 | 0 |  |
| 15 | grad_ds  (tmm)  = [-1.75199732 -0.21552416  1.00729645 -2.08209951] | 1.0073 | 1.0073 | 0 |  |
| 15 | grad_ds  (tmm)  = [-1.75199732 -0.21552416  1.00729645 -2.08209951] | -2.0821 | -2.0821 | 0 |  |
| 17 | grad_eps (tmm)  = [-0.61534061  0.03063751 -0.45200988 -0.64506151] | -0.615341 | -0.615341 | 0 |  |
| 17 | grad_eps (tmm)  = [-0.61534061  0.03063751 -0.45200988 -0.64506151] | 0.0306375 | 0.0306375 | 0 |  |
| 17 | grad_eps (tmm)  = [-0.61534061  0.03063751 -0.45200988 -0.64506151] | -0.45201 | -0.45201 | 0 |  |
| 17 | grad_eps (tmm)  = [-0.61534061  0.03063751 -0.45200988 -0.64506151] | -0.645062 | -0.645062 | 0 |  |
| 17 | grad_eps (FDTD)  = [-0.61546297  0.03088108 -0.45142613 -0.64534188] | -0.615463 | -0.615471 | 1.36e-03% |  |
| 17 | grad_eps (FDTD)  = [-0.61546297  0.03088108 -0.45142613 -0.64534188] | 0.0308811 | 0.0309231 | 0.14% |  |
| 17 | grad_eps (FDTD)  = [-0.61546297  0.03088108 -0.45142613 -0.64534188] | -0.451426 | -0.451391 | 7.70e-03% |  |
| 17 | grad_eps (FDTD)  = [-0.61546297  0.03088108 -0.45142613 -0.64534188] | -0.645342 | -0.645356 | 2.22e-03% |  |
| 17 | grad_ds  (tmm)  = [-0.60214521 -0.07407365  0.34619844 -0.71559827] | -0.602145 | -0.602145 | 0 |  |
| 17 | grad_ds  (tmm)  = [-0.60214521 -0.07407365  0.34619844 -0.71559827] | -0.0740737 | -0.0740737 | 0 |  |
| 17 | grad_ds  (tmm)  = [-0.60214521 -0.07407365  0.34619844 -0.71559827] | 0.346198 | 0.346198 | 0 |  |
| 17 | grad_ds  (tmm)  = [-0.60214521 -0.07407365  0.34619844 -0.71559827] | -0.715598 | -0.715598 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 7 | 0 | 5.00% |
