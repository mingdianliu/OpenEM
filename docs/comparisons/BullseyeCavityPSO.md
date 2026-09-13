# BullseyeCavityPSO

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/BullseyeCavityPSO/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 20 on the reference side, 20 on ours, 1 that did not line up |
| Numbers compared | 12 |
| Largest relative difference | 0 |
| Over 5% / over 20% | 0 / 0 |
| Figures | 7 on the reference side, 7 on ours |
| Largest pixel difference | 62.26% |

Machine-readable form of everything below: [`data/BullseyeCavityPSO.json`](data/BullseyeCavityPSO.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 10 | Best Parameters: r: 0.7463241440689826 p: 0.3812907960183133 w: | 0.746324 | 0.746324 | 0 |  |
| 10 | Best Parameters: r: 0.7463241440689826 p: 0.3812907960183133 w: | 0.381291 | 0.381291 | 0 |  |
| 10 | 0.1245569125769858 h: 0.26602630907617625 t_sio2: | 0.124557 | 0.124557 | 0 |  |
| 10 | 0.1245569125769858 h: 0.26602630907617625 t_sio2: | 0.266026 | 0.266026 | 0 |  |
| 10 | 0.6492437787750074 dcf: 1.0869886995014149 | 0.649244 | 0.649244 | 0 |  |
| 10 | 0.6492437787750074 dcf: 1.0869886995014149 | 1.08699 | 1.08699 | 0 |  |
| 11 | r_cav = 0.746 um | 0.746 | 0.746 | 0 |  |
| 11 | p_bragg = 0.381 um | 0.381 | 0.381 | 0 |  |
| 11 | w_bragg = 0.125 um | 0.125 | 0.125 | 0 |  |
| 11 | h = 0.266 | 0.266 | 0.266 | 0 |  |
| 11 | t_sio2 = 0.649 um | 0.649 | 0.649 | 0 |  |
| 11 | d_cf = 1.087 um | 1.087 | 1.087 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 4 | 0 | 6.94% |
| 6 | 0 | 62.26% |
| 11 | 0 | 3.96% |
| 12 | 0 | 5.97% |
| 14 | 0 | 4.31% |
| 17 | 0 | 4.09% |
| 19 | 0 | 3.15% |
