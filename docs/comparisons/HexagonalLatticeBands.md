# HexagonalLatticeBands

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/HexagonalLatticeBands/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 12 on the reference side, 12 on ours, 0 that did not line up |
| Numbers compared | 20 |
| Largest relative difference | 0 |
| Over 5% / over 20% | 0 / 0 |
| Figures | 5 on the reference side, 5 on ours |
| Largest pixel difference | 8.06% |

Machine-readable form of everything below: [`data/HexagonalLatticeBands.json`](data/HexagonalLatticeBands.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 4 | value, which will set it to 1 internally. A value of 1 means that | 1 | 1 | 0 |  |
| 4 | value, which will set it to 1 internally. A value of 1 means that | 1 | 1 | 0 |  |
| 4 | 'interval > 1' or by choosing alternative 'start' and 'stop' values | 1 | 1 | 0 |  |
| 4 | explicitly setting 'interval=1' in the monitor. | 1 | 1 | 0 |  |
| 4 | value, which will set it to 1 internally. A value of 1 means that | 1 | 1 | 0 |  |
| 4 | value, which will set it to 1 internally. A value of 1 means that | 1 | 1 | 0 |  |
| 4 | 'interval > 1' or by choosing alternative 'start' and 'stop' values | 1 | 1 | 0 |  |
| 4 | explicitly setting 'interval=1' in the monitor. | 1 | 1 | 0 |  |
| 4 | value, which will set it to 1 internally. A value of 1 means that | 1 | 1 | 0 |  |
| 4 | value, which will set it to 1 internally. A value of 1 means that | 1 | 1 | 0 |  |
| 4 | 'interval > 1' or by choosing alternative 'start' and 'stop' values | 1 | 1 | 0 |  |
| 4 | explicitly setting 'interval=1' in the monitor. | 1 | 1 | 0 |  |
| 4 | value, which will set it to 1 internally. A value of 1 means that | 1 | 1 | 0 |  |
| 4 | value, which will set it to 1 internally. A value of 1 means that | 1 | 1 | 0 |  |
| 4 | 'interval > 1' or by choosing alternative 'start' and 'stop' values | 1 | 1 | 0 |  |
| 4 | explicitly setting 'interval=1' in the monitor. | 1 | 1 | 0 |  |
| 4 | value, which will set it to 1 internally. A value of 1 means that | 1 | 1 | 0 |  |
| 4 | value, which will set it to 1 internally. A value of 1 means that | 1 | 1 | 0 |  |
| 4 | 'interval > 1' or by choosing alternative 'start' and 'stop' values | 1 | 1 | 0 |  |
| 4 | explicitly setting 'interval=1' in the monitor. | 1 | 1 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 1 | 0 | 8.06% |
| 6 | 0 | 5.09% |
| 7 | 0 | 2.40% |
| 9 | 0 | 2.59% |
| 11 | 0 | 6.51% |
