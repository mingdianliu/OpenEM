# PhotonicSpinSelector

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/PhotonicSpinSelector/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **trajectory divergence** |
| Reference output | reference rerun |
| Cells | 17 on the reference side, 17 on ours, 2 that did not line up |
| Numbers compared | 54 of 155 paired; the other 101 are near-zero, see below |
| Largest relative difference | 1,265% |
| Over 5% / over 20% | 5 / 4 |
| Figures | 4 on the reference side, 4 on ours |
| Largest pixel difference | 34.84% |

Machine-readable form of everything below: [`data/PhotonicSpinSelector.json`](data/PhotonicSpinSelector.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 5 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 4 | The total numbre of pixels is 435. | 435 | 435 | 0 |  |
| 14 | Starting DBS. Parameters: 435, Initial objective: 0.3477 | 435 | 435 | 0 |  |
| 14 | Starting DBS. Parameters: 435, Initial objective: 0.3477 | 0.3477 | 0.35 | 0.05% | near-zero |
| 14 | [Pass 1] Improved! Idx 319: 1->0 \| Highest objective value: 0.4293 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 319: 1->0 \| Highest objective value: 0.4293 | 319 | 319 | 0 |  |
| 14 | [Pass 1] Improved! Idx 319: 1->0 \| Highest objective value: 0.4293 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 319: 1->0 \| Highest objective value: 0.4293 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 319: 1->0 \| Highest objective value: 0.4293 | 0.4293 | 0.4308 | 0.05% | near-zero |
| 14 | [Pass 1] Improved! Idx 373: 0->1 \| Highest objective value: 0.4333 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 373: 0->1 \| Highest objective value: 0.4333 | 373 | 373 | 0 |  |
| 14 | [Pass 1] Improved! Idx 373: 0->1 \| Highest objective value: 0.4333 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 373: 0->1 \| Highest objective value: 0.4333 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 373: 0->1 \| Highest objective value: 0.4333 | 0.4333 | 0.4353 | 0.05% | near-zero |
| 14 | [Pass 1] Improved! Idx 189: 0->1 \| Highest objective value: 0.4983 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 189: 0->1 \| Highest objective value: 0.4983 | 189 | 189 | 0 |  |
| 14 | [Pass 1] Improved! Idx 189: 0->1 \| Highest objective value: 0.4983 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 189: 0->1 \| Highest objective value: 0.4983 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 189: 0->1 \| Highest objective value: 0.4983 | 0.4983 | 0.5003 | 0.11% | near-zero |
| 14 | [Pass 1] Improved! Idx 212: 1->0 \| Highest objective value: 0.5031 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 212: 1->0 \| Highest objective value: 0.5031 | 212 | 212 | 0 |  |
| 14 | [Pass 1] Improved! Idx 212: 1->0 \| Highest objective value: 0.5031 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 212: 1->0 \| Highest objective value: 0.5031 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 212: 1->0 \| Highest objective value: 0.5031 | 0.5031 | 0.5048 | 0.08% | near-zero |
| 14 | [Pass 1] Improved! Idx 232: 1->0 \| Highest objective value: 0.5213 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 232: 1->0 \| Highest objective value: 0.5213 | 232 | 232 | 0 |  |
| 14 | [Pass 1] Improved! Idx 232: 1->0 \| Highest objective value: 0.5213 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 232: 1->0 \| Highest objective value: 0.5213 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 232: 1->0 \| Highest objective value: 0.5213 | 0.5213 | 0.5223 | 0.04% | near-zero |
| 14 | [Pass 1] Improved! Idx 185: 0->1 \| Highest objective value: 0.6802 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 185: 0->1 \| Highest objective value: 0.6802 | 185 | 185 | 0 |  |
| 14 | [Pass 1] Improved! Idx 185: 0->1 \| Highest objective value: 0.6802 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 185: 0->1 \| Highest objective value: 0.6802 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 185: 0->1 \| Highest objective value: 0.6802 | 0.6802 | 0.6811 | 0.05% | near-zero |
| 14 | [Pass 1] Improved! Idx 75: 0->1 \| Highest objective value: 0.6975 | 1 | 1 | 0 |  |
| 14 | [Pass 1] Improved! Idx 75: 0->1 \| Highest objective value: 0.6975 | 75 | 75 | 0 |  |
| 14 | [Pass 1] Improved! Idx 75: 0->1 \| Highest objective value: 0.6975 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 75: 0->1 \| Highest objective value: 0.6975 | 1 | 1 | 0 |  |
| 14 | [Pass 1] Improved! Idx 75: 0->1 \| Highest objective value: 0.6975 | 0.6975 | 0.698 | 0.07% | near-zero |
| 14 | [Pass 1] Improved! Idx 114: 1->0 \| Highest objective value: 0.7368 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 114: 1->0 \| Highest objective value: 0.7368 | 114 | 114 | 0 |  |
| 14 | [Pass 1] Improved! Idx 114: 1->0 \| Highest objective value: 0.7368 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 114: 1->0 \| Highest objective value: 0.7368 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 114: 1->0 \| Highest objective value: 0.7368 | 0.7368 | 0.7372 | 0.04% | near-zero |
| 14 | [Pass 1] Improved! Idx 36: 0->1 \| Highest objective value: 0.8089 | 1 | 1 | 0 |  |
| 14 | [Pass 1] Improved! Idx 36: 0->1 \| Highest objective value: 0.8089 | 36 | 36 | 0 |  |
| 14 | [Pass 1] Improved! Idx 36: 0->1 \| Highest objective value: 0.8089 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 36: 0->1 \| Highest objective value: 0.8089 | 1 | 1 | 0 |  |
| 14 | [Pass 1] Improved! Idx 36: 0->1 \| Highest objective value: 0.8089 | 0.8089 | 0.8104 | 0.19% |  |
| 14 | [Pass 1] Improved! Idx 276: 1->0 \| Highest objective value: 0.8198 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 276: 1->0 \| Highest objective value: 0.8198 | 276 | 276 | 0 |  |
| 14 | [Pass 1] Improved! Idx 276: 1->0 \| Highest objective value: 0.8198 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 276: 1->0 \| Highest objective value: 0.8198 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 276: 1->0 \| Highest objective value: 0.8198 | 0.8198 | 0.8213 | 0.05% | near-zero |
| 14 | [Pass 1] Improved! Idx 9: 1->0 \| Highest objective value: 0.8335 | 1 | 1 | 0 |  |
| 14 | [Pass 1] Improved! Idx 9: 1->0 \| Highest objective value: 0.8335 | 9 | 9 | 0 |  |
| 14 | [Pass 1] Improved! Idx 9: 1->0 \| Highest objective value: 0.8335 | 1 | 1 | 0 |  |
| 14 | [Pass 1] Improved! Idx 9: 1->0 \| Highest objective value: 0.8335 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 9: 1->0 \| Highest objective value: 0.8335 | 0.8335 | 0.8347 | 0.14% |  |
| 14 | [Pass 1] Improved! Idx 242: 0->1 \| Highest objective value: 0.8394 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 242: 0->1 \| Highest objective value: 0.8394 | 242 | 242 | 0 |  |
| 14 | [Pass 1] Improved! Idx 242: 0->1 \| Highest objective value: 0.8394 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 242: 0->1 \| Highest objective value: 0.8394 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 242: 0->1 \| Highest objective value: 0.8394 | 0.8394 | 0.8407 | 0.05% | near-zero |
| 14 | [Pass 1] Improved! Idx 287: 1->0 \| Highest objective value: 0.8430 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 287: 1->0 \| Highest objective value: 0.8430 | 287 | 287 | 0 |  |
| 14 | [Pass 1] Improved! Idx 287: 1->0 \| Highest objective value: 0.8430 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 287: 1->0 \| Highest objective value: 0.8430 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 287: 1->0 \| Highest objective value: 0.8430 | 0.843 | 0.8438 | 0.03% | near-zero |
| 14 | [Pass 1] Improved! Idx 275: 0->1 \| Highest objective value: 0.8520 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 275: 0->1 \| Highest objective value: 0.8520 | 275 | 275 | 0 |  |
| 14 | [Pass 1] Improved! Idx 275: 0->1 \| Highest objective value: 0.8520 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 275: 0->1 \| Highest objective value: 0.8520 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 275: 0->1 \| Highest objective value: 0.8520 | 0.852 | 0.8526 | 0.02% | near-zero |
| 14 | [Pass 1] Improved! Idx 155: 1->0 \| Highest objective value: 0.8946 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 155: 1->0 \| Highest objective value: 0.8946 | 155 | 155 | 0 |  |
| 14 | [Pass 1] Improved! Idx 155: 1->0 \| Highest objective value: 0.8946 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 155: 1->0 \| Highest objective value: 0.8946 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 155: 1->0 \| Highest objective value: 0.8946 | 0.8946 | 0.8949 | 0.02% | near-zero |
| 14 | [Pass 1] Improved! Idx 171: 1->0 \| Highest objective value: 0.8975 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 171: 1->0 \| Highest objective value: 0.8975 | 171 | 171 | 0 |  |
| 14 | [Pass 1] Improved! Idx 171: 1->0 \| Highest objective value: 0.8975 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 171: 1->0 \| Highest objective value: 0.8975 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 171: 1->0 \| Highest objective value: 0.8975 | 0.8975 | 0.8979 | 0.02% | near-zero |
| 14 | [Pass 1] Improved! Idx 234: 1->0 \| Highest objective value: 0.9143 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 234: 1->0 \| Highest objective value: 0.9143 | 234 | 234 | 0 |  |
| 14 | [Pass 1] Improved! Idx 234: 1->0 \| Highest objective value: 0.9143 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 234: 1->0 \| Highest objective value: 0.9143 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 234: 1->0 \| Highest objective value: 0.9143 | 0.9143 | 0.9146 | 0.01% | near-zero |
| 14 | [Pass 1] Improved! Idx 48: 1->0 \| Highest objective value: 0.9218 | 1 | 1 | 0 |  |
| 14 | [Pass 1] Improved! Idx 48: 1->0 \| Highest objective value: 0.9218 | 48 | 48 | 0 |  |
| 14 | [Pass 1] Improved! Idx 48: 1->0 \| Highest objective value: 0.9218 | 1 | 1 | 0 |  |
| 14 | [Pass 1] Improved! Idx 48: 1->0 \| Highest objective value: 0.9218 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 48: 1->0 \| Highest objective value: 0.9218 | 0.9218 | 0.9219 | 0.01% |  |
| 14 | [Pass 1] Improved! Idx 41: 0->1 \| Highest objective value: 0.9434 | 1 | 1 | 0 |  |
| 14 | [Pass 1] Improved! Idx 41: 0->1 \| Highest objective value: 0.9434 | 41 | 41 | 0 |  |
| 14 | [Pass 1] Improved! Idx 41: 0->1 \| Highest objective value: 0.9434 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 41: 0->1 \| Highest objective value: 0.9434 | 1 | 1 | 0 |  |
| 14 | [Pass 1] Improved! Idx 41: 0->1 \| Highest objective value: 0.9434 | 0.9434 | 0.9438 | 0.04% |  |
| 14 | [Pass 1] Improved! Idx 175: 0->1 \| Highest objective value: 0.9505 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 175: 0->1 \| Highest objective value: 0.9505 | 175 | 175 | 0 |  |
| 14 | [Pass 1] Improved! Idx 175: 0->1 \| Highest objective value: 0.9505 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 175: 0->1 \| Highest objective value: 0.9505 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 175: 0->1 \| Highest objective value: 0.9505 | 0.9505 | 0.951 | 0.03% | near-zero |
| 14 | [Pass 1] Improved! Idx 400: 0->1 \| Highest objective value: 0.9567 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 400: 0->1 \| Highest objective value: 0.9567 | 400 | 400 | 0 |  |
| 14 | [Pass 1] Improved! Idx 400: 0->1 \| Highest objective value: 0.9567 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 400: 0->1 \| Highest objective value: 0.9567 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 400: 0->1 \| Highest objective value: 0.9567 | 0.9567 | 0.9573 | 0.02% | near-zero |
| 14 | [Pass 1] Improved! Idx 363: 1->0 \| Highest objective value: 0.9583 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 363: 1->0 \| Highest objective value: 0.9583 | 363 | 363 | 0 |  |
| 14 | [Pass 1] Improved! Idx 363: 1->0 \| Highest objective value: 0.9583 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 363: 1->0 \| Highest objective value: 0.9583 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 363: 1->0 \| Highest objective value: 0.9583 | 0.9583 | 0.959 | 0.02% | near-zero |
| 14 | [Pass 1] Improved! Idx 95: 0->1 \| Highest objective value: 0.9671 | 1 | 1 | 0 |  |
| 14 | [Pass 1] Improved! Idx 95: 0->1 \| Highest objective value: 0.9671 | 95 | 95 | 0 |  |
| 14 | [Pass 1] Improved! Idx 95: 0->1 \| Highest objective value: 0.9671 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 95: 0->1 \| Highest objective value: 0.9671 | 1 | 1 | 0 |  |
| 14 | [Pass 1] Improved! Idx 95: 0->1 \| Highest objective value: 0.9671 | 0.9671 | 0.9675 | 0.04% |  |
| 14 | [Pass 1] Improved! Idx 390: 1->0 \| Highest objective value: 0.9692 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 390: 1->0 \| Highest objective value: 0.9692 | 390 | 390 | 0 |  |
| 14 | [Pass 1] Improved! Idx 390: 1->0 \| Highest objective value: 0.9692 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 390: 1->0 \| Highest objective value: 0.9692 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 390: 1->0 \| Highest objective value: 0.9692 | 0.9692 | 0.9694 | 5.13e-03% | near-zero |
| 14 | [Pass 1] Improved! Idx 415: 1->0 \| Highest objective value: 0.9720 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 415: 1->0 \| Highest objective value: 0.9720 | 415 | 252 | 39.28% |  |
| 14 | [Pass 1] Improved! Idx 415: 1->0 \| Highest objective value: 0.9720 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 415: 1->0 \| Highest objective value: 0.9720 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 415: 1->0 \| Highest objective value: 0.9720 | 0.972 | 0.9694 | 0.06% | near-zero |
| 14 | [Pass 1] Improved! Idx 388: 1->0 \| Highest objective value: 0.9731 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 388: 1->0 \| Highest objective value: 0.9731 | 388 | 415 | 6.96% |  |
| 14 | [Pass 1] Improved! Idx 388: 1->0 \| Highest objective value: 0.9731 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 388: 1->0 \| Highest objective value: 0.9731 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 388: 1->0 \| Highest objective value: 0.9731 | 0.9731 | 0.9721 | 0.03% | near-zero |
| 14 | [Pass 1] Improved! Idx 107: 0->1 \| Highest objective value: 0.9819 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 107: 0->1 \| Highest objective value: 0.9819 | 107 | 388 | 263% |  |
| 14 | [Pass 1] Improved! Idx 107: 0->1 \| Highest objective value: 0.9819 | 0 | 1 | 93.46% | near-zero |
| 14 | [Pass 1] Improved! Idx 107: 0->1 \| Highest objective value: 0.9819 | 1 | 0 | 93.46% | near-zero |
| 14 | [Pass 1] Improved! Idx 107: 0->1 \| Highest objective value: 0.9819 | 0.9819 | 0.9734 | 0.79% | near-zero |
| 14 | [Pass 1] Improved! Idx 102: 1->0 \| Highest objective value: 0.9842 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 102: 1->0 \| Highest objective value: 0.9842 | 102 | 107 | 4.90% |  |
| 14 | [Pass 1] Improved! Idx 102: 1->0 \| Highest objective value: 0.9842 | 1 | 0 | 98.04% | near-zero |
| 14 | [Pass 1] Improved! Idx 102: 1->0 \| Highest objective value: 0.9842 | 0 | 1 | 98.04% | near-zero |
| 14 | [Pass 1] Improved! Idx 102: 1->0 \| Highest objective value: 0.9842 | 0.9842 | 0.9801 | 0.40% | near-zero |
| 14 | [Pass 1] Improved! Idx 423: 1->0 \| Highest objective value: 0.9868 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 423: 1->0 \| Highest objective value: 0.9868 | 423 | 326 | 22.93% |  |
| 14 | [Pass 1] Improved! Idx 423: 1->0 \| Highest objective value: 0.9868 | 1 | 1 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 423: 1->0 \| Highest objective value: 0.9868 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 423: 1->0 \| Highest objective value: 0.9868 | 0.9868 | 0.9833 | 0.08% | near-zero |
| 14 | [Pass 1] Improved! Idx 31: 1->0 \| Highest objective value: 0.9961 | 1 | 1 | 0 |  |
| 14 | [Pass 1] Improved! Idx 31: 1->0 \| Highest objective value: 0.9961 | 31 | 423 | 1,265% |  |
| 14 | [Pass 1] Improved! Idx 31: 1->0 \| Highest objective value: 0.9961 | 1 | 1 | 0 |  |
| 14 | [Pass 1] Improved! Idx 31: 1->0 \| Highest objective value: 0.9961 | 0 | 0 | 0 | near-zero |
| 14 | [Pass 1] Improved! Idx 31: 1->0 \| Highest objective value: 0.9961 | 0.9961 | 0.9868 | 0.93% |  |
| 14 | Termination condition reached (objective=0.9961 ≥ 0.99). Stopping optimization. | 0.9961 | 0.9903 | 0.58% |  |
| 14 | Termination condition reached (objective=0.9961 ≥ 0.99). Stopping optimization. | 0.99 | 0.99 | 0 |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 5 | 0 | 6.23% |
| 6 | 0 | 24.05% |
| 9 | 0 | 34.84% |
| 15 | 0 | 6.27% |
