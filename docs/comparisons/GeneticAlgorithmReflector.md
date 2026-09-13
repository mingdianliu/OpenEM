# GeneticAlgorithmReflector

[This notebook on the Tidy3D example library](https://www.flexcompute.com/tidy3d/examples/notebooks/GeneticAlgorithmReflector/). The notebook itself is Flexcompute's and is **not** redistributed here; what follows is our comparison of its output against OpenEM's, number by number and figure by figure.

| | |
|---|---|
| Category | Example Library |
| Verdict | **agree** |
| Reference output | reference rerun |
| Cells | 17 on the reference side, 17 on ours, 2 that did not line up |
| Numbers compared | 511 of 802 paired; the other 291 are near-zero, see below |
| Largest relative difference | 425% |
| Over 5% / over 20% | 97 / 39 |
| Figures | 5 on the reference side, 5 on ours |
| Largest pixel difference | 9.87% |

Machine-readable form of everything below: [`data/GeneticAlgorithmReflector.json`](data/GeneticAlgorithmReflector.json). The verdict and the one-line note for this notebook are in the summary table, [`../validation.md`](../validation.md); what the fields here mean is in [`README.md`](README.md).

> 97 of the scraped numbers below differ by more than 5%. The verdict for this notebook is not derived from them: it rests on the conclusion quantity named in its row of [`../validation.md`](../validation.md), plus a visual comparison of the figures. [What produces a large scraped difference](README.md#why-a-scraped-number-can-differ-wildly-while-the-notebook-still-agrees) explains the three causes.

## Numbers

One row per number the two sides printed in the same place. `near-zero` marks a value at least 100 times smaller than the largest number on its own printed line; those are excluded from the counts above, because a relative difference on them is not meaningful.

| Cell | Printed line | Reference | OpenEM | Rel. diff | |
|---:|---|---:|---:|---:|---|
| 9 | solutions_per_pop: 30 | 30 | 30 | 0 |  |
| 9 | n_generations: 25 | 25 | 25 | 0 |  |
| 9 | n_parents_mating: 10 | 10 | 10 | 0 |  |
| 9 | stop_criteria_number: 6.0 | 6 | 6 | 0 |  |
| 9 | keep_parents: -1 | -1 | -1 | 0 |  |
| 9 | keep_elitism: 1 | 1 | 1 | 0 |  |
| 9 | crossover_prob: 0.7 | 0.7 | 0.7 | 0 |  |
| 9 | No. of Parameters: 162 | 162 | 162 | 0 |  |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 0 | 0 | 0 | near-zero |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 1 | 1 | 0 |  |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 2 | 2 | 0 |  |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 3 | 3 | 0 |  |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 4 | 4 | 0 |  |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 5 | 5 | 0 |  |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 6 | 6 | 0 |  |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 7 | 7 | 0 |  |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 8 | 8 | 0 |  |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 9 | 9 | 0 |  |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 10 | 10 | 0 |  |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 11 | 11 | 0 |  |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 12 | 12 | 0 |  |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 13 | 13 | 0 |  |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 14 | 14 | 0 |  |
| 9 | Parameters: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, | 15 | 15 | 0 |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 16 | 84 | 425% |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 17 | 85 | 400% |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 18 | 86 | 378% |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 19 | 87 | 358% |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 20 | 88 | 340% |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 21 | 89 | 324% |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 22 | 90 | 309% |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 23 | 91 | 296% |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 24 | 92 | 283% |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 25 | 93 | 272% |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 26 | 94 | 262% |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 27 | 95 | 252% |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 28 | 96 | 243% |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 29 | 97 | 234% |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 30 | 98 | 227% |  |
| 9 | 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, | 31 | 99 | 219% |  |
| 9 | 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, | 110 | 100 | 9.09% |  |
| 9 | 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, | 111 | 101 | 9.01% |  |
| 9 | 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, | 112 | 102 | 8.93% |  |
| 9 | 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, | 113 | 103 | 8.85% |  |
| 9 | 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, | 114 | 104 | 8.77% |  |
| 9 | 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, | 115 | 105 | 8.70% |  |
| 9 | 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, | 116 | 106 | 8.62% |  |
| 9 | 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, | 117 | 107 | 8.55% |  |
| 9 | 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, | 118 | 108 | 8.47% |  |
| 9 | 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, | 119 | 109 | 8.40% |  |
| 9 | 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, | 120 | 110 | 8.33% |  |
| 9 | 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, | 121 | 111 | 8.26% |  |
| 9 | 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, | 122 | 112 | 8.20% |  |
| 9 | 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, | 123 | 113 | 8.13% |  |
| 9 | 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, | 124 | 114 | 8.06% |  |
| 9 | 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, | 125 | 115 | 8.00% |  |
| 9 | 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, | 126 | 116 | 7.94% |  |
| 9 | 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, | 127 | 117 | 7.87% |  |
| 9 | 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, | 128 | 118 | 7.81% |  |
| 9 | 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, | 129 | 119 | 7.75% |  |
| 9 | 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, | 130 | 120 | 7.69% |  |
| 9 | 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, | 131 | 121 | 7.63% |  |
| 9 | 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, | 132 | 122 | 7.58% |  |
| 9 | 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, | 133 | 123 | 7.52% |  |
| 9 | 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, | 134 | 124 | 7.46% |  |
| 9 | 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, | 135 | 125 | 7.41% |  |
| 9 | 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, | 136 | 126 | 7.35% |  |
| 9 | 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, | 137 | 127 | 7.30% |  |
| 9 | 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, | 138 | 128 | 7.25% |  |
| 9 | 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, | 139 | 129 | 7.19% |  |
| 9 | 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, | 140 | 130 | 7.14% |  |
| 9 | 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, | 141 | 131 | 7.09% |  |
| 9 | 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, | 142 | 132 | 7.04% |  |
| 9 | 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, | 143 | 133 | 6.99% |  |
| 9 | 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, | 144 | 134 | 6.94% |  |
| 9 | 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, | 145 | 135 | 6.90% |  |
| 9 | 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, | 146 | 136 | 6.85% |  |
| 9 | 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, | 147 | 137 | 6.80% |  |
| 9 | 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, | 148 | 138 | 6.76% |  |
| 9 | 0: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 0: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 0: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 1: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 1: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 1: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 2: ParameterInt (0, 1) | 2 | 2 | 0 |  |
| 9 | 2: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 2: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 3: ParameterInt (0, 1) | 3 | 3 | 0 |  |
| 9 | 3: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 3: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 4: ParameterInt (0, 1) | 4 | 4 | 0 |  |
| 9 | 4: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 4: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 5: ParameterInt (0, 1) | 5 | 5 | 0 |  |
| 9 | 5: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 5: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 6: ParameterInt (0, 1) | 6 | 6 | 0 |  |
| 9 | 6: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 6: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 7: ParameterInt (0, 1) | 7 | 7 | 0 |  |
| 9 | 7: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 7: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 8: ParameterInt (0, 1) | 8 | 8 | 0 |  |
| 9 | 8: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 8: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 9: ParameterInt (0, 1) | 9 | 9 | 0 |  |
| 9 | 9: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 9: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 10: ParameterInt (0, 1) | 10 | 10 | 0 |  |
| 9 | 10: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 10: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 11: ParameterInt (0, 1) | 11 | 11 | 0 |  |
| 9 | 11: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 11: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 12: ParameterInt (0, 1) | 12 | 12 | 0 |  |
| 9 | 12: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 12: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 13: ParameterInt (0, 1) | 13 | 13 | 0 |  |
| 9 | 13: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 13: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 14: ParameterInt (0, 1) | 14 | 14 | 0 |  |
| 9 | 14: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 14: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 15: ParameterInt (0, 1) | 15 | 15 | 0 |  |
| 9 | 15: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 15: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 16: ParameterInt (0, 1) | 16 | 16 | 0 |  |
| 9 | 16: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 16: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 17: ParameterInt (0, 1) | 17 | 17 | 0 |  |
| 9 | 17: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 17: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 18: ParameterInt (0, 1) | 18 | 18 | 0 |  |
| 9 | 18: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 18: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 19: ParameterInt (0, 1) | 19 | 19 | 0 |  |
| 9 | 19: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 19: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 20: ParameterInt (0, 1) | 20 | 20 | 0 |  |
| 9 | 20: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 20: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 21: ParameterInt (0, 1) | 21 | 21 | 0 |  |
| 9 | 21: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 21: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 22: ParameterInt (0, 1) | 22 | 22 | 0 |  |
| 9 | 22: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 22: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 23: ParameterInt (0, 1) | 23 | 23 | 0 |  |
| 9 | 23: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 23: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 24: ParameterInt (0, 1) | 24 | 24 | 0 |  |
| 9 | 24: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 24: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 25: ParameterInt (0, 1) | 25 | 25 | 0 |  |
| 9 | 25: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 25: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 26: ParameterInt (0, 1) | 26 | 26 | 0 |  |
| 9 | 26: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 26: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 27: ParameterInt (0, 1) | 27 | 27 | 0 |  |
| 9 | 27: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 27: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 28: ParameterInt (0, 1) | 28 | 28 | 0 |  |
| 9 | 28: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 28: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 29: ParameterInt (0, 1) | 29 | 29 | 0 |  |
| 9 | 29: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 29: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 30: ParameterInt (0, 1) | 30 | 30 | 0 |  |
| 9 | 30: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 30: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 31: ParameterInt (0, 1) | 31 | 31 | 0 |  |
| 9 | 31: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 31: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 32: ParameterInt (0, 1) | 32 | 32 | 0 |  |
| 9 | 32: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 32: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 33: ParameterInt (0, 1) | 33 | 33 | 0 |  |
| 9 | 33: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 33: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 34: ParameterInt (0, 1) | 34 | 34 | 0 |  |
| 9 | 34: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 34: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 35: ParameterInt (0, 1) | 35 | 35 | 0 |  |
| 9 | 35: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 35: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 36: ParameterInt (0, 1) | 36 | 36 | 0 |  |
| 9 | 36: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 36: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 37: ParameterInt (0, 1) | 37 | 37 | 0 |  |
| 9 | 37: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 37: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 38: ParameterInt (0, 1) | 38 | 38 | 0 |  |
| 9 | 38: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 38: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 39: ParameterInt (0, 1) | 39 | 39 | 0 |  |
| 9 | 39: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 39: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 40: ParameterInt (0, 1) | 40 | 40 | 0 |  |
| 9 | 40: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 40: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 41: ParameterInt (0, 1) | 41 | 41 | 0 |  |
| 9 | 41: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 41: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 42: ParameterInt (0, 1) | 42 | 42 | 0 |  |
| 9 | 42: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 42: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 43: ParameterInt (0, 1) | 43 | 43 | 0 |  |
| 9 | 43: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 43: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 44: ParameterInt (0, 1) | 44 | 44 | 0 |  |
| 9 | 44: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 44: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 45: ParameterInt (0, 1) | 45 | 45 | 0 |  |
| 9 | 45: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 45: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 46: ParameterInt (0, 1) | 46 | 46 | 0 |  |
| 9 | 46: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 46: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 47: ParameterInt (0, 1) | 47 | 47 | 0 |  |
| 9 | 47: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 47: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 48: ParameterInt (0, 1) | 48 | 48 | 0 |  |
| 9 | 48: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 48: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 49: ParameterInt (0, 1) | 49 | 49 | 0 |  |
| 9 | 49: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 49: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 50: ParameterInt (0, 1) | 50 | 50 | 0 |  |
| 9 | 50: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 50: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 51: ParameterInt (0, 1) | 51 | 51 | 0 |  |
| 9 | 51: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 51: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 52: ParameterInt (0, 1) | 52 | 52 | 0 |  |
| 9 | 52: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 52: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 53: ParameterInt (0, 1) | 53 | 53 | 0 |  |
| 9 | 53: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 53: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 54: ParameterInt (0, 1) | 54 | 54 | 0 |  |
| 9 | 54: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 54: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 55: ParameterInt (0, 1) | 55 | 55 | 0 |  |
| 9 | 55: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 55: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 56: ParameterInt (0, 1) | 56 | 56 | 0 |  |
| 9 | 56: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 56: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 57: ParameterInt (0, 1) | 57 | 57 | 0 |  |
| 9 | 57: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 57: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 58: ParameterInt (0, 1) | 58 | 58 | 0 |  |
| 9 | 58: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 58: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 59: ParameterInt (0, 1) | 59 | 59 | 0 |  |
| 9 | 59: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 59: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 60: ParameterInt (0, 1) | 60 | 60 | 0 |  |
| 9 | 60: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 60: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 61: ParameterInt (0, 1) | 61 | 61 | 0 |  |
| 9 | 61: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 61: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 62: ParameterInt (0, 1) | 62 | 62 | 0 |  |
| 9 | 62: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 62: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 63: ParameterInt (0, 1) | 63 | 63 | 0 |  |
| 9 | 63: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 63: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 64: ParameterInt (0, 1) | 64 | 64 | 0 |  |
| 9 | 64: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 64: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 65: ParameterInt (0, 1) | 65 | 65 | 0 |  |
| 9 | 65: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 65: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 66: ParameterInt (0, 1) | 66 | 66 | 0 |  |
| 9 | 66: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 66: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 67: ParameterInt (0, 1) | 67 | 67 | 0 |  |
| 9 | 67: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 67: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 68: ParameterInt (0, 1) | 68 | 68 | 0 |  |
| 9 | 68: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 68: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 69: ParameterInt (0, 1) | 69 | 69 | 0 |  |
| 9 | 69: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 69: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 70: ParameterInt (0, 1) | 70 | 70 | 0 |  |
| 9 | 70: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 70: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 71: ParameterInt (0, 1) | 71 | 71 | 0 |  |
| 9 | 71: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 71: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 72: ParameterInt (0, 1) | 72 | 72 | 0 |  |
| 9 | 72: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 72: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 73: ParameterInt (0, 1) | 73 | 73 | 0 |  |
| 9 | 73: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 73: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 74: ParameterInt (0, 1) | 74 | 74 | 0 |  |
| 9 | 74: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 74: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 75: ParameterInt (0, 1) | 75 | 75 | 0 |  |
| 9 | 75: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 75: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 76: ParameterInt (0, 1) | 76 | 76 | 0 |  |
| 9 | 76: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 76: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 77: ParameterInt (0, 1) | 77 | 77 | 0 |  |
| 9 | 77: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 77: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 78: ParameterInt (0, 1) | 78 | 78 | 0 |  |
| 9 | 78: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 78: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 79: ParameterInt (0, 1) | 79 | 79 | 0 |  |
| 9 | 79: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 79: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 80: ParameterInt (0, 1) | 80 | 80 | 0 |  |
| 9 | 80: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 80: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 81: ParameterInt (0, 1) | 81 | 81 | 0 |  |
| 9 | 81: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 81: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 82: ParameterInt (0, 1) | 82 | 82 | 0 |  |
| 9 | 82: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 82: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 83: ParameterInt (0, 1) | 83 | 83 | 0 |  |
| 9 | 83: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 83: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 84: ParameterInt (0, 1) | 84 | 84 | 0 |  |
| 9 | 84: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 84: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 85: ParameterInt (0, 1) | 85 | 85 | 0 |  |
| 9 | 85: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 85: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 86: ParameterInt (0, 1) | 86 | 86 | 0 |  |
| 9 | 86: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 86: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 87: ParameterInt (0, 1) | 87 | 87 | 0 |  |
| 9 | 87: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 87: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 88: ParameterInt (0, 1) | 88 | 88 | 0 |  |
| 9 | 88: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 88: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 89: ParameterInt (0, 1) | 89 | 89 | 0 |  |
| 9 | 89: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 89: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 90: ParameterInt (0, 1) | 90 | 90 | 0 |  |
| 9 | 90: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 90: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 91: ParameterInt (0, 1) | 91 | 91 | 0 |  |
| 9 | 91: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 91: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 92: ParameterInt (0, 1) | 92 | 92 | 0 |  |
| 9 | 92: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 92: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 93: ParameterInt (0, 1) | 93 | 93 | 0 |  |
| 9 | 93: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 93: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 94: ParameterInt (0, 1) | 94 | 94 | 0 |  |
| 9 | 94: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 94: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 95: ParameterInt (0, 1) | 95 | 95 | 0 |  |
| 9 | 95: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 95: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 96: ParameterInt (0, 1) | 96 | 96 | 0 |  |
| 9 | 96: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 96: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 97: ParameterInt (0, 1) | 97 | 97 | 0 |  |
| 9 | 97: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 97: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 98: ParameterInt (0, 1) | 98 | 98 | 0 |  |
| 9 | 98: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 98: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 99: ParameterInt (0, 1) | 99 | 99 | 0 |  |
| 9 | 99: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 99: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 100: ParameterInt (0, 1) | 100 | 100 | 0 |  |
| 9 | 100: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 100: ParameterInt (0, 1) | 1 | 1 | 0 |  |
| 9 | 101: ParameterInt (0, 1) | 101 | 101 | 0 |  |
| 9 | 101: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 101: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 102: ParameterInt (0, 1) | 102 | 102 | 0 |  |
| 9 | 102: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 102: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 103: ParameterInt (0, 1) | 103 | 103 | 0 |  |
| 9 | 103: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 103: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 104: ParameterInt (0, 1) | 104 | 104 | 0 |  |
| 9 | 104: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 104: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 105: ParameterInt (0, 1) | 105 | 105 | 0 |  |
| 9 | 105: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 105: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 106: ParameterInt (0, 1) | 106 | 106 | 0 |  |
| 9 | 106: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 106: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 107: ParameterInt (0, 1) | 107 | 107 | 0 |  |
| 9 | 107: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 107: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 108: ParameterInt (0, 1) | 108 | 108 | 0 |  |
| 9 | 108: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 108: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 109: ParameterInt (0, 1) | 109 | 109 | 0 |  |
| 9 | 109: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 109: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 110: ParameterInt (0, 1) | 110 | 110 | 0 |  |
| 9 | 110: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 110: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 111: ParameterInt (0, 1) | 111 | 111 | 0 |  |
| 9 | 111: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 111: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 112: ParameterInt (0, 1) | 112 | 112 | 0 |  |
| 9 | 112: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 112: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 113: ParameterInt (0, 1) | 113 | 113 | 0 |  |
| 9 | 113: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 113: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 114: ParameterInt (0, 1) | 114 | 114 | 0 |  |
| 9 | 114: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 114: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 115: ParameterInt (0, 1) | 115 | 115 | 0 |  |
| 9 | 115: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 115: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 116: ParameterInt (0, 1) | 116 | 116 | 0 |  |
| 9 | 116: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 116: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 117: ParameterInt (0, 1) | 117 | 117 | 0 |  |
| 9 | 117: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 117: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 118: ParameterInt (0, 1) | 118 | 118 | 0 |  |
| 9 | 118: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 118: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 119: ParameterInt (0, 1) | 119 | 119 | 0 |  |
| 9 | 119: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 119: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 120: ParameterInt (0, 1) | 120 | 120 | 0 |  |
| 9 | 120: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 120: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 121: ParameterInt (0, 1) | 121 | 121 | 0 |  |
| 9 | 121: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 121: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 122: ParameterInt (0, 1) | 122 | 122 | 0 |  |
| 9 | 122: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 122: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 123: ParameterInt (0, 1) | 123 | 123 | 0 |  |
| 9 | 123: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 123: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 124: ParameterInt (0, 1) | 124 | 124 | 0 |  |
| 9 | 124: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 124: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 125: ParameterInt (0, 1) | 125 | 125 | 0 |  |
| 9 | 125: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 125: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 126: ParameterInt (0, 1) | 126 | 126 | 0 |  |
| 9 | 126: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 126: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 127: ParameterInt (0, 1) | 127 | 127 | 0 |  |
| 9 | 127: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 127: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 128: ParameterInt (0, 1) | 128 | 128 | 0 |  |
| 9 | 128: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 128: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 129: ParameterInt (0, 1) | 129 | 129 | 0 |  |
| 9 | 129: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 129: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 130: ParameterInt (0, 1) | 130 | 130 | 0 |  |
| 9 | 130: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 130: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 131: ParameterInt (0, 1) | 131 | 131 | 0 |  |
| 9 | 131: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 131: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 132: ParameterInt (0, 1) | 132 | 132 | 0 |  |
| 9 | 132: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 132: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 133: ParameterInt (0, 1) | 133 | 133 | 0 |  |
| 9 | 133: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 133: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 134: ParameterInt (0, 1) | 134 | 134 | 0 |  |
| 9 | 134: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 134: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 135: ParameterInt (0, 1) | 135 | 135 | 0 |  |
| 9 | 135: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 135: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 136: ParameterInt (0, 1) | 136 | 136 | 0 |  |
| 9 | 136: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 136: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 137: ParameterInt (0, 1) | 137 | 137 | 0 |  |
| 9 | 137: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 137: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 138: ParameterInt (0, 1) | 138 | 138 | 0 |  |
| 9 | 138: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 138: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 139: ParameterInt (0, 1) | 139 | 139 | 0 |  |
| 9 | 139: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 139: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 140: ParameterInt (0, 1) | 140 | 140 | 0 |  |
| 9 | 140: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 140: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 141: ParameterInt (0, 1) | 141 | 141 | 0 |  |
| 9 | 141: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 141: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 142: ParameterInt (0, 1) | 142 | 142 | 0 |  |
| 9 | 142: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 142: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 143: ParameterInt (0, 1) | 143 | 143 | 0 |  |
| 9 | 143: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 143: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 144: ParameterInt (0, 1) | 144 | 144 | 0 |  |
| 9 | 144: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 144: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 145: ParameterInt (0, 1) | 145 | 145 | 0 |  |
| 9 | 145: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 145: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 146: ParameterInt (0, 1) | 146 | 146 | 0 |  |
| 9 | 146: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 146: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 147: ParameterInt (0, 1) | 147 | 147 | 0 |  |
| 9 | 147: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 147: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 148: ParameterInt (0, 1) | 148 | 148 | 0 |  |
| 9 | 148: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 148: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 149: ParameterInt (0, 1) | 149 | 149 | 0 |  |
| 9 | 149: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 149: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 150: ParameterInt (0, 1) | 150 | 150 | 0 |  |
| 9 | 150: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 150: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 151: ParameterInt (0, 1) | 151 | 151 | 0 |  |
| 9 | 151: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 151: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 152: ParameterInt (0, 1) | 152 | 152 | 0 |  |
| 9 | 152: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 152: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 153: ParameterInt (0, 1) | 153 | 153 | 0 |  |
| 9 | 153: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 153: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 154: ParameterInt (0, 1) | 154 | 154 | 0 |  |
| 9 | 154: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 154: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 155: ParameterInt (0, 1) | 155 | 155 | 0 |  |
| 9 | 155: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 155: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 156: ParameterInt (0, 1) | 156 | 156 | 0 |  |
| 9 | 156: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 156: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 157: ParameterInt (0, 1) | 157 | 157 | 0 |  |
| 9 | 157: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 157: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 158: ParameterInt (0, 1) | 158 | 158 | 0 |  |
| 9 | 158: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 158: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 159: ParameterInt (0, 1) | 159 | 159 | 0 |  |
| 9 | 159: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 159: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 160: ParameterInt (0, 1) | 160 | 160 | 0 |  |
| 9 | 160: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 160: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | 161: ParameterInt (0, 1) | 161 | 161 | 0 |  |
| 9 | 161: ParameterInt (0, 1) | 0 | 0 | 0 | near-zero |
| 9 | 161: ParameterInt (0, 1) | 1 | 1 | 0 | near-zero |
| 9 | Maximum Run Count: 780 | 780 | 780 | 0 |  |
| 11 | Best Result: 0.8665353911783293 | 0.866535 | 0.850119 | 1.89% |  |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 0 | 0 | 0 | near-zero |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 0 | 1 | 1,111% | near-zero |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 1 | 1 | 0 |  |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 1 | 1 | 0 |  |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 2 | 2 | 0 |  |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 1 | 1 | 0 |  |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 3 | 3 | 0 |  |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 0 | 1 | 1,111% | near-zero |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 4 | 4 | 0 |  |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 0 | 1 | 1,111% | near-zero |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 5 | 5 | 0 |  |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 1 | 1 | 0 |  |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 6 | 6 | 0 |  |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 0 | 0 | 0 | near-zero |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 7 | 7 | 0 |  |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 1 | 1 | 0 |  |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 8 | 8 | 0 |  |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 1 | 1 | 0 |  |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 9 | 9 | 0 |  |
| 11 | Best Parameters: 0: 0 1: 1 2: 1 3: 0 4: 0 5: 1 6: 0 7: 1 8: 1 9: 1 | 1 | 1 | 0 |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 10 | 10 | 0 |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 1 | 0 | 100% |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 11 | 11 | 0 |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 1 | 0 | 100% |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 12 | 12 | 0 |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 0 | 0 | 0 | near-zero |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 13 | 13 | 0 |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 0 | 1 | 500% | near-zero |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 14 | 14 | 0 |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 0 | 1 | 500% | near-zero |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 15 | 15 | 0 |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 1 | 1 | 0 |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 16 | 16 | 0 |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 1 | 0 | 100% |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 17 | 17 | 0 |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 0 | 1 | 500% | near-zero |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 18 | 18 | 0 |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 1 | 0 | 100% |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 19 | 19 | 0 |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 0 | 0 | 0 | near-zero |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 20 | 20 | 0 |  |
| 11 | 10: 1 11: 1 12: 0 13: 0 14: 0 15: 1 16: 1 17: 0 18: 1 19: 0 20: 1 | 1 | 0 | 100% |  |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 21 | 21 | 0 |  |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 1 | 1 | 0 |  |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 22 | 22 | 0 |  |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 0 | 0 | 0 | near-zero |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 23 | 23 | 0 |  |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 0 | 0 | 0 | near-zero |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 24 | 24 | 0 |  |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 1 | 1 | 0 |  |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 25 | 25 | 0 |  |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 0 | 1 | 323% | near-zero |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 26 | 26 | 0 |  |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 0 | 0 | 0 | near-zero |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 27 | 27 | 0 |  |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 0 | 1 | 323% | near-zero |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 28 | 28 | 0 |  |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 0 | 0 | 0 | near-zero |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 29 | 29 | 0 |  |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 0 | 0 | 0 | near-zero |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 30 | 30 | 0 |  |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 0 | 1 | 323% | near-zero |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 31 | 31 | 0 |  |
| 11 | 21: 1 22: 0 23: 0 24: 1 25: 0 26: 0 27: 0 28: 0 29: 0 30: 0 31: 1 | 1 | 1 | 0 |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 32 | 32 | 0 |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 0 | 0 | 0 | near-zero |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 33 | 33 | 0 |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 0 | 1 | 238% | near-zero |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 34 | 34 | 0 |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 1 | 1 | 0 |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 35 | 35 | 0 |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 1 | 0 | 100% |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 36 | 36 | 0 |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 0 | 0 | 0 | near-zero |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 37 | 37 | 0 |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 0 | 0 | 0 | near-zero |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 38 | 38 | 0 |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 1 | 1 | 0 |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 39 | 39 | 0 |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 0 | 0 | 0 | near-zero |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 40 | 40 | 0 |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 1 | 0 | 100% |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 41 | 41 | 0 |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 1 | 1 | 0 |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 42 | 42 | 0 |  |
| 11 | 32: 0 33: 0 34: 1 35: 1 36: 0 37: 0 38: 1 39: 0 40: 1 41: 1 42: 0 | 0 | 0 | 0 | near-zero |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 43 | 43 | 0 |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 0 | 1 | 189% | near-zero |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 44 | 44 | 0 |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 1 | 1 | 0 |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 45 | 45 | 0 |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 1 | 0 | 100% |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 46 | 46 | 0 |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 1 | 1 | 0 |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 47 | 47 | 0 |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 0 | 0 | 0 | near-zero |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 48 | 48 | 0 |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 1 | 0 | 100% |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 49 | 49 | 0 |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 1 | 0 | 100% |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 50 | 50 | 0 |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 0 | 0 | 0 | near-zero |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 51 | 51 | 0 |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 0 | 0 | 0 | near-zero |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 52 | 52 | 0 |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 0 | 1 | 189% | near-zero |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 53 | 53 | 0 |  |
| 11 | 43: 0 44: 1 45: 1 46: 1 47: 0 48: 1 49: 1 50: 0 51: 0 52: 0 53: 1 | 1 | 0 | 100% |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 54 | 54 | 0 |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 1 | 1 | 0 |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 55 | 55 | 0 |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 1 | 0 | 100% |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 56 | 56 | 0 |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 0 | 1 | 156% | near-zero |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 57 | 57 | 0 |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 1 | 0 | 100% |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 58 | 58 | 0 |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 1 | 0 | 100% |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 59 | 59 | 0 |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 1 | 1 | 0 |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 60 | 60 | 0 |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 0 | 0 | 0 | near-zero |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 61 | 61 | 0 |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 0 | 0 | 0 | near-zero |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 62 | 62 | 0 |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 1 | 1 | 0 |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 63 | 63 | 0 |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 0 | 1 | 156% | near-zero |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 64 | 64 | 0 |  |
| 11 | 54: 1 55: 1 56: 0 57: 1 58: 1 59: 1 60: 0 61: 0 62: 1 63: 0 64: 1 | 1 | 0 | 100% |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 65 | 65 | 0 |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 1 | 1 | 0 |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 66 | 66 | 0 |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 1 | 0 | 100% |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 67 | 67 | 0 |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 0 | 0 | 0 | near-zero |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 68 | 68 | 0 |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 0 | 1 | 133% | near-zero |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 69 | 69 | 0 |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 1 | 0 | 100% |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 70 | 70 | 0 |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 0 | 1 | 133% | near-zero |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 71 | 71 | 0 |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 0 | 1 | 133% | near-zero |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 72 | 72 | 0 |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 0 | 1 | 133% | near-zero |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 73 | 73 | 0 |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 0 | 0 | 0 | near-zero |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 74 | 74 | 0 |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 1 | 0 | 100% |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 75 | 75 | 0 |  |
| 11 | 65: 1 66: 1 67: 0 68: 0 69: 1 70: 0 71: 0 72: 0 73: 0 74: 1 75: 0 | 0 | 1 | 133% | near-zero |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 76 | 76 | 0 |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 1 | 0 | 100% |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 77 | 77 | 0 |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 1 | 0 | 100% |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 78 | 78 | 0 |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 1 | 1 | 0 |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 79 | 79 | 0 |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 1 | 1 | 0 |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 80 | 80 | 0 |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 0 | 0 | 0 | near-zero |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 81 | 81 | 0 |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 1 | 1 | 0 |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 82 | 82 | 0 |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 0 | 1 | 116% | near-zero |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 83 | 83 | 0 |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 1 | 0 | 100% |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 84 | 84 | 0 |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 0 | 1 | 116% | near-zero |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 85 | 85 | 0 |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 0 | 0 | 0 | near-zero |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 86 | 86 | 0 |  |
| 11 | 76: 1 77: 1 78: 1 79: 1 80: 0 81: 1 82: 0 83: 1 84: 0 85: 0 86: 1 | 1 | 0 | 100% |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 87 | 87 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 1 | 1 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 88 | 88 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 0 | 0 | 0 | near-zero |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 89 | 89 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 0 | 1 | 103% | near-zero |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 90 | 90 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 1 | 1 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 91 | 91 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 1 | 0 | 100% |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 92 | 92 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 1 | 1 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 93 | 93 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 1 | 1 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 94 | 94 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 1 | 1 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 95 | 95 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 0 | 1 | 103% | near-zero |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 96 | 96 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 1 | 1 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 97 | 97 | 0 |  |
| 11 | 87: 1 88: 0 89: 0 90: 1 91: 1 92: 1 93: 1 94: 1 95: 0 96: 1 97: 1 | 1 | 1 | 0 |  |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 98 | 108 | 10.20% |  |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 0 | 1 | 93.46% | near-zero |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 99 | 109 | 10.10% |  |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 0 | 1 | 93.46% | near-zero |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 100 | 110 | 10.00% |  |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 0 | 1 | 93.46% | near-zero |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 101 | 111 | 9.90% |  |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 1 | 1 | 0 | near-zero |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 102 | 112 | 9.80% |  |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 0 | 1 | 93.46% | near-zero |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 103 | 113 | 9.71% |  |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 1 | 0 | 93.46% | near-zero |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 104 | 114 | 9.62% |  |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 1 | 0 | 93.46% | near-zero |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 105 | 115 | 9.52% |  |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 0 | 1 | 93.46% | near-zero |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 106 | 116 | 9.43% |  |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 0 | 1 | 93.46% | near-zero |
| 11 | 98: 0 99: 0 100: 0 101: 1 102: 0 103: 1 104: 1 105: 0 106: 0 107: | 107 | 117 | 9.35% |  |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 1 | 1 | 0 | near-zero |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 108 | 118 | 9.26% |  |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 1 | 1 | 0 | near-zero |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 109 | 119 | 9.17% |  |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 0 | 0 | 0 | near-zero |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 110 | 120 | 9.09% |  |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 1 | 0 | 86.21% | near-zero |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 111 | 121 | 9.01% |  |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 1 | 1 | 0 | near-zero |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 112 | 122 | 8.93% |  |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 0 | 0 | 0 | near-zero |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 113 | 123 | 8.85% |  |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 0 | 1 | 86.21% | near-zero |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 114 | 124 | 8.77% |  |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 0 | 0 | 0 | near-zero |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 115 | 125 | 8.70% |  |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 1 | 1 | 0 | near-zero |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 116 | 126 | 8.62% |  |
| 11 | 1 108: 1 109: 0 110: 1 111: 1 112: 0 113: 0 114: 0 115: 1 116: 0 | 0 | 1 | 86.21% | near-zero |
| 12 | Fitness value of the best solution = 0.867 | 0.867 | 0.85 | 1.96% |  |

## Figures

Fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. This is a screening number, not a physical error: an autoscaled colour bar, antialiasing, or a shifted tick label all register. Figures were also compared by eye, and the verdict reflects that.

| Cell | Figure | Pixel difference |
|---:|---:|---:|
| 5 | 0 | 4.20% |
| 12 | 0 | 5.03% |
| 12 | 1 | 6.49% |
| 14 | 0 | 9.87% |
| 15 | 0 | 3.65% |
