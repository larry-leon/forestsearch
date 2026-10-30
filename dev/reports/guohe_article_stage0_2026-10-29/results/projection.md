# Projection (pure arithmetic from measured single-core seconds)

Per-replicate seconds (family + Guo-He + MR, means over 3 replicates): rung 1 72.2, rung 2 19.1, rung 3 8.5.
Truth curve beta2 = 0.3, c_step = 1: n_big 2e5 5.9 s; n_big 5e4 1.2 s. Added once for the single beta2 > 0 cell.

Cells {beta2 = 0, 0.3}; 2 cells x n_rep replicates. Truth at n_big = 2e5 included.

| n_rep | rung | 1 core | 2 workers (assumed linear scaling) | 1 core, truth at 5e4 |
|---:|---|---:|---:|---:|
| 50 | 1 | 2.0 h (121 min) | 1.0 h (60 min) | 2.0 h (120 min) |
| 50 | 2 | 0.5 h (32 min) | 0.3 h (16 min) | 0.5 h (32 min) |
| 50 | 3 | 0.2 h (14 min) | 0.1 h (7 min) | 0.2 h (14 min) |
| 100 | 1 | 4.0 h (241 min) | 2.0 h (121 min) | 4.0 h (241 min) |
| 100 | 2 | 1.1 h (64 min) | 0.5 h (32 min) | 1.1 h (64 min) |
| 100 | 3 | 0.5 h (28 min) | 0.2 h (14 min) | 0.5 h (28 min) |
| 200 | 1 | 8.0 h (482 min) | 4.0 h (241 min) | 8.0 h (482 min) |
| 200 | 2 | 2.1 h (127 min) | 1.1 h (64 min) | 2.1 h (127 min) |
| 200 | 3 | 0.9 h (57 min) | 0.5 h (28 min) | 0.9 h (57 min) |

2-worker column: ASSUMPTION -- replicate-level parallelism scaling linearly; the truth curve is not divided (one serial call).

Reference, not an article size: the manuscript's shape, 6 cells x 2000 replicates at rung 1, single core:
240.8 h (14450 min) (truth curves excluded; the manuscript's caches were n_big = 2e6).
