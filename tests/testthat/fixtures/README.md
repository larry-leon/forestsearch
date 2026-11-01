# Test fixtures

## `grf_disjunction.rds`

Used by `test-betaHhat-contract.R` (not described here; pre-existing).

## `mrct_analytic/`

Records of the mrct-consistency repository (`../mrct-consistency`, read at
`4710058`), copied byte for byte (rows selected by line, header kept) for
`test-ahr-analytic.R` and `test-mrct-dgm-args.R`
(TASK_fs-mrct-analytic-layer_2026-11-01, section 5). Each is named with the
mrct-consistency commit that last changed its source file.

| Fixture | Source (mrct-consistency) | Rows kept | Source commit |
|---|---|---|---|
| `predictions_S2_rmlog4.csv` | `results/ident_rlog4_full/analytic/predictions.csv` | header + the 61 rows of cell `S2_rmlog4` | `568a9bd` (2026-10-05) |
| `grid_b055_h2_n1000_rmlog2.csv` | `results/analytic_maps/grid.csv` | header + line 526 (r = -log 2, b = 0.55, h = 2, n = 1,000) | `6c978ee` (2026-10-05) |
| `roots.csv` | `results/regional_entry/roots.csv` | all (the two roots and two reference rows) | `3deb2e6` (2026-11-01) |
| `build_S2_rmlog4.csv` | `results/ident_rlog4_full/summary/build.csv` | header + the row of `S2_rmlog4` | `fca4c1e` (2026-10-05) |
