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

## `mrct_case_study/`

The case example of the mrct-consistency report, `results/s2_example/` of
`../mrct-consistency` (read at `2341f96`; last changed at `85377ec`,
2026-10-06), copied byte for byte (whole files) for
`test-mrct-case-study.R` (TASK_fs-mrct-article-r2_2026-11-01, section 2):
replicate 5 of the S2 cell, its four sets, its H-hat and its 2 x 2 tables.
The same three files ship with the article
(`vignettes/articles/mrct-consistency/`) and are identical to these.

| Fixture | Source (mrct-consistency) | SHA-256 (12) |
|---|---|---|
| `s2_example_four_sets.csv` | `results/s2_example/s2_example_four_sets.csv` | `f101a1f7b9d9` |
| `s2_example_trials.csv` | `results/s2_example/s2_example_trials.csv` | `4d57e6df5a60` |
| `s2_example_definitions.md` | `results/s2_example/s2_example_definitions.md` | `aa6b37b93ac0` |
