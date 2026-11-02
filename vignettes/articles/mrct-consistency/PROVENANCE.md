# Provenance of the article's data files

Records of the mrct-consistency repository, read at commit `2341f96`
(2026-11-01, the commit that closed `TASK_mr-field-campaign_2026-11-01`),
for the article `vignettes/articles/mrct-consistency.qmd`
(TASK_fs-mrct-article_2026-11-01). The article reads only the files of this
directory; it never reads mrct-consistency.

Every file except two is a byte-for-byte copy of its source (checked with
`cmp` when copied). The source commit is the last mrct-consistency commit that
changed the source path; the hash is the first 12 hex digits of the shipped
file's SHA-256.

| File | Source (mrct-consistency) | Source commit | SHA-256 (12) |
|---|---|---|---|
| `build_full_summary.csv` | `results/ident_rlog4_full/summary/build.csv` | fca4c1e (2026-10-05) | `0203e941848f` |
| `core_n1000.csv` | `results/ident_rlog4_S1/summary/core_n1000.csv` | 671e44b (2026-10-05) | `79e95b9024ec` |
| `core_full_numeric.csv` | `results/ident_rlog4_full/summary/core.csv` | fca4c1e (2026-10-05) | `9da02a15b1b9` |
| `core_S1_numeric.csv` | `results/ident_rlog4_S1/summary/core_S1.csv` | 671e44b (2026-10-05) | `91a933d95467` |
| `entry_cells_core.csv` | `results/regional_entry_cells/summary/core.csv` | 41c6bfa (2026-11-01) | `d99466c3a028` |
| `entry_cells_core_numeric.csv` | `results/regional_entry_cells/summary/core_numeric.csv` | 41c6bfa (2026-11-01) | `11c8aea3587c` |
| `predictions_full.csv` | `results/ident_rlog4_full/analytic/predictions.csv` | 568a9bd (2026-10-05) | `3140741570c9` |
| `predictions_S1.csv` | `results/ident_rlog4_S1/analytic/predictions.csv` | 3134494 (2026-10-05) | `c18eb7b15483` |
| `predictions_n700.csv` | `results/ident_rlog4_n700/analytic/predictions.csv` | 3716174 (2026-10-05) | `0bb42d7933d5` |
| `entry_roots.csv` | `results/regional_entry/roots.csv` | 3deb2e6 (2026-11-01) | `f85e4716e857` |
| `mr_field_core.csv` | `results/mr_field_campaign/summary/mr_field_core.csv` | 23420aa (2026-11-01) | `3be609c0f002` |
| `s2_example_four_sets.csv` | `results/s2_example/s2_example_four_sets.csv` | 85377ec (2026-10-06) | `f101a1f7b9d9` |
| `s2_example_trials.csv` | `results/s2_example/s2_example_trials.csv` | 85377ec (2026-10-06) | `4d57e6df5a60` |
| `s2_example_definitions.md` | `results/s2_example/s2_example_definitions.md` | 85377ec (2026-10-06) | `aa6b37b93ac0` |
| `timing_full.csv` | `results/ident_rlog4_full/summary/timing.csv` | fca4c1e (2026-10-05) | `b444edcdd670` |
| `timing_S1.csv` | `results/ident_rlog4_S1/summary/timing.csv` | 671e44b (2026-10-05) | `fbf51c8b7d8e` |
| `timing_entry_cells.csv` | `results/regional_entry_cells/summary/timing.csv` | 41c6bfa (2026-11-01) | `607a93a30a62` |
| `timing_mr_field_N0.csv` | `results/mr_field_campaign/N0/timing.csv` | 23420aa (2026-11-01) | `27e1ddd4d6fa` |
| `timing_mr_field_N1.csv` | `results/mr_field_campaign/N1/timing.csv` | 23420aa (2026-11-01) | `784b2e4b6c00` |
| `timing_mr_field_S1.csv` | `results/mr_field_campaign/S1/timing.csv` | 23420aa (2026-11-01) | `da57707e9c81` |
| `timing_mr_field_S2.csv` | `results/mr_field_campaign/S2/timing.csv` | 23420aa (2026-11-01) | `d5d6b039b434` |
| `timing_mr_field_S3.csv` | `results/mr_field_campaign/S3/timing.csv` | 23420aa (2026-11-01) | `ef9c60583699` |
| `timing_mr_field_S4.csv` | `results/mr_field_campaign/S4/timing.csv` | 23420aa (2026-11-01) | `48719143ba97` |
| `timing_mr_field_S2_D1.csv` | `results/mr_field_campaign/S2_D1/timing.csv` | 23420aa (2026-11-01) | `e262039aa3fd` |
| `timing_mr_field_S2_D2.csv` | `results/mr_field_campaign/S2_D2/timing.csv` | 23420aa (2026-11-01) | `5ce59a839927` |
| `mr_field_S2_rep5.csv` | `results/mr_field_campaign/S2/mr_field.csv` (header + line 2, replicate 5) | 23420aa (2026-11-01) | `8cf7a9d3cd7c` |
| `latent_R.csv` | written from `build_full("S2")$mrct$latent_R` (`sims/ident_rlog4_full_design.R`, run at `2341f96` with forestsearch 0.4.1.9000 installed from `5db25a27`) | `2341f96` (2026-11-01) | `2f959b361809` |

The two files that are not whole copies:

- `latent_R.csv`: the 12 x 12 latent correlation of the copula expansion of
  the S2 build of record, one row per covariate (`var` gives the row name),
  each entry written with `sprintf("%.17g")` so that it reads back exactly.
  It is identical (to 0) to the correlation that a fresh latentcor fit gives
  in `tests/testthat/helper-mrct-analytic.R`, and a build with it supplied as
  `inject_mrct_structure(latent_R = )` gives a `df_super` identical to the
  fitted build. The correlation does not depend on the cell, so one file
  serves all eight.
- `mr_field_S2_rep5.csv`: the header and the row of replicate 5 of
  `results/mr_field_campaign/S2/mr_field.csv` (lines 1-2), the MR field
  campaign's record of the article's case example; the article checks its
  own MR call against it.

What each file is used for:

| File | Used for |
|---|---|
| `latent_R.csv` | the builds of every cell, in the collapsed "Reproducing this report's build exactly" chunk (the reader's path fits the correlation live) |
| `build_full_summary.csv` | the build check (the `S2_rmlog4` row); the cells' `h` and `r` |
| `predictions_full.csv`, `predictions_S1.csv`, `predictions_n700.csv` | checks of the live analytic layer (truths and design probabilities at n = 1,000 and 700) |
| `entry_roots.csv` | checks of the three roots (r at no delay; D1, D2) |
| `s2_example_four_sets.csv`, `s2_example_trials.csv`, `s2_example_definitions.md` | the replicate-5 reproduction check (four HRs, H-hat, the 2 x 2 tables); replicates 1-5 for the live demonstration's check |

The three `s2_example_*` files also ship, byte-identical, as the fixture of
`tests/testthat/test-mrct-case-study.R` (`tests/testthat/fixtures/mrct_case_study/`,
TASK_fs-mrct-article-r2_2026-11-01), which checks that `mrct_select_case()` and
`mrct_case_study()` reproduce replicate 5.
| `mr_field_S2_rep5.csv` | the check of the replicate-5 MR call |
| `core_n1000.csv` | the core OC table at n = 1,000 (displayed) |
| `core_full_numeric.csv`, `core_S1_numeric.csv` | the recount of the core table's counts; observed trigger and oracle-joint frequencies |
| `entry_cells_core.csv`, `entry_cells_core_numeric.csv` | the later-enrolling cells' table (displayed, recounted), observed frequencies |
| `mr_field_core.csv` | the MR field rows (displayed) |
| `timing_*.csv` | the cells' specifications (seeds, h, r, a, replicates) and the wall-clock of the committed runs |
