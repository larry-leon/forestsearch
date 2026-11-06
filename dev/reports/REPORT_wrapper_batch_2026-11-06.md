# REPORT — wrapper batch: P3, P4, P6, P7, P8 (2026-11-06)

Task: `dev/tasks/TASK_wrapper_batch_2026-11-06.md`.
**Outcome: every gate passes, with 0 errors, 0 warnings and 0 notes from `check()`, the same as the reused baseline.**
The work is on branch `fs-wrappers-batch`, unmerged and unpushed, for your review. The final tree is `efc785af`; this report's commit sits on top of it.

## Baseline reuse (no baseline run, per the testing cadence)

`git diff 84afa853..7e955ba7 -- R tests man vignettes inst DESCRIPTION NAMESPACE .Rbuildignore` is **empty**. The files changed in that range are only these:

```
dev/reports/REPORT_mrct_trial_analysis_2026-11-06.md
dev/tasks/TASK_mrct_trial_analysis_2026-11-06.md
dev/tasks/TASK_wrapper_batch_2026-11-06.md
```

So `REPORT_mrct_trial_analysis_2026-11-06.md`'s final figures are the baseline: `check()` 0/0/0, status OK; suite 0 failures, 3 skips, 38 warnings, 7,962 expectations.

## Commits

```
28681558 feat(B1): mrct_core_summary(), the core rows of an mrct_region_sims() campaign
963512f2 feat(B2): ahr_truths_table() and ahr_design_table(), the analytic layer over cell
94bc34e3 feat(B3): mr_method_comparison() and mr_cost_table()
24c3fb05 fix(B1): mrct_core_summary()'s example reaches the trigger; a real campaign in it
3bb782f9 feat(B4): fs_declaration_table(), the declaration calibration as one tidy table
02d1bfc0 fix(B4): fs_declaration_table()'s example runs in a plain check
62d21fc4 feat(B5): mrct_record(), fs_in_subgroup() and mrct_build_dgm()
3e304126 fix(B1): mrct_core_summary() reads its columns directly
efc785af fix(B1): the en dash in mrct_core_summary() as an escape
```

**The doc asks for one commit per item, and this is nine.** There are the five item commits, plus four labelled fixes I made afterwards: three to B1 (its example, then two small code edits) and one to B4's example. Each fix is labelled with its item, so commits still map to items. The history was not rewritten; that is your call.

## Testing cadence: what ran

| Cadence rule | What happened |
|---|---|
| 1. No baseline run | None. The git diff above, and the P1+P2 figures reused |
| 2. Inner loop: new test files and each item's fast fixture only | Followed. No full suite or `check()` until the final stage. Lightweight extras: the new functions' roxygen examples, and `codetools` on the new functions only |
| 3. Each gate's fixture once | Each of B1-B5 ran once |
| 4. Full suite, `check()`, `check_pkgdown()` once, at the final stage | Once each, at `efc785af` |

## Gate 0: Stage 0 field notes (PASS)

Every field the doc names was verified, or its absence is reported here. Lines are at `7e955ba7`.

- **B1 (P3).**
  - `mrct_region_sims()` (`R/mrct_simulation.R:260`) returns a plain `data.table` (`:706`, `:742`), with **no class of its own**; the gate printed `data.table data.frame`. By D-B1.1 that means an exported function, `mrct_core_summary(res)`, not a `summary()` method.
  - The `"estimators"` attribute is set at `:721`, under `record = "full"` only. Its keys are `itt_all`, `itt_nonregion`, `itt_region`, `region_Ghat`, and, with `flag_harm`, `region_Gstar` and `Hstar_nonregion`.
  - The fields used are `sim`, `estimator`, `n`, `hr`, `lo`, `hi`. The per-trial columns are `sim`, `any_found` and the `TP_/FP_/FN_/TN_{nonregion,region}` counts.
- **B2 (P4).**
  - `ahr_population()` is at `R/ahr_analytic.R:179`. `ahr_limit()`'s set forms are named strings or a logical vector of length `pop$N` (`:492`), and its fields `hr`, `marghr_inf`, `ahr`, `cde` are at `:511`. `ahr_design()`'s return is at `:595`.
  - `$prob` holds exactly the five named probabilities (`:662-665`). **P(G* < 1 | trigger) is not among them**, so the table computes it as the ratio, as the article does.
- **B3 (P6).**
  - The article reads `H_estimates$H2` and `Hc_estimates$H2` (`R/bootstrap_dofuture_main.R:93-94`, `:574`). `H2` is stored on the fit's own `est.scale` (`:559`), so `est.scale = "1/hr"` inverts only when needed.
  - **The object stores wall time**: `attr(boot$results, "timing")$total_minutes` (`R/bootstrap_analysis_dofuture.R:1004`; read at `R/bootstrap_dofuture_main.R:853`). So `boot_seconds = NULL` reads it. It covers the replicate loop, not the whole call: 3.7 s against 3.9 s measured in the gate.
- **B4 (P7).** `fs_declaration_calibration()` is at `R/fs_declaration_calibration.R:514`. It returns `c0$table` (with `kappa_hat`, `fw_size`), `c0$source`, `kappa_hat`, `fw_size`, `T_hat` and `admitted_calibrated`. `mr_inference$declaration_field` is built at `R/fs_mr_inference.R:1116`, with `family_id` at `:1120` and `meta$selected_index` at `:1124`, and needs `keep_declaration_field = TRUE`.
- **B5 (P8).** `get_dfpred(df.predict, sg.harm, version = 1)` is at `R/forestsearch_helpers.R:90`. The records are in `inst/extdata/mrct-consistency/`, with `PROVENANCE.md`; the dev-tree and installed resolution is in the articles' `packages` chunk.

**Core-table SE definitions, confirmed from the record's own script** (`mrct-consistency/sims/124_ident_rlog4_full_tables.R`): binomial on counts for the proportions, `sd / sqrt(n)` over the trials for the Sens and Spec averages.

## B1 (P3): `mrct_core_summary()` (Gate B1: PASS)

A real 40-trial `record = "full"` campaign, with the `demo-summary` chunk transplanted. Every cell, and the unrounded values, are identical.

```
campaign: 40 trials, n = 400, 10 s; class of the result: data.table data.frame

The article's rows beside the function's:
 row                                    demonstration mrct_core_summary
 trial power (ITT ub < 1.0)             22.5          22.5             
 trials reaching the trigger (count)    4 of 40       4 of 40          
 subgroup declared                      100.0         100.0            
 Sens / Spec of H-hat vs H*, non-region 98.7 / 99.1   98.7 / 99.1      
 Sens / Spec of H-hat vs H*, region     99.0 / 99.4   99.0 / 99.4      
 declared, region(G-hat) < 1.0          50.0          50.0             
 region(G*) < 1.0                       50.0          50.0             
 declared, region(G-hat) >= 1.0         50.0          50.0             

every cell identical: TRUE 
unrounded values identical: TRUE 

The whole table:
 quantity          estimate se      count denominator display
 power              22.50    6.6026  9    40          22.5   
 trigger            10.00        NA  4    40          4 of 40
 declared          100.00    0.0000  4     4          100.0  
 sens_nonregion     98.73    0.8840 NA     4          98.7   
 spec_nonregion     99.09    0.6458 NA     4          99.1   
 sens_region        98.96    1.0417 NA     4          99.0   
 spec_region        99.38    0.3591 NA     4          99.4   
 declared_Ghat_lt1  50.00   25.0000  2     4          50.0   
 Gstar_lt1          50.00   25.0000  2     4          50.0   
 declared_Ghat_ge1  50.00   25.0000  2     4          50.0   
 Ghat_share         75.71        NA NA     4          75.7   

GATE B1: PASS
```

Only 4 of the 40 trials reach the trigger in this fixture, which is enough for equality but makes the numbers themselves uninformative. **The gate ran before two later edits to the function:** `with()` became explicit `demo$` columns, and a literal en dash became an escape. Both leave the arithmetic and strings unchanged, and per the cadence the gate fixture was not re-run. The committed tests, which include the transplanted chunk on a hand-built record and a real campaign, ran after those edits, and in the final suite.

## B2 (P4): `ahr_truths_table()`, `ahr_design_table()` (Gate B2: PASS)

The six `r = -log 4` cells on the shipped `latent_R`, with the article's check code transplanted, against the committed records.

```
 cell    h
   N0 0.55
   N1 0.75
   S1 1.00
   S2 2.00
   S3 3.00
   S4 4.00
six builds: 2.8 s
ahr_truths_table(): 24 rows, 4.4 s
truths against predictions_full.csv + predictions_S1.csv: 96 values, near(1e-9) TRUE; largest |difference| 4.9e-15
 cell          set horizon marghr_40 marghr_inf      ahr     cde
   S2          all      40  0.743115   0.710761 0.673490 1.05863
   S2    nonregion      40  0.723115   0.666126 0.671724 1.05571
   S2       region      40  0.868301   0.678806 0.680574 1.09671
   S2 region_Gstar      40  0.559664   0.616606 0.550000 0.55000

ahr_design_table(), twice: 0.2 s
n = 1,000 against predictions_full.csv + predictions_S1.csv: 5 x 6 values below 1e-9; largest |difference| 5.6e-16
n = 700 against predictions_n700.csv + predictions_S1.csv: 5 x 6 values below 1e-9; largest |difference| 5.6e-16
P(G* < 1 | trigger) against the record: largest |difference| 3e-15
 cell p_itt_ub_lt1 p_region_itt_gt1 p_region_Gstar_lt1 p_trigger
   N0       0.9987          0.05482             0.9135   0.05442
   N1       0.9937          0.08584             0.9135   0.08370
   S1       0.9780          0.13139             0.9135   0.12303
   S2       0.8437          0.32612             0.9135   0.24124
   S3       0.7185          0.47314             0.9135   0.28938
   S4       0.6350          0.56818             0.9135   0.30541
 p_trigger_and_Gstar_lt1 p_Gstar_lt1_given_trigger
                 0.01686                    0.3098
                 0.03649                    0.4359
                 0.06827                    0.5549
                 0.18741                    0.7769
                 0.24547                    0.8483
                 0.26841                    0.8788
 row                       S2   
 r1_itt_ub_lt1             0.844
 r2_region_itt_gt1         0.326
 r3_consistency            0.602
 r4_non_consistency        0.241
 r5_oracle_consistency     0.187
 r6_oracle_non_consistency 0.054
 r7_consistency_or_oracle  0.790
the seven S2 rows print the values of the manuscript's Table 2: TRUE

GATE B2: PASS
```

## B3 (P6): `mr_method_comparison()`, `mr_cost_table()` (Gate B3: PASS)

A real `forestsearch()` fit with MR and a real two-replicate bootstrap, with `three_columns()` and the cost chunk transplanted.

```
fit: 8.8 s; H-hat {x1 <= -0.7}; est.scale hr; complement de-biased TRUE
Bootstrap parallel config: Using user-provided settings
System max cores available: 128
Bootstrap will use:  workers with 'sequential' plan
Bootstrap iteration 1/2
Bootstrap iteration 2/2
bootstrap, B = 2: 3.9 s; stored loop time 3.7 s; nb_boots field 2; H2 3.5900, Hc H2 0.5849

The article's three_columns():
 method region n         naive             bootstrap multiplier       
 FS     H      134 (27%) 3.56 (2.45, 5.16) 3.59      3.62 (1.50, 8.74)
 FS     Hc     366 (73%) 0.54 (0.43, 0.67) 0.58      0.52 (0.29, 0.93)

mr_method_comparison():
 method region n         naive             bootstrap multiplier       
 FS     H      134 (27%) 3.56 (2.45, 5.16) 3.59      3.62 (1.50, 8.74)
 FS     Hc     366 (73%) 0.54 (0.43, 0.67) 0.58      0.52 (0.29, 0.93)

identical, cell for cell: TRUE 

The article's cost chunk:
 method mr_seconds bootstrap_B bootstrap_seconds seconds_per_replicate
     FS      0.015           2             3.921                1.9605
 projected_B_1000
             1960

mr_cost_table(boot_seconds = ):
 method mr_seconds bootstrap_B bootstrap_seconds seconds_per_replicate
     FS      0.015           2             3.921                1.9605
 projected_B_1000
             1960
Warning messages:
1: In identical(as.numeric(unname(cost[[cn]])), as.numeric(ours_cost[[cn]])) :
  NAs introduced by coercion
2: In identical(as.numeric(unname(cost[[cn]])), as.numeric(ours_cost[[cn]])) :
  NAs introduced by coercion
               method            mr_seconds           bootstrap_B 
                 TRUE                  TRUE                  TRUE 
    bootstrap_seconds seconds_per_replicate      projected_B_1000 
                 TRUE                  TRUE                  TRUE 

identical, cell for cell: TRUE

mr_cost_table() reading the stored wall time:
 method mr_seconds bootstrap_B bootstrap_seconds seconds_per_replicate
     FS      0.015           2             3.743                1.8715
 projected_B_1000
             1872

GATE B3: PASS
```

The two coercion warnings above come from my gate script's comparison: it coerced the `method` label to numeric, so that one column's comparison was vacuous there. The printed labels match, and the committed test compares `method` directly and passes.

## B4 (P7): `fs_declaration_table()` (Gate B4: PASS)

A small fit with the declaration field and a three-level `c0` grid, with `s17_column()` transplanted.

```
fit: H-hat {x1 <= -0.7} | family 14 | declaration field kept: TRUE 

fs_declaration_table():
             row   c0 primary alpha n_region T_region T_max kappa_hat
              c0 0.60   FALSE   0.1      134    6.664 6.664 -0.807685
              c0 0.70    TRUE   0.1      134    6.664 6.664 -0.006649
              c0 0.75   FALSE   0.1      134    6.664 6.664  0.360383
 claim_threshold   NA   FALSE   0.1      134    6.664 6.664  2.178527
 implied_level fw_size n_admitted admits_none  source
     -0.580728   0.000          1       FALSE capture
     -0.005305   0.000          1       FALSE capture
      0.281439   0.000          1       FALSE capture
      0.970633   0.272          1       FALSE    <NA>

The article's s17_column() beside the entries read from the table:
 row                                               s17_column        
 Region declared at the default level, N           134               
 T of the region / largest T in the family         6.664 / 6.664     
 Primary protected level c0                        HR 0.70           
 Calibrated cutoff, kappa-hat(c0)                  -0.007            
 Implied level, 2 Phi(kappa-hat) - 1               -0.005            
 Diagnostic at the default screen's cutoff         0.000             
 Calibrated cutoff at the claim threshold          2.179             
 Diagnostic at the claim threshold                 0.272             
 Calibrated rule, every c0 and the claim threshold admits a candidate
 table             
 134               
 6.664 / 6.664     
 HR 0.70           
 -0.007            
 -0.005            
 0.000             
 2.179             
 0.272             
 admits a candidate

every entry identical: TRUE 

GATE B4: PASS
```

On this fixture `kappa_hat` is below 0 at some levels, so `implied_level` is negative. That is the article's formula as given, not an error (OPEN ITEM 3).

## B5 (P8): `mrct_record()`, `fs_in_subgroup()`, `mrct_build_dgm()` (Gate B5: PASS)

```
---- 1. mrct_record(): against the articles' rd()
resolved to: ../../inst/extdata/mrct-consistency 
rep5_trial.csv           1000 rows, 18 columns; identical: TRUE
mr_field_S2_rep5.csv        1 rows, 37 columns; identical: TRUE
predictions_full.csv      305 rows, 13 columns; identical: TRUE
absent directory:   The MRCT records directory (extdata/mrct-consistency) was not found; is forestsearch installed in full? 
absent provenance:  The records directory <dir> has no PROVENANCE.md, so its records cannot be trusted as the shipped ones. 
absent file:        There is no record 'nope.csv' in ../../inst/extdata/mrct-consistency. 

---- 2. fs_in_subgroup(): on a small fit, region as the test set
H-hat: {x1 <= -0.8} | members: 123 of 500 
equal to the raw idiom: TRUE 
equal to the fit's own df.est / df.test classification: TRUE

---- 3. mrct_build_dgm(): against the direct calls
the rep5 build, known to warn: direct 3 warnings; wrapper 3 ; same messages: TRUE 
  Ran out of iterations and did not converge
output equal, attribute aside: TRUE | df_super identical: TRUE 
the S2 build, m = 3, clean: wrapper warnings 0 ; output equal, attribute aside: TRUE | df_super identical: TRUE | model_params identical: TRUE 

GATE B5: PASS
```

## Final stage

**1. Diff confinement (PASS).** The changed existing files are exactly `NAMESPACE`, `NEWS.md` and `_pkgdown.yml`; everything else is new.

```
M	NAMESPACE
M	NEWS.md
A	R/ahr_tables.R
A	R/fs_declaration_table.R
A	R/mr_method_tables.R
A	R/mrct_core_summary.R
A	R/mrct_minor_helpers.R
M	_pkgdown.yml
A	man/ahr_tables.Rd
A	man/fs_declaration_table.Rd
A	man/fs_in_subgroup.Rd
A	man/mr_method_tables.Rd
A	man/mrct_build_dgm.Rd
A	man/mrct_core_summary.Rd
A	man/mrct_record.Rd
A	tests/testthat/test-ahr-tables.R
A	tests/testthat/test-fs-declaration-table.R
A	tests/testthat/test-mr-method-tables.R
A	tests/testthat/test-mrct-core-summary.R
A	tests/testthat/test-mrct-minor-helpers.R
```

```
 NAMESPACE                                  |   9 ++
 NEWS.md                                    |  14 +++
 R/ahr_tables.R                             | 163 ++++++++++++++++++++++++++++
 R/fs_declaration_table.R                   |  91 ++++++++++++++++
 R/mr_method_tables.R                       | 165 +++++++++++++++++++++++++++++
 R/mrct_core_summary.R                      | 140 ++++++++++++++++++++++++
 R/mrct_minor_helpers.R                     | 130 +++++++++++++++++++++++
 _pkgdown.yml                               |   7 ++
 man/ahr_tables.Rd                          |  94 ++++++++++++++++
 man/fs_declaration_table.Rd                |  72 +++++++++++++
 man/fs_in_subgroup.Rd                      |  41 +++++++
 man/mr_method_tables.Rd                    |  84 +++++++++++++++
 man/mrct_build_dgm.Rd                      |  47 ++++++++
 man/mrct_core_summary.Rd                   |  87 +++++++++++++++
 man/mrct_record.Rd                         |  30 ++++++
 tests/testthat/test-ahr-tables.R           |  89 ++++++++++++++++
 tests/testthat/test-fs-declaration-table.R |  97 +++++++++++++++++
 tests/testthat/test-mr-method-tables.R     | 132 +++++++++++++++++++++++
 tests/testthat/test-mrct-core-summary.R    | 141 ++++++++++++++++++++++++
 tests/testthat/test-mrct-minor-helpers.R   | 123 +++++++++++++++++++++
 20 files changed, 1756 insertions(+)
```

**2. Full suite (PASS).** All 573 existing blocks match the reused baseline on expectations, failures, skips and warnings. The baseline's per-block results are kept in the session scratchpad from the P1+P2 run.

```
                        blocks files expectations failed errors skipped_blocks
reused_baseline            573    73         7962      0      0              3
final                      593    78         8126      0      0              3
final_without_new_files    573    73         7962      0      0              3
                        warnings
reused_baseline               38
final                         38
final_without_new_files       38

blocks only in baseline: 0 | only in final, excluding the new files: 0 
blocks compared: 573 | differing: 0 
skip list unchanged: TRUE 

new test files:
  test-mrct-core-summary.R         blocks 4, expectations  32, failed 0, skipped 0, 3.4 s
  test-ahr-tables.R                blocks 3, expectations  46, failed 0, skipped 0, 0.8 s
  test-mr-method-tables.R          blocks 6, expectations  36, failed 0, skipped 0, 11.7 s
  test-fs-declaration-table.R      blocks 3, expectations  23, failed 0, skipped 0, 0.3 s
  test-mrct-minor-helpers.R        blocks 4, expectations  27, failed 0, skipped 0, 1.7 s
```

Skip list, verbatim:

```
1. training_fs_hr for continuous (MD) is identity-scale, not exponentiated ('test-cross-outcome-parity.R:262:3') - Reason: primary fit did not identif
2. fold_numeric_summary row labels are outcome-agnostic ('test-cross-outcome-parity.R:358:3') - Reason: primary fit did not identify; test N/A
3. multisession and sequential tenfold agree on the no-subgroup path ('test-cv-no-subgroup-edges.R:250:5') - Reason: multisession worker could not res
```

**3. `check()` (PASS).**

| | Reused baseline (`84afa853`) | Final (`efc785af`) |
|---|---|---|
| errors | 0 | 0 |
| warnings | 0 | 0 |
| notes | 0 | 0 |
| status | OK | OK |

Verbatim from the examples step on:

```
* checking examples ... [34s/34s] OK
* checking examples with --run-donttest ... [59s/60s] OK
* checking for unstated dependencies in ‘tests’ ... OK
* checking tests ...
  Running ‘testthat.R’ [552s/585s]
 [552s/585s] OK
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes ... OK
* checking re-building of vignette outputs ... [75s/74s] OK
* checking for non-standard things in the check directory ... OK
* checking for detritus in the temp directory ... OK
* DONE
Status: OK
```

**`pkgdown::check_pkgdown()` (PASS):**

```
✔ No problems found.
```

## OPEN ITEMS (non-blocking)

1. **Nine commits, not five.** See Commits. I can squash the four fixes into their items on request; the final tree would be unchanged.
2. **Two problems the single `check()` would have caught, caught earlier.** `with()` in B1 would have raised 8 code-usage NOTEs, and the usual fix is an entry in `R/globals.R`, which this task can't edit. A literal en dash in B1 would have raised a WARNING. Both were found with cheap checks on the new files before the final run, so the one `check()` wasn't spent on them.
3. **`implied_level` can be negative** when `kappa_hat` is below 0, as on the B4 fixture. It is the formula the article uses; the roxygen does not call this out.
4. **Choices where the spec left room.**
   - B1 is one row per quantity (11 rows), so the core table's combined "Sens / Spec" rows become two rows each. Pasting `display` restores them, as the gate does.
   - `Ghat_share` and `trigger` have no SE, because the core table gives none and a new one would be a new statistic.
   - In B2, a single unnamed `x` labels its cell `"cell"`. `analysis_time` and `max_entry` are named arguments on `ahr_design_table()`, beyond the doc's sketch, because a DGM needs them.
   - B3's `multiplier` interval reads `NA` on `wald` or winner-variant fits, following `mr_values()` (D2). The article's `three_columns()` would print those stored limits. Default fits are identical to it, as the gate shows.
   - B3 stops if a fit has no de-biased complement, since the table's rows need one.
   - B4's table also accepts a bare `fs_mr_inference()` result; `source` is a column rather than the article's `stopifnot`, so a `field_matrix` read is reported rather than refused.
5. **Not done, as out of scope:** no article or vignette moved to these functions, no version bump, and the installed library copy was not reinstalled.
