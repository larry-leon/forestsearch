# REPORT — `mrct_trial_analysis()` and `mrct_four_sets()` on a shared core (2026-11-06)

Task: `dev/tasks/TASK_mrct_trial_analysis_2026-11-06.md` (review proposals P1+P2, Option B).
**Outcome: Gates 0, 1, 2 and 3 pass, with 0 errors, 0 warnings and 0 notes from `check()`, the same as baseline.**
The work is on branch `fs-mrct-trial`, unmerged and unpushed, for your review.

## Branch, HEAD and commits

| | |
|---|---|
| Start | `master` at `73a0dff7`, even with `origin/master` |
| Task commit | `168e6a1e` on `master`, the task doc only |
| Work branch | `fs-mrct-trial`, created from `168e6a1e` |
| End | the top of `fs-mrct-trial`: this report's commit, on top of `84afa853` |

Commits on `fs-mrct-trial` beyond `master`:

- `d9fa0255` refactor: `mrct_case_study()` on a shared analysis core (Stage 1)
- `84afa853` feat: `mrct_trial_analysis()` and `mrct_four_sets()` (Stage 2)
- this report (Stage 4)

## Stage 0: source map and baseline (Gate 0: PASS)

The analysis code separates cleanly from the simulation bookkeeping, so the task's STOP condition did not apply.

### Line-partition map of `R/mrct_case_study.R` at `168e6a1e`

**Simulation-specific: stays in `mrct_case_study()`.**

| Lines | What |
|---|---|
| 247-257 | Argument checks (`dgm`, `seed`, `k`, `region_var` in `dgm$df_super`, `n` default) |
| 259-268 | Saving the caller's generator and stream; restoring them on exit |
| 270-279 | The draw: `sim_call_args`, `L'Ecuyer-CMRG`, `simulate_from_dgm(seed = seed + k)`, generator restored |
| 285-293 | The identifier's defaults, with `confounders.name` read from `dgm$df_super` |
| 294-298 | The MR seed: `mr_inference_args$seed`, else `seed + k` |
| 353-354 | The "no flag_harm column" message |
| 364-367 | `call_args` |
| 386-391 | The print header (trial, seed, H-hat, trigger) |
| 447, 479 | `match.arg(which)`; `analysis_time` read from `call_args` for the plot |

**Analysis-generic: moved to `R/mrct_trial_helpers.R`.**

| Lines | What | Now in |
|---|---|---|
| 281-282 | Split into non-region and region | `.mrct_analyze_trial()` |
| 299-308 | The `forestsearch()` call's arguments | `.mrct_analyze_trial()` |
| 309-319 | The call, with warnings muffled and kept | `.mrct_analyze_trial()` |
| 322-328 | The trial as a data frame; `declared`; `in_Hhat` | `.mrct_analyze_trial()` |
| 329-345 | Set indices, labels, the Cox summaries | `.mrct_four_set_table()` |
| 347-352 | The confusion tables | `.mrct_four_set_table()` |
| 357-358 | The trigger | `.mrct_trigger()` |
| 385 | The number formatter | `.mrct_fm()` |
| 392-398 | Printing the set table | `.mrct_print_sets()` |
| 399-403 | Printing the confusion tables | `.mrct_print_confusion()` |
| 404-423 | Printing the MR block | `.mrct_print_mr()` |
| 448-502 | Forest and Kaplan-Meier plots | `.mrct_plot_sets()` |
| 507-514 | `.mrct_case_set()` | `.mrct_set_rows()` |

The object assembly (360-369) stays in `mrct_case_study()`, reading from the core's result.

### Facts the task asked to record

- **Class vector:** the single string `"mrct_case_study"` (line 369).
- **Seeding on the case-study route:**
  - the draw runs under `L'Ecuyer-CMRG` with seed `seed + k`, and the caller's generator is then restored (lines 276-279);
  - the identifier seeds its own random steps;
  - MR runs afterwards under the caller's generator, with `mr_inference_args$seed` or else `seed + k` (lines 294-298).
- **Row set, order and labels** (lines 332-336): `itt_all` "All patients, ITT"; `itt_region` "Region, ITT"; `region_Gstar` "Region minus H* (G*)", only when the trial has `flag_harm`; `region_Ghat` "Region minus H-hat (G-hat)", the whole region when nothing is declared.
- **Conditional components:**
  - `mr` is `NULL` without MR, and `mr_inference_args` is then `list()`;
  - `confusion` and the G* row need `flag_harm`;
  - the MR print block needs `mr$debiased`;
  - the Kaplan-Meier figure needs `weightedsurv`, and is skipped with a message otherwise.

### Gate data and captures

All present in `inst/extdata/mrct-consistency/`: `rep5_trial.csv`, `mr_field_S2_rep5.csv`, `s2_example_four_sets.csv`, `s2_example_trials.csv`, `s2_example_definitions.md`, `latent_R.csv`, `PROVENANCE.md`. The S2 build ran on the shipped `latent_R` in under a second, with no live latentcor fit, using the `reader-data` and `build-fn` chunks of `vignettes/articles/mrct-consistency-simulations.qmd`.

**Baseline captures,** saved in the session scratchpad at `168e6a1e`, each with its object, print text, messages and both plots as PDF:

- **FX-fast:** a synthetic cell (500 seed patients, bootstrap expansion, no latentcor), one case study with `mr_inference = FALSE`: about 1.3 s.
- **FX-rep5:** the S2 build and the article's `find-case` call for replicate 5 with `mr_inference = TRUE` and seed 613800000 + 5: about 60 s.

**A control I added, because `fs` holds run times.** Each capture was made twice, and the two baseline runs were compared component by component. Only timing fields differ: 2 in FX-fast (`fs$find.grps$time_search`, `fs$minutes_all`) and those plus 6 MR timing fields in FX-rep5. So Gate 1 could be strict: **identical everywhere except those named fields.**

## Stage 1: the refactor (Gate 1: PASS)

These outputs are from the final commit `84afa853`. Every compared leaf is `identical()`, which is stricter than the task's `all.equal(tolerance = 1e-12)`; the listed rows are the timing fields, which differ between any two runs. The plot PDFs match as well, beyond the task's "execute cleanly".

```
================ GATE 1: mrct_case_study() against the baseline capture ================

--- FX-fast: every component, baseline vs refactored
class: mrct_case_study | identical to baseline: TRUE 
fields identical: TRUE | leaves identical: 546 
timing fields excluded (they differ between two baseline runs too): 2 
 path                       kind  detail               
 x$fs$find.grps$time_search float max abs diff 0.000117
 x$fs$minutes_all           float max abs diff 0.000417
differences outside the timing fields: 0 
print text byte-identical: TRUE | messages identical: TRUE 
plot(which = "forest"): runs, PDF content identical (dates removed): TRUE
plot(which = "km"): runs, PDF content identical (dates removed): TRUE
caller's generator and stream restored: TRUE 

--- FX-rep5: every component, baseline vs refactored
class: mrct_case_study | identical to baseline: TRUE 
fields identical: TRUE | leaves identical: 1192 
timing fields excluded (they differ between two baseline runs too): 8 
 path                                              kind  detail             
 x$fs$find.grps$time_search                        float max abs diff 0.0001
 x$fs$minutes_all                                  float max abs diff 0.0018
 x$fs$mr_inference$timing_seconds                  float max abs diff 0.399 
 x$fs$mr_inference$field$timing_seconds            float max abs diff 0.233 
 x$fs$mr_inference$field$complement$timing_seconds float max abs diff 0.379 
 x$mr$timing_seconds                               float max abs diff 0.399 
 x$mr$field$timing_seconds                         float max abs diff 0.233 
 x$mr$field$complement$timing_seconds              float max abs diff 0.379 
differences outside the timing fields: 0 
print text byte-identical: TRUE | messages identical: TRUE 
plot(which = "forest"): runs, PDF content identical (dates removed): TRUE
plot(which = "km"): runs, PDF content identical (dates removed): TRUE

--- FX-rep5: the redrawn trial against rep5_trial.csv: 18 columns, per-column all.equal at 1e-12: TRUE; largest relative difference 0

--- FX-rep5, mrct_case_study() (refactored) against the committed record
          set reproduced committed abs_difference
      itt_all      0.808     0.808        2.2e-16
   itt_region      1.036     1.036        1.6e-15
 region_Gstar      0.736     0.736        2.2e-16
  region_Ghat      0.844     0.844        4.4e-16
H-hat: {z_cd40 <= 236} & !{z_preanti <= 196} | record: {z_cd40 <= 236} & !{z_preanti <= 196} 
2 x 2 non-region: 65 8 65 668 | record: 65 8 65 668 
2 x 2 region:     18 2 17 157 | record: 18 2 17 157 
            quantity      ours    record abs_difference
             naive H 2.6613736 2.6613736        8.9e-16
          debiased H 1.6285027 1.6285027        2.9e-15
          IJ lower H 0.8163578 0.8163578        1.1e-16
          IJ upper H 3.2486012 3.2486012        4.4e-16
       field lower H 0.8498440 0.8498440        4.4e-16
            naive Hc 0.7249723 0.7249723        0.0e+00
         debiased Hc 0.7526280 0.7526280        3.3e-16
      field upper Hc 0.9039452 0.9039452        4.4e-16
  Bonferroni lower H 0.7198488 0.7198488        4.4e-16
 Bonferroni upper Hc 0.9297492 0.9297492        3.3e-16
largest |ours - record|: 2.9e-15 | p_hat: 0.186 record: 0.186 

GATE 1: PASS
```

## Stage 2: the new entry points (Gate 2: PASS)

**Route equivalence, verbatim.** `rep5_trial.csv` as an observed trial, with the application article's identifier arguments and MR seed:

```
================ GATE 2: mrct_trial_analysis() on rep5_trial.csv as an observed trial ================
run: 58.1 s; class mrct_trial_analysis; stage complete; trigger TRUE

MRCT trial analysis: n = 1000, region 'z_region'
Stage: complete
H-hat: {z_cd40 <= 236} & !{z_preanti <= 196}
Trigger (ITT ub < 1.0 and region ITT > 1.0): TRUE

 Set                        n    exp n/events ctl n/events HR    95% CI      
 All patients, ITT          1000 500/187      500/218      0.808 0.665, 0.983
 Region, ITT                 194 104/22       90/18        1.036 0.556, 1.932
 Region minus H* (G*)        159 83/12        76/14        0.736 0.340, 1.592
 Region minus H-hat (G-hat)  174 92/16        82/16        0.844 0.422, 1.688

H-hat against H* (TP, FP, FN, TN):
  nonregion  65, 8, 65, 668
  region     18, 2, 17, 157

Multiplier resampling (field, 2000 draws):
 set                           naive debiased IJ 95% CI    field 1-sided 95%
 H-hat                         2.661 1.629    0.816, 3.249 >= 0.850         
 H-hat complement (non-region) 0.725 0.753    0.482, 1.175 <= 0.904         
 field-s 1-sided 95% Bonferroni Bonferroni, field-s
 -                   >= 0.720   >= 0.720           
 <= 0.908            <= 0.930   <= 0.933           

--- real route against the simulated route, same patients
H-hat identical: TRUE | {z_cd40 <= 236} & !{z_preanti <= 196} 
four sets: rows, labels and counts identical: TRUE | largest |HR, CI or SE difference|: 0 
confusion identical: TRUE 
in_Hhat identical: TRUE 
mr_values() across the two routes: largest |difference|, timing excluded: 0 

--- mrct_trial_analysis() against the committed record
          set reproduced committed abs_difference
      itt_all      0.808     0.808        2.2e-16
   itt_region      1.036     1.036        1.6e-15
 region_Gstar      0.736     0.736        2.2e-16
  region_Ghat      0.844     0.844        4.4e-16
H-hat: {z_cd40 <= 236} & !{z_preanti <= 196} | record: {z_cd40 <= 236} & !{z_preanti <= 196} 
2 x 2 non-region: 65 8 65 668 | record: 65 8 65 668 
2 x 2 region:     18 2 17 157 | record: 18 2 17 157 
            quantity      ours    record abs_difference
             naive H 2.6613736 2.6613736        8.9e-16
          debiased H 1.6285027 1.6285027        2.9e-15
          IJ lower H 0.8163578 0.8163578        1.1e-16
          IJ upper H 3.2486012 3.2486012        4.4e-16
       field lower H 0.8498440 0.8498440        4.4e-16
            naive Hc 0.7249723 0.7249723        0.0e+00
         debiased Hc 0.7526280 0.7526280        3.3e-16
      field upper Hc 0.9039452 0.9039452        4.4e-16
  Bonferroni lower H 0.7198488 0.7198488        4.4e-16
 Bonferroni upper Hc 0.9297492 0.9297492        3.3e-16
largest |ours - record|: 2.9e-15 | p_hat: 0.186 record: 0.186 

GATE 2 (route equivalence): PASS
```

**D9 tests, verbatim.** Eight blocks: the FX-fast shape and route agreement, one engineered trial per stage, `proceed_without_trigger`, truth absent and present, the MR print, input checks, and `mrct_four_sets()` against direct `coxph()` at 1e-12. Only the 9-second MR print block carries `skip_on_cran()`.

```
ℹ Testing forestsearch
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 131 ]
 test                                                                            
 FX-fast: mrct_case_study() keeps its fields, class, sets and print layout       
 FX-fast: the three entries agree on the same patients                           
 each stage on an engineered trial, and never an error                           
 proceed_without_trigger runs the identifier and still reports the trigger       
 truth absent: no G* row, no confusion, no error, no message; truth present: both
 with MR, print reports the evaluated construction set, by location              
 the input is checked                                                            
 mrct_four_sets() on a hand-built data frame equals direct coxph() summaries     
 nb failed error skipped warning
 17 0      FALSE FALSE   0      
  8 0      FALSE FALSE   0      
 32 0      FALSE FALSE   0      
  7 0      FALSE FALSE   0      
 10 0      FALSE FALSE   0      
 10 0      FALSE FALSE   0      
  8 0      FALSE FALSE   0      
 39 0      FALSE FALSE   0      
TOTAL expectations 131 | failed 0 | errors 0 | skipped 0 | warnings 0 
```

## Stage 3: preservation (Gate 3: PASS)

**1. Diff confinement.** The only existing `R/` file changed is `R/mrct_case_study.R`; `man/mrct_case_study.Rd` is unchanged.

```
M	NAMESPACE
M	NEWS.md
M	R/mrct_case_study.R
A	R/mrct_four_sets.R
A	R/mrct_trial_analysis.R
A	R/mrct_trial_helpers.R
M	_pkgdown.yml
A	man/mrct_four_sets.Rd
A	man/mrct_trial_analysis.Rd
A	tests/testthat/test-mrct-trial-analysis.R
```

```
 NAMESPACE                                 |   7 +
 NEWS.md                                   |  15 ++
 R/mrct_case_study.R                       | 187 +++-----------------
 R/mrct_four_sets.R                        | 133 ++++++++++++++
 R/mrct_trial_analysis.R                   | 171 ++++++++++++++++++
 R/mrct_trial_helpers.R                    | 281 ++++++++++++++++++++++++++++++
 _pkgdown.yml                              |   2 +
 man/mrct_four_sets.Rd                     | 105 +++++++++++
 man/mrct_trial_analysis.Rd                | 151 ++++++++++++++++
 tests/testthat/test-mrct-trial-analysis.R | 271 ++++++++++++++++++++++++++++
 10 files changed, 1161 insertions(+), 162 deletions(-)
```

**2. Test suite.** All 565 existing blocks match the baseline on expectations, failures, skips and warnings, and the skip list is unchanged.

```
                       blocks files expectations failed errors skipped_blocks
baseline                  565    72         7831      0      0              3
final                     573    73         7962      0      0              3
final_without_new_file    565    72         7831      0      0              3
                       warnings
baseline                     38
final                        38
final_without_new_file       38

blocks only in baseline: 0 | only in final, excluding the new file: 0 
blocks compared: 565 | differing: 0 

skipped, baseline:
                        file
 test-cross-outcome-parity.R
 test-cross-outcome-parity.R
 test-cv-no-subgroup-edges.R
                                                                    test
 training_fs_hr for continuous (MD) is identity-scale, not exponentiated
                    fold_numeric_summary row labels are outcome-agnostic
       multisession and sequential tenfold agree on the no-subgroup path
skipped, final:
                        file
 test-cross-outcome-parity.R
 test-cross-outcome-parity.R
 test-cv-no-subgroup-edges.R
                                                                    test
 training_fs_hr for continuous (MD) is identity-scale, not exponentiated
                    fold_numeric_summary row labels are outcome-agnostic
       multisession and sequential tenfold agree on the no-subgroup path
skip list unchanged: TRUE 

new file: blocks 8 expectations 131 failed 0 seconds 10.4 
existing test-mrct-case-study.R: blocks 7 expectations 64 failed 0 skipped 0 
```

**3. `check()`, side by side.**

| | Baseline (`168e6a1e`) | Final (`84afa853`) |
|---|---|---|
| errors | 0 | 0 |
| warnings | 0 | 0 |
| notes | 0 | 0 |
| status | OK | OK |

Baseline, verbatim from the examples step on:

```
* checking examples ... [31s/31s] OK
* checking examples with --run-donttest ... [44s/45s] OK
* checking for unstated dependencies in ‘tests’ ... OK
* checking tests ...
  Running ‘testthat.R’ [531s/564s]
 [531s/564s] OK
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes ... OK
* checking re-building of vignette outputs ... [74s/73s] OK
* checking for non-standard things in the check directory ... OK
* checking for detritus in the temp directory ... OK
* DONE

Status: OK
```

Final, verbatim from the examples step on:

```
* checking examples ... [32s/32s] OK
* checking examples with --run-donttest ... [46s/46s] OK
* checking for unstated dependencies in ‘tests’ ... OK
* checking tests ...
  Running ‘testthat.R’ [548s/583s]
 [548s/583s] OK
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes ... OK
* checking re-building of vignette outputs ... [77s/76s] OK
* checking for non-standard things in the check directory ... OK
* checking for detritus in the temp directory ... OK
* DONE
Status: OK
```

**4. `pkgdown::check_pkgdown()`:**

```
✔ No problems found.
```

## Before and after: `mrct_case_study()`'s structure

| | Before | After |
|---|---|---|
| `R/mrct_case_study.R` | 514 lines | 377 lines |
| `mrct_case_study()` | 133 lines | 83 lines |
| `print.mrct_case_study()` | 42 lines | 13 lines |
| `plot.mrct_case_study()` | 59 lines | 9 lines |

```
mrct_case_study()                       mrct_trial_analysis()
  checks; generator saved                 checks (columns, 0/1, no NA)
  the draw (seed + k)                     ITT and region fits; stage
  defaults; MR seed + k                   stops here unless triggered or told to proceed
        \                                 /
         .mrct_analyze_trial()  (R/mrct_trial_helpers.R)
           split -> forestsearch() -> H-hat -> sets, confusion, trigger

print and plot methods of mrct_case_study, mrct_trial_analysis, mrct_four_sets
  -> .mrct_print_sets(), .mrct_print_confusion(), .mrct_print_mr(), .mrct_plot_sets()
```

`R/mrct_trial_helpers.R` is 281 lines, `R/mrct_trial_analysis.R` 171 and `R/mrct_four_sets.R` 133.

## OPEN ITEMS (non-blocking)

1. **`id.name` is required in `mrct_trial_analysis()`.** D3 lists the variable-name arguments without saying this. `forestsearch()`'s `id.name` defaults to `"id"`, and without such a column it fails with a cryptic "undefined columns selected". My first draft had `id.name = NULL` and the roxygen example caught it. Now a missing or wrong column gets a clear error.
2. **Defaults that don't mirror the case study.** `confounders.name` is required, because the case study's default reads `dgm$df_super`, which doesn't exist here. `conf_force` defaults to `character(0)`, because the case study's default names simulation columns (`z_age`, `z_bm`). `sg_focus`, `maxk` and the thresholds use the case study's defaults.
3. **Reading of D4 with `proceed_without_trigger = TRUE`.** A declared subgroup gives stage `"complete"` even though the trial has no trigger; `$trigger` still says `FALSE`. "Never `stop()`" is applied to the stages: invalid input still errors, and an error inside `forestsearch()` is not caught.
4. **One addition beyond the spec.** The new plot methods take `analysis_time`, defaulting to the longest follow-up, because an observed trial has no recorded analysis time for the Kaplan-Meier axis.
5. **Two MR print blocks, on purpose.** The case study's block is frozen (D2) and shows the plain field bound and the `joint` pair. It still reads `debiased$lower` and `$upper` as stored, including on `wald` or winner-variant fits, as it always has. The new class's block follows D6 through `mr_values()`, adding the field-s bound and the `joint_s` pair. With 8 columns it wraps at the default console width.
6. **MR seed.** Without `mr_inference_args$seed`, `mrct_trial_analysis()` takes `forestsearch()`'s default seed, since there is no `seed + k` for an observed trial. The roxygen says to give one for a reproducible result.
7. **D9's FX-fast test checks shape and agreement, not absolute numbers.** The baseline capture is not in the tree, so the committed test asserts the frozen fields, class, sets and print layout, and the agreement of the three entries on the same patients. Absolute identity against baseline is Gate 1, run once.
8. **Not done, as out of scope:** no article or vignette edit (the application article still draws its plots by hand), no version bump, and the installed library copy was not reinstalled.
