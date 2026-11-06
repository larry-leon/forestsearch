# REPORT — `mr_values()`, a tidy accessor for `fit$mr_inference` (2026-11-06)

Task: `dev/tasks/TASK_mr_values_2026-11-06.md` (review proposal P5).
**Outcome: Gates 0, 1 and 2 pass, with 0 errors, 0 warnings and 0 notes from `check()`, the same as baseline.**
The work is on branch `fs-mr-values`, unmerged and unpushed, for your review. One thing needs your decision before a pkgdown build: OPEN ITEM 1.

## Branch, HEAD and commits

| | |
|---|---|
| Start | `master` at `e1ed6f1f` (even with `origin/master`) |
| Task commit | `92792bf0`, the task doc only, on `master` |
| Work branch | `fs-mr-values`, created from `92792bf0` before any `R/` change |
| End | the top of `fs-mr-values`: the commit that adds this report, on top of `90eed9b9` |

Commits on `fs-mr-values` beyond `master`:

- `90eed9b9` feat: `mr_values()` (Stage 1)
- the report commit (Stage 3)

The task doc says to work on the current checkout's branch, which was `master`. You asked during the run for R changes to be on a branch and for the task-doc commit to stay on `master`, so `master` is one commit ahead of origin with the task doc only.

## Stage 0: the field map (Gate 0: PASS)

Every D4 quantity has a verified source path, so nothing had to be declared absent under all configurations. Lines are in `R/fs_mr_inference.R` unless stated; paths are inside `mr_inference`.

| Column | Row | Source path | Exists when | Stored scale | Lines |
|---|---|---|---|---|---|
| `n` | H | `n_selected` | full return | count | 1094 |
| `n` | Hc | `complement$n` | complement fit | count | 883 |
| `naive_est`, `naive_lower`, `naive_upper` | H | `naive$est`, `$lower`, `$upper` | full return | effect | 1066-1068 |
| same | Hc | `complement$naive$...` | complement fit | effect | 857-859 |
| `debiased_est` | H | `debiased$est` | full return | effect | 1069 |
| same | Hc | `complement$debiased$est` | complement fit | effect | 860 |
| `ij_lower`, `ij_upper` | H | `debiased$lower`, `$upper` | `ci_method` is `field` or `ij`, and `ij_residual = "two_term"` (SE chosen at 766-770) | effect | 1070-1071 |
| same | Hc | `complement$debiased$lower`, `$upper` | same, and complement fit (SE chosen at 852-854) | effect | 861-862 |
| `ij_se` | H | `debiased$se_ij_two_term` | full return, every `ci_method` and `ij_residual` | **working** | 1078 (755) |
| same | Hc | `complement$debiased$se_ij_two_term` | complement fit | **working** | 869 (846) |
| `field_est2` | H | `field$est2` | `ci_method = "field"` with at least 2 usable draws | effect | 979 |
| same | Hc | `field$complement$est2` | complement field ran (`field_complement` and `include_complement`, 938) | effect | 1382 |
| `field_est2_s` | Hc | `field$complement$est2_s` | same, and `field_scale_complement = "selected"` (1258, 1393) | effect | 1395 |
| `field_lower_1s` | H | `field$lower_1s` | as `field_est2` on H | effect | 980 |
| `field_upper_1s` | Hc | `field$complement$upper_1s` | as `field_est2` on Hc | effect | 1384 |
| `field_upper_1s_s` | Hc | `field$complement$upper_1s_s` | as `field_est2_s` | effect | 1396 |
| `bonf_lower_H`, `bonf_upper_Hc` | H, Hc | `field$joint$bonf_lower_H`, `$bonf_upper_Hc` | complement field ran (attached at 1008, built at 1373) | effect | 1448 |
| `bonf_lower_H_s`, `bonf_upper_Hc_s` | H, Hc | `field$joint_s$...` | same, and `field_scale_complement = "selected"` (attached at 1009, built at 1375) | effect | 1448 |
| `selection_bias` | H | `selection_bias` | full return | **working** | 738, 1089 |
| `p_hat` | H | `reselection$p_hat[[selected_index]]` | `return_reselection = TRUE` | probability | 1096-1099, 1063 |
| `n_family` | H | `n_family` | always, short return included | count | 1094, 641 |
| `n_selected` | H | `n_selected` | full return | count | 1094 |
| `timing_seconds` | H | `timing_seconds` | full return; the main multiplier stage only | seconds | 715, 1095 |
| `field_timing_seconds` | H | `field$timing_seconds` | as `field_est2` on H | seconds | 987 |
| `field_complement_timing_seconds` | H | `field$complement$timing_seconds` | as `field_est2` on Hc | seconds | 1391 |

**Points the task asked to verify:**

- **Scales.** Point estimates, limits and bounds are on the effect scale; `selection_bias` and the SEs are on the working scale (log for ratio measures). This matches the function's own note at lines 500-503 and the code (`to_eff()` is applied at 1066-1072; `selection_bias` at 738 is not transformed).
- **`joint` and `joint_s`, `upper_1s` and `upper_1s_s`.** The current source produces both of each under the recommended defaults. The `_s` members need `field_scale_complement = "selected"`. They are separate columns, as D4 specifies.
- **Which measures are ratio measures.** `log_scale` (1064) comes from `.fs_mr_assemble()` (line 104), set by the pieces: `TRUE` for Cox (`R/consistency_resample.R:90`) and for `OR`, `RR`, `IRR` (`R/consistency_resample.R:323`), `FALSE` for `MD` and `RD`. `scale = "log"` reads this flag.
- **Absent shapes that must not produce values:**
  - the short return (637-642);
  - a note-only complement (885);
  - a note-only field (1047);
  - a note-only complement field (1261, 1334);
  - a note-only joint (1431).

  `list(note = ...)$n` partially matches `note`, so the accessor reads with exact matching (`[[`) throughout; a test covers this on a real note-only fit.

**The `forestsearch()` pass-through.** `R/forestsearch_main.R:3844` calls `fs_mr_inference()` on the consistency path, inside the `tryCatch` that starts at 3801. `:3955` stores the result unchanged as `mr_inference = mr_out`, which is `NULL` when MR is off or failed. The DINA and GRF paths go through `.fs_apply_mr()` (`:2860`, `:3045`; `R/fs_mr_inference_methods.R:141`), which defaults to `ci_method = "ij"`.

**`mr_estimates_table()`.** Its `"mr_values"` attribute is built at `R/mr_estimates_table.R:247-262` and attached at `:336`. Column mapping: `naive_*`, `debiased_est`, `ij_lower` and `ij_upper` have the same names. `field_bound` is `field_lower_1s` on H and `field_upper_1s` on Hc. `bonferroni_bound` is `bonf_lower_H` on H and `bonf_upper_Hc` on Hc, read from `joint`, not `joint_s`.

## Stage 1: implementation (Gate 1: PASS)

New: `R/mr_values.R`, `tests/testthat/test-mr_values.R`, `man/mr_values.Rd`. Changed: `NAMESPACE` (one `export()` line), `NEWS.md` (one entry under 0.4.1.9002). `DESCRIPTION` and every other `man/` file are unchanged; `document()` did not touch them.

**Fixtures**, all standalone `fs_mr_inference()` calls on 300 synthetic survival patients with a three-candidate hand-built family, `reselection = "maxeff"`, 400 draws, `field_R_out = 300`, `field_R_in = 100`, `seed = 4242`:

- **F1 (full coverage):** `ci_method = "field"`, `include_complement = TRUE`, `field_complement = TRUE`, `field_scale_complement = "selected"`, `return_reselection = TRUE`.
- **F2 (NA-fill):** `ci_method = "ij"`, `include_complement = FALSE`.
- **F3 was not needed:** the standalone call produces every mapped component.
- **Variants, one per remaining absence condition:**
  - `field_complement = FALSE`;
  - `field_scale_complement = "none"`;
  - `return_reselection = FALSE`;
  - `ci_method = "wald"`;
  - `ij_residual = "winner"`;
  - a real note-only complement;
  - a real short return;
  - an `OR` fit;
  - an `MD` fit.

The raw-path reads in the test file are written out separately from the implementation. Each row assertion checks every column: a listed column equals its raw value at tolerance 1e-12, and every other column is `NA`, with the absence of the component asserted directly. The whole file takes about 2 seconds, so no `skip_on_cran()` was needed.

**Verbatim identity-test output:**

```
ℹ Testing forestsearch
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 633 ]
 test                                                                       
 F1 (full coverage): every cell equals its raw mr_inference path            
 F2 (NA-fill): field, joint and complement are NA because they are absent   
 each remaining absence condition gives NA, the rest stays identical        
 only the two-term IJ construction is exposed                               
 a complement that could not be fit: an all-NA Hc row and no spurious values
 the short return gives an H row that is NA apart from n_family             
 scale = 'log' transforms the effect columns of ratio measures only         
 MD: effect returns as stored, log is refused                               
 input detection: a fit, a bare result, and informative errors otherwise    
 the overlapping quantities equal mr_estimates_table()'s attribute          
 nb  failed error skipped warning
 117 0      FALSE FALSE   0      
  55 0      FALSE FALSE   0      
 210 0      FALSE FALSE   0      
  11 0      FALSE FALSE   0      
  61 0      FALSE FALSE   0      
  31 0      FALSE FALSE   0      
   9 0      FALSE FALSE   0      
 116 0      FALSE FALSE   0      
   6 0      FALSE FALSE   0      
  17 0      FALSE FALSE   0      

TOTAL: expectations 633 | failed 0 | errors 0 | skipped 0 | warnings 0 
```

**Verbatim `git diff --stat 92792bf0 90eed9b9`:**

```
 NAMESPACE                       |   1 +
 NEWS.md                         |  10 ++
 R/mr_values.R                   | 257 ++++++++++++++++++++++++++++++++++
 man/mr_values.Rd                | 133 ++++++++++++++++++
 tests/testthat/test-mr_values.R | 300 ++++++++++++++++++++++++++++++++++++++++
 5 files changed, 701 insertions(+)
```

```
M	NAMESPACE
M	NEWS.md
A	R/mr_values.R
A	man/mr_values.Rd
A	tests/testthat/test-mr_values.R
```

No line of an existing `R/` file changed.

## Stage 2: preservation (Gate 2: PASS)

| | Baseline (`92792bf0`) | Final (`90eed9b9`) |
|---|---|---|
| `check()` errors | 0 | 0 |
| `check()` warnings | 0 | 0 |
| `check()` notes | 0 | 0 |
| `check()` status | OK | OK |
| Test failures | 0 | 0 |
| Skipped tests | 3 | 3 (the same three) |
| Test warnings | 38 | 38 |
| Expectations | see below | 7,831 (7,198 existing + 633 new) |

**Verbatim `check()`, baseline, from the examples step on:**

```
* checking examples ... [31s/31s] OK
* checking examples with --run-donttest ... [43s/44s] OK
* checking for unstated dependencies in ‘tests’ ... OK
* checking tests ...
  Running ‘testthat.R’ [517s/549s]
 [517s/549s] OK
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes ... OK
* checking re-building of vignette outputs ... [74s/72s] OK
* checking for non-standard things in the check directory ... OK
* checking for detritus in the temp directory ... OK
* DONE
Status: OK
```

**Verbatim `check()`, final, from the examples step on:**

```
* checking examples ... [31s/31s] OK
* checking examples with --run-donttest ... [44s/45s] OK
* checking for unstated dependencies in ‘tests’ ... OK
* checking tests ...
  Running ‘testthat.R’ [526s/560s]
 [526s/560s] OK
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes ... OK
* checking re-building of vignette outputs ... [75s/74s] OK
* checking for non-standard things in the check directory ... OK
* checking for detritus in the temp directory ... OK
* DONE

Status: OK
```

**Verbatim skip list, baseline:**

```
1. training_fs_hr for continuous (MD) is identity-scale, not exponentiated ('test-cross-outcome-parity.R:262:3') - Reason: primary fit did not identif
2. fold_numeric_summary row labels are outcome-agnostic ('test-cross-outcome-parity.R:358:3') - Reason: primary fit did not identify; test N/A
3. multisession and sequential tenfold agree on the no-subgroup path ('test-cv-no-subgroup-edges.R:250:5') - Reason: multisession worker could not res
```

**Verbatim skip list, final:**

```
1. training_fs_hr for continuous (MD) is identity-scale, not exponentiated ('test-cross-outcome-parity.R:262:3') - Reason: primary fit did not identif
2. fold_numeric_summary row labels are outcome-agnostic ('test-cross-outcome-parity.R:358:3') - Reason: primary fit did not identify; test N/A
3. multisession and sequential tenfold agree on the no-subgroup path ('test-cv-no-subgroup-edges.R:250:5') - Reason: multisession worker could not res
```

**A block-by-block comparison, beyond what the gate asks.** The in-repo baseline run recorded failures, skips and warnings but not its total expectation count, because my wrapper script failed to save its results (see OPEN ITEM 4). To compare block by block, I reran the suite on a clean export of `92792bf0` containing `R`, `tests`, `inst`, `man`, `DESCRIPTION`, `NAMESPACE` and `.Rbuildignore`.

Of the 555 existing test blocks, 547 match exactly on expectations, failures, skips and warnings. The other 8 were skipped in the export because their reference files or modules are not in it ("dev/ not shipped", "simulation modules not present", and similar). They were not skipped in the in-repo baseline or in the final run, so the export can't compare those 8 by count; their evidence is the in-repo baseline's 0 failures and identical skip list. No block appears or disappears.

```
                  blocks files expectations failed errors skipped_blocks
baseline_rerun       555    71         7154      0      0             11
final                565    72         7831      0      0              3
final_without_new    555    71         7198      0      0              3
                  warnings
baseline_rerun          38
final                   38
final_without_new       38

block keys duplicated: baseline 0 final 0 
blocks only in baseline:
character(0)
blocks only in final (excluding the new file):
character(0)

blocks compared: 555 | differing: 8 
                                                                                                   k
             test-betaHhat-contract.R :: T11: frames are bitwise identical to the shims they replace
                             test-betaHhat-contract.R :: T13: theta-dagger matches the shims exactly
             test-betaHhat-contract.R :: T14 negative control: the wrong outcome_type must NOT agree
    test-binary-default-or.R :: the campaign template's call re-resolves to OR on its own thresholds
                                  test-directive-c.R :: Part 5 changes no resolved value on any path
 test-fs-oc-grid.R :: fs_oc_predict() is identical to the 0.2.4 pre-refactor reference on both gates
                           test-fs-oc-null.R :: the alternative path is untouched by the null branch
     test-threshold-pair-directive-a.R :: the probe's validate = FALSE run is the committed baseline
 nb.b failed.b error.b skipped.b warning.b nb.f failed.f error.f skipped.f
    1        0   FALSE      TRUE         0    4        0   FALSE     FALSE
    1        0   FALSE      TRUE         0    3        0   FALSE     FALSE
    3        0   FALSE      TRUE         0    3        0   FALSE     FALSE
    1        0   FALSE      TRUE         0    7        0   FALSE     FALSE
    1        0   FALSE      TRUE         0   28        0   FALSE     FALSE
    1        0   FALSE      TRUE         0    5        0   FALSE     FALSE
    1        0   FALSE      TRUE         0    3        0   FALSE     FALSE
    1        0   FALSE      TRUE         0    1        0   FALSE     FALSE
 warning.f
         0
         0
         0
         0
         0
         0
         0
         0

failed or errored, final:
[1] file   test   failed error 
<0 rows> (or 0-length row.names)
failed or errored, baseline rerun:
[1] file   test   failed error 
<0 rows> (or 0-length row.names)

skipped, final:
                        file
 test-cross-outcome-parity.R
 test-cross-outcome-parity.R
 test-cv-no-subgroup-edges.R
                                                                    test
 training_fs_hr for continuous (MD) is identity-scale, not exponentiated
                    fold_numeric_summary row labels are outcome-agnostic
       multisession and sequential tenfold agree on the no-subgroup path
skipped, baseline rerun:
                              file
          test-betaHhat-contract.R
          test-betaHhat-contract.R
          test-betaHhat-contract.R
          test-binary-default-or.R
       test-cross-outcome-parity.R
       test-cross-outcome-parity.R
       test-cv-no-subgroup-edges.R
                test-directive-c.R
                 test-fs-oc-grid.R
                 test-fs-oc-null.R
 test-threshold-pair-directive-a.R
                                                                           test
                    T11: frames are bitwise identical to the shims they replace
                                    T13: theta-dagger matches the shims exactly
                    T14 negative control: the wrong outcome_type must NOT agree
           the campaign template's call re-resolves to OR on its own thresholds
        training_fs_hr for continuous (MD) is identity-scale, not exponentiated
                           fold_numeric_summary row labels are outcome-agnostic
              multisession and sequential tenfold agree on the no-subgroup path
                                   Part 5 changes no resolved value on any path
 fs_oc_predict() is identical to the 0.2.4 pre-refactor reference on both gates
                           the alternative path is untouched by the null branch
                     the probe's validate = FALSE run is the committed baseline

new file: blocks 10 expectations 633 failed 0 seconds 1.8 
```

## Improvement demonstration (non-gating; no article edited)

The raw-path version is the MR part of `s16_column()` in `vignettes/articles/actg175-glm-outcomes.qmd`. Its `in_region()` and event counts are removed, because they read the data rather than `mr_inference`, and `mr_values()` does not replace them. On F1 the two give identical output.

The size saving is modest: 19 lines against 7 as written, and 704 against 600 characters. The real gains are that a reader no longer needs to know the paths or which of `joint` and `joint_s` is which, the columns are fixed, and absent components are `NA` rather than `NULL` (a `NULL` would make the `sprintf()` calls return an empty result).

The same output also includes an end-to-end check on a real `forestsearch()` fit: replicate 5 of S2, with the application article's `case-fs` call. `mr_values()` agrees with `mr_estimates_table()`'s attribute at a difference of exactly 0, and with the bare result by `identical()`.

```
=== F1: selected g1 
                                [,1]       [,2]     
row                             H          Hc       
n                               151        149      
naive_est                       0.9279521  0.7760111
naive_lower                     0.6438566  0.5378543
naive_upper                     1.337402   1.119621 
debiased_est                    0.8592249  0.8688268
ij_lower                        0.4085225  0.4625036
ij_upper                        1.807164   1.632117 
ij_se                           0.3793353  0.3216842
field_est2                      0.8523576  0.8855541
field_est2_s                    <NA>       0.8859312
field_lower_1s                  0.6274062  <NA>     
field_upper_1s                  <NA>       1.196823 
field_upper_1s_s                <NA>       1.194639 
bonf_lower_H                    0.5985389  <NA>     
bonf_upper_Hc                   <NA>       1.289082 
bonf_lower_H_s                  0.5985389  <NA>     
bonf_upper_Hc_s                 <NA>       1.325338 
selection_bias                  0.08277333 <NA>     
p_hat                           0.4575     <NA>     
n_family                         3         <NA>     
n_selected                      151        <NA>     
timing_seconds                  0.018      <NA>     
field_timing_seconds            0.022      <NA>     
field_complement_timing_seconds 0.051      <NA>     

=== demo
 row                                   raw_paths           
 H: N                                  151                 
 H: unadjusted, 95% CI                 0.928 (0.644, 1.337)
 H: corrected (two-term), IJ interval  0.859 (0.409, 1.807)
 H: field lower bound                  0.627               
 Hc: N                                 149                 
 Hc: unadjusted, 95% CI                0.776 (0.538, 1.120)
 Hc: corrected (two-term), IJ interval 0.869 (0.463, 1.632)
 Hc: field-s upper bound               1.195               
 Bonferroni pair (joint_s)             0.599 / 1.325       
 Selection bias                        0.083               
 Re-selection family                   3                   
 mr_values           
 151                 
 0.928 (0.644, 1.337)
 0.859 (0.409, 1.807)
 0.627               
 149                 
 0.776 (0.538, 1.120)
 0.869 (0.463, 1.632)
 1.195               
 0.599 / 1.325       
 0.083               
 3                   
identical(): TRUE 
body lines: raw_path 14 | via_mr_values 5 (as written above)

=== replicate 5 of S2: real forestsearch() fit; forestsearch ; H-hat {z_cd40 <= 236} & !{z_preanti <= 196} 
                                [,1]      [,2]     
row                             H         Hc       
n                                73       733      
naive_est                       2.6613736 0.7249723
naive_lower                     1.4818101 0.5772907
naive_upper                     4.7799037 0.9104335
debiased_est                    1.628503  0.752628 
ij_lower                        0.8163578 0.4820973
ij_upper                        3.248601  1.174968 
ij_se                           0.3523348 0.2272618
field_est2                      1.5592652 0.7533035
field_est2_s                    <NA>      0.7529483
field_lower_1s                  0.849844  <NA>     
field_upper_1s                  <NA>      0.9039452
field_upper_1s_s                <NA>      0.907819 
bonf_lower_H                    0.7198488 <NA>     
bonf_upper_Hc                   <NA>      0.9297492
bonf_lower_H_s                  0.7198488 <NA>     
bonf_upper_Hc_s                 <NA>      0.9329443
selection_bias                  0.4877466 <NA>     
p_hat                           0.186     <NA>     
n_family                        3003      <NA>     
n_selected                      73        <NA>     
timing_seconds                  4.21      <NA>     
field_timing_seconds            27.455    <NA>     
field_complement_timing_seconds 5.592     <NA>     
attribute overlap: shared columns max |diff| 0 | field 0 | Bonferroni 0 
identical to the bare result: TRUE 
```

## OPEN ITEMS

1. **`pkgdown::check_pkgdown()` now fails, and this needs your decision.** It reports `In _pkgdown.yml, 1 topic missing from index: "mr_values".` It passed before this change. `_pkgdown.yml` is outside the task's permitted change set (Gate 1), so I did not edit it. The fix is one line, `- mr_values`, after `- mr_estimates_table` at `_pkgdown.yml:504`. `R CMD check` is unaffected.
2. **Behaviour you may want to review: `ij_lower` and `ij_upper` are `NA` under `ci_method = "wald"` or `ij_residual != "two_term"`.** On those fits `debiased$lower` and `$upper` hold a robust-SE or winner-variant interval, which the task excludes, and the two-term interval isn't stored; computing it would break "nothing new is computed". `ij_se` still reports the two-term SE. `mr_estimates_table()`'s attribute reads those limits regardless of these settings, so the two differ on such fits. Default fits are unaffected.
3. **Choices where the spec left room.**
   - D2 says the `"Hc"` row appears "when the complement block exists". I took that literally, so a note-only complement gives an all-`NA` `"Hc"` row.
   - D4 says "timing". The source has three separate timings, so there are three columns, and `timing_seconds` covers the main multiplier stage only, not the whole MR step.
   - D2 places `selection_bias` on the H row. The complement's own `complement$selection_bias` (line 882) is therefore not exposed.
   - `n_selected` repeats `n` on the H row, because D4 lists both.
4. **My error in the baseline run.** I passed the output variable to the wrapper wrongly, so saving the results failed after each run had finished. The test and check logs were complete, which is where the baseline figures here come from, but the suite's total expectation count wasn't captured. The clean-export rerun above fills that gap. The final run used a corrected wrapper.
5. **Not done, as out of scope:**
   - no version bump (NEWS is under 0.4.1.9002);
   - the installed library copy of forestsearch was not reinstalled, so `library(forestsearch)` doesn't yet have `mr_values()`;
   - no article uses it yet;
   - `mr_estimates_table()` does not delegate to it.
