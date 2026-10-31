# Package facts — installed forestsearch 0.4.1 (= `R/` at this branch; `git diff c2967362 HEAD -- R/` empty)

## 1. Where each Table S16 cell lives in a GLM fit (OR and MD paths alike)

With `g <- fit$mr_inference`, `cc <- g$complement`, `fb <- g$field`, `fc <- fb$complement`, `js <- fb$joint_s`
(the paths of the 2026b documents' `extract` chunk, `calls_2026b/*_assembly.R`). Every path below returned the
Table S16 value on M-bin and M-cont (`rungs/M-bin_compare.csv`, `rungs/M-cont_compare.csv`). All are on the natural
scale (OR, or MD in CD4 cells); `to_eff()` maps the log/identity scale back.

| Table S16 cell | element | written at |
|---|---|---|
| Region | `fit$sg.harm` (pasted with `" & "`); membership `fit$df.est$treat.recommend == 0` | — |
| N (events) | `sum(in_H)`; events `sum(y[in_H])` (binary only) | — |
| Unadjusted, 95% CI | `g$naive$est`, `$lower`, `$upper` | `R/fs_mr_inference.R:1066` |
| Corrected (two-term), IJ interval | `g$debiased$est`, `$lower`, `$upper` | `R/fs_mr_inference.R:1069` |
| Field, lower bound | `fb$lower_1s` | `R/fs_mr_inference.R:980` |
| Complement unadjusted / corrected | `cc$naive$*`, `cc$debiased$*` | `R/fs_mr_inference.R:857`, `:860` |
| Field-s, upper bound | `fc$upper_1s_s` (needs `field_scale_complement = "selected"`) | `R/fs_mr_inference.R:1396` |
| Bonferroni pair | `js$bonf_lower_H`, `js$bonf_upper_Hc` | `R/fs_mr_inference.R:1448` |
| Selection bias | `g$selection_bias` (log scale for OR, cells for MD) | `R/fs_mr_inference.R:1089` |
| Re-selection family | `g$n_family` | `R/fs_mr_inference.R:1094` |
| Whole trial | not in the fit: the documents compute it (`reference-effects` chunk; `glm(binomial)` Wald for OR, `lm` for MD) | `calls_2026b/*_data_prep.R` |

## 2. What `print()` and `summary()` show on these fits

Captured on the M-bin and M-cont fits: `logs/print_summary_glm.txt`. Code: `R/forestsearch_methods.R`
(`print.forestsearch` 244, `summary.forestsearch` 367, the shared block `.fs_print_mr_products` 158–197,
fed by `.fs_mr_products` 82–125).

- **Shown:** the field lower bound in Ĥ, the field-s upper bound in Ĥᶜ (tagged `[field-s]`), the Bonferroni pair
  with its γ, the two-sided IJ interval of Ĥ (lower, upper only), p̂(H); `summary()` adds the field, field-s and IJ
  standard errors, the top-3 re-selection mass, and the family size (`n_family`). The bound lines carry the
  measure's name (`OR`, `MD`) from `g$measure`.
- **Not shown:** the corrected point estimate `g$debiased$est`, the unadjusted estimate and interval in either
  region, the complement's corrected estimate and IJ interval, and `g$selection_bias`. An article that reports
  those reads them from `fit$mr_inference`.
- **Finding:** the selected-subgroup header prints the consistency table's `hr` column under the label "HR"
  (`print`, line 277: `"  HR:"`) and "Hazard ratio" (`summary`, line 499). On the GLM paths that column is the
  log odds ratio or the mean difference (the 2026b binary document says so), so the binary fit prints
  `HR: 1.277` (= log 3.584) and the continuous fit `HR: 71.096` (an MD in cells). Not an `R/` edit here; recorded.

## 3. `fs_declaration_calibration()` on the MD path

`c0` is given on the natural scale of `c2` and mapped by `.fs_decl_c0_cmp()` (`R/fs_declaration_calibration.R:124`):
`log()` on a ratio path (refusing `c0 <= 0` there, line 138), **the identity otherwise** (line 145), with `c0 <= c2`
required (line 150). A negative `c0` on the difference scale is therefore accepted as is: the D-cont rung ran the
document's `fs_declaration_calibration(fit, alpha = 0.10, c0 = c(-5, -10, -15))` and its `c0 = NULL` boundary call
without error and reproduced every Table S17 MD cell (`logs/compare_D.txt`). The capture must be on at the fit
(`mr_inference_args$keep_declaration_field = TRUE`, `declaration_c0 = …`); the documents assert `cal$c0$source == "capture"`.

## 4. `forestsearch_bootstrap_dofuture()` and `mr_estimates_table()` on a binary fit

Read from the code only: the bootstrap rungs (B-S3-*) did not run (Part 2 stopped at Gate C).
**v2:** the B-S3 rungs ran (report §11.3): on all three binary fits the call returns `H_estimates` with `H2`, `H2_lower`,
`H2_upper`, `results$tmins_iteration` per replicate, and the `timing` attribute (`avg_seconds_per_boot`). Defaults in force:
`seed = 8316951`, `mr_in_replicates = FALSE`. `mr_estimates_table()` was not called.

- `forestsearch_bootstrap_dofuture(fs.est, nb_boots, seed = 8316951L, details, show_three, parallel_args, digits,
  mr_in_replicates)` — no `est.scale` argument; the replicate re-runs the search from `fs.est$args_call_all`.
  Returns `SG_CIs`, `FSsg_tab`, `H_estimates`, `Hc_estimates` (roxygen `@return`, `R/bootstrap_dofuture_main.R` 89–94); the
  bias-corrected H estimate is `H_estimates$H2` with `H2_lower`/`H2_upper` (line 656) — the fields the 2026a
  payload's `fb_*` columns read (`calls_2026a/table_assembly.R`, `.fb_triple`). Per-replicate time is
  `results$tmins_iteration` (`R/bootstrap_analysis_dofuture.R:780`, documented at 141–142); the
  `timing` attribute carries the total and the mean per replicate (lines 1004–1010).
- `mr_estimates_table(fs, boot_results, est.scale = "hr", digits = 2)` (`R/mr_estimates_table.R:43`):
  `est.scale` is documented as `"hr"` or `"1/hr"` only (line 18) and passed through to `format_bootstrap_table()`
  (line 109); it needs `boot_results$FSsg_tab` and `$SG_CIs` as a template. The bias-corrected column is found by
  pattern `(HR|IRR|OR|RR|RD|IRD|MD)\*` (line 60), so a binary fit's `OR*` column is matched without an OR-specific
  `est.scale`. The 2026a document notes the same for `summarize_bootstrap_results()` (the OR column "is
  auto-detected from FSsg_tab; the est.scale argument only affects the legacy HR-flavored caption helper").
  Whether `"1/hr"` has any meaning on an OR fit was not exercised.
