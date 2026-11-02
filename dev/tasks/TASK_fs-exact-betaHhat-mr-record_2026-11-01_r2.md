# TASK: forestsearch — the exact survival β(Ĥ) option, and MR bounds in the per-replicate record (r2)

Date: 2026-11-01 (r2 supersedes the unlaunched r1: adds §2b)
Repository: **forestsearch**, branch `mrct-analytic`, from `43441a42` (the end of `TASK_fs-mrct-article_2026-11-01`). mrct-consistency read only. CC does not push.

**`R/` category: changes behaviour, additive, defaults unchanged.** Two new options; no existing call returns anything different unless it asks. Nothing changes the method: the exact route is an option beside the Monte Carlo evaluation, not a replacement.

## 0. First action

Copy this file into `dev/tasks/` and commit it on its own, on `mrct-analytic`.

## 1. The exact survival β(Ĥ) (`R/betaHhat_truth.R`)

Today the survival target is evaluated by `fs_build_eval_frame()`: one `simulate_from_dgm()` draw of the whole super-population (`replace = FALSE`, `analysis_time`, `cens_adjust`, `eval_seed`), then `fs_betaHhat_one()` fits `coxph()` on the realised Ĥ and its complement. The analytic layer (`ahr_limit()`) gives the same population quantity exactly when the censoring is administrative only.

Add `method = c("frame", "exact")` to `fs_attach_betaHhat()` (and to whatever entry points share its signature; read the file), default `"frame"`:

- `"frame"`: today's path, untouched.
- `"exact"`: survival only; `ahr_population(dgm, analysis_time, max_entry = <the frame's entry maximum>)` once, then for each realised rule `ahr_limit(pop, set = <membership of Ĥ in df_super>, psi = pop$psi)` and the same on the complement; the returned columns are those the frame route returns (β(Ĥ), β(Ĥᶜ) and whatever accompanies them), so downstream consumers see no difference in shape. Membership is evaluated from the rule exactly as `fs_betaHhat_one()` does (`get_dfpred()` route), so a rule the frame route cannot evaluate is NA on both routes.
- **Censoring guard.** `"exact"` with `cens_adjust` other than the value meaning no random censoring (read from `simulate_from_dgm()`'s source what that value is, and whether the frame route's default `log(1.5)` adds random censoring) stops with a message saying the exact route covers administrative censoring only. Do not approximate.
- **Entry.** The layer's entry distribution is U(0, `max_entry`) per group; the frame route's draw uses `simulate_from_dgm()`'s `max_entry`. Pass the same value; document it.

Roxygen: what "exact" computes, when it applies, and that the two routes estimate the same population quantity (the frame route with Monte Carlo error from the single evaluation draw).

## 2. MR bounds in `mrct_region_sims(record = "full")` (`R/mrct_simulation.R`)

When `fs_args` carries `mr_inference = TRUE`, the per-trial `results` table under `"full"` gains, NA otherwise: `mr_ci_method`, `mr_debiased_H`, `mr_ij_lower_H`, `mr_ij_upper_H`, `mr_field_lower_1s_H`, `mr_debiased_Hc`, `mr_ij_lower_Hc`, `mr_ij_upper_Hc`, `mr_field_upper_1s_Hc`, `mr_field_upper_1s_s_Hc`, `mr_bonf_lower_H`, `mr_bonf_upper_Hc`, `mr_p_hat_H`, `mr_selection_rate`. Element names as the notes of `TASK_mr-field-rep5-timing_2026-11-01` list them; complement and field columns NA when those options were not used. Under `"point"` nothing changes. The columns come after the existing ones.

## 2b. Two small corrections left by the article task

- `R/forestsearch_main.R`, the Example 1 comment: `mr_estimates_table()` no longer needs a bootstrap template; say it is optional.
- When a template is supplied, the bias-corrected footnote `format_bootstrap_table()` writes reads "Bootstrap-adjusted estimate using infinitesimal jackknife method (draws iterations)" over MR values. Make the footnote name the method actually used (multiplier resampling with the IJ variance, with the number of draws) on the `mr_estimates_table()` template route only; the bootstrap tables keep their footnote. Text only; no number changes.

## 3. Tests

1. **Exact vs frame.** On the S2 build (`helper-mrct-analytic.R`) with a handful of realised rules (the committed replicate-5 Ĥ among them) and administrative censoring only: `method = "exact"` and `method = "frame"` agree within the frame route's Monte Carlo error; quantify that error by running the frame route at three `eval_seed`s and checking the exact value lies within 3 SD of their spread, or, if the package has a cheaper standard for this comparison, use it and say so. Also: `"exact"` with random censoring stops with the documented message; a rule with no members gives NA on both routes.
2. **Exact against the campaign.** For replicate 5 of S2, the exact β(Ĥ) on the non-region patients' population equals `ahr_limit()` called directly on the same membership (identity), and `exp(β(Ĥᶜ))` for the all-patient complement lies where the committed case example puts the population quantities (read `results/s2_example/` in mrct-consistency only if a committed value of this exact quantity exists; otherwise skip this part and say so).
3. **MR columns.** `mrct_region_sims(record = "full", fs_args = list(mr_inference = TRUE, mr_inference_args = list(ci_method = "field", include_complement = TRUE, field_complement = TRUE, draws = 200, field_R_out = 50, field_R_in = 25, seed = 1)))` on 2 replicates, 2 workers: the new columns are filled where declared, NA where not; `record = "point"` with the same `fs_args` returns exactly today's columns; a `"full"` call without `mr_inference` has the columns all NA.

## 4. Gates (stop on failure)

`devtools::document()` clean; `devtools::install()`; export list unchanged; the three new test files/blocks pass; the existing `test-betaHhat-contract.R` and `test-mrct-region-sims-record.R` still pass. No full check.

## 5. NEWS

Add both options to the 0.4.1.9000 entry.

## 6. Notes file

`~/Downloads/NOTES_fs-exact-betaHhat-mr-record_2026-11-01_r2.md`: status, commits, the exact-vs-frame differences with the frame route's SD, test counts, any disagreement with this task.

## 7. Commits (on `mrct-analytic`)

Two: this task; the code, tests, docs and NEWS (§2b included in the second).

## Recorded, not in this task

- The censoring factor in the analytic layer (Ḡ(u) per arm), which would make the exact route general: a derivation, implementation and validation of its own.
- The exact per-candidate effects and SEs, and an exact correlation for overlapping candidate sets, in `fs_oc_family_enumerate()` / `fs_oc_predict()`: a proposal of its own; the current `glm_dgm`-based route is untouched.
