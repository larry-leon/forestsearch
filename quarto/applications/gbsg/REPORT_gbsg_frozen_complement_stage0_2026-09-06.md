# REPORT — GBSG frozen-family illustration, complement rows and the joint pair: Stage 0 (Discovery)

**Task:** `dev/tasks/TASK_gbsg_frozen_complement_2026-09-06.md` (eac66f3b). L1–L4 at defaults. Sequencing satisfied: `TASK_mr_field_complement` (Stage 3 at 815debf1) and `TASK_complement_refinements` (Stage 3 at e2e0d747) both complete; the installed package is 8fbb3bdc's (`field_complement`, `ij_residual`, `field$joint` all present).
**Date:** 2026-09-06. No compute; no `R/` changes.

## Element names (verified from source, `R/fs_mr_inference.R` at 8fbb3bdc)

- `fs$mr_inference$complement`: `naive$est/lower/upper`; `debiased$est/lower/upper/lower_1s/se/se_ij/se_wald/var_ij/ij_source/ij_draws` plus the winner-only elements `se_ij_two_term`, `se_ij_winner`, `ij_source_winner`, `se_ij_winner_floor`, `ij_source_winner_floor`, `lower_w/upper_w/lower_1s_w/upper_1s_w`, `lower_wf/upper_wf/lower_1s_wf/upper_1s_wf`; `selection_bias`, `fixed_bias`, `n`. (Harm-side `debiased` carries the same winner elements without `upper_1s_w/_wf`.)
- `fs$mr_inference$field$complement` (needs `field_complement = TRUE`): `lambda_mean`, `lambda_sd`, `q05…q975`, `n_in_used_mean`, `est2`, **`upper_1s`** (primary), `lower_1s`, `lower_2s/upper_2s`, `se_field`, `lower_se/upper_se`, `n_out_used`, `n_out_dropped_unfit`, `n_complement_fits`, `n_new_fits`, `share_draws_new_fit`, `R_out`, `R_in`, `timing_seconds`.
- `fs$mr_inference$field$joint`: `gamma`, `joint_prob`, `alpha`, `lower_H`, `upper_Hc`, `bonf_gamma`, `bonf_lower_H`, `bonf_upper_Hc`, `bonf_joint_prob`, `corr`, `n_joint_draws`, `grid_gamma`, `grid_joint_prob`.
- `fs$mr_inference$ij_residual` (top level; `"two_term"` here).

## The document and its committed payload (`quarto/applications/gbsg/`)

- **Gate call** (`analysis_gbsg_survival_frozen_family.qmd:284–286`): `mr_inference_args = list(draws = mr_draws, ci_method = "field", return_reselection = TRUE, field_uniform = TRUE, field_M_cap = 1000L)`. `include_complement` is not stated and resolves to `TRUE` through `forestsearch()`'s `.g_mr(mr_inference_args$include_complement, TRUE)` (`R/forestsearch_main.R:3392`), so stating it is a no-op and **`field_complement = TRUE` is the only behavioural change**; `ij_residual` stays at its default, so every reported harm-side number is computed exactly as before (the refinements task's J1 identities).
- **FB complement:** `fs_bc$Hc_estimates` exists with `H0/sdH0/H0_lower/H0_upper`, `H1/…`, `H2/sdH2/H2_lower/H2_upper` (the committed payload stores it under `extras$fb$Hc_estimates`) → the FB row of Table 3 is available under L1 (`H2` = the two-source bias-corrected complement estimate; one-sided upper `exp(log H2 + 1.645·sdH2/H2)`, the delta-method convention Table 1 uses).
- **Committed payload (identity anchor):** `_payloads/analysis_gbsg_survival_frozen_family/analysis_gbsg_survival_frozen_family_payload.rds` (ca2a2f93/b94478c9, 2026-09-06 10:37), `extras$mr = {naive, debiased, settings, n_family, selection_bias, timing_seconds, field (with uniform), reselection$p_hat}`, `extras$fb = {H_estimates, Hc_estimates}`, `extras$gh = {naive_hr, debiased_hr, bound_hr, selected, n_family, B, r}`, `labels$sg_harm = "{pgr <= 32.5} {er <= 0}"`; harm-side anchors: MR (IJ) est 1.753, field est2 1.79 / lower_1s 0.88, κ* 1.59, FB H2 1.96, G&H 1.36 / 0.96.
- **Payload location:** `_payloads/<dirout>/<dirout>_payload.rds` with `dirout <- NULL → the .qmd stem` (`:79–80`, `:905–911`). Rendering as-is would overwrite the committed payload above. To keep committed payloads read-only, this render sets `dirout <- "analysis_gbsg_survival_frozen_family_complement"` (a document knob, documented in place), so the new payload lands beside the committed one; the tracked HTML is refreshed per the directory's convention.
- **LOO cache:** `loo_cache <- Sys.getenv("LOO_CACHE", unset = file.path(gh_dir, "cv_out_<doc>_<focus>_<rule>.rds"))` (`:576`); the committed cache `_payloads_2026-09-01_complete/analysis_gbsg_survival_frozen_family/cv_out_analysis_gbsg_survival_frozen_family_maxeff_neighborhood.rds` (the one the previous render used) is passed through `LOO_CACHE`. Confirmed present (1,190 bytes, 2026-09-01).
- **Timing table** (`:868–882`): rows Fit + MR gate, FB, MR gate, field block, uniform sweep, G&H, 10-fold CV, LOO; gains "MR complement field block (within gate)".
- **`@sec-intervals`** (`:722–866`): Table 1 (`intervals-table`), Table 2 (`reselection-table`), the callout with the generated reading (`intervals-reading`). Tables 3–4 and the complement reading go after the callout, before `timing`.

## Plan (Stage 1)

Edit per the task's Changes 1–4; render with `LOO_CACHE` set; identities against the committed payload (naive / FB `H_estimates` and `Hc_estimates` / MR (IJ) / MR field harm block / κ / G&H / selected subgroup, ≤ 1e-12); complement blocks finite, bound identities, `gamma ∈ [0.025, 0.05]`; STOP on any harm-side difference.
