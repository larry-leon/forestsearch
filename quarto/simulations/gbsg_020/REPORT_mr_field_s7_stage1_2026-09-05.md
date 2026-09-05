# REPORT — MR (field) on Section 7 cells: Stage 1 (Template, identities, projection) — PASS after fix

**Task:** `dev/tasks/TASK_mr_field_section7_2026-09-05.md` (bb516569); Stage 0 record ca2996b3.
**Date:** 2026-09-05. Executor: Claude Code.
**History:** the first Stage 1 pass ended in a GATE 1 STOP (committed 8293a3ce; see that revision of this file for the full failure record): under `ci_method = "field"` the complement's MR interval was built from the Wald SE instead of the IJ SE (`R/fs_mr_inference.R:607` tested `ci_method == "ij"` only), failing 4 of 48 identities (`mr_Hc_lo`/`mr_Hc_hi`, 0.21–0.27 rel, both cells) while `mr_Hc_est`/`se_ij` were exact. Larry resolved the STOP by approving the one-token fix.

---

## GATE 1: PASS — fix applied and guarded, 48/48 identities at zero relative difference, projection far under the 12 h ceiling. Stage 2 started under the F7 pre-authorization.

## The fix (commit a6702fd8)

`R/fs_mr_inference.R:607`: `if (ci_method == "ij")` → `if (ci_method %in% c("ij", "field"))`.
**Classification (Larry's, for the record):** R/ change on the `"field"` path only — complement interval Wald → IJ, restoring the documented intent that the de-biased element is IJ-computed under `"field"`; default `"ij"` path untouched; no method change. Roxygen `@return` gained one clarifying sentence (the complement CI follows the winner's SE convention); `man/fs_mr_inference.Rd` regenerated.

### Post-R/ guards (all green before the fix was committed)

- **Default-path byte-identity, three Guo–He Stage 1 fixed-seed cases** (F0 of `quarto/GuoHe/mr_field_stage1_checks.R`, run against the pre-change 736088e3 reference objects): Sec 5.1 `t35_beta2_00` m=1, Sec 5.2 `t7_beta2_00` m=1, and K=1 — all `identical()` component-for-component, `timing_seconds` excluded (the asserted sole exclusion). Note this reference predates the *field feature itself*, so the chain also re-confirms the feature's no-op on `"ij"`.
- **Complement guard (added):** `ij` and `wald` runs with `include_complement = TRUE` (fixed seed 424242, K=1 case) captured with the pre-fix installed package, then re-run post-fix: both `identical()` (timing excluded). **Positive control:** the `field` run's complement `lower/upper/lower_1s/se` now equal the `ij` run's exactly.
- **Test suite:** 0 fail / 4941 pass / 3 skip / 32 warn — the standing 0.3.5 parity.
- Installed via `R CMD INSTALL --preclean` (workers use the installed package).

## 1a — Template (committed with this report)

`quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`, derived from the committed `sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_1_1000.qmd` with knob deltas only:

- `ci_method = "field"`; recorder gains `fld_H_*`: `est2`, two-sided quantile interval (`lo2s`/`hi2s`, primary per F4), one-sided `lo1s`, SE-type interval (`lo_se`/`hi_se`), `se` (= λ-SD), `lam_mean`, seven λ quantiles, draw-usage (`nout`, `nin_mean`), the field's own `secs`, a degenerate-case `note`. Ĥ only (F3; Ĥᶜ stays on IJ).
- `fb_mode` ∈ {`run`, `join`, `none`}; `join` reads `results/fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds`, joins `fb_*` by `sim_id` with a per-row naive-identity assertion (1e-12), and records sim_ids outside the bundle as FB-absent.
- Collision-guarded stem (Gate 0 finding 2): `fs_maxeffCons_fb_mr_field_m1_h<3-digit>_knoise<k>_n<n>_<campaign_tag>`; three-digit h-tag from `round(100·hr)` and a mandatory `campaign_tag`, both also gating combine-mode poolability.
- All knobs overridable per render via `FS_S7_*` env vars (echoed in the render, stored in `meta`) — every batch runs from this one committed file.
- Tables (mean + median) gain **MR (field)** rows for every target (Ĥ block); FB rows render from data presence; new `# MR (field) diagnostics` section (Wilson-interval two-sided coverage, one-sided supplementary coverage, retained bias in log-HR and SD units, λ-SD/SD and IJ-SE/SD, draw-usage/cost); agreement panel gains **est2 vs FB**.

## 1b — Identities after the fix: 48/48 PASS

Smoke renders re-run from scratch under the fixed installed package (sims 1–5, committed seeds, production knobs, campaign `smoke`; h100 `fb_mode=join`, h175 `fb_mode=none`). Every check from the STOP report's table now passes, **all numeric identities at exactly 0 relative difference** (tolerance was 1e-12), including the previously failing `mr_Hc_lo`/`mr_Hc_hi` on both cells:

- Detection flags (h100: 11001; h175: 11111), status, memberships (n_sel, n_harm, n_true, sens, spec, ppv, npv): exact.
- Naive and oracle est/lo/hi/se, both blocks: exact. MR (IJ) est/lo/hi/se_ij, both blocks: exact. β(Ĥ)/β(Ĥᶜ): exact.
- MR (field): finite on all detected reps (3 of 5 / 5 of 5), no degenerate notes; bound↔λ-quantile identities (`lo2s = exp(β̃−q975)`, `hi2s = exp(β̃−q025)`, `lo1s = exp(β̃−q95)`, `est2 = exp(β̃−λ̄)`, `lo_se = exp(log est2 − 1.96·λsd)`) ≤ 3.3e-16; two-sided and one-sided coverage indicators recomputed from stored bounds agree with the document's `.cover` logic (h100: 1.000 / 0.667; h175: 0.800 / 1.000 at 3 / 5 detections — indicative only at this n).
- FB join (h100): all 5 sim_ids matched, all 8 `fb_*` columns exactly equal to the committed FB bundle, naive cross-bundle identity exact.
- `sg_def` strings: identical on h100, spelling-drift on h175 with exact memberships (the agreed comparison).

Smoke evidence bundles committed beside the results: `..._h100_knoise0_n500_smoke_res_1_5.rds`, `..._h175_knoise0_n500_smoke_res_1_5.rds`. Check script: session scratchpad `stage1_identity_checks.R` (full 48-line console record in the session log).

## 1c — Projection and worker decision (measured at target concurrency)

Larry asked for a ~20-replicate burst at 60 and 100 workers; a 20-replicate burst occupies only 20 workers and cannot separate the two settings, so each burst was sized to **two full worker-waves** (deviation noted): h175 cell (the heavier one), `fb_mode = none`, campaigns `burst60`/`burst100`, sim_id from 1.

| Setting | Reps | Run-loop wall | Throughput | fit+MR /rep (incl. field) | field /rep | Peak host mem (baseline ~35 GB) |
|---|---|---|---|---|---|---|
| 60 workers | 120 | 74.2 s | **97 reps/min** | mean 28.3 s | 13.2 s | 81.5 GB |
| 100 workers | 200 | 110.4 s | **109 reps/min** | mean 38.9 s (q90 50.1) | 15.3 s | 108.2 GB |

Memory budget: 70% of 251 GB ≈ 176 GB; both settings are far inside it. **Decision rule → Stage 2 at 100 workers** (higher throughput, ~12%). Field Monte Carlo stayed healthy under load (min n_out 766/1000; zero NA `est2` on 188 detections). Burst bundles committed as the measurement record.

**Projection (100 workers):** 1,000-replicate batch ≈ 1000 / 109 ≈ 9.2 min run-loop + ~1 min fixed overhead; 4 batches + 2 combine renders ≈ **~45 min wall, ~45 core-h** for 2 cells × 2,000 replicates. Against the F7 conditions: identities 48/48 and projection ≪ 12 h → **Stage 2 started immediately** under the 14 h hard timeout, no mid-run changes, fail-fast per batch.

## Deviations

- Burst size 120/200 instead of "~20" (reasoned above; both measurements recorded).
- The complement-guard reference objects were captured minutes before the fix with the pre-fix installed package (the original Guo–He F0 reference predates the field feature and cannot exercise `include_complement`); both are session-scratchpad artifacts quoted here.
