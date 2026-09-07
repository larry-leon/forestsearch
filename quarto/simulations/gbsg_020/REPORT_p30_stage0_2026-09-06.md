# REPORT — Prevalence ~30%: Stage 0 (Discovery) — stopped at Gate 0 for Larry's review

**Task:** `dev/tasks/TASK_prevalence30_2026-09-06.md` (c17e31fb). H-P1–H-P4 at defaults; winner-only / winner-floor variants excluded from every table, figure and report line (H-P3).
**Date:** 2026-09-06. No `R/` changes; no template edits yet (Stage 1). The only compute was the DGM builds and calibrations quoted in 0b (~1 min total). Source: `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` at 37f38540, `R/sim_aft_gbsg.R`, `R/setup_gbsg_dgm.R`, `R/oc_analyses.R`, `R/betaHhat_truth.R`.

---

## GATE 0: PASS on both stop conditions — the truth attachment follows the rule, and a same-form rule reaches [0.28, 0.32]. Stopped here per instruction: the proposed rule, its realized prevalence and the calibrated targets are below for approval before Stage 1.

## 0a — How the template defines the harm subgroup, its prevalence, the calibration and the truth targets

**The rule (template `@sec-truth`, `:78–83`):** "$Q = \{er \le 8\} \cap \{meno = 0\}$ — low estrogen receptor **and** premenopausal. The `er` threshold is the 25th percentile of the GBSG covariate template … The region covers **12.2%** of the template (84 of 686 patients), and the DGM plants the harm effect there via `setup_gbsg_dgm()`, which encodes it as `z1 = 1 & z3 = 1`."

**Where the rule lives (engine `.create_gbsg_dgm_()`, `R/sim_aft_gbsg.R:279–290`):**

```r
er_threshold <- stats::quantile(dfa$er, probs = z1_quantile)
dfa$z1 <- ifelse(dfa$er <= er_threshold, 1L, 0L)
dfa$z3 <- ifelse(dfa$meno == 0, 1L, 0L)  # Premenopausal
dfa$zh <- dfa$treat * dfa$z1 * dfa$z3
dfa$flag.harm <- ifelse(dfa$z1 == 1 & dfa$z3 == 1, 1L, 0L)
```

`z1_quantile` is a formal of both `setup_gbsg_dgm()` (`R/setup_gbsg_dgm.R:93`, default `0.25`) and the engine (`:241`, default `0.25`); the threshold is the `z1_quantile` quantile of `er` on the 686-row GBSG template, and prevalence ≈ 12.2% *arises* from P(er ≤ 8 & meno = 0) on that template (12.42% on the 100k super-population, which resamples the template). The prevalence is therefore a consequence of the rule, not a separate parameter, and the rule's form — a single `er` cut in conjunction with the raw `meno` binary — is fixed by the engine.

**How `k_inter` acts (`:416–420`, documented at `:148–159`):** the engine fits an AFT model to GBSG with `zh` as a covariate and then scales the fitted interaction, `gamma["zh"] <- k_inter * gamma["zh"]`; "HR(H) = exp(−gamma[treat]/sigma − gamma[zh]/sigma)". So `k_inter` is a *multiplier of a data-fitted coefficient whose value depends on the rule* (see 0b).

**The template's DGM build (`:571–576`):**

```r
k_inter <- calibrate_k_inter(target_hr_harm = target_hr_harm, model = dgm_model, use_ahr = FALSE)
dgm <- setup_gbsg_dgm(model = dgm_model, k_inter = k_inter, n_super = n_super, seed = seed_base)
dgm <- compute_dgm_cde(dgm)   # attaches CDE (theta-ddagger) targets
```

`calibrate_k_inter()` (`R/sim_aft_gbsg.R:1004–1050`) root-finds `k_inter` so that the engine's `hr_H_true` (the Cox HR in `flag.harm == 1` on the super-population) equals `target_hr_harm`, and **forwards `...` to `.create_gbsg_dgm_()`**, so a `z1_quantile` passed to it is honoured — the calibration follows whatever rule the document passes, provided the same value is passed to both calls.

**Truth attachment follows the rule (Gate 0 condition 1):**
- θ† (marginal): `truth$marg_H/marg_Hc <- dgm$hr_H_true / hr_Hc_true` (`:581–582`), computed by the engine from `flag.harm` (`R/generate_aft_dgm_helpers.R:1112–1138`, Cox fits within `flag_harm == 1 / 0`).
- θ‡ (CDE): `compute_dgm_cde()` (`R/oc_analyses.R:131–165`) auto-detects `flag.harm`/`flag_harm` in `df_super` and takes `mean(exp(theta_1[in_H])) / mean(exp(theta_0[in_H]))` and its complement — rule-agnostic.
- β(Ĥ)/β(Ĥᶜ): `fs_build_eval_frame(dgm, …)` (`R/betaHhat_truth.R:639`) evaluates the full super-population once; `fs_attach_betaHhat(results, eval_df, focus = "harm", …)` (template `:1174`) evaluates each replicate's *realized* rule string on that frame; the θ†-on-frame cross-check reads `harm.name = "flag_harm"` (`:715–728`). Nothing hard-codes `er ≤ 8`: the per-subject `flag_harm` carried in `df_super` is the only channel, and it is set by the rule.
- Search side: `fs_conf_force <- c("meno == 0", "er <= 0", "pgr <= 0")` and `conf.cont_jcuts = list(er = 10)` (`:465–466`) — the forced cuts and the J-quantile `er` grid `0, 3, 9, 17, 30, 44, 69.8, 100, 173.5, 293.7` (`:97–99`) are rule-independent; the design property "the true cut is not on the grid" must be re-checked for the new threshold (it holds, 0b).

## 0b — Proposed ~30% rule (same form), realized prevalence, calibrated `k_inter`, truth targets

**Rule:** the same form, `{er ≤ q_{0.60}(er)} ∩ {meno = 0}`, i.e. `z1_quantile = 0.60` → **er ≤ 59** (the 60th percentile of the GBSG template; er median 36). Realized prevalence on the 100k super-population (seed 8316951): **0.3065**. The grid over `z1_quantile` on the same super-population, for the record:

| `z1_quantile` | er threshold | P(er ≤ thr) | prevalence of `er ≤ thr & meno = 0` |
|---|---|---|---|
| 0.25 (M1) | 8 | 0.264 | 0.1242 |
| 0.50 | 36 | 0.502 | 0.2542 |
| 0.55 | 45 | 0.554 | 0.2915 |
| **0.60** | **59** | 0.602 | **0.3065** |
| 0.65 | 75 | 0.654 | 0.3313 |
| 0.70 | 89 | 0.703 | 0.3472 |

P(meno = 0) = 0.420 caps the conjunction at 42%; 0.55 (0.2915) and 0.60 (0.3065) are both inside [0.28, 0.32]; 0.60 is proposed as the nearer to 30% and the further from the 0.28 edge under trial-level sampling variation (at n = 500 the realized prevalence has SD ≈ 0.02). The threshold 59 sits between the grid's 44 and 69.8, so the exact boundary is again **not** a candidate cut — the design property is preserved (the nearest grid cut, 69.8, over-covers; 44 under-covers).

**Calibrated `k_inter` and truth targets at `z1_quantile = 0.60`** (`calibrate_k_inter(target_hr_harm, model = "alt", use_ahr = FALSE, z1_quantile = 0.60)`, then `setup_gbsg_dgm(…, z1_quantile = 0.60, n_super = 100000L, seed = 8316951L)` + `compute_dgm_cde()`), with the M1 values beside them:

| Cell target HR(Ĥ) | `k_inter` (0.60) | θ†(Ĥ) | θ‡(Ĥ) = AHR(Ĥ) | θ†(Ĥᶜ) | θ‡(Ĥᶜ) | overall causal HR | M1 (0.25): `k_inter`; θ†(Ĥ)/θ†(Ĥᶜ); θ‡(Ĥ)/θ‡(Ĥᶜ); overall |
|---|---|---|---|---|---|---|---|
| 1.00 | −8.579 | 0.9999 | 1.0000 | **0.7206** | **0.6564** | 0.7922 | 0.568; 1.0005/0.6569; 1.0000/0.5848; 0.6847 |
| 1.50 | −19.497 | 1.4990 | 1.7087 | 0.7206 | 0.6564 | 0.8694 | 1.112; 1.5086/0.6569; 1.6713/0.5848; 0.7041 |
| 1.75 | −23.583 | 1.7462 | 2.0881 | 0.7206 | 0.6564 | 0.8944 | 1.320; 1.7691/0.6569; 2.0355/0.5848; (0.72) |

Every target is hit within the calibrator's tolerance (θ†(Ĥ) 0.9999 / 1.4990 / 1.7462). Cell 4 (HR 1.00, n = 1000) shares cell 1's DGM.

**Two things Larry should see before approving (not stop conditions):**

1. **The `k_inter` values are large and negative because the data-fitted interaction changes with the rule.** `k_inter` multiplies the GBSG-fitted AFT coefficient of `zh`; at `z1_quantile = 0.25` that coefficient is −0.674 (log-time; a genuine harm interaction in the ER-low premenopausal corner of the real data), so M1 needs `k_inter` ∈ [0.57, 1.32]. At 0.60 the fitted `zh` coefficient is **+0.035** — essentially zero and of the opposite sign (the real GBSG data show no interaction in the wider ER ≤ 59 region) — so harm must be manufactured by a large negative multiplier. The *effective* planted interaction, `k_inter × gamma_zh`, is of the same magnitude as M1's: −0.30 / −0.68 / −0.83 at 0.60 vs −0.38 / −0.75 / −0.89 at 0.25 (log-time units; sigma 0.713 in both). The DGM is well-defined and the targets reproduce; the point is that the 30% harm subgroup is no longer "where the data say harm is", only where it is planted. The Cox warnings ("coefficient may be infinite") appeared during the root search's bracket extension at extreme trial values, not at the solution.
2. **The complement's benefit is weaker under the 0.60 rule:** θ†(Ĥᶜ) = 0.721 and θ‡(Ĥᶜ) = 0.656 (vs 0.657 / 0.585 at M1), and the fitted treatment main effect shifts (AFT `gamma[treat]` 0.300 vs 0.383) because the AFT refit absorbs the different `zh`. The complement is also smaller (69% vs 88% of the trial). Both are consequences of the same-form rule on the GBSG covariate structure, not choices — but they change what "a benefit of at most U" means in the p30 cells, and the reading aids (0.80 / 0.85) sit closer to θ†(Ĥᶜ) = 0.72 than they did to 0.66.

If Larry prefers the harm region to remain one the real data support, that is a change of form (H-P2) — e.g. a different covariate pair — and is his choice at this gate; the same-form default is deliverable as proposed.

## 0c — The knob (proposed lines, add-only, default reproduces M1 exactly)

In the DGM knob block (after `n_super`, template `:403`):

```r
# ── Harm-subgroup prevalence knob (TASK_prevalence30_2026-09-06) ─────────────
# z1_quantile is the engine's own formal: the harm rule is
#   {er <= quantile(er, z1_quantile)} & {meno == 0}
# (the same form at every value).  0.25 is M1 (er <= 8, prevalence 12.4%) and
# the package default, so a render without FS_S7_Z1Q is byte-identical to the
# committed s7/s7c/map1/map1c bundles; 0.60 is the ~30% cell (er <= 59,
# prevalence 0.3065 on the super-population).  Passed to BOTH the calibration
# and the DGM build, so k_inter is calibrated under the same rule it plants.
harm_z1_quantile <- .env_num("FS_S7_Z1Q", 0.25)
stopifnot(harm_z1_quantile > 0, harm_z1_quantile < 1)
```

In `build-dgm` (`:571–575`):

```r
k_inter <- calibrate_k_inter(target_hr_harm = target_hr_harm,
                             model = dgm_model, use_ahr = FALSE,
                             z1_quantile = harm_z1_quantile)
dgm <- setup_gbsg_dgm(model = dgm_model, k_inter = k_inter,
                      z1_quantile = harm_z1_quantile,
                      n_super = n_super, seed = seed_base)
```

Byte-identity argument: `calibrate_k_inter()` passes `z1_quantile` through `...` to `.create_gbsg_dgm_()`, whose default is 0.25; `setup_gbsg_dgm()`'s own default is 0.25. Passing 0.25 explicitly is the same call as omitting it (`match.arg`-free numeric formal), so the DGM object — and every downstream column — is identical; Stage 1's knob-inert identity against the s7c bundle (5 replicates, every column ≤ 1e-12) is the machine check.

Stem, meta and pooling (`:350–352`, `:1280`, `:1317–1318`): the stem gains a prevalence tag only when the knob is non-default — `..._n%d[_z1q%02d]_%s` with `round(100 * harm_z1_quantile)` (so p30 stems read `..._n500_z1q60_p30`; M1 stems are unchanged, preserving every committed stem); `meta` gains `harm_z1_quantile` and `harm_prevalence_super` (= `mean(dgm$df_super$flag_harm)`); `harm_z1_quantile` joins the combine poolability gate (a 0.25 batch and a 0.60 batch must never pool). The knobs echo line gains `z1q=%.2f`. The realized per-replicate prevalence is already recorded (`n_true / n_sample`); the truth block records `prevalence`.

**Batch-save guard (task protocol; the cell-7 incident), proposed lines before `saveRDS(…, rds_path)` at `:1250` and before the pooled save at `:1364`:**

```r
# Refuse to overwrite a git-tracked bundle (committed bundles are read-only;
# guard from the cell-7 campaign-tag incident, REPORT_complement_refinements_gate2).
.refuse_if_tracked <- function(path) {
  tracked <- tryCatch(
    system2("git", c("ls-files", "--error-unmatch", shQuote(path)),
            stdout = FALSE, stderr = FALSE) == 0L,
    error = function(e) FALSE)
  if (isTRUE(tracked))
    stop("refusing to overwrite the git-tracked bundle '", path,
         "': committed bundles are read-only -- use a new campaign_tag.",
         call. = FALSE)
  invisible(path)
}
.refuse_if_tracked(rds_path)
```

(`system2()` returns the exit status when `stdout = FALSE`; `git ls-files --error-unmatch` exits 1 for an untracked path and non-zero outside a repository, both of which fall through to "not tracked". Add-only; a batch whose path is untracked behaves exactly as today. The `FS_S7_SAVE_COMBINED` guard already protects committed pooled bundles from combine renders; this one closes the batch path.)

## 0d — Cost anchors (`REPORT_mr_field_complement_gate2_2026-09-06.md`, `REPORT_complement_refinements_gate2_2026-09-06.md`)

At 100 workers with the complement field on: h100 n500 16 min, h150 n500 19 min, h175 n500 19–20 min, h100 n1000 26 min per 2,000 replicates; fit+MR 39–46 s per replicate at n = 500 (complement block 2.8–3.0 s, 350–600 complement fits per replicate), 80 s at n = 1000. At 30% prevalence the harm subgroup is ~150 of 500, so the candidate family is larger and complements differ more (more distinct complement fits — the task's 1c expectation); the per-replicate cost may rise by the complement-fit share (a few seconds). **Projection for the four H-P1 cells ≈ 16 + 19 + 20 + 26 ≈ 1.4 h, allowing up to ~2 h** — within the ~2.5 h expectation stated under H-P5.

## What is needed from Larry at this gate

- Approve the same-form rule `z1_quantile = 0.60` (er ≤ 59, prevalence 0.3065) with the calibrated targets above — or choose a different form (H-P2) given item 1 of 0b (the harm region at 30% is planted, not data-supported) and item 2 (θ†(Ĥᶜ) = 0.72).
- Approve the knob name `FS_S7_Z1Q` / variable `harm_z1_quantile`, the conditional stem tag `_z1q60`, and the batch-save guard as quoted, for Stage 1.
