# REPORT — MR (field) on Section 7 cells: Stage 1 (Template, identities, projection)

**Task:** `dev/tasks/TASK_mr_field_section7_2026-09-05.md` (bb516569); Stage 0 record `REPORT_mr_field_s7_stage0_2026-09-05.md` (ca2996b3).
**Date:** 2026-09-05. Executor: Claude Code, under Larry's Gate 0 acknowledgment (findings 1–4 accepted, F1–F5 as ruled, F7 pre-authorization conditional on "every Stage 1 identity passes").

---

## GATE 1: STOP — one identity fails, and its only fix touches the gate's arithmetic

**44 of 48 machine-checked identities pass at exactly 0 relative difference** (not merely within the 1e-12 tolerance). The four failures are one defect seen through four columns: under `ci_method = "field"`, the gate builds the **complement's MR interval from the Wald SE instead of the IJ SE**, so `mr_Hc_lo` / `mr_Hc_hi` differ from the committed IJ bundles by 0.21–0.27 relative on both cells — while `mr_Hc_est` and `mr_Hc_se_ij` are exactly identical.

This is a STOP under two of Larry's rules simultaneously:
- Task 1b requires "MR (IJ) est / SE / **bounds** identical"; Gate 1 requires 1b to pass; the F7 pre-authorization required "every Stage 1 identity passes". Stage 2 is therefore **not started**.
- The one-token candidate fix (below) changes which SE constructs the complement interval — that is the gate's interval arithmetic, and "anything touching the gate's arithmetic is a STOP" (task protocol). No `R/` edit was made.

### Root cause (quoted)

The H block deliberately keeps the IJ interval under `"field"` — `R/fs_mr_inference.R:539–542`:

```r
# "field" keeps the debiased element on the IJ interval (identical to "ij");
# only "wald" switches it to the robust SE.
se    <- if (ci_method == "wald") se_wald else se_ij$se
```

The complement block tests for `"ij"` alone — `R/fs_mr_inference.R:607`:

```r
sec_used <- if (ci_method == "ij") se_ijc$se else sec
```

so under `ci_method = "field"` the complement's `lower/upper/lower_1s` (`:613–615`) are built from `sec` (Wald) while `se_ij` is still stored (`:616`) — exactly the observed signature: est exact, `se_ij` exact, bounds off by the IJ↔Wald SE ratio. The branch predates "field" (it was written when `"ij"`/`"wald"` were the only options) and was not updated when `87880a24` added the third value.

### Candidate fix (NOT applied — Larry's decision)

`R/fs_mr_inference.R:607`: `if (ci_method == "ij")` → `if (ci_method %in% c("ij", "field"))`, aligning the complement with the H block's stated convention ("field keeps the IJ interval"). One token; behavior under `"ij"` and `"wald"` unchanged; under `"field"` it changes only the three complement bound fields. If approved, the standing post-`R/`-edit guards apply (test suite, `fidelity_fs_oc_predict_2026-08-28.R`, `prerefactor_reference_2026-08-29.R check`), then Stage 1's identity render re-runs from scratch. The alternative is to accept complement-on-Wald under field runs and rewrite the 1b expectation — but that leaves the template's Ĥᶜ "MR (IJ)" rows not reproducing the committed IJ intervals, which contradicts decision F3's "complement stays on IJ" as acknowledged.

### Scope of the defect

- **Ĥ block: unaffected.** All H-block MR (IJ) columns (est/SE/bounds) are exactly identical to the committed bundles; the field block itself is finite, internally consistent, and healthy on every detected replicate.
- **Committed Section 7 bundles: unaffected** (they ran under `ci_method = "ij"`, where the branch is correct).
- **The completed MR-field-vs-Guo–He arc (`1f771bf8..8c3fd6c0`): unaffected.** Its adapter (`quarto/GuoHe/mr_vs_guohe_sim.R`) never passes `include_complement` (formal default `FALSE`, `R/fs_mr_inference.R:398`) and its bundles record no complement columns; "complement" appears there only in the V8a return-names contract.
- **Affected:** any future `include_complement = TRUE` run under `ci_method = "field"` — i.e., precisely this task's template.

---

## 1a — Template (delivered, uncommitted pending the STOP resolution)

`quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`, derived from the committed `sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_1_1000.qmd` with knob deltas only (left **untracked** in the working tree per stop semantics — commit is Larry's call with the STOP resolution):

- `ci_method = "field"`; the recorder gains the `fld_H_*` columns (est2; two-sided quantile interval `lo2s/hi2s`; one-sided `lo1s`; SE-type interval `lo_se/hi_se`; `se` = λ-SD; `lam_mean`; the seven λ quantiles; draw-usage `nout`/`nin_mean`; the field's own `secs`; a degenerate-case `note`). Field is recorded for Ĥ only (F3).
- `fb_mode` knob: `"run"` | `"join"` (reads `results/fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds`, joins `fb_*` by `sim_id`, per-row naive-identity assertion at 1e-12, sim_ids outside the bundle recorded FB-absent) | `"none"`.
- `campaign_tag` mandatory in the stem, plus a `fb_mr_field` token and a three-digit h-tag (`h100`/`h175` from `round(100*hr)`) — closing the Gate 0 `h18` collision (finding 2) three ways.
- Every knob overridable per render via `FS_S7_*` environment variables (echoed in the render and stored in `meta`), so all batches run from one committed file.
- Estimation/coverage tables (mean + median companions) gain **MR (field)** rows (est2 with the two-sided quantile interval, F4) for every target, Ĥ block; FB rows render from data presence (run or joined). New `# MR (field) diagnostics` section: two-sided coverage with **Wilson intervals** (primary) and one-sided lower coverage (supplementary) for every estimator × target × block; retained bias of β̃ and est2 vs β(Ĥ) in log-HR and SD units; λ-SD/SD and IJ-SE/SD ratios; draw-usage/cost table. The FB↔MR agreement panel gains an **est2 vs FB** comparison.
- The committed documents are untouched; template purl parses clean (145 top-level expressions).

## 1b — Identities (measured values)

Two smoke renders of the template, sims 1–5 at the committed seeds, production knobs (`mr_draws` 5000, `n` 500), campaign `smoke`, 5 workers: h100 with `fb_mode = join`, h175 with `fb_mode = none`. Compared against `fs_maxeffCons_mr_m1_h10_knoise0_n500_res_1_1000.rds`, `fs_maxeffCons_mr_m1_h18_knoise0_n500_res_1_1000.rds`, `fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds`. Check script: scratchpad `stage1_identity_checks.R`; full console output preserved in this session's log.

| Identity (per cell) | h100 (key) | h175 (separation) |
|---|---|---|
| Rows present, sim_id 1–5 | PASS | PASS |
| Detection flag vs committed | PASS (11001) | PASS (11111) |
| Status trichotomy | PASS | PASS |
| Membership: n_sel, n_harm, n_true, sens, spec, ppv, npv | all 0 | all 0 |
| Naive est/lo/hi/se, both blocks (8 cols) | all 0 | all 0 |
| Oracle est/lo/hi/se, both blocks (8 cols) | all 0 | all 0 |
| MR (IJ) Ĥ: est / lo / hi / se_ij | all 0 | all 0 |
| MR (IJ) Ĥᶜ: est / se_ij | 0 / 0 | 0 / 0 |
| **MR (IJ) Ĥᶜ: lo / hi** | **FAIL 0.259 / 0.206** | **FAIL 0.272 / 0.214** |
| betaHhat_H / betaHhat_Hc | 0 / 0 | 0 / 0 |
| MR (field) finite on all detected reps (3 / 5 detections) | PASS | PASS |
| Field identities: lo2s = exp(β̃−q975), hi2s = exp(β̃−q025), lo1s = exp(β̃−q95), est2 = exp(β̃−λ̄), lo_se = exp(log est2 − 1.96·λsd) | ≤ 3.3e-16 | ≤ 1.3e-16 |
| 2-sided / 1-sided field coverage indicators recomputed from stored bounds | PASS (1.000 / 0.667) | PASS (0.800 / 1.000) |
| FB join: all 5 sim_ids matched; 8 fb_* columns = committed FB bundle; naive cross-bundle identity | all 0 | n/a (no FB) |

`sg_def` strings: identical on h100; spelling-different on h175 (expected 0.2.0-vintage drift; memberships exact — the agreed comparison).

## 1c — Projection (~60 workers), indicative

Smoke field cost (5 workers, lightly loaded): **8.2–9.0 s per detected replicate**, with healthy draw usage (n_out 974–1000 of 1000; n_in mean 496–500 of 500); fit+MR 19.5–23.1 s. Anchors at load: committed 115-worker runs measured fit+MR (IJ) at mean 28.3 s (h10) / 35.0 s (h175) per replicate. Allowing ×2–4 memory-bandwidth inflation on the field's matmul (host is bandwidth-bound), per-replicate ≈ 45–70 s; at 60 workers a 1,000-replicate batch ≈ ceil(1000/60) × per-rep ≈ **13–20 min**, so 2 cells × 2 batches + combines ≈ **1–1.5 h wall, ~50–80 core-h** — comfortably inside the 12 h go-condition. Cost was never the blocker; the identity was.

---

## State left behind (all untracked; no `R/` edits; nothing committed but this report)

- Template: `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`
- Smoke evidence bundles: `results/fs_maxeffCons_fb_mr_field_m1_h100_knoise0_n500_smoke_res_1_5.rds`, `results/fs_maxeffCons_fb_mr_field_m1_h175_knoise0_n500_smoke_res_1_5.rds`
- Smoke renders: `smoke_field_h100.html`, `smoke_field_h175.html` (same directory)

## Deviations

- None from the Gate 0-acknowledged plan up to the STOP; Stage 2/3 not entered (F7 condition unmet). The template commit, prescribed by 1a, is deferred to the STOP resolution because stop semantics commit only the report.

**Stopped at Gate 1.** Decision needed from Larry: approve the one-token complement-SE fix in `R/fs_mr_inference.R:607` (then guards + full Stage 1 re-run), or accept complement-on-Wald under field with a rewritten 1b expectation.
