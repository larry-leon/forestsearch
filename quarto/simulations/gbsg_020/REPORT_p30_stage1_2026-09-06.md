# REPORT — Prevalence ~30%: Stage 1 (Knob, guard, identities, smoke, projection)

**Task:** `dev/tasks/TASK_prevalence30_2026-09-06.md` (c17e31fb); Stage 0 record a969c56a (Gate 0 PASS; H-P2 approved by Larry as proposed: `{er ≤ q₀.₆₀(er)} ∩ {meno = 0}`, `harm_z1_quantile = 0.60`, prevalence 0.3065; caution 1 recorded in the template comment: the 30% region is planted rather than data-supported, with interaction sizes comparable to M1's).
**Date:** 2026-09-06. Executor: Claude Code, unattended per Larry's pre-authorization. H-P1–H-P4 at defaults; winner variants excluded (H-P3). No `R/` changes.

---

## GATE 1: PASS — knob-inert identity exact (every column and the truth object, vs s7w and s7c), all four 30% smokes pass, suite green (0 fail / 5025 pass / 3 skip / 32 warn), projection ≈ 1.4 h at 100 workers (< 3 h ceiling, 4 h hard timeout). Compute go: all four H-P1 cells run, none deferred.

## 1a — Template changes (document-level, add-only; `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`)

- **Knob** `harm_z1_quantile <- .env_num("FS_S7_Z1Q", 0.25)` (placed beside `target_hr_harm`/`n_sample`, before the stem), passed to **both** `calibrate_k_inter(…, z1_quantile = harm_z1_quantile)` and `setup_gbsg_dgm(…, z1_quantile = harm_z1_quantile)`; `harm_prevalence_super <- mean(dgm$df_super$flag_harm)` printed and stored; the knobs echo line gains `z1q=`.
- **Stem tag** `_z1q60` only when the knob is non-default (`z1q_tag`), so no committed M1 stem changes; `meta` gains `harm_z1_quantile` and `harm_prevalence_super` (carried into the pooled meta); `harm_z1_quantile` joins the combine poolability gate.
- **Save guard** `.refuse_if_tracked()` (`git ls-files --error-unmatch`; untracked and out-of-repo both fall through) before the batch save and before the pooled save. Unit check from the results directory: the tracked `…_s7c_res_1_1000.rds` → "refusing to overwrite the git-tracked bundle"; an untracked path → passes.
- **Winner rows** (H-P3): `show_winner_rows <- identical(.env_chr("FS_S7_WINNER_ROWS", "FALSE"), "TRUE")` gates the MR (IJ, winner) / (IJ, winner-floor) rows out of every estimation and coverage table and both displays by default; their recorder columns are still written (the task permits) but never summarised.

## 1b — Identities and smoke (all PASS; `p30_smoke_check.R`, session scratchpad)

**Knob-inert identity** (campaign `z1qinert`: h100 n500, sims 1–5 at the committed seeds, default knob, `FS_S7_FIELD_COMPLEMENT=TRUE`, FB none): all **114** shared columns identical to the committed `s7w` bundle and all **93** to `s7c` (worst relative difference 0.0; `fb_*`/wall-clock excluded), the `truth` object `identical()`; meta records `harm_z1_quantile = 0.25`, `harm_prevalence_super = 0.1242`; stem untagged. The knob is inert at its default by construction (explicit 0.25 = the engine's and wrapper's defaults) and now by machine check.

**30% smoke** (campaign `p30smoke`, `FS_S7_Z1Q=0.60`, 5 replicates per H-P1 cell at the committed seeds, 5 workers, rendered sequentially):

| Cell | realized prevalence per replicate (n_true/n) | truth (θ†H/θ†Hc; θ‡H/θ‡Hc) | detected | fields finite | γ in range | bound ids | β(Ĥ)/β(Ĥᶜ) attached | complement fits (share new) | mean |Ĥ| | fit+MR s/rep (light) |
|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.00 n500 | 0.308, 0.344, 0.324, 0.342, 0.310 | 1.000/0.721; 1.000/0.656 | 4/5 | yes | yes | 5.6e-17 | 4/4 | 428 (0.061) | 66 | 20.3 |
| HR 1.50 n500 | same draws | 1.499/0.721; 1.709/0.656 | 5/5 | yes | yes | 1.1e-16 | 5/5 | 390 (0.101) | 70 | 27.1 |
| HR 1.75 n500 | same draws | 1.746/0.721; 2.088/0.656 | 5/5 | yes | yes | 1.1e-16 | 5/5 | 372 (0.109) | 70 | 28.2 |
| HR 1.00 n1000 | 0.304, 0.316, 0.313, 0.328, 0.300 | 1.000/0.721; 1.000/0.656 | 5/5 | yes | yes | 1.1e-16 | 5/5 | 414 (0.034) | 125 | 34.1 |

All in [0.24, 0.36]; super-population prevalence 0.3065 in every bundle's meta; `k_inter` calibrated under the 0.60 rule (targets reproduced to the third decimal). One early observation for Stage 3, not a check: at n = 500 the identified Ĥ averages 66–70 subjects against a true harm subgroup of ~155 — the search returns a sub-region of the planted 30% (the `er ≤ 59` boundary is not on its grid; 44 under-covers), so β(Ĥ) and the realized rules will differ from the truth targets more than at M1. Smoke bundles and HTMLs committed as evidence (`…_z1qinert_res_1_5.rds`, `…_z1q60_p30smoke_res_1_5.rds` × 4).

**Suite:** 0 fail / 5025 pass / 3 skip / 32 warn (no `R/` change; run for the record).

## 1c — Projection (100 workers)

Light-load fit+MR per replicate (20–34 s) is at or under the M1 smokes' (23–25 s at n = 500); complement fits per replicate (370–430) and the new-fit share (3–11%) are in the s7c range. Anchors (s7c/map1c walls): h100 n500 16 min, h150 n500 19, h175 n500 20, h100 n1000 26 → **≈ 1.4 h for the four cells**, allowing up to ~2 h; inside the 3 h ceiling, 4 h hard timeout. Order: HR 1.00 n500 → 1.50 n500 → 1.75 n500 → 1.00 n1000; campaign `p30`, stems `…_z1q60_p30`; two seed-disjoint batches then combine; per-cell fail-fast (`stage2P_driver.sh`, session scratchpad; every render is the committed template driven by `FS_S7_*` env only). The save guard is live on every batch and pooled save.
