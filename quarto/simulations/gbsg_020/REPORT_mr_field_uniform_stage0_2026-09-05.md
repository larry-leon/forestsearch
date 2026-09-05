# REPORT — Field interval, uniform calibration: Stage 0 (Discovery)

**Task:** `dev/tasks/TASK_mr_field_uniform_2026-09-05.md` (committed `1313f597` with the four reference files: `limit_sweep.py`, `limit_sweep_K.py`, `limit_sweep_K2_rules_2026-09-05.csv`, `limit_sweep_K10_kappa_2026-09-05.csv`). `TASK_mr_field_adaptive_2026-09-05.md` is withdrawn and was not started (no copy exists in the tree or `~/Downloads`).
**Date:** 2026-09-05. Executor: Claude Code. Stage 0 only; stopped at Gate 0 per instruction. Decisions H1–H5 at defaults; H6 at Gate 1; H7 after Stage 3.
**Source state:** forestsearch 0.3.5 + `a6702fd8` (installed); s7 campaign bundles at `36881ba7`.

---

## GATE 0: PASS — κ is computable without altering the field's own draws or the re-selection routine

- **Draws untouched:** the field block seeds its stream once, `R/fs_mr_inference.R:653` — `if (!is.null(seed)) set.seed(as.integer(seed) + 900000L)` — and consumes it entirely inside lines 654–656 (`Zo`/`Zi` built up front) before the selection loop. A κ block placed **after** the field block, seeded under its own derived offset (e.g. `seed + 910000L`, mirroring the `900000L` pattern that already protects the `"ij"`/`"wald"` streams), leaves every existing draw byte-identical — the same architecture the field feature itself used to leave `"ij"` untouched.
- **Re-selection routine untouched:** the κ sweep needs only to *call* the selection map, and every piece is in scope as a local function — `.admit` (`:463–472`), `.zcons` (`:486`), `.fs_mr_select` (definition `:132`), and the field's own `sel_one` closure (`:659–664`, including the `fast` argmax path `:658`, `:670–675`). Nothing needs modification.

Two add-only pass-through requirements surfaced (0b), both anticipated by the task's method-change paragraph; neither is a stop.

## 0a — In-scope inputs after the field block (`R/fs_mr_inference.R`, current tree)

Every input the κ computation (task §method, steps 1–6) names is a live local at the point where the field block ends (`:706`) and before the return assembly (`:713`):

| κ input | Object | Where bound |
|---|---|---|
| B_eff (n × K; Σ̂ = B_effᵀB_eff) | `B` | `:440` `B <- asm$B` (assembled `:102`); the field already draws ζ = Bᵀξ exactly this way, `:655–656` `Zo <- crossprod(B, matrix(stats::rnorm(Np * field_R_out), Np, field_R_out))` |
| β̂ (candidate effects) | `bh` | `:440` (`asm$beta_hat`, filled `:99`) |
| σ̂ per candidate (σ̂_Ĥ for the winner-profile) | `sdv`; `sdv[sel]` (= `se_wald`, `:530`) | `:440` |
| Candidate sizes | `sz` | `:440` |
| Ĥ (winner index) | `sel` | `:427` `sel <- match(sel_lab, asm$names)` |
| β̃(Ĥ) | `beta_deb` | `:529` `beta_deb <- beta_naive - selection_bias - fb` |
| Shrunk field w | `w` | `:657` `w <- bh; w[sel] <- beta_deb` (field-block local) |
| p̂ (re-selection frequencies) | from `winner` | `winner` filled per draw `:493/:500`; p̂ = `tabulate(winner[!is.na(winner)], nbins = length(asm$names)) / draws` — the exact expression already used at `:735` (inside `return_reselection`, `:734`), computable unconditionally |
| S (thresholds, sizes, rule) | `.admit`/`.zcons`/`.fs_mr_select`/`sel_one` + `t_g`, `reselection`, `effect_neighborhood`, `selection_rule`, `log_scale` | `:463–486`, `:132`, `:659–664` |

**Return-list assembly (the D6 pattern):** `out <- list(...)` at `:713–732`; conditional attachments `if (isTRUE(return_reselection)) out$reselection <- ...` `:734–738` and `if (!is.null(field)) out$field <- field` (final line before `out`). The uniform block attaches as `field$uniform` (a sub-element of the existing conditional attachment), so the default return names are unchanged — the Guo–He V8a return-names contract (`quarto/GuoHe/mr_vs_guohe_sim.R:360–371`) stays satisfied under defaults.

## 0b — Pass-through and recorder sites

**`fs_mr_inference()` signature:** `field_R_out = 1000L, field_R_in = 500L` sit at `:406–407` after `return_reselection` (`:405`); `field_uniform = FALSE` joins them add-only.

**The wrapper (flagged, add-only):** `forestsearch()` builds the gate call with a fixed argument list — `R/forestsearch_main.R:3385` `ci_method = .g_mr(mr_inference_args$ci_method, "ij"),` — and, per the Stage 0 finding of the s7 arc, does **not** forward `field_R_out`/`field_R_in`; `field_uniform` has the same gap. The template runs through `forestsearch()`, so the method change needs one add-only forwarding line in that call (e.g. `field_uniform = .g_mr(mr_inference_args$field_uniform, FALSE)`). Same class as the documented `field_R_*` gap; called out here for the Stage 1 record.

**Template (`quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`):**
- `mr_ci_method <- "field"` `:485`; `mr_inference_args <- c(list(ci_method = mr_ci_method, draws = mr_draws, include_complement = TRUE, confirm_rule = mr_confirm_rule), ...)` `:491–494` — the `field_uniform` knob joins this list.
- Recorder: `.na_record` field columns `:690–701` (`fld_H_est2` … `fld_H_note`); extraction `if (is.list(g$field)) { f <- g$field; ... }` `:831–849`. The eight new columns (task 1a: `fld_H_kappa`, `fld_H_M`, `fld_H_mass`, `fld_H_minC1`, `fld_H_lo2u`, `fld_H_hi2u`, `fld_H_kappa_mcse`, `fld_H_uniform_secs`) extend both sites; batch/combine machinery carries new columns without change (columns travel in `results`).

**Guo–He harness (`quarto/GuoHe/mr_vs_guohe_sim.R`):** `mv_mr()` `:124–138` calls `forestsearch:::fs_mr_inference` directly with `return_reselection = TRUE` (`:136` — p̂ already exercised there) and appends `field_R_out`/`field_R_in` when `ci_method == "field"` (`:137–138`); `field_uniform` passes by extending that same `args <- c(args, list(...))` line. No wrapper in the path.

## 0c — s7 seed scheme regenerates replicates bit-identically

Two independent proofs, both from committed artifacts:
1. **Same-vintage regeneration:** the committed `smoke`-campaign bundles (sims 1–5, separate renders) vs the `s7` combined bundles at the same sim_ids: **64 of 64 shared columns `identical()`** on both cells (excluded only: wall-clock columns and the FB-join columns, which differ by `fb_mode`/`join_skip` configuration, not by computation). Every per-replicate seed is `seed_base + sim_id` (template `:466/:482/:502` lineage; `seed_base = 8316951`), so κ can be evaluated on the same replicates by regeneration.
2. **Cross-vintage bound (Gate 2 record, 2026-09-05):** fresh full runs matched the 0.2.x-committed bundles exactly on all gated columns except the four enumerated vintage flips — under the *current fixed* package (the only one Stage 1+ will use) regeneration is exact, period.

## 0d — Cost anchors (s7 campaign, 100 workers, loaded host)

| Cell | fit_mr_secs (search + MR-IJ + field) mean / q50 / q90 / q99 | field block alone, mean / q90 | detections |
|---|---|---|---|
| h100 | 28.3 / 34.2 / 47.0 / 54.8 s | 14.7 / 18.7 s | 1361 / 2000 |
| h175 | 41.7 / 43.8 / 52.4 / 60.4 s | 16.3 / 20.1 s | 1900 / 2000 |

(The task's "~39 s" anchor is the mean over both cells.) **Cost-risk note for Stage 1d:** the current field block costs ~15 s at `R_out × R_in = 1000 × 500` selection evaluations on the full enumerated family. The κ sweep at H4 defaults is, per δ, R_rep × (1 + R_in + R_out·R_in-equivalent) selection work on the **reduced** family — the top-M cap (H3, M ≤ 12) and the vectorized `fast` argmax path are what the task's "seconds per δ / minutes per analysis" estimate leans on; whether the maxeffCons selection map on the reduced family can take a vectorized path (or M-sized `sel_one` calls are cheap enough) is exactly what Stage 1d must measure before H6. At face value, "minutes per analysis" × ~1,300–1,900 detections / 100 workers ≈ 0.5–1.3 h per cell — feasible but to be measured, not assumed.

## Deviations

None. No `R/` edits, nothing re-run, committed bundles read-only throughout. **Stopped at Gate 0** per instruction; Stage 1 (implementation + identities against the two reference CSVs + projection) awaits Larry's go and the H6 compute decision at Gate 1.
