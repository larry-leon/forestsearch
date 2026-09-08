# REPORT — Studentized complement field E1: Stage 0 (discovery), Stage 1 (the R1 "field-s" construction; identity gates G2a–G2c) and Gate 1 (projection)

**Task:** `dev/tasks/TASK_field_studentize_e1_2026-09-08.md` (9a6e21cc). Proposal §5 variant **R1** (decided on record by the E0 review). Predecessors: Stage 1 record 7245e898, E0 record 9d71f736. Decisions E-1–E-5 as stated.
**Date:** 2026-09-08. Executor: Claude Code (Linux), unattended. `R/` changes in exactly the two named files, classified **adds code; byte-identical defaults; new outputs only on the enabled path**. Winner-only and winner-floor excluded from every table and line.

---

## GATE 2a/2b/2c: PASS — G2a (both knobs off): 131 / 131 pre-existing non-timing columns `identical()` to the committed nb20-A HR 1.50 rows 1–5, `truth` identical, all 22 new columns present and NA. G2b (decompose + scale on): 131 / 131 identical; every `_s` and `joint_s` column finite on 5 / 5 detected rows; interval invariants hold on the `_s`, unscaled and harm blocks; `se_field_s > 0`; `joint_s` γ in [0.025, 0.05]; `se_field_s / (ρᶜ · se_field)` 0.9925–1.0012 (mean 0.997). G2c (e0stud config, decompose on, scale off): 135 / 135 committed columns identical to the committed e0stud rows 1–5, the four scale columns included — the ensure-fit / `s` relocation changed nothing. **GATE 1: compute go** — projection ≈ 2 h 52 m at 100 workers (≤ 3.5 h ceiling), all six cells, none deferred at the start.

## Stage 0 — Discovery (quotes from HEAD 9a6e21cc = the 7245e898 tree for these files; line numbers of that tree)

**The Stage-1 decomposition block with the `sel` ensure-fit (`R/fs_mr_inference.R:1064–1086`):**
```r
  decomp_fields <- NULL
  if (isTRUE(field_decompose)) {
    if (!fit_ok[sel]) {                      # ensure the selected complement's column is populated
      comp_idx <- setdiff(seq_len(Nall), kept[[sel]])
      if (length(comp_idx) >= 6L) {
        pcc <- tryCatch(.fs_mr_pieces(df[comp_idx, , drop = FALSE], spec), error = function(e) NULL)
        if (!is.null(pcc) && length(pcc$dfbeta) == length(comp_idx)) {
          Bc[comp_idx, sel] <- pcc$dfbeta; bh_c[sel] <- pcc$beta_hat; fit_ok[sel] <- TRUE
        }
      }
    }
    s  <- sqrt(colSums(Bc * Bc))             # per-candidate complement noise scale; 0 for unfit
    gG <- G_out[ok_c]; s_win <- s[gG]
    zg <- Zo_c[cbind(gG, ok_c)]; mi <- zg - lf
    decomp_fields <- list(
      scale_sel      = if (fit_ok[sel]) s[sel] else NA_real_,
      scale_win_mean = mean(s_win),
      scale_win_cv   = stats::sd(s_win) / mean(s_win),
      scale_ratio_c  = if (fit_ok[sel]) s[sel] / mean(s_win) else NA_real_,
      var_zeta_G = stats::var(zg), var_m_in = stats::var(mi), cov_zeta_m = stats::cov(zg, mi))
  }
```
**The `lam_c` loop (1026–1041):**
```r
  fit_ok <- is.finite(bh_c)
  Zo_c <- crossprod(Bc, Xo)                  # Ncol x R_out : zeta^c (outer)
  Zi_c <- crossprod(Bc, Xi_f)                # Ncol x R_in  : zeta'^c (inner)
  lam_c  <- rep(NA_real_, R_out)
  n_in_c <- rep(NA_real_, R_out)
  n_drop_unfit <- 0L
  for (r in which(!is.na(G_out))) {
    G <- G_out[r]
    if (!fit_ok[G]) { n_drop_unfit <- n_drop_unfit + 1L; next }
    wi <- W_in[r, ]
    ok_in <- which(!is.na(wi))
    ok_in <- ok_in[fit_ok[wi[ok_in]]]
    if (!length(ok_in)) { n_drop_unfit <- n_drop_unfit + 1L; next }
    lam_c[r]  <- Zo_c[G, r] - mean(Zi_c[cbind(wi[ok_in], ok_in)])
    n_in_c[r] <- length(ok_in)
  }
```
**Bound assembly (1052–1056, 1090–1106) and the joint call (1088–1089):**
```r
  lf <- lam_c[ok_c]
  qs <- stats::quantile(lf, c(.05, .25, .50, .75, .95, .025, .975),
                        names = FALSE, type = 7)
  est2_w <- bdc - mean(lf)
  sd_c <- stats::sd(lf)
  ...
  joint <- if (!is.null(lam_H) && is.finite(beta_deb))
    .fs_mr_field_joint(lam_H[ok_c], lf, beta_deb, bdc, to_eff, alpha) else NULL
  complement <- c(list(
    lambda_mean = mean(lf), lambda_sd = sd_c,
    q05 = qs[1], q25 = qs[2], q50 = qs[3], q75 = qs[4], q95 = qs[5],
    q025 = qs[6], q975 = qs[7],
    n_in_used_mean = mean(n_in_c[ok_c]),
    est2 = to_eff(est2_w),
    # Primary: the one-sided 95% UPPER bound (benefit claim "at most U").
    upper_1s = to_eff(bdc - qs[1]),
    lower_1s = to_eff(bdc - qs[5]),
    lower_2s = to_eff(bdc - qs[7]), upper_2s = to_eff(bdc - qs[6]),
    se_field = sd_c,
    lower_se = to_eff(est2_w - z975 * sd_c),
    upper_se = to_eff(est2_w + z975 * sd_c)),
    counts,
    list(timing_seconds = as.numeric((proc.time() - t0c)["elapsed"])))
  if (!is.null(decomp_fields)) complement$decomp_fields <- decomp_fields   # add-only; absent when FALSE
  list(complement = complement, joint = joint)
```
**`.fs_mr_field_joint()` signature (1133–1134):** `.fs_mr_field_joint <- function(lam_H, lam_c, beta_deb, bdc, to_eff, alpha = 0.05)` — grid of γ from α down to α/2 (step 0.001), the largest γ with joint probability ≥ 1 − α, else the α/2 fallback with its achieved probability.
**`forestsearch_main.R` forwarding (3416–3419):**
```r
        field_complement = .g_mr(mr_inference_args$field_complement, FALSE),
        # Add-only pass-through (TASK_field_studentize_stage1_e0_2026-09-08):
        # the complement field's scale diagnostics; FALSE is the gate's default.
        field_decompose = .g_mr(mr_inference_args$field_decompose, FALSE),
```
**Template:** knob `mr_field_decompose <- identical(.env_chr("FS_S7_FIELD_DECOMP", "FALSE"), "TRUE")` (`:595`), forwarded at `:608`; recorder `fld_Hc_se = NA_real_` (`:838`), the four scale columns `fld_Hc_scale_sel / _win / _cv / _ratio` (`:844–845`); the existing complement bound / estimate / SE columns (`:833–838`, `:846`): `fld_Hc_est2, fld_Hc_up1s, fld_Hc_lo1s, fld_Hc_lo2s, fld_Hc_hi2s, fld_Hc_lo_se, fld_Hc_hi_se, fld_Hc_se, fld_Hc_lam_mean`; the joint columns (`:867–871`): `fld_joint_gamma, fld_joint_prob, fld_joint_loH, fld_joint_upHc, fld_joint_bonf_loH, fld_joint_bonf_upHc, fld_joint_bonf_prob, fld_joint_corr, fld_joint_n`. The `_s` companions mirror these names verbatim.
**Comparators (all `_combined_1_2000.rds`, 2000 rows, 136 columns, seeds `8316951 + sim_id`, `field_complement = TRUE`, `ij_residual = two_term`, `fb_mode = none`, J = 10, z1q 0.60, n 500):** `fs_effMaxSG_…_nb20_p30sgnb20` (ε 0.20; HR 1.50 / 1.75), `fs_effMaxSG_…_nb30_banddial` (ε 0.30; HR 1.50 / 1.75; 1999 detected each), `fs_maxSG_…_z1q60_banddial` and `fs_minSG_…_z1q60_banddial` (HR 1.75; 1999 detected each).

## Stage 1 — The R1 construction (diff: 5 files, 156 insertions, 22 deletions; the deletions are the relocated ensure-fit block and re-formed lines)

| File | Change |
|---|---|
| `R/fs_mr_inference.R` | `fs_mr_inference()` gains `field_scale_complement = c("none", "selected")` (`match.arg`; roxygen), forwarded at the complement call site; `field$joint_s` attached beside `field$joint` when returned. `.fs_mr_field_complement()` gains the argument; `scale_on <- identical(field_scale_complement, "selected")`. The ensure-`sel`-fit is relocated unchanged to immediately after the original lazy-fit loop and `s <- sqrt(colSums(Bc * Bc))` computed there under `scale_on || field_decompose` (invariance comment: `sel ∈ winset`, so the multiplier stage already attempted its fit; the block can only populate `Bc[, sel]` when that fit failed, and then `bdc` is NA and the function has already returned). `lam_cs <- rep(NA_real_, R_out)`; after the `lam_c[r]` line, under `scale_on`: `lam_cs[r] <- (s[sel] / s[G]) * Zo_c[G, r] - mean((s[sel] / s[wi[ok_in]]) * Zi_c[cbind(wi[ok_in], ok_in)])`. After `sd_c`: `lfs <- lam_cs[ok_c]`, quantiles (same probs, type 7), `est2s_w`, `sd_cs`. The decompose block reuses `s`. `joint_s <- .fs_mr_field_joint(lam_H[ok_c], lfs, beta_deb, bdc, to_eff, alpha)` under `scale_on`. The `complement` list gains `lambda_mean_s, se_field_s, est2_s, upper_1s_s, lower_1s_s, lower_2s_s, upper_2s_s, lower_se_s, upper_se_s` (each the exact analogue of its twin, inverted around the same `bdc`); the return list carries `joint_s` (NULL when off). |
| `R/forestsearch_main.R` | `field_scale_complement = .g_mr(mr_inference_args$field_scale_complement, "none")` (one line + comment). |
| `man/fs_mr_inference.Rd` | regenerated. |
| `NEWS.md` | one development bullet. |
| Template | knob `FS_S7_FIELD_SCALEC` (default `none`; `stopifnot` in {none, selected}) → `mr_inference_args$field_scale_complement`; echoed in the knob audit line and the MR-settings readout; recorder gains the nine `fld_Hc_*_s` companions and the nine `fld_joint_s_*` mirrors, filled `%||% NA_real_`; batch and combined meta gain `field_decompose` / `field_scale_complement` (record only, not in the poolability keys). |

`devtools::document()`; `devtools::install(dependencies = FALSE)`: installed 0.3.5; `deparse()` of the installed `fs_mr_inference` and `.fs_mr_field_complement` identical to source; installed `forestsearch()` carries the pass-through.

## Identity gates (template-driven; 5 workers; seeds `8316951 + sim_id`, sim_id 1–5; `gate_e1.R`, session scratchpad; the five timing columns excluded as always)

**G2a — both knobs off** (nb20-A HR 1.50 config, tag `e1inert`, 67 s): **PASS.** 131 / 131 pre-existing non-timing columns `identical()` to `..._nb20_p30sgnb20_res_1_1000.rds` rows 1–5; `truth` identical; fresh bundle 158 columns = 136 + 4 scale + 9 `_s` + 9 `joint_s`, all NA.

**G2b — `FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_SCALEC=selected`** (tag `e1scale`, 69 s): **PASS.** 131 / 131 identical; `truth` identical; all 18 `_s` / `joint_s` columns finite on 5 / 5 detected rows; invariants `lo2s_s ≤ lo1s_s ≤ est2_s ≤ up1s_s ≤ hi2s_s` TRUE (and on the unscaled and harm blocks); `se_field_s > 0`; `joint_s` γ ∈ [0.025, 0.05] with achieved probability ≥ 0.95 off the α/2 fallback (3 of 5 rows at the fallback, as 4 of 5 for the unscaled joint on these rows).

| sim_id | \|Ĥ\| | p̂ | `nv_Hc_se` | `fld_Hc_se` | `fld_Hc_se_s` | ρᶜ | `fld_Hc_up1s` | `fld_Hc_up1s_s` | joint γ / joint_s γ | `se_field_s / (ρᶜ·se_field)` |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 161 | 0.013 | 0.1509 | 0.1340 | 0.1529 | 1.143 | 0.882 | 0.910 | 0.025 / 0.025 | 0.998 |
| 2 | 69 | 0.186 | 0.1261 | 0.1264 | 0.1220 | 0.973 | 1.243 | 1.233 | 0.025 / 0.025 | 0.993 |
| 3 | 113 | 0.064 | 0.1394 | 0.1335 | 0.1402 | 1.053 | 0.923 | 0.934 | 0.025 / 0.025 | 0.998 |
| 4 | 180 | 0.037 | 0.1485 | 0.1400 | 0.1542 | 1.105 | 1.077 | 1.098 | 0.025 / 0.026 | 0.996 |
| 5 | 87 | 0.159 | 0.1409 | 0.1347 | 0.1373 | 1.018 | 0.921 | 0.928 | 0.026 / 0.026 | 1.001 |

R1 vs the global rescale: 0.9925–1.0012, mean 0.997 (SD 0.003) — at a winner-scale CV of 0.03–0.06 the per-draw form and the global form agree to within 1%, as expected; recorded, not gated. `fld_Hc_se_s` tracks `nv_Hc_se` on every row (0.1529 vs 0.1509, 0.1220 vs 0.1261, …) where `fld_Hc_se` does not (0.134 flat). Light-load fit+MR 36.1 s/rep (nb20-A smoke reference 36.3; `e1inert` 35.7); complement field block 2.0 s vs 1.9 s — the studentized companions cost nothing measurable.

**G2c — regression on Stage 1** (e0stud config: HR 1.75, decompose on, scale off; tag `e1regr`, 69 s): **PASS.** 135 / 135 committed non-timing columns `identical()` to the committed `..._nb20_e0stud_res_1_200.rds` rows 1–5 — the four scale columns included (`fld_Hc_scale_ratio` 1.1568, 1.0039, 1.1682, 1.0642, 1.0126 on both sides); `truth` identical.

Gate bundles and renders committed beside the campaign's (`e1inert`, `e1scale`, `e1regr`, `_res_1_5.rds` / `_batch_1_5.html`).

## Gate 1 — Projection (100 workers): compute go

Light-load fit+MR per replicate under field-s: 36.1 s (ε 0.20 HR 1.50) against the nb20-A smoke's 36.3 s and the banddial smoke's 36.8 s (ε 0.30) — ×1.0. Scaling the reference walls (band cells 31–32 min per 2,000, size-rule cells 23 min): 4 × 32 + 2 × 23 ≈ **2 h 34 m of renders, ≈ 2 h 52 m with combines and gates** — inside the 3.5 h ceiling (5 h hard timeout). All six cells proceed, none deferred at the start; the driver (`e1_campaign.sh`, session scratchpad) re-projects each cell from the realized wall of the same kind and defers, listing it, any cell whose projected finish would cross the ceiling; every render is the committed template driven by `FS_S7_*` env only; the save guard is live. Order: ε 0.20 (HR 1.50 → 1.75), ε 0.30 (1.50 → 1.75), `maxSG`, `minSG`; two seed-disjoint batches then combine per cell; Gate 2 per cell after each combine; stop-on-failure per cell, the next proceeds.
