# REPORT — Field block for the complement Ĥᶜ: Stage 0 (Discovery)

**Task:** `dev/tasks/TASK_mr_field_complement_2026-09-06.md` (2c519337). H-C1–H-C4 at defaults; H-C5 unattended pre-authorization (6 h projection ceiling at 100 workers, 8 h hard timeout).
**Date:** 2026-09-06. Executor: Claude Code. No compute, no `R/` edits in this stage. Source read at tip dadec301 (`R/fs_mr_inference.R` last changed 72bd714a/a6702fd8; the installed package's `fs_mr_inference` deparses identically to the source, so the installed 0.3.5 is a valid pre-change reference for the Stage 1 byte-identities).

---

## GATE 0: PASS — the complement field can reuse the harm field's ξ draws and winners without changing them (argument under 0a.4).

## 0a — The gate's complement path and the field block (quoted from `R/fs_mr_inference.R`)

**1. Complement fits for draw winners, the complement influence matrix (`:584–609`).** Runs only under `include_complement = TRUE`, *before* the field block:

```r
winset <- sort(unique(c(winner[!is.na(winner)], sel)))
Bc   <- matrix(0, Nall, Ncol)
bh_c <- rep(NA_real_, Ncol); sdv_c <- rep(NA_real_, Ncol)
for (w in winset) {
  comp_idx <- setdiff(seq_len(Nall), kept[[w]])
  if (length(comp_idx) < 6L) next
  pcc <- tryCatch(.fs_mr_pieces(df[comp_idx, , drop = FALSE], spec), error = function(e) NULL)
  if (is.null(pcc) || length(pcc$dfbeta) != length(comp_idx)) next
  Bc[comp_idx, w] <- pcc$dfbeta
  bh_c[w] <- pcc$beta_hat; sdv_c[w] <- pcc$sigma_D
}
Pc <- crossprod(Bc, Xi)                      # Ncol x draws : D_{complement}(b)
```

So the complement influence matrix `Bc` is already allocated at full family width (`Nall × Ncol`, columns aligned with `asm$names`), populated lazily for the multiplier-draw winners plus the selected candidate; an unfit column is all-zero with `bh_c[w] = NA`. That is the cache the task's step 2 asks to extend: the same loop body, applied to any further candidate index, is the "lazy fit keyed by candidate". `kept <- candidates[asm$keep]` (`:587`) supplies the membership for any column.

**2. `complement$debiased` and its bias terms (`:626–652`).** `use_c <- which(is.finite(selb_c))`; `selbias_c = mean(selb_c[use_c])` (= `bias_sel_c`), `fixed_c = mean(Pc[sel, use_c])` (= `bias_fix_c`); `bnc <- bh_c[sel]`; `bdc <- bnc − sbc − fcc` is **β̃ᶜ**, the two-term de-biased complement estimate the task fixes as the centre of the complement field's inversion; `sec_used` is the IJ SE under `"ij"`/`"field"` (fix a6702fd8). The element is `complement = list(naive, debiased = list(est, lower, upper, lower_1s, se, se_ij, se_wald, var_ij, ij_source, ij_draws), selection_bias, fixed_bias, n)`; `complement <- list(note = "complement subgroup could not be fit")` when `bnc` is not finite.

**3. The field block's outer/inner structure and seeds (`:677–731`).**

```r
if (!is.null(seed)) set.seed(as.integer(seed) + 900000L)
Np <- nrow(B)
Zo <- crossprod(B, matrix(stats::rnorm(Np * field_R_out), Np, field_R_out))
Zi <- crossprod(B, matrix(stats::rnorm(Np * field_R_in), Np, field_R_in))
w <- bh; w[sel] <- beta_deb
fast <- identical(reselection, "maxeff") && is.null(t_g)
...
for (r in seq_len(field_R_out)) {
  v <- w + Zo[, r]
  G <- if (fast) which.max(v) else sel_one(v)
  if (is.na(G)) next
  if (fast) {
    win <- max.col(t(v + Zi), ties.method = "first")
    lam[r] <- Zo[G, r] - mean(Zi[cbind(win, ii)]);  n_in_used[r] <- field_R_in
  } else {
    wi <- vapply(ii, function(j) sel_one(v + Zi[, j]), integer(1))
    ok_in <- which(!is.na(wi)); if (!length(ok_in)) next
    lam[r] <- Zo[G, r] - mean(Zi[cbind(wi[ok_in], ok_in)]);  n_in_used[r] <- length(ok_in)
  }
}
```

Seed: `seed + 900000L`, consumed by exactly two `rnorm` calls (`Np·R_out` then `Np·R_in` variates) — the ξ matrices are drawn as plain `N(0, I)` matrices and projected at once through `crossprod(B, ·)`; the raw ξ are *not retained*. The inner batch `Zi` is shared across all outer draws (one `R_in`-column matrix). The winners `G` (outer) and `win`/`wi` (inner, per outer draw) are consumed inside the loop and not retained. The uniform sweep (`field_uniform`) draws under its own `+ 910000L` *after* this stream (`:737–751`).

**4. Why Gate 0 passes.** The complement field needs (i) the same ξ_r and ξ'_j and (ii) the same winners G_r and G(v_r + ζ'_j). Both are available add-only:

- (i) Keep the two `rnorm` matrices in local variables *before* projecting them (`Xo <- matrix(rnorm(...)); Zo <- crossprod(B, Xo)`) — the RNG stream and `Zo`/`Zi` are bit-identical to today's, because the same two `rnorm` calls are made in the same order with the same sizes. Under `field_complement = TRUE` the complement perturbations are then `Zo_c <- crossprod(Bc, Xo)`, `Zi_c <- crossprod(Bc, Xi_f)` — task step 3 (same multipliers, so the complement noise carries its correct cross-candidate correlation `B_effᵀ B_effᶜ`). Memory: `Np × (R_out + R_in)` doubles = 6 MB at n = 500, 18 MB at n = 1500.
- (ii) Record `G_r` (length `R_out`) and the inner winners (an `R_out × R_in` integer matrix, `NA` where the inner draw had no winner) *as the harm loop computes them*, under the flag. The harm loop's arithmetic, its `lam`, `n_in_used`, and every `field$*` output are untouched — the recording is a side assignment. This is the design choice over the alternative of a second pass re-deriving the winners from `Zo`/`Zi` (which would double the ~15 s selection cost for no gain).
- The complement block then runs *after* the harm field's assembly and *before* the uniform sweep's seed reset (the uniform sweep, when requested, still sees an identical RNG state because the complement block draws nothing). Task step 4: `Λ*ᶜ_r = Zo_c[G_r, r] − mean_j Zi_c[win_{r,j}, j]` over the inner draws with a winner *whose complement is fit*; outer draws dropped when the harm block dropped them, when `G_r`'s complement is unfit, or when no inner draw survives — counted.
- Candidates needing a complement fit beyond `winset`: the distinct values of `G_r` and of the inner-winner matrix not already in `winset`; fit through the same loop body (0a.1), recording `n_complement_fits` (distinct fits, `winset` included) and `share_draws_new_fit` (share of outer + inner draw-winner readings whose candidate was outside the multiplier-draw `winset`).

**5. Return assembly (`:775–802`).** `out$field <- field` is attached last; the complement field goes under `field$complement` (task step 5), so the `"ij"`/`"wald"` returns and `field`'s existing top-level names are untouched. `field$uniform` (when present) stays where it is.

**6. Forwarding.** `forestsearch()` forwards `field_uniform`, `return_reselection`, `field_M_cap` as add-only pass-throughs (`R/forestsearch_main.R:3392–3401`); `field_complement` takes the same one-line form there, with a roxygen `\item` in the `mr_inference_args` list (`:797–812`). The DINA/GRF caller `.fs_mr_run_generic()` (`R/fs_mr_inference_methods.R:127–142`) forwards neither `field_uniform` nor `field_R_*` and is left alone, as the uniform task left it.

## 0b — The template's complement recorder and the display

**Recorder (`sim_fs_maxeffCons_fb_mr_field_m1_template.qmd:683–714, 830–872`).** Complement columns: `betaHhat_Hc` (`:683`), `or_Hc_est/lo/hi/se` (`:686`, from the oracle fit at `:938`), `nv_Hc_est/lo/hi/se` and `mr_Hc_est/lo/hi/se_ij` (`:688, :692`, filled at `:830–839` from `g$complement$naive` / `g$complement$debiased`, `nv_Hc_se = debiased$se_wald`), `fb_Hc_*` (`:690`, FB join), `cde_Hc`/`marg_Hc` in `truth` (`:562–564`). The harm field block is recorded at `:844–872` from `g$field` with `fld_H_note` for the degenerate case; the uniform sub-block at `:859–870`. The new `fld_Hc_*` columns follow the same shape, read from `g$field$complement`.

**Coverage helpers.** `.cover(target, lo, hi)` (`:1518`), `.cov_n` / `.cov1_n` (`:1950–1960`; one-sided = `target >= lower bound`), `.wilson` (`:1944`); `est_keys_for()` (`:1549`) drops `MR (field)` unless `suffix == "H"`, `.est_cols()` (`:1557`) and `.get_se()` (`:1537`) read `fld_H_*` only; the one-sided bound for the field row is `fld_H_lo1s` (`:1976`). The interval invariant `.ci_check` (`:1025–1051`) lists field pairs explicitly. The diagnostics section (`:1930–2072`) is harm-only (`betaHhat_H`, `fld_H_*`), and the display chunk `field-bias-coverage` (`:2006–2022`) calls `fs_sim_bias_coverage(results, block = "H", estimators = c("mr", "fld"))`.

**`fs_sim_bias_coverage()` `block = "Hc"` branch (`R/fs_bias_coverage.R:56–60`):** drops `"fld"` with a message ("the field block exists for the harm block only"); `cols("fld")` (`:76–79`) is hard-wired to `fld_H_*`; the one-sided coverage is `tgt >= lo1` (`:104`) and the Gaussian reference is `pnorm(z1·r − b)` (`:117`). The `side = c("lower", "upper")` extension: under `"upper"`, `cov1 = mean(tgt <= up1)` with `up1` = `fld_<blk>_up1s` for the field row and `exp(log e + z1·se)` for naive/IJ, and `cov1_ref = pnorm(z1·r + b)` (the mirror: positive bias *helps* an upper bound).

## 0c — Cost anchors

| Cell (committed bundle) | Detected / 2000 | fit+MR s/rep (mean / q50 / q90, 100 workers, loaded) | field s/rep (mean / q90) | complement IJ SE/SD |
|---|---|---|---|---|
| h100 n500 (s7) | 1361 (68.0%) | 39.0 / 39.6 / 48.7 | 14.7 / 18.7 | 0.255 / 0.137 = 1.86 |
| h175 n500 (s7) | 1900 (95.0%) | 43.6 / 44.1 / 52.6 | 16.3 / 20.1 | 0.257 / 0.143 = 1.80 |
| h150 n500 (map1) | 1822 (91.1%) | 42.4 / 43.0 / 51.5 | 15.9 / 19.7 | 1.82 |
| h150 n1500 (map1) | 1976 (98.8%) | 147.3 / 148.6 / 172.8 | 36.0 / 46.6 | 1.84 |
| h075 n500 (map1) | 1042 (52.1%) | 36.9 / 37.5 / 46.1 | 13.9 / 18.0 | 1.85 |
| h100 n1000 (map1) | 1319 (66.0%) | 74.0 / 76.5 / 94.7 | 20.9 / 28.8 | 1.85 |
| h175 knoise3 (map1) | 1945 (97.2%) | 85.8 / 86.4 / 105.5 | 24.9 / 31.4 | 1.83 |

(`fit_mr_secs` includes the field; the s7 Stage 1 anchor under load was 28–40 s with the field at 13–15 s. `mr_Hc_est` is finite on every detected replicate in all seven cells.)

**Complement fit per candidate:** one `.consistency_cox_pieces()` call on a treatment-only Cox model is **2.4 ms at n = 400, 3.5 ms at n = 1,200** (200-call microbenchmark, this host, light load). Even several hundred new complement fits per replicate cost well under 2 s. The complement block's other costs are two `crossprod(Bc, ·)` products (`Ncol × Np × R`, ≈ 10⁸ flops at n = 500) and the `R_out × R_in` gather — sub-second. **Expected incremental cost per replicate: ≈ 1–3 s** (≤ 10% of the field block), so the Stage 2 projection is the s7/map1 wall (35–53 min per cell at 100 workers) plus a small increment.

## Stage 1 plan (for the record)

- `fs_mr_inference(field_complement = FALSE)`: guard `isTRUE(field_complement) && ci_method == "field" && include_complement`; retain `Xo`/`Xi_f` and the winners under the flag; complement block after the harm field's `field <- list(...)`, before the uniform sweep; `field$complement` per task step 5, or `list(note = ...)` when the selected complement is unfit or fewer than 2 outer draws survive.
- Roxygen: orientation (benefit claim, the exposed limit is the **upper** bound; under-correction vulnerable, over-correction immune — the mirror of the harm side); `@param field_complement`; `@return` paragraph. `man/` regenerated; `forestsearch()` forwarding line + `\item`; two tests (default-path identity with the flag off; K = 1 identities and bound identities with the flag on).
- `fs_sim_bias_coverage(side = c("lower", "upper"))`, `block = "Hc"` no longer drops `"fld"` when `fld_Hc_est2` is present.
- Template: `FS_S7_FIELD_COMPLEMENT` knob (default FALSE), `fld_Hc_*` recorder columns, meta `field_complement`, `.ci_check` pairs, `est_keys_for`/`.est_cols`/`.get_se` complement-aware, one-sided **upper** coverage for Ĥᶜ in the Wilson table, complement λ-SD/SD and retained-bias rows, `fs_sim_bias_coverage(block = "Hc", side = "upper")` chunk.
- Identities per task 1b, on the installed 0.3.5 as the pre-change reference (Guo–He three fixed-seed cases via `quarto/GuoHe/mr_vs_guohe_sim.R`'s `mv_mr()`, which gains a `field_complement` pass-through as it gained `field_uniform`), then smoke at campaign `s7csmoke` (5 replicates per s7 cell, `FS_S7_FIELD_COMPLEMENT=TRUE`) against the committed s7 bundles.
