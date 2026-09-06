# REPORT — Complement inference refinements: Stage 0 (Discovery)

**Task:** `dev/tasks/TASK_complement_refinements_2026-09-06.md` (95e6c8c5). K-1–K-3 at defaults; K-4 unattended pre-authorization (5 h ceiling at 100 workers, 7 h hard timeout).
**Date:** 2026-09-06. Executor: Claude Code (same session as the complement task per Larry's queue; the document's "fresh session" line is superseded). No compute, no `R/` edits. Source read at tip 815debf1 (`R/fs_mr_inference.R` as committed at ef3e609a; the installed 0.3.5 deparses `fs_mr_inference` identically, so it is the pre-change reference for Stage 1).

---

## GATE 0: PASS — the joint draws `(Λ*_r, Λ*ᶜ_r)` are already aligned by outer-draw index; nothing in either block needs to change (0b).

## 0a — The IJ residual, the finite-B correction, the fallbacks, and where the SEs feed the bounds (`R/fs_mr_inference.R`)

**Harm block (`:601–610`):**

```r
r_H   <- (selection_bias + fb) - sel_bias - P[sel, ]
mean_r <- mean(r_H[ok_H])
ijH   <- .fs_mr_ij_var(Xi, r_H, ok_H)
se_ij <- .fs_mr_se_from_ij(ijH, se_wald)
se    <- if (ci_method == "wald") se_wald else se_ij$se
```

`sel_bias[b] = P[winner_b, b] = D_{Ĥ*_b}(b)` (the re-selected winner's perturbation, `:537`), `P[sel, ]` = `D_Ĥ(b)` (the same-draws term), both over `ok_H` (draws with a winner). The two-term residual is exactly the task's `r_b`; the winner-only residual is `r_b^w = selection_bias − sel_bias` (drop `fb` and `P[sel, ]`).

**Complement block (`:675–690`):**

```r
r_c   <- (sbc + fcc) - selb_c - Pc[sel, ]
mean_r_c <- mean(r_c[use_c])
ijC   <- .fs_mr_ij_var(Xi, r_c, use_c)
se_ijc <- .fs_mr_se_from_ij(ijC, sec)
sec_used <- if (ci_method %in% c("ij", "field")) se_ijc$se else sec
complement <- list(naive = ..., debiased = list(est = to_eff(bdc),
    lower = to_eff(bdc - z975 * sec_used), upper = to_eff(bdc + z975 * sec_used),
    lower_1s = to_eff(bdc - stats::qnorm(0.95) * sec_used),
    se = sec_used, se_ij = se_ijc$se, se_wald = sec, var_ij = se_ijc$var,
    ij_source = se_ijc$source, ij_draws = ijC$B_ok), ...)
```

`selb_c[b] = Pc[winner_b, b]` (`:645`), `Pc[sel, ]` the complement's same-draws term, over `use_c` (draws with a winner whose complement is fit). Winner-only complement residual: `sbc − selb_c`. `sec = sdv_c[sel]` is the complement's naive (robust) SE — the floor for `winner_floor`.

**Finite-B correction and fallbacks (`:195–231`):** `.fs_mr_ij_var()` centres the multipliers over the used draws, `cov_i = (1/B_ok) Σ_b (K*_bi − K̄_i) r_b`, `tilde_V = Σ_i cov_i²`, `hat_V = tilde_V − (N/B_ok)·mean(r_b²)` (Wager 2014). `.fs_mr_se_from_ij()` takes `hat_V` when positive (`ij_source = "ij"`), else `tilde_V` (`"ij_raw"`), else the robust SE (`"wald_fallback"`). The harm `debiased` element (`:853–859`) reads `se`, `se_ij`, `var_ij`, `ij_source`, `ij_draws` from `se_ij`/`ijH`; the harm one-sided bound `ci_lo_1s = beta_deb − z₀.₉₅·se` (`:613`) also feeds the `"ci"` confirm rule.

**Design for A.** Compute, beside the existing `ijH`/`ijC`, `ijH_w = .fs_mr_ij_var(Xi, selection_bias − sel_bias, ok_H)` and `ijC_w = .fs_mr_ij_var(Xi, sbc − selb_c, use_c)` — same `Xi`, same draw sets, one term dropped — resolve each through `.fs_mr_se_from_ij()` (same fallback ladder, floor = naive SE for `winner_floor` as `sqrt(max(V, σ̂²))`), and attach `se_ij_winner`, `se_ij_winner_floor`, `lower_w/upper_w/lower_1s_w`, `lower_wf/upper_wf/lower_1s_wf` (complement: `upper_1s_w/_wf` too, its exposed side) as **additional** elements of `debiased` / `complement$debiased`. `ij_residual` then selects which variance populates the *existing* `se`/`se_ij`/`var_ij`/`ij_source`/`lower`/`upper`/`lower_1s` (and `sec_used` → the complement's) — under the default `"two_term"` those are computed exactly as today, so the default is byte-identical with the new elements excluded. The harm flag under `confirm_rule = "ci"` follows the reported `se`, as documented. Point estimates untouched. The uniform sweep and the field block read `beta_deb`/`sdv[sel]` only (not `se`), so they are unaffected by `ij_residual`.

## 0b — Alignment of the harm field's `lam` and the complement field's `Λ*ᶜ`

**Harm field (`:745–768`):** `lam <- rep(NA_real_, field_R_out)`; per outer draw `r`, `lam[r] <- Zo[G, r] − mean(Zi[cbind(win, ii)])` (fast path `:756`) or `... − mean(Zi[cbind(wi[ok_in], ok_in)])` (`:763`); a draw with no outer winner (`next` at `:751`) or no inner winner (`:761`) leaves `lam[r]` NA. `ok_f <- which(is.finite(lam))` (`:768`), `lf <- lam[ok_f]`; the quantiles and `est2` are taken over `ok_f`. Under `field_complement`, the same iteration stores `G_out[r]` and `W_in[r, ]` (`:758`, `:765`) — written only on iterations that also wrote `lam[r]`.

**Complement field (`.fs_mr_field_complement()`, `:927–948`):** `lam_c <- rep(NA_real_, R_out)` with `R_out = length(G_out)`; `for (r in which(!is.na(G_out)))`: `G <- G_out[r]`; skip (counted in `n_drop_unfit`) if `G`'s complement is unfit or no inner winner's complement is fit; else `lam_c[r] <- Zo_c[G, r] − mean(Zi_c[cbind(wi[ok_in], ok_in)])` (`:935`) with `wi <- W_in[r, ]`. `ok_c <- which(is.finite(lam_c))` (`:938`).

So `lam` and `lam_c` are two length-`R_out` vectors indexed by the **same** outer draw `r`, built from the same `Xo[, r]` column and the same winner `G_out[r] = G`; `ok_c ⊆ ok_f` by construction (the complement writes only where the harm field wrote), and in every Stage 2 cell `n_out_dropped_unfit = 0`, so `ok_c = ok_f`. The joint draws for B are `(lam[ok_c], lam_c[ok_c])` — aligned without re-drawing and without touching either block. The complement helper does not currently return `lam_c`; B needs the per-draw vectors. **Chosen design:** the helper gains arguments `lam_H = lam`, `beta_deb` and `alpha`, computes the joint quantities inside (where both vectors are in scope) and returns `list(complement = <unchanged list>, joint = <B>)`; the caller attaches `field$complement` and `field$joint`. `field$complement`'s contents are then unchanged element-for-element (the s7c/map1c anchor), and no per-draw vector is added to the return.

**Joint level γ (B):** over `ok_c`, with `Λ = lam[ok_c]`, `Λᶜ = lam_c[ok_c]`, grid `γ ∈ {α, α − 0.001, …, α/2}`, find the largest γ with `mean(Λ ≤ q_{1−γ}(Λ) & Λᶜ ≥ q_γ(Λᶜ)) ≥ 1 − α` (type-7 quantiles as the field uses); return `gamma`, `joint_prob` (achieved), `lower_H = to_eff(beta_deb − q_{1−γ}(Λ))`, `upper_Hc = to_eff(bdc − q_γ(Λᶜ))`, the Bonferroni pair at γ = α/2, `corr = cor(Λ, Λᶜ)`, `n_joint_draws = length(ok_c)`; the marginal one-sided bounds (`field$lower_1s` = `beta_deb − q₀.₉₅(Λ)` over `ok_f`, `field$complement$upper_1s` = `bdc − q₀.₀₅(Λᶜ)` over `ok_c`) are untouched. Note the orientation: the harm lower bound inverts the *upper* tail of `Λ*` and the complement upper bound the *lower* tail of `Λ*ᶜ`, so the joint event is `Λ ≤ q_{1−γ}(Λ)` and `Λᶜ ≥ q_γ(Λᶜ)` — the task's statement. Under independence `(1−γ)² ≥ 1−α` gives γ = 1 − √(0.95) = 0.0253.

## 0c — Cost anchors (`REPORT_mr_field_complement_gate2_2026-09-06.md`)

Cell walls at 100 workers with the complement field on: h100 n500 16 min, h175 n500 19, h150 n500 19, h150 n1500 57, h075 14, h100 n1000 26, h175 knoise3 35 — **3 h 09 m for all seven**. Per replicate: fit+MR 39–46 s (n = 500), 80 s (n = 1000), 91 s (knoise3), 162 s (n = 1500). A and B add two `.fs_mr_ij_var()` calls per block (each an `N × B_ok` matrix-vector product, ≈ 10⁶–10⁷ flops — milliseconds) and a 26-point quantile grid on two length-≤1000 vectors: no measurable cost. **Projection for Stage 2 ≈ 3.2 h**, inside the 5 h ceiling → all seven cells.

## Note on the cited PoC file

`dev/tasks/POC_mr_interval_alternatives_2026-09-05.md` is present (the "bag"/"bagfloor" constructions: bag SE/SD 0.51–1.00 and one-sided coverage 0.63–0.95 across the PoC scenarios — exact when one candidate dominates, under-covering at ties; bagfloor 0.97–1.31 / 0.925–0.951). The companion `poc_ci_results_2026-09-05.csv` the task names is **not** in `dev/tasks/` (only `poc2_results_2026-09-05.csv`, the tie-sign PoC); the markdown's quoted figures are the reference used.

## Stage 1 plan

- A as designed above; `ij_residual` forwarded from `forestsearch()` (add-only line + `\item`); template knob `FS_S7_IJ_RESIDUAL` (default `two_term`, stored in meta), recorder columns `mr_H_se_w`, `mr_H_lo_w/hi_w`, `mr_H_se_wf`, `mr_H_lo_wf/hi_wf` and the `Hc` twins; B's `fld_joint_gamma`, `fld_joint_loH`, `fld_joint_upHc`, `fld_joint_bonf_loH`, `fld_joint_bonf_upHc`, `fld_joint_corr`, `fld_joint_n`, `fld_joint_prob`; tables: complement and harm blocks gain "MR (IJ, winner)" and "MR (IJ, winner-floor)" rows; a joint-bounds table (three pairs' joint coverage and margins); the complement display gains the two winner rows (as extra `fs_sim_bias_coverage()` estimators `mr_w`, `mr_wf`, add-only).
- Identities per task 1b against the 15 pre-change reference objects captured with the installed ef3e609a package (A/B/C × ij/field × complement, plus `field_complement = TRUE`), the K = 1 4σ² identity, an exchangeable K = 10 tie simulation, the joint-γ checks (range, achieved probability, the independence construction), then smoke at campaign `s7wsmoke` against the s7c bundles.
