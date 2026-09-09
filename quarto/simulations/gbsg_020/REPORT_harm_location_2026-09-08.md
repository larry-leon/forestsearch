# REPORT — Harm-side location by selection stratum on `e1stud` (zero compute): does the stable-pick under-correction sit on both blocks?

Date: 2026-09-08. Task: `dev/tasks/TASK_harm_location_2026-09-08.md` (committed as received, 3a23fd5f). Executor: Claude Code (Linux). **Report-and-wait: no repair proposal, no recommendation.** No compute; no `R/` change; no bundle written (asserted at the end).

Predecessors: `REPORT_complement_location_2026-09-08.md` (A5: the complement's uncorrected fraction ≈ 0 / 0.13–0.22 / 0.41–0.44 across p̂ tertiles in the four band cells), `REPORT_field_studentize_e1_2026-09-08.md` (constructions tables, harm block: field lower coverage 0.936–0.980, IJ two-term 0.971–0.999 at cell level), `summary_complement_variance.qmd` (A5 machinery; extended here with A6).

## Stage 0 — quotes for the record (HEAD 3a23fd5f)

**Line-number note.** The task cites `R/fs_mr_inference.R` lines 776–782 for the shrunk-field comment. In HEAD those lines are inside the complement's `debiased = list(...)` (the winner-only `lower_w` / `upper_w` / `lower_1s_w` entries); the comment the task describes is at **lines 800–812**, and `w` is assembled at **line 827**. Quoted from HEAD:

```r
  # Shrunk field w = beta_hat with the winner's entry replaced by the two-term
  # de-biased estimate (E1).  Gaussian-multiplier perturbations zeta = B' xi,
  # xi ~ N(0, I_n) -- no Cholesky, no explicit Sigma.  Per outer draw r:
  #   v_r = w + zeta*_r;  G_r = S(v_r);
  #   m-hat(v_r) = mean_j zeta'_{j, S(v_r + zeta'_j)}  (shared inner draws,
  #                draws with no winner skipped, the bias_sel convention);
  #   Lambda*_r  = zeta*_{r, G_r} - m-hat(v_r).
  # The interval inverts Lambda* around beta_deb (basic-bootstrap form).
  # S is the gate's own re-selection: under maxeff with no admission floor it
  # is a plain argmax (vectorized over inner draws; ties.method = "first"
  # matches which.max); any other configuration goes through .fs_mr_select
  # per draw, identically to the main loop above.
```

Line 827: `    w <- bh; w[sel] <- beta_deb`

The harm field's estimate and one-sided lower bound (lines 868, 875–876): `est2_w <- beta_deb - mean(lf)` … `est2 = to_eff(est2_w)`, `lower_1s = to_eff(beta_deb - qs[5])` (the 95% quantile of Λ*, the basic-bootstrap form; **not** an SE form); `lambda_mean = mean(lf)`, `se_field = sd_f`. The harm block's `debiased$est` is `to_eff(beta_deb)` (line 957).

**`.fs_mr_select()`, `effMaxSG` branch** (line 173, inside the `switch(rule, ...)` at lines 170–178):

```r
    effMaxSG = { b <- .inband(); b[which.max(sizes[b])] },
```

with `.inband()` (lines 135–167) taking the passers' natural-scale effects `eff` and sizes `sz`, the band `ib <- .compute_inclusion_band(hr_vec = eff, n_vec = sz, selection_rule = selection_rule, effect_neighborhood = nbhd) == 1L`, the empty-band fallback `if (!any(ib)) ib <- rep(TRUE, length(passers))`, and returning `passers[ib]`.

**Template harm-block recorder columns** (`sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`, names verbatim). Declared (lines 805, 809, 813, 820–825, 906):

```r
  betaHhat_H = NA_real_, betaHhat_Hc = NA_real_,
  nv_H_est = NA, nv_H_lo = NA, nv_H_hi = NA, nv_H_se = NA,
  mr_H_est = NA, mr_H_lo = NA, mr_H_hi = NA, mr_H_se_ij = NA,
  fld_H_est2 = NA_real_,                                   # two-term shrunk-field estimate
  fld_H_lo2s = NA_real_, fld_H_hi2s = NA_real_,            # two-sided 95% quantile interval (primary, F4)
  fld_H_lo1s = NA_real_,                                   # one-sided 95% lower bound (supplementary)
  fld_H_lo_se = NA_real_, fld_H_hi_se = NA_real_,          # SE-type interval: est2 -+ z * lambda_sd
  fld_H_se = NA_real_,                                     # se_field = lambda_sd (log-HR scale)
  fld_H_lam_mean = NA_real_,                               # mean Lambda* (log-HR scale)
  p_hat_H = NA_real_, p_hat_sum = NA_real_,
```

Assigned (lines 1044–1046, 1078–1083, 1013–1014): `rec$nv_H_est <- g$naive$est; … rec$nv_H_se <- g$debiased$se_wald %||% NA_real_`; `rec$mr_H_est <- g$debiased$est; rec$mr_H_lo <- g$debiased$lower; rec$mr_H_hi <- g$debiased$upper`; `rec$mr_H_se_ij <- g$debiased$se_ij %||% NA_real_`; `rec$fld_H_est2 <- f$est2`; `rec$fld_H_lo1s <- f$lower_1s`; `rec$fld_H_se <- f$se_field`; `rec$fld_H_lam_mean <- f$lambda_mean`; `rec$p_hat_H <- if (!is.na(rec$label) && rec$label %in% names(ph)) unname(ph[[rec$label]]) else NA_real_` with `ph <- g$reselection$p_hat`.

**Which harm-side analogues of the complement quantities exist as columns.** All of them: `nv_H_est` / `nv_H_se` (naive), `mr_H_est` / `mr_H_se_ij` (two-term, IJ), `fld_H_est2` / `fld_H_se` / `fld_H_lo1s` (field), **`fld_H_lam_mean` is recorded** (so `lam_H` is read directly and the identity `ef_H = e_H − lam_H` is checkable), `betaHhat_H`, `p_hat_H`. Not present on the harm block: any `_s` companion (`field_scale_complement` acts on the complement only; the harm field has no studentized variant), and a stored one-sided IJ lower bound (built as `exp(log(mr_H_est) − z₀.₉₅·mr_H_se_ij)`, the `e1stud_findings.R` construction, line 81). The harm block's `selection_bias` / `fixed_bias` split is not a recorded column either; `c_H` is read whole, as `c` was in A5.

## A6 — Harm-side location by stratum

**Data.** The six committed `e1stud` pooled bundles (`results/*_e1stud_combined_1_2000.rds`, 2,000 rows each, 158 columns), read by `summary_complement_variance.qmd` with `FS_SUMCV_GLOBS="results/*_e1stud_combined_1_2000.rds"`. A6 takes A5's rows (`loc()`, detected with every complement input finite) and requires the harm-block inputs finite as well: **1999 per cell in both** (the identity table below carries n detected / n A5 / n analysed). Winner-only and winner-floor excluded. The rendered `summary_complement_variance.html` (committed) is the document of record; every number below is its "Numbers for the record (A6)" section verbatim. The A5 tables re-rendered in the same pass reproduce `REPORT_complement_location_2026-09-08.md` (checked on the eps 0.20 HR 1.50 T3 rows: identical).

**Definitions (per replicate, log scale).** `a_H = log(nv_H_est) − log(betaHhat_H)` (naive error on Ĥ; positive = harm over-estimated), `c_H = log(nv_H_est) − log(mr_H_est)` (two-term correction), `e_H = a_H − c_H`, `lam_H = fld_H_lam_mean`, `ef_H = log(fld_H_est2) − log(betaHhat_H)`. From the Stage 0 quotes, `ef_H = e_H − lam_H` on the working scale; asserted `<= 1e-12` in the document, realized 2.2e-16 to 3.3e-16 (table below). **Observed coverage** is that of the harm-side one-sided **lower** bounds: field `mean(betaHhat_H >= fld_H_lo1s)` (the stored quantile-form bound) and IJ two-term `mean(betaHhat_H >= exp(log(mr_H_est) − z₀.₉₅·mr_H_se_ij))`; Wilson 95% limits. **Gaussian-implied field coverage** is Φ(z₀.₉₅·r − b) with b = mean(ef_H)/SD(ef_H), r = mean(fld_H_se)/SD(ef_H) from the stratum's own numbers (the error SD; the marginal SD of `log fld_H_est2` is shown beside it). The Gaussian-implied value is the SE-form analogue of a bound that is realized in quantile form. Tertiles within the cell as in A5.

**Reproduction of the E1 constructions table (harm block, cell level).** The "all" rows below reproduce `REPORT_field_studentize_e1_2026-09-08.md`'s `Hhat (lower)` rows: field one-sided coverage 0.974 / 0.970 / 0.980 / 0.976 / 0.947 / 0.936 and IJ two-term 0.985 / 0.985 / 0.995 / 0.994 / 0.999 / 0.971 (cells in the order effMaxSG ε 0.20 HR 1.50, HR 1.75, ε 0.30 HR 1.50, HR 1.75, maxSG, minSG); the field's bias (log) −0.056 / −0.078 / −0.093 / −0.097 / −0.005 / 0.022 = mean `ef_H` −0.0555 / −0.0783 / −0.0930 / −0.0969 / −0.0048 / 0.0221; the naive bias 0.430 / 0.365 / 0.299 / 0.254 / 0.046 / 0.326 = mean `a_H`; the IJ bias −0.003 / −0.037 / −0.068 / −0.082 / 0.002 / 0.074 = mean `e_H`; error SD 0.314 and marginal SD 0.299 (first cell) = 0.3145 / 0.2994; Wilson limits identical.

### Tables (verbatim from the rendered `summary_complement_variance.html`, section "Numbers for the record (A6)")

Cells labelled as the document labels them: `HR x n500 | e1stud (focus [eps], J 10, prev 31%)`. Order: effMaxSG ε 0.20 (HR 1.50, 1.75) → ε 0.30 (HR 1.50, 1.75) → maxSG → minSG. Mean (SE) on the log scale; "all" is the whole cell. In the H3 fraction table the complement rows (`a`, `ef_s`, `ef_s / a`) are A5's L3 values on the same replicates.



**Identity check (harm block)** (max |ef_H - (e_H - lam_H)|; n detected / n in A5 / n analysed):

| cell | regime | detected | n A5 | n | max abs ef_H |
|----|----|----|----|----|----|
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | 1999 | 1999 | 1999 | 3.3e-16 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | 1999 | 1999 | 1999 | 3.3e-16 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | 1999 | 1999 | 1999 | 2.8e-16 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | 1999 | 1999 | 1999 | 3.3e-16 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | 1999 | 1999 | 1999 | 2.2e-16 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | 1999 | 1999 | 1999 | 2.8e-16 |

**H1 -- harm-block location by p-hat tertile** (mean (SE) on the log scale):

| cell | regime | stratum | n | mean p-hat | mean Hhat/H | a_H | c_H | e_H | lam_H | ef_H |
|----|----|----|----|----|----|----|----|----|----|----|
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat all | 1999 | 0.103 | 0.834 | 0.4305 (0.0063) | 0.4336 (0.0017) | -0.0031 (0.0068) | 0.0524 (0.0009) | -0.0555 (0.0070) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T1 [0.000, 0.043] | 669 | 0.022 | 1.009 | 0.3115 (0.0098) | 0.4717 (0.0024) | -0.1602 (0.0103) | 0.0685 (0.0013) | -0.2287 (0.0105) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T2 [0.044, 0.109] | 664 | 0.071 | 0.814 | 0.3960 (0.0097) | 0.4387 (0.0026) | -0.0427 (0.0099) | 0.0541 (0.0014) | -0.0968 (0.0101) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T3 [0.109, 0.818] | 666 | 0.215 | 0.677 | 0.5844 (0.0103) | 0.3902 (0.0029) | 0.1942 (0.0106) | 0.0346 (0.0016) | 0.1596 (0.0108) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat all | 1999 | 0.119 | 0.846 | 0.3647 (0.0066) | 0.4022 (0.0018) | -0.0374 (0.0071) | 0.0408 (0.0009) | -0.0783 (0.0074) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T1 [0.000, 0.051] | 667 | 0.026 | 1.009 | 0.2511 (0.0108) | 0.4489 (0.0025) | -0.1978 (0.0112) | 0.0612 (0.0013) | -0.2591 (0.0114) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T2 [0.051, 0.129] | 666 | 0.085 | 0.828 | 0.3291 (0.0097) | 0.4038 (0.0026) | -0.0747 (0.0098) | 0.0409 (0.0014) | -0.1156 (0.0099) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T3 [0.129, 0.928] | 666 | 0.245 | 0.700 | 0.5141 (0.0112) | 0.3537 (0.0029) | 0.1604 (0.0115) | 0.0203 (0.0016) | 0.1401 (0.0116) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat all | 1999 | 0.064 | 1.107 | 0.2989 (0.0055) | 0.3667 (0.0017) | -0.0678 (0.0058) | 0.0252 (0.0009) | -0.0930 (0.0060) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T1 [0.000, 0.027] | 671 | 0.013 | 1.315 | 0.2344 (0.0080) | 0.4007 (0.0026) | -0.1663 (0.0081) | 0.0386 (0.0014) | -0.2049 (0.0083) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T2 [0.027, 0.065] | 662 | 0.043 | 1.078 | 0.2591 (0.0088) | 0.3671 (0.0027) | -0.1080 (0.0091) | 0.0251 (0.0014) | -0.1331 (0.0094) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T3 [0.066, 0.687] | 666 | 0.135 | 0.926 | 0.4036 (0.0101) | 0.3321 (0.0027) | 0.0715 (0.0103) | 0.0118 (0.0014) | 0.0597 (0.0104) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat all | 1999 | 0.080 | 1.094 | 0.2543 (0.0057) | 0.3365 (0.0017) | -0.0822 (0.0060) | 0.0148 (0.0009) | -0.0969 (0.0062) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T1 [0.000, 0.035] | 667 | 0.017 | 1.291 | 0.1844 (0.0085) | 0.3762 (0.0028) | -0.1918 (0.0087) | 0.0304 (0.0014) | -0.2223 (0.0089) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T2 [0.035, 0.083] | 666 | 0.056 | 1.063 | 0.2279 (0.0093) | 0.3381 (0.0027) | -0.1102 (0.0094) | 0.0141 (0.0014) | -0.1243 (0.0095) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T3 [0.084, 0.872] | 666 | 0.167 | 0.927 | 0.3508 (0.0104) | 0.2950 (0.0027) | 0.0557 (0.0106) | -0.0002 (0.0014) | 0.0560 (0.0107) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat all | 1999 | 0.415 | 2.633 | 0.0465 (0.0021) | 0.0445 (0.0008) | 0.0020 (0.0027) | 0.0068 (0.0003) | -0.0048 (0.0029) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat T1 [0.005, 0.123] | 667 | 0.061 | 2.197 | -0.0171 (0.0036) | 0.0791 (0.0012) | -0.0962 (0.0043) | 0.0175 (0.0005) | -0.1138 (0.0045) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat T2 [0.123, 0.644] | 666 | 0.365 | 2.759 | 0.0238 (0.0017) | 0.0421 (0.0008) | -0.0183 (0.0022) | 0.0050 (0.0004) | -0.0233 (0.0024) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat T3 [0.644, 1.000] | 666 | 0.821 | 2.943 | 0.1329 (0.0026) | 0.0122 (0.0003) | 0.1206 (0.0029) | -0.0021 (0.0002) | 0.1227 (0.0029) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat all | 1999 | 0.086 | 0.406 | 0.3262 (0.0061) | 0.2519 (0.0022) | 0.0743 (0.0068) | 0.0521 (0.0010) | 0.0221 (0.0073) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat T1 [0.000, 0.014] | 669 | 0.005 | 0.405 | 0.3750 (0.0101) | 0.1891 (0.0036) | 0.1859 (0.0112) | 0.0268 (0.0016) | 0.1591 (0.0119) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat T2 [0.014, 0.078] | 664 | 0.037 | 0.407 | 0.3172 (0.0108) | 0.2686 (0.0035) | 0.0486 (0.0117) | 0.0579 (0.0017) | -0.0093 (0.0124) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat T3 [0.078, 0.867] | 666 | 0.216 | 0.406 | 0.2861 (0.0105) | 0.2983 (0.0028) | -0.0123 (0.0112) | 0.0719 (0.0015) | -0.0841 (0.0117) |

**H1 -- harm-block coverage by p-hat tertile** (one-sided lower bounds; Gaussian-implied from the stratum’s own mean ef_H, SD ef_H, mean fld_H_se):

| cell | regime | stratum | n | mean naive SE | mean fld_H_se | mean se_ij | SD ef_H (error) | SD log fld_H_est2 (marginal) | mean ef_H | b | r | field cov [Wilson] | IJ cov [Wilson] | Gaussian-implied (field) |
|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat all | 1999 | 0.2497 | 0.3237 | 0.3600 | 0.3145 | 0.2994 | -0.0555 | -0.177 | 1.029 | 0.974 [0.967, 0.981] | 0.985 [0.979, 0.990] | 0.969 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T1 [0.000, 0.043] | 669 | 0.2253 | 0.3311 | 0.3465 | 0.2727 | 0.2251 | -0.2287 | -0.839 | 1.214 | 0.999 [0.992, 1.000] | 1.000 [0.994, 1.000] | 0.998 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T2 [0.044, 0.109] | 664 | 0.2477 | 0.3190 | 0.3561 | 0.2592 | 0.2220 | -0.0968 | -0.373 | 1.230 | 0.995 [0.987, 0.998] | 0.998 [0.992, 1.000] | 0.992 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T3 [0.109, 0.818] | 666 | 0.2760 | 0.3208 | 0.3775 | 0.2779 | 0.2625 | 0.1596 | 0.574 | 1.155 | 0.929 [0.907, 0.947] | 0.958 [0.940, 0.971] | 0.907 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat all | 1999 | 0.2425 | 0.3229 | 0.3565 | 0.3286 | 0.3119 | -0.0783 | -0.238 | 0.983 | 0.970 [0.962, 0.977] | 0.985 [0.979, 0.989] | 0.968 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T1 [0.000, 0.051] | 667 | 0.2213 | 0.3320 | 0.3453 | 0.2956 | 0.2341 | -0.2591 | -0.876 | 1.123 | 0.999 [0.992, 1.000] | 1.000 [0.994, 1.000] | 0.997 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T2 [0.051, 0.129] | 666 | 0.2405 | 0.3164 | 0.3525 | 0.2563 | 0.2189 | -0.1156 | -0.451 | 1.235 | 0.992 [0.983, 0.997] | 0.995 [0.987, 0.998] | 0.993 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T3 [0.129, 0.928] | 666 | 0.2658 | 0.3203 | 0.3716 | 0.2988 | 0.2792 | 0.1401 | 0.469 | 1.072 | 0.919 [0.896, 0.937] | 0.959 [0.942, 0.972] | 0.902 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat all | 1999 | 0.2130 | 0.2979 | 0.3336 | 0.2670 | 0.2701 | -0.0930 | -0.348 | 1.116 | 0.980 [0.973, 0.986] | 0.995 [0.991, 0.997] | 0.985 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T1 [0.000, 0.027] | 671 | 0.1944 | 0.3067 | 0.3241 | 0.2162 | 0.2006 | -0.2049 | -0.948 | 1.419 | 1.000 [0.994, 1.000] | 1.000 [0.994, 1.000] | 0.999 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T2 [0.027, 0.065] | 662 | 0.2129 | 0.2943 | 0.3310 | 0.2414 | 0.2194 | -0.1331 | -0.551 | 1.219 | 0.992 [0.982, 0.997] | 1.000 [0.994, 1.000] | 0.995 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T3 [0.066, 0.687] | 666 | 0.2317 | 0.2927 | 0.3457 | 0.2676 | 0.2519 | 0.0597 | 0.223 | 1.094 | 0.949 [0.930, 0.963] | 0.985 [0.973, 0.992] | 0.942 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat all | 1999 | 0.2098 | 0.2970 | 0.3333 | 0.2765 | 0.2862 | -0.0969 | -0.351 | 1.074 | 0.976 [0.968, 0.982] | 0.994 [0.990, 0.997] | 0.983 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T1 [0.000, 0.035] | 667 | 0.1925 | 0.3063 | 0.3243 | 0.2295 | 0.2092 | -0.2223 | -0.968 | 1.334 | 0.997 [0.989, 0.999] | 1.000 [0.994, 1.000] | 0.999 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T2 [0.035, 0.083] | 666 | 0.2104 | 0.2929 | 0.3308 | 0.2462 | 0.2081 | -0.1243 | -0.505 | 1.190 | 0.992 [0.983, 0.997] | 0.998 [0.992, 1.000] | 0.993 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T3 [0.084, 0.872] | 666 | 0.2264 | 0.2919 | 0.3448 | 0.2764 | 0.2645 | 0.0560 | 0.203 | 1.056 | 0.938 [0.918, 0.954] | 0.983 [0.971, 0.991] | 0.938 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat all | 1999 | 0.1333 | 0.1261 | 0.2207 | 0.1312 | 0.0867 | -0.0048 | -0.037 | 0.961 | 0.947 [0.936, 0.956] | 0.999 [0.997, 1.000] | 0.947 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat T1 [0.005, 0.123] | 667 | 0.1461 | 0.1298 | 0.2163 | 0.1171 | 0.0353 | -0.1138 | -0.971 | 1.108 | 0.999 [0.992, 1.000] | 1.000 [0.994, 1.000] | 0.997 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat T2 [0.123, 0.644] | 666 | 0.1298 | 0.1254 | 0.2163 | 0.0607 | 0.0403 | -0.0233 | -0.383 | 2.066 | 1.000 [0.994, 1.000] | 1.000 [0.994, 1.000] | 1.000 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat T3 [0.644, 1.000] | 666 | 0.1239 | 0.1230 | 0.2296 | 0.0759 | 0.0755 | 0.1227 | 1.618 | 1.621 | 0.842 [0.813, 0.868] | 0.998 [0.992, 1.000] | 0.853 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat all | 1999 | 0.3200 | 0.3151 | 0.3907 | 0.3257 | 0.2838 | 0.0221 | 0.068 | 0.967 | 0.936 [0.924, 0.946] | 0.971 [0.963, 0.977] | 0.936 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat T1 [0.000, 0.014] | 669 | 0.3213 | 0.3238 | 0.4026 | 0.3067 | 0.2723 | 0.1591 | 0.519 | 1.056 | 0.892 [0.867, 0.914] | 0.955 [0.937, 0.968] | 0.888 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat T2 [0.014, 0.078] | 664 | 0.3198 | 0.3077 | 0.3738 | 0.3202 | 0.2850 | -0.0093 | -0.029 | 0.961 | 0.952 [0.933, 0.966] | 0.973 [0.958, 0.983] | 0.946 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat T3 [0.078, 0.867] | 666 | 0.3187 | 0.3137 | 0.3957 | 0.3016 | 0.2670 | -0.0841 | -0.279 | 1.040 | 0.964 [0.947, 0.976] | 0.985 [0.973, 0.992] | 0.977 |

**H2 -- harm-block location by |Hhat|/|H| tertile** (mean (SE) on the log scale):

| cell | regime | stratum | n | mean p-hat | mean Hhat/H | a_H | c_H | e_H | lam_H | ef_H |
|----|----|----|----|----|----|----|----|----|----|----|
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H all | 1999 | 0.103 | 0.834 | 0.4305 (0.0063) | 0.4336 (0.0017) | -0.0031 (0.0068) | 0.0524 (0.0009) | -0.0555 (0.0070) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T1 [0.372, 0.690] | 667 | 0.162 | 0.569 | 0.5851 (0.0105) | 0.4295 (0.0030) | 0.1556 (0.0114) | 0.0550 (0.0015) | 0.1006 (0.0120) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T2 [0.690, 0.946] | 666 | 0.093 | 0.806 | 0.4390 (0.0103) | 0.4311 (0.0030) | 0.0079 (0.0111) | 0.0516 (0.0015) | -0.0437 (0.0115) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T3 [0.947, 2.085] | 666 | 0.052 | 1.127 | 0.2672 (0.0080) | 0.4402 (0.0028) | -0.1730 (0.0089) | 0.0507 (0.0015) | -0.2237 (0.0096) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H all | 1999 | 0.119 | 0.846 | 0.3647 (0.0066) | 0.4022 (0.0018) | -0.0374 (0.0071) | 0.0408 (0.0009) | -0.0783 (0.0074) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T1 [0.350, 0.705] | 667 | 0.184 | 0.578 | 0.5347 (0.0110) | 0.3961 (0.0031) | 0.1386 (0.0121) | 0.0432 (0.0015) | 0.0954 (0.0126) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T2 [0.705, 0.972] | 666 | 0.110 | 0.833 | 0.3515 (0.0111) | 0.3985 (0.0032) | -0.0470 (0.0115) | 0.0384 (0.0017) | -0.0854 (0.0119) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T3 [0.972, 1.929] | 666 | 0.062 | 1.127 | 0.2077 (0.0080) | 0.4119 (0.0029) | -0.2042 (0.0092) | 0.0408 (0.0016) | -0.2451 (0.0100) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H all | 1999 | 0.064 | 1.107 | 0.2989 (0.0055) | 0.3667 (0.0017) | -0.0678 (0.0058) | 0.0252 (0.0009) | -0.0930 (0.0060) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T1 [0.372, 0.968] | 667 | 0.096 | 0.772 | 0.4313 (0.0108) | 0.3726 (0.0028) | 0.0587 (0.0113) | 0.0336 (0.0013) | 0.0250 (0.0116) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T2 [0.969, 1.235] | 666 | 0.066 | 1.098 | 0.2446 (0.0081) | 0.3604 (0.0029) | -0.1158 (0.0088) | 0.0213 (0.0015) | -0.1371 (0.0093) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T3 [1.236, 2.405] | 666 | 0.029 | 1.451 | 0.2208 (0.0068) | 0.3671 (0.0029) | -0.1463 (0.0076) | 0.0207 (0.0015) | -0.1670 (0.0082) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H all | 1999 | 0.080 | 1.094 | 0.2543 (0.0057) | 0.3365 (0.0017) | -0.0822 (0.0060) | 0.0148 (0.0009) | -0.0969 (0.0062) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T1 [0.393, 0.982] | 667 | 0.116 | 0.783 | 0.3765 (0.0117) | 0.3429 (0.0029) | 0.0336 (0.0122) | 0.0239 (0.0014) | 0.0097 (0.0125) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T2 [0.982, 1.217] | 666 | 0.087 | 1.094 | 0.1906 (0.0083) | 0.3267 (0.0031) | -0.1361 (0.0090) | 0.0092 (0.0015) | -0.1453 (0.0095) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T3 [1.218, 2.354] | 666 | 0.036 | 1.404 | 0.1956 (0.0071) | 0.3398 (0.0030) | -0.1442 (0.0079) | 0.0112 (0.0015) | -0.1554 (0.0085) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H all | 1999 | 0.415 | 2.633 | 0.0465 (0.0021) | 0.0445 (0.0008) | 0.0020 (0.0027) | 0.0068 (0.0003) | -0.0048 (0.0029) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H T1 [0.454, 2.519] | 667 | 0.096 | 2.085 | -0.0225 (0.0037) | 0.0808 (0.0012) | -0.1033 (0.0044) | 0.0174 (0.0005) | -0.1207 (0.0046) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H T2 [2.520, 2.892] | 678 | 0.489 | 2.733 | 0.0734 (0.0030) | 0.0311 (0.0008) | 0.0423 (0.0036) | 0.0032 (0.0003) | 0.0390 (0.0038) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H T3 [2.893, 3.640] | 654 | 0.665 | 3.088 | 0.0889 (0.0027) | 0.0213 (0.0006) | 0.0676 (0.0031) | -0.0002 (0.0003) | 0.0679 (0.0033) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H all | 1999 | 0.086 | 0.406 | 0.3262 (0.0061) | 0.2519 (0.0022) | 0.0743 (0.0068) | 0.0521 (0.0010) | 0.0221 (0.0073) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H T1 [0.330, 0.391] | 679 | 0.083 | 0.375 | 0.3624 (0.0107) | 0.2342 (0.0036) | 0.1281 (0.0118) | 0.0455 (0.0016) | 0.0826 (0.0125) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H T2 [0.391, 0.416] | 654 | 0.088 | 0.404 | 0.3147 (0.0106) | 0.2522 (0.0038) | 0.0625 (0.0120) | 0.0532 (0.0018) | 0.0093 (0.0128) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H T3 [0.416, 1.042] | 666 | 0.087 | 0.440 | 0.3006 (0.0101) | 0.2696 (0.0038) | 0.0310 (0.0114) | 0.0579 (0.0018) | -0.0269 (0.0122) |

**H2 -- harm-block coverage by |Hhat|/|H| tertile** (one-sided lower bounds; Gaussian-implied from the stratum’s own mean ef_H, SD ef_H, mean fld_H_se):

| cell | regime | stratum | n | mean naive SE | mean fld_H_se | mean se_ij | SD ef_H (error) | SD log fld_H_est2 (marginal) | mean ef_H | b | r | field cov [Wilson] | IJ cov [Wilson] | Gaussian-implied (field) |
|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H all | 1999 | 0.2497 | 0.3237 | 0.3600 | 0.3145 | 0.2994 | -0.0555 | -0.177 | 1.029 | 0.974 [0.967, 0.981] | 0.985 [0.979, 0.990] | 0.969 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T1 [0.372, 0.690] | 667 | 0.2916 | 0.3207 | 0.3898 | 0.3087 | 0.3051 | 0.1006 | 0.326 | 1.039 | 0.945 [0.924, 0.959] | 0.963 [0.945, 0.974] | 0.917 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T2 [0.690, 0.946] | 666 | 0.2520 | 0.3275 | 0.3654 | 0.2960 | 0.2731 | -0.0437 | -0.148 | 1.106 | 0.979 [0.965, 0.987] | 0.994 [0.985, 0.998] | 0.975 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T3 [0.947, 2.085] | 666 | 0.2053 | 0.3227 | 0.3248 | 0.2473 | 0.2519 | -0.2237 | -0.904 | 1.305 | 1.000 [0.994, 1.000] | 1.000 [0.994, 1.000] | 0.999 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H all | 1999 | 0.2425 | 0.3229 | 0.3565 | 0.3286 | 0.3119 | -0.0783 | -0.238 | 0.983 | 0.970 [0.962, 0.977] | 0.985 [0.979, 0.989] | 0.968 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T1 [0.350, 0.705] | 667 | 0.2855 | 0.3207 | 0.3886 | 0.3252 | 0.3239 | 0.0954 | 0.293 | 0.986 | 0.931 [0.909, 0.948] | 0.964 [0.947, 0.976] | 0.908 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T2 [0.705, 0.972] | 666 | 0.2422 | 0.3255 | 0.3588 | 0.3062 | 0.2764 | -0.0854 | -0.279 | 1.063 | 0.982 [0.969, 0.990] | 0.994 [0.985, 0.998] | 0.979 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T3 [0.972, 1.929] | 666 | 0.1998 | 0.3225 | 0.3220 | 0.2582 | 0.2616 | -0.2451 | -0.949 | 1.249 | 0.997 [0.989, 0.999] | 0.997 [0.989, 0.999] | 0.999 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H all | 1999 | 0.2130 | 0.2979 | 0.3336 | 0.2670 | 0.2701 | -0.0930 | -0.348 | 1.116 | 0.980 [0.973, 0.986] | 0.995 [0.991, 0.997] | 0.985 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T1 [0.372, 0.968] | 667 | 0.2519 | 0.3009 | 0.3647 | 0.3001 | 0.2820 | 0.0250 | 0.083 | 1.002 | 0.946 [0.926, 0.961] | 0.988 [0.977, 0.994] | 0.941 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T2 [0.969, 1.235] | 666 | 0.2050 | 0.2988 | 0.3261 | 0.2401 | 0.2380 | -0.1371 | -0.571 | 1.244 | 0.995 [0.987, 0.998] | 0.997 [0.989, 0.999] | 0.996 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T3 [1.236, 2.405] | 666 | 0.1819 | 0.2942 | 0.3099 | 0.2122 | 0.2181 | -0.1670 | -0.787 | 1.387 | 1.000 [0.994, 1.000] | 1.000 [0.994, 1.000] | 0.999 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H all | 1999 | 0.2098 | 0.2970 | 0.3333 | 0.2765 | 0.2862 | -0.0969 | -0.351 | 1.074 | 0.976 [0.968, 0.982] | 0.994 [0.990, 0.997] | 0.983 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T1 [0.393, 0.982] | 667 | 0.2457 | 0.2999 | 0.3618 | 0.3233 | 0.3025 | 0.0097 | 0.030 | 0.928 | 0.940 [0.919, 0.956] | 0.987 [0.975, 0.993] | 0.933 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T2 [0.982, 1.217] | 666 | 0.2017 | 0.2962 | 0.3251 | 0.2450 | 0.2474 | -0.1453 | -0.593 | 1.209 | 0.991 [0.980, 0.996] | 0.995 [0.987, 0.998] | 0.995 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T3 [1.218, 2.354] | 666 | 0.1819 | 0.2949 | 0.3130 | 0.2187 | 0.2332 | -0.1554 | -0.710 | 1.348 | 0.997 [0.989, 0.999] | 1.000 [0.994, 1.000] | 0.998 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H all | 1999 | 0.1333 | 0.1261 | 0.2207 | 0.1312 | 0.0867 | -0.0048 | -0.037 | 0.961 | 0.947 [0.936, 0.956] | 0.999 [0.997, 1.000] | 0.947 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H T1 [0.454, 2.519] | 667 | 0.1477 | 0.1334 | 0.2203 | 0.1193 | 0.0482 | -0.1207 | -1.012 | 1.118 | 0.999 [0.992, 1.000] | 1.000 [0.994, 1.000] | 0.998 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H T2 [2.520, 2.892] | 678 | 0.1274 | 0.1224 | 0.2193 | 0.0994 | 0.0924 | 0.0390 | 0.393 | 1.232 | 0.928 [0.906, 0.945] | 0.999 [0.992, 1.000] | 0.949 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H T3 [2.893, 3.640] | 654 | 0.1247 | 0.1224 | 0.2226 | 0.0832 | 0.0832 | 0.0679 | 0.816 | 1.470 | 0.914 [0.890, 0.933] | 1.000 [0.994, 1.000] | 0.946 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H all | 1999 | 0.3200 | 0.3151 | 0.3907 | 0.3257 | 0.2838 | 0.0221 | 0.068 | 0.967 | 0.936 [0.924, 0.946] | 0.971 [0.963, 0.977] | 0.936 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H T1 [0.330, 0.391] | 679 | 0.3213 | 0.3147 | 0.3945 | 0.3261 | 0.2870 | 0.0826 | 0.253 | 0.965 | 0.909 [0.885, 0.928] | 0.959 [0.941, 0.971] | 0.909 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H T2 [0.391, 0.416] | 654 | 0.3190 | 0.3162 | 0.3910 | 0.3273 | 0.2844 | 0.0093 | 0.028 | 0.966 | 0.936 [0.914, 0.952] | 0.977 [0.963, 0.986] | 0.941 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H T3 [0.416, 1.042] | 666 | 0.3196 | 0.3144 | 0.3866 | 0.3145 | 0.2732 | -0.0269 | -0.086 | 1.000 | 0.964 [0.947, 0.976] | 0.977 [0.963, 0.986] | 0.958 |

**H3 -- correlations of p-hat with the naive errors** (Pearson / Spearman):

| cell | regime | n | corr(p-hat, a_H) Pearson | Spearman | corr(p-hat, a) Pearson | Spearman |
|----|----|----|----|----|----|----|
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | 1999 | 0.459 | 0.434 | -0.084 | -0.086 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | 1999 | 0.438 | 0.415 | -0.111 | -0.130 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | 1999 | 0.372 | 0.310 | -0.102 | -0.106 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | 1999 | 0.343 | 0.289 | -0.102 | -0.106 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | 1999 | 0.693 | 0.711 | 0.348 | 0.301 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | 1999 | -0.079 | -0.145 | -0.200 | -0.290 |

**H3 -- uncorrected fractions, p-hat tertiles** (T1 low … T3 high p-hat):

| cell | regime | quantity | all | T1 | T2 | T3 |
|----|----|----|----|----|----|----|
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | a_H (harm) | 0.4305 | 0.3115 | 0.3960 | 0.5844 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | ef_H (harm) | -0.0555 | -0.2287 | -0.0968 | 0.1596 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | ef_H / a_H (harm) | -0.1289 | -0.7340 | -0.2444 | 0.2731 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | a (complement) | -0.1234 | -0.1111 | -0.1228 | -0.1362 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | ef_s (complement) | -0.0308 | -0.0099 | -0.0272 | -0.0555 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | ef_s / a (complement) | 0.2498 | 0.0887 | 0.2215 | 0.4072 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | a_H (harm) | 0.3647 | 0.2511 | 0.3291 | 0.5141 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | ef_H (harm) | -0.0783 | -0.2591 | -0.1156 | 0.1401 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | ef_H / a_H (harm) | -0.2146 | -1.0318 | -0.3511 | 0.2726 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | a (complement) | -0.1021 | -0.0768 | -0.1091 | -0.1205 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | ef_s (complement) | -0.0193 | 0.0162 | -0.0225 | -0.0518 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | ef_s / a (complement) | 0.1894 | -0.2110 | 0.2061 | 0.4300 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | a_H (harm) | 0.2989 | 0.2344 | 0.2591 | 0.4036 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | ef_H (harm) | -0.0930 | -0.2049 | -0.1331 | 0.0597 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | ef_H / a_H (harm) | -0.3110 | -0.8742 | -0.5138 | 0.1480 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | a (complement) | -0.1367 | -0.1168 | -0.1365 | -0.1568 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | ef_s (complement) | -0.0316 | -0.0009 | -0.0285 | -0.0657 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | ef_s / a (complement) | 0.2313 | 0.0075 | 0.2088 | 0.4189 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | a_H (harm) | 0.2543 | 0.1844 | 0.2279 | 0.3508 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | ef_H (harm) | -0.0969 | -0.2223 | -0.1243 | 0.0560 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | ef_H / a_H (harm) | -0.3811 | -1.2056 | -0.5454 | 0.1596 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | a (complement) | -0.1135 | -0.0968 | -0.1094 | -0.1343 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | ef_s (complement) | -0.0217 | 0.0083 | -0.0140 | -0.0595 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | ef_s / a (complement) | 0.1913 | -0.0859 | 0.1277 | 0.4430 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | a_H (harm) | 0.0465 | -0.0171 | 0.0238 | 0.1329 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | ef_H (harm) | -0.0048 | -0.1138 | -0.0233 | 0.1227 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | ef_H / a_H (harm) | -0.1036 | 6.6424 | -0.9770 | 0.9238 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | a (complement) | -0.1914 | -0.3061 | -0.2549 | -0.0131 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | ef_s (complement) | -0.0392 | -0.0829 | -0.0744 | 0.0398 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | ef_s / a (complement) | 0.2048 | 0.2709 | 0.2918 | -3.0399 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | a_H (harm) | 0.3262 | 0.3750 | 0.3172 | 0.2861 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | ef_H (harm) | 0.0221 | 0.1591 | -0.0093 | -0.0841 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | ef_H / a_H (harm) | 0.0679 | 0.4243 | -0.0292 | -0.2940 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | a (complement) | -0.0488 | -0.0027 | -0.0533 | -0.0904 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | ef_s (complement) | -0.0129 | 0.0221 | -0.0148 | -0.0462 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | ef_s / a (complement) | 0.2647 | -8.0868 | 0.2779 | 0.5106 |


## H4 — Reading per cell (in the record, not a task)

The classification asked for: **(A)** competition collapse on both blocks — `a_H` largest in the high-p̂ tertile while `c_H` and `lam_H` shrink, `ef_H` clearly positive there, harm-side coverage lowest there; **(B)** complement-only — `ef_H` ≈ 0 and harm coverage flat across p̂ tertiles while the complement's residual persists; **(C)** something else. Numbers are the H1 / H3 p̂-tertile rows (T1 / T2 / T3).

- **effMaxSG ε 0.20, HR 1.50 — (A).** `a_H` 0.3115 / 0.3960 / 0.5844 (largest in T3); `c_H` 0.4717 / 0.4387 / 0.3902 and `lam_H` 0.0685 / 0.0541 / 0.0346 (both shrink); `ef_H` −0.2287 / −0.0968 / **+0.1596** (SE 0.0108); field lower coverage 0.999 / 0.995 / **0.929** [0.907, 0.947], IJ 1.000 / 0.998 / 0.958. Complement (A5): `ef_s / a` 0.089 / 0.222 / 0.407. Uncorrected fraction on the harm block `ef_H / a_H` in T3: 0.273. corr(p̂, `a_H`) 0.459 Pearson / 0.434 Spearman; corr(p̂, `a`) −0.084 / −0.086 (the complement's `a` is negative; its magnitude also grows with p̂).
- **effMaxSG ε 0.20, HR 1.75 — (A).** `a_H` 0.2511 / 0.3291 / 0.5141; `c_H` 0.4489 / 0.4038 / 0.3537; `lam_H` 0.0612 / 0.0409 / 0.0203; `ef_H` −0.2591 / −0.1156 / **+0.1401** (0.0116); field coverage 0.999 / 0.992 / **0.919** [0.896, 0.937], IJ 1.000 / 0.995 / 0.959. `ef_H / a_H` T3 0.273; complement `ef_s / a` T3 0.430. corr(p̂, `a_H`) 0.438 / 0.415; corr(p̂, `a`) −0.111 / −0.130.
- **effMaxSG ε 0.30, HR 1.50 — (A), smaller amplitude.** `a_H` 0.2344 / 0.2591 / 0.4036; `c_H` 0.4007 / 0.3671 / 0.3321; `lam_H` 0.0386 / 0.0251 / 0.0118; `ef_H` −0.2049 / −0.1331 / **+0.0597** (0.0104); field coverage 1.000 / 0.992 / **0.949** [0.930, 0.963], IJ 1.000 / 1.000 / 0.985. `ef_H / a_H` T3 0.148; complement T3 0.419. corr(p̂, `a_H`) 0.372 / 0.310; corr(p̂, `a`) −0.102 / −0.106.
- **effMaxSG ε 0.30, HR 1.75 — (A), smaller amplitude.** `a_H` 0.1844 / 0.2279 / 0.3508; `c_H` 0.3762 / 0.3381 / 0.2950; `lam_H` 0.0304 / 0.0141 / −0.0002; `ef_H` −0.2223 / −0.1243 / **+0.0560** (0.0107); field coverage 0.997 / 0.992 / **0.938** [0.918, 0.954], IJ 1.000 / 0.998 / 0.983. `ef_H / a_H` T3 0.160; complement T3 0.443. corr(p̂, `a_H`) 0.343 / 0.289; corr(p̂, `a`) −0.102 / −0.106.
- **maxSG, HR 1.75 — (C): the harm block alone shows the high-p̂ collapse; the complement's residual sits in T1–T2.** Harm block: `a_H` −0.0171 / 0.0238 / 0.1329; `c_H` 0.0791 / 0.0421 / 0.0122; `lam_H` 0.0175 / 0.0050 / −0.0021; `ef_H` −0.1138 / −0.0233 / **+0.1227** (0.0029); `ef_H / a_H` T3 0.924 (the correction is nearly absent at high p̂); field coverage 0.999 / 1.000 / **0.842** [0.813, 0.868]; IJ 1.000 / 1.000 / 0.998 (mean `mr_H_se_ij` 0.2296 against an error SD of 0.0759 in T3). Complement (A5): `ef_s` −0.0829 / −0.0744 / +0.0398, over-corrected in T3 (coverage 0.964) with its residual in T1–T2. The harm block satisfies every condition of (A) on its own, but the two blocks' residuals sit in different strata, so the pair is not "both blocks" in the sense of (A). corr(p̂, `a_H`) 0.693 / 0.711; corr(p̂, `a`) 0.348 / 0.301 (both positive here; the complement's `a` is least negative at high p̂).
- **minSG, HR 1.75 — (C): the harm-side residual sits in the LOW-p̂ tertile, the complement's in the high one.** `a_H` 0.3750 / 0.3172 / 0.2861 (largest in T1, not T3); `c_H` 0.1891 / 0.2686 / 0.2983 and `lam_H` 0.0268 / 0.0579 / 0.0719 (both grow with p̂); `ef_H` **+0.1591** / −0.0093 / −0.0841 (SE 0.0119 / 0.0124 / 0.0117); field coverage **0.892** [0.867, 0.914] / 0.952 / 0.964, IJ 0.955 / 0.973 / 0.985 — lowest in T1. Complement (A5): `ef_s / a` 0 → 0.28 → 0.51, residual in T3. corr(p̂, `a_H`) −0.079 / −0.145; corr(p̂, `a`) −0.200 / −0.290. Neither (A) nor (B): an inverted (A) on the harm block.

**Across the four band cells** the harm block shows the same shape as the complement — optimism growing with p̂ while both corrections shrink — with a *larger* T3 residual in absolute terms at ε 0.20 (+0.140 to +0.160 on the harm block against −0.052 to −0.055 on the complement) and a comparable one at ε 0.30 (+0.056 to +0.060 against −0.060 to −0.066). The harm-side uncorrected fraction at T3 (0.15–0.27) is below the complement's (0.41–0.44) because the harm-side optimism `a_H` is 2.2 to 3.6 times the complement's `|a|` at the cell level (2.6 to 4.3 times in T3). The cell-level harm coverage (0.970–0.980) sits above nominal because the T1 stratum is strongly over-corrected (`ef_H` −0.20 to −0.26, coverage 0.997–1.000) and masks the T3 shortfall (0.919–0.949). The IJ two-term lower bound covers ≥ 0.958 in every p̂ stratum of the band cells (its mean SE 0.35–0.38 against an error SD of 0.26–0.30).

**Gaussian-implied vs observed (harm field).** By p̂ tertile, the Gaussian-implied coverage is within 0.022 of the observed in all 24 strata (largest gaps: ε 0.20 HR 1.50 T3 0.907 vs 0.929; ε 0.20 HR 1.75 T3 0.902 vs 0.919; minSG T3 0.977 vs 0.964; maxSG T3 0.853 vs 0.842). By |Ĥ|/|H| tertile it is within 0.032 (largest: maxSG T3 0.946 vs 0.914, ε 0.20 HR 1.50 T1 0.917 vs 0.945, maxSG T2 0.949 vs 0.928, ε 0.20 HR 1.75 T1 0.908 vs 0.931). The Gaussian-implied value is the SE-form analogue; the realized bound is the quantile form (`beta_deb − q95(Λ*)`), which accounts for the residual gaps. It tracks the observed shape in every cell: lowest where the observed is lowest, and the T3 shortfall in the band cells is reproduced from the stratum's own mean and SD of `ef_H`.

**H2 (|Ĥ|/|H| tertiles), briefly.** In the band cells the harm-side residual sits in the low-|Ĥ|/|H| tertile (the pick smaller than H): `ef_H` +0.1006 / −0.0437 / −0.2237 (ε 0.20 HR 1.50), +0.0954 / −0.0854 / −0.2451 (ε 0.20 HR 1.75), +0.0250 / −0.1371 / −0.1670 (ε 0.30 HR 1.50), +0.0097 / −0.1453 / −0.1554 (ε 0.30 HR 1.75), with field coverage 0.945 / 0.931 / 0.946 / 0.940 in T1 against 0.997–1.000 in T3 — the complement's L2 residual `ef_s` in A5 ranged −0.003 to −0.042 across these strata in the band cells with no monotone trend. At maxSG the residual grows with |Ĥ|/|H| (`ef_H` −0.1207 / +0.0390 / +0.0679; coverage 0.999 / 0.928 / 0.914); at minSG it falls (+0.0826 / +0.0093 / −0.0269; coverage 0.909 / 0.936 / 0.964).

## Closing assertions

No compute; no `R/` change. Files changed: `summary_complement_variance.qmd` (A6 appended), `summary_complement_variance.html` (re-rendered on the six `e1stud` bundles), this report. `git status` after the commit shows no new `.rds` under `results/` from this task (the three pre-existing untracked `..._diag_res_*.rds` files and the `diag_h*.html` files under `gbsg_020/`, and the `actg175_extreme_sims_continuous_fixed_10000_payload.rds` under `extreme_subgroups/`, predate this session and were not touched). Branch `feature/glm-extension` left unpushed. Report-and-wait.
