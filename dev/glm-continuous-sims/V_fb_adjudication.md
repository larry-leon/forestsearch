# V — FB adjudication (V1 field provenance, V2 correction mass, V3 family regeneration)

Read-only against `R/`, plus one rerun of replicate 1. No `R/` changes, no DGM
changes. Evidence base: the committed batch bundle
`fb_mr_md_harm/fs_maxeffCons_fb_mr_md40_knoise0_n500_s50_d5000/fs_maxeffCons_fb_mr_md40_knoise0_n500_res_1_50.rds`
(50 replicates, n = 500, nb_boots = 300, mr_draws = 5000) at `993a0586`.

**Headline: V3 falsifies the premise both simulation documents rest on.** The
candidate family is *not* fixed under bootstrap resampling for the consistency
engine. FB and MR therefore do not target the same estimand here, and the FB↔MR
gap reported in the batch is expected behaviour, not a defect in either.

---

## V1 — which `H_estimates` field does `fb_H_est` carry?

**Answer: `H2`, sign-flipped to the oriented scale.** `fb_H_est = -H2`.

The committed qmd, verbatim
(`sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_50.qmd:434-441`):

```r
    if (!is.null(fb) && !is.null(fb$H_estimates)) {
      he <- fb$H_estimates
      .orient <- if (isTRUE(adverse_outcome)) 1 else -1
      rec$fb_H_est <- .orient * (he$H2 %||% NA_real_)
      fb_ci <- sort(.orient * c(he$H2_lower %||% NA_real_,
                                he$H2_upper %||% NA_real_), na.last = TRUE)
      rec$fb_H_lo  <- fb_ci[1]
      rec$fb_H_hi  <- fb_ci[2]
      rec$fb_H_se  <- he$sdH2 %||% NA_real_   # scale-free, NOT negated
    }
```

Source chain for `H2`. `forestsearch_bootstrap_dofuture()` builds `H_estimates`
via `get_dfRes()` (`bootstrap_dofuture_main.R:569-582`), passing
`H2_adj = results$H_biasadj_2`. `get_dfRes()` turns that column into the scalar
`H2` (`bootstrap_calculations_helpers.R:303-316`):

```r
  if (!is.null(H2_adj)) {
    est <- get_targetEst(x = H2_adj, ystar = ystar, cov_method = cov_method, cov_trim = cov_trim)
    q <- est$target_est
    se_new <- est$sehat_new
    cest <- ci_est(x = q, sd = se_new, scale = est.scale, est.loghr = est.loghr)
    H2_lower <- cest$lower
    H2_upper <- cest$upper
    sdH2 <- cest$sd
    H2 <- cest$est
```

and `H_biasadj_2` is the per-resample double correction
(`bootstrap_analysis_dofuture.R:726-728`):

```r
        # Method 1: Simple optimism correction
        H_biasadj_1 <- H_obs - (Hstar_star - Hstar_obs)
        # Method 2: Double correction
        H_biasadj_2 <- 2 * H_obs - (H_star + Hstar_star - Hstar_obs)
```

with the four constituents (`:488`, `:686-697`, `:698-712`):

| term | meaning (source comment) |
|---|---|
| `H_obs` | original subgroup, original data |
| `H_star` | original subgroup H, evaluated on bootstrap data (`:476` "Fit model on H subgroup of bootstrap data") |
| `Hstar_obs` | "New subgroup (from bootstrap) evaluated on ORIGINAL data" (`:686`) |
| `Hstar_star` | "New subgroup (from bootstrap) evaluated on BOOTSTRAP data" (`:698`) |

**`Hstar_*` are the re-selected subgroup.** That is the hinge for V3.

Scale: `est.loghr = FALSE` for MD (`bootstrap_dofuture_main.R:557-563`), so `H2`
and `sdH2` are on the identity MD scale — no log/exp step.

### One replicate's full `H_estimates`, verbatim

Replicate 1, `seed = 8316952`, `nb_boots = 300`, identical config to the batch:

```
--- str(fb$H_estimates) ---
Classes 'data.table' and 'data.frame':	1 obs. of  12 variables:
 $ H0      : num -103
 $ sdH0    : num 31.2
 $ H0_lower: num -164
 $ H0_upper: num -41.9
 $ H1      : num -80.9
 $ sdH1    : num 1.32
 $ H1_lower: num -83.5
 $ H1_upper: num -78.3
 $ H2      : num -81.2
 $ sdH2    : num 41.7
 $ H2_lower: num -163
 $ H2_upper: num 0.626

--- as.data.frame(fb$H_estimates), digits=10 ---
            H0        sdH0     H0_lower     H0_upper           H1        sdH1
1 -103.1366727 31.23977899 -164.3655144 -41.90783098 -80.87327144 1.323439948
      H1_lower    H1_upper          H2        sdH2     H2_lower     H2_upper
1 -83.46716607 -78.2793768 -81.1914644 41.74435021 -163.0088874 0.6259585766

--- names(fb) ---
 [1] "results"        "SG_CIs"         "FSsg_tab"       "Ystar_mat"
 [5] "H_estimates"    "Hc_estimates"   "nb_boots"       "original_sg"
 [9] "outcome_type"   "effect_measure" "est.scale"      "mr_replicates"
```

Cross-checks:

```
gate naive$est (oriented)   = 103.1366726933
H0             (raw)        = -103.1366726933      <- exact negative of naive
H2             (raw)        = -81.1914643994
-H2 (oriented, = fb_H_est)  =  81.1914643994
batch bundle fb_H_est rep 1 =  81.1914644          <- reproduces exactly
```

Note `sdH1 = 1.32` against `sdH2 = 41.74`. `H1` and `H2` are nearly the same
point estimate (−80.87 vs −81.19) with SEs differing by a factor of 32. The
`H1` interval is implausibly tight for n = 500; the reported `fb_H_se` comes
from `H2` and is the defensible one, but the `H1` SE is worth a look
independently of this adjudication.

---

## V2 — per-replicate FB correction, `nv_H_est - fb_H_est`

Both terms oriented, so the correction is positive when FB pulls the naive
estimate toward the null. n = 50, all detected, no `fb_err`.

```
-- deciles --
       0%       10%       20%       30%       40%       50%       60%       70%       80%       90%      100%
10.185757 15.259890 17.143650 18.432914 18.623512 19.803153 20.388193 21.971040 23.138205 24.584733 25.677930

-- summary --
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  10.19   17.63   19.80   19.70   22.80   25.68

mean = 19.70016067   sd = 3.59717924
```

**Near-zero mass: none.**

```
  |correction| <  0.10 :  0 / 50
  |correction| <  0.50 :  0 / 50
  |correction| <  1.00 :  0 / 50
  |correction| <  2.00 :  0 / 50
  |correction| <  5.00 :  0 / 50
  correction <= 0     :  0 / 50
  min |correction|    : 10.18575722
```

Every replicate receives a substantial, strictly positive correction; the
smallest is 10.19. **FB is not silently no-op'ing.** It corrects hard and still
lands short of the target.

For contrast, the MR correction `nv_H_est - mr_H_est` on the same replicates:

```
       0%       10%       20%       30%       40%       50%       60%       70%       80%       90%      100%
30.575811 44.127360 46.384895 49.181222 50.560911 51.870837 53.445743 55.151018 56.990391 58.339464 61.815991
mean = 51.16077950   sd = 6.78982807

mean(FB correction)/mean(MR correction) = 0.385064
```

FB delivers **38.5%** of MR's correction. The two distributions do not overlap:
FB's maximum (25.68) is below MR's minimum (30.58). This is a systematic
difference in what is being corrected, not sampling noise.

---

## V3 — does the candidate family regenerate under resampling?

**Yes. Emphatically.** Replicate 1, `n = 500`, `get_FSdata()` called on the
original sample and on 3 bootstrap resamples with the same arguments
`forestsearch()` uses (`use_lasso = FALSE`, `use_grf = FALSE`,
`cont.cutoff = 4`, `.event_placeholder` per `forestsearch_main.R:2595-2601`).

Cut **definitions** compared, not `confs_names` — the latter are positional
labels (`q1`…`q26`) and comparing them would compare slot numbers, not cuts.

```
ORIGINAL sample (replicate 1, n=500): 26 candidate cut labels
 [1] "age <= 35"      "age <= 34"      "age <= 30"      "age <= 40"
 [5] "preanti <= 350" "preanti <= 132" "preanti <= 0"   "preanti <= 677"
 [9] "wtkg <= 76"     "wtkg <= 75"     "wtkg <= 68"     "wtkg <= 84"
[13] "cd40 <= 346"    "cd40 <= 335"    "cd40 <= 257"    "cd40 <= 415"
[17] "cd80 <= 1022"   "cd80 <= 948"    "cd80 <= 656"    "cd80 <= 1227"
[21] "hemo"           "homo"           "drugs"          "race"
[25] "gender"         "symptom"
```

| | n candidates | shared with original |
|---|---|---|
| original | 26 | — |
| resample 1 | 24 | 12 / 26 |
| resample 2 | 25 | 15 / 26 |
| resample 3 | 26 | 14 / 26 |

Across original + 3 resamples, 60 distinct cut definitions appear.
**10 are present in all four; 50 are not.** The 10 stable ones are exactly the
six categoricals (`hemo`, `homo`, `drugs`, `race`, `gender`, `symptom`) plus
four degenerate-or-coincident continuous cuts (`age <= 30`, `preanti <= 0`,
`wtkg <= 68`, `wtkg <= 76`).

Every quantile-derived continuous cut moves with the resample. `cd40 <= 415`
becomes `<= 420`, `<= 400`, `<= 414`; `cd80 <= 1022` becomes `<= 1009`,
`<= 994`, `<= 1050`; `preanti <= 350` becomes `<= 326`, `<= 330`, `<= 359`.

**The sharpest form of the result.** Replicate 1's selected subgroup was

```
!{cd40 <= 415} & !{cd80 <= 1022}
```

Both constituent cuts appear in the original family and in **none** of the three
resample families (`n_of_4 = 1` for each). The realized rule is not merely
re-estimated under resampling — it is **not expressible** in any of the three
resampled candidate families.

---

## What this means for the FB↔MR gap

Both simulation documents assert the opposite of V3. The batch qmd states:

> For the consistency engine the candidate family is fixed, so FB and MR
> target the same estimand and should agree to leading order

and the survival template draws the same line, treating family regeneration as
the thing that distinguishes `"dina"`/`"grf"` from `"consistency"`
(`sim_fs_maxeffCons_..._batch_1_50.qmd:535-542`). **That distinction does not
hold.** The family is fixed across the *consistency splits* within one fit, but
`forestsearch_bootstrap_dofuture()` re-runs the whole pipeline per resample
(`bootstrap_analysis_dofuture.R`), and `get_FSdata()` recomputes quantile cuts
from each resample's own empirical distribution. So `"consistency"` regenerates
its family under the bootstrap exactly as the other engines do.

Consequences, in the order they bite:

1. **FB and MR target different estimands.** MR conditions on the fixed fitted
   family; FB's `Hstar_obs`/`Hstar_star` terms are computed on subgroups drawn
   from families that differ from the original. The batch's FB↔MR gap is the
   *unconditional-vs-conditional* difference the survival template already
   names for `dina`/`grf` — it was simply thought not to apply here.

2. **The direction is consistent with V2.** FB's correction is diluted, not
   absent: each resample's re-selected subgroup is drawn from a *different*
   family, so `Hstar_star - Hstar_obs` mixes re-selection optimism with
   family-drift noise rather than isolating the optimism the original fit
   incurred. That is a plausible mechanism for a systematically undersized
   correction (38.5% of MR's), though V2 and V3 together establish the
   *association*, not the mechanism.

3. **"FB adjudicates MR" needs re-derivation before it is relied on.** The
   earlier batch summary called FB the adjudicator on the strength of the
   fixed-family premise. With that premise false, the batch's finding that
   "FB is worse than MR" is not evidence that MR is wrong; the two are not
   measuring the same thing.

Not done here, deliberately: no fix, no `R/` change, no re-run of the batch.
The natural next probes, if wanted — (a) pin the cut set across resamples via
`conf_force` and re-run FB on a few replicates to see whether the correction
moves toward MR's; (b) instrument `H_biasadj_2`'s four constituents per resample
to see which term carries the dilution; (c) check whether the `sdH1`/`sdH2`
factor-of-32 discrepancy noted in V1 is related or independent.
