---
title: "Oracle cell suppression and $\\kappa$ — S2 pilot"
bibliography: []
---

# Spec: suppress the oracle's conditional-estimand cells, report $\kappa$

**Supersedes `_implement_betaH_column.md` and `_expected_table_check.md`. Delete both —
they contain a non-existent function name, a target value off by 8%, a column since
dropped, and verification checks that fail on a correct implementation.**

File modified: `fs_t1_t2_m1_h10_knoise0_n500_combine_1_500.qmd` — **that file only**.
Not `betaHhat_truth.R`, not any `*_batch_*.qmd`, not anything under the package `R/`.

---

## Why

$\beta(\hat H)$ is the population effect at the **realized** region. The oracle refits on
the **true** region $H$, so it is not an estimator of $\beta(\hat H)$ and its
$b^{\beta(\hat H)}$ / $\hat C^{\beta(\hat H)}$ cells are neither a bias nor a coverage rate —
they measure the mismatch between the two regions. They are suppressed and that mismatch is
reported directly instead, as

$$\kappa \;=\; \frac{\beta(H)}{\beta(\hat H)},$$

the **identification factor**: how far the effect at the region the search returned sits
from the effect at the true harm region. Measured values for this cell:
$\beta(H) = 1.002890$, $\beta(H^c) = 0.628085$, $\kappa \approx 1.40$ (median $1.42$).

$\kappa$ is worth reporting on its own merits. The paper quantifies identification quality
only on the membership scale (sensitivity $0.39$, PPV $0.31$); $\kappa$ gives it on the
effect scale, which is what bears on inference. §S6.1 already states this qualitatively —
"the region forest search labels harmful is, in the population, mildly protective" —
and $\kappa = 1.40$ is that sentence with a number attached.

---

## Edits

**E1 — compute $\beta(H)$, $\beta(H^c)$.** Combine mode skips the DGM build (line 287), so
add a small chunk after the bundle pool. All constants are already in scope in this file.
About 5 s.

```r
k_inter <- calibrate_k_inter(target_hr_harm = target_hr_harm, model = dgm_model,
                             use_ahr = FALSE)
dgm     <- setup_gbsg_dgm(model = dgm_model, k_inter = k_inter,
                          n_super = n_super, seed = seed_base)
dgm     <- compute_dgm_cde(dgm)
eval_df <- build_eval_frame(dgm, analysis_time = analysis_time,
                            cens_adjust = cens_adjust, eval_seed = 20260628L)

# exp{beta(.)} at the TRUE regions. Same functional as beta(Hhat) -- only the region
# differs -- so the two are comparable on this one evaluation frame.
.bH <- betaHhat_theta_dagger_check(eval_df, harm.name = "flag_harm",
                                   outcome.name = "y_sim", event.name = "event_sim",
                                   treat.name = "treat_sim")
betaH_H  <- unname(.bH[["thetaDagger_H"]])    # name is legacy; the value is exp{beta(H)}
betaH_Hc <- unname(.bH[["thetaDagger_Hc"]])
```

**E2 — $\kappa$.** Insert near the `det_rate` / `cls` diagnostics (lines 930–935).

```r
# Identification factor: effect at the true region relative to the realized one.
# kappa = 1 exactly when the realized region reproduces the true region's effect.
.kap  <- betaH_H / rH$betaHhat_H
.kap  <- .kap[is.finite(.kap)]
kappa_txt <- if (length(.kap))
  sprintf("%.2f (median %.2f)", mean(.kap), stats::median(.kap)) else "\u2014"
```

**E3 — suppress the oracle cells.** In `build_block()` after `blk <- do.call(rbind, rows)`
(line 910), and identically in `build_block_med()` (line ~1028). Applies to **both** the
harm and complement blocks.

```r
  # The oracle refits on the TRUE region H, so it does not estimate beta(Hhat).
  # These two cells are a region mismatch, not a bias or a coverage rate.
  is_or <- blk$Estimator == "oracle"
  blk$b_betaHhat[is_or]   <- NA_real_
  blk$Cov_betaHhat[is_or] <- NA_real_
```

**E4 — render and footnote.** In **both** gt chains — main (937–974) and median (1040–1070):

```r
  sub_missing(columns = c(b_betaHhat, Cov_betaHhat), missing_text = "\u2014") |>
```

and extend the targets footnote:

```r
sprintf(paste("the oracle refits on the true region $H$, so it does not estimate",
              "$\\beta(\\widehat H)$; its $b^{\\beta(\\widehat H)}$ and",
              "$\\hat C^{\\beta(\\widehat H)}$ cells are a region mismatch rather than a",
              "bias or coverage rate and are suppressed. That mismatch is the",
              "identification factor $\\kappa = \\beta(H)/\\beta(\\widehat H)$ = %s, with",
              "$\\beta(H)$ = %.3f."),
        kappa_txt, betaH_H)
```

Check what `cols_hide` at line 958 currently hides before touching the chain.

---

## Verification

| # | Check | Expected |
|---|---|---|
| 1 | $\beta(H)$, $\beta(H^c)$ from E1 | $1.002890$, $0.628085$ |
| 2 | $\kappa$ | mean $\approx 1.40$, median $\approx 1.42$ |
| 3 | Oracle $b^{\beta(\hat H)}$, $\hat C^{\beta(\hat H)}$ | render as "—" in **both** blocks, **both** tables |
| 4 | Oracle $b^{\rm or}$ | still exactly $+0.0$ |
| 5 | naive / FB / MR $b^{\beta(\hat H)}$ | unchanged: $+128.0$, $+71.0$, $+40.0$ |
| 6 | Every other column | unchanged — Avg, SDs, Length, $b^{\rm or}$, $b^\ddagger$, $b^\dagger$ and all coverages |

Check 5 is the one that catches a mis-scoped edit: those three rows must not move.

---

## Manuscript notes (no code)

1. **§S6.1** gains a sentence: the oracle's suppressed cells, and $\kappa = 1.40$ as the
   effect-scale counterpart to the classification metrics.
2. **The SD ratio is free corroboration.** MR's $\widehat{\rm SD}_\beta/{\rm SD}_\beta =
   0.36/0.21 = 1.71$ against the oracle's calibrated $0.33/0.32 = 1.03$. Theorem 2 predicts
   $\hat V \to 4\sigma^2_{D,\hat H}$ under separation, i.e. standard errors inflated by up to
   $2\times$; $1.71$ sits between the local and separated limits exactly as the theory says.
   Both columns are already printed. One sentence turns a table entry into evidence for the
   theorem.
3. **Back-pocket, not for the paper.** The oracle's $b^\dagger = +13.5\%$ decomposes as
   non-collapsibility $1.0024$ $\times$ Jensen $1.0549$ $\times$ detection-conditioning
   $1.0730$. Non-collapsibility contributes essentially nothing at the borderline null. The
   manuscript makes no claim about this number, so nothing needs revising — but it is now
   explicable if a referee asks.

## Optional later

$\kappa$ for the other three DGM cells is four more eval-frame builds. Reported against PPV
across cells, it separates how much of the $b_\theta$ trajectory
($+40 \to +8 \to -1$) is identification improving from how much is the correction improving.
Cheap, but a separate decision from this pilot.
