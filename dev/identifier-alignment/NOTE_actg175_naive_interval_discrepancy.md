---
bibliography: []
---

# Note: the ACTG175 reference payload's naive intervals are not reproducible

**Status:** RESOLVED -- root cause found. See
`NOTE_actg175_naive_interval_rootcause.md`. In short: R 4.6.0 fixed
`stats::dfbeta()` for `glm` objects (it was built on deviance residuals,
which "do not give the one-step approximated dfbeta values"), and the June
reference was built under R 4.5.2. The forestsearch code is correct and
unchanged; the numbers in this note's tables are reproduced exactly by
emulating the pre-4.6.0 formula. No package change is required; GLM-outcome
results produced under R < 4.6.0 need regenerating. Continuous/MD and
survival outcomes are unaffected.

The evidence below is retained as the record of how the discrepancy was
characterised before the cause was known.
**Found:** 2026-08-05, during the identifier-alignment re-analysis (Phase B).
**Affects:** `dev/identifier-alignment/sim_analyses/actg175_table_payload.rds`
(built 2026-06-20, `forestsearch_version = "0.2.0"`), and any manuscript
table drawn from it.
**Does not affect:** `gbsg_table_payload.rds`, which reproduces exactly.

---

## Summary

Re-rendering `analysis_actg175_binary_multimethod.qmd` reproduces the
reference payload's subgroup **labels**, **n**, and naive **point estimates**
exactly, but every naive confidence interval is roughly 13% narrower in the
new run. The difference appears in **all six rows** -- FS, DINA and GRF, both
H and H^c -- so it is not a consequence of any configuration change made
during this task.

Computing the interval independently from the raw `speff2trial::ACTG175`
data shows the **new run is correct** and the reference is not: the new
interval is exactly `beta_naive +/- z_{.975} * sigma_D` with
`sigma_D = sqrt(sum(dfbeta^2))`, which is what
`R/fs_mr_inference.R:454-456` documents the naive interval to be. The
reference's implied standard errors are 12-16% larger and match no standard
estimator computable from the data.

## The evidence

Both runs' intervals are symmetric on the log scale, i.e. genuine Wald
intervals, so the implied standard error can be recovered as
`(log(hi) - log(lo)) / (2 * 1.959964)`.

| region | n | `sqrt(sum(dfbeta^2))` | reference implied SE | new implied SE |
|---|---|---|---|---|
| FS H | 72 | **0.51181** | 0.58019 | **0.51181** |
| DINA H | 98 | **0.42753** | 0.49336 | **0.42753** |
| GRF H | 95 | **0.42789** | 0.49643 | **0.42789** |
| FS H^c | 1011 | **0.13098** | 0.15003 | **0.13098** |

For comparison, on FS H (n = 72): model-based SE 0.49748, HC0 0.49748,
HC1 0.50453. None is 0.58019.

The point estimates agree to `1e-15` in all six rows, so the two runs select
the same subgroup and fit the same model to the same rows. Only the
dispersion differs.

## The inflation is not a constant factor

This is the key structural fact, and it rules out the simplest explanation.

| region | n | ref / new |
|---|---|---|
| FS H | 72 | 1.1336 |
| DINA H | 98 | 1.1540 |
| GRF H | 95 | 1.1602 |
| FS H^c | 1011 | 1.1455 |

A 2.3% spread, uncorrelated with `n`. A `sigma_D` miscomputed by a constant
multiplier (a missing finite-sample correction, a wrong denominator, a
mis-set `z`) would give the same ratio in every row. It does not. Nor is it
a degrees-of-freedom correction: matching FS H would need `k ~ 17` and
FS H^c `k ~ 240`.

The signature instead points to `naive_lo`/`naive_hi` having been populated
in June from a **different quantity** -- one that varies per subgroup in a
way `sqrt(sum(dfbeta^2))` does not.

## What has been ruled out

- **Not `mr_draws`.** Statically: `.fs_mr_assemble(df, candidates, spec)`
  takes no draw count, and `draws` enters only via `.fs_mr_multipliers()`
  -> `P` -> the bias and IJ terms. `se_wald <- sdv[sel]` cannot depend on it.
  (`mr_draws` was raised 2000 -> 5000 in this task to match the reference's
  own `meta$gate_draws`; that change is unrelated.)
- **Not the `use_lasso`/`use_grf` front-end change.** That touches the FS row
  only, and the DINA and GRF rows show the same inflation. Independently,
  the FS MR family size is 2794 in both runs -- the prefilter was inert on
  the full ACTG175 sample.
- **Not the Phase A package changes.** Phase A altered defaults,
  documentation and a validator; it changed no estimation code.
- **Not the platform.** GBSG (survival / Cox) reproduces its reference to
  `1.8e-15` on the same machine and package build, so the run environment is
  not at fault. The effect is confined to the GLM path.
- **Not a change to `sigma_D` itself.** `R/consistency_resample.R`, where
  `sigma_D <- sqrt(sum(dfbeta^2))` is computed, has no commits since before
  the reference was built. (This turned out to be the decisive clue: the
  package code being unchanged *and* reproducing the new value under this R
  is what localised the cause to `stats::dfbeta()` in base R.)

## Downstream consequence

The naive interval is not the only thing affected. The screening threshold is

```
t_g <- pmax(c_screen, c_consistency + z * sdv)
```

so an inflated `sdv` raises the bar every candidate must clear on each
multiplier draw, changing which candidate wins and hence the estimated
selection bias:

| engine | reference `selection_bias` | new |
|---|---|---|
| FS | +0.89932 | +0.69203 |
| DINA | +0.70764 | +0.59230 |
| GRF | +0.45058 | +0.33437 |

Every `mr_*` column in the payload moves accordingly. One root cause
explains all of the ACTG175 differences except the FS `fb_*` block, which
differs by +1.64% for the separate and expected reason that the lasso/GRF
prefilter is inert on the observed data but not inside bootstrap replicates.

## Reproducing this note's numbers

The independent calculation needs no package internals:

```r
library(speff2trial)
d <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
d$treat <- ifelse(d$arms == 1L, 1L, 0L)
d <- d[!is.na(d$cd420), ]
d$y <- 1L - as.integer(d$cd420 > d$cd40)

H <- (d$wtkg > 86) & (d$cd40 > 380)          # the FS subgroup, n = 72
f <- glm(y ~ treat, data = d[H, ], family = binomial)
sqrt(sum(stats::dfbeta(f)[, "treat"]^2))     # 0.51181 = the new run
```

The DINA region is `preanti >= 849.4 & cd40 >= 338` (n = 98); the GRF region
is `wtkg > 84.3696 & cd40 > 368.2` (n = 95).

## A second, independent defect also affects the DINA and GRF rows

Separately from the R 4.6.0 issue above, MR was applying a consistency floor
to DINA and GRF -- engines that have no consistency screen -- because its
admission set was rebuilt from raw parameters rather than resolved once. In
the configuration measured, DINA's multiplier draws passed that floor only
**5.1%** of the time, so its MR correction rested on a twentieth of the
resampling distribution.

The two compound: the offending term is `z * sigma_D`, so an inflated
`sigma_D` makes the floor bite harder. See
`NOTE_mr_admission_set.md`. The `mr_*` columns of the DINA and GRF rows should
be regenerated rather than trusted, on the same footing as the FS row.

## Bearing on the manuscript

Until the source of the June values is identified, the `naive_*` and `mr_*`
columns of `actg175_table_payload.rds` should not be quoted as-is. The
labels, the `n` values, the naive point estimates and the DINA/GRF `fb_*`
columns are unaffected and reproduce.

Whether this is a defect (a wrongly computed `sigma_D`) or a mislabelling
(a different, correct quantity reported under the naive heading) determines
which correction the supplement needs; see
`NOTE_actg175_naive_interval_rootcause.md` for that investigation.
