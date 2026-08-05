---
bibliography: []
---

# Root cause: the ACTG175 naive-interval discrepancy is an upstream R change

**Status:** resolved. No forestsearch code change required.
**Companion:** `NOTE_actg175_naive_interval_discrepancy.md` records the evidence
that the discrepancy exists; this note explains it.
**Answer in one line:** R 4.6.0 fixed `stats::dfbeta()` for `glm` objects, and
the June reference was built under R 4.5.2, before that fix.

---

## Verdict: (a) a defect, not (b) a mislabelling

The quantity was **labelled correctly** and **computed wrongly**.

`sigma_D` is defined, and documented, as the robust subgroup standard error
`sqrt(sum(dfbeta^2))` built from the **one-step approximated dfbeta** of the
within-subgroup treatment effect. That is what the supplement describes and
what `R/consistency_resample.R` asks `stats::dfbeta()` for. In R <= 4.5.x,
`stats::dfbeta()` on a `glm` did not return one-step dfbeta values: it was
built on **deviance** residuals, which R's own NEWS says "do not give the
one-step approximated dfbeta values".

So the June payload's naive interval is not a different quantity wearing the
wrong name. It is the intended quantity, computed from a residual that was
the wrong input. Nothing was mislabelled; an arithmetic input was wrong, and
the error was in base R, not in forestsearch.

**Consequence for the supplement:** the *definition* of the naive interval
needs no correction. What is needed is a reproducibility caveat -- results
computed under R < 4.6.0 carry inflated `sigma_D` on GLM outcomes -- plus
regenerated ACTG175 numbers. The re-run in `dev/identifier-alignment/rerun/`
already supplies those, and they are correct.

The prime candidate considered beforehand -- that `naive_lo`/`naive_hi` had
been populated from the split-matched robust SE of the consistency screen
(Supplement S3.2) -- is **ruled out**. No such substitution occurred.

## The proof

R NEWS, **Changes in R 4.6.0**:

> `weighted.residuals()` now returns weighted working residuals for `"glm"`
> objects. It used to return deviance residuals, but these do not give the
> one-step approximated dfbeta values, as pointed out by Ravi Varadhan.
> `influence.glm()` has been amended so that it now returns an object with
> both a `"wt.res"` and a `"dev.res"` element.
>
> Several influence measures for `"glm"` objects have been updated to use
> Pearson residuals rather than deviance residuals (the latter are not
> unbiased) and to not use leave-one-out estimates of dispersion in models
> with fixed dispersion.

`lm.influence()` computes `dfbeta_i = (X'WX)^-1 x_i sqrt(w_i) * wt.res_i /
(1 - h_i)`. Substituting the two definitions of `wt.res` reproduces both
runs exactly:

| region | n | OLD (deviance residuals) | NEW (weighted working) | June reference | this re-run |
|---|---|---|---|---|---|
| FS H | 72 | **0.58019** | 0.51181 | 0.58019 | 0.51181 |
| DINA H | 98 | **0.49336** | 0.42753 | 0.49336 | 0.42753 |
| GRF H | 95 | **0.49643** | 0.42789 | 0.49643 | 0.42789 |
| FS H^c | 1011 | **0.15003** | 0.13098 | 0.15003 | 0.13098 |

Ratio of emulated-OLD to the June payload: 1.00001, 1.00000, 1.00000,
0.99999. The `NEW` column equals `stats::dfbeta()` in R 4.6.1 exactly, and
equals what this task's re-run produced.

This also explains the structural fact that defeated every constant-factor
explanation: the deviance-to-working residual ratio depends on the fitted
probabilities within each subgroup, so the inflation varies by subgroup
(1.134 to 1.160 here) rather than scaling uniformly.

### Reproducing the emulation

```r
sig <- function(dat) {
  g  <- glm(y ~ treat, data = dat, family = binomial)
  X  <- model.matrix(g); w <- g$weights; h <- stats::hatvalues(g)
  XtWXinv <- solve(crossprod(X * sqrt(w)))
  db <- function(wtres) (XtWXinv %*% t(X * (sqrt(w) * wtres / (1 - h))))["treat", ]
  c(NEW = sqrt(sum(db(sqrt(w) * residuals(g, type = "working"))^2)),   # R >= 4.6.0
    OLD = sqrt(sum(db(residuals(g, type = "deviance"))^2)))            # R <= 4.5.x
}
```

## What was ruled out along the way

Checked statically before touching any history:

- **`mr_draws`** -- `.fs_mr_assemble()` takes no draw count; `draws` reaches
  only the bias and IJ terms, never `se_wald`.
- **Every candidate SE on the observed data** -- model-based, HC0-HC4,
  `sqrt(sum(dfbetas^2))`, jackknife, `sqrt(2) * sigma_D`, and three
  leverage-reweighted variants. None matched; closest ratios were still
  1.10-1.19 and non-constant. (Incidentally, `sqrt(sum(dfbeta^2))` equals
  HC3 exactly for this model class in R 4.6.1.)
- **`se_ij`** -- the infinitesimal-jackknife SE, from a live refit: 0.530817
  for FS H, not 0.58019. The complement's `se_ij` (0.260351) is nowhere near
  its target (0.15003) either.
- **The naive CI formula** -- byte-identical between the June-era commit
  (`e89fb46`) and HEAD: `to_eff(beta_naive +/- z975 * se_wald)`.
- **`.consistency_glm_pieces()`** -- byte-identical June to now
  (`git diff e89fb46 HEAD -- R/consistency_resample.R` is empty).
- **Phase A / Phase B changes** -- Phase A altered no estimation code; the
  front-end change touches the FS row only, yet DINA and GRF show the same
  inflation.

That the package code is unchanged *and* reproduces the new value under this
R is what localised the cause outside the package.

## Scope: which outcomes are affected

The affected call is `stats::dfbeta()` on a `stats::glm()` fit inside
`.consistency_glm_pieces()`. Tested directly:

| `effect_measure` | fit used | affected? | ratio (old/new, example fit) |
|---|---|---|---|
| **OR** (binomial logit) | `glm` | **yes** | 1.165 |
| **RR** (binomial log / Poisson) | `glm` | **yes** | 1.164 |
| **RD** (binomial identity) | `glm` | **yes** | 1.137 |
| **IRR** (Poisson log, offset) | `glm` | **yes** | 1.054 |
| **MD** (continuous) | `stats::lm` | **no** | 1.000 |
| **HR** (survival) | `survival::coxph` | **no** | -- |

**Answer to the question asked: no, the continuous/MD outcome is not
affected.** `.consistency_glm_pieces()` dispatches `MD` to `stats::lm()`, not
`glm()`, and for an `lm` the working, Pearson and deviance residuals coincide
and the weights are 1, so the R 4.6.0 change is a no-op. A `gaussian()` glm
with identity link is likewise unaffected (verified: ratio exactly 1.00000).

So for `fs-glms-interpretable`, which uses both: **binary/OR results computed
under R < 4.6.0 need regenerating; continuous/MD results do not.** Count/rate
(IRR) and the binary RR/RD measures do need regenerating if any were produced
under the old R.

Survival is unaffected on a separate footing: `.consistency_cox_pieces()`
takes `residuals(fit, type = "dfbeta")` from `survival::coxph`, an
independent implementation. This is corroborated empirically -- the GBSG
application reproduced its June reference to 1.8e-15 on the same machine.

## What to do

1. **No forestsearch code change.** The package asks base R for the right
   quantity; base R now returns it.
2. **Regenerate GLM-outcome numbers produced under R < 4.6.0.** For ACTG175
   this is already done and lives in `dev/identifier-alignment/rerun/`.
3. **Add a reproducibility caveat** to the supplement: `sigma_D` on GLM
   outcomes is R-version dependent below 4.6.0, and the naive interval and
   everything downstream of `t_g` inherit it.
4. **Consider recording the R version** in the payload `meta`. The GBSG
   payload already carries `meta$machine$r_version`; the simulation bundles
   record `r_version` too. The ACTG175 payload's `meta` does not, which is
   why the version gap was not visible at comparison time.

## Downstream reach within a run

`sigma_D` is not confined to the naive interval. The screening threshold is
`t_g <- pmax(c_screen, c_consistency + z * sdv)`, so an inflated `sigma_D`
raises the bar every candidate must clear on each multiplier draw. In the
ACTG175 payload this moved the estimated selection bias materially -- FS
+0.899 to +0.692, DINA +0.708 to +0.592, GRF +0.451 to +0.334 -- and with it
every `mr_*` column. Any GLM-outcome analysis regenerated under R >= 4.6.0
should expect its MR block to move, not just its naive interval.
