# Predicted operating characteristics of the search over a candidate family

Monte-Carlo prediction of the maxeffCons search's operating
characteristics – family declaration rate, per-candidate declaration,
the selection distribution, expected selected-subgroup size,
sensitivity, specificity, PPV, NPV, the oriented effect on the selected
rule and the naive selection bias – from the joint normal law of the
candidates' fitted effects. The family is a population enumeration from
[`fs_oc_family_enumerate`](https://larry-leon.github.io/forestsearch/reference/fs_oc_family_enumerate.md)
or a supplied `fs_oc_family` object.

## Usage

``` r
fs_oc_predict(
  dgm = NULL,
  forestsearch_args = list(),
  n,
  c1 = NULL,
  c2 = NULL,
  family = NULL,
  consistency_method = c("resample", "split"),
  pconsistency = NULL,
  draws = 2e+05,
  seed = NULL,
  ...
)
```

## Arguments

- dgm:

  An object of class `"glm_dgm"`; used only when `family = NULL`.

- forestsearch_args:

  Named list of
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  arguments. Supplies the default `c1` (`effect.threshold`) and `c2`
  (`consistency.threshold`) and, when `family = NULL`, the cut
  specification and floors for
  [`fs_oc_family_enumerate`](https://larry-leon.github.io/forestsearch/reference/fs_oc_family_enumerate.md).

- n:

  Integer. Trial size. Overrides any size implied by the arguments; sets
  the size floor and the standard-error scale of an enumerated family
  and converts expected prevalence to expected subjects.

- c1:

  Numeric. Screening floor on the full-sample effect – the `c1` of
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  `effect.threshold` (legacy `hr.threshold`). Default
  `forestsearch_args$effect.threshold`; an explicit value overrides.

- c2:

  Numeric. Consistency floor on each half-sample effect – the `c2` of
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  `consistency.threshold` (legacy `hr.consistency`). An **effect**
  threshold, not a proportion; the proportion is `pconsistency` below.
  Default `forestsearch_args$consistency.threshold`; an explicit value
  overrides.

- family:

  `NULL` (enumerate with
  [`fs_oc_family_enumerate`](https://larry-leon.github.io/forestsearch/reference/fs_oc_family_enumerate.md))
  or an `fs_oc_family` object used as-is.

- consistency_method:

  `"resample"` (the package's production closed-form screen, the
  default) or `"split"` (the analytic document's single-split gate). See
  Details.

- pconsistency:

  Numeric in (0, 1). Consistency-*rate* threshold for the `"resample"`
  gate – the `p*` of
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  `pconsistency.threshold`: the fraction of splits that must clear `c2`,
  not an effect threshold. Default
  `forestsearch_args$pconsistency.threshold`, itself defaulting to
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  default (0.90); an explicit value overrides. Ignored by `"split"`.

- draws:

  Integer. Number of Monte-Carlo draws.

- seed:

  Integer or `NULL`. Passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the draws
  when non-`NULL`.

- ...:

  Passed to
  [`fs_oc_family_enumerate`](https://larry-leon.github.io/forestsearch/reference/fs_oc_family_enumerate.md)
  when `family = NULL` (`max_M`, `verbose`).

## Value

An object of class `c("fs_oc_predict", "list")`:

- `det_rate`, `det_rate_se`:

  Family declaration rate and its MC standard error.

- `P1`, `P1_se`:

  Per-candidate declaration probability.

- `p_sel`, `p_sel_se`:

  Per-candidate selection probability (unconditional).

- `sel_c`:

  Selection distribution given declaration.

- `EnH`:

  Expected selected-subgroup size in subjects, `n * sum(sel_c * Pg)`.

- `Esens`, `Espec`, `Eppv`, `Enpv`:

  Expected classification metrics of the selected rule against Q.

- `EbetaH`:

  Expected oriented true effect on the selected rule.

- `Enaive_bias`:

  Expected naive minus true effect on the selected rule.

- `mass_below`:

  Selection mass on rules whose true mean is below `c1`.

- `M`, `lab`:

  Family size and labels.

- `settings`:

  `n`, `c1`, `c2`, `consistency_method`, `pconsistency` (`NA` for
  `"split"`), `draws`.

- `seed`:

  The seed used.

- `family`:

  The family object.

## Details

The full-sample fitted effects of the M candidates are taken as
\\N(\beta_g, S_g)\\ with \\S_g = \rho \circ (se_g se_g')\\, \\\rho\_{ij}
= P(g_i \cap g_j) / \sqrt{P(g_i) P(g_j)}\\. Two independent half-sample
draws `W1`, `W2`, each with covariance `2 * Sg`, are generated through
[`fs_sym_root`](https://larry-leon.github.io/forestsearch/reference/fs_sym_root.md);
the full-sample effect is `Bhat = (W1 + W2) / 2`.

**The two gates.** Family declaration is any candidate declaring; the
selected rule is the effect maximiser among declaring candidates
(maxeffCons).

- `"resample"` – the package's production screen:

  The consistency stage of
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  (default `consistency_method = "resample"`) represents the random
  50/50 split's half effects as \\\hat\beta \pm D\\ and computes the
  consistency rate in closed form as \\2\Phi((\hat\beta -
  c_2)/\sigma_D) - 1\\, a candidate passing when that rate is at least
  `pconsistency` (`R/consistency_resample.R`;
  `R/subgroup_consistency_helpers.R` drops a candidate when
  `p.consistency < pconsistency.threshold`). Inverting: \\\hat\beta \ge
  c_2 + z_p \sigma_D\\ with \\z_p = \Phi^{-1}((1+p)/2)\\. \\\sigma_D\\
  is `sqrt(sum(dfbeta[i, treat]^2))`, the sandwich standard error of the
  subgroup treatment coefficient; in the draw model above the analogous
  quantity is `se_g`, and on a simulated MD40 trial the two agree within
  a few percent with no prevalence trend, so the wrapper identifies
  \\\sigma_D = se_g\\. The gate is then a single threshold on the
  full-sample draw, `(Bhat >= c1) & (Bhat - c2 >= z_p * se_g)`, and only
  `Bhat ~ N(beta_g, Sg)` is drawn – one matrix, not two.

- `"split"` – the analytic document's historical gate:

  One half-sample pair `W1`, `W2` stands in for the consistency rate:
  `(W1 + W2 >= 2 * c1) & (W1 >= c2) & (W2 >= c2)`, the screening floor
  `c1` on the full-sample effect and the consistency floor `c2` on each
  half. Kept bit-identical to the document's `worked-predictions` chunk.

All expectations are selection-weighted population functionals of the
family, conditional on declaration. Monte-Carlo standard errors are
given for the proportions (`sqrt(p (1 - p) / draws)`).

## See also

[`fs_oc_family_enumerate`](https://larry-leon.github.io/forestsearch/reference/fs_oc_family_enumerate.md),
[`fs_sym_root`](https://larry-leon.github.io/forestsearch/reference/fs_sym_root.md)

## Examples

``` r
# A hand-built three-candidate family
piQ <- 0.34
fam <- structure(list(
  lab = c("Q", "P", "D"), Pg = c(piQ, 0.45, 0.31),
  PQg = c(1, 0.28 / 0.45, 1), sens_g = c(1, 0.28 / piQ, 0.31 / piQ),
  spec_g = c(1, 1 - 0.17 / (1 - piQ), 1),
  ovl = matrix(c(piQ, 0.28, 0.31, 0.28, 0.45, 0.28 * 0.31 / piQ,
                 0.31, 0.28 * 0.31 / piQ, 0.31), 3, 3),
  M = 3L, PQ = piQ), class = c("fs_oc_family", "list"))
fam$beta_g <- 26 + 14 * fam$PQg
fam$se_g   <- 13.7 * sqrt(2) * sqrt(piQ / fam$Pg)
fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
              consistency_method = "resample", draws = 2e4, seed = 1)
#> Predicted operating characteristics (fs_oc_predict)
#>   family      : M = 3;  n = 500;  gate = resample (c1 = 30, c2 = 10, pcons = 0.9);  draws = 20,000;  seed = 1
#>   detection rate (DETECTED share)            : 0.589  (MC SE 0.0035)
#>   E[|Hhat|] given detection                  : 179.5
#>   E[sens] / E[spec] / E[PPV] / E[NPV]        : 0.915 / 0.927 / 0.893 / 0.952
#>   E[beta(Hhat)] oriented                     : 38.51
#>   E[naive - beta(Hhat)]                      : 18.22
#>   selection mass on rules with mean < c1     : 0.000
#>   top selections (P_selected, sel | det):
#>     D                                        0.236  0.401
#>     Q                                        0.187  0.317
#>     P                                        0.166  0.283
fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
              consistency_method = "split", draws = 2e4, seed = 1)
#> Predicted operating characteristics (fs_oc_predict)
#>   family      : M = 3;  n = 500;  gate = split (c1 = 30, c2 = 10);  draws = 20,000;  seed = 1
#>   detection rate (DETECTED share)            : 0.751  (MC SE 0.0031)
#>   E[|Hhat|] given detection                  : 178.8
#>   E[sens] / E[spec] / E[PPV] / E[NPV]        : 0.917 / 0.930 / 0.898 / 0.953
#>   E[beta(Hhat)] oriented                     : 38.57
#>   E[naive - beta(Hhat)]                      : 13.71
#>   selection mass on rules with mean < c1     : 0.000
#>   top selections (P_selected, sel | det):
#>     D                                        0.304  0.405
#>     Q                                        0.244  0.325
#>     P                                        0.203  0.270
```
