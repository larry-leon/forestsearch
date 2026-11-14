# Core operating characteristics of an MRCT simulation campaign

The rows of the regional-consistency core table, computed from one
[`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
result run with `record = "full"`. Everything is read from the result:
its per-trial columns and its `"estimators"` attribute.

## Usage

``` r
mrct_core_summary(res)
```

## Arguments

- res:

  A
  [`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
  result with `record = "full"`. Take the result itself: a `data.table`
  subset of it drops the `"estimators"` attribute.

## Value

A plain data frame with one row per quantity above, in that order:

- `quantity`, `label`: the row's key and its description.

- `estimate`: the value, in percent.

- `se`: its standard error, in percent, as the core table carries it:
  binomial for the proportions, `sd / sqrt(n)` over the trials for the
  Sens and Spec averages, `NA` for `trigger` and `Ghat_share`, which the
  core table gives without one.

- `count`: the trials counted (`NA` for an average).

- `denominator`: the trials the row is taken over.

- `display`: the value as printed, at one decimal (for `trigger`,
  `"K of N"`).

A row with no trial in its denominator has `estimate` `NA` and `display`
an en dash. Without a `flag_harm` column in the simulated trials, the
Sens, Spec and `Gstar_lt1` rows are `NA`.

## Details

The rows, with their denominators:

- `power`: P(ITT ub \< 1.0), of all trials.

- `trigger`: trials reaching the trigger (ITT ub \< 1.0 and region ITT
  \> 1.0), of all trials.

- `declared`: a subgroup declared, of the triggered trials.

- `sens_nonregion`, `spec_nonregion`, `sens_region`, `spec_region`: the
  sensitivity and specificity of \\\hat H\\ against \\H^\*\\, trial by
  trial, averaged over the declared triggered trials. A trial where one
  is undefined does not enter its average.

- `declared_Ghat_lt1`: declared and region(\\\hat G\\) \< 1.0, of the
  triggered trials.

- `Gstar_lt1`: region(\\G^\*\\) \< 1.0, of the triggered trials.

- `declared_Ghat_ge1`: declared and region(\\\hat G\\) not below 1.0, of
  the triggered trials.

- `Ghat_share`: the share of the region's patients kept in \\\hat G\\,
  trial by trial, averaged over the declared triggered trials.

## See also

[`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md),
[`mrct_select_case()`](https://larry-leon.github.io/forestsearch/reference/mrct_select_case.md).

## Examples

``` r
# \donttest{
set.seed(1)
n <- 400
d <- data.frame(age = rnorm(n, 50, 10), lab = rnorm(n, 300, 100),
                f1 = rbinom(n, 1, 0.4), treat = rbinom(n, 1, 0.5))
d$tte <- rexp(n, 0.05 * exp(-0.5 * d$treat))
d$event <- as.integer(d$tte < 40)
d$tte <- pmin(d$tte, 40)
dgm <- suppressWarnings(inject_mrct_structure(
  d, outcome_var = "tte", event_var = "event", treatment_var = "treat",
  continuous_vars = c("age", "lab"), factor_vars = "f1", x_pred = "lab",
  spline_spec = list(knot = 300, zeta = 800, log_hrs = rep(log(0.55), 3)),
  region = list(prevalence = 0.25, or_pred = 1, x1_vars = "f1", or_x1 = 3,
                loghr = -log(3), band = list(cut = 230, or = 1)),
  x3 = list(vars = "lab", cuts = list(lab = 230), loghr = log(3 / 0.55)),
  n_super = 2000, expand = "bootstrap", seed = 7L,
  continuous_vars_cens = c("age", "lab"), factor_vars_cens = "f1"))
dgm <- scale_dgm_hazard(dgm, m = 3)
res <- mrct_region_sims(
  dgm, n_sims = 8, n_sample = 800, region_var = "z_region", sg_focus = "hrMaxSG",
  confounders.name = analysis_covariates(dgm, "observed"), conf_force = character(0),
  fs_args = list(use_grf = FALSE, use_lasso = FALSE),
  sim_args = list(rand_ratio = 1, draw_treatment = TRUE, max_entry = 6, cens_adjust = 20),
  analysis_time = 30, cens_adjust = 20, parallel_args = list(plan = "sequential"),
  details = FALSE, verbose_n_sims = 0, seed = 5000L, record = "full")
mrct_core_summary(res)[, c("quantity", "estimate", "se", "denominator", "display")]
#>             quantity  estimate        se denominator display
#> 1              power  87.50000 11.692679           8    87.5
#> 2            trigger  25.00000        NA           8  2 of 8
#> 3           declared 100.00000  0.000000           2   100.0
#> 4     sens_nonregion  92.81385  2.337662           2    92.8
#> 5     spec_nonregion 100.00000  0.000000           2   100.0
#> 6        sens_region  96.63189  1.549925           2    96.6
#> 7        spec_region 100.00000  0.000000           2   100.0
#> 8  declared_Ghat_lt1 100.00000  0.000000           2   100.0
#> 9          Gstar_lt1 100.00000  0.000000           2   100.0
#> 10 declared_Ghat_ge1   0.00000  0.000000           2     0.0
#> 11        Ghat_share  71.55184        NA           2    71.6
# }
```
