# Bias and coverage summary for one simulation cell

From a combined bundle's per-replicate `results` (the recorder columns
of the `fb_mr_field` template family: `nv_H_*`, `mr_H_*`, `fld_H_*`,
`betaHhat_H`, `detected`, and their `Hc` twins), computes, per estimator
over the detected replicates: the retained bias on the log scale, the
empirical SD of the log estimate, the mean reported SE, their ratios
`b = bias_log / sd_emp` and `r = se_mean / sd_emp`, the observed one-
and two-sided coverage of the target with Wilson intervals, and the
Gaussian reference coverages at (b, r).

## Usage

``` r
fs_sim_bias_coverage(
  results,
  block = c("H", "Hc"),
  estimators = c("naive", "mr", "fld"),
  level = 0.95,
  target = "betaHhat",
  side = c("lower", "upper"),
  scale = c("log", "identity")
)
```

## Arguments

- results:

  Data frame of per-replicate records from a combined bundle
  (`readRDS(...)$results`).

- block:

  `"H"` (harm, default) or `"Hc"` (complement). The `"fld"` estimator is
  dropped with a message when the block's field columns are absent
  (`fld_Hc_*` exist only for bundles run with
  `field_complement = TRUE`).

- estimators:

  Subset of `c("naive", "mr", "fld", "mr_w", "mr_wf")` (default the
  first three). `"mr_w"` / `"mr_wf"` are MR (IJ) with the winner-only /
  winner-floor SE (`mr_<block>_se_w` / `_se_wf` and their `_lo_w/_hi_w`
  / `_lo_wf/_hi_wf` bounds), available only on bundles that recorded
  them.

- level:

  Coverage level (default 0.95).

- target:

  `"betaHhat"` (default; the per-replicate conditional target
  `betaHhat_<block>`), `"oracle"` (the per-replicate oracle estimate
  `or_<block>_est`), or a numeric scalar on the HR scale (e.g. a CDE or
  marginal truth value from the bundle's `truth` element).

- side:

  `"lower"` (default; one-sided coverage of the lower bound, the
  harm-claim orientation) or `"upper"` (one-sided coverage of the upper
  bound, the benefit-claim orientation). Two-sided quantities are
  unaffected.

- scale:

  `"log"` (default) or `"identity"`. The working scale of the estimate
  columns: `"log"` for ratio measures (HR/OR/RR; bias, SD and the
  normal-based one-sided bound `exp(log est -/+ z * SE)` are computed on
  the log scale, and non-positive estimates are dropped), `"identity"`
  for difference measures (MD/RD; bias, SD and the bound
  `est -/+ z * SE` on the natural scale, no positivity guard). Add-only:
  the default reproduces the previous output exactly. Under `"identity"`
  the `bias_log` column holds the identity-scale bias (the column name
  is kept so stacked tables and
  [`fs_plot_bias_coverage()`](https://larry-leon.github.io/forestsearch/reference/fs_plot_bias_coverage.md)
  read either scale unchanged).

## Value

Data frame, one row per estimator: `estimator`, `n`, `bias_log`,
`sd_emp`, `se_mean`, `b`, `r`, `cov1`, `cov1_wilson_lo`,
`cov1_wilson_hi`, `cov2`, `cov2_wilson_lo`, `cov2_wilson_hi`,
`cov1_ref`, `cov2_ref`.

## Details

Conventions: the one-sided lower bound is the stored `fld_<block>_lo1s`
for the field row and `exp(log est - qnorm(level) * SE)` for naive and
MR (IJ); two-sided coverage uses each estimator's stored bounds; the
field estimator is the shrunk-field `est2` with the two-sided
Lambda-quantile interval. Under `side = "upper"` (the complement's
benefit-claim orientation) the one-sided coverage is
`truth <= upper bound`, with the field's stored `fld_<block>_up1s` and
`exp(log est + qnorm(level) * SE)` for naive and MR (IJ), and the
Gaussian reference is `Phi(z_level * r + b)` – a positive retained bias
helps an upper bound where it hurts a lower one.
