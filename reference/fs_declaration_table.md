# Declaration at a protected level, as a table

Reads the two
[`fs_declaration_calibration()`](https://larry-leon.github.io/forestsearch/reference/fs_declaration_calibration.md)
results a report needs from one fit (at the protected levels `c0`, and
at the claim threshold) with the fit's declaration field, into one tidy
table. Nothing is redrawn: the calibration reads the multiplier draws
the fit already made.

## Usage

``` r
fs_declaration_table(fit, alpha = 0.1, c0 = NULL, c0_primary = NULL)
```

## Arguments

- fit:

  A
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  fit (or a bare
  [`fs_mr_inference()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
  result) with the declaration field retained:
  `mr_inference_args = list(keep_declaration_field = TRUE, declaration_c0 = ...)`.
  `declaration_c0` stores the shifted maxima for those `c0` levels;
  `keep_field_matrix = TRUE` lets them be computed for any.

- alpha:

  Target family-wise declaration rate. Default 0.10.

- c0:

  `NULL` (the default, the claim-threshold row only) or a numeric vector
  of protected levels, on the natural scale of the consistency threshold
  (see
  [`fs_declaration_calibration()`](https://larry-leon.github.io/forestsearch/reference/fs_declaration_calibration.md)).

- c0_primary:

  Optional. The level, among `c0`, to mark as primary.

## Value

A plain data frame with one row per `c0`, then the claim-threshold row:

- `row`: `"c0"` or `"claim_threshold"`.

- `c0`: the level (`NA` on the claim-threshold row).

- `primary`: logical: the row is `c0_primary`.

- `alpha`, as given.

- `n_region`: the declared subgroup's size.

- `T_region`, `T_max`: the statistic of the declared subgroup and the
  largest in the family.

- `kappa_hat`: the calibrated cutoff.

- `implied_level`: \\2\Phi(\hat\kappa) - 1\\.

- `fw_size`: see above.

- `n_admitted`: the candidates the calibrated rule admits at this row.

- `admits_none`: logical: the calibrated rule admits nothing at any row.

- `source`: where the shifted maxima came from, `"capture"` or
  `"field_matrix"` (`NA` on the claim-threshold row).

## Details

`fw_size` is the family-wise size of the screen as executed, were every
candidate to carry a uniform effect of `c0`; at the claim threshold it
is a boundary size.

## See also

[`fs_declaration_calibration()`](https://larry-leon.github.io/forestsearch/reference/fs_declaration_calibration.md).

## Examples

``` r
set.seed(1)
n <- 500
d <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
d$time <- rexp(n, 0.1 * ifelse(d$treat == 1, ifelse(d$x1 < -0.6, 3, 0.45), 1))
d$event <- as.integer(d$time < 30)
d$time <- pmin(d$time, 30)
fit <- forestsearch(
  df.analysis = d, outcome.name = "time", event.name = "event", treat.name = "treat",
  id.name = "id", confounders.name = c("x1", "x2"), use_grf = FALSE, use_lasso = FALSE,
  maxk = 1, fs.splits = 200, sg_focus = "hrMaxSG", mr_inference = TRUE,
  mr_inference_args = list(ci_method = "ij", draws = 500, keep_declaration_field = TRUE,
                           declaration_c0 = c(0.6, 0.7)),
  parallel_args = list(plan = "sequential"), quiet = TRUE, details = FALSE,
  plot.sg = FALSE, show_candidate_summary = FALSE)
fs_declaration_table(fit, c0 = c(0.6, 0.7), c0_primary = 0.7)
#>               row  c0 primary alpha n_region T_region    T_max    kappa_hat
#> 1              c0 0.6   FALSE   0.1      134 6.663819 6.663819 -0.807684942
#> 2              c0 0.7    TRUE   0.1      134 6.663819 6.663819 -0.006648818
#> 3 claim_threshold  NA   FALSE   0.1      134 6.663819 6.663819  2.178527462
#>   implied_level fw_size n_admitted admits_none  source
#> 1   -0.58072803   0.000          1       FALSE capture
#> 2   -0.00530495   0.000          1       FALSE capture
#> 3    0.97063321   0.272          1       FALSE    <NA>
```
