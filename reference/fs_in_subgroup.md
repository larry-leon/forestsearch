# Membership of a forestsearch fit's declared subgroup

Applies the fit's declared rule to the rows of `data`, by
[`get_dfpred()`](https://larry-leon.github.io/forestsearch/reference/get_dfpred.md).

## Usage

``` r
fs_in_subgroup(fit, data)
```

## Arguments

- fit:

  A
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  fit.

- data:

  A data frame with the columns the rule names, such as the fit's whole
  trial.

## Value

A logical vector, one value per row: `TRUE` in the declared subgroup
\\\hat H\\. All `FALSE` when nothing was declared.

## See also

[`get_dfpred()`](https://larry-leon.github.io/forestsearch/reference/get_dfpred.md),
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).

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
  maxk = 1, fs.splits = 200, sg_focus = "hrMaxSG",
  parallel_args = list(plan = "sequential"), quiet = TRUE, details = FALSE,
  plot.sg = FALSE, show_candidate_summary = FALSE)
fit$sg.harm
#>           q3.1 
#> "{x1 <= -0.7}" 
table(fs_in_subgroup(fit, d))
#> 
#> FALSE  TRUE 
#>   366   134 
```
