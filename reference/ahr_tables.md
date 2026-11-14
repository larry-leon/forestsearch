# Tables of the MRCT analytic layer, over cells and sets

Tidy data frames from the existing analytic functions, without
quadrature of their own: `ahr_truths_table()` loops over
[`ahr_limit()`](https://larry-leon.github.io/forestsearch/reference/ahr_limit.md)
and `ahr_design_table()` over
[`ahr_design()`](https://larry-leon.github.io/forestsearch/reference/ahr_design.md).

## Usage

``` r
ahr_truths_table(
  x,
  sets = c("all", "nonregion", "region", "region_Gstar"),
  b = NULL,
  h = NULL,
  analysis_time = NULL,
  max_entry = NULL,
  ...
)

ahr_design_table(x, b, h, n, analysis_time = NULL, max_entry = NULL, ...)
```

## Arguments

- x:

  The cell or cells: an
  [`ahr_population()`](https://larry-leon.github.io/forestsearch/reference/ahr_population.md)
  result, a DGM from
  [`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md),
  or a named list of either. The list's names label the cells.

- sets:

  The sets, as
  [`ahr_limit()`](https://larry-leon.github.io/forestsearch/reference/ahr_limit.md)
  takes them: a character vector of set names (`"all"`, `"nonregion"`,
  `"region"`, `"region_Gstar"`, `"Hstar"`, `"Hstar_region"`,
  `"Hstar_nonregion"`), or a list of those and logical vectors of length
  `pop$N`. A name on an element labels the set, and a logical vector
  needs one.

- b, h:

  Hazard ratios outside and inside \\H^\*\\: one number, or a vector
  named by cell. In `ahr_truths_table()`, both `NULL` (the default) use
  the DGM's own per-patient log hazard ratios.

- analysis_time, max_entry:

  Passed to
  [`ahr_population()`](https://larry-leon.github.io/forestsearch/reference/ahr_population.md)
  when `x` holds DGMs, and not used for an
  [`ahr_population()`](https://larry-leon.github.io/forestsearch/reference/ahr_population.md)
  result, which carries its own.

- ...:

  Passed to
  [`ahr_limit()`](https://larry-leon.github.io/forestsearch/reference/ahr_limit.md)
  (for instance `horizon`, `b_R`, `h_R`) or
  [`ahr_design()`](https://larry-leon.github.io/forestsearch/reference/ahr_design.md)
  (`b_R`, `h_R`).

- n:

  The trial size.

## Value

`ahr_truths_table()`: a plain data frame with one row per cell and set:
`cell`, `set`, `horizon` (the analysis time), `marghr` (the marginal
hazard ratio at the horizon), `marghr_inf` (with no censoring), `ahr`,
`cde`.

`ahr_design_table()`: a plain data frame with one row per cell: `cell`,
`n`, `b`, `h`; the five probabilities of
[`ahr_design()`](https://larry-leon.github.io/forestsearch/reference/ahr_design.md)
under their names there (`p_itt_ub_lt1`, `p_region_itt_gt1`,
`p_region_Gstar_lt1`, `p_trigger`, `p_trigger_and_Gstar_lt1`);
`p_Gstar_lt1_given_trigger`, their ratio; and the seven rows above,
`r1_itt_ub_lt1`, `r2_region_itt_gt1`, `r3_consistency`,
`r4_non_consistency`, `r5_oracle_consistency`,
`r6_oracle_non_consistency`, `r7_consistency_or_oracle`.

## Details

`ahr_design_table()`'s `r1_` to `r7_` columns are the seven design
probabilities of a regional-consistency analysis. With \\A\\ = {ITT ub
\< 1.0}, \\B\\ = {region ITT \> 1.0} and \\C^\*\\ = {region(\\G^\*\\) \<
1.0}: (1) \\A\\, (2) \\B\\, (3) consistency \\A \setminus B\\, (4)
non-consistency (the trigger) \\A \cap B\\, (5) oracle consistency \\A
\cap B \cap C^\*\\, (6) oracle non-consistency, and (7) consistency or
oracle consistency. The function stops unless (3) + (4) = (1), (5) + (6)
= (4) and (3) + (5) = (7) at 1e-12.

## See also

[`ahr_population()`](https://larry-leon.github.io/forestsearch/reference/ahr_population.md),
[`ahr_limit()`](https://larry-leon.github.io/forestsearch/reference/ahr_limit.md),
[`ahr_design()`](https://larry-leon.github.io/forestsearch/reference/ahr_design.md).

## Examples

``` r
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
pop <- ahr_population(dgm, analysis_time = 30, max_entry = 6)
ahr_truths_table(pop, sets = c("all", "region", "region_Gstar"), b = 0.55, h = 3)
#>   cell          set horizon    marghr marghr_inf       ahr      cde
#> 1 cell          all      30 0.9027294  0.8078874 0.9079770 1.340451
#> 2 cell       region      30 1.1084587  0.7639976 0.9296125 1.388062
#> 3 cell region_Gstar      30 0.5503904  0.5518979 0.5500000 0.550000
ahr_design_table(pop, b = 0.55, h = 3, n = 800)
#>   cell   n    b h p_itt_ub_lt1 p_region_itt_gt1 p_region_Gstar_lt1  p_trigger
#> 1 cell 800 0.55 3    0.2003976         0.669416          0.9610761 0.08177303
#>   p_trigger_and_Gstar_lt1 p_Gstar_lt1_given_trigger r1_itt_ub_lt1
#> 1              0.07989874                 0.9770793     0.2003976
#>   r2_region_itt_gt1 r3_consistency r4_non_consistency r5_oracle_consistency
#> 1          0.669416      0.1186246         0.08177303            0.07989874
#>   r6_oracle_non_consistency r7_consistency_or_oracle
#> 1               0.001874295                0.1985233
```
