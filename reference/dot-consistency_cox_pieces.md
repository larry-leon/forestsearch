# Per-subject Cox pieces for the consistency approximation

Fits the treatment-only (or covariate-adjusted) Cox model once on a
subgroup and returns the treatment log hazard ratio, the per-subject
treatment `dfbeta` influence contributions, and the split-perturbation
SD.

## Usage

``` r
.consistency_cox_pieces(
  df,
  tte.name = "Y",
  event.name = "Event",
  treat.name = "Treat",
  cox_init = 0,
  adjust_covariates = NULL
)
```

## Arguments

- df:

  data.frame/data.table with the survival columns (and any columns
  referenced by `adjust_covariates`).

- tte.name, event.name, treat.name:

  Character; time, event, and (0/1) treatment column names (defaults
  `"Y"`, `"Event"`, `"Treat"`).

- cox_init:

  Numeric; warm-start for the unadjusted Cox coefficient.

- adjust_covariates:

  Character vector or `NULL`; terms appended after `treat.name`.

## Value

List with `beta_hat`, `dfbeta`, `sigma_D`, `n`, `d`, `log_scale` (TRUE),
`measure` (`"HR"`); or `NULL` if degenerate.
