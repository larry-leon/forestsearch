# Run Single Consistency Split

Performs one random 50/50 split and evaluates whether both halves meet
the HR consistency threshold.

## Usage

``` r
run_single_consistency_split(df.x, N.x, hr.consistency, cox_init = 0)
```

## Arguments

- df.x:

  data.table. Subgroup data with columns Y, Event, Treat.

- N.x:

  Integer. Number of observations in subgroup.

- hr.consistency:

  Numeric. Minimum HR threshold for consistency.

- cox_init:

  Numeric. Initial value for Cox model (log HR).

## Value

Numeric. 1 if both splits meet threshold, 0 if not, NA if error.
