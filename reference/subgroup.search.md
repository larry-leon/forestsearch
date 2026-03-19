# Subgroup Search for Treatment Effect Heterogeneity (Improved, Parallelized)

Searches for subgroups with treatment effect heterogeneity using
combinations of candidate factors. Evaluates subgroups for minimum
prevalence, event counts, and hazard ratio threshold. Parallelizes the
main search loop.

## Usage

``` r
subgroup.search(
  Y,
  Event,
  Treat,
  ID = NULL,
  Z,
  n.min = 30,
  d0.min = 15,
  d1.min = 15,
  hr.threshold = 1,
  max.minutes = 30,
  minp = 0.05,
  rmin = 5,
  details = FALSE,
  maxk = 2,
  parallel_workers = parallel::detectCores()
)
```

## Arguments

- Y:

  Numeric vector of outcome (e.g., time-to-event).

- Event:

  Numeric vector of event indicators (0/1).

- Treat:

  Numeric vector of treatment group indicators (0/1).

- ID:

  Optional vector of subject IDs.

- Z:

  Matrix or data frame of candidate subgroup factors (binary
  indicators).

- n.min:

  Integer. Minimum subgroup size.

- d0.min:

  Integer. Minimum number of events in control.

- d1.min:

  Integer. Minimum number of events in treatment.

- hr.threshold:

  Numeric. Hazard ratio threshold for subgroup selection.

- max.minutes:

  Numeric. Maximum minutes for search.

- minp:

  Numeric. Minimum prevalence rate for each factor.

- rmin:

  Integer. Minimum required reduction in sample size when adding a
  factor.

- details:

  Logical. Print details during execution.

- maxk:

  Integer. Maximum number of factors in a subgroup.

- parallel_workers:

  Integer. Number of parallel workers (default: all available cores).

## Value

List with found subgroups, maximum HR, search time, configuration info,
and filtering statistics.
