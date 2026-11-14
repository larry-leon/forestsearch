# Select a case-study trial from an MRCT simulation

Reads a
[`mrct_region_sims`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
result run with `record = "full"` and returns the trials that meet the
conditions of a regional-consistency case study: the **trigger** (the
upper 95\\ limit of the all-patient hazard ratio below `itt_ub` and the
region's hazard ratio above `region_itt`), a **declared** subgroup
\\\hat H\\, and a Cox hazard ratio below 1.0 in **region(\\\hat G\\)**,
the region without \\\hat H\\.

## Usage

``` r
mrct_select_case(
  res,
  conditions = c("trigger", "declared", "region_Ghat_lt1"),
  which = c("first", "all"),
  itt_ub = 1,
  region_itt = 1
)
```

## Arguments

- res:

  A
  [`mrct_region_sims`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
  result with `record = "full"` (it carries the `"estimators"`
  attribute; take the attribute's owner, not a subset of it).

- conditions:

  Character. The conditions a trial must meet, any of `"trigger"`,
  `"declared"`, `"region_Ghat_lt1"`. Default: all three.

- which:

  `"first"` (default): the first trial by index meeting the conditions;
  `"all"`: every such trial.

- itt_ub:

  Numeric. The trigger's bound on the all-patient upper 95\\ limit
  (strictly below); `Inf` drops the bound. Default 1.0.

- region_itt:

  Numeric. The trigger's bound on the region's hazard ratio (strictly
  above); `-Inf` drops the bound. Default 1.0.

## Value

A data.frame with one row per selected trial (zero rows, with a message,
when none qualifies): `sim`; the flags `trigger`, `declared`,
`region_Ghat_lt1` (logical); and the trial's hazard ratios `hr_itt_all`,
`itt_ub` (its upper 95\\ `hr_itt_region`, `hr_region_Ghat` and, when the
record has the `region_Gstar` row, `hr_region_Gstar`; and `sg_found`.
The attribute `"n_meeting"` is the number of trials meeting the
conditions.

## Details

The flags are evaluated per trial from the `"estimators"` attribute: the
trigger from the `itt_all` row's `hi` and the `itt_region` row's `hr`;
`declared` from the trial's `any_found == 1`; region(\\\hat G\\) \< 1.0
from the `region_Ghat` row's `hr`. A missing hazard ratio never
satisfies a condition. With `which = "first"` the result is the first
trial by index meeting every condition in `conditions`: with the three
default conditions, the rule that chose the case example of the
mrct-consistency report (replicate 5 of its S2 cell).

## See also

[`mrct_case_study`](https://larry-leon.github.io/forestsearch/reference/mrct_case_study.md)
to redraw and examine the selected trial;
[`mrct_region_sims`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md).

## Examples

``` r
if (FALSE) { # \dontrun{
res <- mrct_region_sims(dgm, n_sims = 20, n_sample = 1000,
                        region_var = "z_region", seed = 613800000,
                        record = "full")
mrct_select_case(res)                    # the first qualifying trial
mrct_select_case(res, which = "all")     # every qualifying trial
} # }
```
