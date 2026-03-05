# Evaluate Subgroup Consistency

Evaluates candidate subgroups using split-sample consistency validation.
For each candidate, repeatedly splits the data and checks whether the
treatment effect direction is consistent across splits.

## Usage

``` r
subgroup.consistency(
  df,
  hr.subgroups,
  hr.threshold = 1,
  hr.consistency = 1,
  pconsistency.threshold = 0.9,
  m1.threshold = Inf,
  n.splits = 100,
  details = FALSE,
  by.risk = 12,
  plot.sg = FALSE,
  maxk = 7,
  Lsg,
  confs_labels,
  sg_focus = "hr",
  stop_Kgroups = 10,
  stop_threshold = NULL,
  showten_subgroups = FALSE,
  pconsistency.digits = 2,
  seed = 8316951,
  checking = FALSE,
  use_twostage = FALSE,
  twostage_args = list(),
  parallel_args = list()
)
```

## Arguments

- df:

  Data frame containing the analysis dataset. Must include columns for
  outcome (Y), event indicator (Event), and treatment (Treat).

- hr.subgroups:

  Data.table of candidate subgroups from subgroup search, containing
  columns: HR, n, E, K, d0, d1, m0, m1, grp, and factor indicators.

- hr.threshold:

  Numeric. Minimum hazard ratio threshold for candidates. Default: 1.0

- hr.consistency:

  Numeric. Minimum HR required in each split for consistency. Default:
  1.0

- pconsistency.threshold:

  Numeric. Minimum proportion of splits that must be consistent.
  Default: 0.9

- m1.threshold:

  Numeric. Maximum m1 threshold for filtering. Default: Inf

- n.splits:

  Integer. Number of splits for consistency evaluation. Default: 100

- details:

  Logical. Print progress details. Default: FALSE

- by.risk:

  Numeric. Risk interval for KM plots. Default: 12

- plot.sg:

  Logical. Generate subgroup plots. Default: FALSE

- maxk:

  Integer. Maximum number of factors in subgroup. Default: 7

- Lsg:

  List of subgroup parameters.

- confs_labels:

  Character vector mapping factor names to labels.

- sg_focus:

  Character. Subgroup selection criterion: "hr", "maxSG", or "minSG".
  Default: "hr"

- stop_Kgroups:

  Integer. Maximum number of candidates to evaluate. Default: 10

- stop_threshold:

  Numeric in `[0, 1]` or `NULL`. When a candidate subgroup's consistency
  probability (`Pcons`) meets or exceeds this threshold, evaluation
  stops early — remaining candidates are skipped. Set to `NULL` to
  disable early stopping and evaluate all candidates up to
  `stop_Kgroups`. Default: `NULL`.

  **Note:** Values \> 1.0 are not permitted. To disable early stopping,
  use `stop_threshold = NULL`, not a value above 1.

  **Interaction with `sg_focus`:**

  `"hr"`, `"maxSG"`, `"minSG"`

  :   Early stopping is valid because candidates are sorted by a single
      criterion. The first candidate passing the threshold is optimal
      under that criterion.

  `"hrMaxSG"`, `"hrMinSG"`

  :   Should generally be `NULL`, because these compound criteria
      require comparing HR *and* size across all candidates.
      [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
      automatically resets to `NULL` with a warning for these.

  For parallel execution, early stopping is checked after each batch
  completes, so some additional candidates beyond the first meeting the
  threshold may be evaluated. Use a smaller `batch_size` in
  `parallel_args` for finer-grained early stopping.

- showten_subgroups:

  Logical. If TRUE, prints up to 10 candidate subgroups after sorting by
  sg_focus, showing their rank, HR, sample size, events, and factor
  definitions. Useful for reviewing which candidates will be evaluated
  for consistency. Default: FALSE

- pconsistency.digits:

  Integer. Decimal places for consistency proportion. Default: 2

- seed:

  Integer. Random seed for reproducible consistency splits.
  Default: 8316951. Set to NULL for non-reproducible random splits. The
  seed is used both for sequential execution (via set.seed()) and
  parallel execution (via future.seed).

- checking:

  Logical. Enable additional validation checks. Default: FALSE

- use_twostage:

  Logical. Use two-stage adaptive algorithm. Default: FALSE

- twostage_args:

  List. Parameters for two-stage algorithm:

  n.splits.screen

  :   Splits for Stage 1 screening. Default: 30

  screen.threshold

  :   Consistency threshold for Stage 1. Default: auto

  batch.size

  :   Splits per batch in Stage 2. Default: 20

  conf.level

  :   Confidence level for early stopping. Default: 0.95

  min.valid.screen

  :   Minimum valid Stage 1 splits. Default: 10

- parallel_args:

  List. Parallel processing configuration:

  plan

  :   Future plan: "multisession", "multicore", or "sequential"

  workers

  :   Number of parallel workers

  batch_size

  :   Number of candidates to evaluate per batch. Smaller values provide
      finer-grained early stopping but may increase overhead. Default:
      When stop_threshold is set and sg_focus is "hr" or "minSG",
      defaults to 1 (stop immediately when first candidate passes). For
      other sg_focus values with stop_threshold, defaults to
      min(workers, n_candidates/4). When stop_threshold is NULL,
      defaults to workers\*2 for efficiency.

  show_message

  :   Print parallel config messages

## Value

A list containing:

- out_sg:

  Selected subgroup results

- sg_focus:

  Selection criterion used

- df_flag:

  Data frame with treatment recommendations

- sg.harm:

  Subgroup definition labels

- sg.harm.id:

  Subgroup membership indicator

- algorithm:

  "twostage" or "fixed"

- n_candidates_evaluated:

  Number of candidates actually evaluated

- n_candidates_total:

  Total candidates available

- n_passed:

  Number meeting consistency threshold

- early_stop_triggered:

  Logical indicating if early stop occurred

- early_stop_candidate:

  Index of candidate triggering early stop

- stop_threshold:

  Threshold used for early stopping

- seed:

  Random seed used for reproducibility (NULL if not set)

## Examples

``` r
if (FALSE) { # \dontrun{
# Standard evaluation
result <- subgroup.consistency(
  df = trial_data,
  hr.subgroups = candidates,
  sg_focus = "hr",
  n.splits = 400,
  parallel_args = list(plan = "multisession", workers = 6)
)

# Show top 10 candidates before evaluation
result <- subgroup.consistency(
  df = trial_data,
  hr.subgroups = candidates,
  sg_focus = "hr",
  showten_subgroups = TRUE,  # Display candidates
  n.splits = 400
)

# With early stopping and custom batch size
result <- subgroup.consistency(
  df = trial_data,
  hr.subgroups = candidates,
  sg_focus = "hr",
  stop_threshold = 0.95,
  showten_subgroups = TRUE,
  parallel_args = list(
    plan = "multisession",
    workers = 6,
    batch_size = 2  # Check early stopping after every 2 candidates
  )
)
} # }
```
