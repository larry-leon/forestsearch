# The frame `beta(Hhat)` is scored on

One entry point per outcome family, dispatching on `outcome_type` rather
than exposing a separate function per family. Defaults reproduce the
simulation modules this replaces exactly.

## Usage

``` r
fs_build_eval_frame(
  dgm,
  outcome_type = c("survival", "binary", "continuous", "count"),
  eval_seed = 20260628L,
  analysis_time = 84,
  cens_adjust = log(1.5),
  n_eval = NULL
)
```

## Arguments

- dgm:

  A DGM object carrying `df_super`.

- outcome_type:

  Character. One of `"survival"`, `"binary"`, `"continuous"`, `"count"`.

- eval_seed:

  Integer. Fixes the single realization on the full pool. Ignored for
  `"continuous"` and `"count"`.

- analysis_time, cens_adjust:

  Survival only.

- n_eval:

  Defunct. Non-`NULL` is a hard error: this always evaluates the full
  super-population, so both engines score the identical target.

## Value

A data frame: the evaluation population.

## What each family returns

- **survival** – the entire fixed super-population, every subject
  exactly once (`replace = FALSE`), under the same randomized/censored
  analysis the trials run.

- **binary** – the same construction on the GLM simulator: every subject
  once, under one fixed treatment/outcome realization.

- **continuous / count** – `dgm$df_super` **unchanged**. No simulation
  occurs. The mean difference is collapsible, so the target is an exact
  finite mean over the super-population: the scoring frame *is* the
  population. `eval_seed` is accepted and ignored here so that generic
  harness code need not branch, and the target carries **zero Monte
  Carlo error** – uniformity of the call surface is not sameness of the
  object.

## Rejected arguments

`analysis_time`, `cens_adjust` and `n_eval` are survival-only. Supplying
any of them on a non-survival path is an error rather than a silent
no-op: quietly ignoring an argument that means something on another path
is how conventions drift apart.
