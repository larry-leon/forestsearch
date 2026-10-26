# Collect Results from foreach with Error Handling

Separates successful results from errors in the output of a
[`foreach`](https://rdrr.io/pkg/foreach/man/foreach.html) loop run with
`.errorhandling = "pass"`. Successful results are combined via
[`rbindlist`](https://rdrr.io/pkg/data.table/man/rbindlist.html) (with
an `rbind` fallback if data.table is not available) and failed tasks are
counted, reported, and optionally escalated to a hard error.

## Usage

``` r
collect_results(
  raw_list,
  label = "",
  stop_on_all_fail = TRUE,
  max_error_messages = 3L,
  keep_diagnostics = FALSE
)
```

## Arguments

- raw_list:

  List returned by `foreach(..., .errorhandling = "pass")`. Each element
  is either a data frame, data.table, or row-bindable object (success)
  or a condition object (failure).

- label:

  Character. Optional label included in the warning message to identify
  which simulation phase produced the errors (e.g. `"H1 alt"`,
  `"Bootstrap"`). Default `""`.

- stop_on_all_fail:

  Logical. If `TRUE` (default), raise an error when every task in
  `raw_list` failed. If `FALSE`, return an empty data frame and only
  emit a warning. Default `TRUE`.

- max_error_messages:

  Integer. Maximum number of distinct error messages to print in the
  diagnostic output. Default `3L`.

- keep_diagnostics:

  Logical. If `TRUE`, scan each successful element for a `"diagnostics"`
  attribute and collect them into a list attached to the returned object
  as `attr(result, "diagnostics")`. Used with
  [`run_simulation_analysis`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)'s
  `keep` / `keep_first_n` mechanism to retain per-replicate heavy
  objects (FS or GRF result objects, candidate tables) that
  [`data.table::rbindlist()`](https://rdrr.io/pkg/data.table/man/rbindlist.html)
  would otherwise drop. Default `FALSE` (preserves minimal-memory
  behavior).

## Value

A data frame of row-bound successful results. Always a data frame;
**never** `NULL` – an empty data frame is returned when every task
failed and `stop_on_all_fail = FALSE`. Two attributes are always
attached:

- `n_failed`:

  Integer. Number of tasks that returned errors.

- `n_total`:

  Integer. Total number of tasks (successes + failures).

When `keep_diagnostics = TRUE`, a third attribute `"diagnostics"` may be
attached – a named list of per-replicate diagnostic objects, keyed by
`"sim<ID>"` when the source element carries an `"sim_id"` attribute
([`run_simulation_analysis`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)
sets this automatically).

## Details

This replaces the fragile pattern of calling
[`nrow()`](https://rdrr.io/r/base/nrow.html) or
[`rbind()`](https://rdrr.io/r/base/cbind.html) directly on a raw
`foreach` result list, which fails silently when error objects are mixed
in with data-frame rows.

Design rationale for the reporting strategy:

- A [`warning`](https://rdrr.io/r/base/warning.html) is emitted so that
  calling code (and `tryCatch` / `withCallingHandlers` blocks) can
  capture the failure condition programmatically.

- The error summary is also printed with
  [`cat`](https://rdrr.io/r/base/cat.html) to `stderr`, because
  [`warning`](https://rdrr.io/r/base/warning.html) and
  [`message`](https://rdrr.io/r/base/message.html) output is sometimes
  buffered or suppressed inside future/foreach parallel workers or when
  `options(warn = -1)` is set by user code upstream. The
  [`cat()`](https://rdrr.io/r/base/cat.html) output is the
  visible-no-matter-what fallback.

- When `stop_on_all_fail = TRUE` (the default), a completely empty
  result list is escalated to an error rather than silently returning an
  empty data frame. This avoids downstream code proceeding with zero
  rows and producing misleading summaries. Calibration code that may
  legitimately tolerate some failures (e.g. Monte Carlo FPR estimation)
  can set `stop_on_all_fail = FALSE`.

The returned object is always a data frame – never `NULL` – because
attributes (`n_failed`, `n_total`) attached to `NULL` are silently
discarded by R, whereas attributes on an empty data frame are preserved
and inspectable by the caller.

## Preserving per-replicate diagnostics

[`data.table::rbindlist()`](https://rdrr.io/pkg/data.table/man/rbindlist.html)
drops per-element attributes during row-binding. When the caller has
populated each replicate with heavy diagnostic objects (e.g. via
[`run_simulation_analysis`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)
with a non-empty `keep` argument), set `keep_diagnostics = TRUE` to
retain them. The diagnostics attribute on each input element must carry
`attr(..., "sim_id")` for stable keying; that is set automatically by
[`run_simulation_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md).

## See also

[`reset_workers`](https://larry-leon.github.io/forestsearch/reference/reset_workers.md)
for resetting the future plan and releasing worker memory between
simulation phases.

## Examples

``` r
# \donttest{
library(foreach)
library(doFuture)
#> Loading required package: future
future::plan("sequential")

raw <- foreach(i = 1:5, .errorhandling = "pass") %dofuture% {
  if (i == 3L) stop("simulated failure")
  data.frame(task = i, value = rnorm(1))
}
#> Warning: UNRELIABLE VALUE: At least one of iterations 1-5 of the foreach() %dofuture% { ... }, part of chunk #1 (‘doFuture2-1’), unexpectedly generated random numbers without declaring so. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify foreach() argument '.options.future = list(seed = TRUE)'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, set option 'doFuture.rng.onMisuse' to "ignore". [future ‘doFuture2-1’ (7e617749df86f8c8556890313ff90163-1); on 7e617749df86f8c8556890313ff90163@runnervmtai0g<30337>]
res <- collect_results(raw, label = "demo")
#> Warning: demo: 1 of 5 tasks failed. First error: simulated failure
# Warning: demo: 1 of 5 tasks failed.
nrow(res)
#> [1] 4
# [1] 4
attr(res, "n_failed")
#> [1] 1
# [1] 1
# }
```
