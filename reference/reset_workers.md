# Reset a Parallel Worker Pool

Tears down the current future plan, forces garbage collection,
optionally pauses, and spins up a fresh worker pool. Used to release
worker-side memory between simulation phases (e.g. between alternative
and null Monte Carlo runs) so that accumulated R-session state in
long-running workers does not bloat memory or cause sporadic failures.

## Usage

``` r
reset_workers(
  workers = NULL,
  strategy = "multisession",
  sleep = 1,
  gc = TRUE,
  quiet = FALSE
)
```

## Arguments

- workers:

  Integer. Number of workers for the new pool. Default `NULL` leaves the
  current plan's worker count unchanged (useful if a reset without
  resize is wanted; this requires `future` to be installed).

- strategy:

  Character. The future strategy to reinstate. Default `"multisession"`.
  Other reasonable choices are `"multicore"` (Linux/macOS only) and
  `"sequential"` (to disable parallelism entirely).

- sleep:

  Numeric. Seconds to pause between teardown and setup. Default `1`. Set
  to `0` to skip.

- gc:

  Logical. Whether to force garbage collection between teardown and
  setup. Default `TRUE`.

- quiet:

  Logical. Suppress informational
  [`message()`](https://rdrr.io/r/base/message.html) output. Default
  `FALSE`.

## Value

Invisibly, a list with the previous and current plan class names:
`list(previous = <chr>, current = <chr>)`.

## Details

The worker reset sequence is:

1.  `future::plan("sequential")` – dismantles the existing worker pool;
    R processes terminate and their memory is returned to the OS.

2.  `gc(verbose = FALSE)` – reclaims main-process memory holding
    references to now-defunct futures.

3.  `Sys.sleep(sleep)` – brief pause to let the OS reclaim sockets /
    pipes before the new pool is requested. A small non-zero sleep (the
    default `1`) materially reduces `multisession` setup failures
    observed on busy Linux servers; set to `0` to skip.

4.  `future::plan(strategy, workers = workers)` – starts a fresh pool.

This function is a no-op if `future` is not installed or if the
requested strategy cannot be initialised; a warning is emitted and the
previous plan is left in place.

## See also

[`collect_results`](https://larry-leon.github.io/forestsearch/reference/collect_results.md)
for error-tolerant collection of `foreach` output.

## Examples

``` r
# \donttest{
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("multisession", workers = 2L)
  # ... run alternative simulations ...
  reset_workers(workers = 2L)
  # ... run null simulations ...
  future::plan("sequential")
}
#> reset_workers: plan reset to FutureStrategy with 2 worker(s) (was FutureStrategy).
# }
```
