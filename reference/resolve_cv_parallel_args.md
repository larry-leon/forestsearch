# Resolve Parallel Arguments for Cross-Validation

Helper function to resolve and validate parallel processing arguments,
similar to bootstrap's `resolve_bootstrap_parallel_args`.

## Usage

``` r
resolve_cv_parallel_args(parallel_args, fs_args, details = FALSE)
```

## Arguments

- parallel_args:

  List. User-provided parallel arguments.

- fs_args:

  List. Original ForestSearch call arguments.

- details:

  Logical. Print configuration messages.

## Value

List with resolved plan, workers, show_message.
