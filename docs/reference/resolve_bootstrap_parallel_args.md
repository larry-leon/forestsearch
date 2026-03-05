# Resolve parallel processing arguments for bootstrap

If parallel_args not provided, falls back to forestsearch call's
parallel configuration. Always reports configuration to user.

## Usage

``` r
resolve_bootstrap_parallel_args(parallel_args, forestsearch_call_args)
```

## Arguments

- parallel_args:

  List or empty list

- forestsearch_call_args:

  List from original forestsearch call

## Value

List with plan, workers, show_message
