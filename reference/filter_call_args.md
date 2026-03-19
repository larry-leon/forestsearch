# Filter and merge arguments for function calls

Simplifies the common pattern of filtering arguments from a source list
to match a target function's formal parameters, then adding/overriding
specific arguments.

## Usage

``` r
filter_call_args(source_args, target_func, override_args = NULL)
```

## Arguments

- source_args:

  List of all arguments (typically from
  [`mget()`](https://rdrr.io/r/base/get.html) or a stored args list).

- target_func:

  Function whose formals define the filter criteria.

- override_args:

  List of arguments to add or override (optional).

## Value

List of filtered arguments ready for
[`do.call()`](https://rdrr.io/r/base/do.call.html).

## Details

This function:

1.  Extracts formal parameter names from `target_func`

2.  Keeps only arguments from `source_args` that match those names

3.  Adds or overrides with any `override_args` provided

Reduces boilerplate and improves readability across the codebase.
