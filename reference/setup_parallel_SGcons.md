# Set up parallel processing for subgroup consistency

Sets up parallel processing using the specified approach and number of
workers.

## Usage

``` r
setup_parallel_SGcons(
  parallel_args = list(plan = "multisession", workers = 4, show_message = TRUE)
)
```

## Arguments

- parallel_args:

  List with `plan` (character), `workers` (integer), and `show_message`
  (logical).

## Value

None. Sets up parallel backend as side effect.
