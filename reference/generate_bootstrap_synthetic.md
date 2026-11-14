# Generate Synthetic Data using Bootstrap with Perturbation

Generate Synthetic Data using Bootstrap with Perturbation

## Usage

``` r
generate_bootstrap_synthetic(
  data,
  continuous_vars,
  cat_vars,
  n = NULL,
  seed = 123,
  noise_level = 0.1,
  id_var = NULL,
  cat_flip_prob = NULL,
  preserve_bounds = TRUE,
  ordinal_vars = NULL
)
```

## Arguments

- data:

  Original dataset to bootstrap from

- continuous_vars:

  Character vector of continuous variable names

- cat_vars:

  Character vector of categorical variable names

- n:

  Number of synthetic observations to generate (default: NULL, meaning
  same as original)

- seed:

  Random seed for reproducibility

- noise_level:

  Noise level for perturbation (0 to 1, default 0.1)

- id_var:

  Optional name of ID variable to regenerate (will be numbered 1:n)

- cat_flip_prob:

  Probability of flipping categorical values (default: NULL, giving
  noise_level/2)

- preserve_bounds:

  Logical: should continuous variables stay within original bounds?
  (default: TRUE)

- ordinal_vars:

  Optional character vector naming which of the `cat_vars` are ordinal.
  A numeric or character variable listed in both `cat_vars` and
  `ordinal_vars` is changed to a neighbouring level (in sorted order)
  instead of a random other level. A factor listed in both is set to NA,
  with a warning, at each row selected for change, unless its levels are
  "1", "2", ... in that order. A variable listed only in `ordinal_vars`
  is not perturbed: it is kept as resampled and is named in the 'kept
  as-is' message. Known limitation: for a numeric-coded variable at its
  lowest or highest level, the new value is not always the adjacent
  level (with codes 1, ..., k, a row at 1 that is selected for change
  keeps its value half the time, and a row at k moves to any of 1, ...,
  k - 1 with equal probability).

## Value

A data frame with synthetic data

## Examples

``` r
# \donttest{
# Example 1: Using with GBSG dataset
synth_gbsg <- generate_bootstrap_synthetic(
  data = survival::gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er", "rfstime"),
  cat_vars = c("meno", "hormon", "status"),
  # listed only in ordinal_vars, so it is resampled but not perturbed;
  # list it in cat_vars as well to perturb it
  ordinal_vars = c("grade"),
  id_var = "pid",
  n = 1000,
  seed = 123,
  noise_level = 0.15
)
#> Note: The following variables were not specified as continuous or categorical and will be kept as-is from bootstrap sample:
#>    grade

# Example 2: Using with any dataset
my_data <- data.frame(
  id = 1:100,
  height = rnorm(100, 170, 10),
  weight = rnorm(100, 70, 15),
  age = sample(20:80, 100, replace = TRUE),
  gender = sample(c("M", "F"), 100, replace = TRUE),
  education = sample(1:5, 100, replace = TRUE),
  smoker = sample(0:1, 100, replace = TRUE)
)

synth_data <- generate_bootstrap_synthetic(
  data = my_data,
  continuous_vars = c("height", "weight", "age"),
  cat_vars = c("gender", "smoker"),
  # listed only in ordinal_vars, so it is resampled but not perturbed;
  # list it in cat_vars as well to perturb it
  ordinal_vars = c("education"),
  id_var = "id",
  n = 150,
  seed = 456
)
#> Note: The following variables were not specified as continuous or categorical and will be kept as-is from bootstrap sample:
#>    education
# }
```
