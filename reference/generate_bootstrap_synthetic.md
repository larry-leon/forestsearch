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

  Number of synthetic observations to generate (default: same as
  original)

- seed:

  Random seed for reproducibility

- noise_level:

  Noise level for perturbation (0 to 1, default 0.1)

- id_var:

  Optional name of ID variable to regenerate (will be numbered 1:n)

- cat_flip_prob:

  Probability of flipping categorical values (default: noise_level/2)

- preserve_bounds:

  Logical: should continuous variables stay within original bounds?
  (default: TRUE)

- ordinal_vars:

  Optional character vector of ordinal categorical variables (these will
  be perturbed to adjacent values rather than randomly flipped)

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
  ordinal_vars = c("grade"),
  id_var = "pid",
  n = 1000,
  seed = 123,
  noise_level = 0.15
)
#> Warning: internal error 1 in R_decompress1 with libdeflate
#> Error: lazy-load database '/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/forestsearch/R/forestsearch.rdb' is corrupt

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
  ordinal_vars = c("education"),
  id_var = "id",
  n = 150,
  seed = 456
)
#> Warning: restarting interrupted promise evaluation
#> Warning: internal error 1 in R_decompress1 with libdeflate
#> Error: lazy-load database '/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/forestsearch/R/forestsearch.rdb' is corrupt
# }
```
