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
