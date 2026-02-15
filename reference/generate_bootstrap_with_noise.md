# Generate Bootstrap Sample with Added Noise

Creates a bootstrap sample from a dataset with controlled noise added to
both continuous and categorical variables. This function is useful for
generating synthetic datasets that maintain the general structure of the
original data while introducing controlled variation.

## Usage

``` r
generate_bootstrap_with_noise(
  data,
  n = NULL,
  continuous_vars = NULL,
  cat_vars = NULL,
  id_var = "pid",
  seed = 123,
  noise_level = 0.1
)
```

## Arguments

- data:

  A data frame containing the original dataset to bootstrap from.

- n:

  Integer. Number of observations in the output dataset. If NULL
  (default), uses the same number of rows as the input data.

- continuous_vars:

  Character vector of column names to treat as continuous variables. If
  NULL (default), automatically detects numeric columns.

- cat_vars:

  Character vector of column names to treat as categorical variables. If
  NULL (default), automatically detects factors, logical columns, and
  numeric columns with 10 or fewer unique values.

- id_var:

  Character string specifying the name of the ID variable column. This
  column will be reset to sequential values (1:n) in the output. Default
  is "pid".

- seed:

  Integer. Random seed for reproducibility. Default is 123.

- noise_level:

  Numeric between 0 and 1. Controls the amount of noise added. For
  continuous variables, this is multiplied by the standard deviation to
  determine noise magnitude. For categorical variables, this is divided
  by 2 to determine the probability of value changes. Default is 0.1.

## Value

A data frame with the same structure as the input data, containing
bootstrap sampled observations with added noise.

## Details

The function performs the following operations:

### Bootstrap Sampling

Samples n observations with replacement from the original dataset.

### Continuous Variable Noise

- Adds Gaussian noise with standard deviation = original SD ×
  noise_level

- Constrains values to remain within original variable bounds

- Preserves integer type for variables that appear to be integers

### Categorical Variable Perturbation

- Changes values with probability = noise_level / 2

- Binary variables: flips to opposite value

- Multi-level unordered: randomly selects from other levels

- Ordered factors: weights selection toward adjacent levels

- Preserves factor levels and ordering from original data

## Note

- The function assumes that categorical variables with numeric encoding
  should maintain their numeric type unless they are factors in the
  input

- Missing values (NA) are handled appropriately in calculations but are
  not imputed

- For ordered factors or variables named "grade", the perturbation
  favors transitions to adjacent levels over distant levels

## See also

[`sample`](https://rdrr.io/r/base/sample.html) for bootstrap sampling,
[`rnorm`](https://rdrr.io/r/stats/Normal.html) for noise generation

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example dataset
data(gbsg, package = "survival")

# Basic usage with automatic variable detection
synthetic_data <- generate_bootstrap_with_noise(
  data = gbsg,
  seed = 123
)

# Specify variables explicitly
synthetic_data <- generate_bootstrap_with_noise(
  data = gbsg,
  n = 1000,
  continuous_vars = c("age", "size", "nodes", "pgr", "er", "rfstime"),
  cat_vars = c("meno", "grade", "hormon", "status"),
  id_var = "pid",
  seed = 456,
  noise_level = 0.15
)

# Create multiple synthetic datasets
synthetic_list <- lapply(1:10, function(i) {
  generate_bootstrap_with_noise(data = gbsg, seed = i)
})
} # }
```
