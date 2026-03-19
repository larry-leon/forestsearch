# Automatically Detect Variable Types in a Dataset

Analyzes a data frame to automatically classify variables as continuous
or categorical, and returns a subset of the data with specified
variables excluded.

## Usage

``` r
detect_variable_types(data, max_unique_for_cat = 10, exclude_vars = NULL)
```

## Arguments

- data:

  A data frame to analyze

- max_unique_for_cat:

  Integer. Maximum number of unique values for a numeric variable to be
  considered categorical. Default is 10.

- exclude_vars:

  Character vector of variable names to exclude from both classification
  and the returned dataset (e.g., ID variables, timestamps). Default is
  NULL.

## Value

A list containing:

- continuous_vars:

  Character vector of variable names classified as continuous

- cat_vars:

  Character vector of variable names classified as categorical

- data_subset:

  Data frame with exclude_vars columns removed

## Details

The function classifies variables using the following rules:

- Numeric variables with more than `max_unique_for_cat` unique values
  are classified as continuous

- Numeric variables with `max_unique_for_cat` or fewer unique values are
  classified as categorical

- Factor, character, and logical variables are always classified as
  categorical

- Variables listed in `exclude_vars` are omitted from classification and
  removed from the returned dataset
