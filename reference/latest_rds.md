# Find the Most Recent .rds File Matching a Pattern

Scans `output_dir` for files matching the given prefix pattern and
returns the path to the most recently modified one. Date-stamped
filenames (e.g., `results_20260402.rds`) sort naturally, but this
function uses filesystem modification time for robustness.

## Usage

``` r
latest_rds(prefix, output_dir = "_data", verbose = TRUE)
```

## Arguments

- prefix:

  Character. Filename prefix to match (e.g.,
  `"selection_survival_frontier"`).

- output_dir:

  Character. Directory to scan. Default: `"_data"`.

- verbose:

  Logical. Print the matched file. Default: TRUE.

## Value

Full path to the most recent matching file, or NULL if no match is
found.

## Examples

``` r
if (FALSE) { # \dontrun{
path <- latest_rds("selection_survival_frontier")
res  <- readRDS(path)

path <- latest_rds("calibration_binary_leff_grid")
} # }
```
