# Format ForestSearch Details Output for Beamer Two-Column Display

Captures `forestsearch(details = TRUE)` console output and splits it
into two columns for readable beamer slides. Left column shows variable
selection (GRF, LASSO, candidate factors); right column shows subgroup
search, consistency evaluation, and results.

## Usage

``` r
format_fs_details(
  fs_output,
  split_after = "Candidate factors",
  fontsize = "scriptsize",
  col_widths = c(0.48, 0.52),
  max_width = 48
)
```

## Arguments

- fs_output:

  Character vector of captured output lines from
  `capture.output(forestsearch(..., details = TRUE))`.

- split_after:

  Character string (regex). The output is split after the block matching
  this pattern. Default: `"Candidate factors"` which splits after the
  candidate factor list.

- fontsize:

  Character. LaTeX font size for the output text. One of `"tiny"`,
  `"scriptsize"`, `"footnotesize"`, `"small"`. Default: `"scriptsize"`.

- col_widths:

  Numeric vector of length 2. Column widths as fractions of
  `\textwidth`. Default: `c(0.48, 0.52)`.

- max_width:

  Integer. Maximum character width per line before wrapping. Long lines
  are wrapped at comma or space boundaries with a 4-space continuation
  indent. Default: 48 (suitable for half-slide columns at scriptsize).

## Value

Invisibly returns a list with `left` and `right` character vectors. Side
effect: emits LaTeX via [`cat()`](https://rdrr.io/r/base/cat.html) for
use in a chunk with `results='asis'`.

## Quarto Setup

No special LaTeX packages required. Works in any beamer frame without
the `fragile` option.

## Usage

In a Quarto beamer chunk with `results='asis'` and `echo=FALSE`, first
capture the forestsearch output with
[`capture.output()`](https://rdrr.io/r/utils/capture.output.html), then
call `format_fs_details(fs_output)`.

## Examples

``` r
if (FALSE) { # \dontrun{
fs_output <- capture.output(
  fs_res <- forestsearch(
    df = dat, outcome.name = "time", event.name = "status",
    treat.name = "trt", id.name = "id",
    confounders.name = confs, details = TRUE
  )
)
# Default: split after candidate factors
format_fs_details(fs_output, fontsize = "tiny")

# Custom split: after filtering summary
format_fs_details(fs_output, split_after = "Found.*candidate")
} # }
```
