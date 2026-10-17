# Identification-structure table

The supplement's structure table: detection rate, the five structure
shares, anchor recovery, and the conditional partner-versus-proxy split,
one column per sweep point. Returned as a data frame by default so a
LaTeX document can style it however it likes; `as_gt = TRUE` renders it
with gt.

## Usage

``` r
fs_identification_table(x, as_gt = FALSE, digits = 1)
```

## Arguments

- x:

  An `fs_identification` object.

- as_gt:

  Logical. Render with gt instead of returning the data frame.

- digits:

  Integer. Digits on the percent labels.

## Value

A data frame, or a `gt_tbl`. The accuracy sentence the supplement
carries as a table footnote is attached as the `"accuracy_note"`
attribute of the data frame.
