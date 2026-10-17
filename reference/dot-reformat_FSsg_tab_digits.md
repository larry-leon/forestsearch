# Re-format Numeric Cells in FSsg_tab at a Lower Precision

Parses formatted-string cells in `FSsg_tab` back to numerics and
re-formats at the requested `digits`. Designed for display-time
precision control in
[`summarize_bootstrap_results`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_results.md).

## Usage

``` r
.reformat_FSsg_tab_digits(tab, digits)
```

## Arguments

- tab:

  Data frame or matrix with formatted-string cells.

- digits:

  Integer. Decimal places.

## Value

Reformatted `tab`, same shape.

## Details

Cells are identified by column name pattern, not content sniffing. Three
cell formats are recognized:

1.  Single value with optional sign and percent suffix (Mean(C)/(T),
    Rate(C)/(T), Diff): `"-56.4321"`, `"+12.5"`, `"42.1%"`.

2.  Three-value CI (effect-measure column, `"\{measure\} (95% CI)"`):
    `"82.7831 (14.8067, 150.7585)"`.

3.  Three-value CI without space (bias-corrected `"\{measure\}*"`
    column): `"72.4500 (-61.5400,206.4400)"`.

Columns NOT touched: `Subgroup` (text label), `n` and `n1`
(count-with-percent format like `"87 (8.0%)"` — different semantics).

Robustness contract: cells that don't match any expected pattern are
returned unchanged. Construction-time precision must be at least
`digits` for the reformat to be informative; if construction was at
lower precision, this function will zero-pad to `digits` but cannot
recover information.
