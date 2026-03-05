# Generate Prediction Dataset with Subgroup Treatment Recommendation

Creates a prediction dataset with a treatment recommendation flag based
on the subgroup definition. Supports both label expressions (e.g.,
`"\{er <= 0\}"`) and bare column names (e.g., `"q3.1"`).

## Usage

``` r
get_dfpred(df.predict, sg.harm, version = 1)
```

## Arguments

- df.predict:

  Data frame for prediction (test or validation set).

- sg.harm:

  Character vector of subgroup-defining labels. Values may be wrapped in
  braces and optionally negated, e.g. `"\{er <= 0\}"` or
  `"!\{size <= 35\}"`. Plain column names (e.g., `"q3.1"`) are treated
  as binary indicators that must equal 1.

- version:

  Integer; encoding version (maintained for backward compatibility).
  Default: 1.

## Value

Data frame with treatment recommendation flag (`treat.recommend`): 0 for
harm subgroup, 1 for complement.

## Details

Each element of `sg.harm` is processed as follows:

1.  Outer braces and leading `!` are stripped.

2.  If the result matches `"var op value"` (where `op` is one of `<=`,
    `<`, `>=`, `>`, `==`, `!=`), the comparison is executed directly on
    `df.predict[[var]]`.

3.  Otherwise the expression is treated as a column name and membership
    is `df.predict[[name]] == 1`.

## See also

[`evaluate_comparison`](https://larry-leon.github.io/forestsearch/reference/evaluate_comparison.md)
for the operator-dispatch logic,
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
for the main analysis function.

## Examples

``` r
if (FALSE) { # \dontrun{
# With brace-wrapped label expressions
sg <- c("{er <= 0}", "{size <= 35}")
df_out <- get_dfpred(df.predict = test_data, sg.harm = sg)

# With negation
sg_neg <- c("{er <= 0}", "!{size <= 35}")
df_neg <- get_dfpred(df.predict = test_data, sg.harm = sg_neg)

# With bare column names (binary indicators)
sg_col <- c("q1.1", "q3.1")
df_col <- get_dfpred(df.predict = encoded_data, sg.harm = sg_col)
} # }
```
