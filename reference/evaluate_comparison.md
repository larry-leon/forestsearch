# Evaluate a Comparison Expression Without eval(parse())

Parses a string of the form `"var op value"` and evaluates it directly
against a data frame column using operator dispatch. Falls back to
column-name lookup for bare names.

## Usage

``` r
evaluate_comparison(expr, df)
```

## Arguments

- expr:

  Character. An expression like `"er <= 0"`, `"size > 35"`,
  `"grade3 == 1"`, or a bare column name like `"q3.1"`.

- df:

  Data frame whose columns are referenced by `expr`.

## Value

Logical vector of length `nrow(df)`.

## Details

Supported operators (matched longest-first to avoid partial-match
ambiguity): `<=`, `>=`, `!=`, `==`, `<`, `>`.

If no operator is found, `expr` is treated as a column name and the
result is `df[[expr]] == 1`.

The value on the right-hand side is coerced to numeric when possible,
otherwise kept as character for string comparisons.
