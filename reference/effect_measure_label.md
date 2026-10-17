# Human-Readable Effect Measure Label

Returns a descriptive label for use in table footnotes and narrative
text. Falls back to `"HR"` for survival analyses.

## Usage

``` r
effect_measure_label(effect_measure = NULL, outcome_type = NULL)
```

## Arguments

- effect_measure:

  Character or `NULL`.

- outcome_type:

  Character or `NULL`.

## Value

Character label (e.g., `"Cox HR"`, `"OR"`, `"mean difference"`).
