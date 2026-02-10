# Extract HR from DGM (Backward Compatible)

Extracts hazard ratios from DGM object, supporting both old and new
formats.

## Usage

``` r
get_dgm_hr(dgm, which = "hr_H")
```

## Arguments

- dgm:

  DGM object (gbsg_dgm or aft_dgm_flex)

- which:

  Character. Which HR to extract: "hr_H", "hr_Hc", "ahr_H", "ahr_Hc",
  "hr_overall", "ahr"

## Value

Numeric hazard ratio value
