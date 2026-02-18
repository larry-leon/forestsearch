# Extract HR from DGM (Backward Compatible)

Extracts hazard ratios from DGM object, supporting both old and new
formats. Also supports CDE (controlled direct effect) extraction for
Table 5 of Leon et al. (2024) alignment (theta-ddagger).

## Usage

``` r
get_dgm_hr(dgm, which = "hr_H")
```

## Arguments

- dgm:

  DGM object (gbsg_dgm or aft_dgm_flex)

- which:

  Character. Which HR to extract: `"hr_H"`, `"hr_Hc"`, `"ahr_H"`,
  `"ahr_Hc"`, `"hr_overall"`, `"ahr"`, `"cde_H"`, `"cde_Hc"`, `"cde"`.

## Value

Numeric hazard ratio value
