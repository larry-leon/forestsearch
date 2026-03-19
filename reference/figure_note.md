# Generate Figure Note for Quarto/RMarkdown

Formats the figure note from plot_sg_weighted_km() output for use in
Quarto or RMarkdown documents.

## Usage

``` r
figure_note(
  x,
  prefix = "*Note*: ",
  include_definition = TRUE,
  include_hr_explanation = TRUE,
  custom_text = NULL
)
```

## Arguments

- x:

  Output from plot_sg_weighted_km()

- prefix:

  Character. Prefix for the note. Default uses italic Note.

- include_definition:

  Logical. Include subgroup definition. Default: TRUE

- include_hr_explanation:

  Logical. Include HR(bc) explanation. Default: TRUE

- custom_text:

  Character. Additional custom text to append. Default: NULL

## Value

Character string formatted as a figure note, or NULL if no content
