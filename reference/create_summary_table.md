# Create Enhanced Summary Table for Baseline Characteristics

Generates a formatted summary table comparing baseline characteristics
between treatment arms. Supports continuous, categorical, and binary
variables with p-values, standardized mean differences (SMD), and
missing data summaries.

## Usage

``` r
create_summary_table(
  data,
  treat_var = "treat",
  vars_continuous = NULL,
  vars_categorical = NULL,
  vars_binary = NULL,
  var_labels = NULL,
  digits = 1,
  show_pvalue = TRUE,
  show_smd = TRUE,
  show_missing = TRUE,
  table_title = "Baseline Characteristics by Treatment Arm",
  table_subtitle = NULL,
  source_note = NULL,
  font_size = 12,
  header_font_size = 14,
  footnote_font_size = 10,
  use_alternating_rows = TRUE,
  stripe_color = "#f9f9f9",
  indent_size = 20,
  highlight_pval = 0.05,
  highlight_smd = 0.2,
  highlight_color = "#fff3cd",
  compact_mode = FALSE,
  column_width_var = 200,
  column_width_stats = 120,
  show_column_borders = FALSE,
  custom_css = NULL
)
```

## Arguments

- data:

  Data frame containing the analysis data

- treat_var:

  Character. Name of treatment variable (must have 2 levels)

- vars_continuous:

  Character vector. Names of continuous variables

- vars_categorical:

  Character vector. Names of categorical variables

- vars_binary:

  Character vector. Names of binary (0/1) variables

- var_labels:

  Named list. Custom labels for variables (optional)

- digits:

  Integer. Number of decimal places for continuous variables

- show_pvalue:

  Logical. Include p-values column

- show_smd:

  Logical. Include SMD (effect size) column

- show_missing:

  Logical. Include missing data rows

- table_title:

  Character. Main title for the table

- table_subtitle:

  Character. Subtitle for the table (optional)

- source_note:

  Character. Source note at bottom (optional)

- font_size:

  Numeric. Base font size in pixels (default: 12)

- header_font_size:

  Numeric. Header font size in pixels (default: 14)

- footnote_font_size:

  Numeric. Footnote font size in pixels (default: 10)

- use_alternating_rows:

  Logical. Apply zebra striping (default: TRUE)

- stripe_color:

  Character. Color for alternating rows (default: "#f9f9f9")

- indent_size:

  Numeric. Indentation for sub-levels in pixels (default: 20)

- highlight_pval:

  Numeric. Highlight p-values below this threshold (default: 0.05)

- highlight_smd:

  Numeric. Highlight SMD values above this threshold (default: 0.2)

- highlight_color:

  Character. Color for highlighting (default: "#fff3cd")

- compact_mode:

  Logical. Reduce spacing for compact display (default: FALSE)

- column_width_var:

  Numeric. Width for Variable column in pixels (default: 200)

- column_width_stats:

  Numeric. Width for stat columns in pixels (default: 120)

- show_column_borders:

  Logical. Show vertical column borders (default: FALSE)

- custom_css:

  Character. Additional custom CSS styling (optional)

## Value

A gt table object (or data frame if gt not available)

## Details

Binary variables specified via `vars_binary` display a single row
showing the count and proportion for the "1" level. Categorical
variables specified via `vars_categorical` that happen to be
binary-coded (i.e., have exactly two levels: 0 and 1) are automatically
detected and displayed in the same compact single-row format, showing
only the "1" proportion.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
create_summary_table(
  data = trial_data,
  treat_var = "treatment",
  vars_continuous = c("age", "bmi"),
  vars_categorical = c("sex", "stage")
)

# Customized appearance
create_summary_table(
  data = trial_data,
  treat_var = "treatment",
  vars_continuous = c("age", "bmi"),
  vars_categorical = c("sex", "stage"),
  font_size = 11,
  header_font_size = 13,
  use_alternating_rows = TRUE,
  highlight_pval = 0.05,
  compact_mode = TRUE
)
} # }
```
