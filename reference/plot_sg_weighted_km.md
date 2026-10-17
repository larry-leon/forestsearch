# Plot Weighted Kaplan-Meier Curves for ForestSearch Subgroups

Creates weighted Kaplan-Meier survival curves for the identified
subgroups (H and Hc) using the weightedsurv package, matching the
pattern used in
[`sg_consistency_out()`](https://larry-leon.github.io/forestsearch/reference/sg_consistency_out.md).

## Usage

``` r
plot_sg_weighted_km(
  fs.est,
  fs_bc = NULL,
  outcome.name = "Y",
  event.name = "Event",
  treat.name = "Treat",
  by.risk = NULL,
  sg0_name = NULL,
  sg1_name = NULL,
  conf.int = TRUE,
  show.logrank = TRUE,
  show.cox = TRUE,
  show.cox.bc = TRUE,
  put.legend.lr = "topleft",
  ymax = 1.05,
  xmed.fraction = 0.65,
  hr_bc_position = "bottomright",
  hr_bc_cex = 0.725,
  title = NULL,
  verbose = FALSE
)
```

## Arguments

- fs.est:

  A forestsearch object containing `df.est` with `treat.recommend`
  column, or a data frame directly.

- fs_bc:

  Optional. Bootstrap results from
  [`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
  containing bias-corrected HR estimates. If provided, bias-corrected
  HRs will be annotated on the plots.

- outcome.name:

  Character. Name of time-to-event column. Default: `"Y"`

- event.name:

  Character. Name of event indicator column. Default: `"Event"`

- treat.name:

  Character. Name of treatment column. Default: `"Treat"`

- by.risk:

  Numeric. Risk interval for plotting. Default: `NULL` (auto-calculated
  as max(outcome)/12)

- sg0_name:

  Character. Label for H subgroup (treat.recommend == 0). Default:
  `NULL` (auto-extracted from forestsearch object as "H: {definition}"
  or "Questionable (H)" if not available)

- sg1_name:

  Character. Label for Hc subgroup (treat.recommend == 1). Default:
  `NULL` (auto-generated as "Hc: NOT {definition}" or "Recommend (Hc)"
  if not available)

- conf.int:

  Logical. Show confidence intervals. Default: `TRUE`

- show.logrank:

  Logical. Show log-rank test. Default: `TRUE`

- show.cox:

  Logical. Show unadjusted Cox HR from weightedsurv. Default: `TRUE`

- show.cox.bc:

  Logical. Show bootstrap bias-corrected HR annotation (requires
  `fs_bc`). Default: `TRUE`

- put.legend.lr:

  Character. Legend position. Default: "topleft"

- ymax:

  Numeric. Max y-axis value. Default: 1.05

- xmed.fraction:

  Numeric. Fraction for median lines. Default: 0.65

- hr_bc_position:

  Character. Position for bias-corrected HR annotation. One of
  "bottomright", "bottomleft", "topright", "topleft". Default:
  "bottomright"

- hr_bc_cex:

  Numeric. Character expansion factor for bias-corrected HR annotation
  text. Default: 0.725 (matches weightedsurv cox.cex default)

- title:

  Character. Overall plot title. Default: `NULL`

- verbose:

  Logical. Print diagnostic messages. Default: `FALSE`

## Value

Invisibly returns a list with subgroup data frames and counting data

## Details

This function uses the exact same calling pattern as
[`plot_subgroup()`](https://larry-leon.github.io/forestsearch/reference/plot_subgroup.md)
in the ForestSearch package. Column names are mapped internally to the
standard names (Y, Event, Treat) expected by weightedsurv.

Subgroup definitions are automatically extracted from the forestsearch
object if available:

- `fs$grp.consistency$out_sg$sg.harm_label` - Human-readable labels

- `fs$sg.harm` - Technical factor names (fallback)

HR display options controlled by `show.cox` and `show.cox.bc`:

- Both TRUE (default): Shows unadjusted HR from weightedsurv AND
  bias-corrected HR annotation

- `show.cox = TRUE, show.cox.bc = FALSE`: Shows only unadjusted HR

- `show.cox = FALSE, show.cox.bc = TRUE`: Shows only bias-corrected HR

- Both FALSE: Shows neither HR estimate

## Examples

``` r
if (FALSE) { # \dontrun{
# After running forestsearch - auto-extracts subgroup definition
plot_sg_weighted_km(fs.est = fs)

# With bootstrap bias-corrected estimates
plot_sg_weighted_km(fs.est = fs, fs_bc = fs_bootstrap)

# With custom column names
plot_sg_weighted_km(
  fs.est = fs,
  outcome.name = "time_months",
  event.name = "status",
  treat.name = "hormon"
)
} # }
```
