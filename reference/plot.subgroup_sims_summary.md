# Forest plot of an extreme-subgroups summary

Draws one of the vignettes' six forest panels from a
[`summary.subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/summary.subgroup_sims.md)
object via
[`gg_forest()`](https://larry-leon.github.io/forestsearch/reference/gg_forest.md):
the empirical-CI distribution of the HR point estimate or of the upper
95% bound UB(HR), for the single-variable panel, the combination +
random- benchmark panel, or the high-risk panel (Pr(UB\>=2) filter with
the ITT row anchored first).

## Usage

``` r
# S3 method for class 'subgroup_sims_summary'
plot(
  x,
  metric = c("hr", "ub"),
  panel = c("single", "combo", "highrisk"),
  hr_true = x$hr_true,
  base_size = 14,
  point_size = 2.2,
  line_size = 0.75,
  widths = c(3.5, 5, 1.1, 1.2, 1.2, 1.1),
  clip_marker = "arrow",
  xlim = NULL,
  ticks_at = NULL,
  xlab = NULL,
  footnote = NULL,
  ...
)
```

## Arguments

- x:

  A `"subgroup_sims_summary"` object.

- metric:

  `"hr"` (HR point estimate) or `"ub"` (upper 95% bound).

- panel:

  `"single"`, `"combo"`, or `"highrisk"`.

- hr_true:

  True hazard ratio for the reference/vertical line and footnotes;
  defaults to the value stored on `x` and must be non-`NULL`.

- base_size, point_size, line_size:

  Passed to
  [`gg_forest()`](https://larry-leon.github.io/forestsearch/reference/gg_forest.md);
  defaults match the vignettes' `fp_*` values.

- widths:

  Column width vector passed to
  [`gg_forest()`](https://larry-leon.github.io/forestsearch/reference/gg_forest.md).

- clip_marker:

  Passed to
  [`gg_forest()`](https://larry-leon.github.io/forestsearch/reference/gg_forest.md);
  `"arrow"` marks whiskers clipped at the axis limits.

- xlim, ticks_at, xlab, footnote:

  `NULL` (default) selects the vignette value for the chosen
  `metric`/`panel`.

- ...:

  Further arguments merged over the defaults and passed to
  [`gg_forest()`](https://larry-leon.github.io/forestsearch/reference/gg_forest.md).

## Value

The
[`gg_forest()`](https://larry-leon.github.io/forestsearch/reference/gg_forest.md)
plot object.

## Details

Defaults reproduce the committed vignette figures exactly:
metric-specific axis limits, tick positions, reference/vertical line
styling, the four annotation columns (N plus tail probabilities and the
unconditional median), and the panel-specific footnotes. Any default can
be overridden through the named arguments or `...`, which is merged over
the assembled argument list before the
[`gg_forest()`](https://larry-leon.github.io/forestsearch/reference/gg_forest.md)
call (so e.g. `ref_col = "black"` or `widths = ...` in `...` win).

Summaries carrying effect metadata (the generic path of
[`summary.subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/summary.subgroup_sims.md),
e.g. from
[`subgroup_glm()`](https://larry-leon.github.io/forestsearch/reference/subgroup_glm.md)
fits) take a scale-aware branch instead: identity-scale measures (MD)
plot with `xlog = FALSE` and the null line at 0, ratio-scale measures
(OR) with `xlog = TRUE` and the null line at 1, and both delegate
limits/ticks to
[`gg_forest()`](https://larry-leon.github.io/forestsearch/reference/gg_forest.md)'s
data-driven defaults (the HR panel constants can hide binary medians and
remain available via `xlim=`/`ticks_at=`); annotation columns and
footnotes are built from the resolved labels and threshold strings
(legacy literals preserved when a pair equals the legacy pair), and
columns for disabled (`NA`) thresholds are dropped (the default `widths`
shrink to match; an explicit `widths` is used as given). `hr_true` is
interpreted on the fitter's estimate scale throughout. Legacy summaries
produce byte-identical
[`gg_forest()`](https://larry-leon.github.io/forestsearch/reference/gg_forest.md)
calls.

## See also

[`forest_height()`](https://larry-leon.github.io/forestsearch/reference/forest_height.md),
[`summary.subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/summary.subgroup_sims.md)

## Examples

``` r
if (FALSE) { # \dontrun{
S <- summary(sims_uniform, hr_true = 0.70)
plot(S, metric = "ub", panel = "combo")
} # }
```
