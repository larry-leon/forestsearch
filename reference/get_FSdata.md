# ForestSearch Data Preparation and Feature Selection

Prepares a dataset for ForestSearch, including options for LASSO-based
dimension reduction, GRF cuts, forced cuts, and flexible cut strategies.
Returns a list with the processed data, subgroup factor names, cut
expressions, and LASSO selection results.

## Usage

``` r
get_FSdata(
  df.analysis,
  use_lasso = FALSE,
  use_grf = FALSE,
  grf_cuts = NULL,
  confounders.name,
  cont.cutoff = 4,
  conf_force = NULL,
  conf.cont_medians = NULL,
  conf.cont_medians_force = NULL,
  replace_med_grf = TRUE,
  defaultcut_names = NULL,
  cut_type = "default",
  exclude_cuts = NULL,
  outcome.name = "tte",
  event.name = "event",
  details = TRUE
)
```

## Arguments

- df.analysis:

  Data frame containing the data.

- use_lasso:

  Logical. Whether to use LASSO for dimension reduction.

- use_grf:

  Logical. Whether to use GRF cuts.

- grf_cuts:

  Character vector of GRF cut expressions.

- confounders.name:

  Character vector of confounder variable names.

- cont.cutoff:

  Integer. Cutoff for continuous variable determination.

- conf_force:

  Character vector of forced cut expressions.

- conf.cont_medians:

  Character vector of continuous confounders to cut at median.

- conf.cont_medians_force:

  Character vector of additional continuous confounders to force median
  cut.

- replace_med_grf:

  Logical. If TRUE, removes median cuts that overlap with GRF cuts.

- defaultcut_names:

  Character vector of confounders to force default cuts.

- cut_type:

  Character. "default" or "median" for cut strategy.

- exclude_cuts:

  Character vector of cut expressions to exclude.

- outcome.name:

  Character. Name of outcome variable.

- event.name:

  Character. Name of event indicator variable.

- details:

  Logical. If TRUE, prints details during execution.

## Examples

``` r
# \donttest{
library(survival)
df <- survival::gbsg
df$grade3 <- as.integer(df$grade == "3")
fs_data <- get_FSdata(df.analysis = df,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", event.name = "status",
  use_lasso = FALSE, use_grf = FALSE)
#> # of continuous/categorical characteristics 5 2 
#> Continuous characteristics: age size nodes pgr er 
#> Categorical characteristics: meno grade3 
#> Default cuts included (1st 20) age <= mean(age) age <= median(age) age <= qlow(age) age <= qhigh(age) size <= mean(size) size <= median(size) size <= qlow(size) size <= qhigh(size) nodes <= mean(nodes) nodes <= median(nodes) nodes <= qlow(nodes) nodes <= qhigh(nodes) pgr <= mean(pgr) pgr <= median(pgr) pgr <= qlow(pgr) pgr <= qhigh(pgr) er <= mean(er) er <= median(er) er <= qlow(er) er <= qhigh(er) 
#> Categorical: meno grade3 
#> 
#> ===== CONSOLIDATED CUT EVALUATION (IMPROVED) =====
#> Evaluating 22 cut expressions once and caching...
#> Cut evaluation summary:
#>   Total cuts:  22 
#>   Valid cuts:  22 
#>   Errors:  0 
#> ✓ All 22 factors validated as 0/1
#> ===== END CONSOLIDATED CUT EVALUATION =====
#> 
#> # of candidate subgroup factors= 22 
#>  [1] "age <= 53.1"  "age <= 53"    "age <= 46"    "age <= 61"    "size <= 29.3"
#>  [6] "size <= 25"   "size <= 20"   "size <= 35"   "nodes <= 5"   "nodes <= 3"  
#> [11] "nodes <= 1"   "nodes <= 7"   "pgr <= 110"   "pgr <= 32.5"  "pgr <= 7"    
#> [16] "pgr <= 131.8" "er <= 96.3"   "er <= 36"     "er <= 8"      "er <= 114"   
#> [21] "meno"         "grade3"      
names(fs_data)
#> [1] "df"          "confs_names" "confs"       "lassokeep"   "lassoomit"  
# }
```
