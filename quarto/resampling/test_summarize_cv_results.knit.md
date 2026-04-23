---
title: "Living test: `summarize_cv_results()`"
subtitle: "Validates Phase A CV diagnostic expansion using GBSG"
author: "Larry F. León"
date: today
format:
  html:
    toc: true
    toc-depth: 3
    toc-location: left
    code-tools: true
    code-fold: true
    theme: cosmo
    fontsize: 11pt
    number-sections: true
params:
  quick_mode: true       # TRUE: sims=10 (~1-2 min). FALSE: sims=100 (~15-20 min).
---

# Purpose {#sec-purpose}

This document exercises every branch of
`summarize_cv_results()` against a `forestsearch_tenfold()` run
on the GBSG dataset.  It is a "living" test: re-running it with
an updated function should produce identical pass/fail results.
A cumulative pass/fail table is rendered in
@sec-summary.

Coverage:

1. **Minimal invocation** — `summarize_cv_results(tf)` with no
   optional arguments.  Inspects structure of the returned
   object and confirms every component slot is populated.
2. **Full invocation** — with `original_sg`,
   `original_grf_cuts`, and `create_plots = TRUE`.  Validates
   `original_agreement` and plot components.
3. **Input validation** — confirms the documented error
   conditions fire with expected messages.
4. **Backward compatibility** — simulates a pre-feature-branch
   object by dropping `grf_cuts` from `fold_summary`; confirms
   graceful NULL returns for the GRF-dependent slots.
5. **Cross-reference against manual aggregation** — reproduces
   `grf_cut_summary` using the `.clean_grf` + `.grf_dist`
   helpers lifted directly from `gbsg_10fold_dmin_sensitivity.qmd`,
   then asserts row-for-row equality with the function's output.
6. **Print method** — confirms `print.fs_cv_summary()` emits the
   expected metadata and section headers.
7. **Phase B diagnostic columns** — confirms the new
   `fold_summary` columns (`pconsistency`, `training_fs_hr`,
   `n_candidates_evaluated`) are present, typed correctly, `NA`
   exactly when their semantics require it, and in valid domains.
   Confirms the `pconsistency_distribution` and
   `fold_numeric_summary` slots are populated from the columns
   and produce the expected row counts.  Confirms
   backward-compatible behaviour when individual columns or the
   entire Phase B set are stripped.

The document uses `dmin.grf = 12, frac.tau = 0.8` (the strict
baseline from `gbsg_10fold_dmin_sensitivity.qmd` config A).  This
setting produces a diverse mix of `{er <= 0}`, `(no subgroup)`,
and non-`er` subgroups across folds, which exercises all
diagnostic branches even at small simulation counts.

**Checkpoint behaviour.** This document does not cache CV results;
the run executes every render.  At quick-mode sims, this takes on
the order of 1–2 minutes.


# Setup {#sec-setup}


::: {.cell}

```{.r .cell-code}
library(forestsearch)
library(survival)
library(data.table)
library(gt)
library(ggplot2)
library(future)
library(foreach)
library(doFuture)
```
:::



::: {.cell}

```{.r .cell-code}
# --- Simulation size -------------------------------------------------------
if (isTRUE(params$quick_mode)) {
  Ksims_cv  <- 10L
  Kfolds_cv <- 10L
  cat("Mode: quick (sims = 10).  Set quick_mode = FALSE for sims = 100.\n")
} else {
  Ksims_cv  <- 100L
  Kfolds_cv <- 10L
  cat("Mode: full (sims = 100).  Expect 15-20 min runtime.\n")
}
```

::: {.cell-output .cell-output-stdout}

```
Mode: quick (sims = 10).  Set quick_mode = FALSE for sims = 100.
```


:::

```{.r .cell-code}
# --- GBSG-specific parameters ---------------------------------------------
confounders.name <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
outcome.name     <- "time_months"
event.name       <- "status"
treat.name       <- "hormon"
id.name          <- "id"

# --- forestsearch() base arguments (mirrors sensitivity QMD config A) -----
fs_args_base <- list(
  hr.threshold              = 1.25,
  hr.consistency            = 1.00,
  pconsistency.threshold    = 0.90,
  stop_threshold            = 0.95,
  sg_focus                  = "maxSG",
  max_subgroups_search      = 10L,
  use_twostage              = FALSE,
  use_grf                   = TRUE,
  return_selected_cuts_only = TRUE,
  use_lasso                 = TRUE,
  cut_type                  = "default",
  frac.tau                  = 0.8,
  dmin.grf                  = 12,
  maxk                      = 2L,
  n.min                     = 60L,
  d0.min                    = 10L,
  d1.min                    = 10L,
  fs.splits                 = 400L,
  details                   = FALSE,
  quiet                     = TRUE,
  plot.sg                   = FALSE,
  plot.grf                  = FALSE
)

# --- Parallel configuration ------------------------------------------------
n_workers <- floor(0.95 * max(1L, parallel::detectCores() - 1L))
```
:::



::: {.cell}

```{.r .cell-code}
df              <- survival::gbsg
df$grade3       <- as.integer(df$grade == "3")
df$id           <- seq_len(nrow(df))
df$time_months  <- df$rfstime / 30.4375
```
:::



::: {.cell}

```{.r .cell-code}
# Test-result collector.  Every test case appends a row; the final
# section renders these into a gt summary table.
test_results <- data.frame(
  Section = character(0),
  Test    = character(0),
  Status  = character(0),
  Detail  = character(0),
  stringsAsFactors = FALSE
)

.tc <- function(section, test, pass, detail = "") {
  status <- if (isTRUE(pass)) "PASS" else "FAIL"
  new_row <- data.frame(
    Section = section,
    Test    = test,
    Status  = status,
    Detail  = substr(detail, 1L, 120L),
    stringsAsFactors = FALSE
  )
  # Append to parent scope
  test_results <<- rbind(test_results, new_row)
  invisible(pass)
}

.expect_error <- function(section, test, expr, pattern) {
  err_msg <- NULL
  tryCatch(
    eval(expr),
    error = function(e) { err_msg <<- conditionMessage(e) }
  )
  pass   <- !is.null(err_msg) && grepl(pattern, err_msg, fixed = FALSE)
  detail <- if (!is.null(err_msg)) err_msg else "<no error thrown>"
  .tc(section, test, pass, detail)
}
```
:::



# Baseline runs {#sec-baseline}

## Full-data ForestSearch (always runs)


::: {.cell}

```{.r .cell-code}
fs_full <- do.call(forestsearch, c(list(
  df.analysis      = df,
  confounders.name = confounders.name,
  outcome.name     = outcome.name,
  event.name       = event.name,
  treat.name       = treat.name,
  id.name          = id.name,
  parallel_args    = list(plan         = "multisession",
                          workers      = n_workers,
                          show_message = FALSE)
), fs_args_base))
plan("sequential")

cat("Full-data identified subgroup:",
    paste(fs_full$sg.harm, collapse = " & "), "\n")
```

::: {.cell-output .cell-output-stdout}

```
Full-data identified subgroup: {er <= 0} 
```


:::

```{.r .cell-code}
cat("Full-data GRF cuts:         ",
    paste(fs_full$grf_cuts, collapse = " | "), "\n")
```

::: {.cell-output .cell-output-stdout}

```
Full-data GRF cuts:          er <= 0 
```


:::
:::


## 10-fold CV


::: {.cell}

```{.r .cell-code}
.run_cv <- function() {
  t0 <- proc.time()
  out <- forestsearch_tenfold(
    fs.est        = fs_full,
    sims          = Ksims_cv,
    Kfolds        = Kfolds_cv,
    details       = FALSE,
    keep_resCV    = FALSE,
    parallel_args = list(plan         = "multisession",
                         workers      = n_workers,
                         show_message = FALSE)
  )
  plan("sequential")
  attr(out, "runtime_seconds") <-
    as.numeric((proc.time() - t0)["elapsed"])
  out
}

message(sprintf("Running %d-sim x %d-fold CV ...",
                Ksims_cv, Kfolds_cv))
tf <- .run_cv()

cat(sprintf("CV completed: %d sims x %d folds = %d sim x fold pairs\n",
            tf$sims, tf$Kfolds, nrow(tf$fold_summary)))
```

::: {.cell-output .cell-output-stdout}

```
CV completed: 10 sims x 10 folds = 100 sim x fold pairs
```


:::

```{.r .cell-code}
cat(sprintf("fold_summary columns: %s\n",
            paste(names(tf$fold_summary), collapse = ", ")))
```

::: {.cell-output .cell-output-stdout}

```
fold_summary columns: sim, fold, n_test, sg1, sg2, grf_cuts, pconsistency, training_fs_hr, n_candidates_evaluated, any_found
```


:::

```{.r .cell-code}
.tc("Baseline", "fold_summary has required columns",
    all(c("sim", "fold", "sg1", "sg2", "any_found") %in%
          names(tf$fold_summary)),
    detail = paste(names(tf$fold_summary), collapse = ", "))

.tc("Baseline", "fold_summary$grf_cuts present",
    "grf_cuts" %in% names(tf$fold_summary),
    detail = "Required for feature-branch diagnostic slots")

.tc("Baseline", "fold_summary$pconsistency present",
    "pconsistency" %in% names(tf$fold_summary),
    detail = "Required for Phase B pconsistency_distribution slot")

.tc("Baseline", "fold_summary$training_fs_hr present",
    "training_fs_hr" %in% names(tf$fold_summary),
    detail = "Phase B diagnostic column")

.tc("Baseline", "fold_summary$n_candidates_evaluated present",
    "n_candidates_evaluated" %in% names(tf$fold_summary),
    detail = "Phase B diagnostic column")
```
:::



# Test 1: Minimal invocation {#sec-test1}

Calls `summarize_cv_results(tf)` with no optional arguments.
Every slot should be populated except `original_agreement` and
`plots` (which require opt-in).


::: {.cell}

```{.r .cell-code}
cv_diag <- summarize_cv_results(tf)

# Top-level slots
slot_status <- list(
  identification_summary    = !is.null(cv_diag$identification_summary),
  grf_cut_summary           = !is.null(cv_diag$grf_cut_summary),
  cut_vs_subgroup_xtab      = !is.null(cv_diag$cut_vs_subgroup_xtab),
  no_subgroup_decomposition = !is.null(cv_diag$no_subgroup_decomposition),
  original_agreement        = is.null(cv_diag$original_agreement),
  metrics_table             = !is.null(cv_diag$metrics_table),
  plots                     = is.null(cv_diag$plots)
)

for (nm in names(slot_status)) {
  expected_populated <- !grepl("^(original_agreement|plots)$", nm)
  observed_populated <- if (grepl("^(original_agreement|plots)$", nm)) {
    !slot_status[[nm]]
  } else {
    slot_status[[nm]]
  }
  .tc("Test 1: minimal",
      sprintf("slot %s matches expectation", nm),
      isTRUE(slot_status[[nm]]),
      detail = sprintf("populated = %s", observed_populated))
}

# Metadata
.tc("Test 1: minimal", "n_sims matches tf$sims",
    cv_diag$n_sims == tf$sims,
    detail = sprintf("obs=%d exp=%d", cv_diag$n_sims, tf$sims))
.tc("Test 1: minimal", "n_folds matches tf$Kfolds",
    cv_diag$n_folds == tf$Kfolds,
    detail = sprintf("obs=%d exp=%d", cv_diag$n_folds, tf$Kfolds))
.tc("Test 1: minimal", "total_pairs matches nrow(fold_summary)",
    cv_diag$total_pairs == nrow(tf$fold_summary),
    detail = sprintf("obs=%d exp=%d",
                     cv_diag$total_pairs, nrow(tf$fold_summary)))
.tc("Test 1: minimal", "has_grf_cuts == TRUE (feature branch)",
    isTRUE(cv_diag$has_grf_cuts))
.tc("Test 1: minimal", "class == 'fs_cv_summary'",
    inherits(cv_diag, "fs_cv_summary"))

# Raw data structure
.tc("Test 1: minimal", "$data$identification is data.frame",
    is.data.frame(cv_diag$data$identification))
.tc("Test 1: minimal", "$data$grf_cuts is data.frame",
    is.data.frame(cv_diag$data$grf_cuts))
.tc("Test 1: minimal", "identification df has expected columns",
    all(c("Rank", "Subgroup", "n_folds", "pct_folds") %in%
          names(cv_diag$data$identification)))
.tc("Test 1: minimal", "grf_cuts df has expected columns",
    all(c("Rank", "grf_cut", "n_folds", "pct_folds") %in%
          names(cv_diag$data$grf_cuts)))

# Row-count sanity.  When top_n truncates, an "(N others)" row
# aggregates the tail; the n_folds column must still sum to total_pairs.
total_expected <- nrow(tf$fold_summary)
total_id       <- sum(cv_diag$data$identification$n_folds)
.tc("Test 1: minimal",
    "identification n_folds sums to total pairs (including (N others))",
    total_id == total_expected,
    detail = sprintf("obs=%d exp=%d", total_id, total_expected))

total_grf <- sum(cv_diag$data$grf_cuts$n_folds)
.tc("Test 1: minimal",
    "grf_cuts n_folds sums to total pairs (including (N others))",
    total_grf == total_expected,
    detail = sprintf("obs=%d exp=%d", total_grf, total_expected))

# If the "(N others)" row is present, it should carry the tail count.
id_has_others  <- any(grepl("^\\(\\d+ others\\)$",
                            cv_diag$data$identification$Subgroup))
grf_has_others <- any(grepl("^\\(\\d+ others\\)$",
                            cv_diag$data$grf_cuts$grf_cut))
.tc("Test 1: minimal", "identification_summary row count <= top_n + 1",
    nrow(cv_diag$data$identification) <= 16L,
    detail = sprintf("rows=%d; has (N others) row=%s",
                     nrow(cv_diag$data$identification), id_has_others))
.tc("Test 1: minimal", "grf_cut_summary row count <= top_n + 1",
    nrow(cv_diag$data$grf_cuts) <= 16L,
    detail = sprintf("rows=%d; has (N others) row=%s",
                     nrow(cv_diag$data$grf_cuts), grf_has_others))
```
:::


## Rendered gt tables from minimal invocation

### Identification summary


::: {.cell}

```{.r .cell-code}
cv_diag$identification_summary
```

::: {.cell-output-display}

```{=html}
<div id="jflwfsyzul" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#jflwfsyzul table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#jflwfsyzul thead, #jflwfsyzul tbody, #jflwfsyzul tfoot, #jflwfsyzul tr, #jflwfsyzul td, #jflwfsyzul th {
  border-style: none;
}

#jflwfsyzul p {
  margin: 0;
  padding: 0;
}

#jflwfsyzul .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 13px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#jflwfsyzul .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#jflwfsyzul .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#jflwfsyzul .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#jflwfsyzul .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#jflwfsyzul .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jflwfsyzul .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#jflwfsyzul .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#jflwfsyzul .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#jflwfsyzul .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#jflwfsyzul .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#jflwfsyzul .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#jflwfsyzul .gt_spanner_row {
  border-bottom-style: hidden;
}

#jflwfsyzul .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#jflwfsyzul .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#jflwfsyzul .gt_from_md > :first-child {
  margin-top: 0;
}

#jflwfsyzul .gt_from_md > :last-child {
  margin-bottom: 0;
}

#jflwfsyzul .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#jflwfsyzul .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#jflwfsyzul .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#jflwfsyzul .gt_row_group_first td {
  border-top-width: 2px;
}

#jflwfsyzul .gt_row_group_first th {
  border-top-width: 2px;
}

#jflwfsyzul .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#jflwfsyzul .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#jflwfsyzul .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#jflwfsyzul .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jflwfsyzul .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#jflwfsyzul .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#jflwfsyzul .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#jflwfsyzul .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#jflwfsyzul .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jflwfsyzul .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#jflwfsyzul .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#jflwfsyzul .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#jflwfsyzul .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#jflwfsyzul .gt_left {
  text-align: left;
}

#jflwfsyzul .gt_center {
  text-align: center;
}

#jflwfsyzul .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#jflwfsyzul .gt_font_normal {
  font-weight: normal;
}

#jflwfsyzul .gt_font_bold {
  font-weight: bold;
}

#jflwfsyzul .gt_font_italic {
  font-style: italic;
}

#jflwfsyzul .gt_super {
  font-size: 65%;
}

#jflwfsyzul .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#jflwfsyzul .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#jflwfsyzul .gt_indent_1 {
  text-indent: 5px;
}

#jflwfsyzul .gt_indent_2 {
  text-indent: 10px;
}

#jflwfsyzul .gt_indent_3 {
  text-indent: 15px;
}

#jflwfsyzul .gt_indent_4 {
  text-indent: 20px;
}

#jflwfsyzul .gt_indent_5 {
  text-indent: 25px;
}

#jflwfsyzul .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#jflwfsyzul div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="4" class="gt_heading gt_title gt_font_normal" style><span data-qmd-base64="KipJZGVudGlmaWVkIHN1Ymdyb3VwcyBhY3Jvc3Mgc2ltIHggZm9sZCBwYWlycyoq"><span class='gt_from_md'><strong>Identified subgroups across sim x fold pairs</strong></span></span></td>
    </tr>
    <tr class="gt_heading">
      <td colspan="4" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>10 simulations x 10 folds = 100 total pairs</td>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Rank">Rank</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Subgroup">Identified subgroup</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="n_folds">n folds</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pct_folds">% folds</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Rank" class="gt_row gt_right">1</td>
<td headers="Subgroup" class="gt_row gt_left">{er &lt;= 0}</td>
<td headers="n_folds" class="gt_row gt_right">48</td>
<td headers="pct_folds" class="gt_row gt_right">48.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">2</td>
<td headers="Subgroup" class="gt_row gt_left">(no subgroup)</td>
<td headers="n_folds" class="gt_row gt_right">21</td>
<td headers="pct_folds" class="gt_row gt_right">21.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">3</td>
<td headers="Subgroup" class="gt_row gt_left">!{age &lt;= 43} &amp; {age &lt;= 47}</td>
<td headers="n_folds" class="gt_row gt_right">10</td>
<td headers="pct_folds" class="gt_row gt_right">10.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">4</td>
<td headers="Subgroup" class="gt_row gt_left">{age &lt;= 47} &amp; !{age &lt;= 43}</td>
<td headers="n_folds" class="gt_row gt_right">5</td>
<td headers="pct_folds" class="gt_row gt_right">5.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">5</td>
<td headers="Subgroup" class="gt_row gt_left">{age &lt;= 48} &amp; !{age &lt;= 43}</td>
<td headers="n_folds" class="gt_row gt_right">5</td>
<td headers="pct_folds" class="gt_row gt_right">5.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">6</td>
<td headers="Subgroup" class="gt_row gt_left">!{age &lt;= 43} &amp; {age &lt;= 48}</td>
<td headers="n_folds" class="gt_row gt_right">3</td>
<td headers="pct_folds" class="gt_row gt_right">3.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">7</td>
<td headers="Subgroup" class="gt_row gt_left">{age &lt;= 50} &amp; !{age &lt;= 43}</td>
<td headers="n_folds" class="gt_row gt_right">2</td>
<td headers="pct_folds" class="gt_row gt_right">2.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">8</td>
<td headers="Subgroup" class="gt_row gt_left">!{size &lt;= 25} &amp; {er &lt;= 8}</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">9</td>
<td headers="Subgroup" class="gt_row gt_left">{er &lt;= 0} &amp; {pgr &lt;= 32}</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">10</td>
<td headers="Subgroup" class="gt_row gt_left">{er &lt;= 0} &amp; {pgr &lt;= 34.5}</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">11</td>
<td headers="Subgroup" class="gt_row gt_left">{er &lt;= 0} &amp; {pgr &lt;= 34}</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">12</td>
<td headers="Subgroup" class="gt_row gt_left">{er &lt;= 8} &amp; !{size &lt;= 25}</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">13</td>
<td headers="Subgroup" class="gt_row gt_left">{grade3} &amp; {pgr &lt;= 8}</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
  </tbody>
  
</table>
</div>
```

:::
:::


### GRF cut summary


::: {.cell}

```{.r .cell-code}
cv_diag$grf_cut_summary
```

::: {.cell-output-display}

```{=html}
<div id="jcmpahvikw" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#jcmpahvikw table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#jcmpahvikw thead, #jcmpahvikw tbody, #jcmpahvikw tfoot, #jcmpahvikw tr, #jcmpahvikw td, #jcmpahvikw th {
  border-style: none;
}

#jcmpahvikw p {
  margin: 0;
  padding: 0;
}

#jcmpahvikw .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 13px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#jcmpahvikw .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#jcmpahvikw .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#jcmpahvikw .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#jcmpahvikw .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#jcmpahvikw .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jcmpahvikw .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#jcmpahvikw .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#jcmpahvikw .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#jcmpahvikw .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#jcmpahvikw .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#jcmpahvikw .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#jcmpahvikw .gt_spanner_row {
  border-bottom-style: hidden;
}

#jcmpahvikw .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#jcmpahvikw .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#jcmpahvikw .gt_from_md > :first-child {
  margin-top: 0;
}

#jcmpahvikw .gt_from_md > :last-child {
  margin-bottom: 0;
}

#jcmpahvikw .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#jcmpahvikw .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#jcmpahvikw .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#jcmpahvikw .gt_row_group_first td {
  border-top-width: 2px;
}

#jcmpahvikw .gt_row_group_first th {
  border-top-width: 2px;
}

#jcmpahvikw .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#jcmpahvikw .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#jcmpahvikw .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#jcmpahvikw .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jcmpahvikw .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#jcmpahvikw .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#jcmpahvikw .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#jcmpahvikw .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#jcmpahvikw .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jcmpahvikw .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#jcmpahvikw .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#jcmpahvikw .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#jcmpahvikw .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#jcmpahvikw .gt_left {
  text-align: left;
}

#jcmpahvikw .gt_center {
  text-align: center;
}

#jcmpahvikw .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#jcmpahvikw .gt_font_normal {
  font-weight: normal;
}

#jcmpahvikw .gt_font_bold {
  font-weight: bold;
}

#jcmpahvikw .gt_font_italic {
  font-style: italic;
}

#jcmpahvikw .gt_super {
  font-size: 65%;
}

#jcmpahvikw .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#jcmpahvikw .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#jcmpahvikw .gt_indent_1 {
  text-indent: 5px;
}

#jcmpahvikw .gt_indent_2 {
  text-indent: 10px;
}

#jcmpahvikw .gt_indent_3 {
  text-indent: 15px;
}

#jcmpahvikw .gt_indent_4 {
  text-indent: 20px;
}

#jcmpahvikw .gt_indent_5 {
  text-indent: 25px;
}

#jcmpahvikw .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#jcmpahvikw div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="4" class="gt_heading gt_title gt_font_normal" style><span data-qmd-base64="KipHUkYgY3V0cyBhY3Jvc3MgdHJhaW5pbmcgc2V0cyoq"><span class='gt_from_md'><strong>GRF cuts across training sets</strong></span></span></td>
    </tr>
    <tr class="gt_heading">
      <td colspan="4" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>Raw grf_cuts strings; 10 simulations x 10 folds = 100 total pairs</td>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Rank">Rank</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="grf_cut">GRF cuts</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="n_folds">n folds</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pct_folds">% folds</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Rank" class="gt_row gt_right">1</td>
<td headers="grf_cut" class="gt_row gt_left">er &lt;= 0</td>
<td headers="n_folds" class="gt_row gt_right">41</td>
<td headers="pct_folds" class="gt_row gt_right">41.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">2</td>
<td headers="grf_cut" class="gt_row gt_left">age &lt;= 48 | age &lt;= 43 | er &lt;= 0</td>
<td headers="n_folds" class="gt_row gt_right">15</td>
<td headers="pct_folds" class="gt_row gt_right">15.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">3</td>
<td headers="grf_cut" class="gt_row gt_left">age &lt;= 47 | age &lt;= 43 | er &lt;= 0</td>
<td headers="n_folds" class="gt_row gt_right">12</td>
<td headers="pct_folds" class="gt_row gt_right">12.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">4</td>
<td headers="grf_cut" class="gt_row gt_left">(no GRF cut)</td>
<td headers="n_folds" class="gt_row gt_right">10</td>
<td headers="pct_folds" class="gt_row gt_right">10.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">5</td>
<td headers="grf_cut" class="gt_row gt_left">age &lt;= 50 | age &lt;= 43 | er &lt;= 0</td>
<td headers="n_folds" class="gt_row gt_right">7</td>
<td headers="pct_folds" class="gt_row gt_right">7.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">6</td>
<td headers="grf_cut" class="gt_row gt_left">age &lt;= 47 | age &lt;= 43 | er &lt;= 369</td>
<td headers="n_folds" class="gt_row gt_right">2</td>
<td headers="pct_folds" class="gt_row gt_right">2.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">7</td>
<td headers="grf_cut" class="gt_row gt_left">age &lt;= 48 | age &lt;= 44 | er &lt;= 0</td>
<td headers="n_folds" class="gt_row gt_right">2</td>
<td headers="pct_folds" class="gt_row gt_right">2.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">8</td>
<td headers="grf_cut" class="gt_row gt_left">age &lt;= 48 | age &lt;= 45 | er &lt;= 0</td>
<td headers="n_folds" class="gt_row gt_right">2</td>
<td headers="pct_folds" class="gt_row gt_right">2.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">9</td>
<td headers="grf_cut" class="gt_row gt_left">age &lt;= 50 | age &lt;= 44 | er &lt;= 0</td>
<td headers="n_folds" class="gt_row gt_right">2</td>
<td headers="pct_folds" class="gt_row gt_right">2.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">10</td>
<td headers="grf_cut" class="gt_row gt_left">age &lt;= 47 | age &lt;= 43 | size &lt;= 14</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">11</td>
<td headers="grf_cut" class="gt_row gt_left">age &lt;= 47 | age &lt;= 44 | er &lt;= 0</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">12</td>
<td headers="grf_cut" class="gt_row gt_left">age &lt;= 48 | age &lt;= 43 | er &lt;= 369</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">13</td>
<td headers="grf_cut" class="gt_row gt_left">age &lt;= 48 | age &lt;= 43 | size &lt;= 19</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">14</td>
<td headers="grf_cut" class="gt_row gt_left">age &lt;= 49 | age &lt;= 43 | er &lt;= 0</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">15</td>
<td headers="grf_cut" class="gt_row gt_left">age &lt;= 49 | age &lt;= 45 | er &lt;= 317</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
    <tr><td headers="Rank" class="gt_row gt_right">16</td>
<td headers="grf_cut" class="gt_row gt_left">(1 others)</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
  </tbody>
  
</table>
</div>
```

:::
:::


### Cut × subgroup cross-tab


::: {.cell}

```{.r .cell-code}
cv_diag$cut_vs_subgroup_xtab
```

::: {.cell-output-display}

```{=html}
<div id="bbdjwuojbk" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#bbdjwuojbk table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#bbdjwuojbk thead, #bbdjwuojbk tbody, #bbdjwuojbk tfoot, #bbdjwuojbk tr, #bbdjwuojbk td, #bbdjwuojbk th {
  border-style: none;
}

#bbdjwuojbk p {
  margin: 0;
  padding: 0;
}

#bbdjwuojbk .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 13px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#bbdjwuojbk .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#bbdjwuojbk .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#bbdjwuojbk .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#bbdjwuojbk .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#bbdjwuojbk .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bbdjwuojbk .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#bbdjwuojbk .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#bbdjwuojbk .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#bbdjwuojbk .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#bbdjwuojbk .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#bbdjwuojbk .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#bbdjwuojbk .gt_spanner_row {
  border-bottom-style: hidden;
}

#bbdjwuojbk .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#bbdjwuojbk .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#bbdjwuojbk .gt_from_md > :first-child {
  margin-top: 0;
}

#bbdjwuojbk .gt_from_md > :last-child {
  margin-bottom: 0;
}

#bbdjwuojbk .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#bbdjwuojbk .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#bbdjwuojbk .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#bbdjwuojbk .gt_row_group_first td {
  border-top-width: 2px;
}

#bbdjwuojbk .gt_row_group_first th {
  border-top-width: 2px;
}

#bbdjwuojbk .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#bbdjwuojbk .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#bbdjwuojbk .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#bbdjwuojbk .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bbdjwuojbk .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#bbdjwuojbk .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#bbdjwuojbk .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#bbdjwuojbk .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#bbdjwuojbk .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bbdjwuojbk .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#bbdjwuojbk .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#bbdjwuojbk .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#bbdjwuojbk .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#bbdjwuojbk .gt_left {
  text-align: left;
}

#bbdjwuojbk .gt_center {
  text-align: center;
}

#bbdjwuojbk .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#bbdjwuojbk .gt_font_normal {
  font-weight: normal;
}

#bbdjwuojbk .gt_font_bold {
  font-weight: bold;
}

#bbdjwuojbk .gt_font_italic {
  font-style: italic;
}

#bbdjwuojbk .gt_super {
  font-size: 65%;
}

#bbdjwuojbk .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#bbdjwuojbk .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#bbdjwuojbk .gt_indent_1 {
  text-indent: 5px;
}

#bbdjwuojbk .gt_indent_2 {
  text-indent: 10px;
}

#bbdjwuojbk .gt_indent_3 {
  text-indent: 15px;
}

#bbdjwuojbk .gt_indent_4 {
  text-indent: 20px;
}

#bbdjwuojbk .gt_indent_5 {
  text-indent: 25px;
}

#bbdjwuojbk .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#bbdjwuojbk div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="4" class="gt_heading gt_title gt_font_normal" style><span data-qmd-base64="KipHUkYgY3V0IHggaWRlbnRpZmllZCBzdWJncm91cCoq"><span class='gt_from_md'><strong>GRF cut x identified subgroup</strong></span></span></td>
    </tr>
    <tr class="gt_heading">
      <td colspan="4" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>Cross-tabulation across 100 sim x fold pairs</td>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="grf_cut">GRF cuts</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="subgroup">Identified subgroup</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="n_folds">n folds</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pct_folds">% folds</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="grf_cut" class="gt_row gt_left">er &lt;= 0</td>
<td headers="subgroup" class="gt_row gt_left">{er &lt;= 0}</td>
<td headers="n_folds" class="gt_row gt_right">35</td>
<td headers="pct_folds" class="gt_row gt_right">35.00</td></tr>
    <tr><td headers="grf_cut" class="gt_row gt_left">age &lt;= 47 | age &lt;= 43 | er &lt;= 0</td>
<td headers="subgroup" class="gt_row gt_left">!{age &lt;= 43} &amp; {age &lt;= 47}</td>
<td headers="n_folds" class="gt_row gt_right">9</td>
<td headers="pct_folds" class="gt_row gt_right">9.00</td></tr>
    <tr><td headers="grf_cut" class="gt_row gt_left">(no GRF cut)</td>
<td headers="subgroup" class="gt_row gt_left">(no subgroup)</td>
<td headers="n_folds" class="gt_row gt_right">9</td>
<td headers="pct_folds" class="gt_row gt_right">9.00</td></tr>
    <tr><td headers="grf_cut" class="gt_row gt_left">age &lt;= 48 | age &lt;= 43 | er &lt;= 0</td>
<td headers="subgroup" class="gt_row gt_left">{er &lt;= 0}</td>
<td headers="n_folds" class="gt_row gt_right">5</td>
<td headers="pct_folds" class="gt_row gt_right">5.00</td></tr>
    <tr><td headers="grf_cut" class="gt_row gt_left">er &lt;= 0</td>
<td headers="subgroup" class="gt_row gt_left">(no subgroup)</td>
<td headers="n_folds" class="gt_row gt_right">4</td>
<td headers="pct_folds" class="gt_row gt_right">4.00</td></tr>
    <tr><td headers="grf_cut" class="gt_row gt_left">age &lt;= 48 | age &lt;= 43 | er &lt;= 0</td>
<td headers="subgroup" class="gt_row gt_left">{age &lt;= 48} &amp; !{age &lt;= 43}</td>
<td headers="n_folds" class="gt_row gt_right">4</td>
<td headers="pct_folds" class="gt_row gt_right">4.00</td></tr>
    <tr><td headers="grf_cut" class="gt_row gt_left">age &lt;= 48 | age &lt;= 43 | er &lt;= 0</td>
<td headers="subgroup" class="gt_row gt_left">!{age &lt;= 43} &amp; {age &lt;= 48}</td>
<td headers="n_folds" class="gt_row gt_right">3</td>
<td headers="pct_folds" class="gt_row gt_right">3.00</td></tr>
    <tr><td headers="grf_cut" class="gt_row gt_left">age &lt;= 47 | age &lt;= 43 | er &lt;= 0</td>
<td headers="subgroup" class="gt_row gt_left">{age &lt;= 47} &amp; !{age &lt;= 43}</td>
<td headers="n_folds" class="gt_row gt_right">3</td>
<td headers="pct_folds" class="gt_row gt_right">3.00</td></tr>
    <tr><td headers="grf_cut" class="gt_row gt_left">age &lt;= 48 | age &lt;= 43 | er &lt;= 0</td>
<td headers="subgroup" class="gt_row gt_left">(no subgroup)</td>
<td headers="n_folds" class="gt_row gt_right">2</td>
<td headers="pct_folds" class="gt_row gt_right">2.00</td></tr>
    <tr><td headers="grf_cut" class="gt_row gt_left">age &lt;= 50 | age &lt;= 43 | er &lt;= 0</td>
<td headers="subgroup" class="gt_row gt_left">{age &lt;= 50} &amp; !{age &lt;= 43}</td>
<td headers="n_folds" class="gt_row gt_right">2</td>
<td headers="pct_folds" class="gt_row gt_right">2.00</td></tr>
    <tr><td headers="grf_cut" class="gt_row gt_left">age &lt;= 50 | age &lt;= 43 | er &lt;= 0</td>
<td headers="subgroup" class="gt_row gt_left">{er &lt;= 0}</td>
<td headers="n_folds" class="gt_row gt_right">2</td>
<td headers="pct_folds" class="gt_row gt_right">2.00</td></tr>
    <tr><td headers="grf_cut" class="gt_row gt_left">age &lt;= 47 | age &lt;= 43 | er &lt;= 369</td>
<td headers="subgroup" class="gt_row gt_left">!{age &lt;= 43} &amp; {age &lt;= 47}</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
    <tr><td headers="grf_cut" class="gt_row gt_left">er &lt;= 0</td>
<td headers="subgroup" class="gt_row gt_left">!{size &lt;= 25} &amp; {er &lt;= 8}</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
    <tr><td headers="grf_cut" class="gt_row gt_left">age &lt;= 48 | age &lt;= 43 | size &lt;= 19</td>
<td headers="subgroup" class="gt_row gt_left">(no subgroup)</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
    <tr><td headers="grf_cut" class="gt_row gt_left">age &lt;= 48 | age &lt;= 44 | er &lt;= 0</td>
<td headers="subgroup" class="gt_row gt_left">(no subgroup)</td>
<td headers="n_folds" class="gt_row gt_right">1</td>
<td headers="pct_folds" class="gt_row gt_right">1.00</td></tr>
    <tr><td headers="grf_cut" class="gt_row gt_left">(18 others)</td>
<td headers="subgroup" class="gt_row gt_left">(18 others)</td>
<td headers="n_folds" class="gt_row gt_right">18</td>
<td headers="pct_folds" class="gt_row gt_right">18.00</td></tr>
  </tbody>
  <tfoot>
    <tr class="gt_sourcenotes">
      <td class="gt_sourcenote" colspan="4">Rows where the identified subgroup does not align with GRF's output reveal divergence between candidate discovery and consistency evaluation.</td>
    </tr>
  </tfoot>
</table>
</div>
```

:::
:::


### No-subgroup decomposition


::: {.cell}

```{.r .cell-code}
cv_diag$no_subgroup_decomposition
```

::: {.cell-output-display}

```{=html}
<div id="nmlhnqtnuu" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#nmlhnqtnuu table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#nmlhnqtnuu thead, #nmlhnqtnuu tbody, #nmlhnqtnuu tfoot, #nmlhnqtnuu tr, #nmlhnqtnuu td, #nmlhnqtnuu th {
  border-style: none;
}

#nmlhnqtnuu p {
  margin: 0;
  padding: 0;
}

#nmlhnqtnuu .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 13px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#nmlhnqtnuu .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#nmlhnqtnuu .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#nmlhnqtnuu .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#nmlhnqtnuu .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#nmlhnqtnuu .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#nmlhnqtnuu .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#nmlhnqtnuu .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#nmlhnqtnuu .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#nmlhnqtnuu .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#nmlhnqtnuu .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#nmlhnqtnuu .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#nmlhnqtnuu .gt_spanner_row {
  border-bottom-style: hidden;
}

#nmlhnqtnuu .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#nmlhnqtnuu .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#nmlhnqtnuu .gt_from_md > :first-child {
  margin-top: 0;
}

#nmlhnqtnuu .gt_from_md > :last-child {
  margin-bottom: 0;
}

#nmlhnqtnuu .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#nmlhnqtnuu .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#nmlhnqtnuu .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#nmlhnqtnuu .gt_row_group_first td {
  border-top-width: 2px;
}

#nmlhnqtnuu .gt_row_group_first th {
  border-top-width: 2px;
}

#nmlhnqtnuu .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#nmlhnqtnuu .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#nmlhnqtnuu .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#nmlhnqtnuu .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#nmlhnqtnuu .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#nmlhnqtnuu .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#nmlhnqtnuu .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#nmlhnqtnuu .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#nmlhnqtnuu .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#nmlhnqtnuu .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#nmlhnqtnuu .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#nmlhnqtnuu .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#nmlhnqtnuu .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#nmlhnqtnuu .gt_left {
  text-align: left;
}

#nmlhnqtnuu .gt_center {
  text-align: center;
}

#nmlhnqtnuu .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#nmlhnqtnuu .gt_font_normal {
  font-weight: normal;
}

#nmlhnqtnuu .gt_font_bold {
  font-weight: bold;
}

#nmlhnqtnuu .gt_font_italic {
  font-style: italic;
}

#nmlhnqtnuu .gt_super {
  font-size: 65%;
}

#nmlhnqtnuu .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#nmlhnqtnuu .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#nmlhnqtnuu .gt_indent_1 {
  text-indent: 5px;
}

#nmlhnqtnuu .gt_indent_2 {
  text-indent: 10px;
}

#nmlhnqtnuu .gt_indent_3 {
  text-indent: 15px;
}

#nmlhnqtnuu .gt_indent_4 {
  text-indent: 20px;
}

#nmlhnqtnuu .gt_indent_5 {
  text-indent: 25px;
}

#nmlhnqtnuu .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#nmlhnqtnuu div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="4" class="gt_heading gt_title gt_font_normal" style><span data-qmd-base64="KipEZWNvbXBvc2l0aW9uIG9mIG5vLXN1Ymdyb3VwIGZvbGRzKio="><span class='gt_from_md'><strong>Decomposition of no-subgroup folds</strong></span></span></td>
    </tr>
    <tr class="gt_heading">
      <td colspan="4" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>Out of 100 sim x fold pairs</td>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Category">Category</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="n_folds">n folds</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pct_of_total">% of all pairs</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pct_of_no_sg">% of no-sg folds</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Category" class="gt_row gt_left">Folds with no subgroup identified</td>
<td headers="n_folds" class="gt_row gt_right">21</td>
<td headers="pct_of_total" class="gt_row gt_right">21.00</td>
<td headers="pct_of_no_sg" class="gt_row gt_right">100.00</td></tr>
    <tr><td headers="Category" class="gt_row gt_left">  GRF returned no cut</td>
<td headers="n_folds" class="gt_row gt_right">9</td>
<td headers="pct_of_total" class="gt_row gt_right">9.00</td>
<td headers="pct_of_no_sg" class="gt_row gt_right">42.86</td></tr>
    <tr><td headers="Category" class="gt_row gt_left">  GRF returned a cut; consistency rejected</td>
<td headers="n_folds" class="gt_row gt_right">12</td>
<td headers="pct_of_total" class="gt_row gt_right">12.00</td>
<td headers="pct_of_no_sg" class="gt_row gt_right">57.14</td></tr>
  </tbody>
  <tfoot>
    <tr class="gt_sourcenotes">
      <td class="gt_sourcenote" colspan="4">Guidance: a high 'GRF returned no cut' share argues for relaxing dmin.grf; a high 'consistency rejected' share argues for relaxing pconsistency.threshold.</td>
    </tr>
  </tfoot>
</table>
</div>
```

:::
:::


### Metrics table (delegated to `cv_metrics_tables()`)


::: {.cell}

```{.r .cell-code}
cv_diag$metrics_table
```

::: {.cell-output-display}

```{=html}
<div id="ihiocftcls" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#ihiocftcls table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#ihiocftcls thead, #ihiocftcls tbody, #ihiocftcls tfoot, #ihiocftcls tr, #ihiocftcls td, #ihiocftcls th {
  border-style: none;
}

#ihiocftcls p {
  margin: 0;
  padding: 0;
}

#ihiocftcls .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#ihiocftcls .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#ihiocftcls .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#ihiocftcls .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#ihiocftcls .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#ihiocftcls .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ihiocftcls .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#ihiocftcls .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#ihiocftcls .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#ihiocftcls .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#ihiocftcls .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#ihiocftcls .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#ihiocftcls .gt_spanner_row {
  border-bottom-style: hidden;
}

#ihiocftcls .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#ihiocftcls .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#ihiocftcls .gt_from_md > :first-child {
  margin-top: 0;
}

#ihiocftcls .gt_from_md > :last-child {
  margin-bottom: 0;
}

#ihiocftcls .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#ihiocftcls .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#ihiocftcls .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#ihiocftcls .gt_row_group_first td {
  border-top-width: 2px;
}

#ihiocftcls .gt_row_group_first th {
  border-top-width: 2px;
}

#ihiocftcls .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ihiocftcls .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#ihiocftcls .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#ihiocftcls .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ihiocftcls .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ihiocftcls .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#ihiocftcls .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#ihiocftcls .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#ihiocftcls .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ihiocftcls .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#ihiocftcls .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#ihiocftcls .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#ihiocftcls .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#ihiocftcls .gt_left {
  text-align: left;
}

#ihiocftcls .gt_center {
  text-align: center;
}

#ihiocftcls .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#ihiocftcls .gt_font_normal {
  font-weight: normal;
}

#ihiocftcls .gt_font_bold {
  font-weight: bold;
}

#ihiocftcls .gt_font_italic {
  font-style: italic;
}

#ihiocftcls .gt_super {
  font-size: 65%;
}

#ihiocftcls .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#ihiocftcls .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#ihiocftcls .gt_indent_1 {
  text-indent: 5px;
}

#ihiocftcls .gt_indent_2 {
  text-indent: 10px;
}

#ihiocftcls .gt_indent_3 {
  text-indent: 15px;
}

#ihiocftcls .gt_indent_4 {
  text-indent: 20px;
}

#ihiocftcls .gt_indent_5 {
  text-indent: 25px;
}

#ihiocftcls .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#ihiocftcls div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="3" class="gt_heading gt_title gt_font_normal" style><span data-qmd-base64="KipDcm9zcy1WYWxpZGF0aW9uIE1ldHJpY3MqKg=="><span class='gt_from_md'><strong>Cross-Validation Metrics</strong></span></span></td>
    </tr>
    <tr class="gt_heading">
      <td colspan="3" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>Subgroup: Identified Subgroup</td>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" style="font-weight: bold;" scope="col" id="Metric">Metric</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" style="font-weight: bold;" scope="col" id="Description">Description</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" style="font-weight: bold;" scope="col" id="Value">Value (%)</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" style="font-weight: bold;" scope="colgroup" id="Agreement">Agreement</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Agreement  Metric" class="gt_row gt_left">Sensitivity (H)</td>
<td headers="Agreement  Description" class="gt_row gt_left">Agreement rate for subgroup H</td>
<td headers="Agreement  Value" class="gt_row gt_right">53.7</td></tr>
    <tr><td headers="Agreement  Metric" class="gt_row gt_left">Sensitivity (Hc)</td>
<td headers="Agreement  Description" class="gt_row gt_left">Agreement rate for complement Hc</td>
<td headers="Agreement  Value" class="gt_row gt_right">95.0</td></tr>
    <tr><td headers="Agreement  Metric" class="gt_row gt_left">PPV (H)</td>
<td headers="Agreement  Description" class="gt_row gt_left">Positive predictive value for H</td>
<td headers="Agreement  Value" class="gt_row gt_right">58.1</td></tr>
    <tr><td headers="Agreement  Metric" class="gt_row gt_left">PPV (Hc)</td>
<td headers="Agreement  Description" class="gt_row gt_left">Positive predictive value for Hc</td>
<td headers="Agreement  Value" class="gt_row gt_right">93.8</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" style="font-weight: bold;" scope="colgroup" id="Subgroup Finding">Subgroup Finding</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Subgroup Finding  Metric" class="gt_row gt_left">Any Found</td>
<td headers="Subgroup Finding  Description" class="gt_row gt_left">Any subgroup identified</td>
<td headers="Subgroup Finding  Value" class="gt_row gt_right">80.0</td></tr>
    <tr><td headers="Subgroup Finding  Metric" class="gt_row gt_left">Exact Match</td>
<td headers="Subgroup Finding  Description" class="gt_row gt_left">Exact match on all factors</td>
<td headers="Subgroup Finding  Value" class="gt_row gt_right">50.0</td></tr>
    <tr><td headers="Subgroup Finding  Metric" class="gt_row gt_left">At Least 1</td>
<td headers="Subgroup Finding  Description" class="gt_row gt_left">At least one factor matches</td>
<td headers="Subgroup Finding  Value" class="gt_row gt_right">50.0</td></tr>
    <tr><td headers="Subgroup Finding  Metric" class="gt_row gt_left">Cov1 Any</td>
<td headers="Subgroup Finding  Description" class="gt_row gt_left">First covariate found (any cut)</td>
<td headers="Subgroup Finding  Value" class="gt_row gt_right">50.0</td></tr>
    <tr><td headers="Subgroup Finding  Metric" class="gt_row gt_left">Cov2 Any</td>
<td headers="Subgroup Finding  Description" class="gt_row gt_left">Second covariate found (any cut)</td>
<td headers="Subgroup Finding  Value" class="gt_row gt_right">NA</td></tr>
    <tr><td headers="Subgroup Finding  Metric" class="gt_row gt_left">Cov1 &amp; Cov2</td>
<td headers="Subgroup Finding  Description" class="gt_row gt_left">Both covariates found</td>
<td headers="Subgroup Finding  Value" class="gt_row gt_right">0.0</td></tr>
    <tr><td headers="Subgroup Finding  Metric" class="gt_row gt_left">Cov1 Exact</td>
<td headers="Subgroup Finding  Description" class="gt_row gt_left">First covariate exact match</td>
<td headers="Subgroup Finding  Value" class="gt_row gt_right">50.0</td></tr>
    <tr><td headers="Subgroup Finding  Metric" class="gt_row gt_left">Cov2 Exact</td>
<td headers="Subgroup Finding  Description" class="gt_row gt_left">Second covariate exact match</td>
<td headers="Subgroup Finding  Value" class="gt_row gt_right">NA</td></tr>
  </tbody>
  <tfoot>
    <tr class="gt_footnotes">
      <td class="gt_footnote" colspan="3"> Based on 10 simulation(s) with 10-fold CV. Values are medians shown as percentages.</td>
    </tr>
  </tfoot>
</table>
</div>
```

:::
:::



# Test 2: Full invocation with originals and plots {#sec-test2}


::: {.cell}

```{.r .cell-code}
cv_diag_full <- summarize_cv_results(
  cv_output         = tf,
  original_sg       = fs_full$sg.harm,
  original_grf_cuts = fs_full$grf_cuts,
  create_plots      = TRUE
)

.tc("Test 2: full", "$original_agreement populated",
    !is.null(cv_diag_full$original_agreement))
.tc("Test 2: full", "$plots populated",
    !is.null(cv_diag_full$plots))
.tc("Test 2: full", "$plots$identification is ggplot",
    inherits(cv_diag_full$plots$identification, "ggplot"))
.tc("Test 2: full", "$plots$grf_cut is ggplot (has_grf_cuts TRUE)",
    inherits(cv_diag_full$plots$grf_cut, "ggplot"))

# original_agreement raw data structure
oa_df <- cv_diag_full$data$original_agreement
.tc("Test 2: full", "original_agreement df has Metric, Value",
    all(c("Metric", "Value") %in% names(oa_df)))
.tc("Test 2: full", "original_agreement includes Exact subgroup match row",
    any(grepl("Exact subgroup match", oa_df$Metric)))
.tc("Test 2: full", "original_agreement includes Exact GRF-cut match row",
    any(grepl("Exact GRF-cut match", oa_df$Metric)))
.tc("Test 2: full", "original_agreement includes Both match row",
    any(grepl("Both subgroup and GRF cut match", oa_df$Metric)))
```
:::


## Rendered original agreement table


::: {.cell}

```{.r .cell-code}
cv_diag_full$original_agreement
```

::: {.cell-output-display}

```{=html}
<div id="rqjcseknza" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#rqjcseknza table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#rqjcseknza thead, #rqjcseknza tbody, #rqjcseknza tfoot, #rqjcseknza tr, #rqjcseknza td, #rqjcseknza th {
  border-style: none;
}

#rqjcseknza p {
  margin: 0;
  padding: 0;
}

#rqjcseknza .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 13px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#rqjcseknza .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#rqjcseknza .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#rqjcseknza .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#rqjcseknza .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#rqjcseknza .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rqjcseknza .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#rqjcseknza .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#rqjcseknza .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#rqjcseknza .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#rqjcseknza .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#rqjcseknza .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#rqjcseknza .gt_spanner_row {
  border-bottom-style: hidden;
}

#rqjcseknza .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#rqjcseknza .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#rqjcseknza .gt_from_md > :first-child {
  margin-top: 0;
}

#rqjcseknza .gt_from_md > :last-child {
  margin-bottom: 0;
}

#rqjcseknza .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#rqjcseknza .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#rqjcseknza .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#rqjcseknza .gt_row_group_first td {
  border-top-width: 2px;
}

#rqjcseknza .gt_row_group_first th {
  border-top-width: 2px;
}

#rqjcseknza .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#rqjcseknza .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#rqjcseknza .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#rqjcseknza .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rqjcseknza .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#rqjcseknza .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#rqjcseknza .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#rqjcseknza .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#rqjcseknza .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rqjcseknza .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#rqjcseknza .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#rqjcseknza .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#rqjcseknza .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#rqjcseknza .gt_left {
  text-align: left;
}

#rqjcseknza .gt_center {
  text-align: center;
}

#rqjcseknza .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#rqjcseknza .gt_font_normal {
  font-weight: normal;
}

#rqjcseknza .gt_font_bold {
  font-weight: bold;
}

#rqjcseknza .gt_font_italic {
  font-style: italic;
}

#rqjcseknza .gt_super {
  font-size: 65%;
}

#rqjcseknza .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#rqjcseknza .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#rqjcseknza .gt_indent_1 {
  text-indent: 5px;
}

#rqjcseknza .gt_indent_2 {
  text-indent: 10px;
}

#rqjcseknza .gt_indent_3 {
  text-indent: 15px;
}

#rqjcseknza .gt_indent_4 {
  text-indent: 20px;
}

#rqjcseknza .gt_indent_5 {
  text-indent: 25px;
}

#rqjcseknza .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#rqjcseknza div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="2" class="gt_heading gt_title gt_font_normal" style><span data-qmd-base64="KipBZ3JlZW1lbnQgd2l0aCBvcmlnaW5hbCBmdWxsLWRhdGEgYW5hbHlzaXMqKg=="><span class='gt_from_md'><strong>Agreement with original full-data analysis</strong></span></span></td>
    </tr>
    <tr class="gt_heading">
      <td colspan="2" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>Fraction of sim x fold pairs matching the original result</td>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Metric"></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Value">Value</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Metric" class="gt_row gt_left">Total sim x fold pairs</td>
<td headers="Value" class="gt_row gt_left">100</td></tr>
    <tr><td headers="Metric" class="gt_row gt_left">Original subgroup definition</td>
<td headers="Value" class="gt_row gt_left">{er &lt;= 0}</td></tr>
    <tr><td headers="Metric" class="gt_row gt_left" style="background-color: #E8F4F8; font-weight: bold;">Exact subgroup match</td>
<td headers="Value" class="gt_row gt_left">48 (48.0%)</td></tr>
    <tr><td headers="Metric" class="gt_row gt_left" style="background-color: #E8F4F8; font-weight: bold;">Partial match (&gt;=1 shared factor)</td>
<td headers="Value" class="gt_row gt_left">51 (51.0%)</td></tr>
    <tr><td headers="Metric" class="gt_row gt_left">Original GRF cut(s)</td>
<td headers="Value" class="gt_row gt_left">er &lt;= 0</td></tr>
    <tr><td headers="Metric" class="gt_row gt_left" style="background-color: #E8F4F8; font-weight: bold;">Exact GRF-cut match</td>
<td headers="Value" class="gt_row gt_left">41 (41.0%)</td></tr>
    <tr><td headers="Metric" class="gt_row gt_left" style="background-color: #E8F4F8; font-weight: bold;">Both subgroup and GRF cut match</td>
<td headers="Value" class="gt_row gt_left">35 (35.0%)</td></tr>
  </tbody>
  
</table>
</div>
```

:::
:::


## Rendered plots


::: {.cell}

```{.r .cell-code}
cv_diag_full$plots$identification
```

::: {.cell-output-display}
![](test_summarize_cv_results_files/figure-html/test2-display-plot-id-1.png){width=864}
:::
:::



::: {.cell}

```{.r .cell-code}
cv_diag_full$plots$grf_cut
```

::: {.cell-output-display}
![](test_summarize_cv_results_files/figure-html/test2-display-plot-grf-1.png){width=864}
:::
:::



# Test 3: Input validation {#sec-test3}

Confirms each documented error condition fires with the
expected message pattern.


::: {.cell}

```{.r .cell-code}
# Case 3.1: non-list input
.expect_error(
  "Test 3: errors", "non-list input rejected",
  expr    = quote(summarize_cv_results(cv_output = "not a list")),
  pattern = "must be a list"
)

# Case 3.2: empty fold_summary
tf_empty <- tf
tf_empty$fold_summary <- tf$fold_summary[0, , drop = FALSE]
.expect_error(
  "Test 3: errors", "empty fold_summary rejected",
  expr    = quote(summarize_cv_results(cv_output = tf_empty)),
  pattern = "missing or empty"
)

# Case 3.3: missing required column
tf_bad <- tf
tf_bad$fold_summary$sg1 <- NULL
.expect_error(
  "Test 3: errors", "missing required column rejected",
  expr    = quote(summarize_cv_results(cv_output = tf_bad)),
  pattern = "missing required column"
)

# Case 3.4: bad top_n
.expect_error(
  "Test 3: errors", "non-positive top_n rejected",
  expr    = quote(summarize_cv_results(cv_output = tf, top_n = 0)),
  pattern = "positive integer"
)

# Case 3.5: fs_kfold redirection
# We synthesise a minimal fs_kfold object with the right class; it must
# emit the pointer to cv_summary_tables()/cv_metrics_tables().
tf_kfold_fake <- structure(
  list(fold_summary = tf$fold_summary),
  class = c("fs_kfold", "list")
)
.expect_error(
  "Test 3: errors", "fs_kfold object redirected",
  expr    = quote(summarize_cv_results(cv_output = tf_kfold_fake)),
  pattern = "cv_summary_tables|cv_metrics_tables"
)
```
:::



# Test 4: Backward compatibility {#sec-test4}

Drops `grf_cuts` from `fold_summary` to simulate a
pre-feature-branch object.  Expects graceful NULL for the
GRF-dependent slots and a populated `identification_summary`.


::: {.cell}

```{.r .cell-code}
tf_nogrf <- tf
tf_nogrf$fold_summary <- tf$fold_summary[
  , setdiff(names(tf$fold_summary), "grf_cuts"),
  drop = FALSE
]

cv_diag_nogrf <- summarize_cv_results(tf_nogrf)

.tc("Test 4: backward-compat", "has_grf_cuts FALSE",
    isFALSE(cv_diag_nogrf$has_grf_cuts))
.tc("Test 4: backward-compat", "identification_summary populated",
    !is.null(cv_diag_nogrf$identification_summary))
.tc("Test 4: backward-compat", "grf_cut_summary NULL",
    is.null(cv_diag_nogrf$grf_cut_summary))
.tc("Test 4: backward-compat", "cut_vs_subgroup_xtab NULL",
    is.null(cv_diag_nogrf$cut_vs_subgroup_xtab))
.tc("Test 4: backward-compat", "no_subgroup_decomposition NULL",
    is.null(cv_diag_nogrf$no_subgroup_decomposition))
.tc("Test 4: backward-compat", "$data$grf_cuts NULL",
    is.null(cv_diag_nogrf$data$grf_cuts))
.tc("Test 4: backward-compat", "metrics_table still populated",
    !is.null(cv_diag_nogrf$metrics_table))
```
:::



# Test 5: Cross-reference against manual aggregation {#sec-test5}

Reproduces `grf_cut_summary` and `identification_summary` using
the `.clean_grf` + `.grf_dist` + `.build_sg_full` helpers lifted
directly from `gbsg_10fold_dmin_sensitivity.qmd`, then asserts
row-for-row equality with the function output.  This is the
refactor-safety check: the function should produce identical
results to the manual one-off aggregation the sensitivity QMD
performed.

To isolate the tabulation logic from the top-n truncation layer,
this test calls `summarize_cv_results(..., top_n = Inf)` so the
function returns the full frequency distribution without an
`"(N others)"` aggregation row.  The row-sum invariant under
truncation is covered separately in Test 1.


::: {.cell}

```{.r .cell-code}
# Untruncated variant for byte-exact comparison
cv_diag_full_topn <- summarize_cv_results(tf, top_n = Inf)
```
:::



::: {.cell}

```{.r .cell-code}
# Lifted verbatim from gbsg_10fold_dmin_sensitivity.qmd
.clean_grf_manual <- function(x) {
  if (is.na(x) || x == "") return("(no GRF cut)")
  x
}

.grf_dist_manual <- function(fsum, top_n = Inf) {
  cuts <- vapply(fsum$grf_cuts, .clean_grf_manual, character(1))
  tab  <- as.data.frame(table(cuts = cuts), stringsAsFactors = FALSE)
  tab  <- tab[order(-tab$Freq), , drop = FALSE]
  if (is.finite(top_n) && nrow(tab) > top_n) tab <- utils::head(tab, top_n)
  tab
}

.build_sg_full_manual <- function(fsum) {
  ifelse(
    is.na(fsum$sg1),
    "(no subgroup)",
    ifelse(
      is.na(fsum$sg2),
      fsum$sg1,
      paste(fsum$sg1, fsum$sg2, sep = " & ")
    )
  )
}
```
:::



::: {.cell}

```{.r .cell-code}
# Manual (no truncation)
manual_grf <- .grf_dist_manual(tf$fold_summary)
manual_grf_aligned <- data.frame(
  grf_cut = as.character(manual_grf$cuts),
  n_folds = as.integer(manual_grf$Freq),
  stringsAsFactors = FALSE
)
rownames(manual_grf_aligned) <- NULL

# Function (untruncated)
fn_grf <- cv_diag_full_topn$data$grf_cuts
fn_grf_aligned <- data.frame(
  grf_cut = fn_grf$grf_cut,
  n_folds = as.integer(fn_grf$n_folds),
  stringsAsFactors = FALSE
)
rownames(fn_grf_aligned) <- NULL

# Same number of rows?
.tc("Test 5: cross-ref", "grf_cut_summary row count matches manual",
    nrow(manual_grf_aligned) == nrow(fn_grf_aligned),
    detail = sprintf("manual=%d fn=%d",
                     nrow(manual_grf_aligned),
                     nrow(fn_grf_aligned)))

# Row-for-row equality
eq_ok <- isTRUE(all.equal(manual_grf_aligned, fn_grf_aligned))
.tc("Test 5: cross-ref", "grf_cut_summary rows match manual exactly",
    eq_ok,
    detail = if (eq_ok) "all.equal TRUE" else
      paste(utils::capture.output(
        all.equal(manual_grf_aligned, fn_grf_aligned)
      ), collapse = "; "))
```
:::



::: {.cell}

```{.r .cell-code}
# Manual identification tabulation (no truncation)
manual_sg <- .build_sg_full_manual(tf$fold_summary)
tab <- as.data.frame(table(Subgroup = manual_sg),
                     stringsAsFactors = FALSE)
tab <- tab[order(-tab$Freq), , drop = FALSE]

manual_id_aligned <- data.frame(
  Subgroup = as.character(tab$Subgroup),
  n_folds  = as.integer(tab$Freq),
  stringsAsFactors = FALSE
)
rownames(manual_id_aligned) <- NULL

fn_id <- cv_diag_full_topn$data$identification
fn_id_aligned <- data.frame(
  Subgroup = fn_id$Subgroup,
  n_folds  = as.integer(fn_id$n_folds),
  stringsAsFactors = FALSE
)
rownames(fn_id_aligned) <- NULL

.tc("Test 5: cross-ref", "identification_summary row count matches manual",
    nrow(manual_id_aligned) == nrow(fn_id_aligned),
    detail = sprintf("manual=%d fn=%d",
                     nrow(manual_id_aligned),
                     nrow(fn_id_aligned)))

eq_ok <- isTRUE(all.equal(manual_id_aligned, fn_id_aligned))
.tc("Test 5: cross-ref", "identification_summary rows match manual exactly",
    eq_ok,
    detail = if (eq_ok) "all.equal TRUE" else
      paste(utils::capture.output(
        all.equal(manual_id_aligned, fn_id_aligned)
      ), collapse = "; "))
```
:::



::: {.cell}

```{.r .cell-code}
# Manual no-subgroup count
manual_no_sg <- sum(manual_sg == "(no subgroup)")
# Function
fn_no_sg <- cv_diag$data$no_subgroup$n_folds[
  cv_diag$data$no_subgroup$Category == "Folds with no subgroup identified"
]

.tc("Test 5: cross-ref", "no_subgroup total matches manual",
    isTRUE(manual_no_sg == fn_no_sg),
    detail = sprintf("manual=%d fn=%d", manual_no_sg, fn_no_sg))

# Manual GRF-none-among-no-sg count
manual_grf_none_in_no_sg <- sum(
  (manual_sg == "(no subgroup)") &
    vapply(tf$fold_summary$grf_cuts, .clean_grf_manual,
           character(1)) == "(no GRF cut)"
)
fn_grf_none_in_no_sg <- cv_diag$data$no_subgroup$n_folds[
  cv_diag$data$no_subgroup$Category == "  GRF returned no cut"
]

.tc("Test 5: cross-ref", "decomposition 'GRF none' matches manual",
    isTRUE(manual_grf_none_in_no_sg == fn_grf_none_in_no_sg),
    detail = sprintf("manual=%d fn=%d",
                     manual_grf_none_in_no_sg, fn_grf_none_in_no_sg))
```
:::



# Test 6: Print method {#sec-test6}

Confirms `print.fs_cv_summary()` emits the expected header
lines and section markers.


::: {.cell}

```{.r .cell-code}
out <- utils::capture.output(print(cv_diag))

.tc("Test 6: print", "print output non-empty",
    length(out) > 0L)
.tc("Test 6: print", "print output contains 'Simulations'",
    any(grepl("Simulations:", out)))
.tc("Test 6: print", "print output contains 'Total sim x fold pairs'",
    any(grepl("Total sim x fold pairs", out)))
.tc("Test 6: print", "print output contains 'Top identified subgroups'",
    any(grepl("Top identified subgroups", out)))
.tc("Test 6: print", "print output contains 'Top GRF cuts'",
    any(grepl("Top GRF cuts", out)))
.tc("Test 6: print", "print output contains 'No-subgroup decomposition'",
    any(grepl("No-subgroup decomposition", out)))
```
:::


Full captured output:


```{.r .cell-code}
cat("```\n")
```

```

```{.r .cell-code}
cat(out, sep = "\n")
```

ForestSearch CV Diagnostic Summary
==================================
Simulations:            10
Folds per simulation:   10
Total sim x fold pairs: 100
grf_cuts available:     yes

Top identified subgroups:
  Rank                   Subgroup n_folds pct_folds
1    1                  {er <= 0}      48        48
2    2              (no subgroup)      21        21
3    3 !{age <= 43} & {age <= 47}      10        10
4    4 {age <= 47} & !{age <= 43}       5         5
5    5 {age <= 48} & !{age <= 43}       5         5

Top GRF cuts:
  Rank                         grf_cut n_folds pct_folds
1    1                         er <= 0      41        41
2    2 age <= 48 | age <= 43 | er <= 0      15        15
3    3 age <= 47 | age <= 43 | er <= 0      12        12
4    4                    (no GRF cut)      10        10
5    5 age <= 50 | age <= 43 | er <= 0       7         7

No-subgroup decomposition:
                                    Category n_folds pct_of_total pct_of_no_sg
1          Folds with no subgroup identified      21           21       100.00
2                        GRF returned no cut       9            9        42.86
3   GRF returned a cut; consistency rejected      12           12        57.14

Access gt-formatted tables via:
  $identification_summary, $grf_cut_summary, $cut_vs_subgroup_xtab,
  $no_subgroup_decomposition, $original_agreement, $metrics_table

```{.r .cell-code}
cat("\n```\n")
```

```


# Test 7: Phase B diagnostic columns {#sec-test7}

Validates the three new `fold_summary` columns populated by
`forestsearch_tenfold()` on the feature branch
(`pconsistency`, `training_fs_hr`, `n_candidates_evaluated`)
and the two new summary slots exposed by
`summarize_cv_results()` (`pconsistency_distribution`,
`fold_numeric_summary`).

## 7.1 Column presence and value domain


::: {.cell}

```{.r .cell-code}
has_pcons <- "pconsistency" %in% names(tf$fold_summary)
.tc("Test 7: pconsistency", "fold_summary$pconsistency column present",
    has_pcons,
    detail = if (has_pcons) "yes" else
      "NO -- rerun forestsearch_tenfold() with Phase B branch")

if (has_pcons) {
  pv <- tf$fold_summary$pconsistency

  .tc("Test 7: pconsistency", "pconsistency is numeric",
      is.numeric(pv),
      detail = sprintf("class=%s", class(pv)[1]))

  # NA iff no subgroup identified
  na_aligned <- all(
    is.na(pv) == (tf$fold_summary$any_found == 0L)
  )
  .tc("Test 7: pconsistency",
      "NA pconsistency iff any_found == 0",
      na_aligned,
      detail = sprintf(
        "NA count=%d; any_found==0 count=%d",
        sum(is.na(pv)),
        sum(tf$fold_summary$any_found == 0L)))

  # Domain check: non-NA values should be in [0, 1]
  pv_nona <- pv[!is.na(pv)]
  in_range <- length(pv_nona) == 0L ||
    (all(pv_nona >= 0) && all(pv_nona <= 1))
  .tc("Test 7: pconsistency",
      "non-NA pconsistency values in [0, 1]",
      in_range,
      detail = if (length(pv_nona) > 0L)
        sprintf("range = [%.3f, %.3f], median = %.3f",
                min(pv_nona), max(pv_nona), median(pv_nona))
      else "no identifying folds")

  # Given pconsistency.threshold = 0.90 in config, all identifying
  # Pcons should be >= that threshold (by construction of
  # forestsearch()'s selection rule).
  thresh <- fs_args_base$pconsistency.threshold
  above_thresh <- length(pv_nona) == 0L || all(pv_nona >= thresh)
  .tc("Test 7: pconsistency",
      sprintf("identifying Pcons >= threshold (%.2f)", thresh),
      above_thresh,
      detail = if (length(pv_nona) > 0L)
        sprintf("min = %.3f; threshold = %.2f", min(pv_nona), thresh)
      else "no identifying folds")
}
```
:::


## 7.2 Slot populated when column present


::: {.cell}

```{.r .cell-code}
if (has_pcons) {
  .tc("Test 7: pconsistency",
      "has_pconsistency metadata TRUE",
      isTRUE(cv_diag$has_pconsistency))
  .tc("Test 7: pconsistency",
      "pconsistency_distribution slot populated",
      !is.null(cv_diag$pconsistency_distribution))
  .tc("Test 7: pconsistency",
      "$data$pconsistency is data.frame",
      is.data.frame(cv_diag$data$pconsistency))

  pcd <- cv_diag$data$pconsistency
  .tc("Test 7: pconsistency",
      "pconsistency df has expected columns",
      all(c("Pcons_range", "n_folds", "pct_of_identifying") %in%
            names(pcd)))

  # Summary rows carry median/IQR/total — last three rows
  last_rows <- pcd$Pcons_range[(nrow(pcd) - 2L):nrow(pcd)]
  .tc("Test 7: pconsistency",
      "last three rows include Median/IQR/total summary",
      any(grepl("Median", last_rows)) &&
      any(grepl("IQR",    last_rows)) &&
      any(grepl("Identifying", last_rows)),
      detail = paste(last_rows, collapse = " | "))

  # Bin row counts should sum to number of identifying folds
  n_identifying <- sum(!is.na(tf$fold_summary$pconsistency))
  bin_rows <- pcd[1:6, , drop = FALSE]  # 6 bins from 0 to 1
  bin_sum <- sum(bin_rows$n_folds, na.rm = TRUE)
  .tc("Test 7: pconsistency",
      "bin n_folds sums to identifying-fold count",
      bin_sum == n_identifying,
      detail = sprintf("sum=%d exp=%d", bin_sum, n_identifying))
}
```
:::


## 7.3 Backward compatibility (column stripped)


::: {.cell}

```{.r .cell-code}
if (has_pcons) {
  tf_nopcons <- tf
  tf_nopcons$fold_summary <- tf$fold_summary[
    , setdiff(names(tf$fold_summary), "pconsistency"),
    drop = FALSE
  ]
  cv_diag_nopcons <- summarize_cv_results(tf_nopcons)

  .tc("Test 7: pconsistency",
      "has_pconsistency FALSE when column absent",
      isFALSE(cv_diag_nopcons$has_pconsistency))
  .tc("Test 7: pconsistency",
      "pconsistency_distribution NULL when column absent",
      is.null(cv_diag_nopcons$pconsistency_distribution))
  .tc("Test 7: pconsistency",
      "$data$pconsistency NULL when column absent",
      is.null(cv_diag_nopcons$data$pconsistency))
  .tc("Test 7: pconsistency",
      "other slots still populated when pconsistency absent",
      !is.null(cv_diag_nopcons$identification_summary) &&
      !is.null(cv_diag_nopcons$grf_cut_summary))
}
```
:::


## 7.4 Rendered pconsistency_distribution table


::: {.cell}

```{.r .cell-code}
if (has_pcons) cv_diag$pconsistency_distribution
```

::: {.cell-output-display}

```{=html}
<div id="erhxutjobc" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#erhxutjobc table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#erhxutjobc thead, #erhxutjobc tbody, #erhxutjobc tfoot, #erhxutjobc tr, #erhxutjobc td, #erhxutjobc th {
  border-style: none;
}

#erhxutjobc p {
  margin: 0;
  padding: 0;
}

#erhxutjobc .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 13px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#erhxutjobc .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#erhxutjobc .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#erhxutjobc .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#erhxutjobc .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#erhxutjobc .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#erhxutjobc .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#erhxutjobc .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#erhxutjobc .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#erhxutjobc .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#erhxutjobc .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#erhxutjobc .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#erhxutjobc .gt_spanner_row {
  border-bottom-style: hidden;
}

#erhxutjobc .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#erhxutjobc .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#erhxutjobc .gt_from_md > :first-child {
  margin-top: 0;
}

#erhxutjobc .gt_from_md > :last-child {
  margin-bottom: 0;
}

#erhxutjobc .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#erhxutjobc .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#erhxutjobc .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#erhxutjobc .gt_row_group_first td {
  border-top-width: 2px;
}

#erhxutjobc .gt_row_group_first th {
  border-top-width: 2px;
}

#erhxutjobc .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#erhxutjobc .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#erhxutjobc .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#erhxutjobc .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#erhxutjobc .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#erhxutjobc .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#erhxutjobc .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#erhxutjobc .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#erhxutjobc .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#erhxutjobc .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#erhxutjobc .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#erhxutjobc .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#erhxutjobc .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#erhxutjobc .gt_left {
  text-align: left;
}

#erhxutjobc .gt_center {
  text-align: center;
}

#erhxutjobc .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#erhxutjobc .gt_font_normal {
  font-weight: normal;
}

#erhxutjobc .gt_font_bold {
  font-weight: bold;
}

#erhxutjobc .gt_font_italic {
  font-style: italic;
}

#erhxutjobc .gt_super {
  font-size: 65%;
}

#erhxutjobc .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#erhxutjobc .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#erhxutjobc .gt_indent_1 {
  text-indent: 5px;
}

#erhxutjobc .gt_indent_2 {
  text-indent: 10px;
}

#erhxutjobc .gt_indent_3 {
  text-indent: 15px;
}

#erhxutjobc .gt_indent_4 {
  text-indent: 20px;
}

#erhxutjobc .gt_indent_5 {
  text-indent: 25px;
}

#erhxutjobc .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#erhxutjobc div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="3" class="gt_heading gt_title gt_font_normal" style><span data-qmd-base64="KipQY29ucyBkaXN0cmlidXRpb24gYW1vbmcgaWRlbnRpZnlpbmcgZm9sZHMqKg=="><span class='gt_from_md'><strong>Pcons distribution among identifying folds</strong></span></span></td>
    </tr>
    <tr class="gt_heading">
      <td colspan="3" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>Out of 100 sim x fold pairs</td>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Pcons_range">Pcons range</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="n_folds">n folds</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pct_of_identifying">% of identifying folds</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Pcons_range" class="gt_row gt_left">&lt;0.5</td>
<td headers="n_folds" class="gt_row gt_right">0</td>
<td headers="pct_of_identifying" class="gt_row gt_right">0.00</td></tr>
    <tr><td headers="Pcons_range" class="gt_row gt_left">0.5-0.7</td>
<td headers="n_folds" class="gt_row gt_right">0</td>
<td headers="pct_of_identifying" class="gt_row gt_right">0.00</td></tr>
    <tr><td headers="Pcons_range" class="gt_row gt_left">0.7-0.8</td>
<td headers="n_folds" class="gt_row gt_right">0</td>
<td headers="pct_of_identifying" class="gt_row gt_right">0.00</td></tr>
    <tr><td headers="Pcons_range" class="gt_row gt_left">0.8-0.85</td>
<td headers="n_folds" class="gt_row gt_right">0</td>
<td headers="pct_of_identifying" class="gt_row gt_right">0.00</td></tr>
    <tr><td headers="Pcons_range" class="gt_row gt_left">0.85-0.95</td>
<td headers="n_folds" class="gt_row gt_right">43</td>
<td headers="pct_of_identifying" class="gt_row gt_right">54.43</td></tr>
    <tr><td headers="Pcons_range" class="gt_row gt_left">&gt;=0.95</td>
<td headers="n_folds" class="gt_row gt_right">36</td>
<td headers="pct_of_identifying" class="gt_row gt_right">45.57</td></tr>
    <tr><td headers="Pcons_range" class="gt_row gt_left" style="background-color: #F0F0F0; font-weight: bold;">Median (identifying): 0.940</td>
<td headers="n_folds" class="gt_row gt_right" style="background-color: #F0F0F0; font-weight: bold;">–</td>
<td headers="pct_of_identifying" class="gt_row gt_right" style="background-color: #F0F0F0; font-weight: bold;">–</td></tr>
    <tr><td headers="Pcons_range" class="gt_row gt_left" style="background-color: #F0F0F0; font-weight: bold;">IQR (identifying): 0.920 - 0.970</td>
<td headers="n_folds" class="gt_row gt_right" style="background-color: #F0F0F0; font-weight: bold;">–</td>
<td headers="pct_of_identifying" class="gt_row gt_right" style="background-color: #F0F0F0; font-weight: bold;">–</td></tr>
    <tr><td headers="Pcons_range" class="gt_row gt_left" style="background-color: #F0F0F0; font-weight: bold;">Identifying folds total</td>
<td headers="n_folds" class="gt_row gt_right" style="background-color: #F0F0F0; font-weight: bold;">79</td>
<td headers="pct_of_identifying" class="gt_row gt_right" style="background-color: #F0F0F0; font-weight: bold;">79.00</td></tr>
  </tbody>
  <tfoot>
    <tr class="gt_sourcenotes">
      <td class="gt_sourcenote" colspan="3">Restricted to folds that identified a subgroup.  A concentration near the pconsistency.threshold (typically 0.90) suggests many identifications are marginal; concentration above 0.95 suggests robust identifications.  Values for rejected candidates are not surfaced -- only the selected candidate's Pcons is captured.</td>
    </tr>
  </tfoot>
</table>
</div>
```

:::
:::


## 7.5 `training_fs_hr` column


::: {.cell}

```{.r .cell-code}
has_hr <- "training_fs_hr" %in% names(tf$fold_summary)
.tc("Test 7: pconsistency",
    "fold_summary$training_fs_hr column present",
    has_hr,
    detail = if (has_hr) "yes" else
      "NO -- rerun forestsearch_tenfold() with Phase B branch")

if (has_hr) {
  hr <- tf$fold_summary$training_fs_hr
  .tc("Test 7: pconsistency", "training_fs_hr is numeric",
      is.numeric(hr),
      detail = sprintf("class=%s", class(hr)[1]))

  # NA iff no subgroup identified (same gate as pconsistency)
  na_aligned <- all(is.na(hr) == (tf$fold_summary$any_found == 0L))
  .tc("Test 7: pconsistency",
      "NA training_fs_hr iff any_found == 0",
      na_aligned,
      detail = sprintf(
        "NA count=%d; any_found==0 count=%d",
        sum(is.na(hr)),
        sum(tf$fold_summary$any_found == 0L)))

  # Domain: non-NA HRs should be positive and finite
  hr_nona <- hr[!is.na(hr)]
  in_domain <- length(hr_nona) == 0L ||
    (all(hr_nona > 0) && all(is.finite(hr_nona)))
  .tc("Test 7: pconsistency",
      "non-NA training_fs_hr values positive and finite",
      in_domain,
      detail = if (length(hr_nona) > 0L)
        sprintf("range = [%.3f, %.3f], median = %.3f",
                min(hr_nona), max(hr_nona), median(hr_nona))
      else "no identifying folds")

  # GBSG at sg_focus = "maxSG" with hr.threshold = 1.25 forces
  # identifying HRs above the threshold.  Sanity check.
  hr_thresh <- fs_args_base$hr.threshold
  above_hr_thresh <- length(hr_nona) == 0L || all(hr_nona >= hr_thresh)
  .tc("Test 7: pconsistency",
      sprintf("identifying HRs >= hr.threshold (%.2f)", hr_thresh),
      above_hr_thresh,
      detail = if (length(hr_nona) > 0L)
        sprintf("min = %.3f; threshold = %.2f", min(hr_nona), hr_thresh)
      else "no identifying folds")
}
```
:::


## 7.6 `n_candidates_evaluated` column


::: {.cell}

```{.r .cell-code}
has_cands <- "n_candidates_evaluated" %in% names(tf$fold_summary)
.tc("Test 7: pconsistency",
    "fold_summary$n_candidates_evaluated column present",
    has_cands,
    detail = if (has_cands) "yes" else
      "NO -- rerun forestsearch_tenfold() with Phase B branch")

if (has_cands) {
  nc <- tf$fold_summary$n_candidates_evaluated
  .tc("Test 7: pconsistency", "n_candidates_evaluated is integer",
      is.integer(nc),
      detail = sprintf("class=%s", class(nc)[1]))

  # Every fold that identified a subgroup must have evaluated
  # at least 1 candidate.  The converse is NOT true: a fold can
  # have evaluated many candidates and still identify none.
  identifying_idx <- which(tf$fold_summary$any_found == 1L)
  if (length(identifying_idx) > 0L) {
    nc_identifying <- nc[identifying_idx]
    nc_identifying_ok <- all(
      !is.na(nc_identifying) & nc_identifying >= 1L
    )
    .tc("Test 7: pconsistency",
        "identifying folds have n_candidates_evaluated >= 1",
        nc_identifying_ok,
        detail = sprintf(
          "identifying folds: n=%d, min candidates=%s",
          length(identifying_idx),
          if (all(is.na(nc_identifying))) "all NA"
          else as.character(min(nc_identifying, na.rm = TRUE))))
  }

  # Domain: non-NA values should be non-negative
  nc_nona <- nc[!is.na(nc)]
  in_domain <- length(nc_nona) == 0L ||
    (all(nc_nona >= 0L) && all(is.finite(nc_nona)))
  .tc("Test 7: pconsistency",
      "non-NA n_candidates_evaluated values non-negative",
      in_domain,
      detail = if (length(nc_nona) > 0L)
        sprintf("range = [%d, %d], median = %d",
                as.integer(min(nc_nona)),
                as.integer(max(nc_nona)),
                as.integer(median(nc_nona)))
      else "no populated values")
}
```
:::


## 7.7 `fold_numeric_summary` slot


::: {.cell}

```{.r .cell-code}
.tc("Test 7: pconsistency",
    "fold_numeric_summary slot populated",
    !is.null(cv_diag$fold_numeric_summary))
.tc("Test 7: pconsistency",
    "$data$fold_numeric_summary is data.frame",
    is.data.frame(cv_diag$data$fold_numeric_summary))

fns <- cv_diag$data$fold_numeric_summary
expected_cols <- c("Metric", "n", "median", "Q1", "Q3", "min",
                   "max", "n_na", "Context")
.tc("Test 7: pconsistency",
    "fold_numeric_summary has expected columns",
    all(expected_cols %in% names(fns)),
    detail = paste(setdiff(expected_cols, names(fns)),
                   collapse = ", "))

# Given GBSG + Phase B, all four metric rows should appear.
expected_metrics <- c(
  "Test-fold size (n_test)",
  "Pcons (identifying folds)",
  "Training-fold subgroup effect (identifying folds)",
  "Candidates evaluated per fold"
)
.tc("Test 7: pconsistency",
    "fold_numeric_summary includes all four Phase B metric rows",
    all(expected_metrics %in% fns$Metric),
    detail = paste(setdiff(expected_metrics, fns$Metric),
                   collapse = ", "))

# Row counts must match the per-column non-NA counts.
for (spec in list(
  list(metric = "Pcons (identifying folds)",
       col    = "pconsistency",
       filt   = function(x) !is.na(x)),
  list(metric = "Training-fold subgroup effect (identifying folds)",
       col    = "training_fs_hr",
       filt   = function(x) !is.na(x)),
  list(metric = "Candidates evaluated per fold",
       col    = "n_candidates_evaluated",
       filt   = function(x) !is.na(x)),
  list(metric = "Test-fold size (n_test)",
       col    = "n_test",
       filt   = function(x) !is.na(x))
)) {
  if (!spec$col %in% names(tf$fold_summary)) next
  row    <- fns[fns$Metric == spec$metric, , drop = FALSE]
  actual <- sum(spec$filt(tf$fold_summary[[spec$col]]))
  .tc("Test 7: pconsistency",
      sprintf("%s: n matches non-NA count", spec$metric),
      nrow(row) == 1L && row$n == actual,
      detail = sprintf("row$n=%d actual=%d",
                       if (nrow(row) == 1L) row$n else NA_integer_,
                       actual))
}
```
:::


## 7.8 Rendered fold_numeric_summary table


::: {.cell}

```{.r .cell-code}
cv_diag$fold_numeric_summary
```

::: {.cell-output-display}

```{=html}
<div id="apdwjlpbdy" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#apdwjlpbdy table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#apdwjlpbdy thead, #apdwjlpbdy tbody, #apdwjlpbdy tfoot, #apdwjlpbdy tr, #apdwjlpbdy td, #apdwjlpbdy th {
  border-style: none;
}

#apdwjlpbdy p {
  margin: 0;
  padding: 0;
}

#apdwjlpbdy .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 13px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#apdwjlpbdy .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#apdwjlpbdy .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#apdwjlpbdy .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#apdwjlpbdy .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#apdwjlpbdy .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#apdwjlpbdy .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#apdwjlpbdy .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#apdwjlpbdy .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#apdwjlpbdy .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#apdwjlpbdy .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#apdwjlpbdy .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#apdwjlpbdy .gt_spanner_row {
  border-bottom-style: hidden;
}

#apdwjlpbdy .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#apdwjlpbdy .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#apdwjlpbdy .gt_from_md > :first-child {
  margin-top: 0;
}

#apdwjlpbdy .gt_from_md > :last-child {
  margin-bottom: 0;
}

#apdwjlpbdy .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#apdwjlpbdy .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#apdwjlpbdy .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#apdwjlpbdy .gt_row_group_first td {
  border-top-width: 2px;
}

#apdwjlpbdy .gt_row_group_first th {
  border-top-width: 2px;
}

#apdwjlpbdy .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#apdwjlpbdy .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#apdwjlpbdy .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#apdwjlpbdy .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#apdwjlpbdy .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#apdwjlpbdy .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#apdwjlpbdy .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#apdwjlpbdy .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#apdwjlpbdy .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#apdwjlpbdy .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#apdwjlpbdy .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#apdwjlpbdy .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#apdwjlpbdy .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#apdwjlpbdy .gt_left {
  text-align: left;
}

#apdwjlpbdy .gt_center {
  text-align: center;
}

#apdwjlpbdy .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#apdwjlpbdy .gt_font_normal {
  font-weight: normal;
}

#apdwjlpbdy .gt_font_bold {
  font-weight: bold;
}

#apdwjlpbdy .gt_font_italic {
  font-style: italic;
}

#apdwjlpbdy .gt_super {
  font-size: 65%;
}

#apdwjlpbdy .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#apdwjlpbdy .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#apdwjlpbdy .gt_indent_1 {
  text-indent: 5px;
}

#apdwjlpbdy .gt_indent_2 {
  text-indent: 10px;
}

#apdwjlpbdy .gt_indent_3 {
  text-indent: 15px;
}

#apdwjlpbdy .gt_indent_4 {
  text-indent: 20px;
}

#apdwjlpbdy .gt_indent_5 {
  text-indent: 25px;
}

#apdwjlpbdy .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#apdwjlpbdy div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="9" class="gt_heading gt_title gt_font_normal" style><span data-qmd-base64="KipGb2xkLWxldmVsIG51bWVyaWMgZGlhZ25vc3RpY3MqKg=="><span class='gt_from_md'><strong>Fold-level numeric diagnostics</strong></span></span></td>
    </tr>
    <tr class="gt_heading">
      <td colspan="9" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>Across 100 sim x fold pairs; per-metric denominator shown in Context</td>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Metric">Metric</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="n">n</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="median">Median</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Q1">Q1</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Q3">Q3</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="min">Min</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="max">Max</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="n_na">n NA</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Context">Context</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Metric" class="gt_row gt_left">Test-fold size (n_test)</td>
<td headers="n" class="gt_row gt_right">100</td>
<td headers="median" class="gt_row gt_right">69.000</td>
<td headers="Q1" class="gt_row gt_right">68.000</td>
<td headers="Q3" class="gt_row gt_right">69.000</td>
<td headers="min" class="gt_row gt_right">68.000</td>
<td headers="max" class="gt_row gt_right">69.000</td>
<td headers="n_na" class="gt_row gt_right">0</td>
<td headers="Context" class="gt_row gt_left">All sim x fold pairs</td></tr>
    <tr><td headers="Metric" class="gt_row gt_left">Pcons (identifying folds)</td>
<td headers="n" class="gt_row gt_right">79</td>
<td headers="median" class="gt_row gt_right">0.940</td>
<td headers="Q1" class="gt_row gt_right">0.920</td>
<td headers="Q3" class="gt_row gt_right">0.970</td>
<td headers="min" class="gt_row gt_right">0.900</td>
<td headers="max" class="gt_row gt_right">1.000</td>
<td headers="n_na" class="gt_row gt_right">21</td>
<td headers="Context" class="gt_row gt_left">Folds that identified a subgroup</td></tr>
    <tr><td headers="Metric" class="gt_row gt_left">Training-fold subgroup effect (identifying folds)</td>
<td headers="n" class="gt_row gt_right">79</td>
<td headers="median" class="gt_row gt_right">2.034</td>
<td headers="Q1" class="gt_row gt_right">1.937</td>
<td headers="Q3" class="gt_row gt_right">2.341</td>
<td headers="min" class="gt_row gt_right">1.646</td>
<td headers="max" class="gt_row gt_right">3.075</td>
<td headers="n_na" class="gt_row gt_right">21</td>
<td headers="Context" class="gt_row gt_left">Natural scale; in-sample; optimistically biased</td></tr>
    <tr><td headers="Metric" class="gt_row gt_left">Candidates evaluated per fold</td>
<td headers="n" class="gt_row gt_right">100</td>
<td headers="median" class="gt_row gt_right">4.500</td>
<td headers="Q1" class="gt_row gt_right">2.750</td>
<td headers="Q3" class="gt_row gt_right">7.000</td>
<td headers="min" class="gt_row gt_right">1.000</td>
<td headers="max" class="gt_row gt_right">10.000</td>
<td headers="n_na" class="gt_row gt_right">0</td>
<td headers="Context" class="gt_row gt_left">Folds where consistency stage ran</td></tr>
  </tbody>
  <tfoot>
    <tr class="gt_sourcenotes">
      <td class="gt_sourcenote" colspan="9">Training-fold subgroup effects are in-sample estimates on the effect measure's natural scale (HR for survival; OR, RR, IRR, RD, IRD, or MD for GLM).  They are optimistically biased and shown here for diagnostic comparison across folds only; do not use for inference.</td>
    </tr>
  </tfoot>
</table>
</div>
```

:::
:::


## 7.9 Backward compat: strip ALL Phase B columns

Also confirms that stripping Phase B columns wholesale doesn't
disturb the non-Phase-B slots.


::: {.cell}

```{.r .cell-code}
tf_noB <- tf
tf_noB$fold_summary <- tf$fold_summary[
  , setdiff(names(tf$fold_summary),
            c("pconsistency", "training_fs_hr",
              "n_candidates_evaluated")),
  drop = FALSE
]
cv_diag_noB <- summarize_cv_results(tf_noB)

.tc("Test 7: pconsistency",
    "pconsistency_distribution NULL when all Phase B stripped",
    is.null(cv_diag_noB$pconsistency_distribution))
.tc("Test 7: pconsistency",
    "fold_numeric_summary still populated (n_test always present)",
    !is.null(cv_diag_noB$fold_numeric_summary))
.tc("Test 7: pconsistency",
    "non-Phase-B slots still populated when Phase B stripped",
    !is.null(cv_diag_noB$identification_summary) &&
      !is.null(cv_diag_noB$grf_cut_summary) &&
      !is.null(cv_diag_noB$no_subgroup_decomposition))

# fold_numeric_summary should have exactly 1 row (n_test only)
fns_noB <- cv_diag_noB$data$fold_numeric_summary
.tc("Test 7: pconsistency",
    "fold_numeric_summary has only n_test row when Phase B stripped",
    nrow(fns_noB) == 1L && fns_noB$Metric[1] == "Test-fold size (n_test)",
    detail = sprintf("rows=%d, metrics=%s",
                     nrow(fns_noB),
                     paste(fns_noB$Metric, collapse = ", ")))
```
:::



# Test summary {#sec-summary}


::: {.cell}

```{.r .cell-code}
n_total <- nrow(test_results)
n_pass  <- sum(test_results$Status == "PASS")
n_fail  <- n_total - n_pass

gt::gt(test_results,
       groupname_col = "Section",
       row_group_as_column = TRUE) |>
  gt::tab_header(
    title    = gt::md(sprintf(
      "**Living test results: %d / %d passing**", n_pass, n_total)),
    subtitle = sprintf(
      "GBSG, %d sims x %d folds, dmin.grf = %g, frac.tau = %g",
      Ksims_cv, Kfolds_cv,
      fs_args_base$dmin.grf, fs_args_base$frac.tau)
  ) |>
  gt::tab_style(
    style    = gt::cell_fill(color = "#d4edda"),
    locations = gt::cells_body(
      columns = "Status", rows = Status == "PASS"
    )
  ) |>
  gt::tab_style(
    style    = list(
      gt::cell_fill(color = "#f8d7da"),
      gt::cell_text(weight = "bold")
    ),
    locations = gt::cells_body(
      columns = "Status", rows = Status == "FAIL"
    )
  ) |>
  gt::cols_label(
    Test   = "Test case",
    Status = "Status",
    Detail = "Detail"
  ) |>
  gt::tab_options(table.font.size = gt::px(12))
```

::: {.cell-output-display}

```{=html}
<div id="wmhoiswpch" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#wmhoiswpch table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#wmhoiswpch thead, #wmhoiswpch tbody, #wmhoiswpch tfoot, #wmhoiswpch tr, #wmhoiswpch td, #wmhoiswpch th {
  border-style: none;
}

#wmhoiswpch p {
  margin: 0;
  padding: 0;
}

#wmhoiswpch .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 12px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#wmhoiswpch .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#wmhoiswpch .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#wmhoiswpch .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#wmhoiswpch .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#wmhoiswpch .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#wmhoiswpch .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#wmhoiswpch .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#wmhoiswpch .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#wmhoiswpch .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#wmhoiswpch .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#wmhoiswpch .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#wmhoiswpch .gt_spanner_row {
  border-bottom-style: hidden;
}

#wmhoiswpch .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#wmhoiswpch .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#wmhoiswpch .gt_from_md > :first-child {
  margin-top: 0;
}

#wmhoiswpch .gt_from_md > :last-child {
  margin-bottom: 0;
}

#wmhoiswpch .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#wmhoiswpch .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#wmhoiswpch .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#wmhoiswpch .gt_row_group_first td {
  border-top-width: 2px;
}

#wmhoiswpch .gt_row_group_first th {
  border-top-width: 2px;
}

#wmhoiswpch .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#wmhoiswpch .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#wmhoiswpch .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#wmhoiswpch .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#wmhoiswpch .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#wmhoiswpch .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#wmhoiswpch .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#wmhoiswpch .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#wmhoiswpch .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#wmhoiswpch .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#wmhoiswpch .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#wmhoiswpch .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#wmhoiswpch .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#wmhoiswpch .gt_left {
  text-align: left;
}

#wmhoiswpch .gt_center {
  text-align: center;
}

#wmhoiswpch .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#wmhoiswpch .gt_font_normal {
  font-weight: normal;
}

#wmhoiswpch .gt_font_bold {
  font-weight: bold;
}

#wmhoiswpch .gt_font_italic {
  font-style: italic;
}

#wmhoiswpch .gt_super {
  font-size: 65%;
}

#wmhoiswpch .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#wmhoiswpch .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#wmhoiswpch .gt_indent_1 {
  text-indent: 5px;
}

#wmhoiswpch .gt_indent_2 {
  text-indent: 10px;
}

#wmhoiswpch .gt_indent_3 {
  text-indent: 15px;
}

#wmhoiswpch .gt_indent_4 {
  text-indent: 20px;
}

#wmhoiswpch .gt_indent_5 {
  text-indent: 25px;
}

#wmhoiswpch .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#wmhoiswpch div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="4" class="gt_heading gt_title gt_font_normal" style><span data-qmd-base64="KipMaXZpbmcgdGVzdCByZXN1bHRzOiA5MyAvIDkzIHBhc3NpbmcqKg=="><span class='gt_from_md'><strong>Living test results: 93 / 93 passing</strong></span></span></td>
    </tr>
    <tr class="gt_heading">
      <td colspan="4" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>GBSG, 10 sims x 10 folds, dmin.grf = 12, frac.tau = 0.8</td>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="a::stub"></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Test">Test case</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Status">Status</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Detail">Detail</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_row_group_first"><td headers="Baseline stub_2_1 stub_1" rowspan="5" class="gt_row gt_left gt_stub_row_group">Baseline</td>
<td headers="Baseline stub_2_1 Test" class="gt_row gt_left">fold_summary has required columns</td>
<td headers="Baseline stub_2_1 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Baseline stub_2_1 Detail" class="gt_row gt_left">sim, fold, n_test, sg1, sg2, grf_cuts, pconsistency, training_fs_hr, n_candidates_evaluated, any_found</td></tr>
    <tr><td headers="Baseline stub_2_2 Test" class="gt_row gt_left">fold_summary$grf_cuts present</td>
<td headers="Baseline stub_2_2 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Baseline stub_2_2 Detail" class="gt_row gt_left">Required for feature-branch diagnostic slots</td></tr>
    <tr><td headers="Baseline stub_2_3 Test" class="gt_row gt_left">fold_summary$pconsistency present</td>
<td headers="Baseline stub_2_3 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Baseline stub_2_3 Detail" class="gt_row gt_left">Required for Phase B pconsistency_distribution slot</td></tr>
    <tr><td headers="Baseline stub_2_4 Test" class="gt_row gt_left">fold_summary$training_fs_hr present</td>
<td headers="Baseline stub_2_4 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Baseline stub_2_4 Detail" class="gt_row gt_left">Phase B diagnostic column</td></tr>
    <tr><td headers="Baseline stub_2_5 Test" class="gt_row gt_left">fold_summary$n_candidates_evaluated present</td>
<td headers="Baseline stub_2_5 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Baseline stub_2_5 Detail" class="gt_row gt_left">Phase B diagnostic column</td></tr>
    <tr class="gt_row_group_first"><td headers="Test 1: minimal stub_2_6 stub_1" rowspan="20" class="gt_row gt_left gt_stub_row_group">Test 1: minimal</td>
<td headers="Test 1: minimal stub_2_6 Test" class="gt_row gt_left">slot identification_summary matches expectation</td>
<td headers="Test 1: minimal stub_2_6 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_6 Detail" class="gt_row gt_left">populated = TRUE</td></tr>
    <tr><td headers="Test 1: minimal stub_2_7 Test" class="gt_row gt_left">slot grf_cut_summary matches expectation</td>
<td headers="Test 1: minimal stub_2_7 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_7 Detail" class="gt_row gt_left">populated = TRUE</td></tr>
    <tr><td headers="Test 1: minimal stub_2_8 Test" class="gt_row gt_left">slot cut_vs_subgroup_xtab matches expectation</td>
<td headers="Test 1: minimal stub_2_8 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_8 Detail" class="gt_row gt_left">populated = TRUE</td></tr>
    <tr><td headers="Test 1: minimal stub_2_9 Test" class="gt_row gt_left">slot no_subgroup_decomposition matches expectation</td>
<td headers="Test 1: minimal stub_2_9 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_9 Detail" class="gt_row gt_left">populated = TRUE</td></tr>
    <tr><td headers="Test 1: minimal stub_2_10 Test" class="gt_row gt_left">slot original_agreement matches expectation</td>
<td headers="Test 1: minimal stub_2_10 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_10 Detail" class="gt_row gt_left">populated = FALSE</td></tr>
    <tr><td headers="Test 1: minimal stub_2_11 Test" class="gt_row gt_left">slot metrics_table matches expectation</td>
<td headers="Test 1: minimal stub_2_11 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_11 Detail" class="gt_row gt_left">populated = TRUE</td></tr>
    <tr><td headers="Test 1: minimal stub_2_12 Test" class="gt_row gt_left">slot plots matches expectation</td>
<td headers="Test 1: minimal stub_2_12 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_12 Detail" class="gt_row gt_left">populated = FALSE</td></tr>
    <tr><td headers="Test 1: minimal stub_2_13 Test" class="gt_row gt_left">n_sims matches tf$sims</td>
<td headers="Test 1: minimal stub_2_13 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_13 Detail" class="gt_row gt_left">obs=10 exp=10</td></tr>
    <tr><td headers="Test 1: minimal stub_2_14 Test" class="gt_row gt_left">n_folds matches tf$Kfolds</td>
<td headers="Test 1: minimal stub_2_14 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_14 Detail" class="gt_row gt_left">obs=10 exp=10</td></tr>
    <tr><td headers="Test 1: minimal stub_2_15 Test" class="gt_row gt_left">total_pairs matches nrow(fold_summary)</td>
<td headers="Test 1: minimal stub_2_15 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_15 Detail" class="gt_row gt_left">obs=100 exp=100</td></tr>
    <tr><td headers="Test 1: minimal stub_2_16 Test" class="gt_row gt_left">has_grf_cuts == TRUE (feature branch)</td>
<td headers="Test 1: minimal stub_2_16 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_16 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 1: minimal stub_2_17 Test" class="gt_row gt_left">class == 'fs_cv_summary'</td>
<td headers="Test 1: minimal stub_2_17 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_17 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 1: minimal stub_2_18 Test" class="gt_row gt_left">$data$identification is data.frame</td>
<td headers="Test 1: minimal stub_2_18 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_18 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 1: minimal stub_2_19 Test" class="gt_row gt_left">$data$grf_cuts is data.frame</td>
<td headers="Test 1: minimal stub_2_19 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_19 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 1: minimal stub_2_20 Test" class="gt_row gt_left">identification df has expected columns</td>
<td headers="Test 1: minimal stub_2_20 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_20 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 1: minimal stub_2_21 Test" class="gt_row gt_left">grf_cuts df has expected columns</td>
<td headers="Test 1: minimal stub_2_21 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_21 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 1: minimal stub_2_22 Test" class="gt_row gt_left">identification n_folds sums to total pairs (including (N others))</td>
<td headers="Test 1: minimal stub_2_22 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_22 Detail" class="gt_row gt_left">obs=100 exp=100</td></tr>
    <tr><td headers="Test 1: minimal stub_2_23 Test" class="gt_row gt_left">grf_cuts n_folds sums to total pairs (including (N others))</td>
<td headers="Test 1: minimal stub_2_23 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_23 Detail" class="gt_row gt_left">obs=100 exp=100</td></tr>
    <tr><td headers="Test 1: minimal stub_2_24 Test" class="gt_row gt_left">identification_summary row count &lt;= top_n + 1</td>
<td headers="Test 1: minimal stub_2_24 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_24 Detail" class="gt_row gt_left">rows=13; has (N others) row=FALSE</td></tr>
    <tr><td headers="Test 1: minimal stub_2_25 Test" class="gt_row gt_left">grf_cut_summary row count &lt;= top_n + 1</td>
<td headers="Test 1: minimal stub_2_25 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 1: minimal stub_2_25 Detail" class="gt_row gt_left">rows=16; has (N others) row=TRUE</td></tr>
    <tr class="gt_row_group_first"><td headers="Test 2: full stub_2_26 stub_1" rowspan="8" class="gt_row gt_left gt_stub_row_group">Test 2: full</td>
<td headers="Test 2: full stub_2_26 Test" class="gt_row gt_left">$original_agreement populated</td>
<td headers="Test 2: full stub_2_26 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 2: full stub_2_26 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 2: full stub_2_27 Test" class="gt_row gt_left">$plots populated</td>
<td headers="Test 2: full stub_2_27 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 2: full stub_2_27 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 2: full stub_2_28 Test" class="gt_row gt_left">$plots$identification is ggplot</td>
<td headers="Test 2: full stub_2_28 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 2: full stub_2_28 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 2: full stub_2_29 Test" class="gt_row gt_left">$plots$grf_cut is ggplot (has_grf_cuts TRUE)</td>
<td headers="Test 2: full stub_2_29 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 2: full stub_2_29 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 2: full stub_2_30 Test" class="gt_row gt_left">original_agreement df has Metric, Value</td>
<td headers="Test 2: full stub_2_30 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 2: full stub_2_30 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 2: full stub_2_31 Test" class="gt_row gt_left">original_agreement includes Exact subgroup match row</td>
<td headers="Test 2: full stub_2_31 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 2: full stub_2_31 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 2: full stub_2_32 Test" class="gt_row gt_left">original_agreement includes Exact GRF-cut match row</td>
<td headers="Test 2: full stub_2_32 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 2: full stub_2_32 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 2: full stub_2_33 Test" class="gt_row gt_left">original_agreement includes Both match row</td>
<td headers="Test 2: full stub_2_33 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 2: full stub_2_33 Detail" class="gt_row gt_left"></td></tr>
    <tr class="gt_row_group_first"><td headers="Test 3: errors stub_2_34 stub_1" rowspan="5" class="gt_row gt_left gt_stub_row_group">Test 3: errors</td>
<td headers="Test 3: errors stub_2_34 Test" class="gt_row gt_left">non-list input rejected</td>
<td headers="Test 3: errors stub_2_34 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 3: errors stub_2_34 Detail" class="gt_row gt_left">`cv_output` must be a list returned by forestsearch_tenfold().</td></tr>
    <tr><td headers="Test 3: errors stub_2_35 Test" class="gt_row gt_left">empty fold_summary rejected</td>
<td headers="Test 3: errors stub_2_35 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 3: errors stub_2_35 Detail" class="gt_row gt_left">`cv_output$fold_summary` is missing or empty. summarize_cv_results() requires the per-fold summary produced by forestsea</td></tr>
    <tr><td headers="Test 3: errors stub_2_36 Test" class="gt_row gt_left">missing required column rejected</td>
<td headers="Test 3: errors stub_2_36 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 3: errors stub_2_36 Detail" class="gt_row gt_left">`fold_summary` is missing required column(s): sg1</td></tr>
    <tr><td headers="Test 3: errors stub_2_37 Test" class="gt_row gt_left">non-positive top_n rejected</td>
<td headers="Test 3: errors stub_2_37 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 3: errors stub_2_37 Detail" class="gt_row gt_left">`top_n` must be a positive integer (or `Inf` to disable truncation).</td></tr>
    <tr><td headers="Test 3: errors stub_2_38 Test" class="gt_row gt_left">fs_kfold object redirected</td>
<td headers="Test 3: errors stub_2_38 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 3: errors stub_2_38 Detail" class="gt_row gt_left">`summarize_cv_results()` operates on forestsearch_tenfold() output (class 'fs_tenfold'), which exposes a `fold_summary` </td></tr>
    <tr class="gt_row_group_first"><td headers="Test 4: backward-compat stub_2_39 stub_1" rowspan="7" class="gt_row gt_left gt_stub_row_group">Test 4: backward-compat</td>
<td headers="Test 4: backward-compat stub_2_39 Test" class="gt_row gt_left">has_grf_cuts FALSE</td>
<td headers="Test 4: backward-compat stub_2_39 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 4: backward-compat stub_2_39 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 4: backward-compat stub_2_40 Test" class="gt_row gt_left">identification_summary populated</td>
<td headers="Test 4: backward-compat stub_2_40 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 4: backward-compat stub_2_40 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 4: backward-compat stub_2_41 Test" class="gt_row gt_left">grf_cut_summary NULL</td>
<td headers="Test 4: backward-compat stub_2_41 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 4: backward-compat stub_2_41 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 4: backward-compat stub_2_42 Test" class="gt_row gt_left">cut_vs_subgroup_xtab NULL</td>
<td headers="Test 4: backward-compat stub_2_42 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 4: backward-compat stub_2_42 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 4: backward-compat stub_2_43 Test" class="gt_row gt_left">no_subgroup_decomposition NULL</td>
<td headers="Test 4: backward-compat stub_2_43 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 4: backward-compat stub_2_43 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 4: backward-compat stub_2_44 Test" class="gt_row gt_left">$data$grf_cuts NULL</td>
<td headers="Test 4: backward-compat stub_2_44 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 4: backward-compat stub_2_44 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 4: backward-compat stub_2_45 Test" class="gt_row gt_left">metrics_table still populated</td>
<td headers="Test 4: backward-compat stub_2_45 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 4: backward-compat stub_2_45 Detail" class="gt_row gt_left"></td></tr>
    <tr class="gt_row_group_first"><td headers="Test 5: cross-ref stub_2_46 stub_1" rowspan="6" class="gt_row gt_left gt_stub_row_group">Test 5: cross-ref</td>
<td headers="Test 5: cross-ref stub_2_46 Test" class="gt_row gt_left">grf_cut_summary row count matches manual</td>
<td headers="Test 5: cross-ref stub_2_46 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 5: cross-ref stub_2_46 Detail" class="gt_row gt_left">manual=16 fn=16</td></tr>
    <tr><td headers="Test 5: cross-ref stub_2_47 Test" class="gt_row gt_left">grf_cut_summary rows match manual exactly</td>
<td headers="Test 5: cross-ref stub_2_47 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 5: cross-ref stub_2_47 Detail" class="gt_row gt_left">all.equal TRUE</td></tr>
    <tr><td headers="Test 5: cross-ref stub_2_48 Test" class="gt_row gt_left">identification_summary row count matches manual</td>
<td headers="Test 5: cross-ref stub_2_48 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 5: cross-ref stub_2_48 Detail" class="gt_row gt_left">manual=13 fn=13</td></tr>
    <tr><td headers="Test 5: cross-ref stub_2_49 Test" class="gt_row gt_left">identification_summary rows match manual exactly</td>
<td headers="Test 5: cross-ref stub_2_49 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 5: cross-ref stub_2_49 Detail" class="gt_row gt_left">all.equal TRUE</td></tr>
    <tr><td headers="Test 5: cross-ref stub_2_50 Test" class="gt_row gt_left">no_subgroup total matches manual</td>
<td headers="Test 5: cross-ref stub_2_50 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 5: cross-ref stub_2_50 Detail" class="gt_row gt_left">manual=21 fn=21</td></tr>
    <tr><td headers="Test 5: cross-ref stub_2_51 Test" class="gt_row gt_left">decomposition 'GRF none' matches manual</td>
<td headers="Test 5: cross-ref stub_2_51 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 5: cross-ref stub_2_51 Detail" class="gt_row gt_left">manual=9 fn=9</td></tr>
    <tr class="gt_row_group_first"><td headers="Test 6: print stub_2_52 stub_1" rowspan="6" class="gt_row gt_left gt_stub_row_group">Test 6: print</td>
<td headers="Test 6: print stub_2_52 Test" class="gt_row gt_left">print output non-empty</td>
<td headers="Test 6: print stub_2_52 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 6: print stub_2_52 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 6: print stub_2_53 Test" class="gt_row gt_left">print output contains 'Simulations'</td>
<td headers="Test 6: print stub_2_53 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 6: print stub_2_53 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 6: print stub_2_54 Test" class="gt_row gt_left">print output contains 'Total sim x fold pairs'</td>
<td headers="Test 6: print stub_2_54 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 6: print stub_2_54 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 6: print stub_2_55 Test" class="gt_row gt_left">print output contains 'Top identified subgroups'</td>
<td headers="Test 6: print stub_2_55 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 6: print stub_2_55 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 6: print stub_2_56 Test" class="gt_row gt_left">print output contains 'Top GRF cuts'</td>
<td headers="Test 6: print stub_2_56 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 6: print stub_2_56 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 6: print stub_2_57 Test" class="gt_row gt_left">print output contains 'No-subgroup decomposition'</td>
<td headers="Test 6: print stub_2_57 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 6: print stub_2_57 Detail" class="gt_row gt_left"></td></tr>
    <tr class="gt_row_group_first"><td headers="Test 7: pconsistency stub_2_58 stub_1" rowspan="36" class="gt_row gt_left gt_stub_row_group">Test 7: pconsistency</td>
<td headers="Test 7: pconsistency stub_2_58 Test" class="gt_row gt_left">fold_summary$pconsistency column present</td>
<td headers="Test 7: pconsistency stub_2_58 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_58 Detail" class="gt_row gt_left">yes</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_59 Test" class="gt_row gt_left">pconsistency is numeric</td>
<td headers="Test 7: pconsistency stub_2_59 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_59 Detail" class="gt_row gt_left">class=numeric</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_60 Test" class="gt_row gt_left">NA pconsistency iff any_found == 0</td>
<td headers="Test 7: pconsistency stub_2_60 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_60 Detail" class="gt_row gt_left">NA count=21; any_found==0 count=21</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_61 Test" class="gt_row gt_left">non-NA pconsistency values in [0, 1]</td>
<td headers="Test 7: pconsistency stub_2_61 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_61 Detail" class="gt_row gt_left">range = [0.900, 1.000], median = 0.940</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_62 Test" class="gt_row gt_left">identifying Pcons &gt;= threshold (0.90)</td>
<td headers="Test 7: pconsistency stub_2_62 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_62 Detail" class="gt_row gt_left">min = 0.900; threshold = 0.90</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_63 Test" class="gt_row gt_left">has_pconsistency metadata TRUE</td>
<td headers="Test 7: pconsistency stub_2_63 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_63 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_64 Test" class="gt_row gt_left">pconsistency_distribution slot populated</td>
<td headers="Test 7: pconsistency stub_2_64 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_64 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_65 Test" class="gt_row gt_left">$data$pconsistency is data.frame</td>
<td headers="Test 7: pconsistency stub_2_65 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_65 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_66 Test" class="gt_row gt_left">pconsistency df has expected columns</td>
<td headers="Test 7: pconsistency stub_2_66 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_66 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_67 Test" class="gt_row gt_left">last three rows include Median/IQR/total summary</td>
<td headers="Test 7: pconsistency stub_2_67 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_67 Detail" class="gt_row gt_left">Median (identifying): 0.940 | IQR (identifying): 0.920 - 0.970 | Identifying folds total</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_68 Test" class="gt_row gt_left">bin n_folds sums to identifying-fold count</td>
<td headers="Test 7: pconsistency stub_2_68 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_68 Detail" class="gt_row gt_left">sum=79 exp=79</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_69 Test" class="gt_row gt_left">has_pconsistency FALSE when column absent</td>
<td headers="Test 7: pconsistency stub_2_69 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_69 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_70 Test" class="gt_row gt_left">pconsistency_distribution NULL when column absent</td>
<td headers="Test 7: pconsistency stub_2_70 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_70 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_71 Test" class="gt_row gt_left">$data$pconsistency NULL when column absent</td>
<td headers="Test 7: pconsistency stub_2_71 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_71 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_72 Test" class="gt_row gt_left">other slots still populated when pconsistency absent</td>
<td headers="Test 7: pconsistency stub_2_72 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_72 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_73 Test" class="gt_row gt_left">fold_summary$training_fs_hr column present</td>
<td headers="Test 7: pconsistency stub_2_73 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_73 Detail" class="gt_row gt_left">yes</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_74 Test" class="gt_row gt_left">training_fs_hr is numeric</td>
<td headers="Test 7: pconsistency stub_2_74 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_74 Detail" class="gt_row gt_left">class=numeric</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_75 Test" class="gt_row gt_left">NA training_fs_hr iff any_found == 0</td>
<td headers="Test 7: pconsistency stub_2_75 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_75 Detail" class="gt_row gt_left">NA count=21; any_found==0 count=21</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_76 Test" class="gt_row gt_left">non-NA training_fs_hr values positive and finite</td>
<td headers="Test 7: pconsistency stub_2_76 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_76 Detail" class="gt_row gt_left">range = [1.646, 3.075], median = 2.034</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_77 Test" class="gt_row gt_left">identifying HRs &gt;= hr.threshold (1.25)</td>
<td headers="Test 7: pconsistency stub_2_77 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_77 Detail" class="gt_row gt_left">min = 1.646; threshold = 1.25</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_78 Test" class="gt_row gt_left">fold_summary$n_candidates_evaluated column present</td>
<td headers="Test 7: pconsistency stub_2_78 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_78 Detail" class="gt_row gt_left">yes</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_79 Test" class="gt_row gt_left">n_candidates_evaluated is integer</td>
<td headers="Test 7: pconsistency stub_2_79 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_79 Detail" class="gt_row gt_left">class=integer</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_80 Test" class="gt_row gt_left">identifying folds have n_candidates_evaluated &gt;= 1</td>
<td headers="Test 7: pconsistency stub_2_80 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_80 Detail" class="gt_row gt_left">identifying folds: n=79, min candidates=1</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_81 Test" class="gt_row gt_left">non-NA n_candidates_evaluated values non-negative</td>
<td headers="Test 7: pconsistency stub_2_81 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_81 Detail" class="gt_row gt_left">range = [1, 10], median = 4</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_82 Test" class="gt_row gt_left">fold_numeric_summary slot populated</td>
<td headers="Test 7: pconsistency stub_2_82 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_82 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_83 Test" class="gt_row gt_left">$data$fold_numeric_summary is data.frame</td>
<td headers="Test 7: pconsistency stub_2_83 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_83 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_84 Test" class="gt_row gt_left">fold_numeric_summary has expected columns</td>
<td headers="Test 7: pconsistency stub_2_84 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_84 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_85 Test" class="gt_row gt_left">fold_numeric_summary includes all four Phase B metric rows</td>
<td headers="Test 7: pconsistency stub_2_85 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_85 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_86 Test" class="gt_row gt_left">Pcons (identifying folds): n matches non-NA count</td>
<td headers="Test 7: pconsistency stub_2_86 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_86 Detail" class="gt_row gt_left">row$n=79 actual=79</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_87 Test" class="gt_row gt_left">Training-fold subgroup effect (identifying folds): n matches non-NA count</td>
<td headers="Test 7: pconsistency stub_2_87 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_87 Detail" class="gt_row gt_left">row$n=79 actual=79</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_88 Test" class="gt_row gt_left">Candidates evaluated per fold: n matches non-NA count</td>
<td headers="Test 7: pconsistency stub_2_88 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_88 Detail" class="gt_row gt_left">row$n=100 actual=100</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_89 Test" class="gt_row gt_left">Test-fold size (n_test): n matches non-NA count</td>
<td headers="Test 7: pconsistency stub_2_89 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_89 Detail" class="gt_row gt_left">row$n=100 actual=100</td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_90 Test" class="gt_row gt_left">pconsistency_distribution NULL when all Phase B stripped</td>
<td headers="Test 7: pconsistency stub_2_90 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_90 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_91 Test" class="gt_row gt_left">fold_numeric_summary still populated (n_test always present)</td>
<td headers="Test 7: pconsistency stub_2_91 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_91 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_92 Test" class="gt_row gt_left">non-Phase-B slots still populated when Phase B stripped</td>
<td headers="Test 7: pconsistency stub_2_92 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_92 Detail" class="gt_row gt_left"></td></tr>
    <tr><td headers="Test 7: pconsistency stub_2_93 Test" class="gt_row gt_left">fold_numeric_summary has only n_test row when Phase B stripped</td>
<td headers="Test 7: pconsistency stub_2_93 Status" class="gt_row gt_left" style="background-color: #D4EDDA;">PASS</td>
<td headers="Test 7: pconsistency stub_2_93 Detail" class="gt_row gt_left">rows=1, metrics=Test-fold size (n_test)</td></tr>
  </tbody>
  
</table>
</div>
```

:::

```{.r .cell-code}
# Also emit a one-line headline suitable for grep-based CI.
if (n_fail == 0L) {
  message(sprintf("ALL TESTS PASSED: %d / %d", n_pass, n_total))
} else {
  message(sprintf("TESTS FAILED: %d / %d failing", n_fail, n_total))
  print(test_results[test_results$Status == "FAIL",
                     c("Section", "Test", "Detail")])
}
```

::: {.cell-output .cell-output-stderr}

```
ALL TESTS PASSED: 93 / 93
```


:::
:::



# Session info {#sec-session}


::: {.cell}

```{.r .cell-code}
sessionInfo()
```

::: {.cell-output .cell-output-stdout}

```
R version 4.5.2 (2025-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Tahoe 26.4

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Los_Angeles
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] doFuture_1.2.1      foreach_1.5.2       future_1.70.0      
[4] ggplot2_4.0.2       gt_1.3.0            data.table_1.18.2.1
[7] survival_3.8-6      forestsearch_0.2.0 

loaded via a namespace (and not attached):
 [1] sass_0.4.10          generics_0.1.4       xml2_1.5.2          
 [4] shape_1.4.6.1        lattice_0.22-9       listenv_0.10.1      
 [7] digest_0.6.39        magrittr_2.0.5       grf_2.6.1           
[10] evaluate_1.0.5       grid_4.5.2           RColorBrewer_1.1-3  
[13] iterators_1.0.14     policytree_1.2.4     fastmap_1.2.0       
[16] jsonlite_2.0.0       Matrix_1.7-5         glmnet_4.1-10       
[19] scales_1.4.0         weightedsurv_0.1.0   codetools_0.2-20    
[22] cli_3.6.6            rlang_1.2.0          litedown_0.9        
[25] parallelly_1.46.1    future.apply_1.20.2  commonmark_2.0.0    
[28] splines_4.5.2        base64enc_0.1-6      withr_3.0.2         
[31] yaml_2.3.12          otel_0.2.0           tools_4.5.2         
[34] parallel_4.5.2       dplyr_1.2.1          globals_0.19.1      
[37] vctrs_0.7.3          R6_2.6.1             lifecycle_1.0.5     
[40] randomForest_4.7-1.2 fs_2.0.1             htmlwidgets_1.6.4   
[43] pkgconfig_2.0.3      progressr_0.19.0     pillar_1.11.1       
[46] gtable_0.3.6         glue_1.8.0           Rcpp_1.1.1          
[49] xfun_0.57            tibble_3.3.1         tidyselect_1.2.1    
[52] rstudioapi_0.18.0    knitr_1.51           farver_2.1.2        
[55] patchwork_1.3.2      htmltools_0.5.9      labeling_0.4.3      
[58] rmarkdown_2.31       compiler_4.5.2       S7_0.2.1            
[61] markdown_2.0        
```


:::
:::


