# L_eff Calibration Workflow

## Directory Structure

```
forestsearch/
+-- R/
|   +-- leff_calibration_io.R        # save/load helpers (174 lines)
|   +-- suggest_thresholds.R         # threshold frontier
|   +-- test_hte_crump.R             # HTE tests
|
+-- quarto/
|   +-- calibration_binary_leff.qmd  # PRODUCES binary .rds
|   +-- calibration_survival_leff_grid.qmd  # PRODUCES survival .rds
|   +-- [calibration_count_leff.qmd] # PRODUCES count .rds (MISSING)
|   +-- selection_binary_frontier.qmd     # CONSUMES binary .rds
|   +-- selection_survival_frontier.qmd   # CONSUMES survival .rds
|   +-- selection_glm_guide.qmd           # CONSUMES all .rds
|   +-- validation_glm_simulation_study.qmd # CONSUMES binary .rds
|   +-- _output/
|       +-- calibration_binary_leff.rds
|       +-- calibration_survival_leff.rds
|       +-- calibration_count_leff.rds    (once created)
```


## Data Flow

```
[Calibration doc]                    [.rds file]                [Selection doc]
                                                               
calibration_binary_leff.qmd          _output/                  selection_binary_frontier.qmd
  run 5000 H0 sims at 4 N  ──save──> calibration_binary_       ──load──> C_prior, alpha_prior
  fit power law: C, alpha              leff.rds                           instant grid evaluation
                                                               
calibration_survival_leff_grid.qmd                             selection_survival_frontier.qmd
  run 500 H0 sims at 4 N   ──save──> calibration_survival_     ──load──> C_prior, alpha_prior
  fit power law: C, alpha              leff.rds                           instant grid evaluation
                                                                
                                                               validation_glm_simulation_study.qmd
                                     calibration_binary_        ──load──> L_eff_C_binary, ...
                                      leff.rds                           threshold selection
```


## Pattern A: Calibration Document (Producer)

At the END of calibration_binary_leff.qmd, add:

```r
# ---- Save calibration results ----
save_leff_calibration(
  C = C_fit,
  alpha = alpha_fit,
  n_min = 60,
  outcome_type = "binary",
  cal_data = cal_df,   # data.frame with N, sim_fpr, P1, L_eff
  dgm_description = paste(
    "Binary DGM, 4 confounders (bm1, bm2, age, ecog),",
    "p(event)=0.30, 30% harm prevalence"),
  n_sims_per_N = 5000L
)
```

This writes `quarto/_data/calibration_binary_leff.rds`.


## Pattern B: Selection Document (Consumer)

At the TOP of selection_binary_frontier.qmd, replace the
hardcoded priors:

```r
# ---- Load calibration (replaces hardcoded C_prior/alpha_prior) ----
cal <- load_leff_calibration("binary")
C_prior     <- cal$C
alpha_prior <- cal$alpha
n_min_val   <- cal$n_min

cat(sprintf("Using L_eff from: %s\n", cal$dgm_description))
cat(sprintf("Calibrated: %s\n", format(cal$timestamp)))
```

If the .rds file doesn't exist, `load_leff_calibration()` stops
with a clear error message telling you to run the calibration
document first.


## Pattern C: Validation Document (Consumer)

In validation_glm_simulation_study.qmd, replace the hardcoded
L_eff values:

```r
# ---- Load L_eff calibrations ----
cal_binary <- load_leff_calibration("binary")
L_eff_C_binary     <- cal_binary$C
L_eff_alpha_binary <- cal_binary$alpha

# Count: use binary calibration as approximation (or count-specific
# once calibration_count_leff.qmd has been run)
cal_count <- tryCatch(
  load_leff_calibration("count"),
  error = function(e) {
    message("Count calibration not found; using binary as proxy.")
    cal_binary
  }
)
L_eff_C_count     <- cal_count$C
L_eff_alpha_count <- cal_count$alpha
```

This pattern lets the count pipeline automatically upgrade from
the binary approximation to count-specific values once the count
calibration document is rendered.


## Pattern D: Quick L_eff Lookup

For ad-hoc use or in suggest_thresholds():

```r
# What is L_eff at N = 800 for binary?
L_eff_800 <- get_leff(800, "binary")
cat(sprintf("L_eff at N=800: %.2f\n", L_eff_800))

# Compare across outcomes
for (ot in c("binary", "survival")) {
  L <- get_leff(700, ot)
  cat(sprintf("  %s: L_eff(700) = %.2f\n", ot, L))
}
```


## What the .rds Bundle Contains

```r
cal <- load_leff_calibration("binary")
str(cal, max.level = 1)
# List of 10
#  $ C                   : num 0.22
#  $ alpha               : num 1.3
#  $ n_min               : int 60
#  $ outcome_type        : chr "binary"
#  $ dgm_description     : chr "Binary DGM, 4 confounders ..."
#  $ n_sims_per_N        : int 5000
#  $ timestamp           : POSIXct "2026-03-30 ..."
#  $ R_version           : chr "R version 4.5.2 ..."
#  $ forestsearch_version: chr "0.1.0"
#  $ cal_data            : data.frame (N, sim_fpr, P1, L_eff)
#  $ extra               : list()
```

The `cal_data` data frame preserves the raw calibration points
so the L_eff fit can be re-examined or re-plotted without
re-running the simulations.


## Migration Checklist

To migrate existing documents to the .rds workflow:

1. Place `leff_calibration_io.R` in `R/`

2. In `calibration_binary_leff.qmd`:
   - Add `source(here::here("R", "leff_calibration_io.R"))` to setup
   - Add `save_leff_calibration(...)` at the end
   - Render once to produce the .rds

3. In `calibration_survival_leff_grid.qmd`:
   - Same pattern as step 2, with `outcome_type = "survival"`

4. In `selection_binary_frontier.qmd`:
   - Replace `C_prior <- 0.220; alpha_prior <- 1.298` with
     `cal <- load_leff_calibration("binary")`
   - Keep the old values as comments for reference

5. In `validation_glm_simulation_study.qmd`:
   - Replace `L_eff_C_binary <- 0.220` with
     `cal <- load_leff_calibration("binary")`

6. Test: delete the .rds, try to render a selection document,
   confirm it fails with a clear error message.


## Future: Count Calibration

When `calibration_count_leff.qmd` is created and rendered:

1. It saves `calibration_count_leff.rds` to `_output/`
2. `validation_glm_simulation_study.qmd` automatically picks it up
   via the `tryCatch` pattern in Pattern C above
3. The count thresholds update from the binary approximation to
   count-specific values without any manual edits
