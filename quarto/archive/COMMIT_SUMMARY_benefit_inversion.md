# Git Commit Summary — Benefit-Search Auto-Inversion
# Branch: feature/glm-extension
# Date: 2026-04-04

## Commit message (suggested)

```
feat: auto-invert HRs when subgroup_notation = "benefit"

When subgroup_notation = "benefit", all table and interpretation
functions now automatically invert hazard ratio values from the
switched treatment scale (HR > 1) to the original scale (HR < 1)
before computing bias, formatting labels, and rendering output.

This eliminates the need for manual inversion in user-facing
documents and ensures that benefit-search results are always
displayed on the clinically interpretable scale.

New internal helpers:
- invert_hr_columns(): inverts 14 HR-related columns in results
- invert_dgm_hrs(): inverts hazard_ratios list and legacy fields

Functions updated (auto-inversion + Q/Qc labels + footnotes):
- build_classification_table()
- build_estimation_table()
- interpret_estimation_table()
- format_oc_results()
- render_reference_table()

QC: Synthetic data tests confirm element-wise inversion of all
HR columns, preservation of non-HR columns, correct bias sign
on original scale, and get_dgm_hr() compatibility.
```

## Files changed

### R/simulation_tables.R (+79 lines)

| Change | Location |
|--------|----------|
| New `invert_hr_columns()` | After `get_sg_labels()`, ~L59-87 |
| New `invert_dgm_hrs()` | ~L89-120 |
| Auto-invert in `build_classification_table()` | Inside scenario loop, after `res`/`dgm` extraction |
| Auto-invert in `build_estimation_table()` | After `L <- get_sg_labels()`, before `res <- as.data.table()` |
| Auto-invert in `interpret_estimation_table()` | Same position |
| (Previously applied) Q/Qc labels, notation-aware footnote, bold threshold fix | Throughout |

### R/oc_analyses.R (+4 lines)

| Change | Location |
|--------|----------|
| Auto-invert in `format_oc_results()` | After `L <- get_sg_labels()`, inverts `results` only (no `dgm` param) |

### quarto/actg175_simulations_benefit.qmd (-57 lines)

| Change | Description |
|--------|-------------|
| Removed `invert-to-original` chunk | Manual `invert_hr_columns()` / `invert_dgm_hrs()` definitions and `_orig` objects no longer needed |
| Simplified reporting calls | Pass raw `results_alt`, `dgm_calibrated` directly — package auto-inverts |
| Subtitle HRs | Use `1 / dgm_calibrated$hazard_ratios$...` inline for display only |

## How to revert

```bash
git log --oneline -5
git revert <commit-hash>
```

## Post-install verification

```r
devtools::document()
devtools::install()

# Quick smoke test with synthetic data
library(forestsearch)

# Mock switched-scale results
mock_res <- data.frame(
  analysis  = "FS",
  any.H     = 1,
  hr.H.hat  = 2.66,  # switched scale
  hr.Hc.hat = 1.19,
  hr.H.true = 2.04,
  hr.Hc.true = 1.11,
  sens = 0.80, spec = 0.90, ppv = 0.70, npv = 0.95,
  size.H = 258, size.Hc = 742
)

# Mock DGM
mock_dgm <- list(
  hazard_ratios = list(
    harm_subgroup = 2.04,
    no_harm_subgroup = 1.11,
    overall = 1.30
  ),
  df_super = data.frame(flag_harm = c(rep(1, 130), rep(0, 370)))
)

# Harm notation (default) — values unchanged
cat("=== Harm notation (default) ===\n")
# build_estimation_table uses raw values internally

# Benefit notation — values auto-inverted
cat("\n=== Benefit notation ===\n")
cat("Input:  hr.H.hat = 2.66 (switched)\n")
cat("Expect: hr.H.hat ~ 0.376 (original, < 1)\n")
# The table functions handle this automatically:
# build_estimation_table(mock_res, mock_dgm, subgroup_notation = "benefit")
# interpret_estimation_table(mock_res, mock_dgm, subgroup_notation = "benefit")

# Verify inversion helpers directly
res_inv <- forestsearch:::invert_hr_columns(mock_res)
dgm_inv <- forestsearch:::invert_dgm_hrs(mock_dgm)
cat(sprintf("Inverted: hr.H.hat = %.3f (should be < 1)\n", res_inv$hr.H.hat))
cat(sprintf("Inverted: dgm HR(Q) = %.3f (should be ~ 0.49)\n",
    dgm_inv$hazard_ratios$harm_subgroup))
stopifnot(res_inv$hr.H.hat < 1)
stopifnot(dgm_inv$hazard_ratios$harm_subgroup < 1)
stopifnot(identical(mock_res$sens, res_inv$sens))  # non-HR unchanged
cat("\nAll checks passed.\n")
```
