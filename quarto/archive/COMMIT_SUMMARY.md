# Git Commit Summary
# Branch: feature/glm-extension
# Date: 2026-04-04

## Commit message (suggested)

```
fix: censoring model predict() dimension mismatch + DGM field compatibility

Two bugs that prevent generate_aft_dgm_flex() from working on non-GBSG
datasets when AIC selects a covariate censoring model.

Bug #1 (CRITICAL) — R/generate_aft_dgm_helpers.R
  prepare_censoring_model(): survreg() was fitted with a pre-built matrix
  formula (~ X_cens).  predict.survreg() silently falls back to the
  training data when newdata lacks the matrix variable, producing
  nrow(df_work) predictions instead of nrow(df_super).  Fix: use
  reformulate() to build a proper formula with column names.  Same fix
  applied to prepare_censoring_model_legacy(), which also lacked the
  mu_cens subtraction (double-counting the intercept).  Defensive
  dimension checks added to both functions.

Bug #2 (MEDIUM) — R/simulation_tables.R
  build_classification_table(), build_estimation_table(), and
  interpret_estimation_table() accessed dgm$hr_H_true, dgm$hr_Hc_true,
  dgm$hr_causal, and dgm$df_super_rand$flag.harm directly — fields that
  only exist on gbsg_dgm objects from setup_gbsg_dgm().  Fix: use
  get_dgm_hr() and dgm$df_super_rand %||% dgm$df_super for
  compatibility with aft_dgm_flex objects from generate_aft_dgm_flex().
```

## Files changed

### R/generate_aft_dgm_helpers.R

| Location | What changed |
|----------|-------------|
| `prepare_censoring_model()` ~L146-230 | Replaced `X_cens <- as.matrix(...); survreg(~ X_cens, ...)` with `reformulate()` + `survreg(cens_formula, ...)`. Added 14-line defensive dimension check after predict() calls. |
| `prepare_censoring_model_legacy()` ~L1017-1085 | Same `reformulate()` fix. Added `- mu_cens` to both predict() calls (was missing — double-counted intercept). Added dimension check. |

### R/simulation_tables.R

| Location | What changed |
|----------|-------------|
| `build_classification_table()` ~L112-130 | `dgm$hr_causal` → `get_dgm_hr(dgm, "hr_overall")` with fallback. `dgm$df_super_rand$flag.harm` → `dgm$df_super_rand %||% dgm$df_super` with `flag.harm`/`flag_harm` column auto-detect. `dgm$hr_H_true` / `dgm$hr_Hc_true` → `get_dgm_hr()`. |
| `build_estimation_table()` ~L358-393 | `dgm$hr_H_true` / `dgm$hr_Hc_true` → `get_dgm_hr(dgm, "hr_H")` / `get_dgm_hr(dgm, "hr_Hc")`. CDE fallback broadened: `dgm$df_super_rand` → `dgm$df_super_rand %||% dgm$df_super`. |
| `interpret_estimation_table()` ~L757-762 | Same `get_dgm_hr()` replacement. |

## How to revert

```bash
# From the forestsearch repo root on feature/glm-extension:
git log --oneline -5          # find the commit hash
git revert <commit-hash>      # creates a revert commit
# OR to discard entirely:
git reset --hard <hash-before-this-commit>
```

## Post-install verification

```r
devtools::document()
devtools::install()

# Quick smoke test — should succeed without dimension error:
library(forestsearch)
library(speff2trial)
df <- subset(ACTG175, arms %in% c(1, 3))
df$y     <- df$days
df$event <- df$cens
df$treat <- ifelse(df$arms == 1, 1L, 0L)
df$z1    <- as.factor(ifelse(df$age > 34, 1L, 0L))
df$z2    <- as.factor(ifelse(df$preanti <= 744.5, 1L, 0L))

dgm <- generate_aft_dgm_flex(
  data            = df,
  outcome_var     = "y",
  event_var       = "event",
  treatment_var   = "treat",
  continuous_vars = character(0),
  factor_vars     = c("z1", "z2"),
  subgroup_vars   = c("z1", "z2"),
  subgroup_cuts   = list(z1 = 1L, z2 = 1L),
  model           = "alt",
  n_super         = 5000,
  verbose         = TRUE
)
cat("DGM created successfully. nrow(df_super):", nrow(dgm$df_super), "\n")

# Verify GBSG path still works (regression check):
dgm_gbsg <- setup_gbsg_dgm(model = "alt")
cat("GBSG DGM OK. HR(H):", round(dgm_gbsg$hr_H_true, 3), "\n")
```
