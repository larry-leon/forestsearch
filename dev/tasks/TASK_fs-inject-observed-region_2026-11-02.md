# TASK: `inject_mrct_structure()` — observed-region mode (forestsearch, 2026-11-02)

**Repository:** forestsearch. **Branch:** create `fs-inject-observed-region` from `master` (at `0b1f5ba7` or later). Not on `mrct-article-r3`.

**R/ callout: changes behaviour, additively.** One new option in `inject_mrct_structure()` (`R/mrct_inject_dgm.R`). Defaults unchanged: every existing call builds the same DGM as before. The method is not changed.

## 0. First action

Copy this file to `dev/tasks/TASK_fs-inject-observed-region_2026-11-02.md` and commit it alone ("task: inject_mrct_structure observed-region mode").

## 1. Why

`inject_mrct_structure()` always plants its region with `make_region_model()`. A trial with an observed region (the analyst's own trial; in the article, replicate 5 treated as the analyst's data) cannot be built with the copula expansion, and a DGM built another way (`generate_aft_dgm_flex()`) lacks `dgm$mrct$region$z_name`, which `ahr_population()` requires. The observed-region mode removes both gaps.

## 2. The change

New element of the `region` list: `observed`, the name of a 0/1 column of `seed_data` (a 2-level factor is recoded as the other covariates are). When `region$observed` is given:

1. **Name.** `region$name` is `region$observed`; `z_name` is `paste0("z_", region$name)` as now.
2. **Expansion.** The region column is added to `covs`, so `expand_covariates()` expands it with the other covariates (copula or bootstrap, as `expand` says). `latent_R`, when supplied, must then include the region column (the existing `.check_latent_R()` check enforces this; say so in the `latent_R` documentation).
3. **No region model.** `make_region_model()` is not called. `df_pop[[region$name]]` is the expanded column; `seed_df[[region$name]]` is the observed column. The two-levels guard on the seed is replaced by an error if the observed column does not have both levels.
4. **Prognostic log-HR.** `region$loghr = NULL` (the default in this mode): the region's effect is estimated by the AFT fit, i.e. `z_<region>` is not in `set_beta_spec`. A numeric `region$loghr` fixes it, as in the planted mode. Read `generate_aft_dgm_flex()` for how an empty `set_beta_spec` must be passed when no term is fixed (no `x3`, no `region_treat_loghr`, `loghr = NULL`).
5. **Incompatible arguments.** In this mode, an error if any of `prevalence`, `or_pred`, `x1_vars`, `or_x1`, `scale`, `band` is supplied in `region` (they describe the planted model). Detect "supplied" against the user's `region` list before `modifyList()` fills the defaults.
6. **The record.** `dgm$mrct$region` gains `observed = TRUE` (and `observed = FALSE` in the planted mode); `prevalence_target`, `alpha0`, `coefs`, `scale`, `band` are `NULL` in this mode; `prevalence_pop` as now; `prognostic_loghr` is the value used: the fixed `loghr`, or the estimated one read from `dgm$model_params$b0[z_name]` (check from source that `b0` holds the log-HR scale; the article's build check reads `exp(b0["z_region"])` as a hazard ratio). `dgm$mrct$region_model` is `NULL`. `classes$X1` is `character(0)`.

**Readers of the record.** Grep `R/` for every read of `$mrct$region` fields and `$mrct$region_model` (`ahr_population()` and the other analytic functions, `analysis_covariates()`, `mrct_truth_by_region()`, `plot_region_modifier()`, `scale_dgm_hazard()`, `mrct_region_sims()`, `mrct_case_study()`, anything else). If any reads a field that is `NULL` in this mode, list it in the notes and stop; do not change other functions in this task.

**Documentation.** Roxygen for `region`: the `observed` element, what it does, the error for the planted-model elements, and `loghr = NULL` meaning estimated. One sentence in `@details`: an observed region is carried through the expansion with the covariates; nothing about its membership is modelled. NEWS: one bullet under the development heading.

## 3. Tests

New file `tests/testthat/test-inject-observed-region.R`, on ACTG 175 (`speff2trial`, `skip_if_not_installed`) with an observed region column made in the test from a fixed rule (e.g. `symptom` plus a seeded Bernoulli draw), `n_super` small, `expand = "bootstrap"` so latentcor is not needed (one further test with `expand = "copula"`, `skip_if_not_installed("latentcor")`):

- the build succeeds; `dgm$mrct$region$observed` is TRUE and `z_name` names a 0/1 column of `df_super`;
- with `expand = "bootstrap"`, every `df_super` row's region equals its source row's;
- `loghr = NULL`: `prognostic_loghr` equals `b0[z_name]`; numeric `loghr`: it equals the value given;
- each planted-model element in `region` gives an error; a one-level region column gives an error;
- `ahr_population()` accepts the DGM, and `ahr_limit(pop, "region", ...)` returns finite values.

Existing tests of `inject_mrct_structure()` and of the analytic layer pass unchanged (they are the check that the planted mode is unchanged).

## 4. Gates (stop on failure)

`devtools::document()`; install; the new test file and the existing `inject`/analytic-layer test files; the export list unchanged. No full R CMD check in this task.

## 5. Commits and notes

Commit on `fs-inject-observed-region` (code + docs + NEWS; tests may be their own commit). Do not push. Write `~/Downloads/NOTES_fs-inject-observed-region_2026-11-02.md`: status, commits, the readers-of-the-record list from §2, test counts, any error text.
