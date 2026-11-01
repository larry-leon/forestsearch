# TASK: forestsearch — the analytic layer, the DGM arguments and the latent correlation (Parts A, B, D)

Date: 2026-11-01
Repository: **forestsearch** (`../forestsearch` on the Pop!_OS box), from `19b703b` (0.4.1). mrct-consistency at `4710058` is read only, as the source of the code being moved and of the records used as tests.
Basis: `PROPOSAL_forestsearch-mrct-functions_2026-11-01.md`, approved 2026-11-01. Part C (the per-replicate output of `mrct_region_sims()`) is a later task.

**`R/` categories:** Part A moves existing code (new exports; no existing function touched). Parts B and D change behaviour through new arguments whose defaults reproduce today's behaviour exactly. Nothing changes the method. CC does not push.

## 0. First action

Copy this file and the proposal (`~/Downloads/PROPOSAL_forestsearch-mrct-functions_2026-11-01.md`) into forestsearch `dev/tasks/` and commit them on their own.

## 1. Backward compatibility (the rule for the whole task)

- No existing export is renamed, removed or given a changed signature. Everything new is a new function or a new argument with a `NULL` or identity default.
- No existing returned object gains, loses or reorders a column or element unless the caller asked for it with a new argument.
- The gates in §6 check this against committed records.

## 2. Part D — the latent correlation (`R/mrct_inject_dgm.R`)

1. `expand_covariates(..., latent_R = NULL)` and `inject_mrct_structure(..., latent_R = NULL)`: a supplied matrix (named rows and columns matching the covariates; checked symmetric and positive definite after the existing `nearPD` step) is used as the latent correlation and no fit is attempted. It is stored in `dgm$mrct$latent_R` as the fitted one is now.
2. When `latent_R` is `NULL` and the latentcor fit is not made (package missing or its call failed), issue a warning naming the fallback (the normal-score correlation) and the reason. The computation is unchanged.
3. Document both in the roxygen of the two functions.

## 3. Part A — the analytic layer (new file `R/ahr_analytic.R`, or a few files)

Move the code of mrct-consistency `docs/theory/scripts/sk_limit.R`, `docs/theory/scripts/ahr_setup.R` (the design-layer part) and `sims/regional_entry_design.R` into the package, under the per-patient-entry formulation of `sims/regional_entry_design.R` (with no delay it is the library as it stands). Exports and interfaces as the proposal gives them:

- `ahr_population(dgm, analysis_time, max_entry, entry_delay = NULL)` — from an `aft_dgm_flex` DGM made by `inject_mrct_structure()`; `entry_delay` a per-patient numeric vector or `list(var = <column>, delay = <months>)`. Checks the DGM's `theta_0`, `theta_1`, `loghr_po` consistency as `ahr_pop()` does.
- `ahr_limit(pop, set, b = NULL, h = NULL, psi = NULL, horizon = NULL, variance = TRUE)` — `set` a name among `"all"`, `"nonregion"`, `"region"`, `"region_Gstar"`, `"Hstar"`, `"Hstar_region"`, `"Hstar_nonregion"`, or a logical vector. Returns what `sk_limit()` returns (limit, information, robust and model variances, event probability, expected events per patient, the types table, ω) together with the set's AHR, CDE and margHR at infinity.
- `ahr_design(pop, b, h, n, sets = c("all", "region", "region_Gstar"))` — SEs at n, exact and SE-ratio correlations, and the probabilities `design_prob()` computes. Name the outputs descriptively (`p_itt_ub_lt1`, `p_region_itt_gt1`, `p_trigger`, `p_trigger_and_Gstar_lt1`); do not expose the region ITT > 0.9 variants.
- `ahr_solve(f, target, interval, ...)` — the scan-and-`uniroot()` pattern of `results/analytic_maps/boundary_n1000.R` and `sims/210_regional_entry_analytic.R`: `f` maps the free parameter to P(trigger); the scan brackets every crossing; returns all roots with the number of roots. Building a population per r is the caller's `f`.

Internal: the quadrature grid, the Laplace-sum cache keyed by (population, horizon), the ladder. Roxygen on every export with the limits stated (Weibull AFT DGMs from `inject_mrct_structure()`; uniform entry per group; administrative censoring only; exact correlations for nested sets), with the derivation references pointing to the technical report (`docs/manuscript/techreport/`, Chapter 6) and the theory notes.

Dependencies: `mvtnorm` (already imported); `pmvnorm(..., algorithm = Miwa(steps = 1024))` as now.

## 4. Part B — DGM arguments

1. `scale_dgm_hazard(dgm, m)` (new export): what `ideal_set_m()` in mrct-consistency `sims/ideal_design.R` does — `model_params$mu <- mu_fitted - tau * log(m)` from the stored fitted value — with validation (`m > 0`; the fitted `mu` present), recording `m` in `dgm$mrct$hazard_multiplier` and the attribute `ideal_set_m()` sets, so existing scripts that read that attribute still work. `m = 1` returns the DGM unchanged except for the record.
2. `simulate_from_dgm(..., entry_delay = NULL)`: with `entry_delay = list(var = <column>, delay = a)`, entry = a·1{var = 1} + U(0, `max_entry`), redrawn each call as U(0, `max_entry`) is now; with a numeric vector of length n, used as given. `entry_var` keeps precedence if supplied (document the `entrytime` convention `mrct_region_sims()` uses). `mrct_region_sims()` passes `entry_delay` through `sim_args` as it passes everything else; no change needed there beyond documentation.

## 5. Tests (`tests/testthat/`)

Fixtures copied from mrct-consistency with their commits named in a `fixtures/README`:
- `results/ident_rlog4_full/analytic/predictions.csv` rows for `S2_rmlog4`;
- `results/analytic_maps/grid.csv` row b = 0.55, h = 2, n = 1,000, r = −log 2;
- `results/regional_entry/roots.csv`;
- `results/ident_rlog4_full/summary/build.csv`.

Tests:
1. **Build identity.** The S2 DGM built with `inject_mrct_structure()` as mrct-consistency `sims/ident_rlog4_full_design.R` → `build_full("S2")` builds it (replicate the call with its arguments in the test, using `scale_dgm_hazard()` for m = 3) matches the committed `build.csv` row for `S2_rmlog4` in every numeric column (tolerance 1e-8). Skipped, with a message, when latentcor is not installed.
2. **The layer at no delay.** `ahr_design()` on that population reproduces the committed S2 predictions and the `grid.csv` row (probabilities to 1e-9; events to 1e-9 relative), as `sims/210_regional_entry_analytic.R`'s checks did.
3. **The layer with a delay.** `ahr_design()` at the two committed roots reproduces `roots.csv` (same tolerances); `ahr_limit(..., set = "region")` at delay a equals `ahr_limit()` at horizon 40 − a with no delay (1e-6).
4. **`entry_delay`.** With a fixed seed, `simulate_from_dgm(entry_delay = NULL)` equals the current call element for element; with a delay, region entry times lie in [a, a + `max_entry`] and non-region in [0, `max_entry`].
5. **Part D.** Supplying the stored `latent_R` reproduces the build; the fallback warns once when latentcor is unavailable (test by mocking `requireNamespace`).

## 6. Gates (stop on failure)

1. `devtools::document()` clean; `devtools::install()`; the export list grows by exactly `ahr_population`, `ahr_limit`, `ahr_design`, `ahr_solve`, `scale_dgm_hazard`, and nothing is removed.
2. The new tests pass (`testthat::test_file()` on each new file).
3. **Legacy reproduction, from mrct-consistency with the new package installed:**
   - `source("sims/ident_rlog4_full_design.R"); build_full("S2")` matches the committed `build.csv` row (the design file's own `full_build_record()`), with no code change in mrct-consistency;
   - `Rscript sims/210_regional_entry_analytic.R` into a scratch directory (`RE_OUT`) reproduces the committed `results/regional_entry/roots.csv` byte for byte on its numeric columns (to 1e-9).
4. No full `R CMD check` in this task; it is its own task at the end of the series.

## 7. NEWS

A `NEWS.md` entry under a development version (bump to 0.4.1.9000): the five new exports, the two new arguments, the new warning.

## 8. Notes file

`~/Downloads/NOTES_fs-mrct-analytic-layer_2026-11-01.md`: status, commit hashes, the export diff, test counts, the gate-3 differences (largest), wall-clock, and any place the source disagreed with this task. Do not paste output back.

## 9. Commits

Three:
1. the task document and the proposal (step 0);
2. Part D and Part B;
3. Part A with its tests, NEWS and documentation.
