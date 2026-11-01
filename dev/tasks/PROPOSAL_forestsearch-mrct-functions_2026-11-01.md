# PROPOSAL: forestsearch functions for the MRCT consistency work

Date: 2026-11-01
For Larry's decision. Nothing here is a task. Each part names its `R/` category. None changes the method: identification, selection and inference are untouched.

Purpose: a pkgdown article in forestsearch (`vignettes/articles/`) that reproduces the replicate-5 analysis and the simulation results without reading from the mrct-consistency repository, and the Statistics in Medicine paper's code statement. Four parts; the first and fourth are what the article cannot do without.

Sources read: forestsearch `19b703b` (0.4.1): `R/mrct_inject_dgm.R`, `R/simulate_from_dgm.R`, `R/mrct_simulation.R`, `R/betaHhat_truth.R`, `R/fs_oc_family.R`; mrct-consistency `51d3b6f`: `docs/theory/scripts/sk_limit.R`, `ahr_setup.R`, `sims/regional_entry_design.R`, `sims/ideal_design.R`, `sims/regional_entry_cells_design.R`.

## Part A — the analytic layer

**Category: moves existing code** (from `docs/theory/scripts/sk_limit.R`, `ahr_setup.R` and `sims/regional_entry_design.R`), with new exported interfaces.

What it computes, for any patient set of a Weibull AFT DGM in a finite analysis window with staggered entry: the exact limit of the two-sample Cox estimate (the Struthers–Kalbfleisch root), its information and robust variance, the expected events, the correlation between two sets' estimates, the AHR and CDE of the set, and the normal-approximation probabilities of the trigger and the oracle joint at a trial size n. Validated across 30 campaign cells (trigger counts within about 2.4 SE; one systematic miss, the oracle joint understated by about 0.4 points) and the two later-enrolling-region cells (largest |z| 1.68).

Proposed exports (names open):

| Function | Replaces | Interface |
|---|---|---|
| `ahr_population(dgm, analysis_time, max_entry, entry_delay = NULL)` | `ahr_pop()`, `am_population()` | Builds the population object from an `aft_dgm_flex` DGM made by `inject_mrct_structure()`: per-patient log hazards under each arm, H* and region indicators, the Weibull `mu` and `tau`, and each patient's entry window. `entry_delay` is a per-patient vector or `list(var = "z_region", delay = a)`. The core takes per-patient quantities, so a later split into its own package needs no rewrite. |
| `ahr_limit(pop, set, b = NULL, h = NULL, psi = NULL, horizon = NULL, variance = TRUE)` | `sk_limit()`, `re_exact()` | `set` is a name (`"all"`, `"region"`, `"region_Gstar"`, …) or a logical vector. Returns the limit, information, robust variance, event probability, expected events, AHR, CDE, margHR at the horizon and at infinity. Per-patient entry factors throughout, as `sims/regional_entry_design.R` does; with no delay this is the library as it stands (checked to 1e-13 against the committed records). |
| `ahr_design(pop, b, h, n, sets = c("all", "region", "region_Gstar"))` | `design_cell()`, `design_prob()`, `re_design_cell()` | SEs at n, exact correlations, P(ITT ub < 1.0), P(region ITT > 1.0), P(trigger), P(trigger and region(G*) < 1.0). |
| `ahr_solve(pop_fun, target, over, interval, ...)` | `boundary_n1000.R`, `profile_r…R`, `210_regional_entry_analytic.R` | Root of P(trigger) = target over one of `r`, `entry_delay`, `b`, `h`, the others fixed; scan plus `uniroot()`; reports the number of roots. Needs a population builder for `r`, since r enters the build. |

Internal (not exported): the quadrature grid, the Laplace-sum caches, the ladder.

Limits to state in the documentation: Weibull AFT DGMs from `inject_mrct_structure()`; entry uniform per group; administrative censoring only (the random censoring of `cens_adjust` is not in the layer); correlations exact for nested sets.

Tests: the committed records become regression tests — `results/ident_rlog4_full/analytic/predictions.csv` (S2 at r = −log 4), `results/analytic_maps/grid.csv` (b = 0.55, h = 2, r = −log 2), `results/regional_entry/roots.csv` (the two roots), copied into `tests/testthat/` with their commits named.

Dependencies: `mvtnorm` is already imported.

Relation to existing package routes (for later, not this proposal): the survival β(Ĥ) target in `R/betaHhat_truth.R` is a Cox fit to one censored draw of the super-population, and `fs_oc_family_enumerate()` gives each candidate an approximate effect and a scaled SE. The layer could give both exactly. Making it the default there would **change the method** and is not proposed; offering it as an option would not.

## Part B — DGM arguments

**Category: changes behaviour; defaults unchanged.**

1. **Hazard multiplier m.** `sims/ideal_design.R`'s `ideal_set_m()` sets `dgm$model_params$mu <- attr(dgm, "mu_fitted") - tau * log(m)` and an attribute, reaching into the object. Proposed: `scale_dgm_hazard(dgm, m)` exported, doing the same with validation, and recording `m` in `dgm$mrct`. Default m = 1 means no change to any existing build.
2. **Region entry delay.** `mrct_region_sims()` passes a super-population column `entrytime` to `simulate_from_dgm()` as `entry_var` if it exists (`R/mrct_simulation.R`, undocumented). The campaign cells used that hook with one entry draw per super-population row. Proposed: `simulate_from_dgm(..., entry_delay = NULL)` with `entry_delay = list(var = "z_region", delay = a)`: entry = a·1{var} + U(0, `max_entry`), redrawn each trial like the current U(0, `max_entry`); `mrct_region_sims()` passes it through `sim_args` as now. Document the `entrytime` hook. Default `NULL` is the current behaviour.

## Part C — per-replicate output of `mrct_region_sims()`

**Category: changes behaviour; additive, default unchanged.**

`mrct_region_sims()` returns point estimates (`hr_itt`, `hr_train`, `hr_test`, `hr_sg`, counts). Every campaign therefore redrew all 10,000 trials from their seeds to obtain the rest (`sims/123_ident_rlog4_full_records.R`). Proposed: `mrct_region_sims(..., record = "full")` adding, per trial and for each of ITT, non-region, region, region(Ĝ) and region(G*): the Wald 95% limits, the SE of the log HR, patients and events by arm, and the AHR and CDE on the same patients; plus the Ĥ vs H* counts (TP, FP, FN, TN) on the non-region and on the region side. `record = "point"` (default) returns exactly today's columns. With this, the core OC table is a summary of one returned table.

## Part D — the latent correlation

**Category: changes behaviour.** (The proposal deferred on 2026-10-31.)

`expand_covariates()` in `R/mrct_inject_dgm.R` uses latentcor when installed and otherwise falls back silently to a normal-score correlation. The built DGM differs between the two, which is why replicate 5 reproduced on Pop!_OS and not on the Mac. The fitted matrix is already kept in `dgm$mrct$latent_R`.

Proposed:
1. `inject_mrct_structure(..., latent_R = NULL)`: a supplied matrix is used as the latent correlation and no fit is made.
2. When `latent_R` is `NULL` and latentcor is not installed, a warning names the fallback.
3. The article ships the S2 build's `latent_R` (from the committed record, with its commit) and passes it, so a reader on any platform gets the same super-population, the same replicate 5 and the same Ĥ.

## Order

D and A first: D is a prerequisite for reproducibility, and A is the article's overview. Then C, then B. A, B and D can be one task in the forestsearch repository; C is its own, since it touches the main simulation loop.

## What is not proposed

- The Weibull model-based analyses (`sims/160`–`180`), the soft/bagged Ĝ probes (`190`, `200`), truth by 2,000,000-patient draws, the cell and seed bookkeeping, and `get_dfpred()` re-implementations (`get_dfpred()` is exported already).
- Any change to the identifier, the screen, the selection rule or the inference.

## The article (for the record of decisions, 2026-11-01)

- A pkgdown article in `vignettes/articles/` (excluded from the CRAN build by `.Rbuildignore`), built on the 64-core box.
- Order: overview (trigger predictions and the design specifications, computed live by Part A); replicate 5 (live, with the reproduction check); simulation results; appendix on the development of the analytic tools.
- Scope: the n = 1,000 cells N0–S4 and the two later-enrolling-region cells, with the committed summaries shipped in the article's directory with their mrct-consistency commits, the full campaign code shown and not executed, and a small live run of one cell; n = 700 analytic predictions only, computed live (the campaign's committed n = 700 predictions exist; the later-enrolling cells' n = 700 predictions are one `ahr_design()` call each).
