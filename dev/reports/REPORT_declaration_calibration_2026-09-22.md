# REPORT — calibrated declaration threshold (κ̂_α) and family-wise size diagnostic (FŴ_α)

Task: `dev/tasks/TASK_declaration_calibration_2026-09-22_v2.md` (v2; v1 not executed).
Branch `feature/glm-extension`. No compute beyond unit tests; nothing launched; no push.

## Pins and commits

- HEAD at start: `be3d5a92` (not `c4baf79`; `git diff c4baf79 be3d5a92 -- R/` is empty, so the task's quoted pointers were re-verified from source, below).
- Step 0 commit (the baseline pin): `8a318c10`, `docs(tasks): add declaration-calibration task document v2 (2026-09-22)`.
- Implementation commit: `c5e9f91e`, `feat(declaration): add opt-in calibrated declaration threshold and family-wise size diagnostic`.
- Report commit: the commit that adds this file (see OPEN ITEMS).
- Line numbers below are at the pin `8a318c10`.

## Branch taken: B

- The fitted object does **not** retain the shared multipliers, the per-candidate `dfbeta`, or the per-candidate `beta_hat` / `sigma_D`. `fit$mr_inference` carries only summaries (`selected_index`, `naive`, `debiased`, `settings`, `reselection$p_hat`, ...). Its settings record the law and `draws` but not the seed.
- The seed is recoverable (`mr_inference_args$seed`, default `seedit`), but regenerating `D[b, g]` from it would need the `dfbeta` refits, which D5 forbids. So Branch A is not achievable.
- The pre-reduction family *is* reachable at the capture point without behaviour change. `forestsearch()` hands `fs_mr_inference()` its own enumeration of `Z` (`R/forestsearch_main.R:3713-3726`), which is built independently of `remove_near_duplicate_subgroups()`. The shared multiplier matrix `Xi` is in scope there. No STOP condition was met.
- There are two insertion points:
  - `R/fs_mr_inference.R` gains the formals `keep_declaration_field = FALSE` and `keep_field_matrix = FALSE`, plus one capture block after `out` is complete.
  - `R/forestsearch_main.R` gains two `.g_mr(mr_inference_args$..., FALSE)` pass-throughs (no new formal on `forestsearch()`) and a roxygen item.
- The capture stores `declaration_field = list(Mstar, Zstar, beta_hat, sigma_D, family_id, meta)`. `Zstar` is kept only under `keep_field_matrix = TRUE`.

## Step 1 findings (quoted at the pin)

- **Closed screen.** `R/subgroup_consistency_helpers.R:1463-1465`: `rr <- consistency_resample(` / `df.x, hr.consistency = hr.consistency, method = "closed",` / `adjust_covariates = adjust_covariates, cox_init = cox_init`. The GLM branch is at `:1486`. `consistency_method = "resample"` is the `forestsearch()` default.
- **The closed form.** `R/consistency_resample.R:453`: `out$rate_closed <- max(0, 2 * stats::pnorm(delta / pieces$sigma_D) - 1)`, with `sigma_D <- sqrt(sum(dfbeta^2))` at `:86`.
- **c2 (`c_cons`).** Carried by `hr.consistency` in the signature `consistency_resample <- function(df, hr.consistency = 1.0, consistency_threshold = NULL, comparison_threshold = NULL, ...)` (`:401-403`). It is natural-scale on the survival path and logged inside (`:440`, `log(thr_nat)  # ratio measure: natural -> log`). GLM passes `comparison_threshold`, already on the comparison scale.
- **c1 (`c_screen`).** Carried by `hr.threshold`. The consistency stage filters on it at the natural scale for survival (`R/subgroup_consistency_main.R:557`, `found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold, ]`). The admission set holds both floors on the comparison scale (`R/forestsearch_main.R:2304-2305`: `screening = log(hr.threshold)`, `consistency = log(max(hr.consistency, 0.001))`), resolved at `:2362` (`admission_resolved <- .fs_resolve_admission(`).
- **MR replay: event and truncation knobs.** `R/forestsearch_main.R:3707-3712`: "the per-arm event minima (d0.min/d1.min) and max_subgroups_search are NOT replayed here, so this family is a superset of the one the identifier actually chose among".
- **MR replay: size minimum.** Applied once, at the observed data, when the family is enumerated (`if (length(mem) >= n.min)`).
- **MR replay: fixed `t_g`.** `R/fs_mr_inference.R:628`: `t_g <- pmax(admission$effect_floor, c_cons + z * sdv)`. It is computed once from the observed `sigma_D`. Then `:654-656`: `Xi <- .fs_mr_multipliers(nrow(B), draws, multiplier)` / `P  <- crossprod(B, Xi)` / `beta_star <- bh + P`. So `beta_star` moves per draw while `t_g` stays fixed.
- **`remove_near_duplicate_subgroups()` call site.** `R/subgroup_consistency_main.R:603`: `found.hrs <- remove_near_duplicate_subgroups(found.hrs, details = details)` (non-`maxeff` branch; `maxeff` uses `.maxeff_membership_dedup()` at `:595`). Its position:
  - (i) after candidate enumeration and after the search's own filters;
  - after the c1 filter at `:557`;
  - (ii) **before** the consistency screen (sort `:616`, truncation `:649`, then evaluation);
  - (iii) **outside** the MR replay, which enumerates its own family from `Z` and never sees the reduction.
- **What the fit retains** (item 5):
  - (i) `xi`: no, only the seed.
  - (ii) `db` / `D`: no.
  - (iii) `beta_hat` / `sigma_D`: no.
  - (iv) pre-reduction membership: not stored, but the MR family is built from `Z` at the capture point. The **reduced** family is recoverable post hoc from the retained `find.grps$out.found$hr.subgroups` by replaying the same filter and the same `remove_near_duplicate_subgroups()` call.
- **Multiplier law settable at the capture site: yes.** `multiplier = .g_mr(mr_inference_args$multiplier, "poisson")` (`R/forestsearch_main.R:3765`); `fs_mr_inference()` accepts `"gaussian"`. So Test 8 ran through the production assembly path.
- **Survival vs GLM element names: identical.** `hr.subgroups` carries `grp, K, n, E, d1, m1, m0, HR, L(HR), U(HR)` plus the factor indicators on both paths. The only difference is scale: the survival `HR` column is natural and the GLM column is on the comparison scale. The replay reads `outcome_type` for that. The MR capture is generic in `B` / `bh` / `sdv`. On the continuous GLM fit: replay check TRUE, 0 unmatched, the relabelled rule reproduces the admitted set, family 118 / 118.

## Implementation

- `R/fs_declaration_calibration.R` (new) contains:
  - exported `fs_declaration_calibration(fit, alpha = 0.05, family = c("prereduction", "reduced"), ...)` and `print.fs_declaration_calibration()`;
  - the internal `.fs_decl_field(db, xi, keep_matrix, cor_max = 8)`, the single place the column standardization lives, callable without a fit;
  - `.fs_decl_reduction()`, the post-hoc replay of the near-duplicate reduction;
  - `.fs_decl_key()`, which makes label matching independent of order.
- `κ̂_α` is the empirical `type = 1` quantile. `fw_size` uses the fit's own `p_star`.
- `admitted_current` is read from `out_sg$result`. `admitted_calibrated` and `admitted_pstar` are over the chosen family. `screened` is the set the consistency screen evaluated.
- `...` reads only `p_star` / `c_cons`, and only when the admission set has no consistency floor (`maxeff`). If they are needed and missing, it errors, naming them.
- No new dependency or Suggests (digests in the tests use `serialize()` + `tools::md5sum()`). `devtools::document()` regenerated `NAMESPACE`, `man/fs_declaration_calibration.Rd`, `man/fs_mr_inference.Rd` and `man/forestsearch.Rd`.

## Tests (`tests/testthat/test-declaration-calibration.R`)

- `devtools::load_all()` + `testthat::test_file()` only: **75 expectations, 0 failures, 0 skips, 35.5 s wall** (cap 10 min). No full suite, no `R CMD check`, no vignette build.
- Test fit: GBSG with `sg_focus = "hr"`, `maxk = 2`, `hr.threshold = 1.1`, `hr.consistency = 1.0`, `p* = 0.90`, `n.min = 60`, `d0/d1.min = 12`, `use_twostage = FALSE`, `seedit = 8316951`. MR at `draws = 2000`, `"poisson"`, `ci_method = "ij"`. `hr.threshold = 1.1` was chosen because at 1.25 the reduction removes nothing (10 → 10); at 1.1 it removes one (31 → 30).
- **Test 1 (relabelling exactness): PASS.** The rule `beta_hat >= max(c_screen, c_cons + z_pstar * sigma_D)` over the 30 screened candidates gives `{q4.0 & q18.1}`, the executed screen's admitted set exactly (`T̂ = 1.7995` vs `z = 1.6449`; runner-up 1.5624). The replayed screened set equals the executed count (`replay_check = TRUE`, 0 unmatched).
- **Tests 2-4 (single candidate, independent columns, nested monotonicity): PASS.** These use Gaussian synthetic fields at `B = 20000`, with 4-MC-sd tolerances.
- **Test 5 (conservatism): PASS.** One candidate was removed (`q15.1 & q18.1`). `κ̂(pre) >= κ̂(red)` at α 0.05 and 0.10. The reduced print carries the conditional-on-realized-family caveat.
- **Test 6 (purity): PASS.** The md5 of the fit is unchanged across both family variants.
- **Test 7 (default-off): PASS.**
  - Both formals default `FALSE`, and the unset fit has no `declaration_field`.
  - The field-on MR result, minus `declaration_field`, is `identical()` to the field-off result (timing stripped).
  - The unset call errors naming `keep_declaration_field`, and `family = "reduced"` without the matrix errors naming `keep_field_matrix`.
- **Test 8: PASS.** It runs through both the exported path and the helper on a hand-supplied `db`; details in the Test 8 block.

## Post-conditions (measured values)

- **PC1 (baseline identity): PASS.** Timing-stripped md5s of four fixed-seed objects were identical before (pin `8a318c10`, run twice, deterministic) and after, with the new formals unset:
  - `gbsg_mr_off` `8d2868f0…`;
  - `gbsg_mr_on` `dbbeff93…`;
  - `cont_mr_on` (GLM continuous, MR on) `e033e1b5…`;
  - `mr_direct_B` (`fs_mr_inference()` direct, OLS) `a676d87f…`.
- **PC2 (formals): PASS.** Over the 17 functions in the two touched files, the only change is `fs_mr_inference`, which gains `keep_declaration_field, keep_field_matrix` (both `FALSE`) appended after the existing formals. The existing prefix (names, defaults, positions) is `identical()`. `forestsearch()` formals are unchanged.
- **PC3 (diff scope): PASS.** `git diff --stat 8a318c10` plus new files lists only:
  - the new files `R/fs_declaration_calibration.R`, `tests/testthat/test-declaration-calibration.R` and `man/fs_declaration_calibration.Rd`;
  - the generated `NAMESPACE`, `man/fs_mr_inference.Rd` and `man/forestsearch.Rd`;
  - the insertion points `R/fs_mr_inference.R` (+44) and `R/forestsearch_main.R` (+12);
  - this report.
- **PC4 (scaling): PASS.**
  - GBSG (Poisson, B 2000, G 765): `column_sd` range **[0.9548, 1.0339]**; `mean(Zstar)` **−0.00527** (tolerance 4/√B = 0.0894).
  - OLS (Gaussian, B 200000): `column_sd` **0.99974, 1.00161**; mean **0.000490**.
  - No scaling-convention mismatch.
- **PC5 (shared multipliers): PASS.**
  - `Xi` regenerated from the recorded seed (`set.seed(7)`, `.fs_mr_multipliers`) reproduces the exported path's `Mstar` **identically** through the helper.
  - Element-wise, `Zstar[b, g] == sum(xi[, b] * db[, g]) / sigma_D(g)` holds.
  - `ncol(Zstar)` = 765 = `n_family_prereduction`, and the reduced column count = 764 = `n_family_reduced`.
- **PC6 (σ_D provenance): PASS.**
  - GBSG selected candidate: field 0.3021782247 vs `consistency_resample(method = "closed", cox_init = <screen's warm start>)` 0.3021782247 (difference −2.2e-16; MR uses `cox_init = 0`, the screen a warm start).
  - OLS: g1 0.2213376993 = 0.2213376993, and g2 0.1562950498 = 0.1562950498.
- **PC7 (family size ordering): PASS.** 765 ≥ 764 (GBSG); 118 ≥ 118 (GLM).
- **PC8 (Test 8 closed form): PASS.** See below.

## Test 8 — closed-form two-candidate check

- **Configuration:** supplement S1.7 (`fs_glms_interpretable_supplement.tex`, "The field on a two-candidate example in closed form"), **configuration B**. g1 = {X = 1} (share 1/4) nested in g2 = {X ≤ 2} (share 1/2); OLS mean difference; k = 0 (the tie); σ = 1; n = 400; data seed 20260922. The DGP is copied from `fs-glms-interpretable/dev/verification/field_two_candidate/field_two_candidate_sim.R::gen_data`. Configuration A was not used because its disjoint candidates give ρ = 0.
- **Admission passed to MR:** effect floor NULL, consistency `{c_cons = 0, p_star = 0.90}`, `reselection = "maxeff"`, `seed = 7`, `B = 200000`.
- **ρ̂ = 0.734879**, computed analytically from the dfbeta matrix. The population value is √0.5 = 0.7071. The test is non-degenerate.
- **Closed-form targets:** `rho_hat` 0.734879; `1 − Φ₂(z.95, z.95; ρ̂)` **0.078751**; root of `Φ₂(k, k; ρ̂) = 0.95` **1.867648**.
- **Gaussian multipliers, exported path (the gate):**

  | check | realized | discrepancy | tolerance |
  |---|---|---|---|
  | 8a `field_cor[1,2]` | 0.735304 | +0.000424 | 0.00894 |
  | 8b `fw_size` | 0.078990 | +0.000239 | 0.0027 |
  | 8c `κ̂(0.05)` | 1.869960 | +0.002313 | 0.025 |

- **Gaussian, internal helper on a hand-supplied `db`** (independent `xi`, seed 808): discrepancies +0.000884, +0.000629 and +0.002534. PASS.
- **8d:** `fw_size` is in [0.05, 0.0975] and `κ̂` is in [1.644854, 1.954508]. The closed-form targets decrease strictly in ρ on the grid, and the task's reference table is reproduced to 1e-5 / 1e-6.
- **Production law, recorded as a finding and not gated** (centred Poisson, same data, B 200000): realized cor 0.733366 / fw 0.078330 / κ̂ 1.864597. Discrepancies are **−0.001514 / −0.000421 / −0.003050**.
  - All three are within the Monte Carlo band that applies at this B.
  - At n = 400 the CLT approximation of the Poisson-multiplier field is not separable from MC error.

## κ̂ and fw_size on the fixed-seed test fit (GBSG, Poisson, B 2000)

| family | n family | κ̂(0.05) | κ̂(0.10) | fw_size at p* 0.90 (z 1.6449) |
|---|---|---|---|---|
| prereduction (default) | 765 | **3.623846** | **3.330538** | **0.940000** |
| reduced (conditional on realized family) | 764 | 3.623846 | 3.330538 | 0.940000 |

- The gap is exactly zero on this fit. The one removed near-duplicate never moves the maximum at the quantiles or at 1.645.
- Mstar quantiles (type 1): 50% 2.4886, 90% 3.3305, 95% 3.6238, 99% 4.1202.
- **The p* = 0.90 screen has family-wise size 0.94 on this 765-candidate family.**
- The executed rule admits `{q4.0 & q18.1}` (T̂ 1.7995). The calibrated rule admits **none** of the 765 at α 0.05. Reported only; admission is unchanged.

## Observations for review (not acted on)

- **Relabelling exactness is exact up to `pconsistency.digits`.** The screen compares `round(rate_closed, pconsistency.digits) >= p_star` (`R/subgroup_consistency_helpers.R:1468`; `pconsistency.digits` defaults to 2 in `subgroup.consistency()` and `forestsearch()` does not pass it). A candidate with closed-form rate in [0.895, 0.90) is admitted by the screen but not by `T >= z`. That band did not arise on the test fit or the GLM fit.
- **The pre-reduction family is MR's enumeration.** It excludes nothing on the per-arm event minima (d0/d1 are not replayed). So it is a superset of the family the search's screen saw: 765 enumerated with n ≥ 60 vs 596 passing the search's filters. That superset is the conservative direction D2 asks for, but it is also why the calibrated threshold is high here.
- **`maxeff` fits** have no consistency floor, so `fs_declaration_calibration()` refuses unless `p_star` / `c_cons` are supplied. Their reduced diagnostic reports "not applicable" (membership de-duplication, §5 item 2, out of scope).

## OPEN ITEMS

- The report location was chosen as `dev/reports/`, where the recent package-level R-change reports live (directive A/C, threshold_sync, estimability boundary). `REPORT_*` files are spread over eleven directories, and the MR-field reports of 2026-09-05/06 went to `quarto/simulations/gbsg_020/`. Move it if gbsg_020 is preferred.
- The report cannot quote its own commit SHA. `git log -1 -- dev/reports/REPORT_declaration_calibration_2026-09-22.md` gives it.
- The PC1/PC2 script (`postcond_fits.R`) and the report-values script ran from the session scratchpad and are not committed; PC3 limits the diff to the listed paths. The fit configurations are reproduced in the test file and in this report.
- "Supplement S1.7" was resolved to `fs-glms-interpretable/manuscript/fs_glms_interpretable_supplement.tex` §S1.7 by section numbering. The task named no file.
- The governing handoff (`HANDOFF_declaration_calibration_forestsearch_2026-09-21_v3.md`) was not found locally and was not consulted. The task document is self-contained, as it states.

## Appendix — zero-compute costing read (Step 6)

Sources: `quarto/simulations/gbsg_020/current_status.md` §2.8-2.10; `REPORT_null_gbsg_identification_2026-09-21.md`; `REPORT_null_gbsg_thresholds_2026-09-21.md` (Table 5). Measured values only.

- **`nullid` (c1 0.90 / c2 0.80, p* 0.90):**
  - Grid: 6 cells (uniform HR 0.657 / 0.721 × n 500 / 1000 / 1500) × 3 identifiers × **2,000 replicates** per cell-run, so 18 runs and 36,000 replicate-searches.
  - Machine: Mac-Studio-3.local, **12 workers**, R 4.5.2, forestsearch 0.3.5.9000.
  - **Measured wall: 6,909 s of render (1.919 h)** over the 18 runs.
  - FS runs alone: 316 + 382 + 320 + 399 + 329 + 394 = **2,140 s**.
  - MR off.
- **`nullc125` (c1 1.25 / c2 1.00, p* 0.90):**
  - Same 6 cells × 3 identifiers × **2,000 replicates**, 18 runs.
  - Machine: pop-os, **64 workers**, R 4.6.1, forestsearch 0.3.5.9000.
  - **Measured wall: 4,002 s of render (1.112 h)**; driver 4,005 s end to end.
  - FS runs alone: 186 + 216 + 190 + 200 + 196 + 196 = **1,184 s**.
  - MR off.
- The two walls confound machine, worker count and screen, so they are not a screen comparison (the thresholds report says so).
- For pricing an MR-on rerun: **`nullmr`** (nullid's screen, MR on, same grid, pop-os, 64 workers) measured **19,858 s (5.52 h)**. MR ran on 21,014 / 21,014 declarations, i.e. only on declaring replicates.
- **Is `nullid`'s per-replicate `maxT` sufficient to recompute `Mstar`? No. The multiplier draws must be taken afresh.**
  - `maxT` is `max_g log(HR_g) / se_g` over `out.found$hr.subgroups`, with `se_g` from the Wald interval. It is:
    - an observed statistic, not a draw of the perturbation field;
    - centred at 0, not at `c_cons` (log 0.80 there);
    - scaled by the model-based SE, not the robust `sigma_D`;
    - taken over the floor-cleared family (median 5-100 candidates), not the pre-reduction family (~1,711-1,830 enumerated).
  - `Mstar` needs every candidate's `dfbeta` under one shared multiplier stream per replicate. Neither `nullid` nor `nullmr` retained it (the capture did not exist).
  - The recorded `maxT` is not even the observed side of the calibrated comparison, which needs `T(g) = (β̂ − c_cons) / σ_D` over the family. The retained `declaration_field` supplies both sides.
- **Scope note for that evaluation.** At α = 0.05 and p* = 0.90, `κ̂ ≥ z_0.95` (a maximum dominates any single column), so the calibrated rule admits a subset of the current rule's candidates. MR on declaring replicates then suffices, as `nullmr` ran it.
  - At α = 0.10, `κ̂` can fall below 1.645 on small or highly correlated families. The calibrated rule could then declare where the current one did not, which needs the field on non-declaring replicates too. `forestsearch()` runs MR only when a subgroup is identified.
