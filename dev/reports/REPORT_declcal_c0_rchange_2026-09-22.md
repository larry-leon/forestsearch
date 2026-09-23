# REPORT — clinically specified null level c0 for the calibrated declaration threshold

Task: `dev/tasks/TASK_declcal_c0_rchange_2026-09-22.md`. Branch `feature/glm-extension`. No compute beyond unit tests; nothing launched; no fetch, pull or push.

## Pins and commits

- The task was written against `81752681`, which was HEAD at the start.
- Step 0 commit (the baseline pin): **`f5798c85`**, `docs(tasks): add c0 null-level task document (2026-09-22)`.
- Implementation: `845fea56`, `feat(declaration): add clinically specified null level c0 to the calibrated declaration threshold`.
- Verification scripts: `c11b52b3`, `docs(verification): commit declaration-calibration post-condition scripts`.
- Report: the commit that adds this file (see OPEN ITEMS).
- Source line numbers below are at the pin `f5798c85` unless marked otherwise.

## R/ call-out

- **Moves existing code:** none.
- **Changes behaviour:** none at any default. Every new formal is `NULL`-default and appended. The four baseline digests are unchanged (PC1).
- **Changes the method:** yes, as a generalization. The calibration's null set is `{ beta(g) <= c0_cmp }` for a pre-specified `c0`. At `c0 = c2` it is the current construction bit for bit (Test 1 / PC5).

## Step 1 — source verification (quoted)

1. **Field helper and its maximum.**
   - The signature is `R/fs_declaration_calibration.R:40`: `.fs_decl_field <- function(db, xi, keep_matrix = TRUE, cor_max = 8L) {`.
   - The maximum is taken at `:55`: `m_star <- zs[cbind(seq_len(nrow(zs)), max.col(zs, ties.method = "first"))]`.
   - `Zstar` is dropped only afterwards, in the return list (`Zstar = if (isTRUE(keep_matrix)) zs else NULL`). So the shifted maximum can be taken from `zs` before `keep_matrix` matters, as D1 requires.
2. **How `c2` reaches the comparison scale** (identified unambiguously; no STOP).
   - `consistency_resample()`, survival: `R/consistency_resample.R:427` `thr_nat <- hr.consistency`, then `:440` `log(thr_nat)  # ratio measure: natural -> log`. The identity branch is `:442` `thr_nat  # identity measure`, and a caller-supplied `comparison_threshold` is used as is (`:438`).
   - The value MR actually receives as `c_cons` comes from the admission set.
     - Survival: `R/forestsearch_main.R:2312` `consistency = log(max(hr.consistency, 0.001))`.
     - GLM ratio measures: `:2189` `consistency_threshold <- log(consistency_threshold)`.
     - GLM identity measures (RD / IRD / MD): the threshold is left as is.
   - The threshold is resolved at `:2369-2372` (`admission_resolved <- .fs_resolve_admission(` … `hr.consistency = threshold_config$consistency,`) and stored verbatim at `R/forestsearch_helpers.R:2493` (`consistency <- list(c_cons = as.numeric(hr.consistency), ...`).
   - **`c0` follows the same map.** On a ratio path (`meta$log_scale = TRUE`), `c0_cmp = log(c0)`; on an identity path, `c0_cmp = c0`.
     - At `c0 = c2` this gives `c0_cmp == c_cons` exactly. `log(c2)` and `log(max(c2, 0.001))` coincide for any `c2 >= 0.001`.
     - `delta_g` is therefore exactly `0`, and Test 1 confirms `identical()` on both paths.
3. **`c_cons` and `sigma_D(g)` at the capture block** (`fs_mr_inference()`).
   - `:623` `B <- asm$B; bh <- asm$beta_hat; sdv <- asm$sigma_D; sz <- asm$sizes`, and `:624` `log_scale <- asm$log_scale`.
   - `:642` `c_cons <- if (.has_cons) admission$consistency$c_cons else NULL`.
   - The capture at `:1078` (`fld <- .fs_decl_field(B, Xi, keep_matrix = isTRUE(keep_field_matrix))`) already stores `c_cons` and `log_scale` into `meta` (`:1090-1091`).
   - `delta_g` is formed there with `sdv`, the denominator of `T(g)`. On the GBSG fit, `sdv` is `identical()` to the field's own `sqrt(colSums(db^2))` (`meta$sigma_D_field`).
4. **Per-`c0` cutoff.** `admitted_calibrated` is `names(bh)[bh >= floor_at(kappa_hat)]` (`:328`), with `floor_at <- function(k)` (`:323`) taking the cutoff as an argument. The per-`c0` block calls the same `floor_at` once per `kappa_hat(c0)`, and nothing else changes.

## Implementation (`845fea56`)

- **`.fs_decl_field(..., shift = NULL)`** (D1).
  - `shift` is a length-G vector or a G×K matrix.
  - When supplied, `Mstar_shift` (B×K) is computed from `zs` before `keep_matrix` is consulted.
  - With `NULL`, the return list is unchanged: the element is not present at all.
- **New `@noRd` helpers** in the same file:
  - `.fs_decl_shifted_max()` takes the maximum, looping over the K columns;
  - `.fs_decl_c0_cmp()` does the scale map and the guards;
  - `.fs_decl_c0_shift()` builds the G×K `delta`;
  - `.fs_decl_c0_block()` assembles the per-`c0` result of `fs_declaration_calibration()`.
- **`fs_mr_inference(..., declaration_c0 = NULL)`** (D2). With `keep_declaration_field = TRUE`, the capture stores `declaration_field$Mstar_c0` (B×K, column names `as.character(c0)`) and `meta$c0` / `meta$c0_cmp`. It is ignored when the field is off.
- **`forestsearch()`** gets one line, `declaration_c0 = .g_mr(mr_inference_args$declaration_c0, NULL),`, plus a roxygen `\item`. It has no new formal.
- **`fs_declaration_calibration(fit, alpha, family, ..., c0 = NULL)`** (D3 / D4).
  - `c0` is appended after `...`, so it is named-only.
  - Where the shifted maxima come from:
    - the capture, when it holds every requested `c0` (pre-reduction family, same `c_cons`);
    - otherwise the stored `Zstar`;
    - otherwise an error that names `declaration_c0` and `keep_field_matrix`. It never falls back to the unshifted maximum.
  - `family = "reduced"` always computes from `Zstar` (D5 unchanged).
  - The result gains `c0 = list(table, admitted_calibrated, Mstar_c0, source)`. `table` has one row per `c0` with these columns: `c0`, `c0_cmp`, `kappa_hat`, `fw_size`, `pstar_implied` (`2 * pnorm(kappa_hat) - 1`), `n_admitted_calibrated`, `Mstar_c0_q90/q95/q99`, `is_c2`.
  - Every pre-existing element is unchanged: with `c0` given, the result minus `$c0` is `identical()` to the `c0 = NULL` result.
  - The print method adds one row per `c0` and labels the `c2` row `(= c2, unshifted)`.
- **Scale guards.**
  - `c0 > c2` (natural scale) errors with `c0 = <c0> exceeds c2 = <c2>`.
  - `c0 <= 0` on a ratio path (a log supplied by mistake) errors, pointing to the natural ratio scale.
  - Non-finite or duplicated `c0` errors.
  - `c0` with no consistency floor (`maxeff`) or no recorded `log_scale` errors.
  - On an identity path a wrong-scale `c0` is undetectable; this is documented in the roxygen.
- Roxygen is markdown. `devtools::document()` regenerated `man/fs_declaration_calibration.Rd`, `man/fs_mr_inference.Rd` and `man/forestsearch.Rd`; `NAMESPACE` is unchanged. `tools::checkRd()` is clean on all three. No new dependency or Suggests.

## Tests

- New file `tests/testthat/test-declaration-c0.R`: **71 expectations, 0 failures, 0 skips.**
- Existing file `test-declaration-calibration.R`: **75 expectations, 0 failures, 0 skips.**
- Both files together took **72.7 s wall** (cap 10 min). They were run with `devtools::load_all()` + `testthat::test_file()` only. No full suite, no R CMD check.
- Fixtures:
  - The GBSG fit is the prerequisite's: `sg_focus = "hr"`, `maxk = 2`, `hr.threshold = 1.1`, `hr.consistency = 1.0`, `p* = 0.90`, `n.min = 60`, `d0/d1.min = 12`, `seedit = 8316951`, MR `draws = 2000`, Poisson, `ci_method = "ij"`.
  - The continuous GLM fit is the prerequisite's (`.make_continuous_data(N = 400, MD_harm = 2)`, `hr.threshold = 0.5`, `hr.consistency = 0.25`, MR `draws = 300`).
  - Configuration B is S1.7, data seed 20260922.
- **Test 1 (reduction): PASS.** On GBSG (`c0 = 1.0`) and on the continuous GLM (`c0 = 0.25`):
  - `Mstar_c0[, as.character(c2)]` is `identical()` to `Mstar`;
  - `kappa_hat`, `fw_size` and `admitted_calibrated` at that row are `identical()` to the `c0 = NULL` results;
  - the result minus `$c0` is `identical()` to the `c0 = NULL` result.
  - The field-matrix route and the capture route also agree `identical()`ly, both on `Mstar_c0` and on the whole table.
- **Test 2 (monotonicity): PASS.** At α 0.05 and 0.10 over `c0 = (0.70, 0.75, 0.80, 0.85, 1.0)`, `kappa_hat` and `fw_size` are non-decreasing and `n_admitted_calibrated` is non-increasing.
- **Test 3 (shift correctness): PASS.** The hand-supplied `db` has G = 4, nested supports, B 3000 and `c0 = (0.70, 0.85, 1.0)`.
  - `delta` matches an independent `(c_cons - log(c0)) / sigma_D` construction.
  - `Mstar_shift[, k]` equals `apply(Zstar - rep(delta_k, each = B), 1, max)` element-wise at tolerance 0.
  - A vector shift equals the one-column matrix shift.
  - With no shift, the element `Mstar_shift` is absent.
- **Test 4 (scale guard): PASS.**
  - `c0 = 1.1` errors with "c0 = 1.1 exceeds c2 = 1".
  - `c0 = log(0.75)` and `c0 = 0` give the "natural ratio scale" error.
  - On the identity path, `c0 = 0.3` against `c2 = 0.25` errors.
  - At capture time through `fs_mr_inference()`, `c0 = 0.5` against `c2 = 0` errors.
  - A fit with neither the capture nor the matrix refuses, naming `keep_field_matrix` and `declaration_c0`.
- **Test 5 (shifted closed form): PASS.** Details are in the Test 5 block.
- **Test 6 (default-off): PASS.**
  - `declaration_c0` and `c0` exist and default to `NULL`.
  - The unset fit has no `Mstar_c0` and no `meta$c0` / `c0_cmp`, and its calibration has no `$c0`.
  - The `c0`-set fit is `identical()` to the unset fit (timing stripped) once four things are removed: `Mstar_c0`, `meta$c0`, `meta$c0_cmp`, and the recorded `args_call_all$mr_inference_args$declaration_c0`.
  - Bit-for-bit identity against the pre-change code is PC1.
- **Test 7 (purity): PASS.** The md5 of both fits is unchanged across calls with `c0` through the capture, through the field matrix, and with `family = "reduced"`.

## Post-conditions (measured)

- **PC1 (baseline identity): PASS.**
  - Method: `dev/verification/postcond_fits.R` was run on a `git archive` snapshot of the pin `f5798c85` and again on the changed tree, with every new formal unset.
  - The digests are identical to each other and to the prerequisite's record:
    - `gbsg_mr_off` `8d2868f0e95e79c0124577620d480feb`;
    - `gbsg_mr_on` `dbbeff93dfc85243b92a22aee872ff1c`;
    - `cont_mr_on` `e033e1b5c1f47c8ce05f7e3770625f2b`;
    - `mr_direct_B` `a676d87f80796547ffc3881891b6642e`.
- **PC2 (formals): PASS.** The 17 functions in `R/fs_mr_inference.R` and `R/forestsearch_main.R` were checked, plus all 5 functions of `R/fs_declaration_calibration.R` at the pin. Only three changed, each by appending one `NULL`-default formal with its existing prefix `identical()`:
  - `fs_mr_inference` + `declaration_c0`;
  - `fs_declaration_calibration` + `c0`;
  - `.fs_decl_field` + `shift`.
  - `forestsearch()` formals are unchanged. Four new `@noRd` helpers were added; they are new functions, not changes to existing ones.
- **PC3 (diff scope): PASS.** `git diff --stat f5798c85` lists only:
  - the three `R/` files, at their named insertion points (+12/−0 `forestsearch_main.R`, +29/−2 `fs_mr_inference.R`, +219/−4 `fs_declaration_calibration.R`);
  - the three regenerated `.Rd` files;
  - the new test file;
  - `dev/verification/postcond_fits.R` and `dev/verification/report_values.R` (the prerequisite's scripts, verbatim);
  - this report;
  - the Step 0 task document.
  - No pre-existing untracked file was staged.
- **PC4 (Test 5 closed-form agreement): PASS** at 0.0027 / 0.025 (below).
- **PC5 (reduction identity on GBSG and continuous GLM): PASS.** On the GLM fit, `meta$c_cons` is `identical()` to `0.25`, `log_scale = FALSE`, G = 118, and the `"0.25"` column is `identical()` to `Mstar`.

## Test 5 — shifted closed-form two-candidate check

- **Configuration:** S1.7 configuration B. g1 = {X = 1} is nested in g2 = {X ≤ 2}; OLS mean difference; `c2 = c_cons = 0`; `p* = 0.90`; **`c0 = -0.05`** (identity scale). Gaussian multipliers, `B = 200000`, MR `seed = 7`.
- **`rho_hat = 0.734879`**, analytic, from the dfbeta matrix. `sigma_D = (0.2213376993, 0.1562950498)`.
- **`delta_1 = 0.225899`, `delta_2 = 0.319908`**.
- `Phi2(a, b; rho) = integrate(dnorm(x) * pnorm((b - rho x) / sqrt(1 - rho^2)), -Inf, a)` is base R, with unequal margins and upper limit `a`.
- **Closed-form targets:**
  - `1 - Phi2(z.95 + d1, z.95 + d2; rho_hat)` = **0.045178**;
  - the root of `Phi2(k + d1, k + d2; rho_hat) = 0.95` = **1.597925**.
  - Both are below the unshifted targets (0.078751 and 1.867648), as they must be.

  | path | fw_size(c0) | discrepancy (tol 0.0027) | kappa_hat(c0) | discrepancy (tol 0.025) |
  |---|---|---|---|---|
  | exported (`fs_mr_inference` capture → `fs_declaration_calibration`) | 0.045225 | +0.000047 | 1.597078 | −0.000847 |
  | internal helper, hand-supplied `db`, independent `xi` (seed 808) | 0.045335 | +0.000157 | 1.601449 | +0.003524 |

## GBSG test fit — c0 table

- 765 candidates in the pre-reduction family; `c2 = 1.0` (`c_cons = 0`); `p* = 0.90` (z 1.6449); Poisson multipliers, B 2000; source `capture`.
- The executed screen admits `{q4.0 & q18.1}` (T̂ 1.7995).

| c0 | kappa_hat_05 | kappa_hat_10 | pstar_implied_05 | fw_size at p* 0.90 | n_admitted_cal05 |
|---|---|---|---|---|---|
| 0.70 | 2.661237 | 2.375797 | 0.992215 | 0.392500 | 0 |
| 0.75 | 2.840237 | 2.528865 | 0.995492 | 0.504000 | 0 |
| 0.80 | 2.991732 | 2.694450 | 0.997226 | 0.617000 | 0 |
| 0.85 | 3.170576 | 2.860775 | 0.998479 | 0.733000 | 0 |
| 1.00 (= c2, unshifted) | 3.623846 | 3.330538 | 0.999710 | 0.940000 | 0 |

- The `c0 = 1.0` row reproduces the prerequisite report's κ̂(0.05) 3.623846, κ̂(0.10) 3.330538 and fw 0.940000.
- The shifts span `delta_g ∈ [0.210, 1.238]` at `c0 = 0.85` and `[0.461, 2.717]` at `c0 = 0.70` (σ_D range 0.1313–0.7730). Candidates with small σ_D, which are large or event-rich, move most.
- **Protecting at HR 0.75 instead of HR 1.0 lowers κ̂(0.05) from 3.62 to 2.84.** Measured against the null at c0, the unchanged p* 0.90 screen's family-wise size falls from 0.94 to 0.50. The calibrated rule still admits nothing on this fit at any c0 in the grid, because T̂ max is 1.80.

## Observations for review (not acted on)

- **GLM ratio paths: the task's scale sentence and the source disagree.** §2 says that on GLM paths `c0` "must be on the same scale as `comparison_threshold`". That argument is already logged for ratio measures (`R/consistency_resample.R:438`). But §2 and §4.2 also say `c0` takes "exactly the transform the package applies to `c2`", and the user-facing `c2` for OR / RR / IRR is natural-scale and logged at `R/forestsearch_main.R:2189`.
  - I implemented the second reading: `c0` is natural-scale on every path, logged iff `log_scale`. That is the only reading under which "`c0 <= c2` on the natural scale" and "at `c0 = c2` the construction reduces exactly" both hold.
  - On identity paths (MD / RD) the two readings coincide.
  - A `c0` supplied as a log on a ratio path is always caught when the intended level is `<= 1`, because its log is `<= 0` and trips the natural-ratio-scale guard. When the intended level is above 1 (possible only when `c2 > 1`), its log is positive and usually below `c2`, so it is not detectable.
- **`declaration_c0` errors inside `forestsearch()` cost the whole MR result.** A bad `c0`, for example `c0 > c2` or `maxeff` with no consistency floor, is raised in the capture block. There it falls into `forestsearch()`'s existing `tryCatch` → `warning("mr_inference failed: ...")` and `mr_inference = NULL`. For a campaign it would be safer to validate before the search. That would need `c_cons` / `log_scale` earlier, and I left it out as out of scope.
- `declaration_c0` without `keep_declaration_field = TRUE` is silently ignored, as D2 specifies. This is documented.
- The capture is reused only when `meta$c_cons` is `identical()` to the `c_cons` in force. In the `maxeff` / `...` case there is no capture to reuse, so the post-hoc route computes from `Zstar` or errors.

## OPEN ITEMS

- The report cannot quote its own commit SHA. `git log -1 -- dev/reports/REPORT_declcal_c0_rchange_2026-09-22.md` gives it.
- The values in the c0 table and the Test 5 block came from a session-scratchpad script (`c0_values.R`). It sources the new test file's fixtures, like the prerequisite's `report_values.R`. It is not committed, because §7.3 widens diff scope for the two prerequisite scripts only. Commit it under `dev/verification/` if wanted.
- `dev/verification/postcond_fits.R` is committed verbatim. Its formals capture covers only `R/fs_mr_inference.R` and `R/forestsearch_main.R`, so PC2 for `R/fs_declaration_calibration.R` was run as a separate one-off comparison against the pin snapshot. The script also calls `git rev-parse HEAD`, which warns (status 128) when run outside a git checkout, as on the archive snapshot. That is harmless.
- The Scale sentence for GLM ratio paths (first observation above) should be confirmed before the campaign task is written for OR / RR / IRR designs.
- Supplement S1.7 configuration B was taken from the prerequisite's resolution (`fs-glms-interpretable/manuscript/fs_glms_interpretable_supplement.tex` §S1.7) and its DGP as copied into `test-declaration-calibration.R`.
