# TASK — Clinically specified null level c0 for the calibrated declaration threshold

**Repo:** `larry-leon/forestsearch`
**Branch:** `feature/glm-extension`
**Design written against:** `81752681` (2026-09-22)
**Prerequisite:** the declaration-calibration implementation at `c5e9f91e` (`fs_declaration_calibration()`, the
`keep_declaration_field` / `keep_field_matrix` capture formals on `fs_mr_inference()`, the internal field helper).
**Kind:** add-only R change generalizing the null the calibration protects against, plus its tests and report.
**Compute:** NONE. Unit tests only. The campaign that uses this is a separate task with its own go.
**Authorization:** Larry proposed this generalization (2026-09-22) and approves the R change on this document.

---

## 0. First action

1. Copy this file verbatim to `dev/tasks/TASK_declcal_c0_rchange_2026-09-22.md`; `git add` that explicit path;
   commit `docs(tasks): add c0 null-level task document (2026-09-22)`.
2. Record `git rev-parse HEAD` (the baseline pin) in the report.

No `git fetch` / `pull` / `push`. Every `git add` names explicit paths. Never stage a pre-existing untracked file.

---

## 1. R/ call-out (required classification)

- **Moves existing code:** NONE.
- **Changes behaviour:** NONE at any default. Every new formal defaults to `NULL`, and with `NULL` the executed
  path and the returned object are identical to the baseline — asserted bit-for-bit in §6.
- **Changes the method:** YES, as a generalization with the current method as the special case. The null set
  the calibration controls over becomes `{ beta(g) <= log(c0) for every candidate }` for a pre-specified `c0`,
  instead of `{ beta(g) <= c_cons }`. At `c0 = c2` it is the current construction exactly. This is Larry's
  proposal; the justification is §2.

---

## 2. What is being added — self-contained specification

Three thresholds, not two:

- `c1`, `c2` — **the claim.** Unchanged. A candidate is admitted under the relevance floor and the consistency
  screen at `c2`.
- `c0` — **what the claim is protected against.** A clinically specified, **pre-specified** benefit level
  (for example HR 0.75). The false-declaration rate is controlled at `alpha` whenever every candidate's true
  effect is at least as good as `c0`.

With the existing standardized field `Zstar[b,g] = D_g(b) / sigma_D(g)`:

```
delta_g        = ( c_cons - log(c0) ) / sigma_D(g)          # on the comparison scale; >= 0 when c0 <= c2
Mstar_c0[b]    = max_g { Zstar[b,g] - delta_g }
kappa_hat(c0)  = empirical (1 - alpha) quantile of Mstar_c0
fw_size(c0)    = mean( Mstar_c0[b] > qnorm((1 + p_star)/2) )
```

Admission under the calibrated rule is unchanged in form — `T(g) >= kappa_hat(c0)` with
`T(g) = (beta_hat(g) - c_cons) / sigma_D(g)` — only the constant changes.

Why it is correct: if every candidate's true effect is `log(c0)`, then `T(g)` is centred at `-delta_g`, so the null
law of `max_g T(g)` is the law of the shifted maximum. Strong control over `{beta(g) <= log(c0)}` follows by the
same stochastic-ordering argument as the current Proposition 5, with `log(c0)` in place of `c_cons`. At
`c0 = c2`, `delta_g = 0` and everything reduces to the current construction.

**Scale.** `c0` is supplied on the same scale as `hr.consistency` (natural HR on the survival path) and transformed
to the comparison scale by exactly the transform the package applies to `c2` (§4 item 2 verifies it from source).
On GLM paths `c0` must be on the same scale as `comparison_threshold`. **Require `c0 <= c2` on the natural scale
(`delta_g >= 0`)**; error informatively otherwise — a protected level worse than the harm criterion is not a null.

**`c0` may be a vector.** Then `Mstar_c0` is a `B x K` matrix with column names the `c0` values, and every returned
quantity is indexed by `c0`.

---

## 3. Design decisions (settled)

- **D1 — Where the shift lives.** In the single internal helper that assembles and standardizes the field (the
  one place standardization already lives). New argument `shift = NULL` (numeric length-G vector, or `G x K`
  matrix); `NULL` means no shift and the current behaviour. The shifted maximum is computed **before** any
  `keep_field_matrix` decision, so it is available at capture time without storing the matrix.
- **D2 — Capture-time computation for campaigns.** `fs_mr_inference()` gains `declaration_c0 = NULL`. When
  non-`NULL` and `keep_declaration_field = TRUE`, the capture stores `declaration_field$Mstar_c0` (`B x K`, column
  names = `c0`) alongside the existing unshifted `Mstar`. `forestsearch()` passes it through `mr_inference_args`
  as the existing two formals are; no new formal on `forestsearch()`.
- **D3 — Post-hoc computation for single-trial use.** `fs_declaration_calibration(fit, alpha, family, c0 = NULL,
  ...)`: when `c0` is given, use `declaration_field$Mstar_c0` if the fit carries it for those `c0`; otherwise
  compute from the stored field matrix (needs `keep_field_matrix = TRUE`); otherwise error naming what to set.
  Never silently fall back to the unshifted maximum.
- **D4 — Return shape.** For each `c0`: `kappa_hat`, `fw_size`, `pstar_implied`, `admitted_calibrated`,
  `n_admitted_calibrated`, `Mstar_c0` quantiles at 0.90 / 0.95 / 0.99. The print method shows one row per `c0`,
  with the unshifted (`c0 = c2`) row labelled as such when present.
- **D5 — Family convention unchanged:** the maximum is over the pre-reduction family. `family = "reduced"`
  remains diagnostic-only.
- **D6 — Nothing else changes.** No change to the screen, the correction, the bounds, any default, or the
  rounding convention. `sg_focus` is **not** re-run on the calibrated set (separate pending proposal).

---

## 4. Step 1 — Source verification (read-only, quoted lines in the report)

1. The internal field helper's signature and the line where the maximum is taken.
2. **Exactly how `c2` is transformed to the comparison scale** on the survival path
   (`consistency_resample()`, expected `log(thr_nat)`) and what scale GLM paths use — apply the identical
   transform to `c0`.
3. Where `c_cons` and `sigma_D(g)` are available at the capture block in `fs_mr_inference()`, so `delta_g` can
   be formed there.
4. That `fs_declaration_calibration()`'s existing `admitted_calibrated` logic can take a per-`c0` cutoff.

STOP if the comparison-scale transform for `c2` cannot be identified unambiguously from source.

---

## 5. Step 2 — Implementation

- Helper: `shift` argument as D1.
- `fs_mr_inference()`: `declaration_c0 = NULL` as D2; capture block stores `Mstar_c0`.
- `forestsearch()`: one `.g_mr(mr_inference_args$declaration_c0, NULL)` pass-through and a roxygen item, in the
  same style as the existing two.
- `fs_declaration_calibration()`: `c0 = NULL` as D3/D4; print method extended.
- Roxygen markdown ON; literal `%`, `<`, `>`, `&`; `@section` titles plain. `devtools::document()`. No new
  dependency, no new Suggests. Vectorized; no growing objects.

---

## 6. Step 3 — Tests (`tests/testthat/test-declaration-c0.R`; run ONLY this file and the existing
`test-declaration-calibration.R`; hard cap 10 minutes; no R CMD check, no full suite)

1. **Reduction to current.** With `c0 = c2` (natural scale), `kappa_hat`, `fw_size` and `admitted_calibrated` are
   `identical()` to the `c0 = NULL` results, and `Mstar_c0[, "c2"]` is `identical()` to `Mstar`.
2. **Monotonicity in `c0`.** On the fixed-seed GBSG test fit with `c0 = c(0.70, 0.75, 0.80, 0.85, 1.0)`:
   `kappa_hat` is non-decreasing in `c0` and `fw_size` non-decreasing in `c0`.
3. **Shift correctness.** On a hand-supplied `db`, `Mstar_c0` equals `apply(Zstar - rep(delta, each = B), 1, max)`
   computed independently, element-wise.
4. **Scale guard.** `c0 > c2` errors with a message naming both; `c0` on the wrong path scale is documented and,
   where detectable, errors.
5. **Closed-form two-candidate check, generalized (Test 8 with a shift).** Supplement S1.7 configuration B,
   Gaussian multipliers, `B = 200000`, as in the prerequisite. With shifts `delta_1, delta_2` (from a chosen `c0`
   and the two candidates' `sigma_D`), the null law of the shifted maximum gives:
   - `fw_size(c0)` equals `1 - Phi2(z + delta_1, z + delta_2; rho_hat)` with `z = qnorm(0.95)`;
   - `kappa_hat(c0)` equals the root of `Phi2(k + delta_1, k + delta_2; rho_hat) = 1 - alpha`;
   where `Phi2(a, b; rho) = integrate(function(x) dnorm(x) * pnorm((b - rho * x) / sqrt(1 - rho^2)), -Inf, a)`
   — base R, unequal margins, **upper limit `a`**. Tolerances 0.0027 and 0.025 as in the prerequisite. This
   test is **required**, not optional; it is the sharpest check available.
6. **Default-off.** A fit taken with every new formal unset is `identical()` to the baseline fit
   (timing stripped), and carries no `Mstar_c0`.
7. **Purity.** The fit is unchanged by `fs_declaration_calibration(fit, c0 = ...)`.

---

## 7. Step 4 — Post-conditions (each failure is a STOP)

1. Baseline identity: the four fixed-seed digests of the prerequisite (`gbsg_mr_off`, `gbsg_mr_on`, `cont_mr_on`,
   `mr_direct_B`) are unchanged with all new formals unset.
2. Formals: the only differences are the added `NULL`-default formals, appended; nothing existing changed.
3. Diff scope: the new test file, the touched `R/` files at their named insertion points, regenerated
   `.Rd`/`NAMESPACE`, the report, and **the two verification scripts from the prerequisite task**
   (`postcond_fits.R` and the report-values script), which are to be committed this time under
   `dev/verification/` — the diff-scope rule is widened for exactly those two files and no others.
4. Test 5 closed-form agreement at the stated tolerances.
5. Reduction identity (Test 1) holds on the GBSG fit **and** the continuous GLM fit.

---

## 8. Step 5 — Report

`dev/reports/REPORT_declcal_c0_rchange_2026-09-22.md`, bullets: pins and commits; the quoted transform for
`c2` and how `c0` follows it; every post-condition with measured values; Test 5's `rho_hat`, the two `delta`s, the
closed-form targets, realized values and discrepancies; and on the GBSG test fit the table
`c0 in {0.70, 0.75, 0.80, 0.85, 1.0}` × {`kappa_hat_05`, `kappa_hat_10`, `pstar_implied_05`, `fw_size at p* 0.90`,
`n_admitted_cal05`}. OPEN ITEMS block for documentation gaps.

---

## 9. Commit plan

1. `docs(tasks): add c0 null-level task document (2026-09-22)`
2. `feat(declaration): add clinically specified null level c0 to the calibrated declaration threshold`
3. `docs(verification): commit declaration-calibration post-condition scripts` — the two prerequisite scripts.
4. `docs(report): record c0 null-level implementation and verification`

No push.

---

## 10. Out of scope

- No default changes; no change to the screen, the correction, the bounds, or the rounding.
- No `sg_focus` re-run on the calibrated set; no contrast field; no `pconsistency.digits` pass-through.
- No campaign. The re-run is `TASK_declcal_c0_campaign_2026-09-22.md`, a separate go.
