# TASK v2 — Calibrated declaration threshold (kappa-hat-alpha) and family-wise size diagnostic (FW-hat-alpha)

**Supersedes** `TASK_declaration_calibration_2026-09-21.md` (v1) in full. Do not execute v1. The only
substantive change is the addition of the closed-form two-candidate acceptance check as **Test 8** in §6,
with its supporting notes in §6.1; everything else is reproduced unchanged so this file is the single
authoritative set.

**Repo:** `larry-leon/forestsearch`
**Branch:** `feature/glm-extension`
**Design written against:** `c4baf79` (2026-09-21); nothing under `R/` changed between 2026-09-19 and that pin
**Governing handoff:** `HANDOFF_declaration_calibration_forestsearch_2026-09-21_v3.md` (v3 addendum + v2 body)
**Kind:** add-only R change, plus its tests, its report, and one zero-compute costing read
**Compute:** NONE. No simulation, no campaign, no resampling pass is launched by this task. Unit tests only.
**Authorization:** Larry's approval of the R change is given on THIS document. The handoff does not authorize
it; this document does.

---

## 0. First action (before anything else)

1. `mkdir -p dev/tasks`
2. Copy this file verbatim to `dev/tasks/TASK_declaration_calibration_2026-09-22_v2.md`
3. `git add dev/tasks/TASK_declaration_calibration_2026-09-22_v2.md` and commit with message
   `docs(tasks): add declaration-calibration task document v2 (2026-09-22)`
4. Record `git rev-parse HEAD` (the pre-change baseline pin) into the report. If HEAD is not `c4baf79`,
   that is not a failure — record the actual SHA and re-verify every quoted line number in Step 1 from
   source rather than from this document.

No `git fetch`, `git pull`, or `git push` at any point in this task. Every `git add` names explicit paths.
Never stage a pre-existing untracked file.

---

## 1. R/ call-out (required classification)

This task touches `R/`. Classified as the protocol requires:

- **Moves existing code:** NONE. No function is relocated, renamed, or re-filed.
- **Changes behaviour:** NONE at any default. Under Branch B (Step 3) exactly one new formal is added,
  with default `FALSE`, to the multiplier-resampling entry point; with that formal unset the executed path
  and the returned object are identical to the pre-change baseline, asserted bit-for-bit in Step 5.
  No existing default value is changed. No existing formal is removed, renamed, or reordered.
- **Changes the method:** NO change to any method the package currently executes. The p-star consistency
  screen remains the one and only admission rule the search applies. What is added is a **post-hoc,
  opt-in, reported** quantity computed from draws the package already produces; it does not re-run the
  search, does not mutate the fit, and does not feed back into admission anywhere in the package.
  (The method itself — the calibrated cutoff and the size diagnostic — was raised and approved in the
  governing handoff; this task implements it strictly as an addition.)

---

## 2. What is being added — self-contained specification

Do not consult the project KB or any external document for these definitions. Everything needed is here.

Notation (manuscript Section 4 names in brackets):

- `beta_hat(g)` — the within-candidate coefficient for candidate subgroup `g`, on the scale the
  consistency screen uses it (log-hazard-ratio on the survival path).
- `db[g,i]` — the dfbeta influence of subject `i` in candidate `g`.
- `sigma_D(g)^2 = sum_{i in g} db[g,i]^2` — the robust variance the screen already uses.
- `c_cons` — the consistency threshold; **this is the package's c2**.
- `c_screen` — the relevance floor; **this is the package's c1**.
- `p_star` — the consistency-rate cutoff.
- `alpha` — the target family-wise declaration rate (default `0.05`).

**(a) Standardized statistic — a relabelling, not a new statistic.**

```
T(g) = { beta_hat(g) - c_cons } / sigma_D(g)
```

Because the production screen evaluates the closed form
`max(0, 2 * pnorm((beta_hat - c_cons) / sigma_D) - 1) >= p_star`, admission already *is*
`T(g) >= qnorm((1 + p_star) / 2)` exactly:

```
2 * Phi(T) - 1 >= p_star   <=>   Phi(T) >= (1 + p_star)/2   <=>   T >= z_{(1+p_star)/2}
```

**(b) Standardized perturbation field.** For multiplier draw `b = 1..B` with shared multipliers `xi[b,i]`:

```
Zstar[b,g] = ( sum_{i in g} xi[b,i] * db[g,i] ) / sigma_D(g)
```

i.e. the correction's own perturbation matrix `D[b,g]`, each **column** divided by its own `sigma_D(g)`.
The standardizing denominator is the robust `sigma_D(g)` the screen uses — never a model-based
standard error.

**(c) Family maximum.**

```
Mstar[b] = max over g in FAMILY of Zstar[b,g]
```

The maximum is one-sided and in the harm direction (larger `beta_hat` = more harm). `FAMILY` is defined
in §3 D2 below.

**(d) Calibrated threshold [kappa-hat-alpha].**

```
kappa_hat(alpha) = empirical (1 - alpha) quantile of { Mstar[1], ..., Mstar[B] }
```

**(e) Family-wise size diagnostic [FW-hat-alpha].** The family-wise size the conventional p-star screen
actually has on the realized family:

```
fw_size = mean( Mstar[b] > qnorm((1 + p_star)/2) )
```

Its threshold is the fit's own `p_star`, not `alpha`; document that explicitly in the roxygen.

**(f) Calibrated admission rule — reported, not applied.**

```
beta_hat(g) >= max( c_screen, c_cons + kappa_hat(alpha) * sigma_D(g) )
```

This is the current rule with `z_{(1+p_star)/2}` replaced by `kappa_hat(alpha)`. The function returns the
set of candidates this rule would admit, alongside the set the current rule did admit. It changes nothing.

---

## 3. Design decisions — settled, not CC's to revisit

**D1 — Shape: a new exported function, plus (only if forced) one opt-in capture formal.**
The deliverable is an exported `fs_declaration_calibration()` consuming a fitted object, following the
existing `fs_*` family (`fs_fdr_report()`, `fs_oc_grid()`, `fs_oc_invert()`, `fs_family_report()`).
Not an argument that alters admission; not a silently-populated field on every fit.

**D2 — §5 near-duplicate collapse: the maximum is taken over the family PRIOR to any fitted-quantity
reduction. This is the ruling; implement it, do not re-open it.**

- `remove_near_duplicate_subgroups()` keys on `(K, n, E, d1, m1, m0, HR, L, U)` rounded to 0.001 — i.e. on
  sample-fitted summaries — on every `sg_focus` rule except `maxeff`. A family reduced on estimated
  quantities is outcome-dependent, which is precisely what the declaration theorem cannot tolerate: the
  theorem needs a covariate-measurable family, and covariate-measurability is what separates the
  forest-search guarantee from the conditional-on-family scope GRF and DINA receive.
- Pre-reduction is also the conservative choice: a larger family raises the maximum and hence
  `kappa_hat(alpha)`.
- Therefore: `family = "prereduction"` is the **default** and the only value used for the calibrated rule.
- `family = "reduced"` is implemented as a **diagnostic only** so the gap is measurable and reportable in
  manuscript §4.6. Its return value is labelled conditional-on-realized-family and must carry that label
  in the printed output. It is never used for admission.
- Report both family sizes: `n_family_prereduction` and `n_family_reduced`.
- §5 items 2 (the `maxeff` asymmetry) and 3 (an argument to select or disable the reduction) are
  **out of scope here** — package-level register items, not this task, and not to be bundled in.

**D3 — alpha and reporting.** `alpha = 0.05` default; the report states `kappa_hat` at both `0.05` and
`0.10`. `fw_size` is reported at the fit's own `p_star`.

**D4 — Ride the existing draws; no second resampling pass.** The limit depends on the multiplier law only
through its first two moments, so the MR draws (centred Poisson, B = 2000 by default) serve, as do the
consistency engine's Rademacher draws. The multipliers must be the **shared** ones — the same `xi[b,]`
across all candidates in a draw. Record which law and which `B` were used in the returned object and in
the report. The one exception is Test 8 of §6, which deliberately runs a Gaussian multiplier law; see §6.1.

**D5 — No mutation, no re-fit, no feedback.** `fs_declaration_calibration()` is pure: it reads a fitted
object and returns a new object of class `fs_declaration_calibration`. It does not modify its input, does
not re-run the search, and does not call the consistency engine.

**D6 — Generic across outcome paths.** The quantities are generic in `beta_hat` and `db`. Implement against
whatever element names the fit exposes, dispatching on fit class if the survival and GLM paths differ.
Record in the report whether they differ. Survival (GBSG) is the tested path; the OLS path is exercised by
Test 8.

---

## 4. Step 1 — Source verification and branch selection (read-only)

Verify from source, at the current HEAD, with quoted lines recorded in the report. Do not rely on the line
numbers below; they are from a prior read and are pointers only.

1. Confirm the production consistency screen passes `method = "closed"`
   (expected near `R/subgroup_consistency_helpers.R:1465`). Quote the line.
2. Confirm which formal carries `c_cons` (c2) and which carries `c_screen` (c1) into that screen, and on
   what scale. Quote the signature.
3. Confirm the MR replay structure: that `d0.min`/`d1.min` and `max_subgroups_search` are not replayed
   under resampling, that the size minimum is applied once at the observed data, and that per-draw
   re-evaluation moves `beta_star` with `t_g` fixed at the observed `sigma_D`
   (expected near `R/forestsearch_main.R:3707-3712`). Quote the lines.
4. Locate the call site of `remove_near_duplicate_subgroups()` and record **where it sits** relative to
   (i) candidate enumeration, (ii) the consistency screen, (iii) the MR replay. Quote the call site.
5. Determine what the fitted object already retains:
   - (i) the shared multiplier matrix `xi` (or a seed sufficient to regenerate it exactly),
   - (ii) per-candidate `db[g,i]`, or the assembled `D[b,g]`,
   - (iii) per-candidate `beta_hat(g)` and `sigma_D(g)`,
   - (iv) the **pre-reduction** family membership (or membership counts sufficient to index it).
6. Record whether the multiplier law is settable at the field-capture call site (Test 8 needs Gaussian
   multipliers). If it is not settable, Test 8 runs in its hand-built-field form only; see §6.1.

**Branch on what item 5 returns — this is a source-determined choice, not a question to bring back:**

- **Branch A (preferred; take it if achievable).** Everything in item 5 is already retained on the fit.
  Implement as **one new file under `R/`, with zero edits to any existing file**. This is the ideal
  outcome: add-only in the strictest sense.
- **Branch B.** Something in item 5 is not retained. Add **one** opt-in capture formal,
  `keep_declaration_field = FALSE`, at the multiplier-resampling entry point, at the minimum number of
  insertion points, storing only what is missing into a new list element `declaration_field`:
  `list(Mstar = numeric(B), Zstar = NULL, beta_hat = , sigma_D = , family_id = , meta = )`.
  Store `Mstar` (length `B`) always when on; store the full `Zstar` matrix (`B x G`) only under a second
  formal `keep_field_matrix = FALSE`, since `G` can be large. When the formal is unset, nothing is
  computed and nothing is stored.

**Record which branch was taken, and why, with the quoted evidence.**

**STOP (failure, report and halt) if:** the pre-reduction family is not reachable at the capture point
without changing executed behaviour, or the shared-multiplier draws cannot be recovered or regenerated
exactly. Do not improvise a substitute (no per-candidate re-draw, no model-based standard error, no
post-reduction family standing in for the pre-reduction one). Write the report and stop.

---

## 5. Step 2 — Implementation

New file: `R/fs_declaration_calibration.R`.

**Exported:**

```
fs_declaration_calibration(fit,
                           alpha = 0.05,
                           family = c("prereduction", "reduced"),
                           ...)
```

Returns an object of class `fs_declaration_calibration`, a list with at least:

- `kappa_hat` — the calibrated threshold at `alpha`
- `fw_size` — the family-wise size of the p-star screen on this family
- `alpha`, `p_star`, `z_pstar`, `c_cons`, `c_screen`
- `B`, `multiplier_law`
- `family_source` — `"prereduction"` or `"reduced"`; the latter carries the
  conditional-on-realized-family label
- `n_family_prereduction`, `n_family_reduced`
- `admitted_current` — candidates the executed p-star rule admitted
- `admitted_calibrated` — candidates the calibrated rule of §2(f) would admit
- `Mstar` — the `B` family maxima (kept; it is cheap and the report needs its quantiles)
- `column_sd` — the per-candidate standard deviation of `Zstar[,g]` (the Step 5 scaling assertion needs it)
- `field_cor` — the empirical correlation matrix of `Zstar` when the family has at most 8 candidates,
  `NULL` otherwise (Test 8 needs it; it is not worth storing for a large family)

**Also:** `print.fs_declaration_calibration()` — a compact summary, one quantity per line, that states the
family source and, for `"reduced"`, prints the conditional-on-realized-family caveat.

**Internal (not exported):** a single helper that assembles and column-standardizes the field, so the
standardization lives in exactly one place. It must be callable on a supplied `db` matrix and multiplier
matrix directly, so Test 8 can exercise it without a fit.

Conventions, non-negotiable:

- tidyverse style; roxygen2 with markdown ON — write literal `%`, `<`, `>`, `&` and never Rd-escape;
  `@section` titles stay plain with no markup.
- Vectorized assembly; no growing objects in loops.
- **No new hard dependency, and no new Suggests.** Base R and the package's existing imports only. Test 8's
  bivariate normal integral is done in base R; see §6.1.
- `devtools::document()` to regenerate `NAMESPACE` and the `.Rd` files. Do not hand-edit either.

---

## 6. Step 3 — Tests

New file: `tests/testthat/test-declaration-calibration.R`. Per the standing check policy, run **only this
file** plus the Step 5 post-conditions. Do **not** run `R CMD check`, `rcmdcheck`, vignette builds, or the
full testthat suite. Hard time cap: **10 minutes**; abort and report on overrun.

1. **Relabelling exactness.** With `kappa_hat` replaced by `qnorm((1 + p_star)/2)`, the calibrated
   admission rule of §2(f) reproduces the executed screen's admitted set exactly on a fixed-seed fit.
2. **Single-candidate calibration.** On a synthetic field with `G = 1` and `B >= 2000`,
   `kappa_hat(alpha)` matches `qnorm(1 - alpha)` within Monte-Carlo tolerance.
3. **Independent-column size.** On a synthetic field with `G` independent standard-normal columns,
   `fw_size` matches `1 - (1 - alpha1)^G` within tolerance, where
   `alpha1 = 1 - pnorm(qnorm((1 + p_star)/2))`.
4. **Monotonicity.** `fw_size >= alpha1` for `G >= 1`, and `kappa_hat` is non-decreasing in family size on
   a nested synthetic family.
5. **Conservatism of the pre-reduction family.** On a fit where the reduction removed at least one
   candidate, `kappa_hat(family = "prereduction") >= kappa_hat(family = "reduced")`.
6. **Purity.** The fit object is unchanged by the call — compare a digest of the input before and after.
7. **Default-off.** Under Branch B, a fit taken with the new formal unset is identical to the Branch-B
   baseline; `fs_declaration_calibration()` on it errors with an informative message naming the formal to
   set, rather than silently computing from partial inputs.
8. **Closed-form two-candidate acceptance check (the sharpest test in the suite — see §6.1).** On one
   two-candidate OLS dataset from supplement S1.7 (configuration A or B), with Gaussian multipliers and
   `B = 200000`:
   - `rho_hat = sum_i db[1,i] * db[2,i] / ( sigma_D(1) * sigma_D(2) )`, computed analytically from the
     dfbeta matrix.
   - **8a.** `field_cor[1,2]` equals `rho_hat` within `4 / sqrt(B)` (= 0.00894).
   - **8b.** `fw_size` equals `1 - Phi2(z_0.95, z_0.95; rho_hat)` within `0.0027`
     (4 binomial Monte-Carlo standard deviations at `B = 200000`).
   - **8c.** `kappa_hat(0.05)` equals the root of `Phi2(k, k; rho_hat) = 0.95` within `0.025`
     (4 quantile Monte-Carlo standard deviations at `B = 200000`).
   - **8d (free bracket, holds for any `rho_hat` in [0, 1], no integral needed).**
     `fw_size` lies in `[0.05, 0.0975]` and `kappa_hat(0.05)` lies in `[1.644854, 1.954508]`, with
     `fw_size` decreasing and `kappa_hat` decreasing in `rho_hat`.

---

### 6.1 Notes on Test 8 — read before implementing it

**Why it is the sharpest test available.** With two candidates the family is two columns, so `B` can be
enormous at negligible cost. That collapses the Monte-Carlo tolerance from useless to tight: at
`B = 2000` the 4-sd band on `fw_size` is ±0.027 against a total range of 0.05–0.0975 — no test at all —
while at `B = 200000` it is ±0.0027. Use `B = 200000`. This test is promoted from optional to **required**.

**The Gaussian-multiplier point, which the check depends on.** Given the data, `Zstar[b,]` has mean zero
and covariance exactly `R_hat` for **any** multiplier law with mean 0 and variance 1 — that part is exact
and needs no asymptotics. But `Zstar[b,]` is exactly **normal** only when the multipliers are Gaussian.
Under Rademacher or centred Poisson it is a weighted sum of iid non-Gaussian variables, so its normality
is a CLT approximation in `n`, and the discrepancy from the closed form then has two sources — Monte-Carlo
error in `B` **and** finite-`n` non-normality — not one. So:

- Run **8a–8c with `xi ~ N(0, 1)`**. There the law is exactly `N(0, R_hat)` given the data, the only error
  is Monte-Carlo in `B`, and the tolerances above are correct as stated. This is the equality gate.
- Then run the same dataset **once more under the production law** (whichever §4 item 6 found) at
  `B = 200000`, and **record** the three discrepancies in the report without gating on them. That second
  run measures the CLT approximation at this `n` and is a finding, not a pass/fail.
- If §4 item 6 found the multiplier law is not settable at the capture site, run 8a–8c through the
  internal helper with a supplied Gaussian multiplier matrix (which §5 requires the helper to accept),
  and record that the production-path variant could not be run. Do **not** add a formal to make the law
  settable — that would be a behaviour change outside this task's classification.

**Computing `Phi2` with no new dependency.** Do not add `mvtnorm`, to Imports or to Suggests. Use the
one-dimensional conditional decomposition in base R: for standard bivariate normal with correlation `rho`,

```
Phi2(a, a; rho) = integrate(function(z) dnorm(z) * pnorm((a - rho * z) / sqrt(1 - rho^2)),
                            lower = -Inf, upper = a)$value
```

Note the **upper limit is `a`, not `Inf`** — it is `P(X <= a, Y <= a)` decomposed over `X = z`. Handle the
two degenerate cases explicitly: `rho = 1` gives `pnorm(a)`, `rho = 0` gives `pnorm(a)^2`. Get the root in
8c with `stats::uniroot` on `[1.0, 3.0]`, `tol = 1e-10`.

**Reference values for the bracket and for eyeballing the report** (computed independently of the package;
`z_0.95 = 1.644854`):

| rho_hat | `1 - Phi2(z.95, z.95; rho_hat)` | root of `Phi2(k, k; rho_hat) = 0.95` |
|---|---|---|
| 0.00 | 0.097500 | 1.954508 |
| 0.10 | 0.096287 | 1.950821 |
| 0.30 | 0.092865 | 1.938467 |
| 0.50 | 0.087811 | 1.916332 |
| 0.70 | 0.080401 | 1.877299 |
| 0.90 | 0.068132 | 1.797586 |
| 0.99 | 0.055811 | 1.698675 |
| 1.00 | 0.050000 | 1.644854 |

**Pick the configuration with overlapping candidates.** If configuration A or B has disjoint candidates,
`rho_hat = 0` and Test 8 degenerates into Test 3 at `G = 2` — it then checks nothing new. Report the
`rho_hat` the chosen configuration gives; if it is below 0.05, say so and note that the test ran in its
degenerate form. Do not construct a different dataset to get a better `rho_hat` — record it and move on.

**Prefer the production assembly path.** Run 8a–8c by calling the exported
`fs_declaration_calibration()` on the real two-candidate OLS fit, so the test exercises the actual
assembly. Additionally run the same three assertions through the internal helper on a hand-supplied `db`
matrix, as a unit test of the standardization in isolation. Both are cheap; run both.

---

## 7. Step 4 — Machine-checkable post-conditions

All are assertions in a script CC runs; each failure is a STOP. These are the protection — there is no
stop-to-ask gate anywhere in this task.

1. **Baseline identity (the "every default untouched" guarantee).** Before editing, at the baseline pin,
   compute and save a digest of the returned object of a small fixed-seed fit on each path that Step 1
   found relevant. After the change, recompute with all new formals unset and assert the digests are
   **identical**.
2. **Formals unchanged.** Capture `formals()` of every function in a touched file before the change and
   after; assert the only difference is the addition of the Branch-B formals with default `FALSE`, and that
   no existing default value, name, or position changed.
3. **Diff scope.** `git diff --stat` against the baseline pin lists only: the new `R/` file, the new test
   file, the generated `.Rd` files and `NAMESPACE`, the report, and (Branch B only) the named insertion-point
   file(s). Assert no other path appears.
4. **Scaling assertion (catches a multiplier-convention mismatch).** `mean(Zstar)` is within Monte-Carlo
   tolerance of `0`, and every `column_sd` is within tolerance of `1` (Rademacher and centred Poisson both
   have unit variance, so unit column scale is the correct expectation). If any column's standard deviation
   departs materially from 1, that is a scaling-convention mismatch, not a tolerance miss: **STOP**, record
   the measured values, and do not rescale empirically to make it pass.
5. **Shared multipliers.** Assert the field was assembled from one multiplier vector per draw shared across
   all candidates, and that `ncol(Zstar)` equals the family size for the selected `family`.
6. **sigma_D provenance.** Assert the `sigma_D(g)` used in standardization equals, to within numerical
   tolerance, the value the closed-form screen itself used for the same candidate.
7. **Family size ordering.** `n_family_prereduction >= n_family_reduced`.
8. **Test 8 closed-form agreement.** Assertions 8a–8d of §6 pass at their stated tolerances under Gaussian
   multipliers. A failure here is a STOP: the standardization or the family maximum is wrong, and no other
   result in this task can be trusted.

---

## 8. Step 5 — Report

Write `REPORT_declaration_calibration_2026-09-22.md` **in the directory where the repo's existing
`REPORT_*` files live** — determine it with `git ls-files 'REPORT_*' '*/REPORT_*'`, do not assume a path,
and do not put it in `dev/tasks/`.

Contents, in bullet form (one item per bullet, short):

- The actual HEAD at start, and the final commit SHAs.
- Which branch (A or B) was taken and the quoted source evidence for it.
- The Step 1 findings, with quoted lines: the `method = "closed"` line, the c1/c2 formals and their scale,
  the MR replay lines, the `remove_near_duplicate_subgroups()` call site with its position in the pipeline,
  and whether the multiplier law is settable at the field-capture call site.
- Whether the survival and GLM paths expose the same element names.
- Every post-condition, with its measured value, not just pass/fail — in particular the `column_sd` range
  and `mean(Zstar)`.
- **Test 8, as its own short block:** the configuration used, `rho_hat`, the three closed-form targets, the
  three realized values under Gaussian multipliers with their discrepancies, and the same three
  discrepancies under the production law recorded as a finding (the CLT-approximation measurement).
- `kappa_hat` at `alpha = 0.05` and `0.10`, and `fw_size`, on the fixed-seed test fit, for both
  `family = "prereduction"` and `family = "reduced"`, with the two family sizes.
- An **OPEN ITEMS** block for documentation gaps (missing line numbers, absent citations). Documentation
  gaps do not stop the run; results that cannot be trusted do.

---

## 9. Step 6 — Zero-compute costing read (reporting only, launches nothing)

So the downstream null-cell evaluation can be priced without another round trip:

- Read `quarto/simulations/gbsg_020/current_status.md` and the records of the two structural-null campaigns
  `nullid` (c1 0.90 / c2 0.80) and `nullc125` (c1 1.25 / c2 1.00).
- Record, in an appendix to the same report: the cell count of each campaign, the replicate count per cell,
  the **measured** wall-clock from their own records (measured only — no estimate), the worker count used,
  and whether the per-replicate `max_g T_g` record that `nullid` retains is sufficient to recompute
  `Mstar` without re-running the identification, or whether the multiplier draws must be taken afresh.
- Nothing is launched. This step reads files and writes bullets.

---

## 10. Commit plan

Three commits, each with explicitly named paths:

1. `docs(tasks): add declaration-calibration task document v2 (2026-09-22)` — Step 0.
2. `feat(declaration): add opt-in calibrated declaration threshold and family-wise size diagnostic` —
   the new `R/` file, the test file, regenerated `.Rd`/`NAMESPACE`, and (Branch B only) the insertion-point
   file(s).
3. `docs(report): record declaration-calibration implementation and verification` — the report with its
   costing appendix.

No push. Larry pushes via GitHub Desktop.

---

## 11. Explicitly out of scope

- No default changes anywhere. The p-star screen stays the default admission rule.
- No change to the post-selection correction, the bounds, the field constructions, or any threshold default.
- Not to be bundled with the c1/c2/p-star threshold-specification item, nor with §5 items 2 and 3.
- No re-run of the identifiers campaign; its classification and operating-characteristic results stand as run.
- No simulation, no campaign, no compute. The applications re-run and the null-cell evaluation are separate
  tasks with their own go/no-go.
- No new formal to make the multiplier law settable, even if Test 8 would be tidier with one.
- No `R CMD check`, no vignette build, no full test suite.
- Do not compete with any live campaign found running on the machine; if one is running, this task's unit
  tests still proceed (they are seconds of single-core work) but record that the campaign was observed and
  untouched.
