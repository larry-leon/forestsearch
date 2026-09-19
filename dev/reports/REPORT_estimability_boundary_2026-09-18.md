# REPORT — the estimability boundary (Part A) and `fs_dgm_feasibility()` (Part B)

HEAD at open: `cf97ec70`. Commits: `1719056f` (Part A), `f9b794f6` (Part B).
Branch `feature/glm-extension`. Task: `dev/tasks/TASK_estimability_boundary_2026-09-18.md`.
Predecessor: `dev/reports/REPORT_four_cell_exposure_2026-09-18.md` (`a3374e6f`).
No install, no `R CMD check`, no full suite, not pushed.

---

## 1. Gate A — exposure, read-only, before any edit

**1. Binary — proven, cited.** `a3374e6f`: all 3,313 declared replicates across both committed
`orfs_*` cells, membership size equal to the recorded `n_sel` on every row, **0** violations of the
four-cell condition. The closest approach was a treated-arm non-event count of **4**.

**2. Continuous — invariant by construction.** No existence condition applies to a mean difference:
MD exists on any slice with both arms present. `subgroup.search()` also skips the per-arm events
floor entirely for continuous and count (`R/subgroup_search.R`, Status 3). Nothing in Part A touches
the gaussian path, and the MD micro-fit is byte-identical below.

**3. Survival — two artifact-only instruments, neither regenerating data.**

*(a) The `d0.min`/`d1.min` each committed template passes.* The template of record,
`sim_fs_maxeffCons_fb_mr_field_m1_template.qmd:573-575`, passes `d0_min = 10L; d1_min = 10L`, and FS
declared subgroups under it carry at least 10 events per arm by admission. **But the floors are not
uniform**: 124 template instances — the whole `maxeffCons_mr_coverage_sweep_*` family — pass
`d0_min <- 10L; d1_min <- 0L`. Under those, **the treated arm has no event floor at all** and a
zero-event treated arm reaches `coxph()`. Instrument (a) therefore does not by itself cover the
committed survival record, which is why (b) matters.

*(b) Stored-estimate scan, every committed survival bundle.* 930 bundles under
`quarto/simulations/gbsg{,_020,_redux}`, **568,757 declared replicate rows**, every `*_(H|Hc)_est`,
`_est2`, `_se`, `_se_ij`, `_se_w`, `_se_wf` column plus `betaHhat_H`/`betaHhat_Hc`.

```
non-finite stored estimates or SEs : 0
max |stored log HR|                : 9.69284   (HR ~ 1.6e4)
max  stored SE                     : 6.26196
bundles not readable               : 0
```

**No STOP: zero non-finite values, so no pin-vs-proceed decision is owed.** For the record, 652 rows
carry `|log HR| > 5`, the largest at 9.69, concentrated in `nv_H_est` on the **DINA and GRF**
campaigns at n = 500 — which is exactly where the per-arm floor does not apply, since `d0.min`/`d1.min`
are never forwarded to those collectors. These are finite and so were never caught by anything; they
are reported, not acted on.

A wider first pass over all 994 committed bundles found a maximum stored estimate of 3.0e8. That is
**not** a survival pathology: those rows are the `md*` mean-difference campaigns, where the estimate
is on the natural CD4 scale and a value in the hundreds is ordinary. Recorded so the number is not
later misread.

---

## 2. Part A — the boundary rule

Checked **before** the fit, at the effect-estimator boundary, per estimand and nothing more:

| Estimand | Existence condition | On failure |
|---|---|---|
| OR | all four cells `>= 1` (control/treated x events/non-events) | non-estimable |
| RR, IRR, HR | `>= 1` event in each arm | non-estimable |
| RD, IRD, MD | none | **untouched, byte-identical** |

Non-estimable returns `estimate = NA`, `se = NA`, `converged = FALSE` and a `reason` naming the empty
cell. For the **ratio estimands only**, a fit reporting `converged = FALSE` is also non-estimable
(`"non-convergent fit"`), guarded by `!is.na(estimate)` so the pre-existing error branch returns
unaltered. **RD's tier-3 raw-proportions fallback returns `converged = FALSE` by design and is never
reached by the catch** — asserted in the tests both behaviourally and from source (exactly one
occurrence of the catch, scoped to `c("OR", "RR")`).

**Why an existence check and not a `converged` check.** The baseline digests settle it: on a
zero-cell slice `glm()` returned **19.97 with SE 4809 and `converged = TRUE`**, and `coxph()` on a
zero-event arm returned **-20.88 with SE 23946 and `converged = TRUE`**. Neither raised an error and
neither reported non-convergence, so reading `converged` alone — the audit's §2.4 observation — would
have caught none of these.

### 2.1 Micro-fit digests

| fit | digest before | digest after | verdict |
|---|---|---|---|
| OR, four cells `>= 1` | `a6da7da7b84384b2287ea5ef202a67a1` | `a6da7da7b84384b2287ea5ef202a67a1` | **IDENTICAL** |
| RD, zero cell | `935bca57e6c851310f740359759b0d6e` | `935bca57e6c851310f740359759b0d6e` | **IDENTICAL** |
| RD, normal | `76746b29af4d21440b6d8ff7561d3671` | `76746b29af4d21440b6d8ff7561d3671` | **IDENTICAL** |
| HR, normal | `9c56775fa9d45cf40a34bb9ab9ed7233` | `9c56775fa9d45cf40a34bb9ab9ed7233` | **IDENTICAL** |
| MD, normal | `645e7bd23c9f4267046b4ae826ce5c74` | `645e7bd23c9f4267046b4ae826ce5c74` | **IDENTICAL** |
| OR, zero cell | `bc7321ea07b8f2837e7bc56104d74fd1` | `c6068fea8c56f541dd501b115c15b9ec` | changed, as required |
| RR, zero-event arm | `f765e7fb87d539c9b8c46fbccd599d43` | `afe16f7f0eb3ef5728ebd646fb0da5ec` | changed, as required |
| HR, zero-event arm | `ecac31d4b584582522121bbe1b092817` | `c28347673ac2620659bc1d7f8b381183` | changed, as required |

The three changed fits, before → after:

```
OR zero cell     19.97 (SE 4809),  converged TRUE  ->  NA, "non-estimable: treated non-events = 0"
RR zero-event   -19.69 (SE 6616),  converged TRUE  ->  NA, "non-estimable: treated events = 0"
HR zero-event   -20.88 (SE 23946), converged TRUE  ->  NA, "non-estimable: treated events = 0"
```

Byte-identity is preserved on the estimable path because `reason` is added **only** to non-estimable
returns; a successful fit's list is unchanged in length, order and content.

### 2.2 What did not change

**No admission floor was added anywhere.** `n.min`, `d0.min`/`d1.min` and every DINA/GRF floor are
untouched, nothing is forwarded to DINA or GRF, and no minimum count is imposed. A non-estimable
candidate takes the **existing** fit-failure status (status 5), which an `NA` estimate has always
produced — so no selection logic changes and an `NA` effect simply cannot rank.

### 2.3 Visibility

Counted, never silent: `filter_counts$n_nonestimable` and `filter_counts$nonestimable_reasons`;
printed by `subgroup.search(details = TRUE)` as `Non-estimable (estimand does not exist): k` with a
per-reason breakdown; and a new **"estimability boundary"** stage row in `fs_family_report()` naming
the condition for the run's `effect_measure` and stating that it is not an admission floor.

The survival search path is guarded too, at `evaluate_subgroup_combination()`'s status 5 — this is
what Gate A instrument (a) showed to be necessary, given the 124 templates passing `d1.min = 0`.

---

## 3. Part B — `fs_dgm_feasibility()`

**Interface, established from source.** `fs_oc_family_enumerate()` (`R/fs_oc_family.R:170-182`),
which `fs_oc_predict()` and `fs_oc_grid()` both route through, validates
`inherits(dgm, "glm_dgm")`, `is.data.frame(dgm$df_super)` and the presence of `df_super$flag_harm`.
`fs_dgm_feasibility()` accepts exactly that.

**Survival is refused, and what it needs is reported, not guessed.** `setup_gbsg_dgm()` produces a
different object whose own generator is
`simulate_from_dgm(dgm, n, analysis_time, cens_adjust, seed)` — a different signature carrying two
censoring arguments with no GLM counterpart, and both change the per-arm **event** counts this
function reports. Supporting it needs a generator-dispatch layer plus those two arguments threaded
through. That is stated in the error message and in `?fs_dgm_feasibility`.

**The hard requirement is met.** Replicates are drawn with `simulate_from_glm_dgm()` — the campaign
templates' own function — and there is **no re-implementation and no `RNGkind()` call anywhere in the
file** (both asserted by test). The function seeds with `set.seed()` under whatever generator is
current and restores the caller's `.Random.seed` on exit. The roxygen states the caller's matching
obligation: build and calibrate the DGM *before* any replicate switches the kind.

**Live run, ACTG175 OR 0.75 design, `n_rep = 50`, seed 8316951** (RNG kind verified unchanged
across the call):

```
       n size_mean size_q05 size_min   undeclarable   under_events  non-estimable
     500      48.8     39.5       38          90.0%         42.0%           0.0%
     750      72.2     59.9       56           6.0%          0.0%           0.0%
    1000      96.1     82.0       80           0.0%          0.0%           0.0%
    2000     189.6    165.4      162           0.0%          0.0%           0.0%

  NOT FEASIBLE at n = 500, 750
```

This independently reproduces the admission check's §4.1b, computed there by a different route on 200
replicates: 96.0% vs 90.0% undeclarable at n = 500, 5.0% vs 6.0% at n = 750, 44.5% vs 42.0% under the
events floor at n = 500, and a mean region size of 48.5 vs 48.8.

`$feasible` is `all(share_undeclarable <= tolerance)`. Floors are read and reported, never applied.

---

## 4. The two defaults Larry may override

1. **`tolerance = 0.05`** — the largest undeclarable share that still counts as feasible. Arbitrary,
   and the task names it as Larry's to change. At 0.05 the ACTG175 OR design fails at n = 750 on a
   6.0% share; at 0.10 it would pass there and fail only at n = 500.
2. **Visibility: the candidate takes the existing fit-failure status and cannot rank.** The
   alternative — a distinct status code of its own, separating "non-estimable" from "fit failed" in
   the stage map rather than only in `nonestimable_reasons` — was not taken, per the task's stated
   default. It is a one-status change if wanted.

---

## 5. Findings — no tasks attached

1. **The per-arm floor is not uniform across committed survival templates.** 124 instances pass
   `d1.min = 0`: the treated arm has no event floor, and a zero-event treated arm reached `coxph()`
   in those campaigns. Only the stored-estimate scan, not the floor, establishes that nothing
   diverged to a non-finite value.
2. **`converged` would not have caught this.** All three pathological micro-fits reported
   `converged = TRUE`. Any design that relies on the fitter to report its own failure here is
   unsound; the existence condition is the operative instrument.
3. **The largest committed survival estimates sit where the floor does not apply.** 652 rows with
   `|log HR| > 5`, max 9.69, concentrated in DINA and GRF at n = 500 — the two paths `d0.min`/`d1.min`
   are never forwarded to. Finite, so untouched by Part A, and left alone.
4. **The 3.0e8 maximum over all committed bundles is a scale artefact**, not a pathology: `md*`
   mean-difference campaigns on the natural CD4 scale.
5. **The ACTG175 OR design is infeasible at n = 500 by its own criterion** — 90% of planted regions
   are undeclarable — which is a property of the design, knowable before launch, and now checkable
   in one call.

---

## 6. Tests and scope

`tests/testthat/test-estimability-boundary.R` — **56 pass, 0 fail**.
`tests/testthat/test-fs-dgm-feasibility.R` — **40 pass, 0 fail**.
Run with `devtools::load_all()` + `testthat::test_file()` on these two files only. No `R CMD check`,
no full suite, no install, no push. `devtools::document()` was run (both files carry roxygen);
`NAMESPACE` and `man/` are regenerated, not hand-edited.
