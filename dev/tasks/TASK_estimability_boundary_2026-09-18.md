# CC TASK — the estimability boundary, and `fs_dgm_feasibility()`

**Opened:** 2026-09-18 · **Repository:** forestsearch · **Authorized by:** Larry, 2026-09-18.
**Workstream:** admission floors. **Predecessor:** `dev/reports/REPORT_four_cell_exposure_2026-09-18.md`
(`a3374e6f`) — binary invariance proven on 3,313 replicates; survival not checkable by membership.
**Compute:** micro-fits on tiny synthetic data inside the acceptance tests, plus one small
`fs_dgm_feasibility()` run (`n_rep` ≤ 50). **Hard abort 5 minutes.** No campaign, no bootstrap, no CV.
**Testing:** ONLY this task's own acceptance test files, hard abort 3 minutes. No R CMD check, no full suite.
Every step capped. Report ≤ ~120 lines to `dev/reports/REPORT_estimability_boundary_2026-09-18.md`.
No install. Do not push.

## Larry's governing stance — least restrictive

**No new admission floor anywhere.** `n.min`, `d0.min`/`d1.min`, and every DINA/GRF floor stay exactly as
they are. Nothing is forwarded to DINA or GRF. This task adds only (A) honest handling where the estimand
does not exist, and (B) a design-time feasibility check. Neither restricts what an analysis may admit.

## Part A — the estimator boundary (changes behaviour, narrowly)

**Today:** a candidate on which the estimand does not exist is fitted anyway; `glm()` returns a finite
divergent coefficient, `converged` is computed by every estimator and discarded at
`fit_glm_for_subgroup()`, and the candidate is admitted with a nonsense estimate (audit §2.4).

**Required:** at the effect-estimator boundary — the closures in `glm_effect_estimators.R` and the Cox
fitter, i.e. the place every consumer (FS search, DINA/GRF refit, MR family) obtains an estimate — check the
estimand's **existence condition** before fitting, per estimand and nothing more:

| Estimand | Existence condition | Otherwise |
|---|---|---|
| OR | all four cells ≥ 1 (control/treated × events/non-events) | non-estimable |
| RR, IRR, HR | ≥ 1 event in each arm | non-estimable |
| RD, IRD, MD | none — **untouched, byte-identical behaviour** | — |

A non-estimable candidate returns `estimate = NA`, `se = NA`, `converged = FALSE`, and a **`reason`** string
naming the empty cell (e.g. `"non-estimable: treated non-events = 0"`). Additionally, for the ratio
estimands only, a fit that reports `converged = FALSE` is treated as non-estimable with reason
`"non-convergent fit"` — **RD's tier-3 raw-proportions fallback returns `converged = FALSE` by design and
must not be caught by this; RD is out of scope entirely.**

**Visibility (Larry's default, overridable):** the candidate takes the existing fit-failure status and cannot
rank (an NA effect has no place in selection), but it is **counted with its reason** in `filter_counts` and
surfaced in `fs_family_report()`'s stage map, so a run says "k candidates non-estimable (zero cell)" rather
than dropping them silently. No selection logic changes.

**Exposure gate (Gate A), read-only, before any edit:**
1. Binary: cite `a3374e6f` — proven.
2. Continuous: no rule applies to gaussian — invariant by construction; state it.
3. Survival, two instruments, neither needs data regeneration: (a) grep every committed survival campaign
   template for the `d0.min`/`d1.min` it passes — FS declared subgroups then carry ≥ that many events per
   arm by admission; (b) artifact-only scan of every committed survival bundle's **stored** declared
   estimates and SEs (Ĥ and Ĥᶜ) for non-finite values, which is what a zero-event arm would have produced.
   Report the maximum |log HR| and SE seen, for the record. **Any non-finite stored estimate → STOP and
   present pin-vs-proceed.** Bundles whose stored results cannot be read are listed as not checkable.

**Baseline and gates:** before editing, record digests of micro-fits on tiny synthetic data for: OR with all
cells ≥ 1; OR with one zero cell; RR with a zero-event arm; RD with a zero cell; HR with a zero-event arm; HR
normal; MD normal. After: the all-cells-positive OR, RD, MD and normal-HR fits are **byte-identical** to
baseline; the zero-cell OR, zero-event RR and zero-event HR are NA-with-reason; `fs_family_report()` shows
the count and reason. `converged` is read, not discarded — assert from source.

## Part B — `fs_dgm_feasibility()` (adds code; changes no existing behaviour)

The design-time check Larry asked for, as a second step on a DGM object:

```r
feas <- fs_dgm_feasibility(dgm, n = c(500, 750, 1000, 2000),
                           n.min = 60, d0.min = 10, d1.min = 10,
                           n_rep = 200, tolerance = 0.05)
feas            # loud print, unconditional
feas$feasible   # what a template gates on
```

- **Interface:** accept the `dgm` objects the OC family already accepts (`fs_oc_predict()`, `fs_oc_grid()`)
  — establish that interface from source. If the survival DGM builder produces something the GLM path
  does not share, implement for the GLM `dgm` and **report exactly what survival needs**; do not guess.
- **Per n, over `n_rep` draws of the planted region:** share undeclarable (size ≤ `n.min`, strict as the
  search is), share under the events floor (either arm < `d0.min`/`d1.min`), share non-estimable under
  Part A's estimand-specific condition, plus size and per-arm cell summaries (mean, 5th/95th percentile,
  minimum).
- **`feasible`** = every undeclarable share ≤ `tolerance` (default 0.05 — Larry's to change).
- **The hard requirement:** draw replicates through the DGM's **own** generator path — the same function
  the campaign templates call — never a parallel implementation. The admission check recorded that
  calibrating after an RNG-kind switch yields a different super-population; the roxygen states the
  calibration-before-kind-switch requirement and the function must not itself change the RNG kind.
- Exported, roxygen, `NEWS.md`, and a test on a tiny DGM with small `n_rep` (seconds).

## Not in this task

The oracle helper (keep the legacy pooled 5/5, add the four-cell condition) and the template's Stage 0 call
to `fs_dgm_feasibility()` — study-side, next task. The survival-bundle reproducibility finding — recorded
for the campaign owners, not diagnosed here. Any floor value. Anything under `subgroup_method = "dina"` /
`"grf"` beyond what their refits inherit from the boundary.

## Commits and report

Part A (behaviour) and Part B (addition) as separate commits, tests and NEWS with each. The report carries
Gate A's three instruments, the digests, the boundary rule table, the two defaults Larry may override, and
findings with no tasks attached. **Do not push.**
