# TASK — p\* calibration grid under uniform marginal Cox HR 0.75, GBSG application design

**Repo:** `larry-leon/forestsearch`, branch `feature/glm-extension` — **run CC from the forestsearch clone.**
Nothing is written to `fs-glms-interpretable`. Nothing there is read: the application's settings come from the
committed driver of the 2026-09-23 run, not from the `.qmd`.

**Purpose.** Find the `pconsistency.threshold` (p\*) at which the realized declaration rate under a uniform
marginal Cox HR of 0.75 falls to 10%, on the GBSG application's own design. The 2026-09-23 run measured
0.329 at p\* = 0.90 on a super-population baseline with 1:1 allocation. This task repeats the measurement on
the application's **own covariates and allocation**, and sweeps p\*.

**Kind:** simulation. FS only. No MR. **No `R/` change — nothing in `R/` is read for modification, moved, or
edited.** The driver is a new file; the 2026-09-23 driver is not modified.

**Compute (Larry's go required before Step 5).** Two cells, 5,000 replicates each, preceded by validation and
timing gates of about 400 replicates total. Scaled from the measured 327 s / 1,000 replicates / 14 workers of
the 2026-09-23 run, each cell projects to roughly 8-12 minutes at 48 workers. That projection is an
extrapolation across worker counts and across a changed p\*, not a measured wall clock; Gate C replaces it with
a measured one before the full runs launch.

---

## 0. First action

1. Copy this file to `dev/tasks/TASK_gbsg_pstar_grid_2026-09-23.md`; `git add` that path; commit
   `docs(tasks): add GBSG p* calibration grid task (2026-09-23)`.
2. Record HEAD, installed forestsearch version and build date, R version, platform, worker count.
3. `git status --short`: record pre-existing untracked files; never stage them.

No `git fetch` / `pull` / `push`. Explicit paths on every `git add`.

---

## 1. Settings — inherit, do not re-derive

The application's `forestsearch()` settings were read from source and committed on 2026-09-23. Do **not** re-read
`analysis_gbsg_mr.qmd`.

- Take the settings block verbatim from the committed driver `quarto/simulations/gbsg_app_null/run_gbsg_app_null.R`
  (commit `e6477a22`).
- Carry forward unchanged: `c1 = c2 = 1.00`, `sg_focus = "effMaxSG"`, `selection_rule = "neighborhood"`,
  `effect_neighborhood = 0.20`, `consistency_method = "resample"`, `use_twostage = TRUE`,
  `outcome_type = "survival"`, `effect_measure = NULL`, `maxk = 2`, `n.min = 60`, `d0.min = d1.min = 10`,
  `fs.splits = 1000`, `conf.cont_jcuts = list(er = 10, pgr = 10)`, `cont.cutoff = 4`,
  `conf_force = c("er <= 0", "pgr <= 0")`, `collapse_cuts = TRUE`, `m1.threshold = Inf`, `minp = 0.025`,
  `is.RCT = TRUE`, `est.scale = "hr"`, `use_lasso = use_grf = use_dina = FALSE`.
- **`stop_threshold = NULL` is load-bearing in this task.** At a floor p\* (Step 4) a non-NULL `stop_threshold`
  would halt the scan at the first candidate that clears the floor and destroy the maximum. Assert it is NULL in
  `args_call_all` on every gate replicate.

Only p\* and the simulation call change.

---

## 2. The simulation call — direct, not via the harness

`run_simulation_analysis()` passes only `dgm`, `n`, `seed`, `analysis_time`, `cens_adjust` to
`simulate_from_dgm()` and has no `...`, so it cannot carry the allocation arguments. The driver calls
`simulate_from_dgm()` and `forestsearch()` directly.

Per replicate `b = 1…B`:

```r
df_b <- simulate_from_dgm(
  dgm            = dgm,
  n              = NULL,            # baseline = "fixed" requires NULL or nrow(df_source)
  baseline       = "fixed",         # dgm$df_source: the observed 686, each once
  draw_treatment = <cell>,          # see below
  rand_ratio     = <cell>,
  analysis_time  = Inf,             # NOT the formal default of 48
  cens_adjust    = 0,
  seed           = 8316951L + b
)
```

`analysis_time = Inf` is explicit because `simulate_from_dgm()`'s own default is 48, whereas the 2026-09-23 run
inherited `run_simulation_analysis()`'s `Inf`. Setting it wrong silently truncates follow-up.

**Cell 1 (primary):** `draw_treatment = TRUE`, `rand_ratio = 246/440`. GBSG's covariates fixed, arm
re-randomized each replicate at P(treat) = 0.359.

**Cell 2 (sensitivity):** `draw_treatment = FALSE`. GBSG's covariates and observed arm assignment both fixed;
only outcomes re-drawn.

**Error capture.** The harness's `FS failed:` warning capture is not available on this route. The driver wraps
each fit in `tryCatch` + `withCallingHandlers` and records failures explicitly. A failed fit must never become a
silent non-declaration. The denominator is the full replicate count or the run is void.

**Parallel.** `.options.future = list(seed = TRUE)`, as the programme template does.

---

## 3. Calibration on the fixed baseline (GATE A — stop on failure)

`k_treat = 1.048469` was calibrated to give marginal Cox HR 0.750 on `df_super`. The simulated population is now
`df_source` (686 rows), so the achieved value must be re-verified.

- Build the DGM as the 2026-09-23 driver does (`generate_aft_dgm_flex(model = "null", n_super = 5000,
  seed = 8316951)`).
- Compute the achieved overall marginal Cox HR **on `df_source`** by the same potential-outcome construction the
  prior run used for truth: both potential outcomes per subject under a common extreme-value error, 20 error
  draws per subject, stacked.
- **Gate A passes if the achieved marginal Cox HR is within 1% of 0.750.** If it is not, recalibrate `k_treat`
  against `df_source` (determine the mechanism from `calibrate_k_treat()`'s source) and re-verify. Record both
  the original and any recalibrated `k_treat`.
- Assert `flag_harm` is identically 0.
- Record the patient-level conditional HR for the record. It is not a target.

---

## 4. Candidate family and per-candidate truth on the application's own covariates

With `baseline = "fixed"` the covariates are GBSG's exact 686 rows in every replicate, so cut points and the
candidate family are the application's own, not a population proxy. This replaces the 1,504-candidate
super-population reference family of the 2026-09-23 run and closes that run's per-candidate-null finding on the
real family.

- Enumerate the family on `df_source` by the same route the prior driver used (`get_FSdata()` cuts, `dummy()`
  both directions, combinations up to `maxk = 2`, empty / `minp` / `rmin` / `n.min` floors, identical
  memberships collapsed).
- Record the size and compare to the application's `family_size_prereduction` = 1,744. **Diagnostic, not a hard
  gate:** record the difference and its cause if the counts differ by up to 10%; stop and report if they differ
  by more than 10%, since that would mean the fixed baseline is not reproducing the application's covariates.
- Compute `beta_true(g)` for every candidate, on the marginal Cox HR scale, in **two columns**:
  - **uncensored** (as the prior run did), and
  - **under the DGM's own censoring model**, which is the estimand the screen's partial likelihood actually
    targets and which attenuates less far from 1.
- Report for both columns: min, median, max, and the count and share above HR 0.75.

---

## 5. Validation and timing (GATE B and GATE C — stop on failure)

Run on Cell 1 only.

**Gate B — the floor-p\* equivalence.** Under `consistency_method = "resample"` stages 1 and 2 are bypassed and
the only gate is `Pcons < p*`, so a run at a floor p\* should retain every c1-surviving candidate's `Pcons` in
`grp.consistency$out_sg`, and `I(max Pcons >= 0.90)` should reproduce the p\* = 0.90 declaration indicator
exactly.

- Run 200 replicates at p\* = 0.90 and 200 at **p\* = 0.50** (floor), same seeds.
- Record per replicate at the floor: `max(out_sg$Pcons)` over all rows, or `NA` when `out_sg` is NULL.
- **Gate B passes only on exact agreement across all 200 replicates** between `I(max Pcons >= 0.90)` at the floor
  and the declaration indicator at p\* = 0.90. Not approximate agreement. Report any disagreeing replicate index.
- Also assert `n_candidates_total` matches between the two runs per replicate — the c1 screen is p\*-independent
  and a mismatch would mean the floor run changed the search, not just the threshold.
- **If Gate B fails,** stop and report. Do not fall back to a brute-force grid without a decision; the failure
  itself is the finding.

**Gate C — timing.** 10 replicates at the floor, on the production worker count. Record wall clock and project
5,000. **Hard abort at 2 h for any full cell.** Report the projection and stop for the compute go.

---

## 6. The runs (after Larry's compute go)

- **Cell 1:** 5,000 replicates at p\* = 0.50 (floor), `draw_treatment = TRUE`, `rand_ratio = 246/440`.
- **Cell 2:** 5,000 replicates at p\* = 0.50 (floor), `draw_treatment = FALSE`.
- Per replicate record: `max Pcons`, `n_candidates_total`, treated count, event count and rate, error flag.
- Save the per-replicate frame as `.rds` and `.csv` beside the payload.

---

## 7. Read-out

For each cell, from the per-replicate `max Pcons`:

- The **full empirical CDF** of declaration rate against p\* over [0.50, 1.00], since the rate at any p\* is
  `mean(max Pcons >= p*)`.
- A table at the grid **{0.90, 0.925, 0.95, 0.96, 0.97, 0.98, 0.99}**: declaration rate, Wilson interval at
  alpha = 0.10, and the count.
- **The p\* at which the rate first falls to or below 0.10**, with the Wilson interval at that point and the
  bracketing grid points. State it to the resolution the replicate count supports; do not interpolate beyond it.
- The rate at p\* = 0.90 beside FŴ₀.₁₀(0.75) = 0.651 and beside the 2026-09-23 super-population figure of 0.329.
  This is the first matched-family comparison: same covariates, same allocation, same n as the family FŴ was
  computed on.
- Mean treated fraction (Cell 1) and the assertion that Cell 2's treated count is exactly 246 in every replicate.
- Event rate against GBSG's 0.436.

Present results as one table of cells × grid points followed by a short plain-language reading. No dense prose
summary.

---

## 8. Catalogue and commits

Follow `quarto/simulations/gbsg_app_null/`'s own rule as established on 2026-09-23: the new line goes into
`status_curated.md`, and the regenerated `current_status.md` is committed alone as a child commit. State the
commit the catalogue describes and assert it equals HEAD at commit time.

Commits, explicit paths, in order: task doc; driver; payloads and logs; report; `status_curated.md`;
regenerated `current_status.md`.

---

## POST-CONDITIONS (machine-checkable)

1. Gate A: achieved marginal Cox HR on `df_source` within 1% of 0.750; `flag_harm` identically 0.
2. `args_call_all` on a gate replicate shows `stop_threshold = NULL` and every setting of §1 as received.
3. Gate B: exact agreement on all 200 replicates; `n_candidates_total` matches per replicate.
4. Gate C: measured wall clock recorded; no cell exceeds the 2 h abort.
5. Errors: 0 of 5,000 per cell, or every failure enumerated with its replicate index. Denominator stated
   explicitly and equal to the replicate count.
6. Cell 2: treated count exactly 246 in every replicate. Cell 1: mean treated fraction within 0.01 of 0.359.
7. One replicate re-run sequentially under `RNGkind("L'Ecuyer-CMRG")` reproduces its parallel result exactly.
8. Event rate within 0.02 of 0.436 in both cells.
9. Catalogue pin equals HEAD at commit time.
10. Nothing written outside `quarto/simulations/gbsg_app_null/` and `dev/tasks/`. No `R/` file modified.

---

## OUT OF SCOPE

No MR, no field capture. No κ̂ per replicate. No third cell, no other c0, n or HR. No `R/` change — including the
`modifyList()` NULL-drop in the legacy `run_fs` paths, which stays an open item. No characterisation of declared
subgroups at the selected p\* (the floor run's declared subgroup is not the p\*-specific one; a confirmatory run
at the selected p\* is a separate decision). Nothing written to `fs-glms-interpretable`. Manuscript placement
belongs to the drafting chat.
