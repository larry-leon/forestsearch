# TASK — Realized declaration rate under a uniform marginal Cox HR of 0.75, GBSG application family

**Repo:** `larry-leon/forestsearch` — **run CC from the forestsearch clone, not fs-glms-interpretable's.**
The fs-glms-interpretable clone is read **read-only**, for the application's settings. Nothing is written there.

**Purpose.** The GBSG application reports FŴ₀.₁₀(0.75) = 0.651: the family-wise probability that the p\* = 0.90
screen declares a harm subgroup if every candidate's true effect were a uniform benefit of HR 0.75. That is a
*boundary* quantity computed from the fit's own multiplier draws. This task measures the **realized** rate by
simulation on the same data structure and the same candidate family — replicate trials under a uniform marginal
Cox HR of 0.75, each run through `forestsearch()` at the application's own {c1, c2, p\*}, counting the fraction
that declare.

**Expected direction.** Below 0.651. FŴ is the boundary size with every candidate sitting exactly at 0.75; a real
design sits inside it, and the executed screen prunes candidates the pre-reduction family retains. The finding is
the size of the gap, not its sign. Do not treat a rate below 0.651 as a failure of anything.

**Estimand.** Uniform benefit is calibrated on the **marginal Cox HR**, the reference scale for every simulation in
this programme. Not the AHR, not the patient-level conditional HR.

**Kind:** simulation. FS only. No MR. No R/ change. Nothing existing is modified.
**Compute (Larry's go given 2026-09-23):** one cell, n = 686, 1,000 replicates, preceded by a 10-replicate timing
run. Expected 15–25 minutes on 14 workers.

---

## 0. First action

1. Copy this file to `dev/tasks/TASK_gbsg_app_null_declaration_2026-09-23.md`; `git add` that path; commit
   `docs(tasks): add GBSG uniform-benefit declaration-rate task (2026-09-23)`.
2. Record HEAD, the installed forestsearch version and build date, R version and platform.
3. `git status --short`: record any pre-existing untracked files; never stage them.

No `git fetch` / `pull` / `push`. Explicit paths on every `git add`.

---

## 1. The application's settings — read from source

Read `~/Documents/GitHub/fs-glms-interpretable/quarto/gbsg/analysis_gbsg_mr.qmd`, **read-only**. From its
`forestsearch()` call and its "settings the package received" table, take and quote in the report:

- the analysis frame and its N (expected 686; assert and report if it differs),
- `confounders.name` and the cut specification (`conf.cont_jcuts`, `cut_type`, `cont.cutoff`),
- `maxk`, `n.min`, `d0.min`, `d1.min`, `fs.splits`,
- `effect.threshold` (c1), `consistency.threshold` (c2), `pconsistency.threshold` (p\*),
- `sg_focus`, `selection_rule`, `consistency_method`, `use_twostage`, `is.RCT`, `outcome_type`, `effect_measure`,
- the outcome, event, treatment and id variable names, and the document's seed.

**Gate 1 (STOP on failure):** the file is found and every setting above is read from it. Reason: the whole point is
to mimic this application's family; a guessed setting measures a different family's multiplicity and the result
would not speak to FŴ.

Everything below uses these values. Do not substitute defaults for any of them.

---

## 2. The null DGM — uniform marginal Cox HR 0.75

Build with `generate_aft_dgm_flex()` (not `setup_gbsg_dgm()`, whose covariate construction is fixed), on the same
GBSG frame the application analyses, carrying the application's confounders as its covariates:

```r
base_args <- list(
  data            = <the application's analysis frame>,
  continuous_vars = <the continuous confounders, from §1>,
  factor_vars     = <the factor confounders, from §1>,
  outcome_var     = <from §1>, event_var = <from §1>, treatment_var = <from §1>,
  model           = "null",          # uniform effect, no planted region: flag_harm is zeroed
  n_super         = 5000L,
  seed            = <the document's seed>,
  verbose         = FALSE
)
k <- calibrate_k_treat(target_hr_overall = 0.75, base_args = base_args,
                       use_ahr = FALSE, tol_rel = 1, verbose = TRUE)
dgm <- do.call(generate_aft_dgm_flex, c(base_args, list(k_treat = k)))
```

**Gate 2 (STOP on failure):** `dgm$hazard_ratios$overall` is within 1% of 0.75, and `flag_harm` is identically
zero on `df_super`. Reason: if the achieved HR or the null structure is not what the report claims, every number
below is mislabelled.

**Finding, not a gate — where the per-candidate null actually sits.** FŴ's null set is per-candidate,
{β(g) ≤ log 0.75 for every g}, and calibrating the *overall* HR does not guarantee it: under an AFT generator the
subgroup log-HRs vary. Enumerate the application's candidate family on `df_super` under the §1 cut specification
and floors, compute each candidate's true effect (`fs_betaHhat_table()`, or the nearest available route — say
which was used), and report: the family size, min / median / max of β_true(g) on the HR scale, and how many
candidates exceed HR 0.75. If the maximum sits above 0.75, say so plainly in the report; it does not stop the run,
but it changes how the comparison to FŴ reads.

---

## 3. The timing run

Ten replicates through the full §4 path, on the worker count you intend to use for the full run.

**Gate 3 (STOP on failure):** the projected wall-clock for 1,000 replicates is under 2 hours. Reason: Larry's go
was given against an estimate of 15–25 minutes on 14 workers; a fivefold overrun means something is wrong with the
configuration, and a long run should not launch unattended on that basis. Report the projection either way.

---

## 4. The run

1,000 replicates. Per replicate `b`:

```r
sim  <- simulate_from_dgm(dgm, n = 686, seed = <seed_base + b>)   # the DGM's own censoring model
res  <- run_simulation_analysis(
  sim_id = b, dgm = dgm, n_sample = 686,
  confounders_base = <the application's confounders>,
  fs_params  = <the application's settings from §1>,
  run_fs = TRUE, run_grf = FALSE, run_fs_grf = FALSE,
  seed_base = <fixed>, verbose = FALSE)
```

Record per replicate: the detection flag (`any.H`), and where it declared, the declared subgroup's label and its N.
Save the per-replicate frame as the payload.

Parallelise as the directory's other drivers do. Report any replicate that errored, with its seed; do not silently
drop it — a dropped replicate changes the denominator.

---

## 5. What to report

`dev/reports/REPORT_gbsg_app_null_declaration_2026-09-23.md`, bullets, one item each:

- Pins: HEAD, forestsearch version and build date, R version, platform, seeds.
- The application's settings as read in §1, quoted.
- The calibration: `k_treat`, achieved overall HR, and Gate 2's outcome.
- The per-candidate null: family size, β_true(g) range on the HR scale, count above HR 0.75.
- Realized event rate and censoring rate in the simulated trials against GBSG's own. A finding: if they are far
  apart the candidates' standard errors will not match the application's, which weakens the comparison.
- **The result: the realized declaration rate over 1,000 replicates, with a Wilson interval**, beside
  FŴ₀.₁₀(0.75) = 0.651 and the gap between them.
- Among declaring replicates: the distribution of declared subgroup size, and the most frequent declared rules.
- Timing: the timing run, the projection, and the actual wall-clock.
- OPEN ITEMS.

Present the rate as a declaration rate against the nominal comparison. Do not frame it as a significance
dichotomy at the null, and do not describe a rate below FŴ as a shortfall or a miss.

---

## 6. Placement and commit plan

Script and payload under `quarto/simulations/gbsg_app_null/`. Add one line to
`quarto/simulations/gbsg_020/current_status.md` §2 pointing at the new directory and the report, so the standing
catalogue records that this sibling run exists — this run is on the application-matched GBSG DGM, not gbsg_020's,
and does not belong in gbsg_020's own numbering.

Commits: 1. task doc; 2. script; 3. payload; 4. report; 5. the catalogue line. Explicit paths; no push.

---

## 7. Out of scope

- No MR, no field capture, no calibrated cutoff, no κ̂ per replicate. The capture fires only when the screen
  declares, and under this null most replicates will not; recording max T would need the `declcal` campaign's
  `trace()` workaround and is deliberately left out of this coarse check.
- No change to any forestsearch `R/` file, to any committed payload, or to anything in fs-glms-interpretable.
- No second cell, no other c₀, no other n, no other HR. One cell.
