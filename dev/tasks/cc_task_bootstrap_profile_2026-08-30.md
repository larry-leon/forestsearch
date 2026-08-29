# CC task — profile `forestsearch()` to locate the bootstrap's cost

**File:** `dev/tasks/cc_task_bootstrap_profile_2026-08-30.md` · **Issued:** 2026-08-30 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** Larry places this file in `~/Downloads`; CC copies it into `dev/tasks/` and commits it first (§1).
**Run order:** independent of the OC-wrapper queue, but not concurrently with it.

**Why.** The full bootstrap is the bottleneck in simulations; MR resampling is comparatively fast. `forestsearch_bootstrap_dofuture()` runs a complete `forestsearch()` per replicate, so the cost is inside the search. Chat's reading of the source (`memo_bootstrap_performance_2026-08-30.md`) produced a ranked list of suspected hot spots. **This task measures rather than assumes.** No optimisation is written here; the profile decides what is worth writing.

---

## ⚠ CATEGORY

**No `R/` change of any kind.** No edits to any package file, driver, or document. Writes: one scratch script, its outputs, and the report, all under `dev/glm-continuous-sims/` (or a `dev/profiling/` directory if you prefer — say which). Plus this task document under `dev/tasks/`.

**Compute:** bounded diagnostic. A handful of `forestsearch()` calls and one short bootstrap. Estimate 20–40 minutes. No simulation study, no full bootstrap campaign, no renders.

---

## 1. Provenance and first commit — GATE

Standard block (`hostname; pwd; git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD; git status --porcelain; git log --oneline -4`) plus the installed version. *GATE:* branch `feature/glm-extension`, clean tree. Copy this document into `dev/tasks/` and commit it alone; the `~/Downloads` stem arrives with hyphens stripped, so report both names.

**Record the installed `vi.grf.min` default** — if the queued default task has already landed it is `NULL`, otherwise `-0.2`. Every timing below must state which.

---

## 2. Two source questions, answered before profiling

1. **Are `med0` / `med1` consumed anywhere?** They enter the candidate table as `m1` / `m0` (`create_result_row`). Search every selection and summary path — `sort_subgroups()`, `sort_subgroups_preview()`, the `hrMaxSG` / `hrMinSG` band, `subgroup_consistency_main.R`, the bootstrap summaries, `fs_mr_oc_summary()`, any printing or plotting — and report every read of those columns. **Verdict: display-only, or consumed by selection.** Quote each site.
2. **Does the closed-form consistency ever fall back to `"split"`?** v5 §4 records the fallback exists for candidates where the closed form is inestimable. Establish from source when it fires, and — if cheaply instrumentable — count how often it fires in one search on each fixture. A candidate hitting the fallback costs `fs.splits` refits, so if it is common it dominates everything else.

---

## 3. Profile

Two fixtures, both small enough to profile but large enough to be representative:

- **Continuous / MD** — the MD40 DGM the OC-wrapper scripts rebuild, one trial at `n = 500`, the drivers' `forestsearch_args`.
- **Survival** — a GBSG-based fixture at comparable size, using the gbsg application's arguments. State exactly what you used.

For each fixture, profile **one `forestsearch()` call** with `Rprof(..., line.profiling = TRUE, memory.profiling = FALSE)` or `profvis`, at a sampling interval fine enough to resolve the candidate loop. Report:

- total wall-clock;
- the top 20 functions by **self** time and by **total** time;
- and the share of total time attributable to each of these, named explicitly:

| bucket | what to attribute to it |
|---|---|
| per-candidate effect fit | `fit_glm_for_subgroup` / `fit_cox_for_subgroup` and everything under them, including `lm`/`glm`/`coxph`, `model.frame`, `model.matrix`, `vcov`, `summary.coxph` |
| candidate medians (survival) | `survfit`, `summary.survfit` and below |
| GRF variable importance | `causal_forest` / `causal_survival_forest` / `fit_causal_forest_glm`, `variable_importance` |
| cut construction | `get_FSdata` and below |
| enumeration / floors | `generate_combination_indices`, `evaluate_combination_with_status`, `meets_prevalence_threshold`, `extract_idx_flagredundancy`, excluding the fit |
| consistency screen | `evaluate_subgroup_consistency`, `consistency_resample` and below |
| dedup / selection | `remove_near_duplicate_subgroups`, `.maxeff_membership_dedup`, `sort_subgroups*` |
| everything else | remainder |

Also report, per fixture: the number of candidates reaching a fit, the **mean and median time per candidate fit**, and the number reaching the consistency screen.

**Then repeat the continuous profile with `vi.grf.min` toggled** (`-0.2` and `NULL`) to measure the forest's share directly — this quantifies what the queued default change is already worth to the bootstrap.

---

## 4. Extrapolate to the bootstrap, and check the extrapolation

From the per-call profile, project the wall-clock of a `nb_boots = 1000` bootstrap and state the projection.

Then **check it**: run `forestsearch_bootstrap_dofuture()` at a small `nb_boots` (20–50 is enough; state which) on the continuous fixture, time it, and compare the per-replicate mean against the profiled single call. Report the ratio. If they disagree by more than ~20%, something outside the search dominates and the profile's ranking does not transfer — say so prominently, and report where the extra time goes.

---

## 5. Report

`REPORT_bootstrap_profile_2026-08-30.md`, alongside the other reports:

1. Provenance, including the `vi.grf.min` default in force.
2. §2's two source verdicts, with quoted sites.
3. The bucket table per fixture, with the `vi.grf.min` toggle comparison.
4. Per-candidate fit counts and times.
5. §4's projection and the measured check.
6. **A ranked list of what is worth changing**, each with the measured share it would recover and a one-line note on risk. Rank by measured share, not by chat's memo — and **say explicitly where the measurement contradicts the memo**, which is the most useful thing this report can contain.
7. Ten-line verdict.

Commit the script, its profile outputs and the report by explicit path. **No push.** No `devtools::install()` — nothing in `R/` changed.

---

## 6. Out of scope

- No optimisation, no `R/` edit, no change to any estimator, driver or document. The profile decides what gets written; nothing gets written here.
- No full bootstrap campaign — §4's check is 20–50 replicates, not a production run.
- No renders, no push.
