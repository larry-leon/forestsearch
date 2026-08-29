# CC task — the consistency screen's parallel dispatch: stop treating `plan = "sequential"` as parallel

**File:** `dev/tasks/cc_task_parallel_dispatch_2026-08-30.md` · **Issued:** 2026-08-30 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Follows:** `REPORT_bootstrap_profile_2026-08-30.md`, which measured the effect this task removes.

**What the profile found.** `subgroup.consistency()` accepts `plan = "sequential"` as a **valid parallel configuration**. The bootstrap forces `parallel_args = list(plan = "sequential", workers = 1L)` on the inner search — a caller saying "do not parallelise," twice — and it passes all three validations, sets up a future backend, and batches 749 candidates into **375 two-candidate `future_lapply` calls**, each re-resolving globals. Measured: **29.4 s against 4.5 s** for the identical search, **85.6%** of it globals resolution around **1.6 s** of real consistency evaluation. At B = 1000 that is **490–537 min against ≈ 75 min**.

The driver's `list(plan = "sequential")` avoids this **only by accident** — it trips the "missing required elements" warning at L667 and falls back to the plain loop.

**Larry's decision, 2026-08-30: the general scope.** Fix the dispatch, not the caller.

---

## ⚠ `R/` CALLOUT — **CHANGES BEHAVIOUR**

Same estimand, same results, different execution path — and **that is exactly what §4 must prove.** Warning behaviour also changes; §6 authorises a specific, bounded decrease.

Editable: `R/subgroup_consistency_main.R` (Section 6's validation block, and its roxygen for `parallel_args`), `man/`, `NAMESPACE` as regenerated, `DESCRIPTION`, `NEWS.md`. **Nothing else.**

**Explicitly NOT changed:** `R/bootstrap_dofuture_main.R` (the caller stays as it is — the point of the general scope is that the caller no longer needs to be right); the batch-size logic at L831–858; `valid_plans` itself; any driver, application document, payload, or OC-wrapper file.

**Compute:** §2's reference capture and §4's comparison are a handful of `forestsearch()` calls on two fixtures — tens of seconds each pre-change, single-digit seconds post. §5 is a 20-replicate bootstrap. Estimate **45–75 minutes** including `R CMD check`. No simulation study, no renders.

**Unattended.** Gates stop the task; never ask, never work around.

---

## 1. Provenance and first commit — GATE

Standard block (`hostname; pwd; git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD; git status --porcelain; git log --oneline -4`) plus the installed version (expect **0.3.0**). *GATE:* branch `feature/glm-extension`, clean tree. Copy this document into `dev/tasks/` and commit alone; the `~/Downloads` stem arrives with hyphens stripped, so report both names.

**Then verify chat's source quotations against the working tree** — chat read a 28 Aug snapshot and the tree has moved since. Quote back, from the current file:

- Section 6's validation block (chat has it at L663–683), including `valid_plans`;
- the batch-size block (chat has it at L831–858), specifically the `min(n_workers * 2L, n_candidates)` line;
- how `bootstrap_dofuture_main.R` sets the inner `parallel_args` (chat has `resolve_bootstrap_parallel_args` at L279).

*GATE:* if the dispatch logic differs materially from the above, **stop and report** — the premise has moved and the edit should not be written blind.

---

## 2. Capture the pre-change reference — before any edit

Five configurations. **All five are captured pre-change and re-run post-change.** State the fixtures exactly.

- **Continuous / MD40**, one trial at `n = 500`, the drivers' `forestsearch_args` — the 749-candidate case the profile used.
- **Survival / gbsg-based**, the gbsg application's arguments, comparable size.

| # | fixture | `parallel_args` | consistency | early stopping | what it tests |
|---|---|---|---|---|---|
| C1 | continuous | `list(plan = "sequential", workers = 1L)` | `resample` | off | the main path, batch size 2 |
| C2 | continuous | `list(plan = "sequential", workers = 1L)` | `resample` | `maxeffCons` + `stop_threshold` | batch size 1, and whether the halt lands on the same candidate |
| C3 | continuous | `list(plan = "sequential", workers = 1L)` | `split`, `fs.splits` reduced (state the value; 50 is enough) | off | **the RNG-consuming path** |
| C4 | survival | `list(plan = "sequential", workers = 1L)` | `resample` | off | the other outcome type |
| C5 | continuous | `list(plan = "multisession", workers = 2)` | `resample` | off | **the genuinely parallel path must be untouched** |

Capture, per configuration, at minimum: `sg.harm`; `sort(which(sg.harm.id == 1))`; `find.grps$out.found$hr.subgroups` as a set of rows; `grp.consistency$n_candidates_total`, `n_candidates_evaluated`, `n_passed`, `early_stop_triggered`, `early_stop_candidate`; the fitted effect and CI on the selected subgroup; and wall-clock. Save to `.rds`.

**Compare subject membership, never the definition string** — the smoke-test primitive, reused.

---

## 3. The edit

In Section 6, **before** the `required_parallel` check, add the sequential short-circuit:

- if `parallel_args$plan` is `"sequential"`, **or** `parallel_args$workers` is present and `<= 1`, set `use_parallel <- FALSE` and skip the remaining validations.
- **No `warning()`.** This is the code doing what the caller asked, not an error condition, and a warning here would fire once per bootstrap replicate. A `details`-gated one-line message is fine.

Placing it first has a deliberate second effect: the driver's `list(plan = "sequential")` now takes this branch instead of tripping `"parallel_args missing required elements. Using sequential."` — **the same destination, reached deliberately, without a misleading warning.** §6 handles the parity consequence.

Update the `@param parallel_args` roxygen to state plainly that a sequential plan or a single worker runs the plain loop, that this is a performance property and not a results property, and that `batch_size` is the knob for genuinely parallel runs.

**Do not touch** `valid_plans`, the batch-size logic, or the bootstrap caller.

---

## 4. The identity gates — the whole point of the task

Re-run all five configurations post-change and compare against §2.

**Hard gates — `identical()` on every captured component:**

- **C1, C2, C4** — `resample` consumes no RNG, so there is no mechanism for a difference. Any failure means the edit did more than change which branch runs.
- **C5** — the parallel path must be **bit-identical**; if it is not, the short-circuit is catching a configuration it should not.
- **C2 additionally**: `early_stop_candidate` and `n_candidates_evaluated` must match exactly. Batch size 1 is supposed to halt at the same point as the plain loop; this is where that is verified rather than assumed.

**C3 — characterised, not gated on identity.** `split` consumes RNG, and `future_lapply` lays the stream out differently. Run it with an explicit `set.seed()` immediately before each call. Then:

- if `identical()` — say so, and the concern is closed;
- if not — report the **magnitude**: the distribution of per-candidate differences in `p.consistency`, and whether the selected subgroup changed. A difference consistent with re-drawing at the stated `fs.splits` (i.e. within Monte-Carlo noise for that many splits) is **benign and expected** — proceed, and record it in NEWS as a documented consequence for `consistency_method = "split"` users;
- if the difference is **larger than that**, or the selected subgroup changes, **stop and report.** That would not be an RNG-layout effect.

*GATE:* C1, C2, C4, C5 identical; C3 either identical or explicable as above.

---

## 5. Measure what was bought

Report, per configuration, wall-clock pre and post, and specifically:

- the C1 A/B against the profile's 29.4 s → expected ≈ 4.5 s;
- **the survival case** — the profile recorded gbsg future overhead at 18% but **not which `parallel_args` that profile ran under.** Pin it here: state C4's pre and post wall-clock, so the survival saving is a measured number rather than an inference.

Then re-run the bootstrap check from the profile task: `forestsearch_bootstrap_dofuture()` at the same small `nb_boots` (20–50, state which) on the continuous fixture. Report per-replicate mean against the profile's 32.2 s, and **project B = 1000** against the profile's 490–537 min.

---

## 6. Close-out

- `devtools::document()`; `devtools::test()`.
- **Warning-count parity will legitimately change, downward.** The baseline is 31. Any decrease must be accounted for **test by test**: name each test whose warning count moved, and confirm the warning that disappeared is `"parallel_args missing required elements. Using sequential."` and nothing else. **An increase, or a decrease from any other warning, is a failure — stop and report.** Do not adjust a test to pass without saying so and why.
- `R CMD check` (`RSTUDIO_PANDOC` on Rscript's PATH). Status verbatim; 0/0/0 is the target.
- `DESCRIPTION` → **0.3.1**. If §4's gates all pass, results are unchanged by construction, so this is a patch.
- `NEWS.md`: the dispatch now treats a sequential plan or a single worker as sequential; the measured effect (29.4 s → 4.5 s on the profiled fixture, and whatever C4 shows for survival); no change to results under `resample`; the C3 outcome for `split` users, whatever it was; `parallel_args$batch_size` as the knob for genuinely parallel runs.
- Commit, **explicit paths only**, never `git add -A`. **No push.** `devtools::install()` (`upgrade = FALSE`); confirm `0.3.1`.
- Report `REPORT_parallel_dispatch_2026-08-30.md`: §1's source verification quoted · the diff · §4's five results with the C3 characterisation in full · §5's timings and projection · test and check output raw · the warning-parity accounting · commits · ten-line verdict.

---

## 7. Out of scope

- No edit to `bootstrap_dofuture_main.R`, to `valid_plans`, to the batch-size logic, to any driver, application document, payload, or OC-wrapper file.
- **No re-run or re-verification of anything already committed.** The `gbsg_redux` drivers and their 370 payloads are untouched.
- No optimisation beyond this dispatch condition. The batch-size question (`n_workers * 2L` at realistic worker counts) is **recorded and deliberately not addressed here.**
- No full bootstrap campaign, no simulation study, no renders, no push.
