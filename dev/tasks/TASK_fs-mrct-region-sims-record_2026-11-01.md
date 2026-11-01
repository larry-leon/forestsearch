# TASK: forestsearch — the branch, and the per-replicate record of `mrct_region_sims()` (Part C)

Date: 2026-11-01
Repository: **forestsearch** (`../forestsearch`), from `6b5c85c0` (0.4.1.9000, the end of `TASK_fs-mrct-analytic-layer_2026-11-01`). mrct-consistency at `4710058` is read only.
Basis: Part C of `PROPOSAL_forestsearch-mrct-functions_2026-11-01.md` (in `dev/tasks/`), approved 2026-11-01.

**`R/` category:** changes behaviour, additive only. The default call returns exactly today's columns. Nothing changes the method. CC does not push.

## 0. First action: the branch, then the task document

forestsearch's default branch is `master`. The three commits of the previous task sit on it, unpushed. Move them to a branch before anything else:

1. `git log --oneline 19b703b..master` must list exactly `5d45cd14`, `8c4ef5f6`, `6b5c85c0` and nothing else; `git status --porcelain` must be empty. Stop otherwise.
2. `git branch mrct-analytic`; `git branch -f master 19b703b`; `git checkout mrct-analytic`.
3. Record the three hashes and the result of `git log --oneline -1 master` (must be `19b703b`) in the notes file.

Then copy this file into `dev/tasks/` and commit it on its own, on `mrct-analytic`. Every commit of this task goes on that branch.

## 1. Backward compatibility (as in the previous task)

No existing export renamed, removed or re-signatured; no returned object changes unless a new argument asks for it. `record = "point"` (the default) returns exactly today's `results` table: same columns, same order, same values for the same seed.

## 2. The record

`mrct_region_sims(..., record = c("point", "full"))`. With `"full"`, the returned object carries, in addition to today's `results`, a second table `estimators` with one row per trial and estimator, and today's table gains the Ĥ-vs-H* counts. Match the committed campaign records in mrct-consistency `sims/123_ident_rlog4_full_records.R` (built by redrawing every trial), so that the record is the same thing computed once:

**`estimators`** (one row per trial × estimator), estimators `itt_all`, `itt_nonregion`, `itt_region`, `region_Ghat`, `region_Gstar`, plus `Hstar_nonregion`:
- `sim`, `estimator`, `n`, `n_exp`, `events_exp`, `events_ctl`;
- the unadjusted two-sample Cox fit on the set: `hr`, `lo`, `hi` (Wald 95%), `se_loghr` — the same `safe_coxph()` convention as today's point estimates (a fit with an arm without events is kept, not dropped);
- the targets on the same patients, from the DGM's columns: `ahr` = exp(mean `loghr_po`), `cde` = mean exp(θ₁) / mean exp(θ₀), as the records script defines them (`sims/two_settings_10k_design.R`, `tenk_fit()`);
- `declared` (0/1).

`region_Gstar` and `Hstar_nonregion` need `flag_harm`; when the DGM has none, those rows are omitted with a message, not an error. The margHR columns of the committed records (`marghr_inf`, `marghr_40`) come from the analytic layer and are **not** part of this record; the article computes them with `ahr_limit()`.

**Added to `results`** (today's per-trial table) under `"full"` only: `TP_nonregion`, `FP_nonregion`, `FN_nonregion`, `TN_nonregion`, `TP_region`, `FP_region`, `FN_region`, `TN_region` of Ĥ membership against `flag_harm`, as `hhat_confusion.csv` holds them (counts kept with Ĥ empty in an undeclared trial; NA where `flag_harm` is absent).

**Return shape.** Today `mrct_region_sims()` returns the `results` data.table. Keep that for `"point"`. For `"full"`, return the same `results` table with the `estimators` table attached as an attribute **and** as an element — choose one and document it; the constraint is that code doing `res$hr_itt` or `nrow(res)` on a `"point"` result is unaffected. A short `summaryout_mrct()`-style helper is not needed; the core OC rows are summaries the article computes itself.

## 3. Implementation notes (from source)

- The per-trial work already has `df_nonRegA`, `df_regA` and the full trial frame in the loop (`R/mrct_simulation.R` around the `df_itt`/`hr_test` computations); the sets are formed there from Ĥ (today's `get_dfpred()` route) and `flag_harm`. No redraw.
- `safe_coxph()`'s fit object gives `coef`, `se(coef)` and the 95% limits; take them from one `summary()` per set.
- Workers already return a per-trial `data.table`; add the estimator rows as a list column or a second combined table — CC's choice, as long as the `"point"` path builds exactly what it builds today.
- Also run a one-replicate smoke check of `sim_args = list(entry_delay = list(var = "z_region", delay = 6))` through `mrct_region_sims()` (the pass-through the previous notes flagged), and assert in a test that region entry lies in [6, 17] and non-region in [0, 11].

## 4. Tests (`tests/testthat/test-mrct-region-sims-record.R`)

On the synthetic DGM of `helper-synthetic-dgm.R` (or the small S2 build if faster to reuse), few replicates, fixed seed, 2 workers:
1. **Point path unchanged:** `record = "point"` gives a table `identical()` in names and values to a call without the argument, and to a call at `19b703b`'s code if a fixture is kept (CC may store one small fixture from the pre-change package; name its commit).
2. **Full path consistency:** `estimators$hr` for `itt_all`, `itt_nonregion`, `itt_region`, `region_Ghat` equals `results$hr_itt`, `hr_train`, `hr_test`, `hr_sg` where defined (1e-12).
3. **Counts:** per trial, `n_exp + (n − n_exp) = n`; the Ĥ confusion totals equal the set sizes; `region_Ghat` n = region n − (TP_region + FP_region).
4. **Targets:** `ahr` and `cde` recomputed from the trial frame agree (1e-12).
5. **Entry-delay pass-through** (§3).

## 5. Gate: legacy reproduction (stop on failure)

From mrct-consistency with the new package installed, the campaign call for S2 at a reduced replicate count and the redraw-based records script must agree with the new record: run `IF_REPS=50 IF_OUT=<scratch> Rscript sims/121_ident_rlog4_full_call.R S2_rmlog4` **with `record = "full"` added to the call through a copy of the script in the scratch directory** (mrct-consistency itself is not edited), then `sims/123_ident_rlog4_full_records.R` on the same scratch output, and compare the new `estimators` table to the script's `estimators.csv` on the shared columns (`hr`, `lo`, `hi`, `se_loghr`, `n`, `n_exp`, `events_exp`, `events_ctl`, `ahr`, `cde`) at 1e-9, and the confusion counts exactly. Record the largest difference.

Plus: `devtools::document()` clean; `devtools::install()`; the export list unchanged from `6b5c85c0`; the new test file passes. No full `R CMD check` (end of series).

## 6. NEWS

Add the `record` argument to the 0.4.1.9000 entry.

## 7. Notes file

`~/Downloads/NOTES_fs-mrct-region-sims-record_2026-11-01.md`: the branch step's hashes, status, commits, the gate-5 differences, test counts, wall-clock of the 50-replicate gate, and any disagreement with this task.

## 8. Commits (on `mrct-analytic`)

1. the task document;
2. the code, tests, documentation and NEWS.
