# CC task — `vi.grf.min = -0.2` vs `NULL`: equivalence of results, and whether the ordering buys any efficiency

**File:** `dev/tasks/cc_task_vi_grf_smoke_2026-08-30.md` · **Issued:** 2026-08-30 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Follows:** `REPORT_vi_grf_default_2026-08-30.md`, which stopped at its §2 gate with outcome (B) and made no edit. This task is the evidence that gate asked for, in a form that touches nothing already done.

**Two questions, and the second is the one that decides.**

1. **Equivalence.** Does `vi.grf.min = -0.2` (order only — nothing is filtered, `variable_importance()` being non-negative and the block guarded by `vi_max > 0`) give the same *result* as `NULL` (Section 5 skipped entirely)?
2. **Efficiency.** The ordering was built to make the search faster. **Does it?** Larry does not want it removed as a default until that is assessed, and it has never been measured.

**Nothing historical is touched.** No driver, no application, no payload, no re-verification of anything already committed. Synthetic data only.

---

## ⚠ CATEGORY

**No `R/` change of any kind.** No new export, no default change, no edit to any package file, driver or document. This task produces evidence; the default change is a later, separate decision.

Writes: one scratch script, its outputs, and the report, under `dev/glm-continuous-sims/` (or `dev/profiling/` — say which). Plus this task document.

**Compute:** bounded diagnostic on small synthetic data. Estimate 45–75 minutes. No simulation study, no bootstrap campaign, no renders.

**Unattended-safe.** Gates stop the task; never ask, never work around.

---

## 1. Provenance and first commit — GATE

Standard block plus the installed version and the **current `vi.grf.min` default in force** (it should still be `-0.2`; the default task stopped without editing). Copy this document into `dev/tasks/` and commit alone; the `~/Downloads` stem arrives with hyphens stripped, so report both names.

---

## 2. The comparison primitive — read this before writing anything

**Compare subject membership, never the definition string.** The VI reordering changes the clause order in `sg.harm` (measured in the prior report), so a string comparison flags every case and tells you nothing. The primitives are:

| what | how |
|---|---|
| selected subgroup | `identical(sort(which(sg.harm.id == 1)), ...)` — the subject index set |
| effect | the fitted estimate and CI on the selected subgroup |
| candidate family | `find.grps$out.found$hr.subgroups` compared **as a set of rows**, order-insensitive |
| consistency stage | `grp.consistency$n_candidates_total`, `n_candidates_evaluated`, `n_passed` |
| declaration | `is.null(sg.harm)` — did anything get selected at all |

Classify each paired run as: **identical** · **clause-order only** (membership identical, string differs) · **differs substantively** (membership differs). Report the counts.

---

## 3. Fixtures — designed to break it, not to pass

A test that only exercises benign configurations proves nothing. Each fixture below targets a specific mechanism, and **two of them are expected to differ** — if those two come back identical, the harness is not sensitive and the whole result is void.

Synthetic data throughout; small `n` (200–600) and few covariates, so each `forestsearch()` call is seconds. State the generator.

| # | fixture | mechanism targeted | expectation |
|---|---|---|---|
| F1 | plain continuous DGM, ~8 covariates, one real subgroup | baseline | identical or clause-order only |
| F2 | includes a covariate whose cut yields prevalence just below `(n − rmin)/n` (≈0.99 at n = 500, `rmin = 5`) | `extract_idx_flagredundancy()` tests the **first** factor against the whole sample, so a near-full-sample cut is flagged redundant only when it comes first — the one demonstrated order-dependence | **may differ substantively; this is the sensitive case** |
| F3 | two covariates with (near-)identical membership | tie-breaking and statistics-keyed dedup | may differ |
| F4 | `max_n_confounders` set **low** (e.g. 5) with ~10 covariates | the cap is applied **only inside the VI block**, so it truncates at `-0.2` and is inert at `NULL` | **expected to differ substantively — this is the positive control** |
| F5 | all covariates pure noise (VI ≈ 0 for all) | `vi_max > 0` fails, block skipped; `order()` on an all-equal vector is stable | expected identical |
| F6 | survival outcome with zero events in one arm | the survival guard skips the VI step entirely (`cs.forest <- "zero-events"`) | expected identical |
| F7 | survival DGM, normal | the `causal_survival_forest` path rather than the GLM one | identical or clause-order only |
| F8 | binary outcome | the third VI path | identical or clause-order only |

**≥ 20 seeds per fixture.** One seed proves nothing; report the per-fixture rate of each classification.

*GATE:* if **F4** does not differ substantively, the harness is not measuring what it claims — stop and report.

---

## 4. The efficiency question — the one that decides

For every paired run, at both settings, record:

- **total wall-clock** of the `forestsearch()` call;
- **time attributable to Section 5** (the GRF fit) — time it directly, or difference the two settings on a fixture where nothing else changes;
- `grp.consistency$n_candidates_total`, **`n_candidates_evaluated`**, `n_passed`, `early_stop_triggered`, **`early_stop_candidate`**;
- the number of candidates reaching a model fit in `subgroup.search()`.

**The hypothesis to test, stated so it can fail.** The ordering was intended to reach the consistency loop sooner by putting high-importance covariates first. Chat's reading is that it *cannot*, because `format_search_results()` applies `setorder(hr.out, -HR, K)` and `sort_subgroups_preview()` re-sorts before the consistency loop, so the scan order at the early-stop stage is fitted effect, not VI order. **Verify this from source, quoting it — and then measure it**: if the ordering helps, `n_candidates_evaluated` and `early_stop_candidate` should be systematically *smaller* at `-0.2` than at `NULL`. Report the paired differences and whether they differ from zero.

Use `sg_focus = "maxeffCons"` for the early-stopping arm — the source (`forestsearch_main.R:1574`) says `stop_threshold` is meaningful for that focus alone and is reset to `NULL` for every other. Run at least one fixture under a non-`maxeffCons` focus as well, and note that early stopping cannot fire there.

**Report the net.** Per fixture: mean wall-clock at `-0.2` minus at `NULL`, decomposed into the forest-fit cost and any search-side saving. A negative total means the ordering pays for itself; a positive total means it does not.

---

## 5. Verdict

Three findings, stated separately:

1. **Equivalence** — the classification counts across all fixtures and seeds, with F2 and F4 called out by name.
2. **Efficiency** — does the ordering reduce `n_candidates_evaluated` or wall-clock? Yes with a magnitude, or no with the measured difference and its uncertainty.
3. **The `max_n_confounders` coupling** — from F4, what changing the default would do to a caller relying on the cap. State it plainly; it is the one behaviour that genuinely changes.

No recommendation. The default change is Larry's decision on this evidence.

---

## 6. Report

`REPORT_vi_grf_smoke_2026-08-30.md`: provenance · the source verification of §4's re-sort claim, quoted · fixture generators · the per-fixture classification table · the efficiency table with paired differences · the three verdicts · ten-line summary. Commit the script, its outputs and the report by explicit path. **No push.** No `devtools::install()`.

---

## 7. Out of scope

- No `R/` change, no default change, no edit to any driver, application or document.
- No re-run or re-verification of anything already committed — no application renders, no touching the `gbsg_redux` drivers or their payloads.
- No bootstrap campaign, no simulation study, no renders.
- No recommendation about whether to change the default.
