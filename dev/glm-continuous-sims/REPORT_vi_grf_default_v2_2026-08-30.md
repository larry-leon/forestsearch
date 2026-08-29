# REPORT — `vi.grf.min = NULL` is the default (0.3.0), with a guard for the silent cap

**Task:** `dev/tasks/cc_task_vi_grf_default_v2_2026-08-30.md` (commit `c946def6`), superseding `cc_task_vi_grf_default_2026-08-30.md` (stopped at (B), no edit — that finding stands; its consequence was measured in `REPORT_vi_grf_smoke_2026-08-30.md`).
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Category:** **changes behaviour** — the default of an exported function's formal. Files edited: `R/forestsearch_main.R` (the formal, its roxygen, the `max_n_confounders` roxygen, one `warning()` before Section 5), `man/forestsearch.Rd` (regenerated), `DESCRIPTION` (0.2.7 → **0.3.0**), `NEWS.md`. `NAMESPACE` unchanged (no new export). **Not touched:** `R/run_simulation_analysis.R:68`, any driver, application document, payload or OC-wrapper file; nothing already committed re-run, re-rendered or re-verified.

---

## 1. Provenance (§1 — GATE passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at start e5b15126 (the smoke-test report commit; merge-base check: yes) · git status --porcelain: empty
git log -6: e5b15126 4bc2bc76 01969fa2 695ecd39 81c9568a d5ef8c7c
packageVersion forestsearch 0.2.7 (installed 0.2.7)
vi.grf.min formal before edit: R/forestsearch_main.R:1216 `vi.grf.min = -0.2`
```

`ls ~/Downloads/cc_task_vi_grf_default_v2_2026*`: one match, `cc_task_vi_grf_default_v2_20260830.md` (hyphens stripped) → `dev/tasks/cc_task_vi_grf_default_v2_2026-08-30.md`, committed **`c946def6`**.

Parity baseline `devtools::test()` on the unchanged tree (started before any edit): `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4844 ]`

Pre-change reference captured before any edit (§4): the 300-row synthetic continuous fixture (`.make_continuous_data(N = 300, MD_harm = 1.5, seed = 42)`), `seedit = 1`, `maxk = 2`, `n.min = 30`, `fs.splits = 30`, `pconsistency = 0.5`, `sg_focus = "maxeffCons"`, sequential; three calls — `vi.grf.min = -0.2` explicit, `NULL` explicit, argument omitted — with timings and environment-bearing fields stripped (`minutes_all`, `find.grps$time_search`, `grf_res`, `args_call_all$df.analysis`, `$parallel_args`).

---

## 2. The change (`git diff R/forestsearch_main.R`, code lines)

```diff
-                         vi.grf.min = -0.2,
+                         vi.grf.min = NULL,
```

Roxygen `@param vi.grf.min` rewritten (states: `NULL` skips Section 5, no forest fitted; a numeric value fits the forest, **orders** the cut columns and retains `vi_ratio > vi.grf.min`; at ≤ 0 nothing is filtered — ordering only; the ordering does not reach selection because the table is re-sorted by `(-HR, K)` before the consistency loop, measured effect on `n_candidates_evaluated` / `early_stop_candidate` zero, membership unchanged, only `sg.harm` clause order may differ; pre-0.3.0 default was `-0.2`, pass it explicitly to reproduce; `max_n_confounders` applies only when screening runs). Roxygen `@param max_n_confounders` rewritten from its side (applied only when `vi.grf.min` is numeric; no effect under `NULL`; finite cap below the column count warns).

### The guard (§3)

Inserted immediately before Section 5, after `FSconfounders.name <- FSdata$confs_names` (so the exact number of candidate cut columns the cap would apply to is known — no restructuring needed):

```r
  if (is.null(vi.grf.min) && is.finite(max_n_confounders) &&
      max_n_confounders < length(FSconfounders.name)) {
    warning("max_n_confounders = ", max_n_confounders, " is inert: variable-",
            "importance screening is off (vi.grf.min = NULL), so all ",
            length(FSconfounders.name), " candidate cut columns are searched. ",
            "Set vi.grf.min to a numeric value (e.g. -0.2, the pre-0.3.0 ",
            "default) to apply the cap.", call. = FALSE)
  }
```

One `warning()`, suppressible normally; the cap's placement (now L2858, inside the VI block under `if (vi_max > 0)`) is unchanged.

---

## 3. §4 — characterisation (GATE passed)

Fixture as in §1, 8 candidate cut columns (`age <= 55 | 47 | 62`, `biomarker <= -0.1 | -1 | 0.6`, `biomarker_hi`, `sex`); the post-change reference captured with the identical script after `document()`.

| # | comparison | result |
|---|---|---|
| 1 | `vi.grf.min = -0.2` explicit, post vs pre | **`identical()` = TRUE** |
| 2 | `vi.grf.min = NULL` explicit, post vs pre | **`identical()` = TRUE** |
| 3 | argument omitted post vs `vi.grf.min = NULL` explicit pre | **`identical()` = TRUE** (including `args_call_all$vi.grf.min`, which now records `NULL`) |

*GATE: all three identical — the edit changed only which branch the default routes to.*

**Size of the behaviour change** (argument omitted, post vs pre — i.e. `NULL` vs `-0.2`):

| quantity | omitted, post (`NULL`) | omitted, pre (`-0.2`) |
|---|---|---|
| `sg.harm` | `{age <= 47} {biomarker_hi}` | `{biomarker_hi} {age <= 47}` — **clause order only** |
| selected subjects (`sg.harm.id`) | identical, n = 36 | n = 36 |
| effect (`out_sg$result$hr`) | 1.094112 | 1.094112 |
| candidate table rows | 7 | 7 — same rows, same column set, **different column order** |
| `n_candidates_evaluated` | 1 | 1 |
| wall-clock | 0.21 s | 0.28 s |
| top-level fields that differ | `grp.consistency` (clause order / column order), `find.grps` (column order), `sg.harm` (clause order), `args_call_all` (`vi.grf.min`) | |

This is the smoke report's result reproduced on one fixture: membership, effect and counters identical; clause order and column order differ; ~0.07 s saved.

**Guard behaviour** (same fixture, 8 cut columns; warnings captured with `withCallingHandlers`):

| call | inert-cap warnings |
|---|---:|
| `vi.grf.min = NULL`, `max_n_confounders = 3` | **1** — "max_n_confounders = 3 is inert: variable-importance screening is off (vi.grf.min = NULL), so all 8 candidate cut columns are searched. Set vi.grf.min to a numeric value (e.g. -0.2, the pre-0.3.0 default) to apply the cap." |
| omitted (new default), `max_n_confounders = 3` | **1** (same message) |
| `vi.grf.min = -0.2`, `max_n_confounders = 3` | 0 — the cap applies |
| `NULL`, cap 1000 (default) | 0 |
| `NULL`, cap 8 (= number of columns) | 0 — does not bind |
| `NULL`, cap 7 (< columns) | 1 |
| `NULL`, cap `Inf` | 0 |

(Every call also emits one unrelated warning from the fixture itself — `parallel_args missing required elements. Using sequential.`, from passing `parallel_args = list(plan = "sequential")` — counted separately.)

---

## 4. Close-out

`devtools::document()`: regenerated `man/forestsearch.Rd` only.

`devtools::test()` after the change: `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4844 ]` — parity exact against the baseline (WARN 31 = 31, SKIP 3 = 3, PASS 4844 = 4844). The five test files that call `forestsearch()` without `vi.grf.min` (`test-glm-pipeline-integration.R`, `test-input-validation.R`, `test-mr-inference.R`, `test-return-shape-contract.R`, `test-search-reproducibility.R`) now take the `NULL` path; none of their assertions changed outcome, and no test was adjusted. The new guard fired in no test (no test sets a finite `max_n_confounders`).

`R CMD check` (`devtools::check(document = FALSE, args = "--no-manual")`, `RSTUDIO_PANDOC` on `PATH`, pandoc 3.8.3): 

```
── R CMD check results ───────────────────────────────── forestsearch 0.3.0 ────
Duration: 10m 18.3s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
Status: OK
```


`DESCRIPTION` 0.2.7 → **0.3.0**. `NEWS.md` states all four required items (default changed and screening off unless requested; **reproducing pre-0.3.0 results requires `vi.grf.min = -0.2` explicitly**; the measured effect where the cap does not bind — identical memberships, clause order may differ; `max_n_confounders` has no effect under `NULL`, now with a warning). Staged by explicit path (`R/forestsearch_main.R`, `man/forestsearch.Rd`, `DESCRIPTION`, `NEWS.md`, this report). **No push.** `devtools::install(upgrade = FALSE)`: `installed forestsearch 0.3.0`.

Commits: `c946def6` task doc; the change/report commit and its hash-recording follow-up below.

---

## 5. Verdict (ten lines)

1. `vi.grf.min` default `-0.2` → **`NULL`** (`R/forestsearch_main.R`, one line); version **0.3.0**.
2. Both roxygen entries rewritten to say what each setting does, that the ordering never reaches selection, and that `max_n_confounders` applies only when screening runs.
3. One guard added before Section 5: a finite `max_n_confounders` below the cut-column count under `NULL` warns that it is inert and how to restore it; the cap's placement is unchanged. Verified: fires in exactly the binding cases, silent otherwise.
4. §4 gate: explicit `-0.2`, explicit `NULL`, and omitted-vs-`NULL` all `identical()` pre vs post.
5. Size of the change on the fixture: same 36 subjects, same effect 1.094, same 7 candidates and evaluated count; `sg.harm` clause order and the candidate table's column order differ; 0.07 s faster.
6. `devtools::test()` `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4844 ]` before and after — exact parity, no test adjusted, the guard fired in no test; `R CMD check` 0 / 0 / 0 at 0.3.0.
7. `R/run_simulation_analysis.R:68` (`-0.2`) untouched; drivers, applications, payloads, OC-wrapper files untouched; nothing re-verified.
8. Reproducing pre-0.3.0 output now requires `vi.grf.min = -0.2` explicitly — stated in `NEWS.md` and in the roxygen.
9. The historical record (applications' 0.2.0 reproduction, `gbsg_redux` payloads) stands as a statement about earlier versions.
10. Installed 0.3.0; no push.
