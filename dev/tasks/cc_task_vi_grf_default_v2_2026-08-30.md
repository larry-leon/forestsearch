# CC task (v2) — make `vi.grf.min = NULL` the default, on the smoke-test evidence

**File:** `dev/tasks/cc_task_vi_grf_default_v2_2026-08-30.md` · **Issued:** 2026-08-30 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).

**Supersedes `cc_task_vi_grf_default_2026-08-30.md`**, which stopped at its §2 gate with outcome (B) — reference-producing callers do rely on the default — and made no edit. That finding stands. What has changed is that its *consequence* has now been measured, in `REPORT_vi_grf_smoke_2026-08-30.md`:

- **Equivalence.** 260 paired runs: 132 identical, 116 clause-order only, 12 differing substantively — **all 12 in F4, the `max_n_confounders` positive control**. Candidate families are identical as sets of subject memberships in **220/220** non-F4 pairs; where cut-sets differ, **100% of the differing rows are nested same-variable encodings** (`{age <= 54} & {age <= 63}`) of an identical subject set. F2, the redundancy-walk case, never changed the selected subgroup (0/40).
- **Efficiency: none.** The candidate table is re-sorted by `(-HR, K)` twice before the consistency loop (`format_search_results()` L883, `sort_subgroups_preview()` L660), so the VI order never reaches the scan: `n_candidates_evaluated` and `early_stop_candidate` differ by **exactly 0 in all 240 early-stopping pairs**, including a strict `pconsistency = 0.9` arm. `-0.2` costs **+0.13 s per call** (0.06–0.21, 3–11 SE above zero), essentially all of it the causal forest, with no search-side saving anywhere.
- **The one real hazard.** When `max_n_confounders` binds, the two settings genuinely differ — the cap keeps the top-importance columns and changes the winner in ~30% of seeds — and under `NULL` the cap becomes **silently inert**, invisible in `confounders.evaluated` (pre-cap on both sides). **No repository caller sets it.** §3 below adds a guard rather than leaving it silent.

**Larry's decision, 2026-08-30:** make `NULL` the default. **Nothing already done is touched** — no driver, no application document, no payload, no re-verification, no re-render. The change is forward-looking; the historical record stands as a statement about earlier versions.

---

## ⚠ `R/` CALLOUT — **CHANGES BEHAVIOUR**

Editable: `R/forestsearch_main.R` (the `vi.grf.min` formal, its roxygen, the `max_n_confounders` roxygen, and the §3 guard), `man/`, `NAMESPACE` as regenerated, `DESCRIPTION`, `NEWS.md`. **Nothing else.**

**Explicitly NOT changed:** `R/run_simulation_analysis.R:68` (an explicit caller choice, not a formal default); any driver; any application document; any payload; the OC-wrapper files.

**Compute:** §4's characterisation only — a handful of `forestsearch()` calls on a small synthetic fixture. No simulation study, no renders.

---

## 1. Provenance and first commit — GATE

Standard block plus the installed version. *GATE:* branch `feature/glm-extension`, clean tree, HEAD descends from `e5b15126` (the smoke-test report commit), and the `vi.grf.min` formal still reads `-0.2`. Copy this document into `dev/tasks/` and commit alone; the `~/Downloads` stem arrives with hyphens stripped, so report both names. Run `devtools::test()` for the parity baseline and **capture the pre-change reference** described in §4 before any edit.

---

## 2. The default

`R/forestsearch_main.R`: the formal `vi.grf.min = -0.2` → `vi.grf.min = NULL`.

Roxygen `@param vi.grf.min`, rewritten to state:

- `NULL` (the default) skips GRF variable-importance screening entirely — no forest is fitted.
- A numeric value fits a causal forest, **orders** the candidate covariates by importance, and retains those with `vi_ratio > vi.grf.min`. At values ≤ 0 nothing is retained-out: `grf::variable_importance()` is non-negative and the block is guarded by `vi_max > 0`, so the effect is ordering only.
- The ordering sets the candidate matrix's column order. It does **not** reach candidate selection: the table is re-sorted by `(-HR, K)` before the consistency loop. Measured effect on `n_candidates_evaluated` and `early_stop_candidate`: zero.
- `max_n_confounders` is applied **only when screening runs** — see its own entry.

Roxygen `@param max_n_confounders`: state the same coupling from its side, plainly — the cap has no effect when `vi.grf.min` is `NULL`.

---

## 3. The guard for the silent cap

The one measured hazard is a caller who sets a finite `max_n_confounders` and, under the new default, silently loses it.

Add a single `warning()` in `forestsearch()`, at argument resolution before Section 5: when `vi.grf.min` is `NULL` **and** `max_n_confounders` is finite and less than the number of candidate cut columns that will be searched — or, if that count is not available at that point, simply when it is finite and not the default — warn that the cap is inert because variable-importance screening is off, and that setting `vi.grf.min` to a numeric value restores it.

Keep it to one warning, suppressible in the normal way, and worded so a user can act on it. **Do not change the cap's placement** — moving the truncation out of the VI block would be a second behaviour change and is not authorized here.

If this guard cannot be added without restructuring existing logic, **skip it, say so in the report, and proceed with §2** — the default change is the deliverable; the guard is a mitigation.

---

## 4. Characterise the change — GATE

Before editing, capture the pre-change reference; after editing, compare. On one small synthetic fixture with a fixed seed:

| # | comparison | expected |
|---|---|---|
| 1 | `vi.grf.min = -0.2` explicit, post vs pre | `identical()` — explicit passers unaffected |
| 2 | `vi.grf.min = NULL` explicit, post vs pre | `identical()` — the NULL path unchanged |
| 3 | argument omitted post vs `vi.grf.min = NULL` pre | `identical()` — the default now routes to NULL |

*GATE:* all three `identical()`. Any failure means the edit did more than change which branch runs — stop and report.

Then, for information: argument omitted post vs omitted pre, reporting what differs (`sg.harm` clause order, candidate-table rows, `n_harm`, timing) — the size of the behaviour change as a number.

Also verify §3's guard fires when it should and stays silent when it should not.

---

## 5. Close-out

- `devtools::document()`; `devtools::test()` — counts and **warning-count parity** against §1. A test that relied on the old default will now behave differently: report each explicitly; **do not adjust a test to pass without saying so and why.** The parity baseline is 31 warnings — if the new `warning()` in §3 fires inside any test, parity will break legitimately: say which test and why.
- `R CMD check` (`RSTUDIO_PANDOC` on Rscript's PATH). Status verbatim; 0/0/0 is the target.
- `DESCRIPTION` → **0.3.0**. A default change that can move results is not a patch.
- `NEWS.md`, stating all four of these:
  1. the default changed from `-0.2` to `NULL`; screening is off unless requested;
  2. **reproducing results from versions before 0.3.0 requires passing `vi.grf.min = -0.2` explicitly**;
  3. the measured effect where the cap does not bind — identical subject memberships, clause order in `sg.harm` may differ;
  4. `max_n_confounders` has no effect when `vi.grf.min` is `NULL`.
- Commit, **explicit paths only**, never `git add -A`. **No push.** `devtools::install()` (`upgrade = FALSE`); confirm `0.3.0`.
- Report `REPORT_vi_grf_default_v2_2026-08-30.md`: the diff; §4's three identity results and the measured size of the change; the guard's behaviour; test and check output raw; commits; ten-line verdict.

---

## 6. Out of scope

- No edit to `R/run_simulation_analysis.R`, to any driver, to any application document, or to the OC-wrapper files.
- **No re-verification, re-render or re-run of anything already committed.** The `gbsg_redux` drivers and their 370 payloads are untouched; the applications' 0.2.0 reproduction is a historical statement about earlier versions and is not re-tested here.
- No change to `max_n_confounders`' placement in the code.
- No simulation study, no renders, no push.
