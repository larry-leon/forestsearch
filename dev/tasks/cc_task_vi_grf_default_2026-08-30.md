# CC task — make `vi.grf.min = NULL` the default in `forestsearch()`

**File:** `dev/tasks/cc_task_vi_grf_default_2026-08-30.md` · **Issued:** 2026-08-29 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** Larry places this file in `~/Downloads`; CC copies it into `dev/tasks/` and commits it first (§1).
**Run this only after `cc_task_oc_wrapper_confs_2026-08-30.md` has completed and committed.** Do not run the two concurrently.

**Larry's decision, 2026-08-29:** the GRF variable-importance step was intended to *order* covariates for search efficiency, not to alter the candidate set. At `vi.grf.min = -0.2` it does exactly that — `grf::variable_importance()` is non-negative, the block is guarded by `vi_max > 0`, so `vi_ratio ∈ [0,1]` and `which(vi_ratio > -0.2)` drops nothing. But the ordering it imposes is a per-replicate causal-forest fit, and column order has downstream reach (`extract_idx_flagredundancy()` walks `seq_len(ncol(x))` with the first factor tested against the whole sample). Disabling it by default makes the candidate space deterministic given the data.

---

## ⚠ `R/` CALLOUT — **CHANGES BEHAVIOUR**

**This changes the default of an exported function.** Category: *changes behaviour*, not a move and not a method change — the code paths already exist and neither is altered; only which one runs when the caller is silent.

Files that may be edited: `R/forestsearch_main.R` (the formal at ~L1216 and its roxygen at ~L712), `NAMESPACE`/`man/` as regenerated, `DESCRIPTION`, `NEWS.md`. **Nothing else.**

**Explicitly NOT to be changed in this task:** `R/run_simulation_analysis.R:68`, which sets `vi.grf.min = -0.2` inside `default_fs_params_*()`. That is an explicit caller choice, not a formal default; §2 reports it and Larry decides separately.

**Compute:** a bounded diagnostic — a handful of `forestsearch()` calls on one small fixture (§4). No simulation study, no replicates, no renders.

---

## 1. Provenance and first commit — GATE

```bash
hostname; pwd
git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD
git status --porcelain
git log --oneline -6
Rscript -e 'cat("forestsearch", as.character(packageVersion("forestsearch")), "\n")'
```

*GATE:* branch `feature/glm-extension`; clean tree; the confs task's commits present. Then copy this document into `dev/tasks/` and commit it alone (filename gate as usual — the `~/Downloads` stem arrives with hyphens stripped; copy under the name in the header and report both). Run `devtools::test()` for the parity baseline.

---

## 2. Blast radius — read-only — GATE

Establish, quoting each finding:

1. **Every caller of `forestsearch()` in the repo** — `R/`, `tests/`, `vignettes/`, `man/` examples, and **every `.qmd`/`.R` under `quarto/`** — and for each, whether it passes `vi.grf.min`. Tabulate: path, passes it (yes/no), value if yes.
2. **The applications specifically.** Do `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd` and `quarto/applications/actg175/analysis_actg175_continuous_compare_all.qmd` (and their siblings) pass it? These are the documents whose 0.2.0 reproduction was verified in `REPORT_repro_check_vs_0.2.0_2026-08-27.md`.
3. **Both `-0.2` sites**, verbatim with line numbers: the formal, and `run_simulation_analysis.R:68`.
4. **The `max_n_confounders` coupling.** Confirm from source that L2820's truncation sits inside `if (vi_max > 0)` and therefore does not execute when `vi.grf.min` is `NULL`. State plainly that the default change makes `max_n_confounders` inert by default, and list any caller that sets it to something other than 1000.
5. **RNG.** Does Section 5 consume R's RNG stream (as opposed to grf's internal `seed = seedit`)? If it does, skipping it shifts the stream for everything downstream, which is a far larger change than a reordering. Determine this from source and, if uncertain, by a bounded check: `.Random.seed` before and after a `forestsearch()` call at `-0.2` and at `NULL` on the same fixture.

*GATE — three outcomes:*

- **(A) No reference-producing caller relies on the default** (applications and any committed-payload driver all pass it explicitly) → proceed to §3.
- **(B) Some reference-producing caller relies on the default** → **stop before editing.** Write the report naming which, and what re-verification the change would require (the applications' 0.2.0 reproduction check is ~19 minutes of render). Commit and end; this becomes a compute go/no-go for Larry.
- **(C) §2.5 shows Section 5 consumes R's RNG** → **stop before editing.** That makes the change a stream shift rather than a reordering, and Larry needs to decide with that on the table.

---

## 3. The change — only under (A)

- `R/forestsearch_main.R`: the formal `vi.grf.min = -0.2` → `vi.grf.min = NULL`.
- Roxygen `@param vi.grf.min` (~L712): rewrite. State that `NULL` (the default) skips GRF variable-importance screening entirely; that a numeric value orders candidate covariates by importance and retains those with `vi_ratio > vi.grf.min`; that at values ≤ 0 nothing is filtered and the effect is ordering only; and — because it is not obvious — that `max_n_confounders` is applied **only** when screening runs.
- `@param max_n_confounders`: add the same caveat from its own side.
- `NEWS.md`: a plain entry — the default changed, screening is off unless requested, **results may change for callers that did not set `vi.grf.min`**, and callers wanting the previous behaviour should pass `-0.2` explicitly.
- `DESCRIPTION`: version → `0.3.0`. A default change that can move results is not a patch.

Do not touch `run_simulation_analysis.R`.

---

## 4. Characterise the change — bounded DIAGNOSTIC

On **one** small fixture (a synthetic DGM of the kind `tests/testthat` already uses, or the MD40 build at a small `n_super` — state which), with a fixed seed, run `forestsearch()` three ways and compare against the pre-change package (`git stash` the edit, or install to a temporary lib — say which you used):

| # | comparison | expected |
|---|---|---|
| 1 | `vi.grf.min = -0.2` explicit, post-change vs pre-change | `identical()` — explicit passers are unaffected |
| 2 | `vi.grf.min = NULL` explicit, post-change vs pre-change | `identical()` — the NULL path is unchanged |
| 3 | argument omitted, post-change vs `vi.grf.min = NULL` pre-change | `identical()` — the default now routes to NULL |

*GATE:* all three `identical()`. Any failure means the edit did more than change which branch runs — stop and report.

Then, for information only: argument omitted post-change vs omitted pre-change — report **what actually differs** (`sg.harm`, the candidate table's row count and order, `n_harm`, timing). This is the size of the behaviour change, and Larry should see it as a number rather than a description.

---

## 5. Close-out

- `devtools::document()`; `devtools::test()` — raw counts and **warning-count parity** against §1. Any test that relied on the old default will now behave differently: report each, do not adjust a test to make it pass without saying so explicitly and why.
- `R CMD check` (`RSTUDIO_PANDOC` on Rscript's PATH). Status verbatim; 0/0/0 is the target.
- Commit, **explicit paths only**, never `git add -A`. **No push.** `devtools::install()`; confirm `0.3.0`.
- Report `dev/glm-continuous-sims/REPORT_vi_grf_default_2026-08-30.md`: §2's caller table and gate verdict; the `max_n_confounders` and RNG findings; §3's diff; §4's three identity results and the measured size of the behaviour change; test and check output raw; commits; ten-line verdict.

---

## 6. Out of scope

- No edit to `run_simulation_analysis.R`, to any driver, to any application document, or to the OC-wrapper files.
- No re-run of the applications' 0.2.0 reproduction check — if §2 says it is needed, that is a separate compute go/no-go for Larry.
- No change to `max_n_confounders`' placement. Its coupling to the VI block is **reported as a finding**; moving the truncation out of that block would be a second behaviour change and is not authorized here.
- No simulation study, no replicate runs, no renders.
