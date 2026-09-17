# TASK — Two identifier fixes: GRF membership coding (P1) and DINA proposal-floor orientation (P2)

**File:** `dev/tasks/TASK_grf_dina_fixes_2026-09-16.md` · **Issued:** 2026-09-16 by chat, on Larry's yes to P1 and to P2 under its condition
**Machine:** `pop-os` · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`, HEAD containing `4b27f6f1`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first and alone (§1)
**Directory:** `quarto/simulations/actg175/continuous/`, written `<dir>`
**Record:** `<dir>/REPORT_grf_dina_fixes_2026-09-16.md`
**Evidence:**
- `<dir>/REPORT_md_grf_stage1_2026-09-16.md` §1.6(c) — the GRF membership defect
- `<dir>/REPORT_md_dina_grf_stage0_2026-09-16.md` S0.2 and S0.3 — the DINA floor defect
- `quarto/simulations/gbsg_020/REPORT_grf_factor_exposure_2026-09-16.md` — survival GRF not exposed; the replicate reproduction recipe (RNG kind pinned to L'Ecuyer-CMRG)

## ⚠ CATEGORY — touches `R/`

**P1: moves existing code and changes behaviour; no method change.**
- The problem: `.grf_evaluate_subgroup()` (`R/grf_subgroup_labels.R`) compares raw data columns with numeric cuts. The forest's covariate matrix codes factor covariates first, so candidates on factor covariates get NA membership and are dropped from re-selection and from the MR family.
- The fix: the evaluator codes covariates exactly as the forest's covariate matrix does, before applying cuts.
- If that coding is inline where the forest matrix is built, move it into one internal helper called by both, so the two codings are identical by construction.
- Behaviour changes only where factor covariates reach the evaluator. No committed result is exposed (the exposure record).

**P2: changes behaviour; no method change. Lands only under its condition (§3).**
- The problem: DINA's proposal floor `m_diff` is applied to the raw `tau_hat` (`R/forestsearch_helpers.R:1411`; `R/dina_subgroup.R:726`, `:865`). Under `adverse_outcome = FALSE` it is a benefit floor, while DINA's admission floor is a harm floor.
- The fix: apply the proposal floor to `tau_hat` in the orientation the admission floor uses, transplanting that orientation expression.
- Survival DINA and `adverse_outcome = TRUE` are unchanged by construction.

**Otherwise:**
- New unit tests, `NEWS.md` entries for what lands, and `devtools::document()` only if roxygen text changes.
- One final install into the main library.
- **Compute:** the guard fits of §2 (about 30 identifier fits, one to three replicates each) and two test-suite runs. Ceiling 3 h. Authorized by this kickoff.
- **Unattended.** Gates stop on failure, never to ask.
  - **On a stop:** leave `R/` edits uncommitted if their gate failed, write the record, commit the task document and the record, and stop.
  - A statement here that does not hold is a finding.

## Conventions

1. Verify from source; quote `path:line` at a stated commit.
2. `git add` by explicit path only. Untracked files are never staged. No `fetch`, `pull` or `push`.
3. **Code style:** tidyverse. Internal helpers are documented with roxygen `@noRd`.
4. **Tests** are fast and deterministic, with no network and no parallel backend. Use `skip_if_not_installed()` for suggested packages, and fit no forest where a constructed input suffices.
5. **Environment:** every R process runs with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`, a sequential future plan, and the seeds and RNG kind (L'Ecuyer-CMRG) the source campaign used.
6. **Identity** means `identical()` on the complete returned object, after removing elapsed-time and timing fields only. List the fields removed.
7. **Numbers** are computed and pasted, never typed. The submitted parent paper is not a topic.

---

## 1. Provenance, first commit, baseline — GATE

```bash
cd ~/Documents/GitHub/forestsearch
hostname; git branch --show-current; git rev-parse --short HEAD; git log --oneline -3
git status --porcelain --untracked-files=no
git status --porcelain -- R DESCRIPTION NAMESPACE tests
git merge-base --is-ancestor 4b27f6f1 HEAD && echo "exposure record in HEAD"
git diff --quiet 0071c17e..HEAD -- R/ DESCRIPTION NAMESPACE && echo "R/ unchanged since the install"
ps -eo pid,etime,args | grep -E '[e]xec/R|[R]script|[q]uarto' | head
Rscript -e 'cat(packageDescription("forestsearch")$Built, "\n")'
```

*GATE:* all of the following hold, or stop.
- The host is `pop-os` and the branch is `feature/glm-extension`.
- `4b27f6f1` is in HEAD, and `R/` is unchanged since `0071c17e`.
- The installed `Built` is `2026-09-16 05:57:14 UTC`.
- There are no tracked modifications, nothing pending under `R/`, `DESCRIPTION`, `NAMESPACE` or `tests/`, and no R, Rscript or quarto process running.

Copy this document from `~/Downloads` to `dev/tasks/TASK_grf_dina_fixes_2026-09-16.md` (exact name, else the single match for `*TASK_grf_dina_fixes_2026-09-16*.md`) and commit it alone.

**Baseline tests.** Run `devtools::test()` on HEAD under `timeout 90m`. Record passed, failed and skipped counts, and the name of every failing test. If it times out, run `devtools::test(filter = "grf|dina|forestsearch|mr|consistency|subgroup")` instead, record the fallback as a finding, and use that same filter for the §5 run.

## 2. Stage A — exposure list and before-captures (installed package, before any edit)

### 2.1 Committed DINA runs on binary or continuous outcomes

- **Locate them.** Run `git grep -n -E 'use_dina *= *TRUE|subgroup_method *= *"dina"'` across `quarto/`, `vignettes/`, `inst/`, `tests/` and `dev/`, excluding `dev/tasks/`.
- **For each call, record:**
  - file and line;
  - outcome type;
  - `adverse_outcome`, explicit or default (quote the default's source line);
  - whether it is a document, a simulation template with committed bundles, or a test.
- **Keep for §2.3:** every document or template call on a binary or continuous outcome with `adverse_outcome = FALSE`. Tests are listed separately; they are code, not results.

### 2.2 P1 call sites

- List every call site of `.grf_evaluate_subgroup()` and of the forest-matrix coding it will share.
- State whether any call site lies on GRF's screening-cut path or the standalone tree selection. The exposure record says none does; confirm or contradict.

### 2.3 Before-captures

In a temporary directory outside the repo, using the installed package, run each fit below and save its returned object as RDS.

| Fit | What | MR |
|---|---|---|
| F1 | Survival GRF: `grfmr` cell A124_h150_n500, sim_id 1–3, arguments from the `grfmr` runner (quoted), regenerated by the exposure record's recipe | as `grfmr` ran it |
| F2 | Survival DINA: one committed `dinamr` harm cell (name it), sim_id 1–3, arguments from the `dinamr` runner (quoted) | as `dinamr` ran it |
| F3 | FS on the MD design: md40 n500, sim_id 1–3, `mdsgnb20`'s arguments | on |
| F4 | GRF on the MD design: md40 n500, sim_id 1, the MD template's GRF argument block (commit `894da993`, quoted) | off |
| F5 | DINA on the MD design: md40 n500, sim_id 1, `adverse_outcome = FALSE`, the settings of Stage 0 S0.3 | off |
| F6 | As F5, with the outcome negated and `adverse_outcome = TRUE` | off |
| F7 | Each §2.1 call kept for review, reproduced from its document's setup and data chunks and that call (via `knitr::purl`), or from sim_id 1 of its template | as committed |

- **Warnings.** Capture them for every fit with `withCallingHandlers`.
- **F4:** record the number of candidates enumerated, the number with NA membership, and the number admitted.
- **F5:** record the proposal floor and admission floor as applied, the number of candidates proposed, and the selection.
- **F7:** a call that cannot be reproduced outside its document is recorded as "not reproducible here".

## 3. Implement P1 — GATE

1. **Edit.** Implement P1 as the category box specifies, quoting the forest-matrix coding lines you transplant or move. Add unit tests:
   - a constructed data frame with a binary factor and a numeric covariate;
   - a candidate cut on each;
   - assert membership equals the expected logical vector, with no NA and no warning;
   - a numeric-only case returning what it returned before.
2. **Temporary install.** Install the working tree into a temporary library (`R CMD INSTALL --library=<tmp>/lib_p1 .`). Re-run F1–F7 from fresh R processes loading that library.
3. *GATE P1:*
   - F1, F2, F3, F5, F6 and F7 are identical to their before-captures.
   - F4 has zero factor-comparison warnings and zero NA-membership candidates. Report F4 before and after: enumerated, NA, admitted, and the selection.
   - For every factor covariate of F4's data, the evaluator's coded column is `identical()` to the forest matrix's column.
   - The new tests pass.
4. **Commit** the P1 `R/` change and its tests together, by explicit path.

## 4. Implement P2 — GATE and condition

1. **Edit.** Implement P2 on top of P1, as the category box specifies, quoting the admission floor's orientation lines you transplant. Add unit tests:
   - under `adverse_outcome = FALSE`, a constructed surface where the proposal set must lie on the harm side of `m_diff`;
   - under `adverse_outcome = TRUE`, the set is unchanged.
2. **Temporary install.** Install into `<tmp>/lib_p2`. Re-run F1–F7 from fresh processes.
3. *GATE P2 — correctness:*
   - F1, F2, F3 and F6 are identical to their before-captures, and F4 is identical to its after-P1 capture.
   - F5 proposes a non-empty set, every proposed candidate satisfies oriented `tau_hat ≥ m_diff`, and its proposal and admission floors point the same way.
   - F5's selection (definition and n) equals F6's, since orientation is the only difference between them.
   - If any of these fails, P2 is incomplete. Do not land it: remove the P2 edits from the working tree, keep the diff as `<dir>/P2_dina_floor_orientation.patch`, and continue at §5 with P1 only.
4. **Condition — Larry's.**
   - If every F7 fit is identical to its before-capture, land P2.
   - If any differs, or any is "not reproducible here", do not land P2. Remove the edits, keep the patch as above, and record the list with each call's before and after selection.
   - Tests whose expectations change under P2 are not results. If P2 lands, update them in the P2 commit and quote the old and new expectations.
5. **If P2 lands,** commit the P2 `R/` change and its tests by explicit path.

## 5. Tests, news, final install — GATE

1. Run `devtools::test()` (or the §1 filter) under `timeout 90m`.
   - *GATE:* every test that passed at baseline passes, and every new test passes. Failures present at baseline may remain; list them.
2. Add `NEWS.md` bullets under the development heading, one per landed fix, stating the defect and its effect in one sentence each.
3. Run `devtools::document()` only if roxygen text changed. Commit `NEWS.md` and any `man/` changes by explicit path.
4. **Final install:** `devtools::install(quick = TRUE, upgrade = FALSE)`.
   - Record `Built`, and assert that two doFuture workers report the same value.
   - Re-run F4 (and F5 if P2 landed) against the main library; each is identical to its capture from the temporary library.
5. Remove the temporary directory and libraries.

## 6. Record, catalog, closeout

1. **Record,** at the header path. It contains:
   - provenance and the baseline tests;
   - §2.1 as a table, and the §2.2 call sites;
   - the fit table: F1–F7, before / after-P1 / after-P2, with identity results and the F4/F5 counts;
   - the P1 and P2 diffs, quoted;
   - the condition's outcome; if P2 did not land, the list and the patch path;
   - test results, final `Built`, the landed commits, findings.

   Commit it, and the patch file if one was written.
2. **Catalog.** In `<dir>/status_curated.md`, update the open-work lines:
   - GRF: "membership fix landed (commit); campaign `mdgrf` resumes from its smoke".
   - DINA: "proposal-floor fix landed (commit); campaign pending", or "P2 not landed: see REPORT_grf_dina_fixes_2026-09-16.md".

   Commit, regenerate `current_status.md` as the last commit, and check that `check_current_status.sh --commit` passes.
3. **Post-conditions,** printed in the closing message:
   - `git diff --name-only <§1 HEAD>..HEAD -- R/` lists only the files P1 and a landed P2 changed;
   - the full diff otherwise lists only the task document, tests, `NEWS.md`, any `man/` files, the record, the patch if any, and the two catalog files;
   - there are no tracked modifications;
   - the installed `Built` equals §5's;
   - the temporary directory and libraries are gone.
4. **Closing message:**
   - the commit range to push;
   - P1 landed, with F4 before and after;
   - P2 landed or not, and if not, why, with the list;
   - the guard identities;
   - test results;
   - the final `Built`;
   - findings.

   Then stop.

## Out of scope

- `.fs_apply_mr()`'s IJ fallback, `dmin.grf` behaviour, and any other `R/` change.
- Campaigns; re-rendering any document; version bumps; `R CMD check`; pushing.
