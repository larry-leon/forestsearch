# CC task — deprecate `vi.grf.min` and `max_n_confounders`; remove GRF variable-importance screening from the forest-search path

**File:** `dev/tasks/cc_task_vi_deprecation_2026-10-26.md` · **Issued:** 2026-10-26 by chat, commissioned by Larry (decision: option (3) of the `vi.grf.min` discussion — deprecate at the package level, so that `use_lasso = use_grf = use_dina = FALSE` means no model is fitted anywhere in the forest-search path)
**Repo:** `forestsearch` · **Branch:** new, `feature/vi-deprecation-lean-fits`, cut from `master` in §1. All commits of this task and of the lean-fit (A+F) task that follows go on it. CC never fetches, pulls, pushes or merges; Larry merges in GitHub Desktop.
**Transport:** `~/Downloads` → `dev/tasks/` (this document, committed first and alone). Filenames may arrive with hyphens stripped; locate by pattern.
**Audit basis:** CRAN `forestsearch` 0.4.0. Every line number below is a 0.4.0 line number; `master`'s `R/` files are byte-identical to 0.4.0 (Stage 0b MD5 record: `R/forestsearch_main.R` 11d80e74…, `R/get_fsdata.R` a7f71d65…), so they are the current line numbers until this task's own edits. Verify each named line by content before editing it; if a line has moved, find it by content and say so in the report.

---

## ⚠ CATEGORY — `R/` change: CHANGES BEHAVIOUR (documented API deprecation), decided by Larry

- **What changes.** `forestsearch()` keeps the formals `max_n_confounders = 1000` and `vi.grf.min = NULL` (same positions, same defaults — no existing call breaks), but a non-default value now emits **one deprecation warning and has no effect**. The GRF variable-importance screening/ordering block (0.4.0 "SECTION 5") is **deleted**: no causal forest is fitted in the forest-search path unless `use_grf`, `use_dina` or `use_lasso` asks for one. The formals are removed in a later release.
- **Who sees a change.** Only callers that passed a non-NULL `vi.grf.min` (or a non-default `max_n_confounders`). At `vi.grf.min <= 0` the forest never filtered anything, so for those callers the selected *membership* is unchanged by construction; the clause order inside the printed `sg.harm` label and the redundancy column-walk can differ (the 0.4.0 record says "ordering only"). Callers with a binding `vi.grf.min > 0` lose the screen — the intent. Default callers: bit-identical, verified in §4.
- **Not a method change.** The method's candidate family is the cut grid; this removes a pre-0.3.0 screening step that contradicted the fixed-family statement.
- **Deprecated values are not reset.** `args_call_all` (`mget(names(formals()))`, L1735) records what the caller passed, so a fitted object still says what was asked for; `fs_family_report()` echoes it under status `inert`. Because `do.call(forestsearch, args_call_all)` in the bootstrap and CV paths passes every formal explicitly, replicates of a *default* call never warn; replicates of a call that passed `-0.2` warn inside the existing `suppressWarnings()` (`bootstrap_analysis_dofuture.R` L408; `forestsearch_cross_validation.R` L464/L479/L950/L1001), i.e. silently — acceptable.
- **Writes:** `R/forestsearch_main.R`, `R/fs_family_report.R`, `R/run_simulation_analysis.R`, `R/get_fsdata.R`; regenerated `man/forestsearch.Rd`, `man/fs_family_report.Rd` (and any other `man/*.Rd` roxygen regenerates — list them); `tests/testthat/helper-synthetic-dgm.R`, `tests/testthat/test-fs-family-report.R`, new `tests/testthat/test-vi-deprecation.R`; `vignettes/forestsearch.Rmd`; `NEWS.md`; this document; the battery and report under `dev/reports/`. **Nothing else.** `NAMESPACE` must come out of `devtools::document()` byte-identical (§3 item 6). No `DESCRIPTION` version change (Larry sets the release version). No edits to any application document, driver, template or payload — those are audited and listed, not touched (§6).
- **Installs:** one `devtools::document()` and one `devtools::install()` from the branch working tree after the edits. The Stage-C battery runs against the installed package (parallel workers see only the installed package). The installed package is thereafter the branch build (`0.4.0.9000`); the report says so.
- **Check policy:** no `R CMD check`, no full test suite, no vignette build. CC runs the two touched test files, the new one, and six named dependents of the touched helper (non-gating), each under a hard cap. The full suite is Larry's pre-release gate.
- **Compute:** three fixtures × two sequential `forestsearch()` calls per stage, one `nb = 20` sequential survival bootstrap per stage, nine test files under caps. Hard caps in §2, §4.
- **Unattended.** Gates stop the task; never ask, never work around. Provenance is dirt-tolerant.

---

## 1. Provenance, branch, first commit — GATE

1. `git rev-parse --abbrev-ref HEAD` must be `master`; record HEAD; `git status --porcelain` verbatim (pre-existing untracked dirt is tolerated and never staged). Installed `packageVersion("forestsearch")`: record it with `packageDescription("forestsearch")[c("Version", "Packaged", "Repository", "Built")]` — expected the CRAN 0.4.0 build. If the installed package is *not* the CRAN build, Stage A still runs against it, and the report states exactly what the Stage-A baseline was built from.
2. `git checkout -b feature/vi-deprecation-lean-fits` (from `master` HEAD). **GATE:** the branch is created and checked out; record the base SHA — every post-condition in §7 diffs against it.
3. Copy this document to `dev/tasks/cc_task_vi_deprecation_2026-10-26.md`, `git add` by path, commit alone. Report both filenames.
4. `dev/reports/` is the record directory (the thresholds workstream's `STATUS_thresholds_workstream_2026-09-18.md` lives there). If it does not exist, use `dev/` and say so.

## 2. Stage A — baselines before any edit (installed package = CRAN 0.4.0)

Transplant the three fixtures from the committed Stage-0b bundle — `dev/glm-continuous-sims/stage0b_glm_2026-10-26/stage0b_common.R` (`fixture_survival()`), `fixture_continuous.R`, `fixture_binary.R`, exactly as committed (the application argument blocks; the binary one carries `vi.grf.min = -0.2`). Script: `dev/reports/vi_deprecation_2026-10-26/stageA.R`, sourcing a shared `vi_common.R` beside it that holds the fixtures, the save/compare helpers and this forest-fit counter:

```r
# Counts grf forest fits made anywhere inside expr, via base::trace() on the
# grf namespace (forestsearch calls grf::causal_forest / grf::causal_survival_forest
# through `::`, which resolves in that namespace).  Session-only; untraced on exit.
count_grf_fits <- function(expr) {
  cnt <- new.env(); cnt$n <- 0L
  # The tracer is evaluated in the traced function's frame, so it must carry
  # the counter environment inside the call (a closure named here is not found).
  tracer <- bquote(assign("n", get("n", envir = .(cnt)) + 1L, envir = .(cnt)))
  fns <- c("causal_forest", "causal_survival_forest")
  suppressMessages(for (f in fns) trace(f, tracer = tracer, where = asNamespace("grf"), print = FALSE))
  on.exit(suppressMessages(for (f in fns) untrace(f, where = asNamespace("grf"))), add = TRUE)
  value <- force(expr)
  list(value = value, n_fits = cnt$n)
}
# Verified in chat on grf 2.x: counts 2 for two direct grf::causal_forest() calls,
# 0 for none, 1 through forestsearch:::fit_causal_forest_glm(); untraced after.
# Runs one forestsearch() call, collecting warnings without suppressing the result.
run_collect <- function(args) {
  w <- character(0)
  r <- count_grf_fits(withCallingHandlers(
    suppressMessages(do.call(forestsearch::forestsearch, args)),
    warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") }))
  list(fs = r$value, n_fits = r$n_fits, warnings = w)
}
```

For each fixture, sequentially:

- **(i)** the args as committed;
- **(ii)** the args with `vi.grf.min` set to `NULL` and `max_n_confounders` absent (only the binary fixture differs between (i) and (ii); still run both for all three — the duplicate pair is a free determinism check).

Save per run, as one `.rds` beside the script (`stageA_<fixture>_<i|ii>.rds`): `sg.harm`; the membership (`df.est[[id.name]]` and `df.est$treat.recommend`, sorted by id); `grp.consistency$out_sg$result`; `find.grps$out.found$hr.subgroups`; the candidate-flow counts (`find.grps$filter_counts`, every `grp.consistency$n_candidates_*` element); `n_fits`; the warning count and texts; and `args_call_all` (stored, **excluded** from identity — it echoes the manipulated arguments). Save `sessionInfo()` once.

Survival only: a sequential `nb_boots = 20` bootstrap on the (ii) fit — `forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = 20L, seed = 8316951L, parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE))`, then `future::plan("sequential")`; save `results` minus `tmins_search`/`tmins_iteration` and minus any `timing` attribute, plus `SG_CIs`, `H_estimates`, `Hc_estimates`, `Ystar_mat`, `FSsg_tab` (whichever of these exist on the object — record the names found).

Print a table: fixture × run × (`n_fits`, warnings, `sg.harm`). Expected `n_fits`: 1 for every run whose args carry a non-NULL `vi.grf.min` (binary (i)), 0 otherwise. Commit the script, `vi_common.R` and the outputs by path ("Stage A baselines, vi deprecation"). **GATE:** binary (i) recorded exactly one forest fit and every other run zero; otherwise the premise is wrong — stop and report. Hard cap: 20 min.

## 3. The edits — transplant form, named 0.4.0 lines

Roxygen in this package is **plain Rd markup** (`\code{}`, `\strong{}`, `\link{}`, `\emph{}`), not markdown: use those, and escape a literal percent sign as `\%`.

**`R/forestsearch_main.R`**

1. Formals: leave L1536 (`max_n_confounders = 1000,`) and L1586 (`vi.grf.min = NULL,`) exactly as they are.
2. Deprecation guard: immediately after L1622 (`consistency_method  <- match.arg(consistency_method)`), insert:
   ```r
   # Deprecated arguments (0.4.0.9000): accepted so existing calls do not
   # break, warned about once per call, and ignored.  GRF variable-importance
   # screening/ordering of the candidate cut columns no longer exists: the
   # candidate family is the cut grid alone unless use_lasso, use_grf or
   # use_dina asks for a model-based front end.  Both formals are removed in
   # a later release.  The values are deliberately NOT reset, so that
   # args_call_all records what the caller passed.
   if (!is.null(vi.grf.min)) {
     warning("forestsearch(): 'vi.grf.min' is deprecated and ignored. GRF ",
             "variable-importance screening of the candidate cut columns was ",
             "removed; the candidate family is the cut grid alone unless ",
             "use_lasso, use_grf or use_dina is TRUE. The argument will be ",
             "removed in a future release -- drop it from the call.",
             call. = FALSE)
   }
   if (!missing(max_n_confounders) &&
       !identical(as.numeric(max_n_confounders), 1000)) {
     warning("forestsearch(): 'max_n_confounders' is deprecated and ignored ",
             "(it capped the variable-importance-ordered column list, which ",
             "no longer exists). It will be removed in a future release -- ",
             "drop it from the call.", call. = FALSE)
   }
   ```
3. Delete the inert-cap warning block, L3213–L3225 (from the comment `# The max_n_confounders truncation is applied only inside the variable-` through the closing `}` of the `if (is.null(vi.grf.min) && is.finite(max_n_confounders) && ...)` warning).
4. Delete the whole variable-importance block, L3227–L3313 (the `# SECTION 5: GRF VARIABLE IMPORTANCE SCREENING` banner, the `if (!is.null(vi.grf.min)) {` body and its `} else { conf.screen <- FSconfounders.name }` tail), and put in its place:
   ```r
   # Candidate cut columns are searched as get_FSdata() built them.  The GRF
   # variable-importance screening/ordering that stood here before
   # 0.4.0.9000 is deprecated (see vi.grf.min): no forest is fitted in the
   # forest-search path unless use_lasso, use_grf or use_dina is TRUE.
   conf.screen <- FSconfounders.name
   ```
   Nothing assigned inside the deleted block (`X`, `Y_vi`, `cs.forest`, `vi.cs`, `vi.cs2`, `vi.order`, `vi_max`, `vi_ratio`, `selected.vars`, `tau.rmst`, `tune_arg`, `y_evt_c`, `y_evt_t`) is used after L3313 — only `conf.screen`, consumed at L3325 (`df.confounders <- df[, conf.screen, drop = FALSE]`), which stays. Confirm by grep after the edit.
5. Roxygen (Rd markup):
   - `@param vi.grf.min` (L918–L937) and `@param max_n_confounders` (L507–L513): replace each with — `\strong{Deprecated} (0.4.0.9000): accepted and ignored with a warning. GRF variable-importance screening of the candidate cut columns was removed so that the candidate family is the cut grid alone unless \code{use_lasso}, \code{use_grf} or \code{use_dina} is \code{TRUE}. Will be removed in a future release; drop it from the call.` (For `max_n_confounders` add: `It capped the importance-ordered column list, which no longer exists.`)
   - `@param use_grf`, L387–L392: replace the sentence `This flag governs candidate-cut generation only: GRF variable-importance ordering of the cut columns is a separate step controlled solely by \code{vi.grf.min} and runs, or is skipped, regardless of \code{use_grf}.` with `This is the only switch that fits a GRF in the forest-search path (the former variable-importance screening, \code{vi.grf.min}, is deprecated and inert).`
   - `@param frac.tau`, L557–L559: `Applies to both the GRF variable-importance pre-screen inside \code{forestsearch()} and the downstream call to \code{\link{grf.subg.harm.survival}}.` → `Used by the downstream call to \code{\link{grf.subg.harm.survival}}.`
   - Then `grep -n -i "variable.importance\|section 5\|vi_ratio\|vi_max" R/forestsearch_main.R` and update any remaining roxygen or comment sentence that describes the removed step (L2051's comment `# Survival: not used (causal_survival_forest handles sign correctly)` concerns the GRF-cuts path, not the VI block — read it in context and leave it if so).
6. L1501 `#' @importFrom grf causal_survival_forest variable_importance` **stays as it is** (both symbols are imported by other files too; leaving it keeps `NAMESPACE` byte-identical — `git diff --name-only` must not list `NAMESPACE` after `devtools::document()`; if it does, stop and report the diff).
7. `fit_causal_forest_glm()` (`R/grf_helpers.R` L117) was called only from the deleted block and becomes unreferenced. **Leave it** (out of scope; not an `R CMD check` finding); list it in the report under OPEN ITEMS for Larry's disposition.

**`R/fs_family_report.R`**

8. L268–L283 (`vi_on <- …` through the `add("confounder cap", …)` call): replace the two rows with one:
   ```r
   add("GRF variable-importance screening (deprecated)", c("vi.grf.min", "max_n_confounders"),
       sprintf("vi.grf.min = %s; max_n_confounders = %s", fmt(vi.grf.min), fmt(max_n_confounders)),
       INE,
       "deprecated and ignored since 0.4.0.9000: no forest is fitted in the forest-search path unless use_lasso, use_grf or use_dina is TRUE")
   ```
   Keep L144–L145 (the `.arg()` reads) so the echo works on old and new objects; `n_cut_columns` is used by later rows — keep it. Grep that `vi_on` and `cap_status` are referenced nowhere else.
9. Roxygen: L14 `\code{use_dina = FALSE} and \code{vi.grf.min = NULL} may reasonably believe` → `\code{use_dina = FALSE} may reasonably believe`; L38–L40 `\code{max_n_confounders} is applied only inside the GRF variable-importance block, so it is inert when \code{vi.grf.min} is \code{NULL}; and` → `\code{vi.grf.min} and \code{max_n_confounders} are deprecated and always inert; and`; examples L74 and L87: drop `vi.grf.min = NULL` (and the comma before it) from both calls.

**`R/run_simulation_analysis.R`** — 10. delete L68 (`vi.grf.min               = -0.2,`) from the default argument list.

**`R/get_fsdata.R`** — 11. L598: `"(4) GRF-related thresholds (dmin.grf, frac.tau, vi.grf.min) "` → `"(4) GRF-related thresholds (dmin.grf, frac.tau) "`.

**Tests**

12. `tests/testthat/helper-synthetic-dgm.R`: delete L149 (`vi.grf.min             = -0.2,`) from `.fs_args_for()`. Eighteen test files use that helper; six of them are run in §4 (non-gating), the rest are Larry's pre-release suite.
13. `tests/testthat/test-fs-family-report.R` — exactly these edits, recorded verbatim in the report:
    - L83 title: `"resolution agreement (drift guard): stop_threshold and the max_n_confounders cap match the engine"` → `"resolution agreement (drift guard): stop_threshold matches the engine; the deprecated VI row is inert"`.
    - L104–L107: replace the comment and the two `expect_identical(... "confounder cap" ...)` / `(... "GRF variable-importance ordering" ...)` lines with `expect_null(fit$args_call_all$vi.grf.min)` (keep) and `expect_identical(r[r$stage == "GRF variable-importance screening (deprecated)", "status"], "inert")`.
    - L121–L126: replace the six lines with
      ```r
      # deprecated VI row: inert whatever is passed
      expect_identical(row_status(fs_family_report(.fr_driver_args(vi.grf.min = NULL)), "GRF variable-importance screening (deprecated)"), "inert")
      expect_identical(row_status(fs_family_report(.fr_driver_args(vi.grf.min = -0.2)), "GRF variable-importance screening (deprecated)"), "inert")
      expect_identical(row_status(fs_family_report(.fr_driver_args(max_n_confounders = 5)), "GRF variable-importance screening (deprecated)"), "inert")
      ```
    - Leave L17 (the driver mirror keeps `vi.grf.min = -0.2`: the drivers are not edited in this task), L147–L148, L159 and L177 (that `forestsearch()` call is already inside `suppressWarnings()` and now also exercises that a deprecated value round-trips through `args_call_all`).
14. New `tests/testthat/test-vi-deprecation.R`. Fixture: `.make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)` with the `forestsearch()` call of `test-fs-family-report.R` L88–L95 (`fs.splits = 20, n.min = 30, maxk = 2, use_lasso = FALSE, use_grf = FALSE, sg_focus = "eff", quiet = TRUE, seedit = 1L, parallel_args = list(plan = "sequential")`), with `vi.grf.min` removed from that call and the whole call wrapped in a local `.run <- function(...)` that takes extra arguments. Keep every `forestsearch()` call well under 30 s. Tests:
    - deprecation warnings: `expect_warning(.run(vi.grf.min = -0.2), "deprecated")`; `expect_warning(.run(max_n_confounders = 5), "deprecated")`; a default `.run()` whose warnings are collected with `withCallingHandlers` (the fixture may warn for unrelated reasons) has none matching `"deprecated"`;
    - inert: with warnings suppressed, `.run(vi.grf.min = -0.2)` and `.run(max_n_confounders = 5)` each give `identical()` `sg.harm`, sorted membership, `grp.consistency$out_sg$result` and `find.grps$out.found$hr.subgroups` to `.run()`;
    - no forest: with the `count_grf_fits()` helper from §2 copied into the test file (session-only `trace()`, `untrace()` on exit), `.run(vi.grf.min = -0.2)` makes zero grf fits;
    - echo: `fs_family_report(.run(vi.grf.min = -0.2))` has the `"GRF variable-importance screening (deprecated)"` row at `"inert"` and its `values` contain `"vi.grf.min = -0.2"`.

**Vignette and NEWS**

15. `vignettes/forestsearch.Rmd`: delete L319 (`max_n_confounders          = 1000, …`) and L324 (`vi.grf.min                 = -0.2, …`). No build. `inst/doc/` holds generated copies (`forestsearch.R`, `.Rmd`, `.html`): if it is tracked in the repo, do not edit it (regenerated at release) and say so in the report; if untracked, say that too.
16. `NEWS.md`: the CRAN file starts at `# forestsearch 0.4.0`. If `master` has no development heading above it, add `# forestsearch 0.4.0.9000 (development)`. Under it, a `### Deprecated` entry: both arguments accepted and ignored with a warning; the variable-importance screening block removed; rationale (fixed candidate family: no model is fitted in the forest-search path unless `use_lasso`/`use_grf`/`use_dina` is `TRUE`); effect on callers who passed a non-NULL value (membership unchanged at values `<= 0`; label clause order may differ; a binding value no longer screens); `run_simulation_analysis()`'s default list, `.fs_args_for()` and the vignette no longer pass them; `fs_family_report()` shows one `inert` row for both; formals to be removed in a later release.

Then `devtools::document()` (commit the regenerated `man/*.Rd`, listed by name; `NAMESPACE` must be unchanged) and `devtools::install()` from the working tree; record `packageDescription("forestsearch")[c("Version", "Built")]` afterwards. Commit all edits of this section as **one** commit by explicit path ("Deprecate vi.grf.min and max_n_confounders; remove VI screening"). Stage-set assertion before the commit: `git diff --cached --name-only` ⊆ the allowed list in §7.

## 4. Stage C — verification against the installed branch build

Script `dev/reports/vi_deprecation_2026-10-26/stageC.R`, same `vi_common.R`, same fixtures, same runs and saved components as §2 (`stageC_<fixture>_<i|ii>.rds`). Assertions — each a gate; any failure stops with the component named:

- **survival (i), (ii) and continuous (i), (ii):** every saved component `identical()` to Stage A — bit-identity for default callers; zero deprecation warnings; zero forest fits.
- **binary (ii):** zero warnings; zero fits; `identical()` to Stage A (ii), component by component.
- **binary (i)** (`vi.grf.min = -0.2`, as the application document has it): exactly one warning matching `"deprecated"`; zero forest fits (Stage A had one); every saved component `identical()` to Stage C (ii) and to Stage A (ii). Against Stage A (i) — the old screened run — **report, not gate:** `setequal()` of the `sg.harm` clauses; `identical()` of the sorted membership; equality of the candidate table after sorting by (effect, K, label); whether the clause order of `sg.harm` differs. Expected: membership identical; order may differ.
- **survival `nb = 20` bootstrap payload:** `identical()` to Stage A, component by component (volatile fields excluded as in §2).
- **Tests, gating:** `testthat::test_file()` on `tests/testthat/test-vi-deprecation.R` and `tests/testthat/test-fs-family-report.R`, each under `timeout 10m`; all pass, no skips introduced.
- **Dependents of the edited helper, non-gating** (report pass/fail per file; do not fix anything here): `test-cross-outcome-parity.R`, `test-return-shape-contract.R`, `test-sg-focus-transparency.R`, `test-get-fsdata-diagnostics.R`, `test-dfbeta-lm-influence-path.R`, `test-downstream-consumers.R` — each under `timeout 10m`, total cap 45 min (stop the loop at the cap and list what did not run). Any failure is an OPEN ITEM with its verbatim message.

`stageC.R` writes `stageC_summary.txt`: one line per gate and per reported finding (`PASS`/`FAIL`/`FINDING: …`). Hard cap: 20 min for the identity runs; test caps as above.

## 5. Report — `dev/reports/REPORT_vi_deprecation_2026-10-26.md`

1. Provenance: base SHA on `master`, branch, every commit SHA in order, installed build before (§1) and after (§3), note that the installed package is now the branch build.
2. The edits, file by file, with the final `git diff --stat` against the base SHA and the regenerated `man/` files listed; the verbatim test-file edits of §3 item 13.
3. Stage A vs Stage C: one table per fixture — component × (identical / differs), warning counts, forest-fit counts; the binary (i)-vs-old-(i) findings.
4. Test results, verbatim summary lines, gating and non-gating separately.
5. §6 audit list.
6. Script deviations (none expected; verbatim if any). OPEN ITEMS (at least: `fit_causal_forest_glm()` now unreferenced; the §6 entries). Ten-line verdict.
7. Closing block: the §7 post-condition commands and their output.

## 6. Audit — read-only, listed for Larry's disposition, no edits

`grep -rn "vi\.grf\.min\|max_n_confounders"` over `quarto/`, `dev/`, `inst/`, `R/`, `tests/`, `vignettes/`, `NEWS.md` in this repo, and over `~/GitHub/fs-glms-interpretable/` if it exists (read-only). Tabulate every explicit setting outside this task's edits: file, line, value, and whether the file is a committed application document (with rendered HTML as its record), a driver, a template, or a generated copy (`inst/doc/*`). Expected entries include `quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd` (`vi.grf.min = -0.2`) and the `.fr_driver_args()` mirror in `test-fs-family-report.R` L17. **No edits to any of them in this task:** after this change they emit one warning and behave as `NULL`; whether to clean their source lines (and re-render their records) is Larry's call.

## 7. Commit and post-conditions — GATE

Allowed paths for this task's commits: `dev/tasks/cc_task_vi_deprecation_2026-10-26.md`; `R/forestsearch_main.R`; `R/fs_family_report.R`; `R/run_simulation_analysis.R`; `R/get_fsdata.R`; `man/*.Rd` regenerated by roxygen; `tests/testthat/helper-synthetic-dgm.R`; `tests/testthat/test-fs-family-report.R`; `tests/testthat/test-vi-deprecation.R`; `vignettes/forestsearch.Rmd`; `NEWS.md`; `dev/reports/REPORT_vi_deprecation_2026-10-26.md` and `dev/reports/vi_deprecation_2026-10-26/`. Before every commit `git diff --cached --name-only` must be a subset of that list. After the final commit print in the closing block:

```
git diff --name-only <base SHA> HEAD | grep -v -E '^(dev/|man/.*\.Rd$|NEWS\.md$|vignettes/forestsearch\.Rmd$|R/(forestsearch_main|fs_family_report|run_simulation_analysis|get_fsdata)\.R$|tests/testthat/(helper-synthetic-dgm|test-fs-family-report|test-vi-deprecation)\.R$)'
    # must print nothing (in particular: no NAMESPACE, no DESCRIPTION)
grep -n "variable_importance\|causal_survival_forest\|causal_forest\|fit_causal_forest_glm\|vi_ratio\|vi_max\|vi\.cs\b\|vi\.order" R/forestsearch_main.R | grep -v -E "^[0-9]+:[[:space:]]*#"
    # must print nothing: no code line (comments and roxygen excluded) touches the removed step
grep -n "vi\.grf\.min\|max_n_confounders" R/forestsearch_main.R | grep -v -E "^[0-9]+:[[:space:]]*#"
    # must print exactly: the two formals and the guard's own lines (the two `if (` tests and the strings inside the warning() calls) -- nothing else
grep -rn "vi\.grf\.min\|max_n_confounders" R/ | grep -v -E "^R/(forestsearch_main|fs_family_report)\.R:"
    # must print nothing
git status --porcelain | grep -v '^?? '
    # must print nothing
```

If any command prints something unexpected, stop and report what it prints.

Finally, copy `dev/reports/REPORT_vi_deprecation_2026-10-26.md` and `dev/reports/vi_deprecation_2026-10-26/stageC_summary.txt` to `~/Downloads` and print `ls -la` of both copies.

## 8. Out of scope

- No change to the GRF or DINA identifier paths (`use_grf`, `use_dina`, `subgroup_method = "grf"/"dina"`, `grf.subg.harm.*`, `tune_grf`, `.fs_grf_num_threads`, `fit_causal_forest_glm`).
- No edit to any application document, driver, campaign template, generated `inst/doc/` copy or payload (audit only, §6). No re-render.
- No `DESCRIPTION` version bump; no `NAMESPACE` change; no `R CMD check`; no full test suite; no vignette build; no push; no merge into `master`.
- The lean-fit changes (A and F) are the next task on this branch, not this one.
