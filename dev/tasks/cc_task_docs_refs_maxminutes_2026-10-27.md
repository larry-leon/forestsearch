# CC task — reference update, `max.minutes` deprecation, untested-option notes, manuscript-aligned examples

**Date:** 2026-10-27. **Repo:** forestsearch, branch `master` at `b574c618` (verify with `git rev-parse HEAD` before starting; stop if it differs and report the SHA).
**Authorised by:** Larry, chat of 2026-10-27. **Release:** December batch; nothing here is a CRAN release.
**Task doc path:** copy this file from `~/Downloads/` to `dev/tasks/cc_task_docs_refs_maxminutes_2026-10-27.md` and commit it first and alone.

## 0. Rules in force

- No `R CMD check`, no vignette build, no full test suite inside this task. Larry runs those at release.
- Never edit a legacy analysis document or payload (`quarto/applications/**`, `_payloads/`). `quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd` L≈584 passes `max.minutes = 3`; it stays as is and is recorded in the report (§5).
- Vignettes are a separate workstream: do not touch `vignettes/` except to remove `max.minutes` call-site arguments (§1.4).
- Gates are stop-on-failure. On a failed gate: stop, write the report with the failure, do not work around it.
- Read every file from the working tree before editing; line numbers below are from `b574c618` and are hints, not addresses.
- After `devtools::document()`, restore `Config/roxygen2/version: 8.1.0` in `DESCRIPTION` (roxygen 8.1.1 re-stamps it); the only `DESCRIPTION` lines that may change are the two named in §1.1 and §1.5.
- CC never pushes, fetches, pulls or merges.

## 1. Changes

All edits are to roxygen comments, `NEWS.md`, `DESCRIPTION`, one test file and one forwarding line, plus one new warning. Classification for Larry's record: **(a), (c), (d)** are documentation only; **(b)** adds a warning (behaviour change) and removes one forwarded argument (moves existing code); **(e)** is a `DESCRIPTION` Suggests line. No method change anywhere.

### 1.1 (a) Primary methodology references — both León & Anderson (2026) papers

Reference strings (use verbatim):

- **LA-1 = 2026a:** León, L. F. and Anderson, K. M. (2026a). Inference for standard trial estimands after data-driven subgroup discovery. arXiv:2609.38361 (posted 2026-09-29). `\doi{10.48550/arXiv.2609.38361}`
- **LA-2 = 2026b:** León, L. F. and Anderson, K. M. (2026b). Interpretable subgroup identification across survival and generalized-linear outcomes. arXiv:2610.28379 (posted 2026-10-20). `\doi{10.48550/arXiv.2610.28379}`

Suffixes are by posting date (Larry's decision, 2026-10-27), so LA-1 is always listed before LA-2.

Edits:

1. `R/forestsearch_main.R`, `forestsearch()` `@description` (L≈318–334): after the sentence ending `Leon et al. (2024) \doi{10.1002/sim.10163}.` add:
   ```
   #' Post-selection inference for the reported subgroup (multiplier
   #' resampling, \code{mr_inference = TRUE}) follows León and Anderson
   #' (2026a) \doi{10.48550/arXiv.2609.38361}. The current methodology for
   #' the search itself -- forest search as closed-form stability selection
   #' on splitting consistency, the family-wise false-declaration rate and
   #' the two-subgroup claim -- is given in León and Anderson (2026b)
   #' \doi{10.48550/arXiv.2610.28379}.
   ```
2. `R/forestsearch_main.R`, `forestsearch()` `@references` (L≈1460–1467): make LA-1 (2026a) then LA-2 (2026b) the **first two** `\item`s, each on its own lines in the format `\item León, L. F. and Anderson, K. M. (2026a). Inference ... arXiv:2609.38361. \doi{...}`; keep the three existing items after them; add `\item León et al. (2024) ... \doi{10.1002/sim.10163}` as the fourth item if it is not already there.
3. `DESCRIPTION` `Description:` field: after `Methods are described in Leon et al. (2024) <doi:10.1002/sim.10163>.` and before the existing arXiv sentence (which cites 2609.38361 only — read the field; the handoff places it at L22), leave the existing LA-1 sentence where it is and append after it `The search methodology is given in León and Anderson (2026) <doi:10.48550/arXiv.2610.28379>.` so that the field cites 2024, LA-1, LA-2 in that order. Keep the field within CRAN's formatting (no line over 80 characters, 4-space continuation).
4. Every other `man/` page that cites 2609.38361 (`forestsearch-package`, `print.forestsearch`, `fs_mr_inference`, `fs_betaHhat_one` per the 2026-10-26 handoff — confirm with `grep -rn "2609.38361" R/`): add LA-2 beside LA-1 in the same style as that page uses, **only** where the page's text discusses the search/screen or the method in general. A page that is strictly about MR/IJ inference keeps LA-1 alone. List the decision per file in the report.
5. `README.md` Citation section: add a second `@misc` BibTeX entry for each arXiv paper (`LeonAnderson2026a` = 2609.38361, `LeonAnderson2026b` = 2610.28379) after `Leon2024`; `eprint`, `archivePrefix = {arXiv}`, `primaryClass = {stat.ME}`, `year = {2026}`. Do not touch the rest of the README. (If README.md is generated from a README.Rmd, edit the .Rmd and note it; do not knit.)

### 1.2 (b) `max.minutes` deprecation — both exported formals

Standard: the `vi.grf.min` pattern (formal kept at its position and default; one `warning(call. = FALSE)`; no effect; removal the release after December's).

1. `R/forestsearch_main.R`, `forestsearch()`:
   - Replace the `@param max.minutes` text (L≈884–891) with:
     ```
     #' @param max.minutes \strong{Deprecated} (0.4.0.9000): accepted and
     #'   ignored. It was documented as a wall-clock budget for the
     #'   combination search but was never enforced in the parallel search
     #'   path. Passing it emits one warning. Search scope is governed by
     #'   \code{maxk}, the candidate-factor front ends (\code{use_grf},
     #'   \code{use_lasso}, \code{conf_force}) and the number of workers.
     #'   The formal will be removed in a later release; drop it from calls.
     ```
   - In the deprecated-arguments block (L≈1607–1628), after the `max_n_confounders` warning add:
     ```r
     if (!missing(max.minutes)) {
       warning("forestsearch(): 'max.minutes' is deprecated and ignored. The ",
               "wall-clock budget was never enforced in the parallel search; ",
               "scope is governed by maxk and the candidate-factor controls. ",
               "The argument will be removed in a future release -- drop it ",
               "from the call.", call. = FALSE)
     }
     ```
     Trigger is `!missing()` (any explicitly passed value, default or not), so the warning never fires on a default call and fires exactly once when passed. Do **not** reset the value (`args_call_all` keeps recording it, as for `vi.grf.min`).
   - Remove the forwarding line `max.minutes = max.minutes,` from the internal `subgroup.search(...)` call (L≈3313). `subgroup.search()` then uses its own default; it is inert there, so nothing changes. This is required: otherwise every `forestsearch()` call would trigger `subgroup.search()`'s warning.
   - Update the comment at the top of the deprecated-arguments block to name the three formals.
2. `R/subgroup_search.R`, `subgroup.search()`:
   - Replace the `@param max.minutes` text (L≈39–47) with the same wording as above (substituting `\code{maxk} and the candidate-factor controls`).
   - At the top of the function body (after any `match.arg` lines) add the same `if (!missing(max.minutes)) warning(...)` with prefix `"subgroup.search(): "`.
   - Leave the threading to `search_combinations_parallel()` (L≈171, L≈264) in place: it is internal, never consults the value, and removing the parameter from an internal signature is out of scope.
3. `grep -rn "max.minutes" R/ tests/ vignettes/ inst/` and list every hit in the report. Expected: the two roxygen blocks, the two warnings, the internal threading, `vignettes/forestsearch.Rmd` (§1.4), tests (§3).
4. `NEWS.md`, "Deprecated" section: add one bullet in the style of the `vi.grf.min` bullet: `max.minutes` in `forestsearch()` and `subgroup.search()` deprecated; was documented as inert "scheduled for deprecation in v0.3.0" and never actioned; now warns once when passed; no effect; removed in a later release.

### 1.3 (c) Untested-option notes — `adjust_covariates`, `ps_method`, `ps_adjust_method`, `ps_hat`

Documentation only. In `R/forestsearch_main.R`:

1. Prepend to each of the four `@param` blocks (L≈1077–1120) this first sentence, verbatim:
   `This option is implemented but has not been tested against the manuscript analyses or the simulation studies of León and Anderson (2026a, 2026b), which use the unadjusted treatment-only fits.`
   Keep the rest of each block unchanged. Do not add any "Not validated" label or bold marker (Larry's decision, 2026-10-27).
2. Add one `@section Options not tested against the manuscript analyses:` to `forestsearch()` (place it immediately before `@examples`) with a single paragraph naming the four arguments, repeating that the published analyses and simulations use the unadjusted treatment-only fits, that no coverage or bias results exist for these paths, and that users who need them should check them on their own data (e.g. against a direct `coxph()`/`glm()` fit of the selected subgroup). No claim that they are wrong; only that they are untested.
3. If `subgroup.search()` or any other exported function documents `adjust_covariates` or `ps_*` (grep `@param adjust_covariates` and `@param ps_` across `R/`), add the same first sentence there and list the files in the report.

### 1.4 Remove `max.minutes` from every package call site (decided: Larry, 2026-10-27)

From the §1.2 item 3 grep, delete the `max.minutes = ...` argument (and its trailing comment) from every call in `vignettes/` (known: `vignettes/forestsearch.Rmd` L≈342), `tests/testthat/`, `inst/` and any `R/` wrapper that forwards it (e.g. `run_subgroup_sims.R` / `run_simulation_analysis()` defaults, if present). Numerically nothing changes (the argument is inert). This is the only vignette edit permitted in this task: no other vignette content is touched. `quarto/applications/**` is **not** edited (legacy rule); its hits are recorded in the report only.

### 1.5 (d) `@examples` rewritten to the manuscript analyses

Replace the entire `@examples` block of `forestsearch()` (L≈1393–1456; currently five examples on a non-existent `trial_data`, two of which use argument names that need checking) with three examples mirroring the manuscript analyses. Keep everything inside `\dontrun{}` (each runs minutes with resampling). Argument values are taken from the committed analysis documents and were checked (chat, 2026-10-27) against the lean-build re-renders' resolved-parameter tables; CC verifies each against the current formals (`names(formals(forestsearch))`) and stops if any name does not exist.

```r
#' @examples
#' \dontrun{
#' ## Example 1 -- Survival (GBSG, hormone therapy): the analysis of
#' ## León and Anderson (2026a, Sec. 5) and the survival-post-selection
#' ## vignette. Forest search over the cut grid (no model-based front ends),
#' ## effMaxSG selection in a 10% effect neighbourhood, multiplier-resampling
#' ## inference for the reported subgroup computed inside the fit.
#' library(survival)
#' df_gbsg <- within(gbsg, {
#'   id          <- seq_len(nrow(gbsg))
#'   time_months <- rfstime / 30.4375
#'   grade3      <- ifelse(grade == "3", 1, 0)
#' })
#' fs_surv <- forestsearch(
#'   df.analysis            = df_gbsg,
#'   outcome.name           = "time_months",
#'   event.name             = "status",
#'   treat.name             = "hormon",
#'   id.name                = "id",
#'   confounders.name       = c("age", "meno", "size", "grade3", "nodes",
#'                              "pgr", "er"),
#'   is.RCT                 = TRUE,
#'   use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
#'   conf.cont_jcuts        = list(er = 10, pgr = 10),
#'   conf_force             = c("er <= 0", "pgr <= 0"),
#'   subgroup_method        = "consistency",
#'   sg_focus               = "effMaxSG",
#'   selection_rule         = "neighborhood",
#'   effect_neighborhood    = 0.10,
#'   hr.threshold           = 1.0,
#'   hr.consistency         = 1.0,
#'   pconsistency.threshold = 0.90,
#'   consistency_method     = "resample",
#'   fs.splits              = 1000,
#'   n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
#'   mr_inference           = TRUE,
#'   mr_inference_args      = list(draws = 5000L),
#'   parallel_args          = list(plan = "multisession", workers = 4),
#'   seedit                 = 8316951
#' )
#' print(fs_surv)
#' mr_estimates_table(fs_surv)
#'
#' ## Example 2 -- Binary outcome, odds ratio (ACTG175, ZDV+ddI vs ddI):
#' ## the binary application of León and Anderson (2026b). Adverse-scale
#' ## outcome (1 = no CD4 improvement at week 20), OR thresholds at the null,
#' ## same selection rule and multiplier-resampling gate as Example 1.
#' if (requireNamespace("speff2trial", quietly = TRUE)) {
#'   actg <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
#'   actg <- actg[!is.na(actg$cd420), ]
#'   actg$id           <- seq_len(nrow(actg))
#'   actg$treat        <- as.integer(actg$arms == 1L)
#'   actg$y_binary_neg <- as.integer(!(actg$cd420 > actg$cd40))
#'   bin_vars <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
#'   for (v in bin_vars) actg[[v]] <- as.numeric(actg[[v]])
#'   fs_bin <- forestsearch(
#'     df.analysis            = actg,
#'     outcome.name           = "y_binary_neg",
#'     event.name             = "y_binary_neg",
#'     treat.name             = "treat",
#'     id.name                = "id",
#'     confounders.name       = c("age", "preanti", "wtkg", "karnof", "cd40",
#'                                "cd80", bin_vars),
#'     outcome_type           = "binary",
#'     effect_measure         = "OR",
#'     adverse_outcome        = TRUE,
#'     use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
#'     conf.cont_jcuts        = list(cd40 = 10, wtkg = 10),
#'     subgroup_method        = "consistency",
#'     sg_focus               = "effMaxSG",
#'     selection_rule         = "neighborhood",
#'     effect_neighborhood    = 0.10,
#'     effect.threshold       = 1.0,
#'     consistency.threshold  = 1.0,
#'     pconsistency.threshold = 0.90,
#'     consistency_method     = "resample",
#'     fs.splits              = 500,
#'     n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
#'     mr_inference           = TRUE,
#'     mr_inference_args      = list(draws = 5000L),
#'     parallel_args          = list(plan = "multisession", workers = 4),
#'     seedit                 = 8316951
#'   )
#'   print(fs_bin)
#' }
#'
#' ## Example 3 -- Continuous outcome, mean difference (ACTG175): the
#' ## continuous application of León and Anderson (2026b), with the field
#' ## (calibrated one-sided) bounds for the selected subgroup and its
#' ## complement. Outcome is CD4 decline from baseline to week 20 (adverse
#' ## scale); thresholds are on the mean-difference scale.
#' if (requireNamespace("speff2trial", quietly = TRUE)) {
#'   actg <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
#'   actg <- actg[!is.na(actg$cd420), ]
#'   actg$id        <- seq_len(nrow(actg))
#'   actg$treat     <- as.integer(actg$arms == 1L)
#'   actg$y_decline <- actg$cd40 - actg$cd420
#'   bin_vars <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
#'   for (v in bin_vars) actg[[v]] <- as.factor(actg[[v]])
#'   fs_cont <- forestsearch(
#'     df.analysis            = actg,
#'     outcome.name           = "y_decline",
#'     treat.name             = "treat",
#'     id.name                = "id",
#'     confounders.name       = c("age", "preanti", "wtkg", "karnof", "cd40",
#'                                "cd80", bin_vars),
#'     outcome_type           = "continuous",
#'     effect_measure         = "MD",
#'     adverse_outcome        = TRUE,
#'     use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
#'     conf.cont_jcuts        = list(age = 10, preanti = 10, wtkg = 10,
#'                                   karnof = 10, cd40 = 10, cd80 = 10),
#'     sg_focus               = "effMaxSG",
#'     selection_rule         = "neighborhood",
#'     effect_neighborhood    = 0.20,
#'     effect.threshold       = 10,
#'     consistency.threshold  = 5,
#'     pconsistency.threshold = 0.90,
#'     consistency_method     = "resample",
#'     fs.splits              = 500,
#'     n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
#'     mr_inference           = TRUE,
#'     mr_inference_args      = list(ci_method = "field", draws = 5000L,
#'                                   include_complement = TRUE,
#'                                   field_complement = TRUE,
#'                                   field_scale_complement = "selected",
#'                                   return_reselection = TRUE),
#'     parallel_args          = list(plan = "sequential"),
#'     seedit                 = 8316951
#'   )
#'   print(fs_cont)
#' }
#' }
```

Notes for CC:
- Example 1's `use_twostage`, `minp`, `stop_threshold` are left at the defaults (the manuscript document pre-resolves them to the defaults under effMaxSG); confirm `use_twostage = TRUE`, `minp = 0.025` are the formals' defaults and note if not.
- Example 2 omits `d1.min`/`d0.min`-irrelevant survival arguments; `event.name = outcome.name` is the binary convention used by the committed document.
- If any `mr_inference_args` name in Example 3 is not accepted by `fs_mr_inference()` (`grep -n "mr_inference_args" R/forestsearch_main.R` and the `fs_mr_inference()` formals), drop that element and record it; do not invent a replacement.

### 1.6 (e) `DESCRIPTION` Suggests

`speff2trial` is not in `Suggests` and ACTG175 is not bundled (CC confirms with `ls data/ inst/extdata/ 2>/dev/null`; if a bundled ACTG175 exists, use it in Examples 2–3 instead and skip this item). Add `speff2trial` to `Suggests:` (alphabetical placement is not required; append before `testthat`). This is the second permitted `DESCRIPTION` edit. Examples 2–3 are guarded by `requireNamespace()` as CRAN requires for Suggests.

## 2. Procedure

1. `git rev-parse HEAD` = `b574c618`; working tree clean (`git status --porcelain` empty). Stop otherwise.
2. Copy this task doc into `dev/tasks/` and commit it alone: `dev: task doc — refs, max.minutes deprecation, untested-option notes, examples (2026-10-27)`.
3. (reserved)
4. Record the baseline: `R -q -e 'tools::md5sum(list.files("man", full.names = TRUE))' > dev/reports/docs_refs_2026-10-27/man_md5_before.txt` (create the directory).
5. Apply §1.1–§1.6 in that order, one commit each (messages: `docs: ...`, `deprecate: max.minutes ...`, `docs: untested-option notes ...`, `docs: manuscript-aligned examples ...`, `DESCRIPTION: Suggests speff2trial`; the §1.4 removals go in the `deprecate:` commit).
6. `devtools::document()`; then restore `Config/roxygen2/version: 8.1.0`; `git diff --stat man/ DESCRIPTION NAMESPACE`. Gate: `NAMESPACE` unchanged; `DESCRIPTION` diff is exactly the two field changes (§1.1 item 3, §1.6) — no stamp change.
7. Write the test file §3, run it (gate).
8. Run gates §4. Commit `man/` + test + NEWS: `docs: regenerate man pages; test-max-minutes-deprecation`.
9. Write `dev/reports/REPORT_docs_refs_maxminutes_2026-10-27.md` (§5), commit `dev: report — refs, max.minutes deprecation, examples (2026-10-27)`, and copy the report and `man/forestsearch.Rd` to `~/Downloads/`.
10. Do not reinstall the package unless a gate requires it (§4.2 does); use `devtools::install(upgrade = FALSE)`.

## 3. Test — `tests/testthat/test-max-minutes-deprecation.R`

Pattern of the existing `vi.grf.min` test (find it: `grep -ln "vi.grf.min" tests/testthat/*.R`; reuse its fixture/helper and `skip_on_cran()` convention). Expectations:

1. Default `forestsearch()` call on the smallest fixture emits **no** warning matching `"max.minutes"` (`expect_no_warning` with the `message`/`class` filter the existing test uses, or `expect_silent`-equivalent restricted to that pattern).
2. The same call with `max.minutes = 5` emits **exactly one** warning matching `"'max.minutes' is deprecated"` (`expect_warning(..., regexp = "max.minutes", all = TRUE)` wrapped so that a second warning fails).
3. With a fixed `seedit`, the returned object with `max.minutes = 5` is `identical()` to the default-call object after removing `args_call_all$max.minutes` from both (and any timing fields the existing identity tests already strip — reuse that stripper). The handoff's saved-vs-saved `.rds` pattern applies if a `data.table` is inside: compare `readRDS()` copies.
4. `subgroup.search()` called directly with `max.minutes = 5` emits exactly one warning with prefix `subgroup.search():`; called via `forestsearch()` with default arguments it emits none (covered by 1, but assert the prefix is absent too).
5. `expect_true("max.minutes" %in% names(formals(forestsearch)))` and the same for `subgroup.search` (the formals are kept through this release).

Run: `NOT_CRAN=true Rscript -e 'testthat::test_file("tests/testthat/test-max-minutes-deprecation.R", package = "forestsearch", load_package = "source")'`. Gate: all pass.

## 4. Gates (all stop-on-failure)

- **4.1 Roxygen/Rd.** `R -q -e 'tools::checkRd("man/forestsearch.Rd"); tools::checkRd("man/subgroup.search.Rd")'` clean; `R -q -e 'parse(text = tools::Rd2ex("man/forestsearch.Rd", commentDontrun = FALSE))'` succeeds (the examples parse); every argument name used in the three examples is in `names(formals(forestsearch))` (script this; list the names checked in the report). The Rd Description of `forestsearch` contains both `2610.28379` and `2609.38361`; `R CMD Rd2txt man/forestsearch.Rd | grep -c "León and Anderson"` ≥ 3.
- **4.2 Behaviour.** `test-max-minutes-deprecation.R` green (§3). Also `NOT_CRAN=true` run of the existing `vi.grf.min` deprecation test and `tests/testthat/test-lean-fits.R` — both green (they exercise default calls; this confirms no stray warning changed their expectations). Install first with `devtools::install(upgrade = FALSE)` and use `load_package = "installed"` for the lean-fits file, as in Task 2.
- **4.3 Scope.** `git diff --stat b574c618..HEAD` touches only: `R/forestsearch_main.R`, `R/subgroup_search.R`, any `R/` file from §1.1 item 4 / §1.3 item 3 (listed), `man/*.Rd` regenerated from those, `DESCRIPTION` (two fields), `NEWS.md`, `README.md` (and `README.Rmd` if present), `tests/testthat/test-max-minutes-deprecation.R`, `vignettes/*.Rmd` and `tests/testthat/*.R` only for `max.minutes` removals (§1.4), `dev/tasks/`, `dev/reports/`. Nothing under `quarto/`, `data/`; `inst/` only for `max.minutes` removals (§1.4). `NAMESPACE` unchanged.
- **4.4 Legacy untouched.** `git status --porcelain quarto/` empty at the end.
- **4.5 Line length.** No new line in `R/` or `DESCRIPTION` exceeds 80 characters (`awk 'length > 80' R/forestsearch_main.R R/subgroup_search.R DESCRIPTION` over the changed hunks only; pre-existing long lines are not this task's).

## 5. Report — `dev/reports/REPORT_docs_refs_maxminutes_2026-10-27.md`

Verbatim, not summarised:
- the SHA ladder (task-doc commit → each change commit → man/test commit → report commit);
- §1.1 item 4: per-file decision (LA-2 added / LA-1 only) with the sentence that was added;
- the full `grep -rn "max.minutes"` listing before and after, including the legacy `quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd` hit left in place and every §1.4 removal with file and line;
- the example argument-name check (names verified; any `mr_inference_args` element dropped under §1.5 note 3);
- `git diff` of `DESCRIPTION` and `NEWS.md` in full;
- test output of §3 and §4.2 (pass counts, wall-clock);
- `man_md5_before.txt` vs after: the list of `.Rd` files whose hash changed (expected: `forestsearch.Rd`, `subgroup.search.Rd`, the §1.1 item 4 pages, nothing else);
- anything that surprised CC, with the file and line.

End of task.
