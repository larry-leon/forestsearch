# CC task — `run_null_calibration()` example reduced to an illustrative size (CRAN example-timing NOTE)

**File:** `dev/tasks/cc_task_example_timing_2026-10-29.md` · **Issued:** 2026-10-29 by chat, commissioned by Larry (prompted by the `--as-cran` check on `feature/robustness-0.4.1`: "Examples that took more than 5s: run_null_calibration 79.4 s elapsed, suggest_thresholds 10.2 s", `dev/reports/robustness_041_2026-10-28/check_as_cran_2026-10-29.log`).
**Repo:** `forestsearch` · **Start on `master`** with the robustness branch merged (gate below). CC creates **`feature/example-timing-041`** from it. CC never fetches, pulls, pushes or merges.
**Transport:** `~/Downloads` → `dev/tasks/` (this document, committed first and alone).

## ⚠ CATEGORY — roxygen example text only

- The only `R/` edit is inside the `@examples` block of `run_null_calibration()` in `R/calibrate_null_correction.R` (0.4.0.9000 / master `606a78ea` lines 247–296; verify by content). No function body, formal, default or documentation paragraph changes. `man/run_null_calibration.Rd` regenerates from it (expected `man/` change — this is the one task where a `man/` diff is allowed, and only that file).
- **`suggest_thresholds()` is measured, not changed.** Its example is already small; the 10 s is in the function. The per-call times go in the report for Larry's decision in the release task.
- **Check policy:** no `R CMD check`, no test suite, no vignette build. Acceptance = the example runs and its elapsed time is recorded. The release task's `--as-cran` check is the gate.
- **Compute:** a handful of example runs, each well under a minute. Hard cap 20 min overall.
- **Unattended.** Gates stop the task; never ask.

## 1. Provenance — GATE

1. `git rev-parse --abbrev-ref HEAD` = `master`; `git status --porcelain` has no tracked modifications (untracked dirt tolerated, never staged). Record HEAD.
2. `git merge-base --is-ancestor 1f14b012 HEAD` must exit 0 (the robustness branch is merged). Otherwise stop and report.
3. `git checkout -b feature/example-timing-041`; copy this document to `dev/tasks/cc_task_example_timing_2026-10-29.md`; add by path; commit alone.
4. Record directory `dev/reports/example_timing_2026-10-29/`; report `dev/reports/REPORT_example_timing_2026-10-29.md`.
5. `packageDescription("forestsearch")$Built` — the installed build must be the current tree (the robustness build, 2026-10-29); if in doubt `devtools::install(upgrade = FALSE)` first and say so.

## 2. Baseline timing of the two examples (installed build, before the edit)

Script `dev/reports/example_timing_2026-10-29/time_examples.R`, run with `Rscript` from the repo root, writing `timing_<tag>.txt`:

```r
library(forestsearch)
tag <- commandArgs(trailingOnly = TRUE)[1]
sink(sprintf("dev/reports/example_timing_2026-10-29/timing_%s.txt", tag), split = TRUE)
cat("installed:", packageDescription("forestsearch")$Built, "\n")
t1 <- system.time(example(run_null_calibration, package = "forestsearch", run.donttest = TRUE, echo = FALSE))
cat(sprintf("run_null_calibration example: elapsed %.1f s (user %.1f s)\n", t1[["elapsed"]], t1[["user.self"]]))
# suggest_thresholds: the two calls of its example, timed separately (measurement only)
d_s <- d_eff_survival(n_sg = 60, prop_cens = 0.45)
t2 <- system.time(suggest_thresholds(d_eff = d_s, N = 700))
d_b <- d_eff_binary(n_sg = 150, p_event = 0.30)
t3 <- system.time(suggest_thresholds(d_eff = d_b, N = 500, L_eff_C = 0.220, L_eff_alpha = 1.298))
cat(sprintf("suggest_thresholds survival call: elapsed %.1f s; binary call: elapsed %.1f s\n", t2[["elapsed"]], t3[["elapsed"]]))
sink()
```

Run `timeout 10m Rscript dev/reports/example_timing_2026-10-29/time_examples.R before`. Commit the script and `timing_before.txt` by path.

## 3. The edit

In the `@examples` block of `run_null_calibration()` (`R/calibrate_null_correction.R`), change only these lines, leaving `sim_null()` and everything else as it is:

- in `fs_args`: `use_lasso = TRUE,` → `use_lasso = FALSE,` and `use_grf = TRUE,` → `use_grf = FALSE,`;
- in the call: `N_values = c(200, 500, 1000),` → `N_values = c(200, 500),` and `n_sims = 50,` → `n_sims = 10,`;
- insert, immediately before the line `#' cal <- run_null_calibration(`:
  ```
  #' # Illustrative size only (the function requires n_sims >= 10 and at least
  #' # two N values); a calibration for use takes n_sims of 50 or more over the
  #' # trial sizes of interest.
  ```

`devtools::document()` → `git status --porcelain` must show exactly `R/calibrate_null_correction.R` and `man/run_null_calibration.Rd` modified (plus the `DESCRIPTION` roxygen re-stamp if it appears — `git checkout -- DESCRIPTION`). `devtools::install(upgrade = FALSE)`.

## 4. Timing after the edit — GATE with one fallback

Run `timeout 10m Rscript dev/reports/example_timing_2026-10-29/time_examples.R after`.

- If the `run_null_calibration` example's elapsed time is **< 5.0 s**: done.
- If it is **≥ 5.0 s**: also change `confounders.name = c("bm1", "bm2", "age", "ecog"),` → `confounders.name = c("bm1", "bm2", "age"),` in `fs_args`, `document()`, reinstall, and run the script again with tag `after2`. Report both figures. If still ≥ 5.0 s, stop there (do not reduce below the function's floors, do not touch the function) and report the number — Larry decides between accepting a small NOTE and `\dontrun`.

The example must run without error in every timing (the script errors out otherwise — that is the gate).

## 5. NEWS, report, commit — GATE

- `NEWS.md`, under `### Documentation` of the development version, one bullet: "The `run_null_calibration()` example runs at an illustrative size (two trial sizes, ten simulations, no LASSO/GRF screening); a calibration for use takes more simulations, as the example now says."
- Report `dev/reports/REPORT_example_timing_2026-10-29.md`: provenance; the diff of the example (`git diff <base> HEAD -- R/calibrate_null_correction.R`); the before/after timing lines verbatim; the `suggest_thresholds` per-call times (measurement only, Larry's decision); deviations; closing block.
- Commits by explicit path: `R/calibrate_null_correction.R` + `man/run_null_calibration.Rd` ("Example of run_null_calibration() at illustrative size"); then `NEWS.md`, the record directory and the report.
- Closing block (run after the last commit, output pasted into the report in a final report-only commit):
  ```
  git diff --name-only <base> HEAD | grep -v -E '^(dev/|NEWS\.md$|R/calibrate_null_correction\.R$|man/run_null_calibration\.Rd$)'
      # must print nothing
  git diff <base> HEAD -- R/calibrate_null_correction.R | grep -E '^[-+]' | grep -v -E "^[-+]#'" | grep -v -E '^(\+\+\+|---)'
      # must print nothing: every changed line in the file is a roxygen line
  git status --porcelain | grep -v '^?? '
      # must print nothing
  ```
- Copy the report and `timing_after.txt` (and `timing_after2.txt` if made) to `~/Downloads`; print `ls -la` of the copies.

## 6. Out of scope

`suggest_thresholds()` and its example; any function body; `R CMD check`; push; merge.
