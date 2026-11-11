# REPORT — A-002 magnitude check (bootstrap id merge), 2026-11-10

**Outcome: STOPPED at Step 4.** Both bootstrap arms reached the 30 min wall-clock cap
before finishing; per the task, the runs were stopped and the task ended there. No
bootstrap object was produced, so the results table and the verdict are not available.

## 1. Header

| Item | Value |
|---|---|
| `master` SHA | `217656e2b7142ad62932d12da91698752f174a96` (task-doc commit; `R/` identical to the audited `5b994177`) |
| Version (DESCRIPTION) | 0.4.1.9002 |
| Application document | `quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd` (the named file; present, tracked, one `forestsearch()` call at line 539 and one `forestsearch_bootstrap_dofuture()` call at line 736) |
| `nb_boots` | `NB = 1000L` (chunk `parallel-setup`) |
| `seed` (bootstrap) | not supplied by the document; the formal default `seed = 8316951L` of `forestsearch_bootstrap_dofuture()` applies |
| Fit wall clock | arm A 39.1 s, arm B 38.9 s (Step 3 runs) |
| Arm A wall clock (fit + bootstrap) | 1800.1 s: stopped by `timeout 1800` (exit 124) |
| Arm B wall clock (fit + bootstrap) | 1800.1 s: stopped by `timeout 1800` (exit 124) |

`forestsearch()` arguments as resolved in the script (`sg_focus = "effMaxSG"`):

- `df.analysis`: `survival::gbsg` with `id <- seq_len(nrow)`, `time_months <- rfstime / 30.4375`, `grade3 <- ifelse(grade == "3", 1, 0)` (686 rows).
- `outcome.name = "time_months"`, `event.name = "status"`, `treat.name = "hormon"`, `id.name = "id"`, `confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er")`.
- `is.RCT = TRUE`, `seedit = 8316951`, `est.scale = "hr"`.
- `use_lasso = FALSE`, `use_grf = FALSE`, `use_dina = FALSE`, `conf.cont_jcuts = list(er = 10, pgr = 10)`, `collapse_cuts = TRUE`, `conf_force = c("er <= 0", "pgr <= 0")`.
- `subgroup_method = "consistency"`, `sg_focus = "effMaxSG"`, `selection_rule = "neighborhood"`, `effect_neighborhood = 0.10`.
- `mr_inference = TRUE`, `mr_inference_args = list(draws = 5000L)`.
- `n.min = 60`, `d0.min = 10`, `d1.min = 10`, `maxk = 2`.
- `hr.threshold = 1.0`, `hr.consistency = 1.0`, `pconsistency.threshold = 0.90`, `minp = 0.025`, `use_twostage = TRUE`, `stop_threshold = NULL`, `consistency_method = "resample"`, `fs.splits = 1000`, `max.minutes = 3`.
- `show_candidate_summary = TRUE`, `details = TRUE`, `quiet = FALSE`, `plot.sg = TRUE`.
- `parallel_args = list(plan = "sequential", workers = n_cores, show_message = TRUE)` (document: `"multisession"`).

Bootstrap call: `forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = NB, show_three = FALSE, details = FALSE, parallel_args = list(plan = "sequential", workers = n_cores, show_message = TRUE))` (document: `"multisession"`).

How it ran:

- Scratch directory `/tmp/a002_2026-11-10`. Every run used `Rscript --vanilla` from its own subdirectory, with `pkgload::load_all()` on the repo tree (arm A) or on `/tmp/a002_2026-11-10/fs_patched` (arm B), and the sequential plan.
- The scratch script `/tmp/a002_2026-11-10/app.R` was generated mechanically from the document's chunks `setup`, `focus-config`, `focus-profile`, `parallel-setup`, `data`, `forestsearch-fit` and `bootstrap`. It has three marked changes: `library(forestsearch)` became `pkgload::load_all()` as the task requires, the two `parallel_args` plans became `"sequential"`, and the bootstrap chunk stops after its timing line, dropping the `summarize_bootstrap_results()` call that the two calls do not need. The `.qmd` was not edited.
- Arms A and B ran concurrently, as two separate processes on a 128-core machine, each under `timeout 1800`.
- Not run: the test suite, `check()`, `document()`, `build()`, renders, installs. Nothing under `R/`, `tests/`, `man/`, `vignettes/`, `inst/`, `quarto/`, `DESCRIPTION` or `NAMESPACE` changed.

## 2. Gate checklist

| Gate | Result | Evidence |
|---|---|---|
| Step 0.1: `git status --porcelain` empty on `master` | PASS | 0 lines |
| Step 0.2: task doc committed on `master` | PASS | `217656e2 dev: add A-002 magnitude task` |
| Step 0.3: `fs-a002-check` absent, then created from `master` | PASS | created from `217656e2` |
| Step 1.1: exactly one application candidate | PASS | the named document exists and has both calls; the fallback search was not needed |
| Step 2.2: flag table has one row per row of `df`, same order | PASS | `extract_subgroup()` builds `df_flag <- data.frame(id = df$id, treat.recommend = ...)` (subgroup_consistency_helpers.R:922-950), one row per `df` row in order. The `df` it receives is `df.fs`, passed unchanged by `subgroup.consistency()` to `sg_consistency_out()` (subgroup_consistency_main.R:1141-1142, no reassignment of `df` in the function), which calls `extract_subgroup(df, ...)` at subgroup_consistency_helpers.R:1050. `df.fs` is built row-wise from `forestsearch()`'s `df` (forestsearch_main.R:3390-3391), and that `df` is not reassigned between :3338 and the merge at :3722 |
| Step 3.2: equivalence on original data | PASS | `sg.harm` identical (`{er <= 0} & {size <= 35}`). `df.est`: 686 rows in both, ids identical after sorting, `treat.recommend` identical, 61 in H. `df.predict`: same |
| Step 4: arm A within 30 min | **FAIL** | stopped by the cap at 1800.1 s; no bootstrap object |
| Step 4: arm B within 30 min | **FAIL** | stopped by the cap at 1800.1 s; no bootstrap object |
| Step 4.4: per-replicate selection identical in the two arms | not run | no bootstrap object in either arm |

## 3. The patch (Step 2)

`/tmp/a002_2026-11-10/fs_patched` is `git archive HEAD` of `217656e2`. Only `R/forestsearch_main.R` differs. The `df.predict` line is the DINA branch's line (forestsearch_helpers.R:1666) copied unchanged; `sg.harm` in this block is `grp.consistency$sg.harm` (forestsearch_main.R:3701), the same label vector that `df.test` already passes to `get_dfpred()` at :3738-3742.

```diff
--- R/forestsearch_main.R (repo, 217656e2)
+++ R/forestsearch_main.R (fs_patched)
@@ -3719,7 +3719,8 @@
              "ambiguous merge with internal subgroup-flag tables.",
              call. = FALSE)
       }
-      df.est_out <- merge(df, temp, by.x = id.name, by.y = "id", all.x = TRUE)
+      df.est_out <- df
+      df.est_out$treat.recommend <- temp$treat.recommend
 
       # Return df.predict
       if (!is.null(df.predict)) {
@@ -3730,7 +3731,7 @@
                "ambiguous merge with internal subgroup-flag tables.",
                call. = FALSE)
         }
-        df.predict_out <- merge(df.predict, temp, by.x = id.name, by.y = "id", all.x = TRUE)
+        df.predict_out <- get_dfpred(df.predict, sg.harm, version = 2)
       }
 
       # Return df.test
```

## 4. Results

Not available. Both arms stopped at the cap before `forestsearch_bootstrap_dofuture()` returned, so no H0 or H2 estimate exists for either arm. The number of replicates completed before the stop is not recorded: with `details = FALSE` the bootstrap prints nothing per replicate, and each log ends at the bootstrap's start lines.

What the runs do establish:

- The fit took 39 s in each arm, leaving about 1,761 s for the bootstrap, which did not finish 1,000 replicates in that time. Under these conditions (sequential, two arms running at once), the average replicate therefore took more than about 1.76 s.

## 5. Verdict

Not reached: the task stopped at the Step 4 cap, so neither MATERIAL nor IMMATERIAL can be stated.
