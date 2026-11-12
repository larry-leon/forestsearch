# REPORT — empty identified set: detection flag, PPV, and the CV agreement metrics, 2026-11-11

Task: `dev/tasks/TASK_empty_set_nofind_2026-11-11.md`, committed as received on `fix/npv-nofind-instances` (`1f1af85b`).
It was run as Step 0.3 of `TASK_classification_outside_package_2026-11-11_v2.md`, which requires this task first.
Outcome: **all gates passed.** The branch is not pushed and not merged.

Definition (León et al. 2024, Section 3): with no subgroup identified, the estimated complement is everyone. Under an alternative that replicate has sens = 0, ppv = 0 (0/0 convention), spec = 1, npv = |Hc|/N; under the null spec = 1, npv = 1, sens and ppv NA.

## 1. Provenance (Step 0)

| Item | Value |
|---|---|
| Host | `pop-os`, R 4.6.2 |
| `git status --porcelain` at start | empty (**gate passed**) |
| Working branch | `fix/npv-nofind-instances` (exists; checked out) |
| Start SHA (HEAD before step 4) | `96feee06` (`dev: NPV no-find instances report`) |
| Gate: `grep -cE 'npv +<- +sum\(!true_H\) / nrow\(df\)' R/run_simulation_analysis.R` | 4 (**passed**, ≥ 4) |
| Gate: `grep -c 'hc_share' R/fs_mr_oc_summary.R` | 3 (**passed**, ≥ 1) |
| Task commit | `1f1af85b dev: add empty-set no-find task` (one file, `cmp`-identical to `~/Downloads`) |
| `DESCRIPTION` at start | 0.4.1.9005 (as expected); now 0.4.1.9006 |
| Task file arrival | `~/Downloads/TASK_empty_set_nofind_2026-11-11.md` appeared at 17:24 (the v2 task was received at 16:02 without it; the session waited for it, doing read-only preparation only) |
| Fix commit / report commit | given in the closing message (a report cannot hold its own SHA) |

## 2. Summary for Larry

| Case | Fixed / BLOCKED | Can the real pipeline reach it? (how) | Outputs changed (old → new) |
|---|---|---|---|
| **E1.** FS site 2: result rows, no `M.*` factor recovered (`R/run_simulation_analysis.R:1048-1074` at the start SHA) | **Fixed** (Rule E) | **No.** `forestsearch()` attaches `df.est$treat.recommend` whenever `grp.consistency$out_sg$result` has rows (`forestsearch_main.R:3700`, `:3725-3726`), so the extractor never takes the `M.*` branch from `run_simulation_analysis()`; and `M.1` of a result row is never NA (`subgroup_consistency_helpers.R:1432, :1559, :1904, :2038`). Reached only by direct calls of the internal extractor (as the A-003 tests do) | `any.H` 1 → 0. Every other output already equalled the no-find path |
| **E2.** FS found path, `df.est$treat.recommend` with no harm member while result rows exist | **Fixed** (Rule E) | **No.** Rows exist only on the consistency engine's final return, which always builds `treat.recommend` with ≥ 1 harm member (`subgroup_consistency_helpers.R:934-942`; selected set ≥ 10 members or ≥ n.min + 1 on maxeff rows) | alternative: `any.H` 1 → 0; `ppv` NA → 0; `hr.H.true`, `hr.Hc.true`, `ahr.H.true`, `ahr.Hc.true`, `cde.H.true`, `cde.Hc.true`, `ahr.Hc.hat`, `cde.Hc.hat` value → NA. Null: `any.H` 1 → 0; `hr.Hc.true`, `ahr.Hc.true`, `cde.Hc.true`, `ahr.Hc.hat`, `cde.Hc.hat` value → NA. Edge shapes (§3.4): `spec` NA → 1 when every subject is truly in H; `hr.Hc.hat` NA → `hr.itt` when N ≤ 10; `treat.recommend` all NA or only 1 and NA: **error → no-find row** |
| **E3.** FS found path through `M.*` with an indicator that has no member (no `df.est$treat.recommend`) | **Fixed** (Rule E) | **No** (as E1: the `M.*` branch is not reached from the pipeline; real `M.*` labels are brace labels, for which `create_subgroup_indicator()` returns all 1, not all 0) | as E2 |
| **E4.** GRF found path: `sg.harm.id` set, `data$treat.recommend` with no harm member (`:1409-1430` at the start SHA) | **Fixed** (Rule E) | **Yes, on one non-default path only:** survival, `grf_params$grf_selection = "frontier"`, a factor/character confounder whose labels are numeric strings, and a numeric covariate whose lowest level rounds down at 7 significant digits. `grf_main.R:226` builds GRF's covariate matrix with `apply(..., as.numeric)` (which formats numbers through `as.matrix`), while membership is evaluated on the raw data (`grf_subgroup_labels.R:225`), so a `<=` cut at that level can select nobody in the data. The default pipeline (tree selection) cannot reach it | as E2 (measured on the GRF extractor: identical list) |
| **C1.** CV single-level branch: folds put nobody in H, full-data assignment has both levels (`R/forestsearch_cross_validation.R:1415-1420` at the start SHA) | **Fixed** (Rule C) | **Yes.** A K-fold/tenfold run whose folds identify nothing under an identified primary fit (the "all-ITT sims" shape, `test-cv-mixed-and-fatal.R` row 3) | (`sens_H`, `sens_Hc`, `ppv_H`, `ppv_Hc`), output order: (1, NA, \|orig H\|/N, NA) → (0, 1, 0, \|orig Hc\|/N); test fixture (35 of 100 in H): (1, NA, 0.35, NA) → (0, 1, 0, 0.65) |
| **C2.** CV single-level branch: folds put everybody in H, full data has both levels | **Fixed** (Rule C) | In principle (every held-out subject assigned to H); not seen in a test before this task | (1, NA, \|orig H\|/N, NA) → (1, 0, \|orig H\|/N, 0) |
| **C3.** CV single-level branch, full-data assignment with a single observed level | unchanged, as the task directs | yes | none (pinned by a new expectation) |

- **E1-E3 change nothing that `run_simulation_analysis()` produces today**: the cases are unreachable from it (two independent traces each, §3). The change makes the internal extractor follow Rule E on every input. E4 changes only replicates on the non-default frontier path above.
- **E4 was first classified BLOCKED, then Fixed.** §3.3 gives the reasoning and the reviewer disagreement; it is OPEN ITEM 1 for Larry's decision.

## 3. Part 1 trace

Each case was traced from the source by one read-only agent and checked by a second, independent one (scratch: `scratchpad/trace/{fs_site2,fs_found_empty,grf_found_empty}.md` and `verify_*.md`). Line numbers are at `1f1af85b` (= the start SHA for `R/`).

### 3.1 How the extractors are reached

- FS (`.extract_fs_estimates_gen`, `:921`) is called only from `.run_fs_analysis_gen` (`:842-858`), for the FS, FSlg and `methods=` arms (`:591`, `:616`, `:678`). Its inputs: `fs_res = grp.consistency$out_sg$result` and `fs_full = forestsearch()` result, both only when `has_result` (`:823-831, :844-845`); `forestsearch()` errors become `NULL` (`:818-821`).
- GRF (`.extract_grf_estimates_gen`, `:1306`) is called only from `.run_grf_analysis_gen` (`:1198`, `:1263`), when `methods` is NULL and `run_grf` is TRUE (`:635-653`); `methods = "grf-<focus>"` arms run `forestsearch()` in GRF mode through the FS extractor instead.

### 3.2 FS (E1-E3)

- **Entry.** The found path at `:1041` needs result rows or `.has_treat_rec` (`:1005-1009`). The `M.*` branch (`:1048-1074`) needs, in addition, `fs_full$df.est` without `treat.recommend`.
- **Producer.** `forestsearch()` (`R/forestsearch_main.R`) has these returns: zero-row input (`:1920-1933`), DINA mode (not found: `forestsearch_helpers.R:1646-1654`; found: `:1665, :1675-1687`, `df.est` from `get_dfpred()`, which always assigns `treat.recommend`, `forestsearch_helpers.R:117, :137`), GRF mode (not found `:1967-1975`; found `:1987-1988, :2002-2012`), `get_FSdata` or `subgroup.search` failure (`:3310-3333`, `:3495-3507`), and the consistency path (`:3946-3986`).
  - `out_sg$result` has rows only on the consistency path's final return (`subgroup_consistency_main.R:1194-1208`); every early return there has `out_sg = NULL` (`:459-470, :537-548, :557-568, :576-592, :1057-1093`). In DINA/GRF mode `out_sg` has no element whose name starts with `result`.
  - On that return, `sg.harm <- out_sg$sg.harm_label` is never NULL (`subgroup_consistency_helpers.R:926`), so `forestsearch_main.R:3700` holds and `df.est$treat.recommend` is attached row-aligned (`:3725-3726`), from `ifelse(sg.harm.id == 1, 0, 1)` (`subgroup_consistency_helpers.R:934-942`) with at least one 0: the selected set is re-checked at N.x ≥ 10 on the same data (`:1443-1458`, `:1819-1833`), and maxeff skip rows carry the search's nx > n.min with events in both arms (`subgroup_search.R:699-722`).
  - Result rows always set `M.1` to a label (`subgroup_consistency_helpers.R:1432, :1559, :1904, :2038`); unfilled `M.*` are `""`, never NA.
- **Verdict (both tracers and both verifiers): UNREACHABLE** from `run_simulation_analysis()` for E1, E2 and E3. Reachable by direct calls of the internal extractor. The guard is not triggered: on the only path with result rows, another object (`sg.harm`, `grp.consistency$sg.harm.id`, `df.est$treat.recommend`) always carries the same non-empty set, so if one of these cases were ever reached there it would be a defect elsewhere. The fix ignores `fs_res`'s `M.*` labels when `df.est$treat.recommend` exists; that is consistent only because the case cannot be reached with a non-empty set elsewhere (noted, review F9).

### 3.3 GRF (E4), and why it is Fixed rather than BLOCKED

- **Producer returns** (`grf.subg.harm.survival()`, `R/grf_main.R`; `grf.subg.harm.glm()`, `R/grf_subg_harm_glm.R`): errors → `grf_est = NULL` (no-find, `:1382`); null results (`sg.harm.id = NULL`; the `create_null_result` call at `grf_main.R:235` actually errors at `grf_helpers.R:535` and lands in the same no-find); **tree success** (`grf_main.R:421`): `treat.recommend` is 0 exactly for the selected leaf, predicted from the same tree on the same X, so the set has `Nsg` ≥ max(1, n.min) members — **cannot be empty**; **frontier success, survival** (`grf_main.R:333`): `treat.recommend = .grf_evaluate_subgroup(sg_def, data)` (`:308-309`) — **can be empty** (mechanism in §2, E4; reproduced by both agents with hand-built data, e.g. a selected `zx <= -0.81649660000000002` with 117 members on X and 0 in the data); GLM frontier and tree successes cannot be empty (`.build_grf_X` codes columns before `as.matrix`, `grf_subg_harm_glm.R:891-894`; `|sg_idx|` ≥ n.min, `:590`).
- **Guard.** The task: BLOCKED "if the trace shows that a case can be reached while a non-empty identified set exists for that replicate (membership available from another object)".
  - The tracer read it as BLOCKED: the same result object records `grf.gsub$Nsg` = 117 and the frontier/candidate rows carry that size, and membership can be rebuilt by evaluating the cut on `apply(grf_est$data[, confounders], 2, as.numeric)`.
  - Its verifier, and the independent spec reviewer, read it as REACHABLE-EMPTY: no returned object carries a non-empty **per-subject membership** — `data$treat.recommend` is empty, and `sg.harm.id` evaluated on the returned data gives the same empty set; `Nsg` is a count; the rebuild needs `confounders.name`, which is in neither the GRF object nor the extractor's arguments, and the producer-time `getOption("digits")`.
  - The same producer defect already scores as a no-find on the other GRF route: `methods = "grf-<focus>"` runs `forestsearch()` in GRF mode (frontier by default), whose `df.est$treat.recommend` is then all 1 with `fs_res` NULL, so the FS extractor's no-find path applies (`any.H` 0, `ppv` 0) before and after this task.
  - **Decision: Fixed (Rule E)**, on the literal reading ("membership available") and for consistency with the `forestsearch()` GRF route. The first draft of this change had BLOCKED it; the reversal and both readings are recorded here and in OPEN ITEM 1. The producer defect itself (X coding versus evaluator coding) is out of scope; it also shrinks or grows non-empty sets on the same path (OPEN ITEM 2).

### 3.4 Every output that changes for an empty-set replicate

Measured by running the start-SHA extractors (sourced from `git archive 96feee06`) and the branch's on the same hand-built inputs (N = 200, 50 truly in H, survival with `loghr_po` and `theta_*`); the reviewer's sweep over 1,224 FS input shapes (survival and GLM) gives the same lists.

| Shape | Old → new |
|---|---|
| E1 (site 2) | `any.H` 1 → 0 only |
| E2, E3, E4, alternative | `any.H` 1 → 0; `ppv` NA → 0; `hr.H.true` 1.46439 → NA; `hr.Hc.true` 0.584825 → NA; `ahr.H.true` 1.01831 → NA; `ahr.Hc.true` 0.98807 → NA; `ahr.Hc.hat` 0.995545 → NA; `cde.H.true` 0.853173 → NA; `cde.Hc.true` 1.22825 → NA; `cde.Hc.hat` 1.12071 → NA |
| E2, E3, E4, null | `any.H` 1 → 0; `hr.Hc.true`, `ahr.Hc.true`, `cde.Hc.true`, `ahr.Hc.hat`, `cde.Hc.hat` value → NA (`ppv` NA in both) |
| no harm column | `any.H` 1 → 0; `ahr.Hc.hat`, `cde.Hc.hat` value → NA |
| every subject truly in H | also `spec` NA → 1 |
| N ≤ 10 | also `hr.Hc.hat` NA → `hr.itt` |
| E2/E3 with membership only NA, or NA and "not in H" | **error → the no-find row.** NA membership counts as "no member" (stated interpretation). The start-SHA code stopped with "missing value where TRUE/FALSE needed" (`:1133`) when the harm column was present; without it, `size.Hc` N − #NA → N, `hr.Hc.hat` → `hr.itt`, `attr(, "sg_hat")` NA → 0 |
| zero-row `df.est` with a `treat.recommend` column | error ("replacement has 0 rows") → the no-find row |

- Unchanged in every empty shape: `size.H` (0), `size.Hc` (N), `sens`, `npv`, `hr.itt`, `hr.adj.itt`, `hr.H.hat`, `ahr.H.hat`, `cde.H.hat`, `p.cens`, `taumax`, `attr(, "sg_hat")` (all 0). `hr.Hc.hat` is `identical()` to `hr.itt` in the old found path for N > 10.
- **Non-empty sets: byte-identical.** All 536 non-empty FS shapes of the reviewer's sweep, and the FS (`df.est` and `M.*`) and GRF non-empty cases measured here, give `identical()` columns and `sg_hat`.
- **Not changed (outside the extractors):** the runner columns `sg.def`, `sg.n_factors`, `n_candidates_*`, `n_passed` and `grf_selected_depth` (`:866-885`, `:1279-1288`) come from the search object, not from the extractor's path. On E4 they keep the cut text and its factor count (OPEN ITEM 3).

## 4. The change

### 4.1 `git diff 96feee06 -- R/`, verbatim

```diff
diff --git a/R/forestsearch_cross_validation.R b/R/forestsearch_cross_validation.R
index e885f612..e50d5b2a 100644
--- a/R/forestsearch_cross_validation.R
+++ b/R/forestsearch_cross_validation.R
@@ -1414,10 +1414,26 @@ forestsearch_KfoldOut <- function(res, details = FALSE, outall = FALSE,
       ppvHc <- tabit[2, 2] / sum(tabit[2, ])
     } else {
       tabit <- with(df_CV, table(treat.recommend, treat.recommend.original))
-      sensH <- if (nrow(tabit) > 0 && ncol(tabit) > 0) tabit[1, 1] / sum(tabit[, 1]) else NA
-      sensHc <- NA
-      ppvH <- if (nrow(tabit) > 0) tabit[1, 1] / sum(tabit[1, ]) else NA
-      ppvHc <- NA
+      if (nrow(tabit) == 1L && identical(colnames(tabit), c("0", "1")) &&
+          all(tabit[1, ] > 0)) {
+        # The folds assign a single level while the full-data assignment has
+        # both: the two-level definitions above, with the missing fold level
+        # counted as zero subjects and 0/0 = 0.  Folds identifying nothing
+        # give 0, 1, 0 and the full-data complement share; folds putting
+        # everyone in H give 1, 0, the full-data H share and 0.
+        tab2 <- with(df_CV, table(factor(treat.recommend, levels = c(0, 1)),
+                                  factor(treat.recommend.original, levels = c(0, 1))))
+        .agree <- function(num, den) if (den > 0) num / den else 0
+        sensH  <- .agree(tab2[1, 1], sum(tab2[, 1]))
+        sensHc <- .agree(tab2[2, 2], sum(tab2[, 2]))
+        ppvH   <- .agree(tab2[1, 1], sum(tab2[1, ]))
+        ppvHc  <- .agree(tab2[2, 2], sum(tab2[2, ]))
+      } else {
+        sensH <- if (nrow(tabit) > 0 && ncol(tabit) > 0) tabit[1, 1] / sum(tabit[, 1]) else NA
+        sensHc <- NA
+        ppvH <- if (nrow(tabit) > 0) tabit[1, 1] / sum(tabit[1, ]) else NA
+        ppvHc <- NA
+      }
     }
 
     sens_metrics_original <- c(sensH, sensHc, ppvH, ppvHc)
diff --git a/R/run_simulation_analysis.R b/R/run_simulation_analysis.R
index 766d299e..e3ab7ca9 100644
--- a/R/run_simulation_analysis.R
+++ b/R/run_simulation_analysis.R
@@ -1002,11 +1002,14 @@ run_simulation_analysis <- function(
   # A subgroup is "found" when the consistency result table has rows OR the
   # engine returned a treat.recommend assignment with harm members (GRF / DINA
   # modes, where fs_res is NULL but membership lives in fs_full$df.est).
+  # A treat.recommend assignment with no harm member is an empty identified
+  # set: a no-find, whatever the result table holds.
   .has_treat_rec <- !is.null(fs_full) && !is.null(fs_full$df.est) &&
     "treat.recommend" %in% names(fs_full$df.est) &&
     any(fs_full$df.est$treat.recommend == 0L, na.rm = TRUE)
 
-  if ((is.null(fs_res) || nrow(fs_res) == 0) && !.has_treat_rec) {
+  if ((is.null(fs_res) || nrow(fs_res) == 0 ||
+       "treat.recommend" %in% names(fs_full$df.est)) && !.has_treat_rec) {
     out$hr.Hc.hat <- out$hr.itt
     out$size.H  <- 0L
     out$size.Hc <- nrow(df)
@@ -1054,7 +1057,10 @@ run_simulation_analysis <- function(
     }
     if (length(sg_factors) > 0)
       df$sg_hat <- create_subgroup_indicator(df, sg_factors)
-    else {
+    if (length(sg_factors) == 0 || !any(df$sg_hat == 1L, na.rm = TRUE)) {
+      # No factor recovered, or an indicator with no member: an empty
+      # identified set, scored as the no-find path above scores it.
+      out$any.H <- 0L
       out$hr.Hc.hat <- out$hr.itt
       out$size.H  <- 0L
       out$size.Hc <- nrow(df)
@@ -1406,7 +1412,10 @@ run_simulation_analysis <- function(
   has_tr     <- !is.null(grf_est$data) &&
                 "treat.recommend" %in% names(grf_est$data)
 
-  if (!has_sg || !has_tr) {
+  # A treat.recommend assignment with no harm member is an empty identified
+  # set: a no-find, whatever sg.harm.id holds.
+  if (!has_sg || !has_tr ||
+      !any(grf_est$data$treat.recommend == 0L, na.rm = TRUE)) {
     out$hr.Hc.hat <- out$hr.itt
     out$size.H  <- 0L
     out$size.Hc <- nrow(df)
```

### 4.2 Hunk mapping

| Hunk | Rule | Case |
|---|---|---|
| `R/forestsearch_cross_validation.R @@ -1414,10 +1414,26 @@` | Rule C | C1, C2 (single fold level, both full-data levels observed: `nrow(tabit) == 1`, columns `"0"`, `"1"`, both non-empty); C3 and every other shape keep the old code |
| `R/run_simulation_analysis.R @@ -1002,11 +1002,14 @@` | Rule E | E2: a `df.est$treat.recommend` with no harm member routes to the existing no-find block |
| `R/run_simulation_analysis.R @@ -1054,7 +1057,10 @@` | Rule E | E1 (`any.H <- 0L`) and E3 (an `M.*` indicator with no member takes the same block) |
| `R/run_simulation_analysis.R @@ -1406,7 +1412,10 @@` | Rule E | E4: a GRF `treat.recommend` with no harm member routes to the existing no-find block |

### 4.3 Gates after the edits

| Gate | Result |
|---|---|
| 1. `git diff 96feee06 --stat -- R/` lists only the two files | `R/forestsearch_cross_validation.R` (+20/−4), `R/run_simulation_analysis.R` (+12/−3). **PASS** |
| 2. Every hunk maps to Rule E or Rule C | 4 hunks, all mapped (§4.2). **PASS** |
| 3. `document()` | not run: no roxygen line changed (`git diff 96feee06 -- R/ \| grep -cE "^[-+]#'"` = 0); `git diff 96feee06 --stat -- man/` empty. **PASS** |
| 4. The two Step 0 greps; `grep -c 'spec <- sum(!true_H)' R/run_simulation_analysis.R` | 4, 3; 0. **PASS** |

- `DESCRIPTION` 0.4.1.9005 → 0.4.1.9006 (the only line changed).
- `NEWS.md`: two bullets at the top of `### Behaviour changes` in the development section. Part 1 is the task's sentence verbatim. Part 2 is the task's sentence with the four metric names inserted (`sens_H`, `ppv_H`, `sens_Hc`, `ppv_Hc`, in the order of the task's values) and one added sentence for C2 ("When the folds put everyone in the subgroup, `sens_Hc` and `ppv_Hc` are now 0 (they were NA)"), because the review found the task's sentence silent on C2 and its value order ambiguous against the output order.

## 5. Tests

### 5.1 Updated expectations

| File:line | Old | New |
|---|---|---|
| `tests/testthat/test-cv-mixed-and-fatal.R:142` (now `:143`) | `c(1, NA, 0.35, NA)` | `c(0, 1, 0, 0.65)` (comment above rewritten) |
| `tests/testthat/test-nofind-classification.R:67-71` (A-003 site 2) | `check_any_H = FALSE`, with a comment that site 2 leaves `any.H` at 1 | the flag is now checked (`any.H` 0); comment removed |

### 5.2 New expectations

`tests/testthat/test-nofind-classification.R` (hand-built inputs, no search):

| Block | Checks |
|---|---|
| "empty set: FS site 2 … equals site 1" | Part 3 item 1: alternative and null; `any.H` 0; every column and `sg_hat` `identical()` to the site 1 output for the same data |
| "empty set: GRF found path with no member equals the no-find output" | Part 3 item 2, GRF: `sg.harm.id` set, `treat.recommend` all 1; `any.H` 0; `ppv` 0 (alternative) / NA (null); every column and `sg_hat` `identical()` to the GRF no-find output |
| "empty set: FS found path with no member equals the no-find output" | Part 3 item 2, FS: via `df.est$treat.recommend` all 1, and via an `M.*` factor matching no subject; same checks |
| "non-empty set: FS and GRF found paths keep the found-case formulas" | Part 3 item 3: 60 identified subjects, FS (`df.est`) and GRF: `any.H` 1, `size.H` 60, sens/spec/ppv/npv equal to the formulas computed inline; FS `M.*` path with 50 members: the same |

`tests/testthat/test-cv-mixed-and-fatal.R`, new block "CV agreement: folds identify nothing / everyone, full data has both levels": both rows of the Rule C table (0, 1, 0, 65/100 and 1, 0, 35/100, 0, in output order), the two-level case unchanged from the start SHA (30/35, 55/65, 30/40, 55/60), and a full-data factor with an unused level left as it was (`NaN, NA, 0, NA`).

- **Negative control** (the two files run against `git archive 96feee06`): `test-nofind-classification.R` 20 failed expectations (A-003 site 2 `any.H` 1, site 2 vs site 1 4, GRF empty 5, FS empty 10); `test-cv-mixed-and-fatal.R` 3 (the updated `:143`, the two Rule C rows). The non-empty, two-level and unused-level expectations pass on both.
- No search, no simulation and no parallel plan in any new test.

## 6. Tests and check

| Run | Baseline (`REPORT_npv_nofind_instances_2026-11-11.md` §6, not re-run) | This branch | Result |
|---|---|---|---|
| Changed test files (`load_all` + `test_file`, `NOT_CRAN=true`) | — | `test-nofind-classification.R` 9 blocks, 96 expectations, 0.6 s; `test-cv-mixed-and-fatal.R` 7 blocks, 58 expectations, 14.2 s; **0 failures**, 0 warnings, 0 skips (cap 2 min each) | **PASS** |
| Full suite (`devtools::test()`, once) | 608 blocks, 81 files, 8,218 expectations, 0 failures, 3 skips, 38 warnings; 647 s | **613 blocks, 81 files, 8,270 expectations, 0 failures, 0 errors, 3 skips, 38 warnings**; wall **645 s** (18:04:50-18:15:35; cap 30 min) | **PASS** |
| Difference | — | +5 blocks, +52 expectations, all in the two changed test files; skips the same 3 (`test-cv-no-subgroup-edges.R:250`, `test-cross-outcome-parity.R:262, :358`), warnings unchanged | |
| `devtools::check()` (once) | 0 errors / 0 warnings / 0 notes | **0 errors / 0 warnings / 0 notes**, Status OK; wall **994 s** (R CMD check duration 14 m 46 s; tests 598 s elapsed, vignettes 74 s) (cap 60 min) | **PASS** |

- **Check route:** as `REPORT_npv_nofind_instances_2026-11-11.md` §6: `check(document = FALSE, error_on = "never")` on an export of the working tree (`git ls-files -co --exclude-standard`, 9,626 files) in the session scratchpad; the 7 changed files were `cmp`-identical to the repository's.
- **Timing:** the suite and `check()` ran at the same time (both started 18:04:50), so each wall clock includes contention from the other.
- **Pre-set decision:** no expectation outside the two changed files failed; nothing else needed updating.
- **Commit:** `fix: empty identified set is a no-find on every path; CV agreement when folds find nothing` by explicit path: `DESCRIPTION`, `NEWS.md`, `R/run_simulation_analysis.R`, `R/forestsearch_cross_validation.R`, `tests/testthat/test-nofind-classification.R`, `tests/testthat/test-cv-mixed-and-fatal.R`.

## 7. Independent review

Two read-only reviewers checked the uncommitted change (scratch `scratchpad/es_review/{spec,correct}.md`): one against the task text, one by running the start-SHA and branch functions on 1,224 FS and 1,512 CV input shapes. Neither found a blocking issue; every non-empty FS shape was byte-identical. What changed in response:

| Finding | Change |
|---|---|
| F1 (spec): the GRF BLOCKED call does not follow from the guard as written; the `forestsearch()` GRF route already scores the same defect as a no-find | E4 fixed under Rule E (§3.3); OPEN ITEM 1 |
| S1 (correctness): the Rule C guard tested level names, so a full-data factor with an unused level (one observed level) entered the new branch, which the task says to leave as it is | `&& all(tabit[1, ] > 0)` added; expectation pinned |
| S2 / F7: NEWS value order ambiguous; C2 not stated | metric names inserted, C2 sentence added (§4.3) |
| S3 / F6: the list of changed outputs | §3.4 |
| F3, F4, N1: GRF empty-set test missing; `M.*` non-empty path untested; stale A-003 comment | tests added / updated (§5) |

The first full suite/`check()` pair, run on the version before these changes, was stopped and discarded (suite finished: 612 blocks, 8,257 expectations, 0 failures, 3 skips, 38 warnings, 654 s; `check()` stopped mid-run). Its results are not used.

## OPEN ITEMS

1. **E4 guard reading (Larry's decision).** Fixed here on the literal reading of the guard; the alternative reading (GRF's recorded `Nsg` > 0 is a non-empty set) would have left the GRF extractor unchanged. Reverting is one condition at `R/run_simulation_analysis.R:1417-1418` and the GRF block in `test-nofind-classification.R`. Neither reading changes anything on the default (tree) GRF path.
2. **Producer defect behind E4.** `grf_main.R:226` builds the frontier covariate matrix with `apply(data[, confounders.name], 2, as.numeric)`, which formats numbers to 7 significant digits whenever a confounder is a factor or character column; `.grf_evaluate_subgroup()` evaluates the selected cut on the raw values (`grf_subgroup_labels.R:225`). The GLM path already codes columns first (`grf_subg_harm_glm.R:891-894`). The same mismatch can shrink a `<=` set or grow a `>` set without emptying it. Not changed (no other fix allowed).
3. **Runner columns.** `sg.def`, `sg.n_factors`, `n_candidates_*`, `n_passed`, `grf_selected_depth` (`R/run_simulation_analysis.R:866-885`, `:1279-1288`) are set from the search object whatever the extractor's path, so an E4 replicate (and an existing `forestsearch()` GRF-mode empty selection) shows `any.H` 0 beside a non-empty `sg.def`. Outside the two extractors; not changed.
4. **Pre-existing, unchanged:** a non-empty FS set with NA members and the harm column present still stops at `:1133` (64 shapes in the reviewer's sweep, both versions); CV fold assignments {1, NA} or {0, NA} count NA as a level (`n_sgfound`, `forestsearch_cross_validation.R:1363`) and stop with "subscript out of bounds" at `:1412`, both versions.
5. **CV downstream summaries.** `forestsearch_Kfold()$sens_summary` and `forestsearch_tenfold()`'s `sens_out` and median `sens_summary` (`:616`, `:1201`, `:1229`) now include C1 sims in the `sens_Hc` / `ppv_Hc` medians (`na.rm = TRUE`) and give them `sens_H` 0 instead of 1; readers: `R/cv_summary_tables.R:403, :558, :669, :757`, `R/plot_subgroup_results_forestplot.R:613, :1523`. No stored CV output was regenerated.
6. **NEWS heading** still reads `0.4.1.9002 (development)` (A-003 open item 6).
