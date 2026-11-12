# REPORT detail — CRAN preparation, part 4: per-fix evidence, 2026-11-12

Companion to `REPORT_cran_prep_code_fixes_2026-11-12.md`. One section per fix, in fix order. Each was drafted by a read-only agent from that fix's evidence (implementer notes and runs, verifier verdicts and runs, the main session's red/green runs on the repository) and checked line by line against git (`git show <sha>:<file>`); where a source disagreed with git, git was followed and the section says so under "Corrections to the notes".

- Line numbers: start SHA `76a62a71`, each fix's own commit, and the head of the fixes `61c15088` (`R/` is identical at `88302212`).
- Scratch paths are under `/tmp/cran_prep_code_fixes_2026-11-12/` (`evidence/F*`, `verify/F*`, `main/F*`, `patches/`).
- "Main red/green" = the main session on the repository: only the fix's test file applied (red), then the code (green), at the commit preceding the fix.
- NEWS line numbers cited here are at `88302212`. Four of those bullets (F2, F9, F10, F13) were reworded afterwards in `a577b59d` (main report §7 and OPEN ITEMS 16), so a NEWS wording problem noted in a section below is resolved there.

### F1 — A-147: bootstrap summary bins match their labels

**Status:** PASS. Commit `0dc99947` ("fix: bootstrap summary bins match their labels (A-147)"), parent `76a62a71` (start SHA).

| file (`git show --stat 0dc99947`) | change |
|---|---|
| `R/summarize_bootstrap_subgroups.R` | 6 +- (numstat: 4 insertions, 2 deletions) |
| `tests/testthat/test-bootstrap-summary-bins.R` | new, 100 lines |
| total | 2 files, 104 insertions, 2 deletions |

- No `man/`, `NAMESPACE`, `DESCRIPTION` or `NEWS.md` change in this commit.
- The committed diff equals the gated patch `patches/F1.patch` apart from `index`/`@@` lines. I checked this myself (`diff` empty, exit 0), because `main/F1/` holds no `patch_equal.txt` or `diffstat.txt`.
- `R/summarize_bootstrap_subgroups.R` and the test file are unchanged from `0dc99947` to `61c15088` and to `88302212` (`git diff` empty).

#### Site, before and after

Both `cut()` calls get `right = FALSE`. `include.lowest = TRUE` is kept.

| call | start `76a62a71` | head `61c15088` |
|---|---|---|
| Pcons `cut()` (`sg_found$Pcons_bin <- cut(`) | `R/summarize_bootstrap_subgroups.R:417-422` | `:417-423` (new line `:422`) |
| `size_breaks` / `size_labels` | `:441` / `:442` | `:442` / `:443` |
| N_sg `cut()` (`sg_found$N_sg_bin <- cut(`) | `:444-449` | `:445-451` (new line `:450`) |

Before, `76a62a71`, `R/summarize_bootstrap_subgroups.R:420-422` and `:447-449`:

```r
      labels = c("<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.85", "0.85-0.95", ">=0.95"),
      include.lowest = TRUE
    )
...
      labels = size_labels,
      include.lowest = TRUE
    )
```

After, `61c15088`, `R/summarize_bootstrap_subgroups.R:420-423` and `:448-451`:

```r
      labels = c("<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.85", "0.85-0.95", ">=0.95"),
      include.lowest = TRUE,
      right = FALSE
    )
...
      labels = size_labels,
      include.lowest = TRUE,
      right = FALSE
    )
```

- Before the fix, the bins were right-closed `(a, b]`, so a value on a break was counted in the bin below its label. After the fix they are left-closed `[a, b)`, and the top bins `[0.95, 1]` and `[300, Inf]` stay closed.
- The CV summary already uses this closure: `R/summarize_cv_results.R:764-765` at `76a62a71` has `include.lowest = TRUE, right = FALSE`.

#### Red and green

Filter `bootstrap-summary-bins`. It matches only the new file.

| run | where | result | wall |
|---|---|---|---|
| implementer red | scratch copy of `76a62a71` (`wt/F1`) with only the test file; `evidence/F1/red.txt` | `bootstrap-summary-bins: 12..3456`: 8 expectations, 6 failures, 2 passes, no errors. Failure 1 (`test-bootstrap-summary-bins.R:47:3`): `` `actual`:   "0.85-0.95" "0.85-0.95" "0.8-0.85"  ">=0.95" "0.85-0.95" `` vs `` `expected`: ">=0.95"    "0.85-0.95" "0.85-0.95" ">=0.95" "0.85-0.95" `` | 3.04 s |
| implementer green | `wt/F1` with the patch; `evidence/F1/green.txt` | `bootstrap-summary-bins: ........`, DONE. Silent reporter: 3 blocks, nb 1 + 1 + 6 = 8, failed 0, error FALSE, warning 0, skipped FALSE, passed 8 | 2.84 s (summary reporter), 2.82 s (silent reporter) |
| verifier red | start export plus the new test (`verify/F1/pre`); `verify/F1/test_pre.txt` | `TOTAL nb= 8  failed= 6  passed= 2  error= FALSE  warning= 0`; same six failures as the implementer's | 3.60 s |
| verifier green | `verify/F1/post`; `verify/F1/test_post.txt` | `TOTAL nb= 8  failed= 0  passed= 8  error= FALSE  warning= 0` | 3.08 s |
| main-session red | repository at `76a62a71` with only the patch's `tests/` hunks; `main/F1/red.txt` | `bootstrap-summary-bins: 12..3456`: 6 failures (`:47:3`, `:62:3`, `:89:3`, `:92:3`, `:98:3`, `:99:3`), 2 passes; messages identical to the implementer's | 3.0 s |
| main-session green | repository after the code hunks; `main/F1/green.txt` | `bootstrap-summary-bins: ........`, DONE, no failures | 2.9 s |

- The other five red failures, verbatim and shortened:
  - `:62:3`: `` `actual`:   "<50"   "<50" "50-99"   "200-299" "200-299" `` vs `` `expected`: "50-99" "<50" "100-149" ">=300"   "200-299" ``
  - `:89:3` and `:92:3`: Count `0 0 0 1 3 1` vs `0 0 0 0 3 2`, and `2 1 0 0 2 0` vs `1 1 1 0 1 1`.
  - `:98:3` and `:99:3`: the matching Percent failures.
- Every red failure is an on-break value counted one bin below its label. The two expectations that pass on red are the label-order checks.

#### Expected values

All expected values are in the test file at `0dc99947`.

| expected value (test line) | how it was derived |
|---|---|
| `c(">=0.95", "0.85-0.95", "0.85-0.95", ">=0.95", "0.85-0.95")` for Pcons 0.95, 0.90, 0.85, 1.00, 0.949 (`:43`) | By hand from the labels: "a-b" holds a <= v < b, "<a" holds v < a, ">=a" holds v >= a (header comment `:8-11`). One line per value at `:36-42`. These are the spec's values. |
| `c("50-99", "<50", "100-149", ">=300", "200-299")` for N_sg 50, 49, 100, 300, 299 (`:58`) | Same rule, one line per value at `:52-57`. These are the spec's values. |
| consistency labels `c("<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.85", "0.85-0.95", ">=0.95")` (`:76-79`) | The labels documented in the code (`R/summarize_bootstrap_subgroups.R:420` at `76a62a71`). Comment `:75`. **The verifier flagged this check:** it passes before and after the fix, and its comment does not say so. |
| size labels `c("<50", "50-99", "100-149", "150-199", "200-299", ">=300")` (`:80-83`) | The labels in the code (`:442` at `76a62a71`). Comment `:75`. **Same verifier flag.** |
| consistency Count `c(0L, 0L, 0L, 0L, 3L, 2L)` (`:89`) | Hand tally of the per-value bins: ">=0.95" holds 0.95 and 1.00; "0.85-0.95" holds 0.90, 0.85 and 0.949 (`:85-88`). |
| size Count `c(1L, 1L, 1L, 0L, 1L, 1L)` (`:92`) | Hand tally: 49, 50, 100, none, 299, 300 (`:90-91`). |
| Percent `c(0, 0, 0, 0, 60, 40)` and `c(20, 20, 20, 0, 20, 20)` (`:98`, `:99`) | 100 x Count / n_found with n_found = 5 (`:94-97`). This is the definition in the function body (`:429`, `:456` at `76a62a71`), which the fix does not touch. |

- Inputs only, not expected values: the hold points N_sg = 120 (`:44`, an interior point of "100-149") and Pcons = 0.92 (`:59`, an interior point of "0.85-0.95"), and the helper's `Count == 1L` selector (`:29-30`).
- The verifier re-derived all ten per-value bins, the counts and the percents by hand. It found that every value has a derivation comment and that none needs the fixed code.

#### Unaffected case

- **Spec:** "values not on a bin edge give identical tables."
- **What was run** (`evidence/F1/unaffected.R`, `unaffected.txt`):
  - Inputs built once with `set.seed(20261112)`: three cases of 400 rows, 332 of them found.
    - `consistency_wide`: Pcons 0.301-0.998 on a 0.001 grid, N_sg 21-420.
    - `consistency_default_grid`: Pcons 0.90-0.99 on a 0.01 grid without 0.95, N_sg 21-419.
    - `dina_mode`: Pcons NA, N_sg 20-420.
  - A `stopifnot()` asserts that no Pcons or N_sg lies on any break.
  - Trees: the start export (`/tmp/cran_prep_code_fixes_2026-11-12/start`, `76a62a71`) and the fixed copy `wt/F1`. Each ran in its own `timeout 120 Rscript --vanilla` process. All 11 return slots were compared.
- **Result:** `identical: TRUE` for all three cases, `Differing leaves: none`, and `OVERALL identical(start, fixed): TRUE`.
- **Wall:** start 2.88 s, fixed 2.87 s, driver 6.02 s.
- **Verifier re-run** (`verify/F1/unaffected/`), start export against `verify/F1/post`:
  - The inputs RDS is byte-identical to the implementer's (`cmp`).
  - All cases were identical. Without data.frame normalisation the result is also `raw identical (no normalisation): TRUE`.
  - Wall: workers 3.10 s and 3.09 s, driver 6.46 s.
- **Verifier grid** (`verify/F1/grid/compare.txt`): every Pcons on a 0.001 grid over [0, 1] and every N_sg in 0:450, one row per call.
  - All 994 off-edge Pcons values and all 445 off-edge N_sg values give identical bins on both trees.
  - Only the interior breaks change bin, each moving up one: Pcons 0.50, 0.70, 0.80, 0.85, 0.95 and N_sg 50, 100, 150, 200, 300.
  - Wall: 8.46 s (pre) and 8.65 s (post).

#### Existing tests

- No existing test file was run, because none reaches the changed lines (`evidence/F1/existing.txt`).
  - `grep` of `tests/testthat` for the function names, the slot names and the range columns found no file other than the new one (exit 1). I confirmed with `git grep` at `76a62a71` over `tests` and `inst` (exit 1).
  - The only caller in `R/` is `R/summarize_bootstrap_results.R:394`. `summarize_bootstrap_results()` has no caller in `R/`.
- Existing expectations updated: none.
- The full-suite run at `88302212` was still in progress when this section was written. Its log so far shows `✔ |          8 | bootstrap-summary-bins` (`suite/suite.log:90`).

#### What a user sees

`summarize_bootstrap_subgroups()` on five found rows: Pcons 0.95, 0.90, 0.85, 1.00, 0.949 and N_sg 50, 49, 100, 300, 299 (`evidence/F1/example.R`). `consistency_dist` from `evidence/F1/example.txt`:

```
before (start tree)                          after (fixed tree)
4:          0.8-0.85     1      20           4:          0.8-0.85     0       0
5:         0.85-0.95     3      60           5:         0.85-0.95     3      60
6:            >=0.95     1      20           6:            >=0.95     2      40
```

- Rows 1-3 are 0 in both tables.
- Before the fix, a consistency of exactly 0.95 was shown under "0.85-0.95". After the fix it is under ">=0.95".
- In the same file, `size_dist` moves from `2 1 0 0 2 0` to `1 1 1 0 1 1`: a subgroup of 50 moves from "<50" to "50-99".

#### Verification

- **Verdict:** CONFIRMED. Red, green and the unaffected case were reproduced. There was no repair or re-verify round: the `wf_results.json` entry has only `impl` and `verify`.
- **Problem 1 (minor).** The label-order expectations (test `:76-83`) pass on both trees. The verifier suggested a one-line comment saying so.
  - **Not acted on.** There is no `patches/F1.verified.patch`, and the commit equals `patches/F1.patch`. The comment at `:75` still reads "Labels in the documented order (the rows of each distribution table)."
  - Carried as an open item.
- **Problem 2 (minor).** `R/format_subgroup_summary_tables.R:222` highlights `">= 0.95"` (with a space), which never matches the label `">=0.95"`. The verifier said F1 must not fix it.
  - **Open item.** The line is unchanged at `61c15088:222`.
- **Recipe note (not a problem with the fix).** The verifier recorded that `git init && git apply --index` fails on a fresh export ("does not exist in index"). It succeeds after `git add -A` and a start commit.

#### Departures and interpretations

- **Departures:** none, per the implementer.
- **Extra test checks.** The spec's ten per-value bins form the first two blocks. The third block adds label order, counts and percents, all derived by hand. It does not change what the fix is.
- **One variable at a time.** Each per-value block holds the other variable at an interior point (N_sg 120, Pcons 0.92), so only one bin can move.
- **"Values not on a bin edge"** was read as synthetic tables whose Pcons and N_sg avoid every break, outer breaks included. Three cases were used, one of them DINA mode (Pcons NA, so only `size_dist` is built). Every return slot was compared, not only the two tables.
- **No roxygen or `man/` change.** The fix changes no documented argument or value, so `document()` was not run. The undocumented bin closure is left as an open item.
- **Existing-tests gate.** It was read as met with nothing to run, because no existing test reaches the changed lines.

#### paper_simulations

- No evaluated call reaches the changed lines (implementer; the verifier agreed).
- The article never names `summarize_bootstrap_subgroups()` or `summarize_bootstrap_results()`: `git grep summarize_bootstrap` at `76a62a71` exits 1. Its only `summarize_*` call is `summarize_simulation_results()`.

#### Open items

- **Top consistency bin never highlighted.** At `R/format_subgroup_summary_tables.R:222` (`76a62a71` and `61c15088`), the consistency gt table highlights rows where `` `Consistency Range` %in% c("0.85-0.95", ">= 0.95") ``. The label is `">=0.95"`, with no space, so the top bin is never highlighted green. After F1 more iterations land in that row. The bug predates F1 and is in a file outside its scope.
- **Bin closure not documented.** The roxygen at `R/summarize_bootstrap_subgroups.R:22-23` (`76a62a71` and `61c15088`) says only "Data.table of Pcons distribution by bins" and "Data.table of subgroup size distribution". It does not say the bins are left-closed. This was left out to keep the commit to the code change and its test.
- **Label-order checks pass before and after.** At `tests/testthat/test-bootstrap-summary-bins.R:76-83` (`0dc99947`) the label-order checks pass on both trees. They only fix the row order that the positional Count and Percent checks rely on. The verifier's suggested comment saying this was not added.
- **No test for the callers.** No test covers `summarize_bootstrap_results()` (its call at `R/summarize_bootstrap_results.R:394`) or `format_subgroup_summary_tables()` (`:186-191`). The gt consistency table they build changes for on-edge values but is untested.
- **Out-of-range values (noted, no action).** Pcons < 0 or > 1, or N_sg < 0, gives an NA bin before and after the fix. A proportion or a count cannot produce these values.
- **NEWS (resolved).** The implementer's suggested NEWS entry is in `88302212` (`NEWS.md:295-299`), so it is not open.

#### Corrections to the notes

- **Missing gate files.** `main/F1/` has no `patch_equal.txt` or `diffstat.txt`. The equality check and the `--stat` above come from my own `git show 0dc99947` runs.
- **Code-file line count.** The verifier gives the code file as "+2/-2". `git show --numstat 0dc99947` gives 4 insertions and 2 deletions: two lines gain a trailing comma and two `right = FALSE` lines are added.
- **Where n_found is defined.** The verifier cites n_found at `R/summarize_bootstrap_subgroups.R:107-115`. At `76a62a71` the `sg_found` selection is `:107-114` and `n_found <- nrow(sg_found)` is `:116`.

### F2 — CV fold assignments {1, NA} or {0, NA} no longer stop

**Status:** PASS. Commit `1c2565b3` (parent `0dc99947`, the F1 commit).

| File | Change |
|---|---|
| `R/forestsearch_cross_validation.R` | 10 lines (+7/-3) |
| `tests/testthat/test-cv-fold-na-assignments.R` | new, 78 lines |
| Total | 2 files, 85 insertions, 3 deletions |

No roxygen change, so no `man/` page. `main/F2/patch_equal.txt` is empty: the committed diff equals `patches/F2.patch` apart from index and hunk-header lines (re-checked here: same).

#### Site, before and after

`forestsearch_KfoldOut()` in `R/forestsearch_cross_validation.R`. The file is identical at `28e71ffa` (the spec's pointers) and at the start SHA `76a62a71`.

Unchanged, both trees: `n_sgfound <- length(unique(df_CV$treat.recommend))` at start `:1371`, head `:1401`.

Before, start `76a62a71`, `:1417-1420` (the stop is `:1420`):

```r
    if (n_sgfound == 2) {
      tabit <- with(df_CV, table(treat.recommend, treat.recommend.original))
      sensH <- tabit[1, 1] / sum(tabit[, 1])
      sensHc <- tabit[2, 2] / sum(tabit[, 2])
```

and `:1423-1425` (the duplicate table line `:1424` is removed):

```r
    } else {
      tabit <- with(df_CV, table(treat.recommend, treat.recommend.original))
      if (nrow(tabit) == 1L && identical(colnames(tabit), c("0", "1")) &&
```

After, head `61c15088`, `:1447-1454` (at `1c2565b3`: `:1417-1424`):

```r
    # Fold-by-original table (table() leaves out NA assignments).  Fold
    # assignments {1, NA} or {0, NA} count as two values in n_sgfound but
    # give a one-row table, so the two-level branch is taken only when the
    # table has two rows and two columns; other shapes use the single-level
    # branch.
    tabit <- with(df_CV, table(treat.recommend, treat.recommend.original))
    if (n_sgfound == 2 && nrow(tabit) == 2L && ncol(tabit) == 2L) {
      sensH <- tabit[1, 1] / sum(tabit[, 1])
```

and `:1458-1459` (at `1c2565b3`: `:1428-1429`):

```r
    } else {
      if (nrow(tabit) == 1L && identical(colnames(tabit), c("0", "1")) &&
```

- The single-level branch is otherwise unchanged: start `:1425-1445` = `1c2565b3` `:1429-1449` = head `:1459-1479`.
- At head, the code sits 30 lines lower than at `1c2565b3`. Two later fixes touched this file: F4 (`3612d0ba`) and F3 (`a4cce90a`). `R/` is identical at `61c15088` and `88302212`.

#### Red and green

| Run | Where | Result | Wall |
|---|---|---|---|
| Implementer red | Scratch copy of start `76a62a71` (`wt/F2`), test file only | `.1.2`: 2 passes (premise checks), 2 errors (test `:57:3` and `:70:3`), both "``Error in ``[.default`(tabit, 2, 2)`: subscript out of bounds``" at `R/forestsearch_cross_validation.R:1420:7` | 2.95 s |
| Implementer green | Same scratch copy with the fix | `......`: 6 passes, 0 failures, 0 errors | 2.97 s |
| Verifier red | `verify/F2/pre` (archive of `76a62a71` + test) | `[ FAIL 2 \| WARN 0 \| SKIP 0 \| PASS 2 ]`, same error at `:1420:7` | 3.19 s |
| Verifier green | `verify/F2/post` (archive + patch) | `[ FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 6 ]` | 3.08 s |
| Main red | Repository at `0dc99947` + test file only (this R file is identical to `76a62a71` there) | `.1.2`: 2 passes, 2 errors (test `:57:3`, `:70:3`), same message at `R/forestsearch_cross_validation.R:1420:7` | 3.0 s |
| Main green | Repository with the code applied | `......`: 6 passes, 0 failures | 2.9 s |

Sources: `evidence/F2/red.txt`, `green.txt`; `verify/F2/pre.redgreen.txt`, `post.redgreen.txt`; `main/F2/red.txt`, `green.txt`, `wall.txt`.

#### Expected values

Test file at `1c2565b3`. Data (lines 28-36): 100 subjects. Positions 1-8 are orig 0 and 9-20 orig 1, all with fold NA. Positions 21-47 (27) are orig 0 and 48-100 (53) orig 1, all with fold `lev`. `table()` drops the 20 NA, leaving 80.

| Line | Expected value | How derived |
|---|---|---|
| 55, 68 | `length(unique(.f2_rec(lev)))` is `2L` | By hand: one level and NA. Premise check on the data, not the code. Comment at line 54 (none of its own at 68). |
| 58, 71 | names `c("sens_H", "sens_Hc", "ppv_H", "ppv_Hc")` (line 26) | From the source: `names(sens_metrics_original) <-` at start `:1448` (per notes and verifier; no comment in the test). |
| 64 | `{1, NA}`: `c(0, 1, 0, 53 / 80)` | By hand, from the single-level rule (header lines 15-23; code comment at start `:1427-1431`). sens_H 0/27, sens_Hc 53/53, ppv_H 0/0 -> 0, ppv_Hc 53/80. Comments at lines 59-63. |
| 77 | `{0, NA}`: `c(1, 0, 27 / 80, 0)` | By hand, same rule. sens_H 27/27, sens_Hc 0/53, ppv_H 27/80, ppv_Hc 0/0 -> 0. Comments at lines 72-76. |

- The verifier re-derived both vectors with base-R arithmetic (`verify/F2/hand_arith.txt`: `TRUE`, `TRUE`). It flagged no value as lacking a comment.
- These six are all the expectations in the file.

#### Unaffected case

- **Spec:** assignments {0, 1} and {0, 1, NA} give identical values.
- **Synthetic run.** `forestsearch_KfoldOut()` on hand-built `res` with {0, 1} and {0, 1, NA}. It covers a one- and a two-factor full-data subgroup, `outall = FALSE` and `TRUE`, and `details = TRUE` with the console output captured: 8 cases.
- **Pipeline run.** `forestsearch()` + `forestsearch_Kfold(Kfolds = 5L)` on `.make_survival_data(N = 200L, HR_harm = 3.0, seed = 42L)`. One run on complete data (fold values 0, 1). One with 20 NA ages (fold values 0, 1, NA; 10 NA held out).
- **Result.** `identical()` TRUE for both parts. The implementer compared start vs scratch; the verifier re-ran start vs `post`. {0, 1, NA} keeps NA `sens_Hc`/`ppv_Hc` in both trees: synthetic `0.800000 NA 0.666667 NA`, pipeline `0.431579 NA 0.672131 NA`.
- **Verifier sweep.** 7 fold shapes x 6 original shapes x `outall` FALSE/TRUE = 84 cases. 56 were identical. 28 differ, and all 28 stopped on the start tree. 0 cases that ran before changed (`verify/F2/sweep_compare.txt`).

| Run | Start | Fixed |
|---|---|---|
| Implementer synthetic | 3.12 s | 3.03 s |
| Implementer pipeline | 8.47 s | 8.72 s |
| Verifier synthetic | 3.08 s | 3.08 s |
| Verifier pipeline | 8.33 s | 8.68 s |
| Verifier sweep | 4.24 s | 5.01 s |

#### Existing tests

These files call `forestsearch_KfoldOut()`, `forestsearch_Kfold(` or `forestsearch_tenfold(`. Each run was capped at 120 s. 0 failures everywhere.

| File | Implementer (fixed scratch) | Verifier (`post`) |
|---|---|---|
| test-cv-mixed-and-fatal.R | 58 pass, 17.8 s | 58 pass, 18.6 s |
| test-no-subgroup-cv.R | 72 pass, 29.6 s | 72 pass, 30.7 s |
| test-cv-no-subgroup-edges.R | 44 pass, 1 skip, 34.5 s | 44 pass, 1 skip, 35.4 s |
| test-cv-no-subgroup-dina-grf.R | 48 pass, 14.0 s | 48 pass, 14.7 s |
| test-cross-outcome-parity.R | 35 pass, 12 W, 2 skip, 45.3 s | 35 pass, 12 W, 2 skip, 47.3 s |
| test-downstream-consumers.R | 38 pass, 5 W, 22.9 s | 38 pass, 5 W, 23.0 s |
| test-rng-state-restored.R | 16 pass, 7.0 s | 16 pass, 7.3 s |
| test-input-validation.R | not run | 18 pass, 1 W, 4.5 s |
| test-mr-inference.R, whole file | hit the cap (exit 124), 120.0 s | not run whole |
| test-mr-inference.R, 4 chunks (21 blocks) | 108 / 15 / 12 / 13 pass; 33.2 / 40.0 / 75.4 / 56.9 s | 108 / 15 / 12 / 13 pass; 35.2 / 40.8 / 77.0 / 59.5 s |

- **Warnings and skips predate the fix.** The W and S patterns are identical on a start-SHA copy. Implementer: edges 34.5 s, parity 46.1 s, downstream 22.1 s. Verifier: parity 45.8 s, edges 33.7 s, downstream 22.0 s, input-validation 4.2 s. The skip reasons are "dev-load artifact" and "primary fit did not identify".
- **Existing expectations updated:** none. The commit touches no existing test file.

#### What a user sees

From `evidence/F2/example.txt`, held-out assignments {1, NA}:

```
# start SHA
> forestsearch_KfoldOut(res_1na)$sens_metrics_original   # folds {1, NA}
Error in `[.default`(tabit, 2, 2): subscript out of bounds

# fixed
> forestsearch_KfoldOut(res_1na)$sens_metrics_original   # folds {1, NA}
 sens_H sens_Hc   ppv_H  ppv_Hc
 0.0000  1.0000  0.0000  0.6625
```

Through the callers, inferred from the code and not measured:

- `forestsearch_Kfold()` now returns these values in `sens_summary`.
- `forestsearch_tenfold()` now completes, with these values in the affected simulation's `sens_out` row. Its `sens_summary` is the median across simulations (start `:1237`, head `:1267`).

#### Verification

Verdict: **CONFIRMED**. Red, green and the unaffected case were reproduced. There was no repair. The commit equals the gated patch, and there is no `patches/F2.verified.patch`, so the main session added nothing.

| Problem | Severity | Resolution |
|---|---|---|
| Open item 1 named only the 2x1 sibling. 1x1 tables (-> `1 NA 1 NA`, labelled H) and 1x2 tables with a zero cell (-> `NaN NA 0 NA`) also now reach the position-labelled fallback. | minor | Noted, no code change: Open item 1 below. All these shapes stopped before, and the spec routes them to the existing branch. |
| Open item 2: survival `outall = TRUE` "completes for {1, NA}" only when `treat.recommend.original` has no NA. With NA there, it stops ("incorrect number of dimensions") for every fold shape on both trees. | minor | Noted, no code change: Open item 2 below. Pre-existing and unrelated to F2. |
| The caller-level claim "both report the values above" was not measured. `forestsearch_tenfold()` reports a median across simulations. | minor | Reworded under "What a user sees" above. |

#### Departures and interpretations

- **The guard keeps `n_sgfound == 2`.** The spec says to take the two-level branch "only when the table has two rows and two columns". The commit adds `nrow(tabit) == 2L && ncol(tabit) == 2L` as a further necessary condition and keeps `n_sgfound == 2`. Reasons:
  - `n_sgfound` counts NA as a value: `length(unique(c(0, 1, NA)))` is 3, `length(unique(c(1, NA)))` is 2 (checked in base R). The spec forbids changing that count.
  - At the start SHA, {0, 1, NA} (n_sgfound 3) took the single-level branch. It ended in the fallback (start `:1439-1444`), which gives NA `sens_Hc`/`ppv_Hc`.
  - `table()` drops NA, so {0, 1, NA} gives a 2x2 table. A guard on the table shape alone would send it to the two-level branch and turn those NAs into values (code reading; not run). That breaks the unaffected case the spec names.
  - With both conditions: {0, 1} (n_sgfound 2, 2x2) and {0, 1, NA} (n_sgfound 3) take the same branch as before. {1, NA} and {0, NA} (n_sgfound 2, 1x2) now take the single-level branch.
  - This is the variant both earlier reports recommend: `dev/reports/REPORT_cran_prep_2026-11-11.md:387` and `..._part3_detail.md:366`, `:380`.
- **The duplicate table line is removed.** The `tabit <-` line inside `else` (start `:1424`) was dropped because the table is now built before the branch, as the spec asks.
- **The test has more than the four values.** Besides the four values per shape, it has two premise checks and two names checks.
- **`test-mr-inference.R` was run in four chunks** outside the tree, because the whole file exceeds the 120 s per-run cap. Its content was not changed.

#### paper_simulations

Not reached. `vignettes/articles/paper_simulations.qmd` has no match for `Kfold|tenfold|cross.valid|CV_sgs` (re-checked at `1c2565b3`). `forestsearch_KfoldOut()` is called only from `forestsearch_Kfold()` (start `:609`) and `forestsearch_tenfold()` (start `:1176`), and nothing in `R/` calls those two.

#### Open items

- **Non-2x2 shapes now reach a fallback that labels by position.** Site: `R/forestsearch_cross_validation.R`, head `:1473-1477` (start `:1439-1444`). It takes `tabit[1, 1]` whatever level the only row or column is. These shapes have n_sgfound 2 and stopped before:
  - 2x1, with the full data all Hc: `0.4 NA 1.0 NA` (`evidence/F2/example.txt`).
  - 1x1: `1 NA 1 NA`, labelled H.
  - 1x2 with a zero cell: `NaN NA 0 NA` (`verify/F2/sweep_compare.txt`).
  - Not fixed here: the spec routes them to the existing single-level branch, and the labelling predates F2. The part 3 report judges the 2x1 shape practically unreachable.
- **`outall = TRUE` table assembly still branches on `n_sgfound == 2`.** Head `:1622`, `:1669`, `:1729`, `:1762` (start `:1588`, `:1635`, `:1695`, `:1728`).
  - For {1, NA} with no NA in `treat.recommend.original`, the survival path now completes with an empty `KfA_0` row (`"0 (0.0%)"`, NA estimates).
  - With NA in `treat.recommend.original`, `outall = TRUE` stops with "incorrect number of dimensions" for every fold shape on both trees. That is pre-existing.
  - The GLM `outall = TRUE` path was not exercised with these shapes.
  - Not fixed here: the spec keeps `n_sgfound` as it is, and the stop is outside F2.
- **NA held-out assignments are dropped silently by `table()` in a direct `forestsearch_KfoldOut()` call** (head `:1452`). F3 (`a4cce90a`) added warnings only in `forestsearch_Kfold()` (head `:558`) and `forestsearch_tenfold()` (head `:1228`).
- **`tests/testthat/test-mr-inference.R` exceeds a 120 s per-run cap** as a whole file (120.0 s, exit 124). It passes in four chunks. A full-suite run needs it run without that cap.

#### Corrections to the notes

- Notes, `test-mr-inference.R` chunk B3: the wall is 56.9 s (`existing.txt`: 56.941 s), not 57.0 s.
- Verifier JSON: it puts the derivation comments at test "lines 61-66 and 74-78". At `1c2565b3` they are at lines 59-63 and 72-76, with the expectations at 64 and 77.
- Notes, open item 3: "F3c adds the warning". At head the warning is in the two callers (`:558`, `:1228`), not in `forestsearch_KfoldOut()`.

### F3 — Missing-covariate messages at three sites

**Status:** PASS. Commit `a4cce90a` (parent `685e0c5d`, the F14 commit).

| File | Change |
|---|---|
| `R/forestsearch_cross_validation.R` | 32 lines (+31/-1) |
| `R/forestsearch_helpers.R` | 15 lines (+15) |
| `R/run_simulation_analysis.R` | 50 lines (+50) |
| `tests/testthat/test-missing-membership-messages.R` | new, 363 lines |
| Total | 4 files, 459 insertions, 1 deletion |

- No roxygen change, so no `man/` page.
- `main/F3/patch_equal.txt` is empty: the committed diff equals `patches/F3.patch` apart from index and hunk-header lines.
- The three R files are the same at `61c15088` and `88302212` as at `a4cce90a` (`git diff --stat` empty). Line numbers at head below are at `61c15088`.

#### Site, before and after

**F3a, forest search extractor** (`.extract_fs_estimates_gen()`, `R/run_simulation_analysis.R`).

- Today's stop is unchanged: `out$sens <- if ((tp + fn) > 0) ...` at start `:1134`, `out$spec <- if ((tn + fp) > 0) ...` at start `:1135` (head `:1180`, `:1181`).
- Before, start `:1129-1131`:

```r
    hat_H  <- df$sg_hat == 1L

    tp <- sum(true_H & hat_H);   fp <- sum(!true_H & hat_H)
```

- After, head `:1172-1177` (the new stop is `:1175`):

```r
    hat_H  <- df$sg_hat == 1L
    # NA membership: tp and fn, or tn and fp, below would be NA and the first
    # `if` would stop.  Stop here with the reason (same cases stop).
    if (anyNA(hat_H)) .stop_na_membership(df$sg_hat, sg_def_txt, names(df))

    tp <- sum(true_H & hat_H);   fp <- sum(!true_H & hat_H)
```

- The definition at hand is kept in a local variable on both routes. Start `:1051` / head `:1092-1093` (df.est route):

```r
    df$sg_hat <- as.integer(fs_full$df.est$treat.recommend == 0)
    sg_def_txt <- fs_full$sg.harm  # definition at hand, for the NA message
```

- Start `:1058-1059` / head `:1100-1102` (M.* route):

```r
    }
    sg_def_txt <- sg_factors       # definition at hand, for the NA message
    if (length(sg_factors) > 0)
```

**F3a, GRF extractor** (`.extract_grf_estimates_gen()`, same file).

- Today's stop is unchanged: `if (out$size.H  > 10) out$hr.H.hat  <-` at start `:1456` (GLM) and `:1468` (survival); head `:1506`, `:1518`.
- Before, start `:1441-1444`:

```r
  harm_ind  <- as.integer(grf_data$treat.recommend == 0)
  grf_data$sg_hat <- harm_ind

  out$size.H  <- sum(harm_ind)
```

- After, head `:1487-1494` (the new stop is `:1491-1492`):

```r
  harm_ind  <- as.integer(grf_data$treat.recommend == 0)
  grf_data$sg_hat <- harm_ind
  # NA membership: size.H below would be NA and `if (out$size.H > 10)` would
  # stop.  Stop here with the reason (same cases stop).
  if (anyNA(harm_ind))
    .stop_na_membership(harm_ind, sg_harm_id, names(grf_data))

  out$size.H  <- sum(harm_ind)
```

**F3a, new internal helpers**, head `:917-955`, after `run_simulation_analysis()` ends (start `:914`). `.stop_na_membership()` is `:931-941`; `.sg_def_columns()` is `:947-955`.

```r
.stop_na_membership <- function(sg_hat, sg_def = NULL, cols = character(0)) {
  covs <- tryCatch(.sg_def_columns(sg_def, cols),
                   error = function(e) character(0))
  stop(sprintf(paste0("%d of %d subjects have no subgroup assignment: a ",
                      "covariate that defines the identified subgroup%s has ",
                      "missing values."),
               sum(is.na(sg_hat)), length(sg_hat),
               if (length(covs)) sprintf(" (%s)", paste(covs, collapse = ", "))
               else ""),
       call. = FALSE)
}
```

**F3b, GRF route** (`.forestsearch_grf_select()`, `R/forestsearch_helpers.R`).

- Site confirmed at start `:1990-1991`. The spec's pointer `:1987-1988` (at `28e71ffa`) is a blank line and a comment; the file is unchanged between `28e71ffa` and `76a62a71`.
- Before, start `:1990-1991`:

```r
  tr_est  <- .grf_evaluate_subgroup(sg_def, df)
  df.est  <- df; df.est$treat.recommend <- tr_est
```

- After, head `:1995-2006`. `tr_est` is still `:1990`; `:1991-1994` are a comment.

```r
  if (anyNA(tr_est)) {
    sg_vars <- unique(unlist(lapply(sg_def$conjunctions,
                                    function(cj) as.character(cj$variable)),
                             use.names = FALSE))
    warning(sprintf(paste0("%d of %d subjects have no subgroup assignment ",
                           "(treat.recommend is NA): a covariate that defines ",
                           "the GRF-selected subgroup (%s) has missing values."),
                    sum(is.na(tr_est)), length(tr_est),
                    paste(sg_vars, collapse = ", ")),
            call. = FALSE)
  }
  df.est  <- df; df.est$treat.recommend <- tr_est
```

**F3c, `forestsearch_Kfold()`** (`R/forestsearch_cross_validation.R`).

- The NA comes from `get_dfpred()` inside the fold workers (start `:508`). The fold loop runs under `suppressWarnings()` (start `:466`).
- Before, start `:552-555`:

```r
  resCV <- .collect_cv_results(resCV)

  # ===========================================================================
  # SECTION 7: POST-PROCESSING AND VALIDATION
```

- After, head `:558-564`. It follows `resCV <- .collect_cv_results(resCV)` (head `:552`) and a comment (`:554-557`).

```r
  n_na_cv <- sum(is.na(resCV$treat.recommend))
  if (n_na_cv > 0L)
    warning(sprintf(paste0("%d of %d held-out subjects have no subgroup ",
                           "assignment (treat.recommend is NA): a covariate ",
                           "that defines their fold's subgroup has missing ",
                           "values."), n_na_cv, nrow(resCV)),
            call. = FALSE)
```

**F3c, `forestsearch_tenfold()`** (same file).

- The NA comes from `get_dfpred()` at start `:1037`. The simulation loop runs under `suppressWarnings()` (start `:953`).
- Per-simulation list. Before, start `:1184-1185`:

```r
      mr_replicates         = if (isTRUE(mr_in_replicates)) mr_list else NULL
    )
```

- After, head `:1196-1201`:

```r
      mr_replicates         = if (isTRUE(mr_in_replicates)) mr_list else NULL,
      # NA held-out assignments, counted here and reported by the caller
      # (one warning per call, after collection); not part of the result.
      n_na_assign           = sum(is.na(resCV$treat.recommend)),
      n_assign              = nrow(resCV)
    )
```

- Warning in the calling session. Before, start `:1206-1209`:

```r
    warning(sims - length(valid_results), " simulations failed and were excluded.")
  }

  sens_out <- do.call(rbind, lapply(valid_results, `[[`, "sens_metrics_original"))
```

- After, head `:1228-1237`. A comment (`:1225-1227`) comes first, and `sens_out` follows at `:1239`.

```r
  n_na_cv <- sum(vapply(valid_results, `[[`, integer(1), "n_na_assign"))
  if (n_na_cv > 0L)
    warning(sprintf(paste0("%d of %d held-out assignments across %d ",
                           "simulations have no subgroup assignment ",
                           "(treat.recommend is NA): a covariate that defines ",
                           "the fold's subgroup has missing values."),
                    n_na_cv,
                    sum(vapply(valid_results, `[[`, integer(1), "n_assign")),
                    length(valid_results)),
            call. = FALSE)
```

#### Red and green

Old message in every red run: ``Error in `if ((tp + fn) > 0) tp/(tp + fn) else NA_real_`: missing value where TRUE/FALSE needed``. The GRF extractor gives the same text from `if (out$size.H > 10) ...`. In red, the F3b and F3c tests get no warning: `run$warnings` is empty and `run$has_call` is empty.

| Run | Where | Result | Wall |
|---|---|---|---|
| Implementer red | Scratch copy of start `76a62a71`, test file only | 9 not passed: 3 errors (F3a, test lines 103, 122, 135) and 6 failures (F3b 241, 246; F3c 322, 323, 340, 341); 42 passed | 7.64 s |
| Implementer green | Scratch tree `wt/F3` (start SHA + patch) | 10 tests, 60 expectations, 0 failed, 0 errors, 0 skipped | 7.50 s |
| Verifier red | `verify/F3/pre` (start SHA R/ + test file) | `FAIL 9 \| PASS 42`: 3 errors, 6 failures, the same lines as implementer red | 7.40 s |
| Verifier green | `verify/F3/post` | `FAIL 0 \| PASS 60` | 7.28 s |
| Verifier green, later base | `685e0c5d` + patch (`verify/F3/head`) | tests=10, nb=60, failed=0, err=0 | 7.14 s |
| Main red | Repository at `685e0c5d`, committed test file only | 9 not passed: 3 errors (103, 122, 135) and 6 failures (242, 247, 323, 324, 341, 342); 42 passed | 7.3 s |
| Main green | Repository after the code was applied (`a4cce90a`) | 60 passed, 0 failed | 7.3 s |

- The test lines for F3b and F3c differ by one between the implementer/verifier runs and the main runs. The main session added a comment line at test `:156` before committing (see Verification).
- Sources: `evidence/F3/red.txt`, `green.txt`, `green_counts.txt`; `verify/F3/pre_run.txt`, `post_run.txt`, `head_f3_green.txt`; `main/F3/red.txt`, `green.txt`, `wall.txt`.

#### Expected values

All expected values are worked out from hand-built inputs. Test lines are at `a4cce90a`.

| Expected value (test line) | How it was derived |
|---|---|
| F3a stop, df.est route: `"2 of 40 ... (x1) ..."` (shape A) and `"3 of 40 ... (x1) ..."` (shape B) (106) | By hand. Shape A puts NA on rows 9-10, shape B on rows 20-22 (comment `:65-66`, `:101`). `x1` is named because the definition is `"{x1 <= 0.5}"` (comment `:102`). The wording carries what the spec's F3a bullet requires (comment `:85-87`). |
| F3a stop with no definition at hand: count alone (111) | By hand: `fs_full` has no `sg.harm` (comment `:108`). |
| F3a stop, M.* route: `"... (grp) ..."`, 2 and 3 of 40 (124) | By hand: `create_subgroup_indicator()` reads `"grp.1"` as column `grp`, level `"1"`, and gives NA where `grp` is NA (comment `:116-118`). |
| F3a stop, GRF extractor: 2 and 3 of 40, `(x1)` (138) | By hand: the same shapes; `sg.harm.id = "x1 <= 0.5"` (comment `:134`). |
| `conditionCall(err)` is NULL (107, 139) | From the spec: `call. = FALSE` (comment `:107`; `:139` repeats the check with no comment of its own). |
| `any.H = 1L` (157) | By hand: a subgroup of 12 members is identified. **The verifier flagged this one as lacking a comment.** The main session added the comment at `:156` before the commit. |
| `size.H = 12`, `size.Hc = 28` (158-159) | By counting rows: members are rows 1-8 and 11-14 (comment `:45-51`). |
| sens 8/10, spec 26/30, ppv 8/12, npv 26/28 (146, 161) | By counting rows: tp 8, fn 2, fp 4, tn 26 (comment `:48-50`). |
| `sg_hat = as.integer(rec == 0)` (163) | By hand (comment `:162`). |
| No harm column, NA at rows 9-10: completes; `sg_hat` keeps the NA, `size.H = 12`, `size.Hc = 26` (178-180) | By hand: 40 - 12 - 2 = 26 (comment `:169-171`). |
| `flag_harm` NA, complete membership: R's own `"missing value where TRUE/FALSE needed"` (186-189) | Behaviour at the start SHA: a missing truth label is not a membership problem (comment `:182-183`). |
| F3b: exactly one warning, `"3 of 40 ... (x1, x2) ..."` (242-246) | By hand, three-valued logic on `{x1 <= 0.5} & {x2 > 0.2}`: rows 3, 7 and 11 are NA; rows 1 and 4 are FALSE (comment `:197-206`, `:240-241`). |
| F3b: `has_call` FALSE (247) | From the spec (comment `:247`). |
| F3b: returned `treat.recommend` is 0 at rows 15, 19, ..., 39, NA at rows 3, 7, 11, and 1 elsewhere (250-253); `grp.consistency$sg.harm.id` is the same vector (254-255); 3 NAs (257) | By hand (comment `:205-206`, `:249`, `:256`). |
| F3b, complete covariates: no warning (262); members rows 3, 7, ..., 39 (264-266) | By hand: odd rows with `x2 = 0.9` (comment `:263`). |
| F3c: the full-data fit has no NA (317) | By hand: the definition is on `x2` alone, which is complete (comment `:275-277`). |
| F3c Kfold: exactly one warning, `"3 of 40 held-out subjects ..."` (323); `has_call` FALSE (324) | By hand: `x1` is NA at rows 5, 12, 33; every subject is held out once (comment `:274-280`, `:321-322`, `:324`). |
| F3c Kfold: NA ids 5, 12, 33 (326); the non-NA assignments follow `x1` (330) | By hand (comment `:278-280`, `:327`). |
| F3c tenfold: exactly one warning, `"6 of 80 held-out assignments across 2 simulations ..."` (341); `has_call` FALSE (342) | By hand: 2 x 3 = 6 of 80 (comment `:281`). `:342` has no comment of its own; it repeats the check at `:324`. |
| F3c tenfold: the names of the result equal the 10 documented components (344-347); NA ids 5, 12, 33 in each `resCV_all` element (348-349) | From the documentation: the `@return` list of `forestsearch_tenfold()`, `R/forestsearch_cross_validation.R:726-777` at the start SHA (comment `:343`). The ids are by hand (comment `:280`). |
| F3c, complete `x1`: no warning from Kfold or tenfold, and no NA (357, 358, 362) | By hand, from the data construction (`:282-290`, `with_na = FALSE`). The block has no comment of its own. |

- The message wording is the implementer's design. It carries the content the spec requires (count, total, cause, covariates where known), and the verifier accepted this.
- The verifier re-derived every value by hand. It found no value that could only have come from running the fixed code.

#### Unaffected case

- **What the spec names:** "No returned value changes at any of the three sites". From the Tests bullet: F3a, a 0/1 membership completes with identical values; F3b and F3c, the returned object is identical to the start SHA's. F3 has no separate "Unaffected case" bullet.
- **What was run:** `evidence/F3/unaffected.R`.
  - Trees: the start tree (`git archive 76a62a71`) and the fixed scratch tree `wt/F3`, each in its own `Rscript --vanilla` with `pkgload::load_all()`, a sequential plan and `set.seed` before each call.
  - 15 cases. 11 complete: F3a FS df.est, FS M.* and GRF with 0/1 membership; FS with no harm column and NA membership; F3b with a supplied definition, with and without NA; F3b, a real binary GRF frontier fit on 500 rows; F3c Kfold and tenfold with mocked refits; a real survival fit; F3c, a real survival Kfold on 300 rows with NA in `age` and `bm`. 4 are the F3a stop shapes (A and B, FS and GRF).
- **Result:**
  - All 11 returned values are `identical()` TRUE after dropping timing leaves only (`unaffected.txt`, `unaffected_rawcheck.txt`). Without that, raw `identical()` is FALSE only where timing leaves differ: 5 cases in the implementer run, 6 in the verifier run.
  - All 4 stop cases stop on both trees. Start: `missing value where TRUE/FALSE needed`. Fixed: for example `2 of 40 subjects have no subgroup assignment: a covariate that defines the identified subgroup (x1) has missing values.`
- **Wall clocks:**
  - Implementer: 9.28 s start, 9.22 s fixed, 0.87 s compare.
  - Verifier re-run (`verify/F3/unaff_*`): 8.69 s pre, 8.84 s post. Same results.
- **Verifier extra checks:**
  - F3a property check, 240 random cases (`prop_f3a_cmp.txt`): `problems: 0`. 220 stops get the new message; 40 keep R's message, and these are the truth-label NA cases.
  - Default survival `forestsearch()` on gbsg: `timing-dropped identical: TRUE` (`default_path_cmp.txt`).
  - Survival GRF on NA covariates stops on both trees with `missing values and NaN's not allowed if 'na.rm' is FALSE`: 2.50 s pre, 2.52 s post (`surv_grf_na_out.txt`).

#### Existing tests

- **Implementer** (`evidence/F3/existing.txt`; fixed tree, each file in its own `timeout 120` process): 27 whole files.
  - The files: build_classification_table, nofind-classification, npv-empty-complement, fs_mr_oc_summary, directive-c, cv-no-subgroup-dina-grf, no-subgroup-bootstrap, no-subgroup-dina-grf-plot, no-subgroup-dina-grf-summary, no-subgroup-unified-contract, grf-inclusion-band, grf-survival-matrix, fs-family-report, cross-outcome-parity, cv-mixed-and-fatal, cv-no-subgroup-edges, no-subgroup-cv, downstream-consumers, rng-state-restored, betaHhat-contract, binary-default-or, dina-gaussian-floor, forestplot-xlog, fs-focus-tag, mr-admission-set, sg-focus-transparency and threshold-pair-directive-a.
  - Result: 244 tests, 4351 expectations, 0 failed, 0 errors, 9 skipped.
  - There were 17 uncaught warnings: 12 in cross-outcome-parity and 5 in downstream-consumers. Their texts (`existing_warnings.txt`) are candidate-pool truncation, `nb_boots < 100` and `Could not build metrics_table`. None is an F3 message.
  - `test-mr-inference.R` as one process hit the cap (exit 124, 120.05 s). It was then run in 5 parts (1-6, 7, 8, 9a, 9b): 21 tests, 148 expectations, 0 failed.
- **Verifier** (`verify/F3/existing/`): the same 27 files on its post tree gave the same totals, 244 tests, 4351 expectations, 0 failed. `test-mr-inference.R` ran in 8 parts: 21 tests, 148 expectations, 0 failed.
- **Verifier probe:** the verifier used a copy of the post tree that logs each F3 stop or warning. On it, it ran all 89 existing test files plus the new one, with `test-mr-inference.R` in 8 parts: 97 runs, 662 tests, 8494 expectations, 0 failed, 0 errors. F3 code fired only in the new file: 8 F3a stops, 1 F3b warning, 1 Kfold warning and 1 tenfold warning (`probe_out/missing-membership-messages.probe`).
- **Existing expectations updated: none.** The commit touches no test file other than the new one (`--stat`).

#### What a user sees

From `evidence/F3/example_console.txt`. The call is `forestsearch(subgroup_method = "grf", grf_selection = "frontier", outcome_type = "binary", ...)` on 500 rows with 30 NA in each of `x1`, `x2`, `x3` (`evidence/F3/example_console.R`). The verifier reproduced it exactly (`verify/F3/example_console_out.txt`).

Before (start SHA): no message, and 34 subjects silently have NA membership.

```
[1] "{x1 <= 0.42360689910128713}" "{x2 <= 0.77094318917952487}"

   0    1 <NA> 
 140  326   34 
```

After: the same returned values, plus one warning.

```
Warning message:
34 of 500 subjects have no subgroup assignment (treat.recommend is NA): a covariate that defines the GRF-selected subgroup (x1, x2) has missing values. 
[1] "{x1 <= 0.42360689910128713}" "{x2 <= 0.77094318917952487}"

   0    1 <NA> 
 140  326   34 
```

#### Verification

- **Verdict:** CONFIRMED. The verifier reproduced red, green and the unaffected case.
- **Problems** (3, all minor):

| Problem | Severity | Resolution |
|---|---|---|
| F3c gives the count only. The general bullet ("name the covariates where the function has the subgroup definition at hand") is not discussed as a departure. | minor | Noted. Recorded as departure 7 below; no code change. |
| `.mrct_analyze_trial()` (used by `mrct_trial_analysis()` and `mrct_case_study()`) records every `forestsearch()` warning in `$warnings`. With a supplied `grf_res` and NA covariates, the F3b text now appears there. | minor | Noted. Listed under Open items. |
| `expect_identical(out$any.H, 1L)` has no derivation comment of its own. | minor | Fixed by the main session before the commit, with a comment-only line. `patches/F3.verified.patch` and `patches/F3.patch` differ only by `+    # any.H = 1: a subgroup is identified (12 members, rows 1-8 and 11-14).` (test `:156`), and the test file goes from 362 to 363 lines. |

#### Departures and interpretations

1. **F3a names the covariates.** The F3a bullet asks only for the count and the cause. The general bullet asks to name the covariates where the definition is at hand, and both extractors have it (`fs_full$sg.harm`, the M.* factors, `sg.harm.id`). So the message adds `(x1, ...)` when the definition names data columns, and gives the count alone otherwise. The names are found inside `tryCatch`, so which cases stop never depends on them.
2. **GRF extractor: one check, not two.** The check sits at the first point where the old code stopped (start `:1441-1468`). No check was added in the classification block (start `:1536-1546`): NA membership cannot reach it, because `m_class$sg_hat` is copied from `harm_ind`.
3. **F3b applies to any outcome type.** The spec's "(binary, continuous or count)" was read as describing where the NA arises, not as a filter. Survival reaches this line with NA only through a supplied `grf_res`, because the survival GRF fit stops on NA covariates first.
4. **F3b warns for `df` only.** The spec names "the selected cut applied to every subject", which is the analysis frame. `df.predict` and `df.test` (head `:2007-2012`) do not warn.
5. **F3c tenfold counts inside each simulation.** The spec requires the warning in the calling session, not in a worker. So each simulation returns two internal counts (`n_na_assign`, `n_assign`), and the caller sums them and warns. The counts do not reach the result: the names test and `identical()` both show this.
6. **`forestsearch_KfoldOut()` does not warn.** It is exported, but `forestsearch_Kfold()` and `forestsearch_tenfold()` call it after collection (head `:621`, `:1188`). A warning there would break "one warning per exported call".
7. **F3c gives the count alone** (verifier problem 1). The F3c bullet asks for "the count". Each fold has its own definition, and `resCV` holds only the first two fold labels (`sg1`, `sg2`); this reason is the verifier's suggested wording.

#### paper_simulations

The evaluated chunks `run-alt-sims` (`:908`) and `run-null-sims` (`:966`) run the new F3a lines whenever a subgroup is found. The stop cannot fire there: the gbsg covariates have no NA, and `n_add_noise = 0L`. F3b and F3c are not reached (`run_fs_grf = FALSE`, no CV calls). Article output is unchanged.

#### Open items

- **Repeated F3b warnings in simulation loops.** The warning is at `R/forestsearch_helpers.R:1995-2005` (head). A loop that calls `forestsearch(subgroup_method = "grf")` on data with missing covariates repeats it once per replicate. This covers `run_simulation_analysis(methods = "grf-*")` with a GLM DGM, `fpr_calibration()` and `mrct_simulation()`; the survival MRCT functions reach it only with a supplied `grf_res`. Bootstrap and CV wrap `forestsearch()` in `suppressWarnings()`, so they do not repeat it. Not fixed: de-duplication across calls is outside the spec's "one warning per call".
- **`mrct_trial_analysis()` and `mrct_case_study()` `$warnings` can now hold the F3b text.** This happens when `fs_args` supplies `subgroup_method = "grf"` and a `grf_res`, and the trial has NA covariates. The capture is `.mrct_analyze_trial()`, `R/mrct_trial_helpers.R:45-55` (start and head). Not changed: this is the intended message, passed through as designed.
- **No message for NA in `df.predict` / `df.test` on the GRF route**, `R/forestsearch_helpers.R:2007-2012` (head; start `:1992-1997`). Not fixed: the spec names the analysis frame (departure 4).
- **No message from a direct `forestsearch_KfoldOut()` call** on held-out assignments with NA (`R/forestsearch_cross_validation.R:1341`, head). Not fixed: it would warn twice inside `forestsearch_Kfold()` (departure 6).
- **F3a stop not run end to end** through `run_simulation_analysis(methods = "grf-*")` with a GLM DGM that has missing covariates. It was shown by direct extractor calls on both trees, and the F3b warning was shown on that route's `forestsearch()` call. Not done: the rules allow no simulation campaign.
- **Covariate names are matched approximately.** `.sg_def_columns()` (`R/run_simulation_analysis.R:947-955`, head) matches identifier tokens to column names. It could, for example, name a column called `e` for a cut written `1e-05`. This affects only the message text, never which cases stop. Not fixed: the spec asks only that the covariates be named.
- **Grammar.** With a count of 1, the messages read "1 of N subjects have ..." (`R/run_simulation_analysis.R:934-936`, head, and the F3b and F3c messages). Not fixed: the tests fix the wording.
- **`tests/testthat/test-mr-inference.R` cannot run as one process within the 120 s cap.** It was covered in parts by both the implementer and the verifier. Not F3's file.
- **paper_simulations reaches F3a lines** (`R/run_simulation_analysis.R:1093` or `:1101`, `:1175`, `:1491-1492`, head). The stop cannot fire on gbsg. Listed because the lines are reached; nothing needs re-rendering.

#### Corrections to the notes

- The notes' diffstat (+458/-1, test file +362) and their red test lines for F3b/F3c (241, 246, 322, 323, 340, 341) describe the implementer's patch. The commit has +459/-1 and a 363-line test file, because of the main session's comment at `:156`. Those expectations are at 242, 247, 323, 324, 341, 342 (`main/F3/red.txt`).
- `patches/F3.patch`, which the notes cite, now holds the committed version. The verified version, which matches the implementer's diffstat, is `patches/F3.verified.patch`.
- Notes open item 7 (F2 and F4 share `R/forestsearch_cross_validation.R`) is closed. F3 was committed after F2 (`1c2565b3`) and F4 (`3612d0ba`), on `685e0c5d`, and the committed diff equals the gated patch.
- The verifier cites `R/mrct_trial_helpers.R:44-54`. Git shows the capture at `:45-52` and the message at `:53-55`.
- The verifier writes "all 97 test files ran". These were 97 runs: 89 existing files plus the new one, with `test-mr-inference.R` split into 8 parts.

### F4 — A-103: CV "Exact" means the same subgroup

**Status:** PASS. Commit `3612d0ba` (parent `1c2565b3`, the F2 commit).

| File (`git show --stat 3612d0ba`) | Change |
|---|---|
| `R/forestsearch_cross_validation.R` | +8 / -4 (`CV_sgs()`, one-factor branch only) |
| `tests/testthat/test-cv-exact-match.R` | new, +71 |
| Total | 2 files, 79 insertions, 4 deletions |

No roxygen, `man/`, `NAMESPACE`, `DESCRIPTION` or `NEWS.md` change in this commit.

#### Site, before and after

Before: `R/forestsearch_cross_validation.R:1818-1827` at `76a62a71`.

```r
  } else if (sg_depth == 1) {
    sg1a <- sg_analysis[1]

    # Exact match
    exact_match <- ifelse((sg1 == sg1a | sg2 == sg1a), 1, 0)
    exact_match[is.na(exact_match)] <- 0.0

    one_match <- exact_match
    cov1_match <- exact_match
    cov2_match <- NA
```

After: `R/forestsearch_cross_validation.R:1852-1865` at `61c15088`. The code is the same at `88302212`. In commit `3612d0ba` the block is at :1822-1835, because F2 adds 4 lines above it.

```r
  } else if (sg_depth == 1) {
    sg1a <- sg_analysis[1]

    # Factor appears in the fold's subgroup (either slot)
    contains_match <- ifelse((sg1 == sg1a | sg2 == sg1a), 1, 0)
    contains_match[is.na(contains_match)] <- 0.0

    # Exact match: the fold's subgroup is that factor alone (other slot NA)
    exact_match <- ifelse((sg1 == sg1a & is.na(sg2)) | (sg2 == sg1a & is.na(sg1)), 1, 0)
    exact_match[is.na(exact_match)] <- 0.0

    one_match <- contains_match
    cov1_match <- contains_match
    cov2_match <- NA
```

- `cov1_match` is returned as `cov1_exact`. The line is `cov1_exact = cov1_match,`: :1849 at the start SHA, :1887 at head.
- The depth-2 branch and the no-subgroup branch are unchanged.

#### Red and green

All runs use `devtools::test(filter = "cv-exact-match")`.

| Run | Where | Result | Wall |
|---|---|---|---|
| Implementer red | Scratch copy of `76a62a71` (`wt/F4`), test file only | 2 failures, 4 passes (`cv-exact-match: 1..2..`). `test-cv-exact-match.R:32:3`: `` `actual`: 1.0 1.0 1.0 0.0 0.0 `` vs `` `expected`: 1.0 0.0 0.0 0.0 0.0 ``. `:68:3`: `` `actual`: 0.60 `` vs `` `expected`: 0.20 `` | 2.89 s |
| Implementer green | Same scratch copy, after the fix | 6 passes, 0 failures (`cv-exact-match: ......`) | 2.78 s |
| Verifier red | `verify/F4/pre` (archive of `76a62a71` plus the test file) | FAIL 2, WARN 0, SKIP 0, PASS 4. Same two failures (:32:3, :68:3) | 2.93 s |
| Verifier green | `verify/F4/post` | FAIL 0, WARN 0, SKIP 0, PASS 6 | 2.90 s |
| Main-session red | Repository at `1c2565b3` with only the test file applied | 2 failures, 4 passes. Same two failures, same messages | 3.0 s |
| Main-session green | Repository with the `3612d0ba` code applied | 6 passes, 0 failures | 2.8 s |

Both red failures are the A-103 defect: the containment indicator was used as "Exact". The `one_match`, `cov1_exact`, "At least 1" and "Cov1 exact" expectations pass on both trees. This is by design, because the spec keeps those values.

#### Expected values

These are the six expectations in `tests/testthat/test-cv-exact-match.R` at `3612d0ba`. The folds (`sg1`, `sg2`) are :21-22: (A, NA), (A, B), (B, A), (B, NA), (NA, NA), with A = `{age <= 50}`, B = `{sex >= 1}` and the full-data subgroup {A}.

| Line | Expectation | Expected | How derived |
|---|---|---|---|
| :32 | `out$exact_match` | `c(1, 0, 0, 0, 0)` | By hand from the spec F4 definition: 1 only when the fold's subgroup is {A} alone (comment :28-31). The vector is the one given in the spec's Test line. |
| :36 | `out$one_match` | `c(1, 1, 1, 0, 0)` | By hand: "factor appears in the fold's subgroup", in either slot; none scores 0 (comment :34-35). This is the spec's vector. |
| :37 | `out$cov1_exact` | `c(1, 1, 1, 0, 0)` | Same containment indicator and same comment (:34-35). This is the spec's vector. |
| :68 | `ko$find_metrics["Exact"]` | `1 / 5` | The mean over the 5 folds of the per-fold vector at :32 (comment :65-67). |
| :69 | `ko$find_metrics["At least 1"]` | `3 / 5` | The mean of `c(1, 1, 1, 0, 0)` (comment :65-67). |
| :70 | `ko$find_metrics["Cov1 exact"]` | `3 / 5` | The mean of `c(1, 1, 1, 0, 0)` (comment :65-67). |

- The data frame at :43-62 is constructed input, not an expected value.
- The verifier flagged no missing comment. It re-derived all six values by hand.
- The verifier also checked that the fixed `exact_match` equals the package's canonical-equality measure (`.cv_canonicalize_sg`) on all 16 (sg1, sg2) pairs over {A, B, C, NA}. It did this for three one-factor analyses (`verify/F4/own/own_compare.txt`).

#### Unaffected case

- **Spec:** a two-factor full-data subgroup, and no full-data subgroup, give identical output.
- **Script:** `evidence/F4/unaffected.R`. Each tree ran in its own `Rscript --vanilla` process with `pkgload::load_all()`, a sequential plan and `set.seed()` before every call.

| Object(s) | What was run |
|---|---|
| `cvsgs_two_AB`, `cvsgs_two_BA`, `cvsgs_none_NULL`, `cvsgs_none_chr0` | `CV_sgs()` on all 16 (sg1, sg2) pairs over {A, B, C, NA}, with `sg_analysis` set to `c(A, B)`, `c(B, A)`, `NULL` and `character(0)` |
| `fs2`, `kf2`, `ko2_f`, `ko2_t`, `tf2` | Survival `forestsearch()` with N = 400, HR_harm 3, seed 42. `sg.harm` is `{age <= 63} & !{stage}`; the two-factor premise is checked with `stopifnot`. Then `forestsearch_Kfold(Kfolds = 5)`, `forestsearch_KfoldOut(outall = FALSE / TRUE)` and `forestsearch_tenfold(sims = 2, Kfolds = 3)` |
| `fs0`, `kf0`, `ko0_f`, `ko0_t`, `tf0` | Null survival fit with N = 160, HR_harm 1, seed 7. `sg.harm` is NULL (premise checked with `stopifnot`). Then the same four CV calls |

- **Comparison:** `identical()` after dropping six timing leaves only (`time_search`, `minutes_all`, `timing_seconds`, `timing_minutes`, `tmins_search`, `tmins_iteration`). data.tables were compared as data.frames.
- **Result:** all 14 objects are `identical: TRUE`, and the output reports `ALL IDENTICAL: TRUE`.
  - Implementer (`evidence/F4/unaffected.txt`): start 13.00 s, fixed 12.73 s, compare 0.17 s.
  - Verifier (`verify/F4/unaffected_verify.txt`): pre 12.89 s, post 12.90 s, compare 0.19 s. Also all TRUE.
- **Comparator sanity check** (`evidence/F4/unaffected_sanity.txt`): in a mutated copy, `ko2_f` and `tf2` were reported as `identical: FALSE`, and the output reported `ALL IDENTICAL: FALSE`.
- **Verifier's extra check** (`own_check.R`; 2.48 s pre, 2.47 s post): the two-factor analyses (4 orderings) and the no-subgroup analyses are identical. For one-factor analyses, only `exact_match` differs.

#### Existing tests

All runs were on the fixed tree with `timeout 120` per run. Each cell gives tests/expectations/failed/skipped, then the wall clock.

| File | Implementer (`evidence/F4/existing.txt`) | Verifier (`verify/F4/existing/existing_summary.txt`) |
|---|---|---|
| test-cross-outcome-parity.R | 9/49/0/2, 44.62 s | 9/49/0/2, 45.13 s |
| test-cv-mixed-and-fatal.R | 7/58/0/0, 17.94 s | 7/58/0/0, 17.15 s |
| test-cv-no-subgroup-dina-grf.R | 4/48/0/0, 13.99 s | 4/48/0/0, 13.67 s |
| test-cv-no-subgroup-edges.R | 6/45/0/1, 34.23 s | 6/45/0/1, 33.16 s |
| test-downstream-consumers.R | 11/43/0/0, 22.31 s | 11/43/0/0, 21.88 s |
| test-input-validation.R | 12/19/0/0, 4.04 s | 12/19/0/0, 4.25 s |
| test-no-subgroup-cv.R | 13/72/0/0, 29.76 s | 13/72/0/0, 28.70 s |
| test-rng-state-restored.R | 7/16/0/0, 7.31 s | 7/16/0/0, 6.49 s |
| test-mr-inference.R, whole file | exit 124 at 120.06 s (cap). No failure had been reported when it stopped | not run |
| test-mr-inference.R, one block per process | 21 blocks, 148 expectations, 0 failed. 3.46-67.46 s per block | 6 blocks that reach Kfold or tenfold, 17 expectations, 0 failed. 8.18-31.92 s per block |

- No run had an error.
- The skips have environmental causes, not this fix:
  - cross-outcome-parity, 2 skips: "primary fit did not identify; test N/A".
  - cv-no-subgroup-edges, 1 skip: the multisession worker could not resolve package internals, which is a dev-load artifact.
- **Existing expectations updated:** none.
- The only existing expectation on an Exact value is `test-cv-mixed-and-fatal.R:135-136` at `76a62a71`. It uses a two-factor `sg_analysis` (:58) and passes unchanged.

#### What a user sees

The example is a synthetic survival fit (N = 300, HR_harm 3, seed 31) with full-data `sg.harm: {age <= 55}`. It is run through `forestsearch_Kfold(Kfolds = 5)` and then `forestsearch_KfoldOut(details = TRUE)`. The fold subgroups are the same on both trees:

```
 cvindex         sg1    sg2
       1 {age <= 55} !{sex}
       2 {age <= 63}   <NA>
       3 {age <= 63}   <NA>
       4 {age <= 64}   <NA>
       5 {age <= 55}   <NA>
```

`diff evidence/F4/user_example_start.txt evidence/F4/user_example_fixed.txt` (apart from the tree path and the wall clock):

```
13c13
< Exact match: 0.4
---
> Exact match: 0.2
```

- Before the fix, fold 1 (`{age <= 55} & !{sex}`) counted as exact, giving 2/5. After the fix only fold 5 counts, giving 1/5.
- `At least 1 match: 0.4` and `Cov 1 exact: 0.4` are unchanged.
- The printed `find_metrics` row changes the same way (line 22: Exact `0.4` -> `0.2`).
- The verifier reproduced the example (`verify/F4/user_example_pre.txt` and `user_example_post.txt`).
- Wall clocks: implementer 6.01 s / 6.22 s; verifier 5.91 s / 6.17 s.

#### Verification

**Verdict:** CONFIRMED. The verifier reproduced red, green and the unaffected case. There was no repair round.

`patches/F4.verified.patch` does not exist, so the main session added nothing. `main/F4/patch_equal.txt` is empty. Re-checked here: the committed diff equals `patches/F4.patch`, apart from the index and hunk-header lines. The hunk header reads -1822 instead of -1818 because F2 was committed first.

| Problem (verifier) | Severity | Resolution |
|---|---|---|
| `vignettes/articles/forestsearch-2024.Rmd` chunk `kfold-cv` (`eval = TRUE`) reaches `CV_sgs()`. Its rendered "Exact Match" may fall at the next render | minor | Left as an open item. `vignettes/` is read-only and is not re-rendered |
| The `CV_sgs()` roxygen `@return` and the "Exact match on all factors" label do not state the new meaning | minor | Left as an open item. The NEWS entry states the new meaning (`88302212:NEWS.md:317-322`) |
| The second clause `(sg2 == sg1a & is.na(sg1))`, a fold (NA, A), has no test | minor | Optional and not added; the committed test has five folds. The verifier's own check gives `exact_match` 1 for (NA, A) (`own_compare.txt`, row 4) |

#### Departures and interpretations

- **"The other slot is NA":** read as either slot, so (NA, A) is also exact (the second clause). This is the audit's proposed form. Package-produced folds never give (NA, A), because `sg1` is NA only when nothing was found. The clause matters only for direct callers of the exported `CV_sgs()`.
- **A fold (A, A):** now scores `exact_match` 0, where it scored 1 before. This follows the same wording, since the other slot is not NA. It agrees with canonical equality (`own_compare.txt`, row 1). Package folds do not produce this case.
- **Field names:** the spec names `cov1_match`, but the returned field is `cov1_exact`. The test checks `cov1_exact`, as the spec's Test line does.
- **Extra test block:** the test adds a second block beyond the spec's per-fold check. It calls the exported `forestsearch_KfoldOut()`, and its expected values follow arithmetically from the spec's vectors.
- **Per-block run:** `test-mr-inference.R` was run one block per process, because the whole file exceeds the 120 s per-run cap.
- **Site lines:** the spec's pointer is "near 1818-1826" (at `28e71ffa`). At the start SHA the block is :1818-1827.

#### paper_simulations

Does not reach the changed lines.
- `CV_sgs()` is called only by `forestsearch_KfoldOut()` (:1386 at `76a62a71`).
- `forestsearch_KfoldOut()` is called only by `forestsearch_Kfold()` (:609) and `forestsearch_tenfold()` (:1176).
- Nothing in `R/` calls those two.
- `paper_simulations.qmd` names none of the four functions (0 hits at `76a62a71`) and sources no script.

#### Open items

- **Vignette that reaches the change:** `vignettes/articles/forestsearch-2024.Rmd:598-634` (chunk `kfold-cv`, `eval = TRUE`; same at head) calls `forestsearch_tenfold()` (:606) and `cv_metrics_tables()` (:622), which reach `CV_sgs()`. If that article's `sg.harm` has one factor, its displayed "Exact Match" may fall at the next render. Not re-rendered here because `vignettes/` is read-only.
- **Wording of "Exact":** the `CV_sgs()` roxygen `@return` "exact match" (`R/forestsearch_cross_validation.R:1779-1780` at start, :1813-1814 at head) is unchanged. So are the labels "Exact Match" and "Exact match on all factors" (`R/cv_summary_tables.R:413` and :418, both SHAs). The audit marks this wording change as optional, and the label is in a file outside F4's scope.
- **maxk >= 3:** `CV_sgs()` sees only `sg1` and `sg2`. The fold refits store `fs.train$sg.harm[1]` and `[2]` only (`R/forestsearch_cross_validation.R:505-506` and :1034-1035 at start; :505-506 and :1046-1047 at head). A third factor is dropped, so the depth-2 "Exact" can over-count. Fixing this needs an upstream `sg3` column, which is outside this fix.
- **NaN in "Cov 1 & 2" (pre-existing, not changed):** in the one-factor branch, `cov2_any <- NA` (:1831 at start, :1869 at head). As a result, `find_metrics["Cov 1 & 2"]` is `NaN` when every fold has `cov1_any = 1`; the user example shows this. F4 does not affect it.
- **Untested (NA, A) clause:** optional per the verifier. Package folds cannot produce this case.

#### Corrections to the notes

- **maxk lines:** the notes cite the `sg.harm[1:2]` storage at ":503-504, :1026-1027". Git at `76a62a71` puts the `sg1`/`sg2` assignments at :505-506 and :1034-1035.
- **Vignette chunk:** the verifier gives "about lines 598-635". In git the chunk runs from :598 to :634.
- **Existing Exact expectation:** the verifier cites `test-cv-mixed-and-fatal.R:132-135`. In git the comment is at :132 and the `expect_equal` is at :135-136.
- **"After" line numbers:** the notes' "after" lines (:1818-1831) are for the scratch tree built from the start SHA. In the repository the block is at :1822-1835 in `3612d0ba` and at :1852-1865 at head `61c15088`.
- **NEWS:** the notes leave the NEWS entry to a later commit. It is present at `88302212:NEWS.md:317-322`.

### F5 — A-151: test_hte_from_forestsearch() includes continuous covariates

**Status:** PASS, after one repair and a second verification. Commit `4aaab2f0`.

| File | Change (`git show --stat 4aaab2f0`) |
|---|---|
| `R/test_hte_crump.R` | +22 / -2 |
| `tests/testthat/test-hte-crump-confounders.R` | new, +191 |
| Total | 2 files, 213 insertions, 2 deletions |

- No roxygen, `man/`, `NAMESPACE`, `DESCRIPTION` or `NEWS.md` change.
- Both files are unchanged from `4aaab2f0` to `61c15088` and `88302212` (`git diff` is empty).

#### Site, before and after

Before, at the start SHA `76a62a71`:

```r
# R/test_hte_crump.R:397-400
  if (!is.null(fs_result$confounders.evaluated) &&
      length(fs_result$confounders.evaluated) > 0) {
    covs <- fs_result$confounders.evaluated
    source_desc <- "confounders.evaluated (GRF/LASSO selected)"
```

```r
# R/test_hte_crump.R:409-411 (unchanged by F5; :425-427 at head)
  for (v in covs) {
    col <- df[[v]]
    if (is.null(col)) next
```

```r
# R/test_hte_crump.R:416
  if (length(X_cols) == 0) return(NULL)
```

- A cut label such as `"age <= 55"` is not a column of `df`, so `df[[v]]` is NULL and the label is skipped with no message (:411).
- If nothing remains, :416 returns NULL silently. The no-field branch at :406 warns `"No covariates found."`.

After, at head `61c15088` (same text at `4aaab2f0`):

```r
# R/test_hte_crump.R:397-416
  if (!is.null(fs_result$confounders.evaluated) &&
      length(fs_result$confounders.evaluated) > 0) {
    # Evaluated labels are bare column names (categorical) or cut
    # expressions such as "age <= 55" (continuous): map each label to the
    # df columns it refers to, and use the unique columns.  A label that is
    # exactly a column name is that column (it may not be a syntactic R
    # name, e.g. "ER-pos"); any other label is parsed.
    covs <- character(0)
    for (lab in fs_result$confounders.evaluated) {
      lab_vars <- if (lab %in% names(df)) lab else
        intersect(tryCatch(all.vars(parse(text = lab)),
                           error = function(e) character(0)),
                  names(df))
      if (length(lab_vars) == 0)
        warning("Label '", lab, "' in confounders.evaluated cannot be ",
                "resolved to a column of df; skipped.")
      covs <- c(covs, lab_vars)
    }
    covs <- unique(covs)
    source_desc <- "confounders.evaluated (GRF/LASSO selected)"
```

```r
# R/test_hte_crump.R:432-436
  if (length(X_cols) == 0) {
    warning("No covariate found in df for the ForestSearch covariates; ",
            "returning NULL.")
    return(NULL)
  }
```

- The `source` text is unchanged (start :400, head :416).
- The `confounders.candidate` branch (head :417-420) and the column loop (head :424-431) are unchanged.

#### Red and green

All runs use `devtools::test(filter = "hte-crump-confounders")` with `NOT_CRAN=true`. The test file has 5 blocks and 35 expectations (9, 8, 5, 7, 6).

| Run | Where | Result | Wall |
|---|---|---|---|
| Implementer red, first version | start tree with only the first test file (27 expectations) | 16 of 27 failed, 0 errors. Summary line `.....12.3...456789abcd.ef.g` (`scratch_F5/evidence_prior/red.txt`) | 6.06 s |
| Implementer green, first version | `wt/F5` (scratch copy of `76a62a71`) | 27 of 27 pass | 5.60 s |
| Verifier red / green, round 1 | `verify/F5/pre`, `post` (first version) | pre: 16 of 27 failed; post: 27 of 27 pass. Recorded in `wf_results.json` only; the round-1 folder was removed before round 2 | 6.18 s / 5.61 s |
| Implementer red, after repair | `scratch_F5/red_start`: fresh `76a62a71` export plus only the repaired test file | 19 of 35 failed, 0 errors (`red.txt`, `red_counts_start_tree.txt`). Per block: 3, 5, 5, 5, 1. Key messages: `` `actual`:   "bm" `` vs `` `expected`: "age" "bmi" "bm" `` (:71); `Expected is.null(r$value) to be FALSE.` (:98); `Expected r$warnings to have length 2. Actual length: 0.` (:119, :138); `Expected r2$warnings to have length 1. Actual length: 0.` (:150); block 5 mixed part: `` `actual`:   "bm status" "ER-pos" `` (:190) | 6.15 s (counts 6.07 s) |
| Implementer green, after repair | `wt/F5` (scratch copy of `76a62a71` with the fix) | 35 of 35 pass; failed 0, error 0, skipped 0, warning 0 (`green.txt`, `green_counts.txt`) | 5.44 s (counts 5.47 s) |
| Implementer green, fresh export | fresh `76a62a71` export, staged and committed, regenerated patch applied | 35 of 35 pass (`green_fresh_patch.txt`) | 5.76 s |
| Repaired tests on the prior fix | fresh export with the first-version `R/test_hte_crump.R` | 5 failures, all in block 5, e.g. `` `actual`:   "ER"        "pos" `` vs `` `expected`: "bm status" "ER-pos" `` (`repaired_tests_on_prior_fix.txt`) | 5.68 s |
| Verifier red / green, round 2 | `verify/F5/pre`, `post` | pre: 19 failures, 0 errors (`pre_test.txt`); post: 35 of 35 pass (`post_test.txt`) | 6.13 s / 5.23 s |
| Main-session red | repository at the preceding commit `3612d0ba` plus only the test file | Summary line `.....12.3...456789abcd.efg.hi.....j`: 35 results, 16 pass, 19 fail. Only the first 10 are printed (max 10 reports), all `Failure`, at lines 71, 73, 81, 98, 99, 102, 103, 105, 119, 120 (`main/F5/red.txt`) | 5.9 s |
| Main-session green | repository after the code was applied | 35 of 35 pass, then `DONE` (`main/F5/green.txt`) | 5.4 s |

- `R/test_hte_crump.R` is identical at `76a62a71` and `3612d0ba` (`git diff` is empty), so the main-session red runs the same code as the start tree.
- `main/F5/patch_equal.txt` is empty: the committed diff equals the gated patch `patches/F5.patch` (239 lines, ASCII).

#### Expected values

Test file `tests/testthat/test-hte-crump-confounders.R` at `4aaab2f0`.

| Expected value (line) | How it was derived |
|---|---|
| Block 1 preconditions: labels start with `"age <= "` and `"bmi <= "`; `age`, `bmi` are not labels; `bm` is (:57-60) | get_fsdata() docs and code: a confounder with 4 or more distinct values gets `var <= value` cut labels (process_conf_force_expr); a binary one keeps its bare name. Comment :53-56. |
| Block 1: 0 warnings (:65) | By hand: every label refers to a column (age, bmi, bm). Comment :64. |
| Block 1: `covariates_used == c("age", "bmi", "bm")` (:71) | The spec ("the confounder names passed to `forestsearch()`"). Order from get_fsdata(): `confs = c(conf.cont_cuts, conf.categorical)`, cuts in confounders.name order. Comment :66-70. |
| Block 1: `n_covariates == 3L` (:73) | Three covariates; `ncol()` returns an integer. Comment :72. |
| Block 1: `source == "confounders.evaluated (GRF/LASSO selected)"` (:75) | The start-SHA string; the spec says do not change it. Comment :74. |
| Block 1: `hte_test` equals `test_hte(df$time, df$treat, cbind(age, bmi, bm))` (:81) | The spec ("equals `test_hte()` called directly on those columns"); X built as a numeric matrix, as `test_hte()` documents. Comment :77-78. |
| Block 2 preconditions: `length(ev) > 0`; no label is a column name (:92-93) | Same get_fsdata() rule. Comment :90-91. |
| Block 2: 0 warnings; non-NULL value; class `hte_test` (:97-99) | Spec: "a fit with only continuous confounders returns a test". **No comment** on these three lines (verifier round 2, minor; still open, see below). |
| Block 2: `covariates_used == c("age", "bmi")`; `n_covariates == 2L` (:102-103) | By hand: the two names passed, in the order passed. Comment :100-101. |
| Block 2: `hte_test` equals `test_hte()` on `cbind(age, bmi)` (:105) | The spec, as in block 1. The block-1 comment covers the construction. |
| Block 3, labels `c("age <= 55", "zzz <= 1", "bmi <=", "bm")`: 2 warnings, one naming `"zzz <= 1"`, one naming `"bmi <="`; `covariates_used == c("age", "bm")`; equals `test_hte()` on `cbind(age, bm)` (:119-124) | By hand: `zzz` is not a column; `"bmi <="` does not parse. The spec asks for a warning naming each unresolved label. Comment :116-118. |
| Block 4, label `"zzz <= 1"` alone: NULL; exactly 2 warnings; one names the label; one contains `"No covariate found"` (:137-140) | From the code path: one label warning, then one empty-`X_cols` warning. Comment :130-134. |
| Block 4, `confounders.candidate = c("q1", "q2")`: NULL; exactly 1 warning, containing `"No covariate found"` (:149-151) | From the code path: the candidate branch has no per-label warning; the column loop skips q1, q2; the empty-`X_cols` branch warns once. Comment :142-146. |
| Block 5, labels `c("bm status", "ER-pos")` on a df that also has `ER` and `pos`: 0 warnings; `covariates_used == c("bm status", "ER-pos")`; `n_covariates == 2L`; equals `test_hte()` on those two columns (:176-181) | By hand: each label is exactly a column name, so it is that column, as at the start SHA (the spec's unaffected case). Comment :167-172. |
| Block 5, labels `c("bm status", "age <= 55", "ER-pos")`: 0 warnings; `covariates_used == c("bm status", "age", "ER-pos")` (:189-190) | By hand: two column-name labels are themselves; `"age <= 55"` is parsed to `age`; unique columns in label order. Comment :183-185. |

- No expected value was copied from the fixed code's output. No seed, tolerance or expected value was changed after a run (notes).
- Block 5 data columns are deterministic functions of existing columns; no extra RNG draws (:160-165).
- Both verifiers re-derived every value by hand and agreed.

#### Unaffected case

- **Spec:** "a fit whose evaluated labels are all column names gives identical output."
- **What was run** (`evidence/F5/unaffected.R`): each tree in its own `Rscript --vanilla`, `load_all`, sequential plan, `set.seed` before each call. `forestsearch()` with no LASSO or GRF, then `test_hte_from_forestsearch()` with OLS on `time` and logistic on `event`. Compared: evaluated labels, both results with their warnings, and the fit with timing leaves dropped.
  - s1: syntactic names `bm`, `sex`, `grp`. Labels `bm | sex | grp`.
  - s2: non-syntactic names `"bm status"`, `"ER-pos"`, `grp`. Labels `bm status | ER-pos | grp`. The script asserts all labels are column names and at least one is non-syntactic.
  - s2x: the s2 fit on a df with extra columns `ER` and `pos`.

| Run | Trees | Result | Wall |
|---|---|---|---|
| Implementer, after repair (`unaffected.txt`) | `start` vs `wt/F5` | s1 and s2 (incl. s2x) all `identical()` TRUE; `UNAFFECTED CASE IDENTICAL (all members): TRUE`; 0 warnings everywhere | 4.61 s / 4.48 s / compare 0.18 s |
| Implementer, start vs prior fix (`unaffected_prior_fix_compare.txt`) | `start` vs first version | s1 TRUE; s2 `h_ols`, `h_log`, `hx_ols`, `hx_log` FALSE; overall FALSE. Shows the extended case catches the blocking problem | prior-fix run 4.46 s / compare 0.19 s |
| Verifier round 2, re-run (`verify/F5/unaffected_out.txt`) | `start` vs `verify/F5/post` | all TRUE; overall TRUE | 4.32 s / 4.46 s / compare 0.19 s |
| Verifier round 2, own check (`verify/F5/my_unaffected_out.txt`) | same | 6 minimal result lists (numeric, factor, character, non-syntactic and duplicated column-name labels; candidate branch; no field) x OLS, logistic, Poisson, poly_order 2: 24 of 24 TRUE, `ALL IDENTICAL: TRUE` | 2.51 s / 2.43 s |
| Verifier round 1 (`wf_results.json`) | `start` vs first version | s1 only: all TRUE. Its non-syntactic probe was NOT identical (the blocking problem) | 3.00 s / 3.06 s / 0.16 s |
| Cumulative, at head (`cumulative/F5/unaffected.txt`) | `start` vs `head_fixes` (`61c15088`) | all TRUE; overall TRUE | 4.37 s / 4.44 s / compare 0.19 s |

#### Existing tests

- None were run. `grep` of `tests/testthat` (excluding the new file) for `test_hte_from_forestsearch|test_hte\b|test_hte\(|hte_test|crump|Crump` finds nothing (`existing.txt`).
- `git grep -n test_hte_from_forestsearch 76a62a71 -- R inst tests` finds only the definition, `R/test_hte_crump.R:388`. I re-ran it and got the same.
- Existing expectations updated: none.

#### What a user sees

From `evidence/F5/example.txt`: `forestsearch()` with `age`, `bmi`, `bm`, no LASSO or GRF. `fit$confounders.evaluated` is `age <= 55 | age <= 47 | age <= 63 | bmi <= 27 | bmi <= 25 | bmi <= 30 | bm` on both trees.

Before (start):
```
covariates_used: bm
    Covariates: bm
2. Constant Cond. Ave. TE: chi-sq=  0.11 (df=1) p=0.7367 | z= -0.63 p=0.7347
N0=160, N1=140, K=2
--- only continuous confounders (age, bmi) ---
is.null(result): TRUE
```

After (F5):
```
covariates_used: age bmi bm
    Covariates: age, bmi, bm
2. Constant Cond. Ave. TE: chi-sq=  1.31 (df=3) p=0.7269 | z= -0.69 p=0.7550
N0=160, N1=140, K=4
--- only continuous confounders (age, bmi) ---
is.null(result): FALSE
covariates_used: age bmi
```

- New warnings, where there was silence before (`example_warnings.txt`): `Label 'zzz <= 1' in confounders.evaluated cannot be resolved to a column of df; skipped.` and `No covariate found in df for the ForestSearch covariates; returning NULL.`
- After the repair the example re-ran with identical output apart from timings (`example_repair_rerun.txt`); verifier round 2 also reproduced it (`verify/F5/example_out.txt`).

#### Verification

- **Round 1: PROBLEMS.**

| # | Problem | Severity | Resolution |
|---|---|---|---|
| 1 | The first version parsed every label. A non-syntactic categorical name stored as a bare label was dropped (`"bm status"`) or mapped to unrelated columns (`"ER-pos"` became `ER`, `pos`, with 0 warnings when those columns exist). This broke the spec's unaffected case. | blocking | Fixed in the repair. Code: a label that is exactly a column name is that column (head `R/test_hte_crump.R:406`). Test block 5 added (:154-191). `unaffected.R` gained members s2 and s2x. |
| 2 | Block 4 checked only `>= 1` warning for the candidate case, and never checked that the second warning in the first part is the no-covariate one. | minor | Fixed in the repair: `expect_length(r2$warnings, 1)` (:150) and `grepl("No covariate found")` checks (:140, :151), with code-path comments. |
| 3 | The harness recipe `git init && git apply --index` fails on an unstaged export. | minor (recipe only) | No patch change. On a staged and committed export, `git apply --check`, `--index` and `-R --check` all pass. |

- **Round 2: CONFIRMED.** RED, GREEN and the unaffected case reproduced. A `roxygenise` (rd) run on a scratch copy left every `man/` page unchanged (`verify/F5/doc_check_out.txt`, 11.17 s).

| # | Problem | Severity | Resolution |
|---|---|---|---|
| 1 | Block 2 lines :97-99 (0 warnings, non-NULL, class `hte_test`) have no derivation comment. The values are right. | minor | **Left open.** No comment was added before commit: there is no `patches/F5.verified.patch`, `patch_equal.txt` is empty, and the committed file still has no comment at :97-99. |
| 2 | An open item in the notes cites the round-1 verifier's `'my age'` probe, whose output was deleted with the round-1 folder. | minor | **Noted, not fixed.** An equivalent probe survives at `verify/F5/probe_nonsyn_cont.R` and `probe_nonsyn_cont_out.txt`. See "Corrections to the notes". |

#### Departures and interpretations

- **Column-name labels are resolved before parsing (repair).** The spec's method is `all.vars()` of the parsed label, intersected with `names(df)`. Applied literally, it breaks the spec's own unaffected case for non-syntactic column names, which `forestsearch()` accepts for categorical confounders. A label that is exactly a column name refers to that column, so it is mapped directly. Every other label follows the spec's method. Verifier round 2 judged this acceptable.
- **Mapping only on the `confounders.evaluated` branch.** The spec names only that field. The `confounders.candidate` fallback still looks names up as columns. It now shares the new no-covariate warning, which the spec states generally.
- **One warning per unresolved label, naming it.** A label that fails to parse and one whose variables are not columns get the same message.
- **Warning wording is the implementer's**; the spec gives none. `call.` is left at the default, as in the existing `warning("No covariates found.")` (head :422).
- **"The covariates used are the confounder names passed"** is tested with order: block 1 expects `c("age", "bmi", "bm")` exactly, with continuous confounders listed first on purpose so that get_fsdata()'s cut-then-categorical order matches the order passed.
- **Three extra test blocks.** Blocks 3 and 4 test the spec's two warning clauses. Block 5 tests the column-name member of the unaffected case. The spec's "Test" bullet names only blocks 1 and 2.
- **Roxygen not changed.** The spec does not ask for it.

#### paper_simulations

- Not reached. `vignettes/articles/paper_simulations.qmd` has no call to `test_hte_from_forestsearch()`, `test_hte()` or any Crump function, and nothing in `R/` or `inst/` calls `test_hte_from_forestsearch()`.

#### Open items

- `R/test_hte_crump.R:386` (start and head): the roxygen `@return` still reads "List with hte_test, covariates_used, source." It does not describe how labels map to columns or the two new warnings. Not fixed: the spec asks for no roxygen change.
- `R/test_hte_crump.R:416` (head; :400 at start): the `source` text "confounders.evaluated (GRF/LASSO selected)" is wrong on fits that use neither GRF nor LASSO. Not fixed: the spec says keep it.
- `R/test_hte_crump.R:417-420` (head; :401-404 at start): the `confounders.candidate` fallback holds `q1`, `q2`, ... on the consistency path, which are not `df` columns. That branch now ends in the new warning plus NULL instead of a silent NULL. Not changed: outside the spec's mapping, and normal fits do not reach it.
- `tests/testthat/test-hte-crump-confounders.R:97-99` (head): no derivation comment on three expected values (verifier round 2, minor). Not fixed before commit.
- `quarto/applications/validation_hte_tests_crump.qmd:1366` (chunk `fs-informed-test`, :1360) calls the function when an `fs` object exists; :1313 shows a call in a plain code block. The presentations copy has the same lines. Read-only and not re-rendered; if the chunk runs, its printed covariates would now include continuous ones.
- A cut label on a continuous confounder with a non-syntactic name (e.g. `"my age <= 55"`) would be skipped with a warning. Not reachable at `76a62a71`: `forestsearch()` emits no cut label for `"my age"` (`verify/F5/probe_nonsyn_cont_out.txt`). Not handled here because no such label exists today.
- Observation outside F5: a 2-level factor confounder (`sex`, levels F/M) did not appear in `confounders.evaluated` at the start SHA either. Not investigated.

#### Corrections to the notes

- **paper_simulations grep.** The notes say a case-insensitive `hte` grep matches only "HTEs" at lines 97, 1022, 1028, 1077. `git grep -n -i hte 76a62a71 -- vignettes/articles/paper_simulations.qmd` also matches :117 (`weightedsurv`) and :777 (`plot_weighted_km`), as substrings. The verdict (not reached) is unchanged.
- **Quarto caller lines.** The notes put :1313 and :1366 both in chunk `fs-informed-test`. Per `git show`, :1313 is in a plain ```` ``` ```` block opened at :1308 (not evaluated); only :1366 is in the chunk at :1360.
- **`'my age'` probe.** The notes say the round-1 probe showed "output is identical on both trees". That file no longer exists. The surviving probe also has a continuous `bmi`: on both trees `evaluated` is `bmi <= 27 | bmi <= 25 | bmi <= 30 | bm`, so `"my age"` gets no label. But `covariates_used` is `bm` at start and `bmi bm` after the fix. The difference is the intended A-151 change, not a regression.

### F6 — A-042: LASSO, the intercept is not a selected variable

**Status:** PASS. Commit `f5b65516` ("fix: LASSO intercept is not a selected variable (A-042)"), parent `4aaab2f0`.

| File (from `git show --stat f5b65516`) | Lines |
|---|---|
| `R/forestsearch_main.R` | 9 (+8/-1) |
| `R/get_FSdata_helpers.R` | 20 (+15/-5) |
| `man/lasso_selection.Rd` | 10 (+8/-2) |
| `tests/testthat/test-lasso-intercept.R` (new) | 214 |
| **4 files** | 245 insertions, 8 deletions |

- The files are unchanged between `76a62a71` and the parent `4aaab2f0` (`git diff --stat` is empty).
- `get_FSdata()` (`R/get_fsdata.R`) is not edited. Its existing rule at `R/get_fsdata.R:354` (start), `if (length(lassokeep) > 0)`, now applies to every outcome type.

#### Site, before and after

`R/get_FSdata_helpers.R`, `lasso_selection()`, the selection:

- Start `76a62a71`, lines 333-338:
  ```r
  coefs <- as.vector(coef(fit))
  names(coefs) <- rownames(coef(fit))

  # Determine selected and omitted variables
  selected <- names(coefs)[coefs != 0]
  omitted <- names(coefs)[coefs == 0]
  ```
- Head `61c15088`, lines 339-348:
  ```r
  coefs <- as.vector(coef(fit))
  names(coefs) <- rownames(coef(fit))

  # Determine selected and omitted variables from the confounder rows only.
  # For the GLM families coef(fit) also has an "(Intercept)" row; the
  # intercept is not a candidate variable, so it is kept in `coefficients`
  # but never counted as selected or omitted.
  conf_coefs <- coefs[names(coefs) %in% confounders.name]
  selected <- names(conf_coefs)[conf_coefs != 0]
  omitted <- names(conf_coefs)[conf_coefs == 0]
  ```

`R/get_FSdata_helpers.R`, the roxygen:

| Item | Start `76a62a71` | Head `61c15088` |
|---|---|---|
| Title | :214 `#' LASSO selection for Cox model` | :214 `#' LASSO selection for survival and GLM outcomes` |
| Description | :216 `#' Performs LASSO variable selection using Cox regression.` | :216-219 `#' Performs LASSO variable selection with \code{glmnet}: a Cox model for` ... `#' (Poisson also for a binary outcome when \code{offset.name} is given).` |
| `@return` | :227 `#' @return List with selected, omitted variables, coefficients, lambda, and fits.` | :230 unchanged, plus :231-233 `#'   \code{selected} and \code{omitted} contain confounders only; the` / `#'   \code{"(Intercept)"} of a GLM fit is in neither, but is kept in` / `#'   \code{coefficients}.` |

`man/lasso_selection.Rd` carries the same title, description and value text.

`R/forestsearch_main.R`, the LASSO summary label:

- Start `76a62a71`, lines 3394-3397:
  ```r
  if (details && use_lasso) {
    n_kept <- length(lassokeep)
    n_total <- n_kept + length(lassoomit)
    cat("Cox-LASSO selected:", n_kept, "of", n_total, "candidate factors\n")
  ```
- Head `61c15088`, lines 3393-3403. At `f5b65516` this was 3394-3404; later commits `685e0c5d` and `61c15088` moved it up one line.
  ```r
  if (details && use_lasso) {
    n_kept <- length(lassokeep)
    n_total <- n_kept + length(lassoomit)
    lasso_label <- switch(outcome_type,
      binary     = "Logistic-LASSO",
      continuous = "Linear-LASSO",
      count      = "Poisson-LASSO",
      "Cox-LASSO"
    )
    cat(paste0(lasso_label, " selected:"), n_kept, "of", n_total,
        "candidate factors\n")
  ```
- This is the same switch as `R/get_fsdata.R:343-348` (the same lines at start and head).

#### Red and green

All runs are `devtools::test(filter = "lasso-intercept", reporter = "summary")` with `NOT_CRAN=true`. The progress line has one character per expectation: `.` is a pass, and a digit or letter is a failure or error.

| Run | Where | Result | Wall |
|---|---|---|---|
| Implementer red | `wt/F6` at `76a62a71`, test file only | `.1..23..4..56..7..89......................abc.d.e..`: 37 pass, 14 not (13 failures, 1 error), in tests 1, 3 and 4 | 5.78 s (`red.txt`); a first run with the default 10-report cap took 5.93 s (`red_first10.txt`) |
| Implementer green | `wt/F6` with `F6.patch` | 62 pass, 0 failures, 0 errors | 5.97 s |
| Verifier red | fresh `76a62a71` export (`pre`) + test file | the same progress line, 37 pass / 14 not. `lasso_intercept_counts_pre.txt`: `tests=4 expectations=50 failed=13 errors=1`. Identical to the implementer's `red.txt` except the WALL line | 5.83 s |
| Verifier green | export + `F6.patch` (`post`) | 62 pass. `lasso_intercept_counts_post.txt`: `expectations=62 failed=0 errors=0` | 5.93 s |
| Main red | repository at `4aaab2f0` + only the `tests/` hunks | the same progress line: 37 pass / 14 not. The body lists the first 10 (9 failures + 1 error), then `Maximum number of 10 failures reached` | 5.6 s |
| Main green | repository + the whole patch | `..............................................................`: 62 pass, 0 failures | 5.8 s |

Key failure messages (from `evidence/F6/red.txt`, verbatim, shortened):

- Test 1, at :70 for each outcome: ``Expected `"(Intercept)" %in% res$selected` to be FALSE.``
- Test 1, :82, continuous: `` `actual`:   "(Intercept)" "age" "biomarker" "biomarker_hi" ``
- Test 3, :159 (error): `get_FSdata: no candidate cut expressions could be constructed from the supplied inputs.  Upstream source counts: GRF cuts (0), DINA cuts (0), LASSO selected (1), ...`
- Test 4, :208 for binary: ``Expected `summary_line` to have length 1.`` / `Actual length: 0.` (forestsearch stopped in `get_FSdata`, so no summary line was printed)
- Test 4, :209 for continuous: `x | Cox-LASSO selected: 4 of 5 candidate factors`. For count: `x | Cox-LASSO selected: 2 of 5 candidate factors`.

#### Expected values

Line numbers are those of `tests/testthat/test-lasso-intercept.R` at `f5b65516`. Tests 1-3 loop over binary, continuous and count. Test 4 also covers survival.

| Expected value (line) | How it was derived (the test's comment) |
|---|---|
| `names(res$coefficients) == c("(Intercept)", .f6_confs)` (:65) | `?glmnet`: `intercept` defaults to TRUE for gaussian, binomial and poisson. `coef()` returns `(Intercept)` first, then one row per column of `x`, in the order supplied (:60-64). |
| `"(Intercept)"` not in `selected` (:70) and not in `omitted` (:71) | The F6 spec: the intercept is not a confounder (:68-69). |
| `selected` and `omitted` do not overlap (:76), and together they are the 4 confounders (:77) | The roxygen ("selected, omitted variables") and the spec (:73-75). |
| `selected == .f6_confs[conf_coef != 0]` (:82), `omitted == .f6_confs[conf_coef == 0]` (:83) | Definition of LASSO selection at `lambda.min` (:79-80). |
| Mock test: `lassokeep` has length 0 (:116); `lassoomit == .f6_confs` (:117) | The mock selects nothing, and `get_FSdata()` passes its output through (:89-91, :115). |
| Mock test, full pool: the cut variables are the 4 confounders (:125); `biomarker_hi` and `sex` appear verbatim (:126); some `age <= ` cut (:127); some `biomarker <= ` cut (:128) | `get_FSdata()` cut construction: 2-level factors are categorical (`cont.cutoff = 4`); continuous variables get the default cuts mean, median, qlow, qhigh from `cut_var()` (:119-124). |
| Mock test: `confs` with LASSO on equals `confs` with LASSO off (:134) | The rule "if LASSO selects nothing, keep the full pool"; with `cut_type = "default"` both branches build default cuts (:130-133). |
| Real fit with `coef()` mocked: `lassokeep` has length 0 (:169); `lassoomit == .f6_confs` (:170) | The mock sets every covariate row to 0; the order is the row order of `coef()` (:166-168). |
| Real fit, full pool (:174, :175) | "by the same hand derivation as in the previous test" (:172-173). |
| Label `Logistic-` / `Linear-` / `Poisson-` / `Cox-LASSO` (:209-212) | The labels `get_FSdata()` prints (the switch in `R/get_fsdata.R`) (:181-183). The spec says to mirror that switch. |
| `of 4 candidate factors` and `selected: [0-4]` (:211) | `length(lassokeep) + length(lassoomit)` = the 4 confounders supplied; the intercept is not a candidate factor (:183-185). `[0-4]` is a bound, not a value. |
| Exactly one summary line (:208) | No comment of its own; the header comment at :181-185 covers the line's label and count. The verifier did not flag it. |

- The verifier said: "Every expected value carries a comment saying where it comes from. None was copied from the fixed code's output." It flagged no expected value for a missing comment.
- The verifier also noted that some assertions hold on both trees and are invariant locks, not red evidence (problem 2 below): :65, :71, all of test 2, and the survival case of test 4.

#### Unaffected case

- **Spec:** "the survival path, and a non-survival fit where LASSO selects at least one covariate, give identical candidate cuts."
- **Implementer run** (`evidence/F6/unaffected.R`):
  - One `Rscript --vanilla` per tree, with `pkgload::load_all()`, a sequential plan and `set.seed(2026)` before each call.
  - Data:
    - survival: `.make_survival_data(200, 42)`, confounders age, stage, sex, noise;
    - binary: `.make_binary_data(300, 42)`, with y redrawn with an age effect;
    - continuous: `.make_continuous_data(200, 42)`;
    - count: `.make_count_data(200, 42)`, offset `ftime`.
  - On both trees LASSO kept: survival `stage`; binary `age`; continuous `age biomarker biomarker_hi`; count `age`.
  - Calls: `lasso_selection()`, `get_FSdata(use_lasso = TRUE)` and `forestsearch(use_lasso = TRUE, use_grf = FALSE)` for each type. For survival also `use_grf = TRUE`, and the lines printed by `details = TRUE`.
- **Result** (`evidence/F6/unaffected.txt`):
  - 952 leaves on each side; 6 differ.
  - All 6 are the intended change: `"(Intercept)"` leaves `$binary|continuous|count$lasso_selection$selected` and `$...$get_FSdata$lassokeep`.
  - `identical()` TRUE:
    - every survival leaf, including the printed `Cox-LASSO selected: 1 of 4 candidate factors`;
    - for binary, continuous and count: the candidate cuts (`get_FSdata` `df`, `confs_names`, `confs`), the whole `forestsearch` object without timing leaves, and `coefficients`.
  - `ALL GATES: TRUE`.
  - Wall: start tree 6.31 s, fixed tree 6.35 s (process), compare 0.95 s.
- **Verifier re-run** of the same script against `start` and `post`: the same 952/6 result, `ALL GATES: TRUE`. Wall 6.22 s / 6.30 s / 1.01 s.
- **Verifier's own check** (`verify/F6/my_unaffected.R`):
  - New data: n = 400, 5 confounders, LASSO selecting in all four types.
  - Four `get_FSdata()` configurations (default, default with `collapse_cuts = FALSE`, median, `defaultcut_names` with `conf.cont_medians_force`), plus `survival::gbsg`.
  - Result: candidate cuts, `df` and `lassoomit` identical everywhere; `lassokeep` differs only by the intercept; survival and gbsg are fully identical; `ALL: TRUE`.
  - Wall 5.31 s pre, 5.16 s post (process).

#### Existing tests

No existing expectation was updated. The commit touches only the new file under `tests/`.

| File | Tests | Failed | Errors | Skipped | Warnings | Wall, implementer (fixed) | Wall, verifier (post) |
|---|---|---|---|---|---|---|---|
| test-get-fsdata-diagnostics.R | 19 | 0 | 0 | 0 | 4 | 7.54 s | 7.69 s |
| test-glm-pipeline-integration.R | 4 | 0 | 0 | 0 | 0 | 6.22 s | 6.02 s |
| test-input-validation.R | 12 | 0 | 0 | 0 | 1 | 4.23 s | 4.22 s |
| test-return-shape-contract.R | 10 | 0 | 0 | 0 | 9 | 11.48 s | 10.89 s |
| test-downstream-consumers.R | 11 | 0 | 0 | 0 | 5 | 22.75 s | 22.09 s |
| test-fs-family-report.R | 5 | 0 | 0 | 0 | 0 | 5.14 s | 5.55 s |
| test-cross-outcome-parity.R | 9 | 0 | 0 | 2 | 12 | 45.68 s | 46.53 s |
| test-plan-restored.R (verifier only) | 8 | 0 | 0 | 0 | 0 | - | 20.10 s |

- The verifier's pre-fix run gives the same counts (`verify/F6/existing_pre.txt`).
- The 2 skips are the test's own `primary fit did not identify`.
- The only existing tests with LASSO on are survival tests, at `test-get-fsdata-diagnostics.R:453` and `:477` (start).
- `test-fs-family-report.R:123` passes `use_lasso = TRUE` through `.fr_driver_args()` to `fs_family_report()`. There it is a report flag; no LASSO fit runs.

#### What a user sees

Binary outcome with no covariate signal: `forestsearch(outcome_type = "binary", use_lasso = TRUE, use_grf = FALSE, details = TRUE)` on `.make_binary_data(N = 200, seed = 42)`. Quoted from `evidence/F6/example.txt`:

- **Before:**
  ```
  > lasso_selection(...)$selected
  [1] "(Intercept)"
  > console lines mentioning LASSO selected / Omitted:

  > warnings:
  Error in get_FSdata: get_FSdata: no candidate cut expressions could be constructed from the supplied inputs.  Upstream source counts: GRF cuts (0), DINA cuts (0
  FSdata failure - returning NULL result
  > fs$error_log$stage: get_FSdata
  > number of candidate factors evaluated: 0
  ```
- **After:**
  ```
  > lasso_selection(...)$selected
  character(0)
  > console lines mentioning LASSO selected / Omitted:
  Logistic-LASSO selected: 0 of 4 candidate factors
    Omitted: age, biomarker, biomarker_hi, sex
  > fs$error_log$stage: NULL
  > number of candidate factors evaluated: 8
  ```
- When LASSO keeps a covariate, only the label and count change (continuous, `verify/F6/continuous_label.txt`): `Cox-LASSO selected: 4 of 5 candidate factors` becomes `Linear-LASSO selected: 3 of 4 candidate factors`.

#### Verification

**Verdict: CONFIRMED.** The verifier reproduced red, green and the unaffected case.

- `tools::checkRd` gives `parse ok`.
- `roxygen2` regenerates `man/lasso_selection.Rd` byte-identically (11.78 s).
- No repair or re-verify was run: the F6 entry in `wf_results.json` has only `impl` and `verify`.
- There is no `patches/F6.verified.patch`.
- `main/F6/patch_equal.txt` is empty: the committed diff equals the gated patch. I re-ran that diff against `git diff f5b65516^ f5b65516` and it is also equal.

| # | Severity | Problem | Resolution |
|---|---|---|---|
| 1 | minor | The spec's test "`lasso_selection()` mocked to select nothing, `get_FSdata()` keeps the full pool" (test 2, :88-136) cannot go red. The start tree's guard at `get_fsdata.R:354` was already correct. | Noted, not changed. Test 2 is kept as a lock (Departure 1). The red for the spec's second bullet comes from test 3. The suggested comment above test 2 ("expected to pass on the pre-fix tree") was not added. |
| 2 | minor | :65, :71 and the survival case of test 4 pass on both trees. They lock invariants rather than the defect. | Noted; no code change was needed. The optional comment at :65 and :71 was not added. |
| 3 | minor | The docs go beyond the literal spec: the description and the 3 new `@return` lines (`R/get_FSdata_helpers.R:231-233`). | Accepted as the documented Departure 3. The lines are kept; they describe only the changed contract. |

#### Departures and interpretations

1. **The spec's mocked test passes on the start tree too.**
   - The defect is in `lasso_selection()`, not in `get_FSdata()`'s guard. A mock that replaces `lasso_selection()` cannot show it.
   - Test 2 is kept as specified, as a lock on the rule.
   - The added test 3 runs the real `lasso_selection()` and mocks only `coef()` of the glmnet fit. It gives the red for the downstream symptom.
2. **Added test 4** for the printed label. The label is in the spec's change list, so it gets its own red and green.
3. **The roxygen description and `@return` were corrected as well as the title.** The description repeated the Cox-only claim ("using Cox regression"). `@return` now states the changed contract.
4. **"The existing rule in `get_FSdata()` then applies"** was read as: no edit to `get_FSdata()`. It is unchanged.
5. **"Identical candidate cuts"** was read as: `get_FSdata` `df`, `confs_names` and `confs` are identical, and so is the whole `forestsearch` object. `selected` and `lassokeep` are allowed to lose `"(Intercept)"`, since that is the fix.
6. **The label mirrors the existing switch exactly,** including how it treats a binary outcome with an offset (open item 1).

#### paper_simulations

The implementer's statement: it **reaches the changed lines; the output is unchanged**.

- `run-alt-sims` and `run-null-sims` call `run_simulation_analysis()` with `fs_params$use_lasso = TRUE` on survival DGMs. So `lasso_selection(outcome_type = "survival")` runs the new lines in every simulation, and the label switch runs for sim 1.
- A Cox fit has no intercept row, so `selected`, `omitted` and the printed `Cox-LASSO` line are unchanged.

#### Open items

1. **Binary outcome with an offset.**
   - With `offset.name`, `lasso_selection()` fits `family = "poisson"`: `R/get_FSdata_helpers.R:288-291` at start, `:294-297` at head.
   - Both labels still print `Logistic-LASSO`: `R/get_fsdata.R:343-348`, and `R/forestsearch_main.R:3396-3401` at head.
   - Not fixed here: the spec said to mirror the existing switch.
2. **Comment omits count.** `R/get_fsdata.R:327` (start and head) says "(Cox for survival, GLM for binary/continuous)". Cosmetic and outside F6's files.
3. **`run_simulation_analysis()` roxygen.**
   - `R/run_simulation_analysis.R:341-348` (start and head) says GLM parameters route "not to `\code{forestsearch()}` itself".
   - The A-042 audit notes that the code does pass `outcome_type`.
   - This is a documentation side issue outside F6's files.
4. **Test 3 depends on an import.** It mocks forestsearch's imported `coef` binding (`importFrom(stats, ...)`, `NAMESPACE:497` at start). If that import were removed, the mock would error at setup instead of passing.
5. **Lower counts for GLM outcomes.**
   - The `LASSO selected (%d)` count in `get_FSdata()`'s no-candidates diagnostic is one lower: `R/get_fsdata.R:582-583` at start, `:588-589` at head.
   - `n_kept`/`n_total` in the `forestsearch()` summary are also one lower.
   - This is intended. NEWS covers it at `88302212`, `NEWS.md:332-339` ("counts no longer include the intercept").
6. **Duplicated confounder names.** The verifier noted a theoretical corner case, not raised as a problem: with a duplicated `confounders.name`, `[.data.frame` renames the copy (e.g. `a.1`). That name would now be in neither `selected` nor `omitted`. Not tested.

#### Corrections to the notes

- The notes' "after" line numbers for `R/forestsearch_main.R` (3397-3404) are at `f5b65516`. At head `61c15088` the site is 3396-3403. The `R/get_FSdata_helpers.R` lines are the same at both.
- The notes' open item 7 ("NEWS.md and DESCRIPTION are not edited") is superseded: commit `88302212` adds the A-042 entry at `NEWS.md:332-339`.
- The main-session `red.txt` lists only 10 of the 14 failures and errors, because of testthat's default report cap. Its progress line shows all 14, the same as the implementer's and verifier's.

### F7 — A-107: fpr_calibration() outcome-type fallback

**Status:** PASS. Commit `65831b54` (parent `f5b65516`, start SHA `76a62a71`).

| file (`git show --stat 65831b54`) | lines changed |
|---|---|
| `R/fpr_calibration.R` | 18 (+14 / -4) |
| `man/fpr_calibration.Rd` | 4 (+3 / -1) |
| `tests/testthat/test-fpr-calibration-fallback.R` (new) | 128 (+128) |
| total | 3 files, 145 insertions, 5 deletions |

- Scope: only this fix's files. The committed diff equals the gated patch `patches/F7.patch` apart from index and hunk-header lines (`main/F7/patch_equal.txt` is empty; re-checked with `git show 65831b54`).
- The three files are unchanged from `65831b54` to the head of the fixes `61c15088`. `R/` is identical at `61c15088` and `88302212`.

#### Site, before and after

Before, at `76a62a71`:

`R/fpr_calibration.R:91-92` (end of `@param outcome_type_approx`, 89-92):
```r
#'   for computing the asymptotic \eqn{P_1} approximation.  If \code{NULL}
#'   (default), inferred from \code{fs_params$outcome_type}.
```
`R/fpr_calibration.R:236-237`:
```r
  # ── Infer outcome type for approximation ───────────────────────
  otype <- outcome_type_approx %||% fs_params$outcome_type %||% "binary"
```
`R/fpr_calibration.R:245-247` (end of the `switch()` at 241-247):
```r
      "continuous" = d_eff_continuous(n_min, sigma_y),
      d_eff_binary(n_min, p_event)  # fallback
    )
```
`man/fpr_calibration.Rd:80`: `(default), inferred from \code{fs_params$outcome_type}.}`

After, at `61c15088` (same text as `65831b54`):

`R/fpr_calibration.R:91-94`:
```r
#'   for computing the asymptotic \eqn{P_1} approximation.  If \code{NULL}
#'   (default), inferred from \code{fs_params$outcome_type}.  If neither is
#'   given, \code{"survival"} is used, which is the default
#'   \code{outcome_type} of \code{\link{forestsearch}}.
```
`R/fpr_calibration.R:238-248`:
```r
  # ── Infer outcome type for approximation ───────────────────────
  # Neither given: "survival", forestsearch()'s default outcome_type.
  otype <- outcome_type_approx %||% fs_params$outcome_type %||% "survival"
  otype_ok <- c("survival", "binary", "continuous", "count")
  if (!is.character(otype) || length(otype) != 1L || is.na(otype) ||
      !(otype %in% otype_ok))
    stop(sprintf(paste0(
      "Unrecognised outcome type %s (from outcome_type_approx or ",
      "fs_params$outcome_type); it must be one of %s."),
      paste(deparse(otype), collapse = ""),
      paste0("\"", otype_ok, "\"", collapse = ", ")))
```
`R/fpr_calibration.R:256-257` (end of the `switch()` at 252-257):
```r
      "continuous" = d_eff_continuous(n_min, sigma_y)
    )
```
`man/fpr_calibration.Rd:80-82`: the same added sentence, regenerated by `devtools::document()`. It wrote only `fpr_calibration.Rd` (`evidence/F7/document.txt`).

- The fallback value is forestsearch()'s default. At `76a62a71`: `R/forestsearch_main.R:1746-1747` `outcome_type = c("survival", "binary", "continuous", "count")`, and `:2215` `outcome_type <- match.arg(outcome_type)`.

#### Red and green

| run | where | result | wall |
|---|---|---|---|
| implementer red | scratch tree `wt/F7` (archive of `76a62a71`), test file present, fix stashed | `123456789ab.....`; `[ FAIL 11 \| WARN 0 \| SKIP 0 \| PASS 5 ]` (`red_counts.txt`). Key failures: `` Expected `res$call_args$otype` to be identical to "survival". `` (actual `"binary"`); `` Expected `res$d_eff` to equal `d_eff_survival(n_sg = 60, prop_cens = 0.25)`. `` (actual 15.0, expected 45.0); `` Expected `.fc_call(fs_params = fs_bad)` to throw a error. `` | 3.87 s (`red.txt`); 3.70 s (`red_counts.txt`) |
| implementer green | `wt/F7` with the fix | `................`; `[ FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 16 ]` | 3.42 s (`green.txt`); 3.48 s (`green_counts.txt`) |
| verifier red | `verify/F7/pre` (`76a62a71` + test file; `R/fpr_calibration.R` equal to the start copy) | `123456789ab.....`; 11 failed, 0 errors, 5 passed (blocks 1-3 fail 3, 3, 5; guard block passes 5). Same messages as the implementer's | 3.97 s |
| verifier green | `verify/F7/post` (`76a62a71` + patch) | `................`; 16 passed, 0 failed (3, 3, 5, 5) | 3.40 s |
| main red | repository at `f5b65516` (parent; `R/fpr_calibration.R` and its Rd unchanged since `76a62a71`) + test hunks only | `123456789ab.....`: 11 failures, 5 passes. The summary reporter lists 10, then `... and 1 more` and `Maximum number of 10 failures reached`. First failure: `` Expected `res$call_args$otype` to be identical to "survival". `` | 3.8 s |
| main green | repository, rest of the patch applied | `................`: 16 passes, 0 failures | 3.2 s |

- The 11 red failures are the defect: the binary fallback (d_eff 15 = 60 x 0.5 x 0.5) and the silent acceptance of `"survial"`. The 5 red passes are the guard block, which checks that a given type is still used.
- Counts: the main runs give only the progress line (one character per expectation), so the counts above are read from it.

#### Expected values

Test file at `65831b54`, 16 expectations.

| expected value | test line(s) | how derived |
|---|---|---|
| otype `"survival"` with no type given | 53, 69 | forestsearch()'s `outcome_type` formal lists `"survival"` first, and `match.arg()` returns the first element when the argument is not supplied (`R/forestsearch_main.R:1746-1747`, `:2215`). The spec says the fallback is forestsearch()'s default. |
| `d_eff_survival(n_sg = 60, prop_cens = 0.25)` | 58 | Evaluated in the test, as the spec requires. |
| `45` | 59 | By hand from `?d_eff_survival` (`d_eff = n_H (1 - p_c)`, `R/truefind_asymptotic_glm.R:34`): 60 x (1 - 0.25) = 45. |
| `d_eff_survival(n_sg = 60, prop_cens = 0)` and `60` | 73, 74 | Default `prop_cens = 0` from fpr_calibration()'s formals; 60 x (1 - 0) = 60. |
| error message contains `"survial"` | 85, 94, 98 | The message must name the value given. `"survial"` is not a substring of `"survival"`, so the recognised-value list cannot satisfy it. |
| error message contains `"survival"` | 88, 95 | The message must list the recognised values. Documented values (roxygen of `outcome_type_approx`): "binary", "survival", "count", "continuous". |
| otype `"binary"`, d_eff `15` (guard) | 111, 112 | `?d_eff_binary` `n_H p (1 - p)` at the default `p_event = 0.5`: 60 x 0.5 x 0.5 = 15. |
| otype `"continuous"`, d_eff `80 / 3` (guard) | 122, 123 | Roxygen: `fs_params$outcome_type` is used only when `outcome_type_approx` is NULL. `?d_eff_continuous` `n_H / sigma_y^2` = 60 / 1.5^2 = 80/3. |
| d_eff `30` (guard) | 127 | Roxygen: a supplied `d_eff_approx` is used as given. |

- Each value has a comment in the test saying how it was derived. The verifier re-derived every value and flagged none as lacking a comment.
- No expected value was taken from the fixed code's output.

#### Unaffected case

- Spec: "with a valid type given, output is identical; with no type given, the simulated false-positive rate is identical."
- Implementer (`evidence/F7/unaffected.R`):
  - Setup: start tree and fixed tree, each in its own `Rscript --vanilla` process with `pkgload::load_all()` and a sequential plan.
  - Replicates are real forestsearch() runs, not mocked: synthetic null data, N = 240, `nsims = 20`, `nsims_verify = 10`, `n_min = 40`, `c1 = 1.25`, `c2 = 1.0`, `set.seed(777)` before each call.

| case | compared | identical() | fpr |
|---|---|---|---|
| `outcome_type_approx = "survival"`, `prop_cens = 0.3` | object, print, warnings, `.Random.seed` after | TRUE | 0.55 |
| `fs_params$outcome_type = "survival"`, `prop_cens = 0.3` | same | TRUE | 0.55 |
| `fs_params$outcome_type = "survival"`, `d_eff_approx = 30` | same | TRUE | 0.55 |
| binary data, `fs_params$outcome_type = "binary"`, `p_event = 0.3` | same | TRUE | 0.5 |
| no type given, `prop_cens = 0.3` | FPR fields (`fpr, fpr_ci, n_found, nsims, n_failed, fpr_verify, fpr_verify_ci, nsims_verify, n_verify_failed`), warnings, `.Random.seed` after | TRUE | 0.55 |

- No type given:
  - FPR results: fpr 0.55, CI [0.342082, 0.741805], 11/20 found, fpr_verify 0.7.
  - Leaves that differ, all from the approximation:
    - `P1`: 0.23131 -> 0.200548.
    - `d_eff`: 10 -> 28.
    - `L_eff_implied`: 3.0353717610148196 -> 3.5674933075713406.
    - `call_args$otype`: `"binary"` -> `"survival"`.
- Wall clocks (implementer): start tree 22.35 s (22.16 s inside R); fixed tree 22.07 s (21.86 s inside R); compare 0.18 s. Result line: `UNAFFECTED CASE RESULT: IDENTICAL`.
- Verifier's re-run of the same script (start 21.86 s, post 21.65 s, compare 0.19 s): same result.
- Verifier's own script (`verify/F7/unaffected/my_unaffected.R`; start copy 11.31 s, post 10.64 s, compare 0.18 s):
  - Valid continuous, valid count, and DGM-mode survival: object, print, warnings and RNG all identical() TRUE.
  - DGM mode with no type, and no type with `d_eff_approx = 25`: FPR fields and RNG identical() TRUE.

#### Existing tests

- `git grep` at `76a62a71` finds one existing test file that touches `fpr_calibration`, `outcome_type_approx` or `fpr_calibration_transfer`: `tests/testthat/test-plan-restored.R`.
  - Its block at line 153 calls fpr_calibration() with no type, so it now takes the survival path.
- Implementer, fixed tree: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 20 ]`, 20.15 s (`evidence/F7/existing.txt`).
- Verifier, `verify/F7/post`: 8 tests, 20 passed, 0 failed, 20.09 s (`verify/F7/existing.txt`).
- Existing expectations updated: none.

#### What a user sees

`print()` of `fpr_calibration()` for a survival analysis with no type given (N = 240, `n_min = 40`, `prop_cens = 0.3`; `evidence/F7/unaffected.txt`, case `none`):

```
before:   P1 (per-subgroup): 23.1%  (d_eff = 10.0, uncorrected)
          Implied L_eff:     3.035
after:    P1 (per-subgroup): 20.1%  (d_eff = 28.0, uncorrected)
          Implied L_eff:     3.567
```

The `Sim FPR (cal)` (55.0%, CI 34.2% -- 74.2%), `Corrected FPR` (55.0%) and `Sim FPR (verify)` (70.0%) lines are unchanged.

#### Verification

- Verdict: CONFIRMED. Red, green and the unaffected case were reproduced.
- There was no repair. The main session made no edits: there is no `patches/F7.verified.patch`, and the commit equals `patches/F7.patch`.

| problem | severity | resolution |
|---|---|---|
| Exact matching (`%in%`) rejects abbreviations that forestsearch() accepts through `match.arg()`. `fs_params$outcome_type = "surv"` ran at the start SHA (otype `"surv"`, d_eff 10, fpr 0.4). It now stops with `Unrecognised outcome type "surv" ...` (`verify/F7/unaffected/my_compare.txt`, case `abbrev_fs`). | minor | Code left as specified; listed as an open item. |
| The notes' misspelt-type example shows `Error: Unrecognised outcome type "survial" ...`. That prefix comes from the example script's own `cat("Error:", conditionMessage(e))`. At top level R prints `Error in fpr_calibration(c1 = 1.25, c2 = 1, df.analysis = d, nsims = 2L,  : ` and then the message on the next line (`verify/F7/misspelt/toplevel_post.txt`). The message text is quoted exactly. | minor | Noted. This report's example is the `print()` output, not the error line. No code change. |
| `fpr_calibration_transfer()` drops the nested call's error, so a misspelt `outcome_type_approx` now gives NA silently (`verify/F7/unaffected/transfer_misspelt.txt`). | minor | Left as an open item; the implementer had already disclosed it. |

#### Departures and interpretations

- **Validation is unconditional.** The spec says "An unrecognised type stops with a clear message" with no condition. So the resolved type is checked before the replicates run, even when `d_eff_approx` is supplied. A `stop()` in the switch's default arm would fire only when `d_eff_approx` is NULL; that variant was not used.
- **Matching is exact** (`%in%`), not `match.arg()`. Abbreviations such as `"surv"` stop. A non-character, NA or length > 1 value stops with the same message. The documented values are the four exact strings.
- **The switch's binary fallback arm was deleted**, not replaced. After the up-front check it cannot be reached.
- **One roxygen sentence** states the fallback, as specified. The stop is not described further; the existing "One of ..." list covers it.
- **The test stubs the `%dofuture%` replicate loop** with `local_mocked_bindings()`, the device `test-plan-restored.R` already uses. The type and `d_eff` are resolved before the loop. The unaffected case runs real replicates.
- **The misspelt-type checks were split before green.** The first red run also showed a cascade error (`no applicable method for 'conditionMessage' applied to an object of class "NULL"`). The checks were then split into two `expect_error()` calls each, with the same regexps; no expected value changed.

#### paper_simulations

- No evaluated call reaches the changed lines.
- The changed lines are reachable only through `fpr_calibration()` or `fpr_calibration_transfer()`. `vignettes/articles/paper_simulations.qmd` (2860 lines, 47 R chunks at `76a62a71`) contains no occurrence of "fpr" (case-insensitive, also at `61c15088`).

#### Open items

- `R/fpr_calibration.R:242-243` at `61c15088`: exact matching rejects type abbreviations that forestsearch() accepts through `match.arg()` (`R/forestsearch_main.R:2215`), for example `fs_params$outcome_type = "surv"`. Not fixed here because the spec documents the four exact strings. A later fix could resolve the type with `match.arg(otype, otype_ok)`.
- `R/fpr_calibration.R:569` at `61c15088` (`:559` at `76a62a71`): `fpr_calibration_transfer()` wraps the nested call in `tryCatch(..., error = function(e) e)` and drops the error.
  - A misspelt `outcome_type_approx` passed through `...` gave `fpr_sim` 0, CI [0, 0.657628], `n_sim` 2. It now gives `fpr_sim`, `fpr_sim_lo`, `fpr_sim_hi` NA and `n_sim` 0, with no message (`verify/F7/unaffected/transfer_misspelt.txt`).
  - Not fixed here: outside the spec. A later fix could warn with `conditionMessage()`.
- Same function: for a `cal_result` made with no type and no `d_eff_approx`, `P1_base`, `P1_new` and `fpr_pred` change with `d_eff` (verifier: `P1_base` 0.23131 -> 0.194358, `my_compare.txt`, `transfer_none`). This follows from the fix. It is not tested because the transfer function has no tests.
- `R/fpr_calibration.R:99-101` at `61c15088` (`:97-99` at `76a62a71`): `@param prop_cens` says it is used only when `outcome_type_approx = "survival"`. It is used whenever the resolved type is survival, which now includes the fallback. `@param p_event` (`:102-104`) has the same wording for binary. Not fixed here: the spec allows one roxygen sentence.
- `R/fpr_calibration.R:255` at `61c15088` (`:244` at `76a62a71`): `"count"` uses `d_eff_binary(n_min, p_event)  # approximate` although `d_eff_count()` exists (`R/truefind_asymptotic_glm.R:114`). This predates the fix and is outside A-107.
- `R/fpr_calibration.R` has 16 lines with non-ASCII characters at both SHAs (box-drawing rules, em dashes, `~±`). They flow into `man/fpr_calibration.Rd` (for example lines 52, 57 and 164 at `61c15088`). They predate the fix; the commit adds no non-ASCII line.

#### Corrections to the notes

- The verifier's notes give `R/fpr_calibration.R` as "+13/-3 net (18 lines changed)". `git show --numstat 65831b54` gives +14 / -4 (18 lines).
- The notes' misspelt-type "after" line is `conditionMessage()` output printed by the example script, not R's top-level error form (see Verification).

### F8 — A-047: duplicate-cut test compares names exactly

**Status:** PASS. Commit `ea5f82d9` (parent `65831b54`, start SHA `76a62a71`).

| file (`git show --stat ea5f82d9`) | lines changed |
|---|---|
| `R/get_fsdata.R` | 28 (+17 / -11) |
| `tests/testthat/test-duplicate-cut-names.R` (new) | 194 (+194) |
| total | 2 files, 211 insertions, 11 deletions |

- Scope: only this fix's files. The committed diff equals the gated patch `patches/F8.patch` (sha256 `defe795d...8d43`) apart from index and hunk-header lines: `main/F8/patch_equal.txt` is empty. There is no `patches/F8.verified.patch`.
- No roxygen line changed, so no `man/` page changed. Both files are ASCII only (`grep -P '[^\x00-\x7F]'` on the committed files: 0 hits).
- Both files are unchanged from `ea5f82d9` to `61c15088` and to `88302212` (`git log` on the two paths is empty; `R/get_fsdata.R` has the same md5 at all three SHAs).

#### Site, before and after

All three sites are in `R/get_fsdata.R`, inside `get_FSdata()`. The spec's pointers (near 451, 499, 520) match the start SHA.

| site | start `76a62a71` | head `61c15088` (= `ea5f82d9`) |
|---|---|---|
| median removal (`replace_med_grf`) | changed :448-453 | :448-452 |
| GRF merge | changed :498-499 | :497-502 |
| DINA merge | changed comment :515-516, code :520 | comment :518-521, code :525-526 |

Median removal. Before, `76a62a71` :446-453:
```r
  if (replace_med_grf) {
    if (use_grf && length(conf.cont_medians) > 0 && length(grf_cuts) > 0) {
      # Find which conf.cont_medians are present in any grf_cuts
      to_exclude <- vapply(
        conf.cont_medians,
        function(x) any(grepl(x, grf_cuts)),
        logical(1)
      )
```
After, `61c15088` :446-452:
```r
  if (replace_med_grf) {
    if (use_grf && length(conf.cont_medians) > 0 && length(grf_cuts) > 0) {
      # Find which conf.cont_medians are the variable of a GRF cut.  The
      # variable name of a cut is the text before "<=" or ">=", trimmed;
      # names are compared exactly (not as substrings or regexes).
      grf_cut_names <- trimws(sub("\\s*[<>]=.*$", "", grf_cuts))
      to_exclude <- conf.cont_medians %in% grf_cut_names
```

GRF merge. Before, `76a62a71` :497-500:
```r
    if(length(confs) > 0 & length(grf_cuts) > 0){
      # Vectorized check: keep only grf_cuts not matching any confs
      flag_omit <- sapply(grf_cuts, function(cut) any(sapply(confs, function(x) grepl(x, cut))))
      grf_cuts_keep <- grf_cuts[!flag_omit]
```
After, `61c15088` :496-503:
```r
    if(length(confs) > 0 & length(grf_cuts) > 0){
      # Drop a GRF cut only if its variable name (the text before "<=" or
      # ">=", trimmed) is exactly a categorical confounder, or the whole cut
      # string is exactly in the current pool.  Exact comparison, not a
      # substring or regex match.
      grf_cut_names <- trimws(sub("\\s*[<>]=.*$", "", grf_cuts))
      flag_omit <- (grf_cut_names %in% conf.categorical) | (grf_cuts %in% confs)
      grf_cuts_keep <- grf_cuts[!flag_omit]
```

DINA merge. Before, `76a62a71` :514-521:
```r
  # DINA cuts: merged into the candidate pool exactly like GRF cuts -- keep
  # only those not already represented in confs (substring match), then
  # union-dedup.  dina_cuts are already canonical "var <= value" expressions.
  if(length(dina_cuts) > 0){
    if(details) cat('Factors per DINA:', c(dina_cuts), '\n')
    if(length(confs) > 0){
      flag_omit <- sapply(dina_cuts, function(cut) any(sapply(confs, function(x) grepl(x, cut))))
      dina_cuts_keep <- dina_cuts[!flag_omit]
```
After, `61c15088` :517-527:
```r
  # DINA cuts: merged into the candidate pool exactly like GRF cuts -- keep
  # only those not already represented in confs (exact match: the cut's
  # variable name is a categorical confounder, or the whole cut string is
  # already in confs), then union-dedup.  dina_cuts are canonical
  # "var <= value" or "var >= value" expressions.
  if(length(dina_cuts) > 0){
    if(details) cat('Factors per DINA:', c(dina_cuts), '\n')
    if(length(confs) > 0){
      dina_cut_names <- trimws(sub("\\s*[<>]=.*$", "", dina_cuts))
      flag_omit <- (dina_cut_names %in% conf.categorical) | (dina_cuts %in% confs)
      dina_cuts_keep <- dina_cuts[!flag_omit]
```

- `conf.categorical` at both merges is the vector that built the pool: `confs <- c(conf.categorical, conf.cont_Medcuts)` (start :486, head :485).

#### Red and green

Test file `tests/testthat/test-duplicate-cut-names.R`: 7 blocks, 13 expectations. Red progress line on every tree: `duplicate-cut-names: 12...3456..78` (5 pass, 8 fail).

| run | where | result | wall |
|---|---|---|---|
| implementer red | scratch `wt/F8` (archive of `76a62a71`), test file only | 8 failures, 0 errors, in 5 of 7 blocks: test lines 105, 116, 142, 143, 157, 158, 192, 193. E.g. `Failure ('test-duplicate-cut-names.R:105:3')` … `` `actual`: 0 `` / `` `expected`: 1 ``; `:142:3` … `Absent: "x <= 25.3"`; `:192:3` … `Absent: "bm <= 2.5"` | 3.25 s |
| implementer green | same scratch, full patch | 13 pass, 0 failures, 0 errors | 3.16 s |
| verifier red | `verify/F8/pre` (archive of `76a62a71` + test) | 8 failures, 0 errors, same lines and messages | 3.34 s |
| verifier green | `verify/F8/post` (archive + patch) | 7 blocks, 13 pass, 0 failed, 0 errors | 3.36 s |
| main red | repository at `65831b54` (preceding commit) + test file only | 8 failures, 0 errors, same lines and messages | 3.3 s |
| main green | repository after the code hunks were applied (committed as `ea5f82d9`) | 13 pass, 0 failures, 0 errors | 3.1 s |

- The 5 that pass on red are the guard expectations (lines 127, 128, 129, 168, 177). They pass on both trees by design.
- Every red failure is the substring drop of A-047: `grepl("bm", "bmi <= 25.3")`, or `grepl("x <= 2", "x <= 25.3")`. None is a typo, missing object or error.
- Sources: `evidence/F8/red.txt`, `green.txt`; `verify/F8/red.txt`, `green.txt`; `main/F8/red.txt`, `green.txt`, `wall.txt`.

#### Expected values

All three datasets are built without random numbers. Facts used (test header, lines 15-29): `?is.continuous` (at least `cont.cutoff` = 4 distinct values means continuous); `cut_type = "median"` cuts at `round(median, 2)`; `cut_type = "default"` cuts at mean, median, Q1, Q3 (type-7 quantiles). The helper sets `use_lasso = FALSE`, `collapse_cuts = FALSE` (lines 83-94).

| expected value (test line) | how it was derived |
|---|---|
| `sum(res$confs == "bmi <= 25.3")` = 1L (:105) | By hand, dataset A. The cut's name "bmi" is not exactly "bm". The pool before the merge is {"bm"}. Default bmi cuts are 26, 26, Q1 = 20 + 0.75*4 = 23, Q3 = 28 + 0.25*4 = 29, so none is 25.3. |
| `sum(res$confs == "bmi >= 25.3")` = 1L (:116) | By hand, dataset A. Same reasoning. All default cuts use `<=`. |
| `"bm <= 0"` count 0L, `"bm >= 1"` count 0L (:127, :128) | By hand. The name "bm" is exactly the categorical confounder (2 distinct values < 4, `?is.continuous`). |
| `"bm"` count 1L (:129) | By hand. The bare categorical name stays once (bm == 1 on 50 of 100 rows). |
| setequal {"x <= 2", "x <= 25.3"}, length 2 (:142, :143) | By hand, dataset B. Median of x: the 50th and 51st sorted values are both 2, so the cut is "x <= 2". With `replace_med_grf = FALSE` the pool is {"x <= 2"}. "x" is not categorical and the strings differ, so the GRF cut is kept. `cut_type = "median"` adds no default cuts. |
| setequal {"x <= 2", "x <= 25.3"}, length 2 (:157, :158) | By hand, dataset B. `replace_med_grf = TRUE` and the GRF cut "x <= 2" is on x, so the median cut is removed. The GRF cut is kept, giving {"x <= 2"}. The DINA cut "x <= 25.3" is a different string, so it is kept. |
| `"x <= 2"` (:168) | By hand, dataset B. The GRF cut equals the median string, so it is dropped. |
| `"x <= 25.3"` (:177) | By hand, dataset B. The GRF cut removes the median cut and is kept. The DINA cut is the same string, so it is dropped. |
| setequal {"bm <= 2.5", "bmi <= 25.3"}, length 2 (:192, :193) | By hand, dataset C. bm (1..4) has 4 distinct values, so it is continuous. Its median is (2 + 3)/2 = 2.5. Only bmi's median cut "bmi <= 26" is removed, because the GRF name is exactly "bmi". |

- The verifier re-derived all 13 by hand before looking at the fixed tree's output. It also checked the summary statistics numerically on both trees (`verify/F8/derive_check.txt`: `median(B$x)= 2  median(C$bm)= 2.5  bmi mean/median/q1/q3= 26 26 23 29`).
- The verifier found every expected value has a derivation comment. None is flagged as missing one.

#### Unaffected case

- **Spec:** confounder names with no overlap give identical candidate pools.
- **What was run** (`evidence/F8/unaffected.R`):
  - Synthetic data, n = 300, `set.seed(20261112)`.
  - Confounders `age`, `wt`, `sz` (continuous) and `grp`, `stg` (categorical). `stopifnot()` checks that no name, and no pool string, is a substring of another cut.
  - 7 `get_FSdata()` configurations (c1-c7): default + GRF; GRF + DINA; median with `replace_med_grf` TRUE and FALSE; LASSO + GRF; `conf.cont_medians_force`; `collapse_cuts = FALSE`.
  - 2 `forestsearch()` runs (f1 `use_grf = TRUE`, f2 `use_grf = TRUE, use_dina = TRUE`), `seedit = 42`, `set.seed(1)` before each call.

| comparison | result | walls |
|---|---|---|
| implementer: start vs fixed (`wt/F8`) | `OVERALL: IDENTICAL`, 9 of 9 objects (timing leaves dropped) | 5.98 s / 5.90 s; compare 1.60 s |
| implementer control: start vs start re-run | `OVERALL: IDENTICAL` | 5.74 s; compare 1.56 s |
| verifier: start vs post, pre vs post | `OVERALL: IDENTICAL` for both | 5.97 s / 5.90 s / 6.00 s |
| verifier: plain `identical()` on whole objects | TRUE for c1-c7; FALSE for f1 and f2, which differ only in `find.grps$time_search` and `minutes_all` | — |
| main session, cumulative: start vs `head_fixes` (all 14 fixes; `R/` equal to `61c15088`) | `OVERALL: IDENTICAL`; controls (gated start vs re-run start; gated `wt/F8` vs `head_fixes`) also IDENTICAL | 6.11 s / 5.98 s; compare 1.59 s; controls 1.62 s each |

Sources: `evidence/F8/unaffected.txt`; `verify/F8/unaffected_rerun.txt`, `unaffected_leafdiff.txt`; `cumulative/F8/unaffected.txt`, `control.txt`.

#### Existing tests

- **Implementer** (`evidence/F8/existing.txt`): 27 files on the fixed scratch tree. These are files that call `get_FSdata()`, use `.fs_args_for()`, or set `use_grf`/`use_dina = TRUE`. Each file ran in its own process (120 s cap, 6 at a time).
  - 26 files completed with 0 failures and 0 errors. Per-file walls run from 3.6 s to 46.0 s.
  - `mr-inference` hit the cap twice (exit 124, 120.1 s), once alongside the others and once alone.
- **Verifier** (`verify/F8/existing_summary.txt`): 29 files on `verify/F8/post`, the 27 plus `build_classification_table` and `fs_mr_oc_summary`. Total wall 156.5 s.
  - 28 files completed with 0 failed and 0 errors.
  - `mr-inference` exit 124 at 120.06 s.
  - The 5 files with warnings show the same counts on `pre` and `post` (`warnings_pre_vs_post.txt`).
- `test-mr-inference.R` cannot reach the changed lines. All three of its configurations pass `use_grf = FALSE` (`76a62a71` :28, :304, :450), and none sets `use_dina` or `dina_cuts`.
- **Existing expectations updated:** none.

#### What a user sees

From `evidence/F8/before_after.R` and `before_after.txt`:
- The call is `get_FSdata(df, confounders.name = c("bm", "bmi", "age"), outcome.name = "time", event.name = "event", use_grf = TRUE, grf_cuts = "bmi <= 25.3", dina_cuts = "bmi >= 25.3", details = FALSE)`.
- `bm` is binary, n = 100. Package defaults otherwise, collapse on.

```
== start tree
res$confs:
[1] "bmi <= 26" "bmi <= 23" "bmi <= 29" "age <= 60" "age <= 50" "age <= 70"
[7] "bm"       
== fixed tree
res$confs:
[1] "bmi <= 25.3" "bmi >= 25.3" "bmi <= 26"   "bmi <= 23"   "bmi <= 29"  
[6] "age <= 60"   "age <= 50"   "age <= 70"   "bm"         
```

- Before the fix, both screening cuts on `bmi` were dropped only because "bm" is a prefix of "bmi". After the fix, both are candidate factors.
- Walls were 2.63 s and 2.56 s. The verifier reproduced the output exactly (`verify/F8/before_after_rerun.txt`, 2.30 s and 2.53 s).

#### Verification

- **Verdict: CONFIRMED.** The verifier reproduced red, green and the unaffected case.
- **Mutation checks** (`verify/F8/mutants.txt`):
  - M1 drops the categorical-name clause. It fails :127 and :128.
  - M3 restores the old median `grepl`. It fails :192 and :193.
  - M2 drops the whole-string clause. It passes all 13 (see problem 1).
- **Extra probes** (`verify/F8/probe2.txt`): with categorical `z1` and the GRF cut `z10 <= 3.5`, and with categorical `x.1` (a regex) and the GRF cut `xa1 <= 3.5`, the cut is dropped before the fix and kept after.

| problem | severity | resolution |
|---|---|---|
| The guard "an exact duplicate cut is still dropped" (test :162-178) cannot tell the `cuts %in% confs` clause apart from having no such clause. `confs <- unique(c(confs, *_cuts_keep))` (head :504, :531) already removes exact duplicates, so mutant M2 passes. | minor | Noted, not fixed. The committed diff is the gated patch unchanged. The verifier asked for this caveat in the report, and it is stated here: the spec's "an exact duplicate is still dropped" holds as an observable outcome. The clause's only observable effect is that the duplicate stays out of `grf_cuts_keep`/`dina_cuts_keep`, and so out of `collapse_redundant_cuts()`'s `protect` list (head :576). The optional extra assertion was not added. See open items. |
| In the unaffected evidence, the `forestsearch()` rows print `confs=0`. That is `length(fit$confs)`, and a forestsearch object has no `confs` element. The pool is `fit$confounders.evaluated`, with 11 entries (`verify/F8/unaffected_leafdiff.txt`). | minor | Noted, not fixed. It is an evidence label only, and the whole-object comparison stands. The same script ran unchanged in the cumulative check (`cumulative/F8/unaffected.R` is byte-identical), so `confs=0` still prints there. |

#### Departures and interpretations

- **Operators.** The name is parsed from `<=`/`>=` only, as the spec words it. The audit's sketch used `[<>]=?`. GRF writes only `<=`, and DINA writes only `<=`/`>=`. A hand-made cut with a bare `<` or `>` keeps its whole string as its "name", so it is dropped only as an exact duplicate.
- **"The current pool"** was read as `confs` at each merge. At the DINA merge, `confs` also holds the kept GRF cuts, so exact GRF/DINA duplicates are still dropped.
- **"A categorical confounder"** was read as `conf.categorical` at the merge. This is the vector that built the pool, after any LASSO filtering.
- **Comment.** The DINA comment now says `"var <= value" or "var >= value"`. It sits in the edited lines, and the new parse depends on it. The `replace_med_grf` roxygen ("removes median cuts that overlap with GRF cuts") was left as is, because it still describes the behaviour.
- **Tests beyond the spec's list.**
  - "`x <= 25.3` kept when the pool holds `x <= 2`" is tested twice: at the GRF merge, with a median cut in the pool (:133-144), and at the DINA merge, with a kept GRF cut in the pool (:147-159).
  - A median-removal block (:181-194) covers the third site, which the spec's Change bullet names but its Tests bullet does not.
  - The categorical guard covers both GRF and DINA (:120-130).

#### paper_simulations

- The implementer reports that the evaluated chunks `run-alt-sims` (:908) and `run-null-sims` (:966) reach the changed GRF-merge line through the FS arm (`fs_params$use_grf = TRUE`). They do not reach the median-removal or DINA lines.
- The output is unchanged, both by argument (no name in `confounders_base` is a substring of another) and in a one-trial traced check. That check gave the same result row apart from the elapsed-time column `fs_minutes` (`evidence/F8/paper_sim_check.txt`).

#### Open items

- **paper_simulations reaches the changed GRF-merge line.** The route is `vignettes/articles/paper_simulations.qmd` chunks `run-alt-sims` (:908) and `run-null-sims` (:966) at `76a62a71`. In the FS arm, `fs_params` (:850-879, `use_grf = TRUE` at :856) overrides the arm's `use_grf = FALSE` (`R/run_simulation_analysis.R:587-591`). The result is unchanged by argument and by the one-trial check. The article was not re-rendered, because `vignettes/` is read-only in this task and no campaign is run.
- **`conf.cont_medians_force` is appended to `conf.cont_medians` twice.** This is in `R/get_fsdata.R` at start :434-436 and :470 (head :434-436 and :469). With GRF/DINA off and no `conf_force`, the median cut appears twice in the pool. `evidence/F8/side_issue_medians_force.txt` shows `"wt <= 75.85" "wt <= 75.85" "grp"`, the same on both trees. Not fixed: it is unrelated to A-047 and not in F8's authorized change.
- **`get_cut_name()` uses the same unanchored `grepl(x, thiscut)`.** This is at `R/get_FSdata_helpers.R` start :453-456 (head :463-466). Its only caller, `is_flag_drop()` (start :497, head :507), has no caller in `R/`, `tests/` or `NAMESPACE` at either SHA. Elsewhere it appears only in `man/is_flag_drop.Rd`, `_pkgdown.yml` and a profiling regex string in `quarto/`, never as a call. Not fixed: outside the spec, and no result depends on it.
- **`test-mr-inference.R` did not finish under the 120 s per-run cap** in the implementer's two runs or the verifier's run. It cannot reach the changed lines (`use_grf = FALSE` at :28, :304, :450; no DINA). It passed in the full-suite run at `88302212` (main report, suite section: 0 failures).
- **The exact-duplicate guard does not pin the `%in% confs` clause.** This is `tests/testthat/test-duplicate-cut-names.R:162-178` (verifier problem 1, mutant M2). An assertion through the collapse `protect` path (`collapse_cuts = TRUE`) would pin it, with its expected value worked out by hand. Not added: the gated patch was committed unchanged.

#### Corrections to the notes

- `notes.md` puts the first `conf.cont_medians_force` append at start `R/get_fsdata.R:433-435`. Git shows :434-436; line 433 is blank.
- `notes.md` cites `R/run_simulation_analysis.R:587-590` for the FS-arm `modifyList()`. Git shows the statement spans :587-591; the closing parenthesis is at :591.
- The "164 s total" for the implementer's existing-test run appears only in `notes.md`. No evidence file records it. The per-file walls are in `evidence/F8/existing.txt`.

### F9 — A-125: MRCT default confounders exclude internal spline columns

**Status:** PASS. Commit `a30d6e7c` (parent `ea5f82d9`). Files (`git show --stat a30d6e7c`):

| File | Lines (`--stat`) | `--numstat` |
|---|---|---|
| `R/mrct_case_study.R` | 20 | +18 / -2 |
| `R/mrct_simulation.R` | 22 | +20 / -2 |
| `man/mrct_case_study.Rd` | 7 | +6 / -1 |
| `man/mrct_region_sims.Rd` | 10 | +9 / -1 |
| `tests/testthat/test-mrct-default-confounders.R` (new) | 83 | +83 / -0 |

- 5 files changed, 136 insertions(+), 6 deletions(-).
- The five files are unchanged from `a30d6e7c` to head `61c15088` (`git diff --stat` is empty). `R/` is unchanged from `61c15088` to `88302212`. So head line numbers below hold at `88302212` too.
- The commit diff equals the gated patch `patches/F9.patch`, apart from index and hunk-header lines (`main/F9/patch_equal.txt` is empty; re-checked with `diff`: same).

#### Site, before and after

`mrct_region_sims()`, NULL branch. Start `R/mrct_simulation.R:321-324` (block 315-325):

```r
    confounders.name <- setdiff(
      all_vars[grepl("^z_|^ecog|^strat", all_vars)],
      exclude_vars
    )
```

Head `R/mrct_simulation.R:329-342` (block 323-343; new lines 329-338, changed line 341):

```r
    # A mechanism built with a spline carries three internal columns of the
    # spline variable: <var>_treat, <var>_k and <var>_k_treat. They are not
    # baseline covariates. Exclude exactly these names (never a suffix match,
    # which would also drop a covariate such as z_prior_treat).
    spline_info <- dgm$model_params$spline_info
    spline_helpers <- if (!is.null(spline_info)) {
      paste0(spline_info$var, c("_treat", "_k", "_k_treat"))
    } else {
      character(0)
    }
    confounders.name <- setdiff(
      all_vars[grepl("^z_|^ecog|^strat", all_vars)],
      c(exclude_vars, spline_helpers)
    )
```

`mrct_case_study()`, defaults list. Start `R/mrct_case_study.R:291-293` (block 285-293):

```r
      c("id", "treat", "tte", "event", "y_sim", "treat_sim", "event_sim", "t_true",
        "c_time", "lin_pred_0", "lin_pred_1", "loghr_po", "entrytime", "flag_harm")),
    conf_force = c("z_age <= 65", "z_bm <= 0", "z_bm <= 1", "z_bm <= 2", "z_bm <= 5"))
```

Head `R/mrct_case_study.R:291-300` (new block, same text as in `mrct_region_sims()`, before `defaults <- list(` at 301) and 306-309 (changed 307, new 308):

```r
  spline_info <- dgm$model_params$spline_info
  spline_helpers <- if (!is.null(spline_info)) {
    paste0(spline_info$var, c("_treat", "_k", "_k_treat"))
  } else {
    character(0)
  }
  ...
      c("id", "treat", "tte", "event", "y_sim", "treat_sim", "event_sim", "t_true",
        "c_time", "lin_pred_0", "lin_pred_1", "loghr_po", "entrytime", "flag_harm",
        spline_helpers)),
    conf_force = c("z_age <= 65", "z_bm <= 0", "z_bm <= 1", "z_bm <= 2", "z_bm <= 5"))
```

(Lines 291-294 are the same four comment lines as in `mrct_region_sims()`.)

Roxygen:

- `R/mrct_simulation.R:52` at start: `#'   If NULL, automatically extracted from dgm`. At head, 52-60: `If NULL, automatically extracted from dgm: the columns of \code{dgm$df_super} whose names start with \code{z_}, \code{ecog} or \code{strat}, except that, for a mechanism built with a spline (\code{dgm$model_params$spline_info}), the three internal columns ... are left out by their exact names \code{<var>_treat}, \code{<var>_k} and \code{<var>_k_treat} (for \code{\link{create_dgm_for_mrct}}: \code{z_bm_treat}, \code{z_bm_k}, \code{z_bm_k_treat}). A covariate such as \code{z_prior_treat} and the region indicator stay in the set.`
- `R/mrct_case_study.R:172` (`@param fs_args`, block 166-175 at start, 166-180 at head): after "reproduce its identifier." the head adds 172-177: `The default \code{confounders.name} is that of \code{mrct_region_sims()}: for a mechanism built with a spline, it leaves out the three internal columns of the spline variable \code{<var>} by their exact names ...; a covariate such as \code{z_prior_treat} and the region indicator stay.`
- Rd pages: `man/mrct_region_sims.Rd:63` (start) -> 63-71 (head); `man/mrct_case_study.Rd:43` (start) -> 43-48 (head). Regenerated by `devtools::document()`; `tools::checkRd()` 0 problems on both (`evidence/F9/checkRd.txt`, `verify/F9/checkRd.txt`).

#### Red and green

| Run | Where | Result | Wall |
|---|---|---|---|
| Implementer red | scratch `wt/F9` (archive of `76a62a71`) + test file only | `...1.2..`: 2 failures, 6 passes. `test-mrct-default-confounders.R:65:3`: "Expected `rs_set` to be identical to `expected`", extra `[7] "z_bm_treat"`, `[8] "z_bm_k"`, `[9] "z_bm_k_treat"`. `:76:3`: same for `cs_set` | 7.58 s |
| Implementer green | scratch `wt/F9` with the fix | `........`: 0 failures, 8 passes. Progress reporter: `[ FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 8 ]` (7.32 s) | 7.18 s |
| Verifier red | `verify/F9/pre` + test file | `[ FAIL 2 \| WARN 0 \| SKIP 0 \| PASS 6 ]`; same two failures at `:65:3` and `:76:3` | 7.55 s |
| Verifier green | `verify/F9/post` | `[ FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 8 ]` | 7.22 s |
| Main-session red | repository at `ea5f82d9` + the patch's test hunk only | `...1.2..`: 2 failures, 6 passes; same messages at `:65:3` (`rs_set`) and `:76:3` (`cs_set`) | 6.9 s |
| Main-session green | repository after the rest of the patch | `........`: 0 failures, 8 passes | 6.7 s |

- The two F9 R files are the same at `76a62a71` and `ea5f82d9` (`git diff --stat` between them lists only other fixes' test files).
- Red fails only on the three helper columns, in both functions. The preconditions, `expect_length`, both `z_prior_treat` checks and the identifier check pass on the unfixed tree.

#### Expected values

All from `tests/testthat/test-mrct-default-confounders.R` at `a30d6e7c` (8 expectations).

| Expected value (line) | How derived |
|---|---|
| `"z_bm"` = `spline_info$var` (:37) | `create_dgm_for_mrct()` source passes `spline_spec` with `var = "z_bm"`. Comment :33-36. |
| `z_bm_treat`, `z_bm_k`, `z_bm_k_treat` in `names(df_super)` (:38) | `create_spline_variables()` builds `paste0(spline_var, "_treat" / "_k" / "_k_treat")`. Comment :34-36. |
| `1L` printed `Confounders:` line (:63) | One replicate (`n_sims = 1`, `verbose_n_sims = 1`). **No comment gives this count**; the verifier flagged it (minor). |
| `expected <- c("z_age", "z_bm", "z_male", "z_histology", "z_prior_treat", "z_regA")` (:50; used at :65 and :76) | By hand: `create_dgm_for_mrct()` passes `continuous_vars = c("age", "bm")` and `factor_vars = c("male", "histology", "prior_treat", "regA")`; `process_continuous_vars()` then `process_factor_vars()` write the `z_` columns in that order; a two-level factor is one column; region indicator stays (spec). Comment :40-49. The verifier re-derived it from source and found it correct, but flagged (minor) that the comment does not say why other `df_super` columns are absent (raw copies `age`, `bm`, ...; censoring columns `zcens_age`, `zcens_prior_treat`; none match `^z_`). |
| `TRUE` for `"z_prior_treat" %in%` set (:66, :77) | F9 spec: "it contains `z_prior_treat`". Comment :47-49 and file header :8-9. |
| `cs$Hhat` identical to `as.character(rs$sg_found)` (:82) | `mrct_case_study()` roxygen: a run's direct arguments and its `fs_args` "together reproduce its identifier". Comment :79-81. The verifier flagged (minor) that it passes on both trees, so it is not red/green evidence. |

- No expected value was copied from the fixed code's output (implementer notes; verifier re-derived all from the source at `76a62a71`).

#### Unaffected case

- Spec: "with `confounders.name` given, output is identical."
- Script: `evidence/F9/unaffected.R` (verifier copy `verify/F9/unaffected_verify.R`, which differs only in paths). Each tree in its own `Rscript --vanilla`, `timeout 120`, `pkgload::load_all()`, sequential plan, `set.seed(2026)` before each call. Comparison: `identical()` after dropping timing leaves by name.
- Trees: implementer `start` (export of `76a62a71`) vs scratch `wt/F9`; verifier `start` vs `verify/F9/post`. The verifier checked that the `start` export equals the git archive of `76a62a71`.
- The main session did not run this case on the repository (no such file in `main/F9/`).

| Case | Implementer | Verifier |
|---|---|---|
| `rs_given`: `mrct_region_sims(confounders.name = <six>)`, 3 replicates, `record = "full"` | TRUE | TRUE |
| `rs_given_helpers`: same with the six + 3 helpers | TRUE | TRUE |
| `rs_given_line`: printed `Confounders:` line, `details = TRUE` | TRUE | TRUE |
| `cs_given`: `mrct_case_study(fs_args$confounders.name = <six>)` | TRUE | TRUE |
| `cs_given_helpers`: same with the nine | TRUE | TRUE |
| Extra: mechanism without spline info (`generate_aft_dgm_flex()`, no `spline_spec`), `confounders.name` omitted: `nospline_has_spline_info`, `nospline_rs_default`, `nospline_cs_default` | TRUE (3 of 3) | TRUE (3 of 3) |

- `ALL IDENTICAL: TRUE` in both runs; also TRUE with `ignore.environment = TRUE`.
- Walls (start / fixed / driver): implementer 15.89 / 16.20 / 32.36 s (`evidence/F9/unaffected.txt`); verifier 15.75 / 15.96 / 31.99 s (`verify/F9/unaffected_verify.txt`).

#### Existing tests

- Selection: every test file that references `mrct_region_sims` or `mrct_case_study`, directly or through helpers. Each in its own process, `NOT_CRAN=true`, `timeout 120`.

| File | Implementer (fixed scratch tree) | Verifier (`post`) |
|---|---|---|
| test-mrct-default-confounders.R (new) | PASS 8, FAIL 0, 7.32 s | 8, 7.26 s |
| test-mrct-trial-analysis.R | PASS 131, FAIL 0, 17.24 s | 131, 17.36 s |
| test-mrct-region-sims-mr.R | PASS 71, FAIL 0, 38.60 s | 71, 37.89 s |
| test-mrct-case-study.R | PASS 64, FAIL 0, 91.19 s | 64, 90.80 s |
| test-mrct-core-summary.R | PASS 32, FAIL 0, 8.22 s | 32, 8.56 s |
| test-mrct-region-sims-record.R | PASS 172, FAIL 0, 27.88 s | 172, 28.38 s |
| test-plan-restored.R | PASS 20, FAIL 0, 20.20 s | 20, 19.97 s |

- Implementer: all `[ FAIL 0 | WARN 0 | SKIP 0 ]` (`evidence/F9/existing.txt`).
- Verifier: all expectations with "non-dot=0" (`verify/F9/existing.txt`). Each run exited 1 because of a bug in the verifier's own summary `cat()` after testthat finished (`Summary.data.frame` error in `existing_mrct-default-confounders.log`), not a test failure.
- Existing expectations updated: none. All existing callers pass `confounders.name` explicitly.

#### What a user sees

From `evidence/F9/example.txt` (`mrct_region_sims(..., details = TRUE)`, `create_dgm_for_mrct()` mechanism, `confounders.name` omitted):

- Before (`76a62a71`, 7.17 s):
  `   Confounders:  z_age, z_bm, z_male, z_histology, z_prior_treat, z_regA, z_bm_treat, z_bm_k, z_bm_k_treat`
- After (fixed tree, 6.94 s):
  `   Confounders:  z_age, z_bm, z_male, z_histology, z_prior_treat, z_regA`
- `mrct_case_study()$call_args$fs_args$confounders.name` goes from the same 9 names to the same 6. `sg_found` and `Hhat` are `!{z_male}` on both trees. The verifier reproduced this (`verify/F9/example_verify.txt`, 7.05 / 7.04 s).

#### Verification

- Verdict: CONFIRMED. Red, green and the unaffected case reproduced. No repair or re-verify round.
- Problem 1 (minor): the last expectation (`cs$Hhat` vs `rs$sg_found`, :82) passes on both trees and cannot catch a one-function fix in this trial (the subgroup `!{z_male}` uses no helper). Suggested: a comment calling it a regression guard, or drop it. **Not applied**: open item.
- Problem 2 (minor): the derivation comment (:40-49) does not say why `df_super`'s raw copies and `zcens_age` / `zcens_prior_treat` are absent; `expect_length(line, 1L)` (:63) has no "(one replicate, verbose_n_sims = 1)" note. **Not applied**: open item.
- No main-session comment was added for F9 (no `patches/F9.verified.patch`; committed diff equals `F9.patch`).
- Verifier extra checks:
  - Candidate pool (`verify/F9/pool_check.txt`): before, the case study's `fs$confounders.candidate` and `fs$confounders.evaluated` together hold `z_bm_treat <= 3/0/5`, `z_bm_k <= 2/1/0/4`, `z_bm_k_treat <= 1/0`; after, none ("any helper in candidates: FALSE").
  - S2 build from `inject_mrct_structure()` (`verify/F9/inject_check_{pre,post}.txt`): after the fix the default equals the `z_` columns minus `z_cd40_treat`, `z_cd40_k`, `z_cd40_k_treat`, and equals `analysis_covariates(dgm, "observed")` plus `z_region` as sets (both TRUE; FALSE before).
  - `devtools::document()` on a copy of the fixed tree (roxygen2 8.1.1, 12.85 s) left both Rd pages unchanged.

#### Departures and interpretations

- No departure in the code (implementer).
- "Mechanism has spline information" read as `!is.null(dgm$model_params$spline_info)`; the names are built as `paste0(spline_info$var, c("_treat", "_k", "_k_treat"))`. No suffix match.
- "Identical change in both functions": the same comment and the same `spline_helpers` expression in both. In `mrct_case_study()` the block sits before `defaults <- list(` because the defaults are a list literal.
- "State the exclusion in the roxygen": `mrct_case_study()` has no `@param confounders.name`, so the exclusion is stated in `@param fs_args`. The `mrct_region_sims()` text also names the existing `z_` / `ecog` / `strat` pattern, so that "left out" has a referent.
- "The region indicator stays": not excluded; it is in the default because it matches `^z_` (`z_regA`).
- Test: the `mrct_region_sims()` set is read from the printed `Confounders:` line under `details = TRUE`; the `mrct_case_study()` set from `cs$call_args$fs_args$confounders.name` (documented in its `@return`). The test is `skip_on_cran()`.
- Additions beyond the spec: the identifier expectation (:82) and the no-spline mechanism in the unaffected script. Neither changes the scope.

#### paper_simulations

- Not reached. `paper_simulations.qmd` never names `mrct_region_sims` or `mrct_case_study` (only `create_dgm_for_mrct()` in a prose table, line 457 at `76a62a71`), and no package function calls either one (`git grep` of `R/` finds only messages).

#### Open items

- `R/mrct_simulation.R:543` (start) / `:561` (head): `fs_call_args <- modifyList(fs_args, fs_base_args)` still silently replaces a `confounders.name` passed through `mrct_region_sims(fs_args = ...)` with the default (A-125 "less ordinary path"). Out of F9 scope; after F9 the default no longer carries the helpers.
- `R/mrct_inject_dgm.R:635` (start and head): `analysis_covariates()` still drops columns with `grep("_treat$|_k$", z, value = TRUE)`, which would also drop a covariate such as `z_prior_treat` (A-124). Not part of F9.
- `R/mrct_case_study.R:294` (start) / `:310` (head): `fs_spec <- utils::modifyList(defaults, fs_args)` means `fs_args = list(confounders.name = NULL)` removes the element instead of taking the default. Existing behaviour, unchanged and not exercised by F9.
- `tests/testthat/test-mrct-default-confounders.R:40-49, :63, :82` (head): the verifier's two minor comment requests (why other `df_super` columns are absent; the count in `expect_length(line, 1L)`; :82 is a regression guard that passes before and after) were not applied before commit.
- Test caveat: the multisession workers in `test-mrct-case-study.R` and `test-mrct-region-sims-record.R` load the installed forestsearch 0.4.1.9002, as those files document. The changed lines run in the main process (`load_all()` of the tree), so the caveat does not affect F9's evidence.
- NEWS: resolved in the last commit `88302212` (`NEWS.md:355-359`, "(Code audit A-125.)").

#### Corrections to the notes

- `evidence/F9/notes.md` lists the 8 expectations as "2 preconditions; for each function, set identical to the hand list and contains `z_prior_treat`; and `cs$Hhat == rs$sg_found`", which names 7. git shows an eighth, `expect_length(line, 1L)` at :63.
- The verifier JSON gives "R/mrct_case_study.R +20/-2, R/mrct_simulation.R +22/-2". `git show --numstat a30d6e7c` gives +18/-2 and +20/-2 (20 and 22 are the `--stat` totals).
- The verifier's `expected_values_rederived` says "Every expected value has a comment saying how it was derived", while its own problem 2 says `expect_length(line, 1L)` lacks one. git agrees with problem 2: no comment at :63 gives the count.

### F10 — A-064: count outcomes with GRF are oriented

**Status:** PASS. Commit `f2e8de47` (parent `a30d6e7c`).

| File (from `git show --stat f2e8de47`) | Change |
|---|---|
| `R/grf_subg_harm_glm.R` | 23 lines (+15 / -8) |
| `tests/testthat/test-grf-count-orientation.R` | new, 138 lines |
| total | 2 files changed, 153 insertions(+), 8 deletions(-) |

- No roxygen, `man/`, DESCRIPTION, NAMESPACE or NEWS change in this commit.
- `R/grf_subg_harm_glm.R` is identical at `76a62a71` and `a30d6e7c`. It is also identical at `f2e8de47`, `61c15088` and `88302212`. So start-SHA line numbers are the pre-fix lines, and the head line numbers are the committed lines.

#### Site, before and after

All sites are in `R/grf_subg_harm_glm.R`, `grf.subg.harm.glm()`. The flip block comes after the count transform (start :409-434: `log(Y_raw + 0.5)` or `Y_raw`) and before `causal_forest()`.

**1. Header comment.** Start `76a62a71` :437-441:

```r
  # Adverse outcome flip: align CATE direction with ForestSearch convention.
  # When adverse_outcome = TRUE, Y is flipped (binary: 1-Y, continuous: -Y)
  # so that positive CATE = treatment increases the BAD outcome = harm.
  # This ensures the most-negative-CATE subgroup is the true "harm" group,
  # matching forestsearch()'s internal GRF pre-screening behavior.
```

Head `61c15088` :437-443:

```r
  # Adverse outcome flip: align CATE direction with ForestSearch convention.
  # When adverse_outcome = TRUE, Y is flipped (binary: 1-Y, continuous: -Y,
  # count: minus the transformed count above, -log(Y + 0.5) or -Y) so that a
  # higher Y_grf is better and positive -CATE = treatment increases the BAD
  # outcome = harm.  This ensures the most-negative-CATE subgroup is the true
  # "harm" group, matching forestsearch()'s internal GRF pre-screening
  # behavior.
```

**2. Flip block.** This is the code change. Start :448-452:

```r
      Y_grf <- -Y_grf
      if (verbose) message("[grf.subg.harm.glm] adverse_outcome=TRUE: using -Y for causal_forest()")
    }
    # count: log transform already applied above; flip not needed
  }
```

Head :450-458. The new branch is at :452-457.

```r
      Y_grf <- -Y_grf
      if (verbose) message("[grf.subg.harm.glm] adverse_outcome=TRUE: using -Y for causal_forest()")
    } else if (outcome_type == "count") {
      # count: the log(Y + 0.5) or identity transform above keeps "higher
      # count = higher Y_grf", so negate it: a higher count is then worse,
      # as for binary and continuous.
      Y_grf <- -Y_grf
    }
  }
```

**3. Frontier comment.** Start :515-517:

```r
  # path's score = -CATE convention.  When adverse_outcome = TRUE the Y-flip was
  # applied before causal_forest(), so dr_scores are already harm-oriented and
  # need no special handling here.  grf_selection = "tree" skips this entirely.
```

Head :521-524:

```r
  # path's score = -CATE convention.  When adverse_outcome = TRUE the Y-flip
  # (binary, continuous and count alike) was applied before causal_forest(), so
  # dr_scores are already harm-oriented and need no special handling here.
  # grf_selection = "tree" skips this entirely.
```

#### Red and green

The test has 2 `test_that` blocks with 14 expectations. In every red run, the 4 failures are at `test-grf-count-orientation.R` :110, :114, :134 and :136.

| Run | Where | Result | Wall |
|---|---|---|---|
| Implementer red | Scratch tree `wt/F10` (archive of `76a62a71`), test file only | `....1.2....3.4`: 10 pass, 4 failures, 0 errors. :110 and :134: "Actual comparison: 0.00 < 0.80". :114 and :136: "Expected \`identical(res_adv$data$treat.recommend, res_ben$data$treat.recommend)\` to be FALSE" with "\`actual\`: TRUE" (`evidence/F10/red.txt`) | 4.60 s |
| Implementer green | Same scratch tree (start SHA) with the code applied | `..............`: 14 pass, 0 failures, 0 errors, 0 skips (`evidence/F10/green.txt`) | 4.29 s |
| Verifier red | Verifier's copy `verify/F10/pre` (start + test) | Same 4 failures at the same lines, same messages (`verify/F10/pre_test.txt`) | 4.60 s |
| Verifier green | Verifier's copy `verify/F10/post` (patch applied) | 14 pass, 0 failures (`verify/F10/post_test.txt`) | 4.54 s |
| Main red | Repository at `a30d6e7c` (preceding commit) with only the test file | `....1.2....3.4`: 4 failures at :110, :114, :134 and :136, same messages (`main/F10/red.txt`) | 4.3 s |
| Main green | Repository after the code was applied | `..............`: 14 pass, 0 failures (`main/F10/green.txt`) | 4.2 s |

- Main-session walls are from `main/F10/wall.txt` ("red wall 4.3 s; green wall 4.2 s").
- `main/F10/patch_equal.txt` is empty. The committed diff therefore equals the gated patch `patches/F10.patch`.
- `main/F10/diffstat.txt` matches the commit `--stat`.

#### Expected values

Line numbers refer to `tests/testthat/test-grf-count-orientation.R` at `f2e8de47`.

| Expected value | Lines | How it was derived |
|---|---|---|
| Subgroup members are `treat.recommend == 0L` (definition used by every assertion) | :93-94 | Roxygen `@return`, `R/grf_subg_harm_glm.R` :142-144 ("`0` = in harm/questionable subgroup"). Comment at :93. |
| `sg.harm.id` is not NULL, TRUE and FALSE fits, both blocks | :103-104, :129-130 | By hand: `dmin.grf` defaults to 0 (signature :296), and one region has a positive harm score in either orientation. Header :35-37 says the policy tree "assigns control (the harm action) to one whole region". **Verifier: no inline comment says how (minor).** |
| `sum(members) > 0`, TRUE and FALSE fits, both blocks | :105-106, :131-132 | Same reasoning as the row above. **Verifier: no inline comment (minor).** |
| `adverse_outcome = TRUE`: share of members in `z == 1` is >= 0.8 | :110 (tree + log), :134 (frontier + identity) | Direction from the roxygen of `adverse_outcome` (:100-103): TRUE means harm is "where treatment increases the outcome", which by design is `z == 1` (header :13-16, :18-23). The 0.8 margin is worked out by hand in the header (:28-53, values below). Inline comment at :108-109. |
| `adverse_outcome = FALSE`: share of members in `z == 0` is >= 0.8 | :112, :135 | Same roxygen: with FALSE "higher values are better", so harm is where treatment lowers the count, which is `z == 0` (header :24-25). Comment at :111. |
| TRUE and FALSE memberships differ (`!identical(treat.recommend)`) | :114-115, :136-137 | By hand: the two expected subgroups lie in opposite regions (header :26). Comment at :113. |
| Header margin numbers: E[log(Y + 0.5)] = 2.197 / 1.101 / 0.169 for mu = 9 / 3 / 1; SD at most 0.71; forest-scale CATE +1.10 / -0.93; DR SD about 1.4; SE = 1.4/sqrt(125) = 0.13, about 7 SE; identity scale 6 / -2 / mixed 2 / band floor 0.9 x 6 = 5.4 | :31-53 | By formula, `sum(dpois(y, mu) * log(y + 0.5))`, and Poisson means. The verifier re-derived them independently (`verify/F10/handcalc.txt`: 2.1965, 1.1005, 0.1689; SD 0.3450, 0.6095, 0.7107; CATE 1.096 / -0.932; SE 0.127; 7.32 SE; 6, -2, 2, 5.4). |

- The data constants (n = 1000, mu = 3 / 9 / 1, seeds `20261112L`) are design inputs, not expected values. The header (:55-56) says the seeds were fixed once.
- The verifier robustness check is a diagnostic, not used for tuning (`verify/F10/robust.txt`). It used 15 fresh seeds on both paths.
  - Fixed tree: "all pass: TRUE ; min shareA 0.994 ; min shareB 0.992".
  - Start tree: "all pass: FALSE ; min shareA 0.000", and the memberships are never different.

#### Unaffected case

- **What the spec names:** `adverse_outcome = FALSE` for counts, and binary and continuous outcomes, give identical output.
- **Implementer run** (`evidence/F10/unaffected.R`, `unaffected.txt`):
  - Trees: the start export (`76a62a71`) and the fixed scratch tree, each run in a separate `Rscript --vanilla` process with a sequential plan and `set.seed(4242L)` before each call.
  - Comparison: `identical()` after dropping timing leaves.
  - 12 cases:
    - direct calls `count_F_tree_log`, `count_F_tree_identity`, `count_F_frontier_log`;
    - `binary_T_tree`, `binary_F_tree`, `binary_T_frontier`;
    - `cont_T_tree`, `cont_F_tree`, `cont_F_frontier`;
    - `forestsearch()` cases `fs_binary_use_grf`, `fs_count_F_use_grf`, `fs_cont_grf_method`.
- **Result:** every case `identical = TRUE`; "ALL IDENTICAL: TRUE", exit status 0.
- **Walls:** start 14.94 s, fixed 14.81 s, compare 5.01 s.
- **Verifier re-run** (`verify/F10/unaff/`):
  - Same 12 cases: all `identical = TRUE`, "ALL IDENTICAL: TRUE". Walls: start 15.20 s, post 15.26 s, compare 5.05 s.
  - Independent check (`my_unaff.R`, new data and seeds): 7 cases `identical=TRUE`. The positive control `POSCTRL_cnt_T_tree_log` gives `identical=FALSE` (`my_compare.txt`), which shows the comparison can detect a change. No wall clock is recorded for this check.

#### Existing tests

- **Implementer call trace** (`evidence/F10/existing_trace.log`, 144 traced calls; the 46-file survey is described in the notes). Five files reach `outcome_type = "count", adverse_outcome = TRUE`, all on tree + log:

  | File | Calls |
  |---|---|
  | `test-binary-default-or.R` | 2 |
  | `test-cross-outcome-parity.R` | 13 |
  | `test-glm-pipeline-integration.R` | 1 |
  | `test-no-subgroup-cv.R` | 15 |
  | `test-return-shape-contract.R` | 3 |

  The verifier's trace (`verify/F10/trace/calls.log`) gives the same 5 files with the same counts.
- **Implementer run on the fixed tree** (`evidence/F10/existing.txt`): 16 files, all `rc=0` with 0 failures and 0 errors. Per-file walls run from 4.5 s to 45.6 s; the notes give a 45.66 s batch wall.
  - The five files above.
  - `cv-no-subgroup-dina-grf`, `cv-no-subgroup-edges`, `directive-c`, `get-fsdata-diagnostics`.
  - `no-subgroup-dina-grf-summary`, `threshold-pair-directive-a`, `grf-inclusion-band`, `no-subgroup-bootstrap`.
  - `no-subgroup-dina-grf-plot`, `no-subgroup-unified-contract`, `fs-family-report`.
- **Three files skipped for missing `dev/` or `quarto/` baselines.** They were re-run on a full copy of the start export with only `R/grf_subg_harm_glm.R` patched. All passed with no Skipped block: `binary-default-or` 10.7 s, `directive-c` 8.9 s, `threshold-pair-directive-a` 7.5 s.
- **Verifier run on its post tree** (`verify/F10/existing/summary.txt`): 31 files.
  - 30 files: `failed=0 err=0`.
  - `mr-inference`: "EXIT 124 WALL 120.053307660 s" (the 120 s cap).
  - With baselines added (`summary_withbaselines.txt`): `directive-c` nb=182, `threshold-pair-directive-a` nb=86 and `binary-default-or` nb=66, each with 0 failed and 0 skipped.
- **Existing expectations updated:** none.

#### What a user sees

From `evidence/F10/example.txt`: a direct `grf.subg.harm.glm(outcome_type = "count", adverse_outcome = TRUE, verbose = TRUE)` call on the test data. In this data, treatment raises the count in `z == 1`.

Before (`76a62a71`):

```
[grf.subg.harm.glm] Subgroup found: n(H) = 497 (49.7%) | definition: (x1 <= 2.1027522972419659 & z <= 0) | (x1 > 2.1027522972419659 & x2 <= 0.23570848070085049)
n(H) = 497 | share of H in z == 1 = 0.00 | IRR(H) = 0.290 (GLM)
```

After (fixed tree):

```
[grf.subg.harm.glm] Subgroup found: n(H) = 500 (50.0%) | definition: z > 0
n(H) = 500 | share of H in z == 1 = 1.00 | IRR(H) = 2.815 (GLM)
```

- Before, the "harm" subgroup was the benefit region, with IRR 0.29. After, it is the region where treatment raises the adverse count.
- The verifier reproduced both outputs (`verify/F10/example/example_out.txt`).

#### Verification

- **Verdict:** CONFIRMED.
- **Reproduced:** red, green and the unaffected case.
- **Re-derived:** the expected values, independently.
- **No repair round.**

| Problem (verifier) | Severity | Resolution |
|---|---|---|
| The existence assertions (test :103-106, :129-132) have no inline comment on how they were derived. The reason is in the header and the notes, not at the assertions. | minor | Not changed. There is no `patches/F10.verified.patch`, and `main/F10/patch_equal.txt` is empty, so no comment was added before the commit. Left as an open item. |
| The NEWS hand-off named only the `forestsearch()` routes. It left out the GRF arm of `run_simulation_analysis()`, which resolves `adverse_outcome` to TRUE for count (start `R/run_simulation_analysis.R` :1257-1258) and dispatches to `grf.subg.harm.glm` (start :1196). | minor | Fixed in the NEWS commit `88302212`. `NEWS.md` :361-365 names "the GRF arm of `run_simulation_analysis()`". |
| The GRF CATE footnote in `R/summary_utility_functions.R` :913-918 is labelled "(risk-difference scale)". For count with adverse TRUE its values are now on the negated scale. The label was already wrong for counts at the start SHA. | minor | Pre-existing. Not changed under F10. Left as an open item. |

#### Departures and interpretations

- **Header comment wording.** Besides adding count to the list, "positive CATE = treatment increases the BAD outcome" became "a higher Y_grf is better and positive -CATE = ...". On the flipped scale harm has negative CATE, and the old wording contradicted the next sentence ("most-negative-CATE subgroup is the true harm group"). This is a comment-only change. The verifier judged it correct.
- **"Negate the transformed outcome".** This was read as negating whichever count transform applies, `log(Y + 0.5)` or identity. The negation is placed after the transform (start :409-434) and before `causal_forest()`. The start :451 comment was replaced by a comment inside the new branch (head :453-455).
- **"Mostly".** Read as a share of at least 0.8 of the identified subgroup's members in the expected region. The margin is argued in the test header.
- **"Memberships differ".** Read as `!identical(treat.recommend)` between the TRUE and FALSE fits. Members are `treat.recommend == 0` (roxygen @return).
- **Test scope.** The spec's scenario is run on two paths, tree + log and frontier + identity, instead of one, with the same assertions. The FALSE orientation (mostly `z == 0`) is also checked. This shows that the "differ" check compares two correctly oriented answers.
- **`skip_on_cran()`.** Placed at the top of the test file because the file fits four forests, as `test-grf-inclusion-band.R` does.
- **Not added.** There is no verbose message for the count flip, because the spec did not ask for one. There is no roxygen change, so `document()` was not run.

#### paper_simulations

Not reached. The only evaluated package-calling chunks are `run-alt-sims` and `run-null-sims`. They use survival DGMs, so the GRF arm calls `grf.subg.harm.survival`. The changed line also needs `outcome_type == "count"`.

#### Open items

- **No verbose message for the count flip.** At head `R/grf_subg_harm_glm.R` :452-457 the count flip prints nothing. The binary and continuous flips print one (head :448, :451). With `verbose = TRUE` a count user sees "applying log(Y + 0.5)" but nothing about the sign. Not fixed because the spec did not ask for a message.
- **`grf_y_transform` does not record the sign flip.** The field is set at head :407, :412 and :423 and documented at :157-159. It still reports `"log(Y + 0.5)"` or `"identity"`. This matches binary and continuous, which report `"none"` although `1-Y` or `-Y` is used. Not fixed: outside the spec.
- **Roxygen examples use the default orientation.** The example calls at start :243-256 and :260-272 use the default `adverse_outcome = FALSE`. The audit's optional clarity edit (an explicit `adverse_outcome = TRUE`) was not applied because it is not in the spec.
- **Missing derivation comments in the test.** The existence assertions at `f2e8de47` :103-106 and :129-132 of `tests/testthat/test-grf-count-orientation.R` have no inline comment. The verifier flagged this as minor, and no comment was added before the commit.
- **GRF CATE footnote label.** `R/summary_utility_functions.R` :917 (unchanged between start and head) labels the CATE footnote "(risk-difference scale)". For count fits that label is wrong, and with adverse TRUE the sign of the printed values is now reversed relative to before the fix. The issue is pre-existing and not in F10's scope.
- **NEWS wording.** `NEWS.md` :361-362 at `88302212` attaches "(the default for counts)" to `grf.subg.harm.glm()`. A direct call defaults to `adverse_outcome = FALSE` (`R/grf_subg_harm_glm.R` :304). TRUE is the resolved default only in `forestsearch()` (`R/forestsearch_main.R` :1840-1841) and `run_simulation_analysis()` (head :1303-1304). Wording only; not changed here.
- **`test-mr-inference.R` was not run to completion.** It hit the 120 s cap in both the implementer's survey and the verifier's run (EXIT 124). It uses `use_grf = FALSE` with survival and binary outcomes, so it cannot reach `grf.subg.harm.glm()`.

#### Corrections to the notes

- **Survival branch line.** The notes and the verifier place the survival branch of the `use_grf` route at `R/forestsearch_main.R` :3124. At `76a62a71`, `if (outcome_type == "survival") {` is at :3125; :3124 is a comment. At head `61c15088` it is at :3124, and the `grf.subg.harm.glm` call moves from :3161 to :3160.
- **Existence-assertion lines.** The verifier gives them as ":102-105 and :126-129". At `f2e8de47` they are :103-106 and :129-132.
- **Roxygen example range.** The notes cite the roxygen example as "start :235-256". That range covers the outcome generation and the first call. A second call at :260-272 also uses the default `adverse_outcome`.
- **Transform range.** The notes give the count transform as :409-433. At `76a62a71` the block closes at :434.
- **Batch walls.** The implementer's batch walls (142.47 s trace survey, 45.66 s for 16 files) appear only in `notes.md`, not in `existing.txt`. The per-file walls in `existing.txt` match the notes.

### F11 — DINA orientation sign on the screening route

**Status:** PASS. Commit `61c15088` (start SHA `76a62a71`). R/ and the test file are unchanged from `61c15088` to the final head `88302212` (`git diff --stat 61c15088 88302212 -- R/ tests/testthat/test-dina-screening-tau-sign.R` is empty).

| File (`git show --stat 61c15088`) | Change |
|---|---|
| `R/forestsearch_main.R` | 3 +- |
| `tests/testthat/test-dina-screening-tau-sign.R` (new) | 132 + |
| Total | 2 files, 134 insertions, 1 deletion |

- No roxygen, `man/`, `NAMESPACE`, `DESCRIPTION` or `NEWS.md` change in this commit.
- `main/F11/patch_equal.txt` is empty: the committed diff equals the gated `patches/F11.patch` apart from index and hunk-header lines.

#### Site, before and after

The call sits inside `if (use_dina && is.null(dina_cuts))` (start :3228, head :3227) and `if (isTRUE(da$selected_only))` (start :3273, head :3272). The spec's "near 3254–3265" was taken at `28e71ffa`.

Before, `R/forestsearch_main.R:3299-3303` at `76a62a71` (call opens at :3292):

```r
          sg_focus            = sg_focus,
          selection_rule      = selection_rule,
          effect_neighborhood = effect_neighborhood
        )
        if (isTRUE(sgsel$found)) {
```

After, `R/forestsearch_main.R:3298-3303` at `61c15088` (call opens at :3290; added line :3301):

```r
          sg_focus            = sg_focus,
          selection_rule      = selection_rule,
          effect_neighborhood = effect_neighborhood,
          tau_sign            = .dina_tau_sign(outcome_type, adverse_outcome)
        )
        if (isTRUE(sgsel$found)) {
```

- The route it copies is unchanged: `R/forestsearch_helpers.R:1542` `tau_sign <- .dina_tau_sign(outcome_type, adverse_outcome)` and `:1604` `tau_sign            = tau_sign` (same line numbers at both SHAs).
- The rule is `.dina_tau_sign()` at `R/forestsearch_helpers.R:1413-1420`: -1 when `!isTRUE(adverse_outcome)` and the outcome is continuous or binary, else 1.
- `adverse_outcome` is resolved before this block at `R/forestsearch_main.R:2241-2243` (NULL becomes `outcome_type %in% c("binary", "count")`).

#### Red and green

Command: `devtools::test(filter = "dina-screening-tau-sign", reporter = "summary")` with `NOT_CRAN=true`. The test file has 16 expectations. No run had an error; all failures are expectation failures.

| Run | Where | Result | Wall |
|---|---|---|---|
| Implementer red | `wt/F11`, git archive of `76a62a71` plus the test file | 13 pass, 3 failures. :95:5 twice, "Expected \`run$signs\` to be identical to \`cs$expected\`", actual `1.0`, expected `-1.0` (continuous, `adverse_outcome = NULL`; binary, `adverse_outcome = FALSE`). :122:3, "Expected \`cut_benefit_y\` to match regexp "^x1 <= "", actual `x1 >= -0.481` | 7.17 s |
| Implementer green | `wt/F11` with the fix | 16 pass, 0 failures | 6.92 s |
| Verifier red | `verify/F11/pre` (archive of `76a62a71`) plus the test file | 13 pass, 3 failures; same messages and lines (:95:5 twice, :122:3) | 7.22 s |
| Verifier green | `verify/F11/post` | 16 pass, 0 failures | 6.81 s |
| Main-session red | repository at `a4cce90a` (preceding commit), test hunks only | 13 pass, 3 failures; same messages, at :96:5 twice and :125:3 | 6.8 s |
| Main-session green | repository after the code hunk | 16 pass, 0 failures | 6.7 s |

- Line numbers move from :95/:122 to :96/:125 because the main session added three comment lines to the test before commit (see Verification).
- The main-session red ran with F14 (`685e0c5d`) and F3 (`a4cce90a`) already committed. It gave the same three failures.
- The implementer also had an accidental second red run (first `Edit` refused): same three failures, 7.04 s. Its output was overwritten, so this run is in the notes only.

#### Expected values

Line numbers are for the test file at `61c15088`.

| Expected value (line) | How it was derived (comment in test) |
|---|---|
| continuous, `adverse_outcome = NULL` -> `-1` (:66) | `?forestsearch` `adverse_outcome`: NULL means FALSE for continuous ("higher values indicate benefit"). `?dina_subgroup` `tau_sign`: forestsearch() "passes -1 for continuous and binary outcomes with adverse_outcome = FALSE". Comment :56-63, :65. |
| continuous, `TRUE` -> `1` (:70) | Not in the -1 set, so the default 1. Comment :60-63, :69. |
| binary OR, `FALSE` -> `-1` (:74) | In the -1 set (`?dina_subgroup`). Comment :73. |
| binary OR, `NULL` -> `1` (:78) | `?forestsearch`: NULL means TRUE for binary, so 1. Comment :58-59, :77. |
| count IRR, `FALSE` -> `1` (:82) | Count is never flipped: `.dina_tau_sign()` roxygen (`R/forestsearch_helpers.R:1406-1407`, `@noRd`), "Survival, count, and `adverse_outcome = TRUE` keep `1`". Comment :62-63, :81. |
| survival, `NULL` -> `1` (:86) | Same roxygen sentence. Comment :85. |
| `expect_length(run$signs, 1L)` (:95, once per case) | `subgroup_method = "consistency"`, so `.forestsearch_dina_select()` is not run and the screening call is the only call. Comment :94. |
| A call with no `tau_sign` is read as `1` (:32, wrapper, not an assertion) | `?dina_subgroup`: "`1` (default) or `-1`". Comment :23-25. |
| `expect_length(cut_benefit_y, 1L)` (:124), `expect_length(cut_harm_y, 1L)` (:130) | `dina_args = list(max_depth = 1L)`: one single-covariate cut, so one cut string. Comment :122-123, added by the main session. **The verifier flagged both as lacking a comment.** :130 has no comment of its own; it uses the same `screen_cut()` and the comment at :122-123. |
| `expect_match(cut_benefit_y, "^x1 <= ")` (:125) | By construction (:107): treatment changes y by -2 where `x1 <= 0` and by +2 where `x1 > 0`. With `adverse_outcome = FALSE`, higher y is benefit, so harm is at `x1 <= 0`, and a single x1 cut points at it only if it keeps the low side. Comment :101-102, :118-120. |
| `expect_match(cut_harm_y, "^x1 >= ")` (:131) | Mirror case. With `adverse_outcome = TRUE`, higher y is harm, so harm is at `x1 > 0` and the cut keeps the high side. This is the control and passes before and after. Comment :127-128. |

- The inputs `set.seed(2026L)`, `N = 300L` and `effect.threshold = 1` were fixed before the first run (notes). No seed, tolerance or expected value changed after a run.
- The verifier re-derived every value independently and found "No value appears to come from running the fixed code".

#### Unaffected case

- **Spec names:** survival and count outcomes, binary with `adverse_outcome = TRUE`, and every `subgroup_method = "dina"` fit.
- **Script:** `evidence/F11/unaffected.R`. The verifier re-ran the same script (byte-identical copy).
  - Data: the start export's `helper-synthetic-dgm.R`, N = 300.
  - Each run used a sequential plan and `set.seed(123)` before each call, with GRF and LASSO off.
  - The compare drops only `$minutes_all` and `$find.grps$time_search`, then calls `identical(ignore.environment = TRUE, ignore.bytecode = TRUE, ignore.srcref = TRUE)`.
- **Trees:**
  - Implementer: `/tmp/cran_prep_code_fixes_2026-11-12/start` (export of `76a62a71`) against `wt/F11`.
  - Verifier: `verify/F11/pre` against `post`. Both are archives of `76a62a71`; `post` has the patch.

| Cases | Covers | identical() |
|---|---|---|
| U1-U4: `use_dina`, survival; count default; count `FALSE`; binary OR `TRUE` | named | TRUE (all 4) |
| D1-D9: `subgroup_method = "dina"`, all outcome types and orientations. Includes D7-D9, where a subgroup is found under a flipped orientation (`{age <= 91}`, `{age <= 71}`, `{age <= 70} & {wtkg >= 61.67...}`) | named | TRUE (all 9) |
| X1, X3: `use_dina`, continuous `TRUE`. X2: `use_dina`, continuous, `selected_only = FALSE` | extra | TRUE (all 3) |

- Result: `UNAFFECTED GATE (all 16 cases identical): PASS` in both `evidence/F11/unaffected.txt` and `verify/F11/unaffected_compare.txt`.
- Implementer walls: start 10.24 s, fixed 10.21 s, compare 1.65 s.
  - An earlier 12-case run without D7-D9 and X3 also passed: 7.82 / 7.70 / 1.63 s. These figures are in the notes only.
- Verifier walls: in-R totals 9.73 s (pre) and 9.30 s (post), from `unaffected_run_{pre,post}.txt`. The verifier's JSON gives process walls of 9.95 s and 9.52 s.

#### Existing tests

- **Implementer, fixed scratch tree (`evidence/F11/existing.txt`):** each file ran in its own process with a 120 s cap.

| File | Result | Wall |
|---|---|---|
| test-dina-gaussian-floor.R (reaches the line) | FAIL 0, PASS 31 | 6.35 s |
| test-directive-c.R | FAIL 0, SKIP 1, PASS 154 | 8.38 s |
| test-directive-c.R + dev baseline CSV (Gate C) | FAIL 0, SKIP 0, PASS 182 | 9.39 s |
| test-dina-proposal-orientation.R | FAIL 0, PASS 25 | 3.14 s |
| test-fs-family-report.R | FAIL 0, PASS 105 | 5.43 s |
| test-dina-bootstrap-sg-focus-whitelist.R | FAIL 0, PASS 15 | 3.63 s |
| test-cv-no-subgroup-dina-grf.R | FAIL 0, PASS 48 | 14.09 s |
| test-no-subgroup-dina-grf-summary.R | FAIL 0, PASS 21 | 5.69 s |
| test-no-subgroup-dina-grf-plot.R | FAIL 0, PASS 22 | 4.89 s |

- **Implementer, 39 other files that call `forestsearch()`** (`existing_broad.txt`):
  - 38 ran to the end as whole files: FAIL 0, WARN 31, SKIP 5, PASS 1755. The wall sum over all 39 whole-file runs is 612.99 s.
  - `test-mr-inference.R` hit the cap as a whole file (exit 124, 120.05 s). It never sets `use_dina`.
  - Its 21 `test_that()` blocks were then run one per process (`existing_mr_inference_blocks.txt`): FAIL 0, PASS 148, all exit 0, wall sum 384.10 s, longest block 68.36 s.
- **Verifier, post tree:** 0 failures in every file run.
  - The same 8 files (`existing_post.txt`).
  - Directive C with the dev baseline: 182 (`existing_post_dev.txt`).
  - 12 more files that mention DINA (`existing_post_more.txt`).
  - 9 broad `forestsearch()` files (`existing_post_broad.txt`).
- Only `test-dina-gaussian-floor.R` and `test-directive-c.R` set `use_dina = TRUE` (`git grep` at `76a62a71`). Neither asserts on the sign.
- **Existing expectations updated:** none.
- The broad runs left a stray `tests/testthat/Rplots.pdf` in the scratch tree. It was deleted before the patch and is not in the commit.

#### What a user sees

Made-up continuous data from test 2. The call is `forestsearch(outcome_type = "continuous", effect_measure = "MD", effect.threshold = 1, sg_focus = "maxSG", use_dina = TRUE, use_grf = FALSE, use_lasso = FALSE, dina_args = list(max_depth = 1L), details = TRUE)`, with the default orientation. Output from `evidence/F11/example.txt`:

```
# start (76a62a71), WALL 3.52 s
DINA cuts identified: 1
  Cuts: x1 >= -0.481
fs$dina_cuts: x1 >= -0.481
fs$sg.harm:   {x1 <= 0}

# fixed, WALL 3.34 s
DINA cuts identified: 1
  Cuts: x1 <= 0.352
fs$dina_cuts: x1 <= 0.352
fs$sg.harm:   {x1 <= 0.352}
```

- Before the fix, the proposal pointed at the benefit side. After the fix, it points at the harm region.
- `sg.harm` also changes. Under `maxSG`, the harm-side cut joins the pool and gives a larger consistent subgroup.
- The verifier reproduced both outputs exactly (`verify/F11/my_example_{pre,post}.txt`, `impl_default`).

#### Verification

- **Verdict:** CONFIRMED. The verifier reproduced red, green and the unaffected case.
- **Minor: loose source for the "keeps 1" rule.** The test-1 comment credited it to `?dina_subgroup`, but that sentence is in the `@noRd` roxygen of `.dina_tau_sign()`.
  - Resolved by the main session before commit, comment only. The test now cites both sources (:60-63 at `61c15088`; diff `patches/F11.verified.patch` vs `patches/F11.patch`).
- **Minor: no derivation comment on the two `expect_length(..., 1L)` calls in test 2.**
  - Resolved by the main session before commit, comment only. It added `# max_depth = 1L: dina_subgroup() selects one single-covariate cut, so one` / `# cut string.` above the first one (:122-123). This is what the verifier suggested.
  - :130 shares that comment.
- The two edits are the only differences between the verified and the committed patch, apart from the index and hunk header. The code hunk is identical.
- After the edits, the main session re-ran red and green on the repository (table above).

#### Departures and interpretations

- **Code:** no departure. The commit adds exactly the one argument the spec names. No roxygen changed, so `document()` was not run.
- **Site location:** the spec's "near 3254–3265" was found by content. It is `R/forestsearch_main.R:3292-3303` at `76a62a71` and `:3290-3302` at `61c15088`.
- **How test 1 observes the sign:** it reads the effective sign from a `local_mocked_bindings` wrapper. That is the passed `tau_sign`, or 1 (the documented default) when none is passed. So the four sign-1 cases pass before and after the fix. They are controls against over-flipping, and red fails only where the effective sign is wrong.
- **"Points at that region" (test 2):** read as a single x1 cut (`max_depth = 1L`) keeping the harm side, checked by the operator (`^x1 <= ` / `^x1 >= `), not by a cut value. The `adverse_outcome = TRUE` mirror is the control.
- **"Every `subgroup_method = "dina"` fit":** read as each outcome type and orientation on that route (D1-D9), including cases where a subgroup is found. D7-D9 were added after a first 12-case run. X1-X3 are extras beyond the spec.
- **Survival:** the screening route resolves NULL to FALSE for survival, while the `dina` route passes TRUE (`R/forestsearch_main.R:2816-2817`). `.dina_tau_sign()` returns 1 for survival either way, so the two routes agree.
- **Existing tests:** `test-mr-inference.R` was run block by block after the whole file hit the 120 s cap.

#### paper_simulations

- Does not reach the changed line. The article has no `dina` or `use_dina` (0 case-insensitive matches at `76a62a71`).
- Its evaluated `run_simulation_analysis()` calls (chunks `run-alt-sims`, `run-null-sims`) pass no `methods`. `fs_params` and `default_sim_params()` leave `use_dina = FALSE`, and the outcome is survival (sign 1 anyway).

#### Open items

- **`selected_only = FALSE` pool not oriented.**
  - Where: `R/forestsearch_main.R:3315` at `76a62a71` (`:3314` at `61c15088`), `fr$cut_expr`.
  - What: with `dina_args = list(selected_only = FALSE)`, the screening pool is `dina_frontier()`'s `cut_expr`, which reads the raw tau-hat.
  - Why not fixed here: this is A-020 in the 2026-11-11 report, and F11's spec says nothing else on either route changes.
- **`dina_select_statistic = "effect"` re-selection only on the `dina` route.**
  - Where: `R/forestsearch_helpers.R:1611-1624`.
  - What: the screening route keeps the tau-hat winner. So the roxygen sentence "the same cut `subgroup_method = "dina"` would select" (`R/forestsearch_main.R:440-445`) holds only for the default statistic.
  - Why not fixed here: out of F11's scope.
- **Quarto smoke test not re-rendered.**
  - Where: `quarto/dina/smoke-tests/dina_screening_smoketest.qmd` runs a continuous battery with `use_dina = TRUE` (:112, :126, :134; continuous battery at :169).
  - What: its DINA cuts and pass/fail checks may change.
  - Why not fixed here: `quarto/` is read-only in this task. Per the implementer and the verifier, the other quarto files that use `use_dina` or `fs-*` labels are survival, or binary with `adverse_outcome` TRUE or the default.
- **`sg.harm` can change, not only `dina_cuts`.** This is documented behaviour of `maxSG` (see the example). Addressed in the separate NEWS commit: `NEWS.md:367-372` at `88302212` says "the candidate pool and the selected subgroup can change".
- **`test-mr-inference.R` exceeds the 120 s cap as a whole file.** This was under CPU contention from parallel agents. It does not set `use_dina`, and all 21 blocks passed when run separately. Noted only; not an F11 defect.

#### Corrections to the notes

- The notes' "After" range `:3292-3304` is the scratch tree (start + F11 only). In the repository at `61c15088` the call is `:3290-3302`, with the added line at `:3301`, because the F14 commit (`685e0c5d`) removed one line in each of two hunks above it (`@@ -2818,10 +2818,9 @@`, `@@ -3281,13 +3280,12 @@`).
- The notes' failure lines `:95` / `:122` and diffstat (129 lines; 131 insertions) are for the pre-comment test file. The committed file has 132 lines, the failures at `:96` / `:125`, and 134 insertions.
- The notes' `R/forestsearch_helpers.R:1611-1625`: the `if` block ends at `:1624`; `:1625` is blank.
- The notes' paper_simulations `:908` / `:966` are the chunk headers. The `run_simulation_analysis(` calls are at `:925` / `:983`.
- The verifier's unaffected walls (9.95 s / 9.52 s) are not in a file. The files record in-R totals of 9.73 s / 9.30 s.

### F12 — Three displays use the fit's own orientation

**Status:** PASS. Commit `25ecaa0f` (parent `f2e8de47`). Files (`git show --stat 25ecaa0f`):

| File | Lines (`--stat`) | `--numstat` |
|---|---|---|
| `R/mr_estimates_table.R` | 8 | +8 / -0 |
| `R/plot_sg_glm_outcomes.R` | 9 | +9 / -0 |
| `R/plot_subgroup_results_forestplot.R` | 13 | +13 / -0 |
| `tests/testthat/test-display-orientation.R` (new) | 227 | +227 / -0 |

- 4 files changed, 257 insertions(+). No roxygen, `man/`, `NEWS.md` or `DESCRIPTION` change.
- The four files are unchanged from `25ecaa0f` to head `61c15088` (`git diff --stat` is empty). They are also unchanged from the start SHA `76a62a71` to the parent `f2e8de47`. `R/` is unchanged from `61c15088` to `88302212`, so the head line numbers below also hold at `88302212`.
- The commit diff equals the gated patch `patches/F12.patch`, apart from index and hunk-header lines. `main/F12/patch_equal.txt` is empty, and a re-check with `diff` against `git show 25ecaa0f` found no difference.

#### Site, before and after

Each display built its estimator without `adverse_outcome`, so `make_effect_estimator()`'s formal `adverse_outcome = TRUE` applied (`R/glm_effect_estimators.R:107` at `76a62a71`). The fix applies `forestsearch()`'s rule at `R/forestsearch_main.R:1840-1841` (`76a62a71`), which the GLM bootstrap fix also uses (`R/bootstrap_dofuture_main.R:306-309`):

```r
  if (is.null(adverse_outcome)) {
    adverse_outcome <- (outcome_type %in% c("binary", "count"))
```

**`plot_sg_glm_outcomes()`.** At the start, `R/plot_sg_glm_outcomes.R:163-166` and the call at `:254-262` (`:259-261` shown):

```r
    adjust_covariates <- fs.est$args_call_all$adjust_covariates
  }

  is_log_scale <- effect_measure %in% c("OR", "RR", "IRR")
...
      effect_measure = effect_measure,
      adjust_covariates = adjust_covariates,
      ps_adjust_method = "none"
```

At head, the new block is at `:166-172` (shown `:164-174`), and the call at `:262-271` gains `:268`:

```r
  }

  # Orientation: the fit's resolved adverse_outcome, so the annotated effect
  # is the fit's own estimate, sign included.  NULL only on an object not
  # built by forestsearch(); resolved by forestsearch()'s rule.
  adverse_outcome <- fs.est$args_call_all$adverse_outcome
  if (is.null(adverse_outcome)) {
    adverse_outcome <- outcome_type %in% c("binary", "count")
  }

  is_log_scale <- effect_measure %in% c("OR", "RR", "IRR")
...
      effect_measure = effect_measure,
      adverse_outcome = adverse_outcome,
      adjust_covariates = adjust_covariates,
```

**`plot_subgroup_results_forestplot()`.** At the start, `R/plot_subgroup_results_forestplot.R:309-312`. The call is at `:432-440` (`:437-439` shown), inside `create_glm_row()` (`:392`), which builds every GLM row:

```r
    adjust_covariates <- fs.est$args_call_all$adjust_covariates
  }

  is_glm <- outcome_type != "survival"
...
      effect_measure = effect_measure,
      adjust_covariates = adjust_covariates,
      ps_adjust_method = "none"
```

At head, the new block is at `:312-322` (shown `:310-324`), and the call at `:444-453` gains `:450`. `create_glm_row()` is now at `:404`.

```r
  }

  # Orientation of the GLM rows: the fit's resolved adverse_outcome, so each
  # row's effect is on the fit's own scale, sign included.  NULL only on an
  # object not built by forestsearch(); resolved by forestsearch()'s rule.
  # Without a fit (fs.est = NULL) the estimator's default (TRUE) is kept.
  adverse_outcome <- TRUE
  if (!is.null(fs.est)) {
    adverse_outcome <- fs.est$args_call_all$adverse_outcome
    if (is.null(adverse_outcome)) {
      adverse_outcome <- outcome_type %in% c("binary", "count")
    }
  }

  is_glm <- outcome_type != "survival"
...
      effect_measure = effect_measure,
      adverse_outcome = adverse_outcome,
      adjust_covariates = adjust_covariates,
```

**`mr_estimates_table()` with `boot_results = NULL`.** This route goes through `.mr_table_from_fs()` (`:174`, called at `:81-82`). At the start, `R/mr_estimates_table.R:203-208`:

```r
  } else {
    est_fn <- make_effect_estimator(
      ...
      offset.name = a$offset.name, effect_measure = a$effect_measure,
      adjust_covariates = a$adjust_covariates)
```

At head, the new block is at `:204-210` and the call at `:211-216` gains `:215`:

```r
  } else {
    # Orientation: the fit's resolved adverse_outcome, so the effect cells are
    # the fit's own estimates, sign included.  NULL only on an object not
    # built by forestsearch(); resolved by forestsearch()'s rule.
    adverse_outcome <- a$adverse_outcome
    if (is.null(adverse_outcome)) {
      adverse_outcome <- outcome_type %in% c("binary", "count")
    }
    est_fn <- make_effect_estimator(
      ...
      offset.name = a$offset.name, effect_measure = a$effect_measure,
      adverse_outcome = adverse_outcome,
      adjust_covariates = a$adjust_covariates)
```

Not changed, as the spec says: `R/summary_utility_functions.R`, where the arm means and `Diff` stay raw. The plot's per-arm means also stay raw.

#### Red and green

Test file: `tests/testthat/test-display-orientation.R`, with 6 blocks and 93 expectations.

| Run | Where | Result | Wall |
|---|---|---|---|
| Implementer red | Scratch copy `wt/F12` at `76a62a71`, test file only | 93 expectations, **24 failed**, 0 errors, 0 skipped. By block (nb/failed): 7/0, 9/6, 10/3, 8/3, 28/12, 31/0. Every failure is a sign mismatch. Example (`evidence/F12/red.txt:88-91`): "Expected forest est H to equal `ref[[k]]$est`." `actual`: -1.1, `expected`: 1.1 | 7.28 s |
| Implementer green | Same copy, code patch applied | 93 expectations, 0 failed, 0 errors, 0 skipped. All 6 blocks 0 failed | 7.08 s |
| Verifier red | `verify/F12/pre` (`git archive 76a62a71` + test) | 93, **24 failed**, 0 errors, 0 skipped. 24 "Failure" entries and 0 "Error" entries. Same per-block counts and messages, e.g. `pre_test.txt:108`: "Expected MR sign Hc to be identical to `sign(want[1])`." actual -1.0, expected 1.0 | 7.6 s in R, 7.78 s process |
| Verifier green | `verify/F12/post` | 93, 0 failed, 0 errors, 0 skipped | 6.68 s in R, 6.84 s process |
| Main-session red | Repository at `f2e8de47` (the parent), test file only | The progress line has 93 marks: 69 passes and 24 not passed (marks `1`-`9`, `a`-`o`). It is identical to the implementer's and the verifier's red line. 10 are listed, all "Failure", then "... and 14 more" ("Maximum number of 10 failures reached"). Example (`main/F12/red.txt:88-91`): "Expected forest est H to equal `ref[[k]]$est`." actual -1.1, expected 1.1 | 6.8 s |
| Main-session green | Repository with the code hunks applied (the tree committed as `25ecaa0f`) | 93 passes, "DONE", no failures | 6.4 s |

- The implementer's first red run stopped listing at 10 failures (7.54 s, per the notes). It was re-run identically with `options(testthat.summary.max_reports = Inf)`, and `red.txt` is that second run.
- Blocks 1 (the fixture) and 6 (`adverse_outcome = TRUE`) pass on both trees by design.

#### Expected values

All values are in `tests/testthat/test-display-orientation.R` at `25ecaa0f`. The test's header comment (`:16-30`) gives the derivation.

| Expected value (test lines) | How it was derived |
|---|---|
| Hand MD estimate for ITT, H and Hc: `mean(y \| T = 1) - mean(y \| T = 0)`, negated when `adverse = FALSE` (`.dor_md()`, `:78-84`) | By hand: the unadjusted OLS coefficient of `y` on a 0/1 treatment. Orientation is from `?forestsearch` `@param adverse_outcome` (`R/forestsearch_main.R:1104-1107`: TRUE means higher values indicate harm; FALSE is the "default for continuous"), with the convention that effect > 0 means harm. Comments are at `:16-24` and `:77`. |
| Hand pooled SE `sqrt(s2 * (1/n0 + 1/n1))`, `s2 = (SS0 + SS1)/(n0 + n1 - 2)` (`:81`, `:83`) | By hand (header `:17-19`). The standard error is the same under both orientations (`:24`). |
| Plot label = estimate and limits `est -/+ z*se`, `z = qnorm(0.975)` (`:99`, `:105-106`). Within `5e-4` (`:108`). 3 numbers per label (`:107`) | z is from `conf.level = 0.95`, the documented default (comment `:99`, formal `R/plot_sg_glm_outcomes.R:109`). The tolerance is half a unit at 3 printed decimals (`:27-28`). |
| Label names `ITT`, `Questionable`, `Recommend` (`:102`) | The documented defaults `sg0_name`/`sg1_name` (`R/plot_sg_glm_outcomes.R:50-53`). Commented "default names". |
| Forest plot `est` and `se` of the ITT, Hc and H rows within `1e-10` (`:127-130`). One row per type (`:126`). `length(row_types) == nrow(data)` (`:120`) | By hand, as above. The row-type locators `"itt"`, `"posthoc"` = Hc and `"posthoc_complement"` = H are structural and commented at `:121-123`. **The verifier flagged that the `1e-10` tolerance (unrounded columns) is explained only in `notes.md`, not in the test.** |
| MR table rows `c("Qstnbl", "Recmnd")` (`:138`) | A locator taken from the code (`sg0_name = "Qstnbl", sg1_name = "Recmnd"`, `R/mr_estimates_table.R:202`). It has no comment at the line, and the verifier counted it among the documented defaults. |
| MR table `MD (95% CI)` estimate within `5e-5`, sign, and 3 numbers (`:142-146`). Hc limits `est -/+ 1.96*se` within `1e-4` (`:150-152`) | By hand. 1.96 is the z the table's descriptive builder uses (`analyze_subgroup_glm()`, `R/summary_utility_functions.R:223` `z <- 1.96`). The H row's limits are not checked: `?mr_estimates_table` documents them as `fs$mr_inference$naive`'s robust (sandwich) limits (comment `:147-149`). **The verifier flagged the comment at `:29-30`, which justifies `1e-4` by a 1.96 vs `qnorm(0.975)` difference; the expected limits also use 1.96, so the real bound is double rounding at 4 decimals.** |
| Sign of each displayed estimate equals the sign of the hand estimate (`:109`, `:146`) | Follows from the hand estimate (the invariant, `:12-13`). |
| `fit$args_call_all$adverse_outcome` is `FALSE` (`:161`) | `forestsearch()`'s rule (comment `:161`), as documented in `?forestsearch` ("default for continuous"). |
| `fit$max_sg_est` and `fit$mr_inference$naive$est` equal the hand H estimate within `1e-10` (`:168`, `:170`). For the TRUE fit, `max_sg_est` within `1e-10` (`:220`) | By hand (`.dor_md()`). This ties the hand value to the fit's own estimate, independently of the displays (`:25-26`, `:165-166`). |
| Hand H estimate `> 0` (`:169`) | The planted region raises `y` by 1.5, and `ben` negates `y`, so under FALSE the region reads as harm (comments `:35-38`, `:165-166`). |
| `fit$mr_inference$complement` is NULL (`:164`) | `include_complement = FALSE` in the fixture (comment `:164`, argument `:58`). |
| `sg.harm` and `mr_inference` are not NULL (`:162-163`, `:218`) | Fixture assumptions: a subgroup is found, and `mr_inference = TRUE`. |
| `"adverse_outcome"` is absent from `args_call_all` (`:204`) | Structural: the block removes it at `:203`. |
| TRUE fit `adverse_outcome` is `TRUE` (`:217`) | Passed explicitly (`:71`). |
| `fit_f$sg.harm` is identical to `fit_t$sg.harm` (`:226`) | Mirror twins: the same data with `y` negated and the opposite orientation (comment `:224-225`). |

- The verifier rederived the values in base R with no package code (`verify/F12/hand.R`, `hand.txt`). Oriented MD (95% CI), se:
  - ITT: 0.2287 (-0.0070, 0.4644), se 0.1203
  - H: 1.0941 (0.5232, 1.6651), se 0.2913. `lm` on `-y` in H gives 1.094112.
  - Hc: 0.0958 (-0.1464, 0.3380), se 0.1236
- The verifier found no expected value copied from the fixed code's output.

#### Unaffected case

- **Spec:** survival fits, and non-survival fits with `adverse_outcome = TRUE`, give identical display data.
- **Run:** `evidence/F12/unaffected.R`, one `Rscript --vanilla` process per tree (`pkgload::load_all(tree)`, sequential plan, `set.seed(2026)` before each call). Start export vs fixed tree, plus a second start run as a control. Cases (`evidence/F12/unaffected.txt`):

| Case | Selected subgroup | `adverse_outcome` | MR |
|---|---|---|---|
| `survival` (MR on, complement on) | `{age <= 63} & !{stage}` | FALSE | TRUE |
| `continuous_true` (TRUE) | `{age <= 47} & {biomarker_hi}` | TRUE | TRUE |
| `binary_OR` (unset) | `{age <= 55} & {biomarker_hi}` | TRUE | TRUE |
| `binary_RD` (TRUE) | `{age <= 55} & {biomarker_hi}` | TRUE | TRUE |
| `count` (IRR with offset, unset) | `{age <= 47} & {biomarker_hi}` | TRUE | TRUE |
| `forest_no_fit` (`fs.est = NULL`) | - | - | - |

- **Saved per case:**
  - the fit, without timing leaves;
  - the forest plot's `data`, `row_types` and `cv_metrics`;
  - `plot_sg_glm_outcomes()`: every field except `plot`, plus `ggplot_build()` layer data and the labels. For survival this is the error message;
  - `mr_estimates_table()`: `_data`, `_boxhead`, `_heading`, `_footnotes`, `_source_notes`, `_spanners`, the `mr_values` attribute and the raw HTML.
- **Result, start vs fixed:**
  - `survival` and `forest_no_fit`: `identical()` = TRUE.
  - The four non-survival cases: `identical()` = FALSE. Each differs only at `<case>$glm_plot$labels (attributes)`.
  - The control (start vs start) reports exactly the same four leaves. The leaf is the `S7_class` attribute of ggplot2's labels object, a constructor closure whose environment differs after each `readRDS()`.
  - With that attribute dropped, **all six cases are `identical()` = TRUE**, overall TRUE, in both comparisons.
- **Verifier:** reproduced the same output (`verify/F12/unaff/cmp_start_post.txt`, `cmp_control.txt`).
- **Wall:**

| Run | Start | Fixed | Start (control) | Compare start vs fixed | Compare control |
|---|---|---|---|---|---|
| Implementer | 10.26 s | 10.55 s | 10.21 s | 0.75 s | 0.81 s |
| Verifier (process) | 10.42 s | 10.18 s | 10.27 s | - | - |

#### Existing tests

These are the five test files that call the three displays (`git grep` at `76a62a71`). No `R/` function calls them. Each file ran on the fixed tree in its own process, with `NOT_CRAN=true` and a sequential plan (`evidence/F12/existing.txt`, `verify/F12/existing_post.txt`).

| File | Blocks | Expectations | Failed / errors / skipped | Wall (implementer) | Wall (verifier, process) |
|---|---|---|---|---|---|
| `test-forestplot-xlog.R` | 3 | 7 | 0 / 0 / 0 | 3.30 s | 3.39 s |
| `test-mrct-case-study.R` | 7 | 64 | 0 / 0 / 0 | 88.85 s | 90.58 s |
| `test-mr-estimates-table-nullboot.R` | 4 | 33 | 0 / 0 / 0 | 16.50 s | 16.87 s |
| `test-mr_values.R` | 10 | 633 | 0 / 0 / 0 | 5.32 s | 5.47 s |
| `test-no-subgroup-dina-grf-plot.R` | 4 | 22 | 0 / 0 / 0 | 4.87 s | 4.72 s |

- **Existing expectations updated: none.**
- One of these files leaves `tests/testthat/Rplots.pdf` when run through `test_file()`. It is git-ignored, was deleted from the scratch copy before the scope check, and is not in the patch. The new test opens `grDevices::pdf(NULL)` and leaves no file.

#### What a user sees

The example is a continuous fit at the default orientation: `sg.harm = {age <= 47} & {biomarker_hi}`, `adverse_outcome = FALSE`, `max_sg_est = 1.09411`. The output below is `mr_estimates_table(fit)` with `boot_results = NULL`, from `evidence/F12/example.txt`. The verifier's re-run (`verify/F12/example/pre.txt`, `post.txt`) matches it line for line.

Before (start SHA):

```
 Subgroup Mean(C) Mean(T)  Diff         MD (95% CI)        MR_debiased
   Qstnbl   -0.33   -1.42 -1.09   1.09 (0.50, 1.68) 0.77 (-0.13, 1.66)
   Recmnd    0.04   -0.05 -0.10 -0.10 (-0.34, 0.15)                  —
```

After:

```
 Subgroup Mean(C) Mean(T)  Diff        MD (95% CI)        MR_debiased
   Qstnbl   -0.33   -1.42 -1.09  1.09 (0.50, 1.68) 0.77 (-0.13, 1.66)
   Recmnd    0.04   -0.05 -0.10 0.10 (-0.15, 0.34)                  —
```

- Before the fix, the table mixed the two orientations. `Qstnbl` (from `fs$mr_inference$naive`) was +1.09, while `Recmnd` (the table's own estimator) was -0.10.
- After the fix, both rows are on the fit's scale.
- `Mean(C)`, `Mean(T)` and `Diff` stay raw, so for this fit `Diff = -MD` in the same row.
- The same file shows the plot annotation for Questionable changing from `MD = -1.094 (-1.665, -0.523)` to `MD = 1.094 (0.523, 1.665)`.

#### Verification

**Verdict: CONFIRMED.** The verifier reproduced the red run, the green run and the unaffected case. It found 2 problems, both minor. There was no repair or re-verify round for F12, and the commit equals the gated patch: there is no `patches/F12.verified.patch`, and `main/F12/patch_equal.txt` is empty.

| # | Severity | Problem | Resolution |
|---|---|---|---|
| 1 | minor | Other displays still build their effect with the estimator's default `TRUE`. These are `sg_tables()` (`R/summary_utility_functions.R:829-839`, `:982-992`) and the CV summary tables (`R/forestsearch_cross_validation.R:1475-1483`, `:1531-1541` at `76a62a71`). After F12, a continuous fit at the default shows opposite signs: `sg_tables(fit)` reads Questionable `-1.094 (-1.665, -0.523)`, ITT `-0.229` and Recommend `-0.096` (`verify/F12/side_sgtables_post.txt`), against `max_sg_est = 1.094112`. The verifier also named `compute_frontier_cis()`, whose fallback is `%||% TRUE`. | Left as an open item (below). No code changed under F12, since these sites are outside its named scope. |
| 2 | minor | Two tolerance comments in the test header are muddled or incomplete. (a) `:29-30` justify the `1e-4` MR-limit tolerance by "the literal 1.96 rather than qnorm(0.975)", but the expected limits also use 1.96. (b) The forest-plot `1e-10` (`:127-130`, unrounded columns) is not explained in the test, although the header says the tolerances are "the display rounding". | Not changed: the committed test equals the gated patch, and no comment was added before commit. Left as an open item. The verifier notes that neither point affects whether a tolerance separates the trees, because the sign flips move values by 0.19 to 2.2. |

#### Departures and interpretations

- **Forest plot without a fit.** With `fs_results$fs.est = NULL` there is no fit, so no effect is built "for a `forestsearch()` fit". The estimator's default `TRUE` is kept and is now passed explicitly (`R/plot_subgroup_results_forestplot.R:316` at head). The unaffected case `forest_no_fit` is identical to the start SHA. Applying the rule there would flip continuous rows drawn without a fit.
- **All GLM rows of the forest plot.** When a fit is given, the stored value also orients the `subgroup_list` reference and post-hoc rows, not only H and Hc. The reason is that `create_glm_row()` builds every GLM row. The verifier judged this to match the spec's intent.
- **Fallback rule input.** The fallback `outcome_type %in% c("binary", "count")` uses each display's resolved `outcome_type`. This follows the bootstrap pattern at `R/bootstrap_dofuture_main.R:306-309`.
- **"The fit's own estimate".** The fit stores its own estimate only for H (`max_sg_est`, `mr_inference$naive$est`). The test therefore checks every set against a hand value on the fit's orientation, and ties the hand H value to `max_sg_est` within `1e-10`.
- **MR table H-row limits.** These limits are documented as robust (sandwich) limits, so only the H estimate and sign are checked. This narrowing was made before the first red run. No expected value, seed or tolerance changed after the red or green runs.
- **Unaffected case extended.** Beyond the spec's cases, the run adds binary and count fits with `adverse_outcome` unset (resolving to TRUE) and `forest_no_fit`. "Identical" was judged after dropping ggplot2's `S7_class` attribute, which differs between two runs of the same tree (control).
- **No roxygen change.** The spec names a code change only, so no `man/` page changed.

#### paper_simulations

The change does not reach it. `vignettes/articles/paper_simulations.qmd` (2,860 lines at `76a62a71`) names none of the three displays, nor `outcome_type`: `git grep` finds no match. No `R/` function calls the displays, and the article is survival-only.

#### Open items

- **Other displays keep the estimator's default orientation.**
  - `sg_tables()` builds its effect without `adverse_outcome` at `R/summary_utility_functions.R:829-839` (GRF branch) and `:982-992` (main branch). These lines are the same at `76a62a71` and `61c15088`.
  - The CV summary tables do the same at `R/forestsearch_cross_validation.R:1475-1483` and `:1531-1541` (`76a62a71`), which are `:1509-1517` and `:1565-1575` at `61c15088`.
  - For a fit whose resolved value is FALSE, their effect column now has the opposite sign of the fit and of the three F12 displays. Example: `sg_tables()` Questionable MD -1.094 vs `max_sg_est` 1.094.
  - Not fixed here: these sites are outside F12's named sites, and the rule is one fix per commit. The same pattern would apply.
- **`compute_frontier_cis()` fallback.** `R/frontier_cis.R:150` (`76a62a71` and `61c15088`) uses `adverse_outcome = args$adverse_outcome %||% TRUE` instead of `forestsearch()`'s rule. It only matters for an object without a stored value. Not in F12's sites.
- **Test comments (verifier problem 2).** In `tests/testthat/test-display-orientation.R` (`25ecaa0f`), `:27-30` explain the `1e-4` MR-limit tolerance by a z difference that the test does not have, and `:127-130` use `1e-10` without a comment. A comment-only fix was suggested and not applied, because the commit was kept equal to the gated patch.
- **Display docs.** None of the three displays' roxygen says that the effect follows the fit's orientation while the arm means and `Diff` stay raw. For a continuous fit at the default, `Diff = -MD` in the same row. Not fixed: the spec asked for code only.
- **Propensity adjustment.** The displays pass `ps_adjust_method = "none"`. This is explicit in both plots: `R/plot_sg_glm_outcomes.R:270` and `R/plot_subgroup_results_forestplot.R:452` at `61c15088` (`:261`, `:439` at `76a62a71`). `mr_estimates_table()` uses it by default. For an `iptw` or `dr_gcomp` fit the displayed effect is unweighted while the fit's is weighted: the sign is right, but the value differs. Not in scope.
- **Estimator doc.** The `make_effect_estimator()` `@param adverse_outcome` doc (`R/glm_effect_estimators.R:54-60`, both SHAs) says it is "ignored for survival, continuous, and count". However, `.make_lm_estimator()` and `.make_lm_estimator_fast()` negate `y` (`:945-947` and `:1082-1084` at `76a62a71`; `:956-958` and `:1093-1095` at `61c15088`). Not in scope.
- **Binary `effect_measure` default.** `plot_sg_glm_outcomes()` defaults a missing binary `effect_measure` to `"RD"` (`R/plot_sg_glm_outcomes.R:149-152`, both SHAs). `forestsearch()` and `make_effect_estimator()` use `"OR"`. Only an object without `effect_measure` reaches this. Not in scope.
- **Overridden `outcome_type`.** A caller who overrides `outcome_type` to differ from the fit's still gets the fit's stored `adverse_outcome`. A survival fit, for example, stores FALSE. This edge case was left as is.

#### Corrections to the notes

- `notes.md` and the implementer's open items cite `make_effect_estimator()` `@param adverse_outcome` at `R/glm_effect_estimators.R:54-59`. The quoted phrase "ignored for survival, continuous, and count" ends on `:60` (`git show 76a62a71`), so the range is `:54-60`.
- The verifier cites `compute_frontier_cis()` at `R/frontier_cis.R:149`. The `%||% TRUE` line is `:150` at both `76a62a71` and `61c15088`; `:149` is `effect_measure = effect_measure,`.
- The verifier's CV-table line numbers (`:1475-1483`, `:1531-1541`) hold at the start SHA only. At head `61c15088` they are `:1509-1517` and `:1565-1575`.

### F13 — A-055: risk-difference standard error, unadjusted case (change of method)

**Status:** PASS. Commit `89ad9c19` ("fix: risk-difference standard error in the unadjusted case (A-055)"), parent `25ecaa0f`. Start SHA `76a62a71`.

| file (`git show --stat 89ad9c19`) | change |
|---|---|
| `R/glm_effect_estimators.R` | 13 +- (numstat: 12 insertions, 1 deletion) |
| `tests/testthat/test-rd-standard-error.R` | new, 82 lines |
| total | 2 files, 94 insertions, 1 deletion |

- No roxygen or `man/` change: `.estimate_rd` is `@noRd` (start `:456`).
- The committed diff equals the gated patch `patches/F13.patch` apart from `index`/`@@` lines (`main/F13/patch_equal.txt` is empty; I re-ran the `diff` and it is empty too). There is no `patches/F13.verified.patch`.
- `R/glm_effect_estimators.R` is byte-identical at `89ad9c19`, `61c15088` and `88302212` (same md5). It is unchanged from `76a62a71` to `25ecaa0f`. The test file is unchanged from `89ad9c19` to `88302212`.

#### Site, before and after

`R/glm_effect_estimators.R`, `.estimate_rd()` (start `:457`), tier 2 (logistic g-computation, start `:514-548`, where the spec pointer said). Lines 1-535 are the same at start and head.

Before, start `76a62a71`, `:534-536`:

```r
    # Convert logistic SE to RD scale using average marginal effect
    avg_deriv <- mean(pred1 * (1 - pred1))
    se_rd <- abs(avg_deriv) * se_treat
```

After, head `61c15088`, `:534-547` (new branch `:536-547`; the new SE is `:544`, the old expression is kept at `:546`):

```r
    # Convert logistic SE to RD scale using average marginal effect
    avg_deriv <- mean(pred1 * (1 - pred1))
    se_rd <- if (n_adj == 0L) {
      # Unadjusted (A-055): the two-arm logistic model is saturated, so the
      # arm means of pred0 and pred1 are the observed proportions p0 and p1,
      # and the exact delta-method SE of rd_est is the two-sample binomial
      # SE.  The scaled-coefficient form below drops the control arm's
      # variance and degenerates when an arm is at 0% or 100%.  No weights
      # reach this tier: the IPTW path (the only weighted one) returns
      # above, and neither glm() call here takes a weights argument.
      sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
    } else {
      abs(avg_deriv) * se_treat
    }
```

- `p0`, `p1`: start `:461-462`, `mean(y0, na.rm = TRUE)` and `mean(y1, na.rm = TRUE)`. `n0`, `n1`: arm sizes from the caller, start `:394-395`. `n_adj`: start `:480`.
- Unchanged: `converged = TRUE` (start `:541`, head `:552`), `measure = "RD"` (`:544` / `:555`), `method_used = if (n_adj > 0) "logistic_gcomp_adj" else "logistic_margins"` (`:545` / `:556`), the `error = function(e) NULL)` handler (`:548` / `:559`).
- `se_treat` (`:533`) and `avg_deriv` (`:535`) are still computed in both branches. So the slices that fall to tier 3 are the same as before. Nothing is routed to the raw-proportion tier.

**How weights reach `.estimate_rd`** (lines at the start SHA; all lie before `:536`, so they are the same at head unless marked):

| path | what reaches `.estimate_rd` | weighted fit? | reaches the new branch? |
|---|---|---|---|
| Signature `:457-458` | `data_slice, fmla, treat.name, outcome.name, n0, n1, adjust_covariates, ps_adjust_method`. No weights argument. | — | — |
| Only caller: `.make_glm_binary_estimator()` (`:366`), RD arm `:414-419` | `fmla` from `.build_adjusted_formula()` (`:1236-1243`). Covariates are passed only when `ps_adjust_method == "none"`, else `NULL` (`:416`, `:418`). | — | — |
| `ps_adjust_method = "iptw"` with an `sw` column (column made by `estimate_propensity_scores()`, `:1444-1451`, or by `forestsearch()` for a user `ps_hat`, `R/forestsearch_main.R:2720-2725`) | weighted means at `:465-478`; returns early with `method_used = "iptw_rd"` | yes | no (returns before tier 1) |
| Tier 1 `glm()` `:485-490`; tier 2 `glm()` `:516-520` | `stats::glm(fmla, data = data_slice, family = ...[, start = ...])`; no `weights =` | no | — |
| `"none"`, no covariates | `fmla = y ~ treat`, `n_adj = 0` | no | yes |
| `"dr_gcomp"`, no covariates | same as `"none"` (RD does not use `ips_covar`) | no | yes |
| `"iptw"` without an `sw` column | same as `"none"` | no | yes |
| **Verifier's addition:** `adjust_covariates` with `"dr_gcomp"`, or with `"iptw"` and no `sw` column | covariates dropped at `:416`/`:418`, so `n_adj = 0` and the new SE applies. Only direct `make_effect_estimator()` calls: `forestsearch()`'s PS rebuild (`R/forestsearch_main.R:2748-2756`) passes no `adjust_covariates`. | no | yes |
| `"none"` with covariates | `n_adj > 0` | no | no: keeps `abs(avg_deriv) * se_treat` (head `:546`) |

- `make_effect_estimator()` (start `:100-112`) has no weights argument. Its 12 call sites at the start SHA (`git grep`): `bootstrap_dofuture_main.R:310`, `forestsearch_cross_validation.R:1476, 1532`, `forestsearch_main.R:2263, 2748`, `frontier_cis.R:143`, `mr_estimates_table.R:204`, `plot_sg_glm_outcomes.R:254`, `plot_subgroup_results_forestplot.R:432`, `run_subgroup_sims.R:321`, `summary_utility_functions.R:830, 983`. Weights exist only as the `sw` column.
- Conclusion: inside tier 2, "no adjustment covariates and no weights" is exactly `n_adj == 0L`. The BLOCKED condition was not met.
- Run evidence (`evidence/F13/example.txt`, section 3; control 30%, treated 100%): `"iptw"` with `sw` gives `se 0.072457 ... method iptw_rd` on both trees. The other five paths go from `se 5.4089e-06` to `se 0.072457`, `method logistic_margins`.
- Verifier probe (`verify/F13/probe_drgcomp_adj.txt`), `adjust_covariates = "x"`, same slice: `dr_gcomp` and `iptw` (no `sw`) go from `se 5.40885e-06` to `se 0.0724569`, `method logistic_margins`. `"none"` with `"x"` stays `se 5.34438e-06`, `method logistic_gcomp_adj`.

#### Red and green

Test file `tests/testthat/test-rd-standard-error.R`, filter `rd-standard-error`.

| run | where | result | wall |
|---|---|---|---|
| implementer red | scratch copy of the start SHA (`wt/F13`), test file only | 14 expectations, 4 failed, 0 errors (`evidence/F13/red.txt`, `red_counts.txt`) | 3.10 s |
| implementer green | same scratch copy, patch applied | 14 expectations, 0 failed, 0 errors, 0 warnings (`green.txt`, `green_counts.txt`) | 2.89 s |
| verifier red | start export, test file only (`verify/F13/pre`) | 14 expectations, 4 failed, 0 errors, 0 skipped (`verify/F13/pre_rdse.txt`) | 2.77 s in R, 3.01 s process |
| verifier green | start export with the patch (`verify/F13/post`) | 14 expectations, 0 failed, 0 errors, 0 skipped (`post_rdse.txt`) | 2.76 s in R, 3.01 s process |
| main-session red | repository at `25ecaa0f`, test file only | `rd-standard-error: 12.....34.....`: 14 expectations, 4 failures (`main/F13/red.txt`) | 2.8 s |
| main-session green | repository, code applied | `rd-standard-error: ..............`: 14 expectations, 0 failures (`main/F13/green.txt`) | 2.7 s |

- The four failures are the four SE checks, the same in all three red runs:
  - `:50` "Expected `res$se` to equal `sqrt(0.3 * 0.7/40 + 1 * 0/40)`." `actual: 0.000`, `expected: 0.072`.
  - `:51` "Expected `round(res$se, 5)` to equal 0.07246." Same actual.
  - `:72` "Expected `res$se` to equal `sqrt(0.8 * 0.2/40 + 0 * 1/40)`." `actual: 448.5`, `expected: 0.1`.
  - `:73` "Expected `round(res$se, 5)` to equal 0.06325." Same actual.
- Unrounded old values (`example.txt`, section 1): `se 5.4089e-06` and `se 448.55`.
- The estimate, `converged`, `method_used`, `measure` and `n0`/`n1` checks pass in red as well. They guard the spec's "do not change" clause.

#### Expected values

Test file at `89ad9c19`. The slices hold 40 patients per arm (`:21-27`).

| expected value (line) | how it was derived |
|---|---|
| `sqrt(0.3 * 0.7 / 40 + 1.0 * 0.0 / 40)` = 0.0724569 (`:50`) | By hand from the spec's formula, p0 = 12/40, p1 = 40/40. Comment `:47-49` (0.21/40 = 0.00525; 0.0724569^2 = 0.0052500). |
| `0.07246` (`:51`) | Spec value sqrt(0.21/40), compared with `round(res$se, 5)`. |
| `1.0 - 0.3` (`:54`) | Difference of the arm proportions (spec). Comment `:53`. |
| `sqrt(0.8 * 0.2 / 40 + 0.0 * 1.0 / 40)` = 0.0632456 (`:72`) | By hand, p0 = 0/40, p1 = 32/40. Comment `:69-71` (0.16/40 = 0.004; 0.0632456^2 = 0.0040000). |
| `0.06325` (`:73`) | Spec value sqrt(0.16/40), compared with `round(res$se, 5)`. |
| `0.8 - 0.0` (`:76`) | Difference of the arm proportions (spec). Comment `:75`. |
| estimate tolerance `.rdse_est_tol <- 1e-6` (`:41`, used at `:54`, `:76`) | Set before any run. glm() stops at a relative deviance change below 1e-8 under quasi-separation. Comment `:34-40`. The verifier measured deviations of about 1e-9 to 1e-8. The SE checks use the default tolerance. |
| `converged` TRUE (`:59`, `:78`) | Spec: unchanged. Tier 2 sets it at start `:541`. Comment `:56-58`. |
| `method_used` `"logistic_margins"` (`:60`, `:79`) | Unadjusted tier-2 label, start `:545`. Tier 3 would give `"raw_means"` with `converged = FALSE`. Comment `:56-58`, header `:14-17`. |
| `measure` `"RD"` (`:61`, `:80`) | Start `:544`. |
| `c(res$n0, res$n1)` = `c(40L, 40L)` (`:62`, `:81`) | Slice size (`:21`). Integer because `n0`/`n1` come from `sum()` at start `:394-395`. |

- **Flag (verifier, minor):** in the second block (`:78-81`), the `converged`, `method_used`, `measure` and `n0`/`n1` expectations have no comment of their own. The first block's comment (`:56-58`) and the file header cover them. The suggested one-line comment was not added before the commit.
- Also, `measure` `"RD"` has no comment naming its source in either block.

#### Unaffected case

- **Spec:** the odds-ratio and risk-ratio measures, and a risk difference with adjustment covariates, give identical output.
- **Run** (`evidence/F13/unaffected.R`, `unaffected.txt`). Each tree runs in its own `Rscript --vanilla` with `pkgload::load_all()`, a sequential plan and `set.seed` before each call. Trees: start export vs scratch `wt/F13`.
  - Part A: `make_effect_estimator(outcome_type = "binary")`, 14 configurations x 8 slices = 112 results.
    - Configurations: OR unadjusted; OR with `c("x","z")`; RR unadjusted; RR with `c("x","z")`; RR with `robust_se = FALSE`; RD with `"x"`; RD with `c("x","z")`. Each with `adverse_outcome` TRUE and FALSE.
    - Slices: 3 interior (n = 80, 200, 24) and 5 with a degenerate arm (control/treated 30/100, 0/80, 50/0, 100/50, 0/100 %).
    - The 32 adjusted-RD results used `identity_binomial_adj` (12) and `logistic_gcomp_adj` (20). So the adjusted tier 2, which keeps the old expression, was exercised.
  - Part B: `forestsearch()` on a synthetic binary outcome (N = 300). Once with `effect_measure` unset (resolves to OR), and once with `effect_measure = "RD", adjust_covariates = "age"`.
  - Comparison: `identical()`, after dropping timing leaves and converting data.tables to data.frames.
- **Result:** `Part A identical(): TRUE`, both Part B runs `TRUE`, `Differing leaves: none`, `OVERALL identical: TRUE`.
- **Wall:** start 4.75 s, fixed 4.85 s, compare 0.25 s (process walls).
- **Verifier:** re-ran the script on the start export and on `verify/F13/post`. `OVERALL identical: TRUE` (`verify/F13/unaffected_rerun.txt`); process walls 4.82 s and 4.87 s.
  - Own probe (`verify/F13/my_unaffected_compare.txt`): OR/RR/RD x no/`x`/`x,z` covariates x `ps` none/iptw-with-`sw` x `adverse_outcome` x 9 slices = 324 results.
  - 312 were identical. All 12 that differ are RD, unadjusted, `ps` none, `logistic_margins`, and only `se` changed.
  - The 288 spec-named results were all identical (`Spec unaffected subset: 288 results, identical: TRUE`). No wall recorded.
- **Scope note:** "RD with adjustment covariates" is identical for `ps_adjust_method = "none"` (and for `"iptw"` with an `sw` column). With `"dr_gcomp"`, or `"iptw"` without `sw`, the covariates are dropped and the SE changes (see the weights table).

#### Existing tests

Files chosen by grepping `tests/testthat` for `"RD"`, `.estimate_rd`, `make_effect_estimator` and the RD consumers. The implementer compared the fixed tree with `R/` at the start SHA (`git stash`; `evidence/F13/existing.txt`, `existing_start.txt`). The verifier compared post with pre (`verify/F13/existing_post.txt`, `existing_pre.txt`).

| file | blocks / expectations | failed / errors | skipped | start vs fixed | wall, fixed (implementer, process) |
|---|---|---|---|---|---|
| estimability-boundary | 13 / 60 | 0 / 0 | 0 | identical | 3.1 s |
| threshold-pair-directive-a | 17 / 86 | 0 / 0 | 1 | identical | 6.9 s |
| binary-default-or | 9 / 60 | 0 / 0 | 1 | identical | 9.9 s |
| binary-default-or-entry-points | 7 / 25 | 0 / 0 | 0 | identical | 2.9 s |
| lean-fits | 4 / 362 | 0 / 0 | 0 | identical | 6.6 s |
| fs_dgm_scale | 21 / 45 | 0 / 0 | 0 | identical | 3.6 s |
| fs-dgm-feasibility | 8 / 40 | 0 / 0 | 0 | identical | 3.3 s |
| dina-gaussian-floor | 5 / 31 | 0 / 0 | 0 | identical | 6.1 s |
| iptw-continuous | 2 / 7 | 0 / 0 | 0 | identical | 3.7 s |
| directive-c | 31 / 155 | 0 / 0 | 1 | identical | 8.7 s |
| no-subgroup-dina-grf-plot | 4 / 22 | 0 / 0 | 0 | identical | 5.0 s |
| cross-outcome-parity | 9 / 49 | 0 / 0 | 2 | identical (12 warnings on both) | 45.9 s |
| glm-pipeline-integration | 4 / 10 | 0 / 0 | 0 | identical | 5.8 s |
| mr-inference | not completed | none printed | — | capped on both trees after the same 133 passing expectations (exit 124) | 120.1 s |

- **Existing expectations updated: none.**
- An instrumented copy (`verify/F13/instr_hits.txt`) found that only estimability-boundary reaches the new branch among existing tests (2 hits).
  - Its RD block (`tests/testthat/test-estimability-boundary.R:75-81`) asserts only that the estimate is not NA and `reason` is NULL.
  - Its SE on those frames moves from `1.5299e-05` to `0.21909` = sqrt(0.4*0.6/5) (`example.txt`, section 5).
- mr-inference uses binary OR. "RD" appears only in the constant lookup `.fs_mr_confirm_null("RD")`. It reaches the new branch 0 times (capped run).

#### What a user sees

`make_effect_estimator("binary", effect_measure = "RD")` on 40 per arm, control 30%, treated 100%. The table is the ITT row of `SG_tab_estimates_glm()` (`evidence/F13/example.txt`, sections 1 and 2; the verifier's re-run `verify/F13/example_rerun.txt` matches).

Before (start SHA):

```
control  30%, treated 100%: estimate 0.7000  se 5.4089e-06  95% CI (0.7000, 0.7000)  converged TRUE  method logistic_margins
             RD (95% CI)
 0.7000 (0.7000, 0.7000)
```

After (F13):

```
control  30%, treated 100%: estimate 0.7000  se 0.072457  95% CI (0.5580, 0.8420)  converged TRUE  method logistic_margins
             RD (95% CI)
 0.7000 (0.5580, 0.8420)
```

#### Verification

- **Verdict:** CONFIRMED. Red, green and the unaffected case were reproduced, and the expected values were re-derived by hand.
- No repair round: the `wf_results.json` entry holds only `impl` and `verify`.
- Problems, all minor:

| # | problem | resolution |
|---|---|---|
| 1 | The notes' weights trace left out one route: `adjust_covariates` with `dr_gcomp`, or with `iptw` and no `sw`. There the covariates are dropped at `:416`/`:418`, so the new SE applies. | Row added to the weights table in this report. No code change. |
| 2 | The comment kept above the change (start/head `:531-532`, "this is conservative for subgroup estimation") contradicts the new comment (head `:540-541`), which says that form "drops the control arm's variance". | Left as an open item. Only the adjusted case still uses that form. |
| 3 | The second test block (`:78-81`) has no comment of its own for `converged`, `method_used`, `measure`, `n0`/`n1`. | Not added. The main session added no comment for F13: there is no `F13.verified.patch`, and the commit equals `F13.patch`. Noted under "Expected values" and in the open items. |

- Setup note from the verifier: the harness command `git init && git apply --index` fails on a fresh repo, because the archive's files are not yet in the index. The verifier committed the archive first, and the patch then applied cleanly. This is a harness flaw, not a patch flaw.

#### Departures and interpretations

- **"No adjustment covariates and no weights"** is read as `n_adj == 0L` inside tier 2. The reason: no weighted fit reaches tier 2 (weights table).
- **Covariates dropped before the call.** With `dr_gcomp`, or with `iptw` and no `sw`, the closure drops `adjust_covariates` before `.estimate_rd`. The fitted model then has no covariates and no weights, so the new SE applies. This is read as within the spec. As a result, the unaffected clause holds for RD with covariates under `ps_adjust_method = "none"` only (plus `iptw` with `sw`, which returns before tier 2). `forestsearch()` never builds the other combination.
- **"Do not route these cases to the raw-proportion tier."** `se_treat` and `avg_deriv` are still computed in the unadjusted branch, so failures route exactly as before. The tests assert `method_used = "logistic_margins"` and `converged = TRUE`.
- **Extra test assertions.** The tests go beyond the spec's two values: `round(se, 5)`, `method_used`, `measure`, `n0`/`n1`. All come from the spec or the start-SHA source.
- **Estimate tolerance.** The 1e-6 tolerance was set before any run, from glm's 1e-8 stopping rule.
- **mr-inference** did not finish within the 2-minute cap on either tree. It was judged out of reach of `.estimate_rd`: binary OR only, 0 hits in the instrumented run.

#### paper_simulations

- The implementer reports that F13 does not reach the changed lines. The article's evaluated chunks call `run_simulation_analysis()` with no `outcome_type` or `effect_measure` (survival), so `is_glm` is FALSE. `.estimate_rd` is reached only with `outcome_type = "binary"` and `effect_measure = "RD"`, and neither appears in an evaluated chunk.

#### Open items

- **Adjusted tier-2 SE (OPEN per spec):** with `n_adj > 0`, `se_rd` is still `abs(avg_deriv) * se_treat` (head `R/glm_effect_estimators.R:546`, start `:536`). It drops the control arm's variance and degenerates at 0%/100% arms. Not fixed because the spec says to record it, not fix it.
- **Both arms degenerate (OPEN per spec):** with control 0% and treated 100%, the new SE (head `:544`) is exactly 0 (was `2.3099e-07`), and the 95% CI is `(1.0000, 1.0000)` (`example.txt`, section 1). Not fixed because the spec says to record it.
- **Stale comment:** "this is conservative for subgroup estimation" (head `:531-532`, same at start) does not hold for the adjusted form that remains, and it contradicts the new comment at `:540-541`. Comment only. It belongs with the adjusted-case fix, which is out of scope here.
- **Test comment:** `tests/testthat/test-rd-standard-error.R:78-81` (at `89ad9c19`) has no comment of its own for `converged`, `method_used`, `measure`, `n0`/`n1`. The verifier rated this minor, and nothing was added before the commit.
- **Empty arm:** with `n0 = 0` (or `n1 = 0`) the slice reaches tier 2 through a rank-deficient fit (estimate 0, `converged` TRUE, as at the start SHA). The SE at head `:544` changes from `NA` to `NaN` (`is.na()` TRUE for both; `example.txt`, section 4). Pre-existing routing, not in the spec.
- **Missing outcomes:** `p0`, `p1` use `na.rm = TRUE` (start `:461-462`), while `n0`, `n1` count every row of each arm (start `:394-395`). With missing outcomes, the new SE uses arm sizes that include those rows. Tier 3 does the same (start `:554-556`). Found by reading the source; not exercised.
- **IRD shortcut:** `se = abs(avg_deriv) * se_treat` in `.make_poisson_rate_estimator()`'s `dr_gcomp` IRD branch (start `:845`, head `:856`). Not part of A-055.
- **Tier-1 warning handler:** it calls `invokeRestart("muffleWarning")` from an exiting `tryCatch` handler (start and head `:505-510`). Carried over from the A-055 audit; not touched, not exercised.
- **mr-inference:** `tests/testthat/test-mr-inference.R` does not finish within the 120 s cap on either tree (exit 124, 133 expectations passed, no failure). It does not reach `.estimate_rd`.
- **Non-ASCII:** "León" in roxygen at start `:65` and `:71` was already there. This fix adds no non-ASCII characters.

#### Corrections to the notes

- Tier-1 `glm()` call: `git show 76a62a71` puts it at `:485-490`. The notes say `:487-492`.
- IRD shortcut: it is at start `:845` (head `:856`). The notes and the JSON open item say "start :856".
- IPTW early return: the block is `:465-478`; `:479` is blank. The notes say `:465-479`.
- User `ps_hat` in `forestsearch()`: `sw` is assigned at `R/forestsearch_main.R:2723-2725`. The notes say `:2720-2723`.

### F14 — DINA floor: fit, record and replay agree in the forced-gaussian corner

**Status:** PASS. Commit `685e0c5d` ("fix: DINA floor: fit, record and replay agree in the forced-gaussian corner"), parent `89ad9c19` (F13).

- `git show --stat 685e0c5d`:
  - `R/forestsearch_main.R` | 16 (+7 −9)
  - `tests/testthat/test-dina-gaussian-floor.R` | 84 (+77 −7)
  - 2 files changed, 84 insertions(+), 16 deletions(-).
- Scope gate: only the fix's code file and its test file. No roxygen, `man/`, `NEWS.md` or `DESCRIPTION` change.
- `main/F14/patch_equal.txt` is empty: the committed diff equals the gated `patches/F14.patch` (index and hunk-header lines aside). There is no `F14.verified.patch`.
- Setting (spec): binary RD or count IRD, `dina_args = list(family = "gaussian")`, and a supplied ratio-scale threshold that the package warns about and remaps (RD to 0.05, IRD to 0.01).
- Defect at the start SHA: the fit used the raw supplied value (for example 1.5) as the DINA floor. The object recorded the remapped value, and a bootstrap or CV replay of `args_call_all` used that remapped value.

#### Site, before and after

Both sites are in `R/forestsearch_main.R`. The spec pointers at `28e71ffa` were "near 2784–2785 and 3251–3253". I located both sites by content. Later commits (`a4cce90a` F3, `61c15088` F11) do not move these lines. In `R/forestsearch_main.R`, the only later change is F11's hunk at head `:3297-3301`. So the after-lines are the same at `685e0c5d` and at head `61c15088`.

**Site 1: `subgroup_method = "dina"`.** This is the last argument of the `.forestsearch_dina_select(` call that starts at `:2806` on both SHAs.

Start `76a62a71`, `R/forestsearch_main.R:2820-2824`:
```r
      admission = admission_resolved,
      # Gaussian floor: the resolved threshold when neither threshold was
      # supplied (A-001); a supplied threshold is used as given.
      effect.threshold = if (user_set_threshold) hr.threshold
                         else effect.threshold)
```
Head `61c15088`, `R/forestsearch_main.R:2820-2823`:
```r
      admission = admission_resolved,
      # Gaussian floor: the resolved threshold (A-001), supplied or not, so
      # the fit, the recorded value and a bootstrap / CV replay agree.
      effect.threshold = effect.threshold)
```

**Site 2: `use_dina = TRUE` selected-cut screening (`m_diff_sel`).**

Start `76a62a71`, `R/forestsearch_main.R:3283-3291`:
```r
        # continues with no DINA cuts.
        # With neither threshold supplied, the gaussian floor is the resolved
        # effect.threshold (SECTION 2B-ii), the value the returned object
        # records and a replay re-uses (A-001); a supplied threshold is used
        # as given.
        .dina_assert_ratio_estimand(da$fit$family, effect_measure)
        m_diff_sel <- if (identical(da$fit$family, "gaussian")) {
                        if (user_set_threshold) hr.threshold else effect.threshold
                      } else log(hr.threshold)
```
Head `61c15088`, `R/forestsearch_main.R:3282-3289`:
```r
        # continues with no DINA cuts.
        # The gaussian floor is the resolved effect.threshold (SECTION
        # 2B-ii), supplied or not: the value the returned object records and
        # a replay re-uses (A-001).
        .dina_assert_ratio_estimand(da$fit$family, effect_measure)
        m_diff_sel <- if (identical(da$fit$family, "gaussian")) {
                        effect.threshold
                      } else log(hr.threshold)
```

**Why only the corner changes.** I read this at `76a62a71`:

- `effect_threshold <- hr.threshold` (`:2318`).
- The only thing that changes a supplied value is the RD/IRD remap: `if (user_set_threshold && effect_threshold > 1.0)` (`:2334`), then `switch(effect_measure, RD = 0.05, IRD = 0.01)` (`:2348-2349`).
- SECTION 2B-ii sets `effect.threshold <- if (.fs_identity_scale) effect_threshold else hr.threshold` (`:2573-2575`).
- So with a supplied threshold, `effect.threshold` and `hr.threshold` differ only after that remap. With no threshold supplied, both sites already used `effect.threshold`.

**Test-file changes** (`tests/testthat/test-dina-gaussian-floor.R` at `685e0c5d`):

- `.gf_replay()` gains a defaulted `df`:
  - `:63` `.gf_replay <- function(fs, df = .gf_df()) {`
  - `:65` `a$df.analysis <- df`
  - Its comment is extended at `:61-62`.
- The pinning block is updated (`:124-153`, was `:123-147`).
- A new corner test is added (`:155-217`).

#### Red and green

All runs use `devtools::test(filter = "dina-gaussian-floor")` with `NOT_CRAN=true`, under `timeout 120`. The two main-session runs went through `..._apply_fix.sh` and ran under `timeout 300`. The summary line has 55 positions. Failures are counted from the `Failure (` blocks.

| run | where | result | wall |
|---|---|---|---|
| implementer red | scratch copy `wt/F14` at start SHA `76a62a71`, new test file only (`evidence/F14/red.txt`, `red_full.txt`) | 55 expectations: 45 pass, **10 failures**, 0 errors. Summary `.........................1.2....3.4...5.6...7.8...9.a..`. `red_full.txt` (`max_reports = 100`) confirms 10, none cut off. 2 failures at `:150:5` (`dina`, `use_dina`): "Expected `unique(fit$floors)` to be identical to 0.05." `actual`: 1.5. 4 at `:209:7` (RD/IRD × dina/use_dina): "Expected `fit$floors` to be identical to `fs$args_call_all$effect.threshold`." `actual`: 1.5. 4 at `:214:7`: "Expected `rep$floors` to be identical to `fit$floors`." `expected`: 1.5. | 11.36 s (`red.time`), 11.29 s (`red_full.time`) |
| implementer green | `wt/F14` = start SHA + `F14.patch` (`evidence/F14/green.txt`) | 55 expectations, 0 failures, 0 errors, 0 skips | 14.55 s (`green.time`) |
| verifier red | `verify/F14/pre` (`git archive 76a62a71` + new test file) (`verify/F14/red.txt`) | 55 expectations: 45 pass, **10 failures**, 0 errors. The file is byte-identical to the implementer's `red_full.txt`. | 11.73 s (`verify/F14/red.time`) |
| verifier green | `verify/F14/post` (`76a62a71` + `F14.patch`) (`verify/F14/green.txt`) | `nb 55, failed 0, skipped 0, error 0` | 14.97 s (`verify/F14/green.time`) |
| verifier: old test file on fixed code | `verify/F14/post_oldtest` (`verify/F14/oldtest_on_post.txt`) | 2 failures, both at `:144:5` (`dina`, `use_dina`): "Expected `unique(fit$floors)` to be identical to 1.5." `actual`: 0.1 (waldo's display of 0.05). This is the approved old-behaviour pin. | 8.06 s (verifier JSON; no time file) |
| main-session red | repository at `89ad9c19` (parent of `685e0c5d`) with only the test file applied (`main/F14/red.txt`) | 55 expectations: 45 pass, **10 failures**, 0 errors. Same summary line and failures as the implementer's red (the files are byte-identical). | 10.8 s (`main/F14/wall.txt`) |
| main-session green | repository at `685e0c5d` (`main/F14/green.txt`) | 55 expectations, 0 failures, 0 errors | 14.0 s (`main/F14/wall.txt`) |

- waldo prints 0.05 as `0.1` and 0.01 as `0.0` at its display precision. The `:209` failures for IRD show `expected`: 0.0.
- `red.txt` and `main/F14/red.txt` end with "Maximum number of 10 failures reached, some test results may be missing." The summary line still marks all 55 positions with exactly 10 failures, and `red_full.txt` confirms that no failure was cut off.

#### Expected values

Line numbers are for `tests/testthat/test-dina-gaussian-floor.R` at `685e0c5d`. The file is unchanged through `88302212`.

| line(s) | expected value | how it was derived | F14 change |
|---|---|---|---|
| `:150` | `unique(fit$floors)` = `0.05` (binary RD, gaussian forced, supplied 1.5; both routes) | Comment `:147-149`: the RD row of `?forestsearch` "Resolved defaults, per estimand" (`man/forestsearch.Rd:1207`, c1 = 0.05) and the remap warning text "Remapping to default: 0.05". The floor rule is in comment `:127-129`, quoting `?forestsearch` `dina_args` (`man/forestsearch.Rd:233-234`: "identity for \code{gaussian} -- from the resolved \code{effect.threshold}"). | updated: was `1.5` at `:144` |
| `:151` | `args_call_all$effect.threshold` = `0.05` | The same comment `:147-149`. The value is the same as at the start SHA (`:145`) and passes on both trees. | moved only (`:145` → `:151`) |
| `:173` → used at `:206`, `:207` | RD: `args_call_all$effect.threshold` and `threshold_config$screening` = `0.05` | Comment `:159-161`: the table's RD c1 = 0.05. RD is on the identity scale, so `screening` holds the value itself, not its log. | new |
| `:179` → used at `:206`, `:207` | IRD: the same two fields = `0.01` | Comment `:159-161`: the table's IRD c1 = 0.01 (`man/forestsearch.Rd:1208`). **Verifier flag (minor):** this table covers the case where no threshold is supplied. A supplied out-of-range value is remapped only per the warning text, and for IRD that text says 0.05 (see Open items). The value is right, but the cited derivation does not literally cover the corner. | new |
| `:204`, `:213` | `expect_length(..., 1L)` | Comment `:203`: one `dina_subgroup()` call per fit, because each route has one derivation site (header `:5-7`). | new |
| `:209` | fit floor == recorded `args_call_all$effect.threshold` (relational) | Comments `:162-165` and `:208`: the `?forestsearch` `dina_args` floor rule. This is the spec's own test. | new |
| `:214` | replay floor == fit floor (relational) | Comments `:164-165` and `:211`. The replay follows the bootstrap's `args_call_all` pattern (comment `:59-62`, `R/bootstrap_analysis_dofuture.R`). | new |
| `:85`, `:92`, `:100`, `:106` | `expect_length(..., 1L)` | Header `:5-7`, `:14-15`. | unchanged |
| `:86-88`, `:101-102` | `0` (continuous MD, both thresholds unset) | Header comment `:8-12` ("MD resolves the screening threshold to 0 and records 0"). This matches `man/forestsearch.Rd:1209` (MD c1 = 0.0). | unchanged |
| `:89`, `:93`, `:103`, `:107` | relational (floor == recorded; replay == fit) | Header `:8-12`. | unchanged |
| `:116-118` | `0.4` (supplied, either spelling, both routes) | The supplied value itself. The block title is at `:110`. | unchanged |
| `:234` | `0.3` | The supplied `hr.threshold`. The helper's default is `effect.threshold = hr.threshold` (comment `:220-221`). | unchanged |
| `:235` | `0` | The supplied `effect.threshold = 0`. | unchanged |

- The 1.5 in the corner calls is an input, not an expected value.
- The verifier re-derived every expected value. It flagged no value as copied from the fixed code's output and no expected value as lacking a comment. Its only flag is the one in the IRD row.

#### Unaffected case

**What the spec names.** For every other threshold configuration the existing test file covers, and for a continuous outcome with a supplied threshold, the fit's floor is unchanged from the start SHA. If any row other than the corner differs, the fix is FAILED.

**What was run** (`evidence/F14/unaffected.R`):

- Trees: the start export `/tmp/cran_prep_code_fixes_2026-11-12/start` (`76a62a71`) and `wt/F14` (start + `F14.patch`). Each tree ran in its own `Rscript --vanilla` process with `pkgload::load_all(tree)`, a sequential plan and `set.seed(20261112)` before every call.
- Data:
  - continuous: `.make_continuous_data(N = 200L, MD_harm = 1.5, seed = 42L)`
  - binary: `set.seed(7); N <- 300L`, the same data as the test file
  - count: `.make_count_data(N = 300L, IRR_harm = 2.0, seed = 7L)`
  - survival: `.make_survival_data(N = 300L, HR_harm = 2.0, seed = 42L)`
- The floor was captured by wrapping `dina_subgroup()`, the only consumer of `m_diff` at both sites.
- Each `forestsearch()` row also replays the parent's `args_call_all`.
- Whole fit and replay objects were compared with `identical()`, after dropping timing leaves and turning data.tables into data.frames. The start-tree objects contain no function, environment or external-pointer leaf.
- Row groups:
  - R01-R10: the configurations the test file covers at the start SHA.
  - R11-R12: the spec's "continuous outcome with a supplied threshold" (MD with 1.5, which MD never remaps).
  - R13-R26: further checks.

**Result** (`evidence/F14/unaffected.txt`; the notes' table agrees row for row):

| row | configuration | fit floor, start SHA | fit floor, fixed | identical? | replay floor start / fixed | whole object identical (timing dropped)? |
|---|---|---|---|---|---|---|
| R01 | continuous MD, both unset, dina | 0 | 0 | yes | 0 / 0 | yes |
| R02 | continuous MD, both unset, use_dina | 0 | 0 | yes | 0 / 0 | yes |
| R03 | continuous MD, effect.threshold = 0.4, dina | 0.4 | 0.4 | yes | 0.4 / 0.4 | yes |
| R04 | continuous MD, hr.threshold = 0.4, dina | 0.4 | 0.4 | yes | 0.4 / 0.4 | yes |
| R05 | continuous MD, effect.threshold = 0.4, use_dina | 0.4 | 0.4 | yes | 0.4 / 0.4 | yes |
| R06 | continuous MD, hr.threshold = 0.4, use_dina | 0.4 | 0.4 | yes | 0.4 / 0.4 | yes |
| R07 | **CORNER** binary RD, gaussian forced, effect.threshold = 1.5, dina | 1.5 | 0.05 | **no** (intended) | 0.05 / 0.05 | **no** |
| R08 | **CORNER** binary RD, gaussian forced, effect.threshold = 1.5, use_dina | 1.5 | 0.05 | **no** (intended) | 0.05 / 0.05 | **no** |
| R09 | `.forestsearch_dina_select()`, hr.threshold = 0.3, no effect.threshold | 0.3 | 0.3 | yes | NA / NA | yes |
| R10 | `.forestsearch_dina_select()`, hr.threshold = 0.3, effect.threshold = 0 | 0 | 0 | yes | NA / NA | yes |
| R11 | continuous MD, effect.threshold = 1.5 (never remapped for MD), dina | 1.5 | 1.5 | yes | 1.5 / 1.5 | yes |
| R12 | continuous MD, effect.threshold = 1.5 (never remapped for MD), use_dina | 1.5 | 1.5 | yes | 1.5 / 1.5 | yes |
| R13 | **CORNER** count IRD, gaussian forced, effect.threshold = 1.5, dina | 1.5 | 0.01 | **no** (intended) | 0.01 / 0.01 | **no** |
| R14 | **CORNER** count IRD, gaussian forced, effect.threshold = 1.5, use_dina | 1.5 | 0.01 | **no** (intended) | 0.01 / 0.01 | **no** |
| R15 | **CORNER** binary RD, gaussian forced, hr.threshold = 1.5 (legacy spelling), dina | 1.5 | 0.05 | **no** (intended) | 0.05 / 0.05 | **no** |
| R16 | **CORNER** binary RD, gaussian forced, hr.threshold = 1.5 (legacy spelling), use_dina | 1.5 | 0.05 | **no** (intended) | 0.05 / 0.05 | **no** |
| R17 | binary RD, gaussian forced, effect.threshold = 0.07 (in range), dina | 0.07 | 0.07 | yes | 0.07 / 0.07 | yes |
| R18 | binary RD, gaussian forced, effect.threshold = 0.07 (in range), use_dina | 0.07 | 0.07 | yes | 0.07 / 0.07 | yes |
| R19 | binary RD, gaussian forced, both unset, dina | 0.05 | 0.05 | yes | 0.05 / 0.05 | yes |
| R20 | binary RD, gaussian forced, both unset, use_dina | 0.05 | 0.05 | yes | 0.05 / 0.05 | yes |
| R21 | count IRD, gaussian forced, effect.threshold = 0.02 (in range), dina | 0.02 | 0.02 | yes | 0.02 / 0.02 | yes |
| R22 | count IRD, gaussian forced, effect.threshold = 0.02 (in range), use_dina | 0.02 | 0.02 | yes | 0.02 / 0.02 | yes |
| R23 | binary OR, default binomial family, effect.threshold = 1.5, dina | 0.405465108108164 (= log 1.5) | 0.405465108108164 | yes | 0.405465108108164 / 0.405465108108164 | yes |
| R24 | binary OR, default binomial family, effect.threshold = 1.5, use_dina | 0.405465108108164 | 0.405465108108164 | yes | 0.405465108108164 / 0.405465108108164 | yes |
| R25 | survival HR, default cox family, both unset, dina | 0.22314355131421 (= log 1.25) | 0.22314355131421 | yes | 0.22314355131421 / 0.22314355131421 | yes |
| R26 | survival HR, default cox family, both unset, use_dina | 0.22314355131421 | 0.22314355131421 | yes | 0.22314355131421 / 0.22314355131421 | yes |

- **`identical()` is TRUE for all 20 non-corner rows**: fit floor, replay floor, whole fit object and whole replay object. The script prints `Non-corner rows that differ: none` and `UNAFFECTED GATE: PASS`.
- **The 6 corner rows** (R07, R08, R13-R16) differ only in the floor and what follows from it:
  - `dina` route: `floors`, `fit$grp.consistency`, `fit$df.est`, `fit$sg.harm`.
  - `use_dina` route: `floors`, `fit$dina_cuts`.
  - The recorded `args_call_all$effect.threshold` and the replay floor are the same on both trees (0.05 or 0.01).
- **Wall clocks:** start 60.58 s, fixed 64.04 s, compare 1.65 s. A trial start run of 61.03 s was discarded (notes).
- **Verifier re-run** (`verify/F14/unaff_compare.txt`), on its own `startcopy` and `post` trees:
  - Same 26-row table, `UNAFFECTED GATE: PASS`.
  - Wall clocks: start 60.84 s, post 64.29 s (`unaff_run_start.txt`, `unaff_run_post.txt`).
- **Verifier extra rows** (`verify/F14/extra_compare.txt`), each on both routes: 18 rows, walls 18.46 s (start) and 20.31 s (post).
  - The floor is unchanged in X1-X7:
    - binary OR and count IRR with gaussian forced, 1.5;
    - RD at the 1.0 boundary;
    - RD at −0.02;
    - MD with `hr.threshold = 1.25` explicit;
    - survival with gaussian forced;
    - survival on cox with `hr.threshold = 1.5`.
  - Only the extra corner rows change:
    - X8, IRD with `hr.threshold = 1.25` explicit: 1.25 → 0.01;
    - X9, RD with 1.5 and `consistency.threshold = 0`: 1.5 → 0.05.
- **Cumulative re-run.** The main session ran this after the fixes. It is outside F14's per-fix gate, and its files are in `cumulative/F14/`.
  - Against `head_fixes` (`.sha` `61c15088`, all 14 fixes), every non-corner fit floor and replay floor is unchanged. The script printed `UNAFFECTED GATE: FAIL`, because the whole objects of R02, R05 and R06 differ (`unaffected_cmp.log`).
  - Those three are the continuous MD `use_dina` rows, and they differ in `dina_cuts` and the search that follows from them.
  - With F11 removed from that tree (`cumulative/F14/minus_F11`), every row is identical to the gated fixed tree (`minusF11_vs_gated_fixed.txt`) and the gate prints PASS (`unaffected_cmp_minus_F11.log`).
  - So the R02/R05/R06 differences come from F11's `tau_sign` on the `use_dina` route (head `R/forestsearch_main.R:3301`), not from F14.
  - Wall clocks: head 63.84 s, start re-run 62.49 s, minus-F11 65.04 s, each compare about 1.6 s.

#### Existing tests

**Implementer** (`evidence/F14/existing.txt`):

- Which files: 52 files on `wt/F14`, one process each, `timeout 120`, 8 at a time. These are every file that mentions DINA, `test-threshold-sync.R`, and every other file that calls `forestsearch(`.
- 51 files finished with rc=0 and 0 failures/errors.
- `test-mr-inference.R` hit the 120 s cap (rc=124, wall 120.1 s, no failure before the cap). It was re-run one `test_that()` block per process: 21 of 21 blocks passed with 0 failures, the longest taking 69.2 s.
- 5 blocks skipped, in 4 files:
  - `test-binary-default-or.R:188`, `test-directive-c.R:667` and `test-threshold-pair-directive-a.R:148` need a baseline under `dev/` or a template under `quarto/`. On a scratch copy with those files copied from the start export, all three ran with 0 skips and 0 failures (10.3 s, 8.8 s, 7.4 s).
  - `test-cross-outcome-parity.R:262` and `:358` skip with the reason "primary fit did not identify; test N/A". They were not re-run.

**Verifier** (`verify/F14/existing/summary_post.txt`):

- 51 files on `post`, 12 at a time, 57 s for the batch.
- Totals: nb 5456, failed 0, error 0, skipped 5.
- `directive-c` (182 expectations), `threshold-pair-directive-a` (86) and `binary-default-or` (66) passed with 0 skips on a copy that held the baselines.
- `test-mr-inference.R` was not run (no DINA route).

**Existing expectations updated** (only this one block, `tests/testthat/test-dina-gaussian-floor.R`):

- `:144` → `:150`: `expect_identical(unique(fit$floors), 1.5, info = route)` → `expect_identical(unique(fit$floors), 0.05, info = route)`.
- Title `:123` → `:124`: "... that is remapped for RD is still the floor as given" → "... is the floor as remapped".
- Comment `:124-126` ("a supplied threshold leaves the floor unchanged from before A-001") → `:125-129`, plus a new derivation comment at `:147-149`.
- Nothing else. The verifier's old-file run on the fixed code fails only at `:144` (2 failures).

#### What a user sees

The call is `forestsearch(..., outcome_type = "binary", effect_measure = "RD", subgroup_method = "dina", dina_args = list(family = "gaussian"), effect.threshold = 1.5, details = TRUE)` on the test file's binary data (`set.seed(7)`, N = 300).

- Both trees first print the same two warnings (`verify/F14/console_example.txt`). Their first lines are:
  - "effect.threshold = 1.50 but effect_measure = 'RD' (identity scale, range [-1, 1])."
  - "effect.threshold = 1.50 appears to be on a ratio scale, ...". This warning goes on to say "Remapping to default: 0.05."
- The blocks below are an excerpt of the `details = TRUE` output plus two printed fields, as captured by `evidence/F14/console_example.R` (`evidence/F14/console_example.txt`).

Start SHA:
```
DINA selection: no subgroup met the criteria
  Harm floor:          m_diff = 1.5000
fs$sg.harm: NULL
fs$args_call_all$effect.threshold: 0.05
```
Fixed tree:
```
DINA-selected subgroup: {x1 <= 2.75}
  Harm floor:          m_diff = 0.0500
fs$sg.harm: [1] "{x1 <= 2.75}"
fs$args_call_all$effect.threshold: 0.05
```

- Before the fix, the fit used 1.5, which a risk difference cannot reach, so it found nothing. Its bootstrap would re-select with 0.05.
- After the fix, the fit selects the same subgroup as the call with both thresholds unset (`unaffected.txt` R19: `sg.harm [{x1 <= 2.75}]`).
- Wall clocks: start 2.99 s, fixed 4.82 s. The verifier's re-run with warnings captured took 3.05 s and 4.76 s.

#### Verification

**Verdict: CONFIRMED.** The verifier reproduced red, green and the unaffected case. There was no repair round. It raised three problems, all minor:

| # | severity | problem | resolution |
|---|---|---|---|
| 1 | minor | The new test's comment (`:159-161`) cites only the "Resolved defaults, per estimand" table for RD 0.05 and IRD 0.01. That table is stated for "when neither threshold is supplied". The remap of a supplied out-of-range value is documented only in the warning text, which for IRD names 0.05. | **Left as noted.** The committed test equals the gated patch (`patch_equal.txt` empty), and the comment was not extended. Listed under Open items. |
| 2a | minor | `?forestsearch` (`man/forestsearch.Rd:1212-1214`; roxygen `R/forestsearch_main.R:1313-1315`) says a supplied threshold "suppresses the corresponding remap and the supplied value is used". This is inexact in the corner. | **Open.** `man/forestsearch.Rd` is unchanged from `76a62a71` to `88302212`. Listed under Open items. |
| 2b | minor | The part-3 NEWS bullet for A-001 ("A supplied threshold, and every non-gaussian DINA family, give the same floor as before") becomes inexact. | **Fixed in the NEWS commit `88302212`.** The F14 entry at `NEWS.md:390-396` says "This narrows the A-001 entry above: every other supplied threshold gives the same floor as before." The A-001 bullet is at `NEWS.md:172-179` in the same "Behaviour changes" section. |
| 3 | minor | The notes say "the same warning". There are two, and the console block is a grep excerpt. | **Fixed in this section:** both warnings are quoted under "What a user sees", and the block is labelled as an excerpt. |

#### Departures and interpretations

- **Test helper.** `.gf_replay()` gained `df = .gf_df()` (`:63`, `:65`), so the corner test can replay binary and count fits. Existing callers pass no `df` and behave as before. This is test-only.
- **Where the new test lives.** It was added to `tests/testthat/test-dina-gaussian-floor.R`, the file the spec names for the pin update.
- **"Bootstrap and CV replays".** These were checked with the bootstrap's own `do.call(forestsearch, args_call_all)` pattern, as the existing test file already does. `forestsearch_bootstrap_dofuture()` and `forestsearch_Kfold()` themselves were not run.
- **"Every other threshold configuration the existing test file covers".** This was read as the configurations in the file at the start SHA: R01-R10, including the two direct `.forestsearch_dina_select()` rows (R09, R10).
- **"Continuous outcome with a supplied threshold".** This was read as MD with a supplied 1.5 (R11-R12), a value MD never remaps. R03-R06 already cover a supplied 0.4.
- **Corner membership.** The legacy spelling `hr.threshold = 1.5` (R15, R16) and the verifier's X8/X9 are counted as corner rows, because the same RD/IRD remap fires there.
- **No roxygen or `man/` change.** `?forestsearch` `dina_args` (`man/forestsearch.Rd:233-234`) already said the gaussian floor comes "from the resolved \code{effect.threshold}". The fix brings the code in line with that text. `devtools::document()` was not run for this fix.
- **`test-mr-inference.R`.** It ran one block per process instead of one process, because the spec's 2-minute cap could not be met for the whole file.
- **Skip-gated blocks.** They ran on a scratch copy holding three read-only files from the start export. Nothing was added to the working tree.

#### paper_simulations

- Does not reach the changed lines. The only evaluated chunks that call into `forestsearch()` are `run-alt-sims` (`:908`) and `run-null-sims` (`:966`) of `vignettes/articles/paper_simulations.qmd`. They call `run_simulation_analysis()` with no `methods` and with `fs_params` that set none of `use_dina`, `subgroup_method` or `dina_args`.
- So `forestsearch()` runs survival with `use_dina = FALSE` and `subgroup_method = "consistency"`. The article contains no "dina" text (`git grep` count 0 at `61c15088`).

#### Open items

- **IRD remap warning names the wrong default.** `R/forestsearch_main.R:2334-2350` (start `76a62a71`, same lines at head `61c15088`) passes a literal `0.05` and `"OR"` into the "Remapping to default: %.2f" warning (`:2346`) for both RD and IRD, but IRD is remapped to 0.01 (`:2349`). The earlier validator warning (`:301-309`) also says "A risk difference cannot exceed 1.0" and "Did you mean effect_measure = 'OR'" for IRD. This is message text only, and the value is correct. It is outside F14, whose spec changes only the floor.
- **`?forestsearch` threshold-vocabulary sentence.** `man/forestsearch.Rd:1212-1214` (roxygen `R/forestsearch_main.R:1313-1315`, start and head) says a supplied threshold "suppresses the corresponding remap and the supplied value is used". For an out-of-range RD/IRD value (> 1), the value is instead remapped with a warning, and since F14 the remapped value is also the DINA floor. The sentence was not changed because F14 changed no roxygen. It needs a roxygen edit and `devtools::document()`.
- **Test comment derivation.** `tests/testthat/test-dina-gaussian-floor.R:159-161` (at `685e0c5d`) cites only the no-threshold defaults table for RD 0.05 and IRD 0.01. It does not spell out that the remap goes "to default", or that the IRD warning misnames 0.05. The values are correct (verifier problem 1). The comment was not changed because the committed patch is the gated one.
- **Helper default kept.** `.forestsearch_dina_select()` keeps `effect.threshold = hr.threshold` as its default (`R/forestsearch_helpers.R:1489`, start and head) for direct callers such as `quarto/dina/smoke-tests/validate_screen3.R`. `forestsearch()` always passes the resolved value (R09/R10 unchanged). This was outside the spec.
- **Gaussian forced on ratio estimands.** Survival with `dina_args = list(family = "gaussian")` uses the raw 1.25 as an identity-scale floor, and OR/IRR with gaussian forced use the raw 1.5 (verifier `extra_compare.txt` X1, X2, X6). This is the same on both trees and was not in F14's scope.
- **Real resampling not run.** `forestsearch_bootstrap_dofuture()` and `forestsearch_Kfold()` were not run in the corner. The replay check uses their `args_call_all` pattern.
- **`test-mr-inference.R` exceeds 120 s.** The file cannot finish in one 120 s process (rc=124). It passes block by block (21/21). Its run time is outside F14.

#### Corrections to the notes

- The notes say the spec pointers at `28e71ffa` were "2784-2785 and 3250-3252". The spec says "near 2784–2785 and 3251–3253", and at `28e71ffa` the `m_diff_sel` block is at `R/forestsearch_main.R:3251-3253` (guard line `:3252`).
- The notes give the pinning block as `:124-150` after. In git at `685e0c5d` the block is `:124-153`, with the updated expectation at `:150`.
- The notes place the `.gf_replay()` change at `:62-66`. In git the changed lines are `:61-63` and `:65`.
- The notes say three files skip one block each. `existing.txt` also shows `test-cross-outcome-parity.R` skipping 2 blocks (`:262`, `:358`, "primary fit did not identify"), so there are 5 skipped blocks in 4 files.
- The notes say "Both trees first print the same warning". There are two warnings (`verify/F14/console_example.txt`).

