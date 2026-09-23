# REPORT -- `R CMD check --as-cran` cleanup (2026-09-23)

Task: `dev/tasks/TASK_cran_check_cleanup_2026-09-23.md` (committed `4388ea80`).
Outcome: **2 ERRORs, 1 WARNING, 2 NOTEs -> 1 NOTE.** No new ERROR, WARNING or NOTE.
The one remaining NOTE (CRAN incoming: `Version contains large components
(0.3.5.9000)`) is **deliberately left open** (Larry, 2026-09-23): it is the standard
marker for an unreleased development version and is expected. It clears at
submission with a release version bump, not in this task. It is not an unresolved item.

## Environment

| item | value |
|---|---|
| HEAD at start | `4388ea80` (task doc commit; parent `114dcba8`, the digits task's last commit) |
| installed forestsearch | 0.3.5.9000, built 2026-09-23 02:32 UTC (predates the digits commits; not used -- both checks build and install from a clean export) |
| R | 4.6.1 (2026-06-24), x86_64-pc-linux-gnu |
| machine | pop-os, Linux 7.1.5 |
| HTML Tidy | 5.6.0 (`/usr/bin/tidy`) |
| `DESCRIPTION` `Encoding:` | `UTF-8` (already declared) |
| pre-existing untracked (never staged) | `quarto/simulations/actg175/binary_020/mr_or_harm/fs_effMaxSG_mr_field_or075_n500_nb20_{redes,relaunch}_d5000/`, `quarto/simulations/actg175/binary_020/smoke_{redes,relaunch}.html`, `quarto/simulations/gbsg_020/scripts_dinamr/logs/nullmr_findings.err` |

Check surface: `rcmdcheck::rcmdcheck(args = "--as-cran")` with
`RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64`,
on a clean export (`git archive HEAD` for the baseline; the tracked files of the
modified tree for the post-change run). The PDF manual was built in both runs.

## Amendment to the committed task doc (from Larry, 2026-09-23, at launch)

This is not in the committed task doc.

1. **tidy is now installed.** The HTML tidy NOTE (§2 item 6) no longer applies.
2. **Fresh baseline.** §1's baseline was taken fresh at current HEAD with tidy
   present. The digits task's saved check output, which predates both, was
   **not** reused.
3. **§2 item 6 replaced:** report whatever the HTML validation now says. Any
   finding it surfaces is in scope under §3's boundary rule, since it would be a
   roxygen/Rd issue, not function logic.

## §1 Baseline (HEAD `4388ea80`, run 14:03-14:14)

`Status: 2 ERRORs, 1 WARNING, 2 NOTEs`. Full item list verbatim:

```
checking examples with --run-donttest ... [19s/17s] ERROR
Running examples in 'forestsearch-Ex.R' failed
The error most likely occurred in:

> ### Name: fs_dgm_feasibility
> ### Title: Design-time feasibility of a DGM's planted region
> ### Aliases: fs_dgm_feasibility print.fs_dgm_feasibility
>
> ### ** Examples
>
> ## No test:
> dgm <- generate_glm_dgm(n_super = 2000, outcome_type = "binary",
+                         effect_measure = "OR", seed = 8316951)
Error in generate_glm_dgm(n_super = 2000, outcome_type = "binary", effect_measure = "OR",  :
  argument "data" is missing, with no default
Calls: generate_glm_dgm -> stopifnot -> is.data.frame
Execution halted
```

```
checking tests ...
  Running 'testthat.R' [348s/294s]
 [348s/295s] ERROR
Running the tests in 'tests/testthat.R' failed.
Last 13 lines of output:
      ▆
   1. ├─base::paste(...) at test-fs-dgm-feasibility.R:103:3
   2. └─base::readLines(...)
   3.   └─base::file(con, "r")
  ── Failure ('test-fs-family-report.R:115:3'): coverage guard: every forestsearch() formal is classified or explicitly out of scope ──
  Expected `intersect(classified, .FR_OUT_OF_SCOPE)` to be identical to `character(0)`.
  Differences:
  `actual`:   "effect_measure"
  `expected`:

  [ FAIL 2 | WARN 22 | SKIP 77 | PASS 5299 ]
  Error:
  ! Test failures.
  Execution halted
```

```
checking code files for non-ASCII characters ... WARNING
Found the following file with non-ASCII characters:
  R/fs_bias_coverage.R
Portable packages must use only ASCII characters in their R code and
NAMESPACE directives, except perhaps in comments.
Use \uxxxx escapes for other characters.
Function 'tools::showNonASCIIfile' can help in finding non-ASCII
characters in files.
```

```
checking CRAN incoming feasibility ... [6s/15s] NOTE
Maintainer: 'Larry Leon <larry.leon.05@post.harvard.edu>'

Version contains large components (0.3.5.9000)
```

```
checking R code for possible problems ... [41s/41s] NOTE
fs_plot_bias_coverage: no visible binding for global variable 'b'
fs_plot_bias_coverage: no visible binding for global variable 'cov1'
fs_plot_bias_coverage: no visible binding for global variable 'cov'
fs_plot_bias_coverage: no visible binding for global variable 'r'
fs_plot_bias_coverage: no visible binding for global variable
  'estimator'
fs_plot_bias_coverage: no visible binding for global variable 'cell'
fs_plot_bias_coverage: no visible binding for global variable 'cov2'
fs_plot_bias_coverage: no visible binding for global variable 'ref'
fs_plot_bias_coverage: no visible binding for global variable 'obs'
Undefined global functions or variables:
  b cell cov cov1 cov2 estimator obs r ref
Consider adding
  importFrom("stats", "cov")
to your NAMESPACE file.
```

HTML validation (amended item 6):

```
* checking PDF version of manual ... OK
* checking HTML version of manual ... [15s/15s] OK
```

The task doc expected 3 NOTEs, the third being HTML tidy. With tidy installed
there is no HTML NOTE and no finding.

## §2 Diagnosis

1. **ERROR: `fs_dgm_feasibility` example.** Root cause in the **example**
   (`R/fs_dgm_feasibility.R` roxygen `@examples`). It called
   `generate_glm_dgm()` without its four required arguments `data`,
   `factor_vars`, `outcome_var` and `treatment_var` (`R/generate_glm_dgm.R:179`).
   Because the example sits inside `\donttest{}`, only `--run-donttest` reaches it.
2. **ERROR: tests.**
   - `test-fs-dgm-feasibility.R:103`: the test reads
     `test_path("..", "..", "R", "fs_dgm_feasibility.R")`. Under `R CMD check`
     the tests run against the installed package, and the `R/` source tree is
     not on disk, so `file()` fails. It passes under `devtools::test()`. Root
     cause in the **test**.
   - `test-fs-family-report.R:115`: `fs_family_report()` **does** classify
     `effect_measure`: `R/fs_family_report.R:315` adds the "estimability
     boundary" row with `arguments = "effect_measure"` (commit `1719056f`). The
     test's `.FR_OUT_OF_SCOPE` list still carried it. Root cause in the
     **test's list**.
3. **WARNING: non-ASCII in `R/fs_bias_coverage.R`.** All occurrences are
   **inside string literals that reach user-visible output** (the `caption =` of
   `fs_plot_bias_coverage()`'s panels 1 and 2). None are in comments or roxygen.

   | line | characters |
   |---|---|
   | 238 | `Φ` U+03A6, `·` U+00B7, `∈` U+2208 |
   | 239 | `−` U+2212 |
   | 251 | `Φ` U+03A6 ×2, `·` U+00B7 ×2, `−` U+2212 ×4 |

   `DESCRIPTION` already declares `Encoding: UTF-8`. That does not exempt R
   code from the ASCII rule.
4. **NOTE: CRAN incoming feasibility.** Sub-items: (a) the `Maintainer:` line
   is informational and always printed, so nothing needs fixing. (b) `Version
   contains large components (0.3.5.9000)` is fixable in the package only by
   changing the version to a release number. That is a `DESCRIPTION` version
   bump, which needs Larry's approval, and it is expected on a development
   version. **Deliberately left open** (Larry, 2026-09-23): the `.9000` suffix is the
   standard unreleased-development marker; the NOTE clears at submission with the
   release bump, not now. The "new submission" and "days since last
   update" sub-items the task doc anticipated do not appear in this output.
5. **NOTE: globals in `fs_plot_bias_coverage`.** Variables: `b cell cov cov1
   cov2 estimator obs r ref`. These are ggplot2 `aes()` columns of `tbl`,
   `cur1`/`cur2` and `long`. `cov` is the data-frame column, not
   `stats::cov`, so the check's `importFrom("stats", "cov")` suggestion is wrong
   for this case. Convention: `R/globals.R` (consolidated; CLAUDE.md).
6. **HTML validation.** `OK`, with no findings (see amendment).

**`devtools::document()` error at `R/fs_bias_coverage.R:17`:**

```
✖ fs_bias_coverage.R:17: @description failed to evaluate inline markdown code.
Caused by error:
! Failed to parse the inline R code: `r = se_mean / sd_emp`
```

Line 17 is the start of the block. The offending text is on line 22,
`` `r = se_mean / sd_emp` ``: roxygen markdown treats an inline code span
beginning with `r ` as inline R to evaluate. It is **not** related to item 3
(encoding). It had kept `man/` from regenerating for the whole package. The
committed `man/fs_sim_bias_coverage.Rd` shows the intended rendering,
`\code{r = se_mean / sd_emp}`.

## §3 Fixes

No exported function's logic changed. No item's root cause lay in function logic.

| file | change |
|---|---|
| `R/fs_bias_coverage.R:22` | `` `r = se_mean / sd_emp` `` -> `\code{r = se_mean / sd_emp}` (Rd markup directly; rendered Rd byte-identical, `fs_sim_bias_coverage.Rd` unchanged by `document()`) |
| `R/fs_bias_coverage.R:238,239,251` | non-ASCII -> `\u` escapes, below |
| `R/globals.R` | new section: `"b", "r", "cell", "estimator", "cov", "cov1", "cov2", "obs", "ref"` |
| `R/fs_dgm_feasibility.R` `@examples` | builds a small synthetic binary data frame (N = 400) and passes it to `generate_glm_dgm()` with `factor_vars`, `continuous_vars`, `outcome_var`, `treatment_var`, `subgroup_vars`/`subgroup_cuts` (age > q55) and `k_inter = 1`; the `n_super = 2000`, `effect_measure = "OR"` and `seed = 8316951` of the original are kept. Still `\donttest{}`; runs in ~0.3 s |
| `tests/testthat/test-fs-dgm-feasibility.R:103` | `readLines(test_path("..","..","R",...))` -> `deparse(fs_dgm_feasibility)`. Same assertions; now scoped to the function rather than the whole file. Verified against the installed package (no srcref) |
| `tests/testthat/test-fs-family-report.R:35` | `"effect_measure"` removed from `.FR_OUT_OF_SCOPE` |
| `man/fs_dgm_feasibility.Rd` | regenerated |
| `NEWS.md` | one bullet |

**Non-ASCII substitutions (post-condition 6):**

| line | char | escape | occurrences |
|---|---|---|---|
| 238 | `Φ` | `Φ` | 1 |
| 238 | `·` | `·` | 1 |
| 238 | `∈` | `∈` | 1 |
| 239 | `−` | `−` | 1 |
| 251 | `Φ` | `Φ` | 2 |
| 251 | `·` | `·` | 2 |
| 251 | `−` | `−` | 4 |

Rendered text unchanged: all 97 string literals in the file, parsed and
evaluated, are `identical()` between HEAD and the modified file. No instance
stopped under the "rendered text would change" rule.

## §4 Verification

Post-change check (modified tree, run 14:15-14:26): `Status: 1 NOTE`.

| # | item | baseline | after | disposition |
|---|---|---|---|---|
| 1 | examples `--run-donttest` (`fs_dgm_feasibility`) | ERROR | OK `[215s/102s]` | **resolved** |
| 2 | tests (`test-fs-dgm-feasibility.R:103`, `test-fs-family-report.R:115`) | ERROR `[ FAIL 2 \| WARN 22 \| SKIP 77 \| PASS 5299 ]` | OK `[ FAIL 0 \| WARN 21 \| SKIP 77 \| PASS 5302 ]` | **resolved** |
| 3 | non-ASCII in `R/fs_bias_coverage.R` | WARNING | OK | **resolved** |
| 4 | CRAN incoming feasibility | NOTE | NOTE (identical text) | **deliberately left open**: standard `.9000` development-version marker; clears with the release bump at submission |
| 5 | globals in `fs_plot_bias_coverage` | NOTE | OK | **resolved** |
| 6 | HTML version of manual | OK | OK | n/a (amended: no finding) |
| - | `devtools::document()` inline-code error | error | clean | **resolved** |
| - | new items | - | none | - |

The same 56 check steps ran in both. The donttest examples take 215 s against
19 s because the baseline halted at the first failing example; every
`\donttest` example after `fs_dgm_feasibility` now runs, and all pass.

`devtools::test()` (full suite, source tree):

| | suite line |
|---|---|
| before (HEAD `4388ea80`) | `[ FAIL 1 \| WARN 32 \| SKIP 3 \| PASS 5710 ]` (the `effect_measure` coverage guard) |
| after | `[ FAIL 0 \| WARN 32 \| SKIP 3 \| PASS 5711 ]` |

No test that passed before now fails (matched by file × test name).
`test-fs-dgm-feasibility.R:103` passes in both under `devtools::test()`; its
failure showed only under `R CMD check`.

## Post-conditions

1. Six items, verbatim text and root cause: **met** (item 6 per amendment).
2. No new ERROR/WARNING/NOTE: **met**.
3. Before/after table with dispositions: **met**. None skipped under the 30-minute cap.
4. Suite line before/after; no regressions: **met**.
5. No exported function's logic changed: **met**.
6. Every non-ASCII substitution listed: **met**.
7. `man/` regenerated where roxygen changed: **met** (`fs_dgm_feasibility.Rd`; `fs_sim_bias_coverage.Rd` regenerated byte-identical).
8. Pre-existing untracked files never staged: **met**.
