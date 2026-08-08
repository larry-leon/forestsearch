# REPORT: blockers 1–3 — eval/theta closeout, MR-guard trace, 0/50 diagnosis

Task document: `CC_TASK_md_mr_blockers_1-3.md`.
Authority: `dev/glm-continuous-sims/HANDOFF_md_mr_harness_session.md` (v2),
`dev/betaHhat-consolidation/SPEC_eval_theta_functions.md`,
`dev/glm-continuous-sims/NOTE_threshold_sign_md.md`.

**A Task 1 stop condition fired.** `R CMD check --as-cran` disagreed with the
certified line. Per Task 1 step 3, nothing from the closeout was committed, no
debugging toward green was attempted, no test was touched, and Tasks 2–3 were
halted. This report contains what exists. Blocker 4 (re-pilot / STOP A) did not
start. No `R/` file and no harness file was modified.

---

## Header — Task 0 preconditions

`git fetch origin` ran clean.

```
$ git rev-parse HEAD
594346907857fa6466d675eeabd81afeef5728e2

$ git rev-parse origin/feature/mr-in-replicates
594346907857fa6466d675eeabd81afeef5728e2

$ git status --porcelain
 M NAMESPACE
 M R/betaHhat_truth.R
 M tests/testthat/test-betaHhat-contract.R
?? dev/version2-preparation/
?? man/fs_betaHhat_theta_dagger_check.Rd
?? man/fs_build_eval_frame.Rd
?? quarto/simulations/actg175/continuous/mr_coverage_sweep_md_harm.qmd
?? quarto/simulations/actg175/continuous/mr_sweep_md_harm/
```

Local `==` remote (both `5943469`). MODIFIED files are exactly the three the
ground truth expects and no others. **No Task 0 stop condition fired.**

Untracked scratch beyond the harness pair, listed and **not touched**:

- `dev/version2-preparation/` — one file, `RELEASE_STATUS.md`, 7005 bytes,
  mtime `2026-08-06 21:28`.

The two new `man/` pages and the harness pair
(`mr_coverage_sweep_md_harm.qmd`, `mr_sweep_md_harm/`) are the expected
untracked entries.

---

## 1 — Eval/theta closeout (blocker 1) — STOPPED

### 1.1 Tests: PASS

`devtools::test(filter = "betaHhat-contract")`, verbatim:

```
ℹ Testing forestsearch
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 159 ]
```

Re-run with a per-test reporter so T11–T14 are individually identifiable
(`nb` = assertions):

```
                                                           test nb failed skipped error warning
    T1: nH_eval + nHc_eval == N for every rule shape and family 30      0   FALSE FALSE       0
       T1 negative control: partial-NA membership FAILS the sum  3      0   FALSE FALSE       0
 T2: a rule naming a missing column gives an all-NA unresolved  18      0   FALSE FALSE       0
   T2 negative control: the unguarded path yields NA membership  4      0   FALSE FALSE       0
     T3: a realized disjunction and its structured sg_def agree  6      0   FALSE FALSE       0
      T3 negative control: split-first FAILS on the same string  3      0   FALSE FALSE       0
 T4: negation is equivalent to the flipped comparison, and part  3      0   FALSE FALSE       0
         T5: no subgroup gives nHc_eval == N and the ITT effect 12      0   FALSE FALSE       0
 T5 (attach): a results frame with no sg_def column takes the I  4      0   FALSE FALSE       0
 T5 negative control: a real subgroup must NOT give nHc_eval ==  2      0   FALSE FALSE       0
   T8: the same rule gives identical membership across families 10      0   FALSE FALSE       0
  T9: repeated evaluation of a resolved target is bit-identical  6      0   FALSE FALSE       0
           focus is required, validated, and inverts the region  4      0   FALSE FALSE       0
 .fs_rule_columns finds the referenced columns in every rule sh  4      0   FALSE FALSE       0
 T6: counters match the unresolvable rules in a synthetic bundl 20      0   FALSE FALSE       0
    T6 negative control: a clean bundle reports zero unresolved  3      0   FALSE FALSE       0
                   T7: the parity guard passes on a clean table  2      0   FALSE FALSE       0
              T7: the guard FIRES on an injected dropped target  2      0   FALSE FALSE       0
  T7: the guard FIRES when a C_betaHhat row is missing entirely  1      0   FALSE FALSE       0
                      T7: strict = FALSE warns without stopping  2      0   FALSE FALSE       0
 T7: a table without the required columns is a no-op, not an er  1      0   FALSE FALSE       0
 T7 companion: the report separates unresolvable from degenerat  5      0   FALSE FALSE       0
    T11: frames are bitwise identical to the shims they replace  4      0   FALSE FALSE       0
  T12: the n_eval trap fires; survival-only args error off-path  4      0   FALSE FALSE       0
                    T13: theta-dagger matches the shims exactly  3      0   FALSE FALSE       0
    T14 negative control: the wrong outcome_type must NOT agree  3      0   FALSE FALSE       0

TOTALS: tests=26 assertions=159 failed=0 skipped=0 error=0 warning=0
```

T11 (4 assertions), T12 (4), T13 (3), T14 (3) all pass, individually
identified. The test names match `SPEC_eval_theta_functions.md`'s T11–T14
definitions (frame identity; the `n_eval` trap and off-path argument rejection;
theta identity to the shims; the wrong-`outcome_type` negative control).

### 1.2 Check: DISAGREES with the certified line — STOP

Invocation, the first option the task document offers:
`rcmdcheck::rcmdcheck(args = "--as-cran", error_on = "never")`.

Final status line, **verbatim**:

```
0 errors ✔ | 1 warning ✖ | 2 notes ✖
```

Expected per the refined ground truth:

```
0 errors ✔ | 0 warnings ✔ | 0 notes ✔
```

**This is a hard stop.** Per Task 1 step 3: nothing from the closeout was
committed; no fix was attempted; no test was touched.

Verbatim results header:

```
── R CMD check results ───────────────────────────────── forestsearch 0.2.0 ────
Duration: 11m 15.6s
```

The one WARNING, verbatim (the LaTeX error block repeats; 62 `LaTeX Error`
lines in total, over exactly two distinct characters):

```
❯ checking PDF version of manual ... [12s/12s] WARNING
  LaTeX errors when creating PDF version.
  This typically indicates Rd problems.
  LaTeX errors found:
  ! LaTeX Error: Unicode character ─ (U+2500)
                 not set up for use with LaTeX.

  See the LaTeX manual or LaTeX Companion for explanation.
  Type  H <return>  for immediate help.
```

The two distinct offending characters, verbatim:

```
Unicode character ᶜ (U+1D9C)
Unicode character ─ (U+2500)
```

The two NOTEs, verbatim:

```
❯ checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found.
  Please obtain a recent version of HTML Tidy by downloading a binary
  release or compiling the source code from <https://www.html-tidy.org/>.

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ‘forestsearch-manual.tex’
```

The test stage inside the check passed:

```
* checking tests ...
  Running ‘testthat.R’ [237s/187s]
 [237s/187s] OK
```

### 1.3 What the divergence is and is not — observations only, nothing acted on

Recorded because the task document asks for the divergence as a finding. **No
change was made in response to any of this.**

1. **None of the three findings involves any file under closeout.** The
   offending characters trace to tracked, unmodified `man/` pages:

   ```
   $ grep -rl "─\|ᶜ" man/
   man/fpr_calibration.Rd
   man/sg_tables.Rd

   $ git log --oneline -1 -- man/fpr_calibration.Rd man/sg_tables.Rd
   8a7096c extreme subgroups vignette draft

   $ git status --porcelain man/fpr_calibration.Rd man/sg_tables.Rd
   (empty — tracked and unmodified)
   ```

   The two new man pages are pure ASCII (`grep -c '[^ -~]'` returns `0` for
   both `man/fs_build_eval_frame.Rd` and
   `man/fs_betaHhat_theta_dagger_check.Rd`).

2. **The two invocations the task document offers as equivalent are not
   equivalent on exactly these items.** `formals(devtools::check)$manual` is
   `FALSE` (devtools 2.5.2), so `devtools::check(args = "--as-cran")` passes
   `--no-manual` and skips the PDF-manual check entirely — which would suppress
   the WARNING, and with it the `forestsearch-manual.tex` byproduct that
   produces the second NOTE. `rcmdcheck::rcmdcheck(args = "--as-cran")`
   (rcmdcheck 1.4.0) does not skip it. Which invocation produced `check5.log`
   is not recorded in the handoff.

3. `tidy` is not on `PATH` (`Sys.which("tidy")` is empty); `pdflatex` is
   (`/home/larryleon/bin/pdflatex` → TinyTeX, symlink dated 2026-05-15, so the
   LaTeX toolchain predates this work).

4. Timings differ substantially from `check5.log` (`Duration: 14m 46.8s`,
   `checking tests ... [654s/484s]`) — here `11m 15.6s` and `[237s/187s]`.
   Machine load, not a package property.

**Not established:** whether `check5.log`'s green line came from a `--no-manual`
invocation. Confirming that would require re-running the check a second way,
which is debugging toward green, so it was not done. Deciding this is the
reviewer's call.

### 1.4 Closeout state

- Commit made: **none**. `R/betaHhat_truth.R`,
  `tests/testthat/test-betaHhat-contract.R`, `NAMESPACE`,
  `man/fs_build_eval_frame.Rd`, `man/fs_betaHhat_theta_dagger_check.Rd` remain
  uncommitted, **unmodified by this session**.
- Man page filenames, for the record when the closeout does commit:
  `man/fs_build_eval_frame.Rd`, `man/fs_betaHhat_theta_dagger_check.Rd`.
- `devtools::install()` (Task 1 step 6) was **not** run — it is conditional on
  the commit.
- The engine re-point stays deferred and queued. Not started, not prepared.

---

## 2 — The MR guard, from source (blocker 2)

Task 2 is investigation-only and was completed from source before the Task 1
stop was known. It required no execution beyond triggering the guard, and no
code was changed. Line numbers re-verified on the local tree; the working tree
differs from remote `5943469` only in the three Task-1 files, none of which is
`R/forestsearch_main.R` or `R/forestsearch_helpers.R`.

### 2.a The guard, verbatim

**Gate — `R/forestsearch_main.R:1292–1295`** (matches the task document's
pointer exactly), placed before any fitting:

```r
  if (isTRUE(mr_inference)) {
    .validate_mr_configuration(args_call_all,
                               context = "forestsearch(mr_inference = TRUE)")
  }
```

Its preceding comment (`R/forestsearch_main.R:1280–1291`) states the intent:
"Multiplier resampling linearizes the selection event, which is valid only when
selection ranks on the inferential coefficient and the candidate family is
fixed. Checked here, immediately after the capture, so an unusable
configuration fails before any fitting."

**Validator — `.validate_mr_configuration()`, `R/forestsearch_helpers.R:1861`**
(matches). Its `fail()` at **`:1873–1876`** (matches):

```r
  fail <- function(msg) {
    stop(sprintf("%s requires an identifier aligned with the reported effect.\n  %s",
                 context, msg), call. = FALSE)
  }
```

**Condition 2 — `R/forestsearch_helpers.R:1909–1927`** (task document said
`~:1910–1927`; the section banner starts at `:1909`):

```r
  # ---------------------------------------------------------------------------
  # Condition 2: fixed candidate family (consistency / FS only)
  # ---------------------------------------------------------------------------
  if (identical(method, "consistency")) {
    front_ends <- c(use_lasso = isTRUE(.first("use_lasso", FALSE)),
                    use_grf   = isTRUE(.first("use_grf",   FALSE)),
                    use_dina  = isTRUE(.first("use_dina",  FALSE)))
    on <- names(front_ends)[front_ends]
    if (length(on)) {
      fail(sprintf(paste0(
        "The model-based front end%s %s set TRUE, so the candidate family is ",
        "estimated from the same data the search runs on and is not fixed.\n",
        "  Set %s and re-run."),
        if (length(on) > 1L) "s" else "",
        paste(paste(on, collapse = " and "),
              if (length(on) > 1L) "are" else "is"),
        paste(sprintf("%s = FALSE", on), collapse = ", ")))
    }
  }
```

**Exact `stop()` text**, produced by throwaway `forestsearch()` calls
(`subgroup_method = "consistency"`, `outcome_type = "continuous"`,
`effect_measure = "MD"`, `mr_inference = TRUE`, n = 60 synthetic rows). The
guard fires before any fitting, so each call returns essentially instantly.

`use_lasso = TRUE`:

```
forestsearch(mr_inference = TRUE) requires an identifier aligned with the reported effect.
  The model-based front end use_lasso is set TRUE, so the candidate family is estimated from the same data the search runs on and is not fixed.
  Set use_lasso = FALSE and re-run.
```

`use_grf = TRUE`:

```
forestsearch(mr_inference = TRUE) requires an identifier aligned with the reported effect.
  The model-based front end use_grf is set TRUE, so the candidate family is estimated from the same data the search runs on and is not fixed.
  Set use_grf = FALSE and re-run.
```

Both TRUE (recorded because this is the harness's actual configuration —
see §3):

```
forestsearch(mr_inference = TRUE) requires an identifier aligned with the reported effect.
  The model-based front ends use_lasso and use_grf are set TRUE, so the candidate family is estimated from the same data the search runs on and is not fixed.
  Set use_lasso = FALSE, use_grf = FALSE and re-run.
```

### 2.b Classification: STRUCTURAL

Evidence from the MR machinery, not the validator's comments. The chain from
family construction to consumption:

**The family is enumerated once, before MR is called.**
`R/forestsearch_main.R:2769` builds the cut matrix from the screened
confounders:

```r
  Z <- as.matrix(df.confounders)
```

and the MR block enumerates the candidate family from that fixed `Z` at
`R/forestsearch_main.R:3182–3197`, a single `for (kk in seq_len(tot_mr))` loop
producing the named list `fam`, before `fs_mr_inference()` is called at
`R/forestsearch_main.R:3218`. Its own comment (`:3176–3181`) states the family
is built "with the search's OWN helpers, so the family is identical to the
space `subgroup.search()` ranked over -- with two known gaps: the per-arm event
minima (d0.min/d1.min) and max_subgroups_search are NOT replayed here, so this
family is a superset of the one the identifier actually chose among, which
inflates the selection-bias term." The superset caveat is orthogonal to the
question here — superset or not, the family is enumerated once — but it is
quoted in full rather than trimmed to the convenient clause.

**Each candidate is fit once, into a fixed influence matrix.**
`.fs_mr_assemble()`, `R/fs_mr_inference.R:88–108`, loops over candidates once
and stores per-candidate quantities:

```r
    B[idx, g] <- pc$dfbeta
    bh[g] <- pc$beta_hat; sdv[g] <- pc$sigma_D; sz[g] <- length(idx)
```

Called once, at `R/fs_mr_inference.R:381`:

```r
  asm <- .fs_mr_assemble(df, candidates, spec)
```

**Every draw applies weights to those fixed quantities.**
`R/fs_mr_inference.R:445–447`:

```r
  Xi <- .fs_mr_multipliers(nrow(B), draws, multiplier)
  P  <- crossprod(B, Xi)                 # S x draws : D_g(b)
  beta_star <- bh + P
```

All `draws` perturbations are produced in a single matrix product against the
one fixed `B`. The draw loop, `R/fs_mr_inference.R:450–457`, only re-runs
admission and the selection rule over the same `S` columns:

```r
  for (b in seq_len(draws)) {
    bs <- beta_star[, b]
    pass <- .admit(bs)
    if (!length(pass)) next
    s <- .fs_mr_select(bs, .zcons(bs), sz, pass, reselection,
                       effect_neighborhood, selection_rule, log_scale)
```

The unrestricted admission branch indexes the same fixed family
(`R/fs_mr_inference.R:428`, `.admit <- function(bs) seq_along(bs)`), and the
complement block reuses the identical multiplier matrix `Xi`
(`R/fs_mr_inference.R:529`, `Pc <- crossprod(Bc, Xi)`).

**There is no re-enumeration anywhere in the loop.** Weights are applied to
fixed per-candidate quantities. By the task document's own criterion — "Reuse
of a fixed family = structural" — the incompatibility is **structural**, not
incidental. A model-based front end makes the enumerated family a function of
the observed data, so the single `B` that every draw reuses would no longer
correspond to the family the identifier could have selected from on a
perturbed sample, and the linearized selection event MR corrects for would not
be the selection event that occurred.

### 2.c The documented nuance at `family_status` — verbatim

**Location correction:** the task document cites
`R/forestsearch_helpers.R, ~:1988–1995`. On the local tree (and, since the file
is unmodified, on remote `5943469`) the passage is at
**`R/forestsearch_helpers.R:1939–1945`**. Line 1988 is inside an unrelated
`sg_focus` dispatch roxygen block. Verbatim:

```r
#'   \item{\code{"no-front-end"}}{\code{subgroup_method = "consistency"} with
#'     all three model-based front ends off, so no fitted model shapes the
#'     family on the observed data.  Deliberately NOT called "fixed": the
#'     manuscript's Section 2.1 fixed family additionally requires the cut
#'     locations to be resample-invariant, which quantile-derived cuts
#'     (\code{cut_type = "default"}) are not.  This level reports the weaker,
#'     checkable property.}
```

**Does the MR code depend on cut invariance, or only on family fixedness?**

As written, **only on family fixedness.** The MR path never derives, re-derives
or inspects a cut: `grep` for `cut_type`, `quantile` and `invarian` over
`R/fs_mr_inference.R`, `R/fs_mr_inference_methods.R` and
`R/consistency_resample.R` returns exactly one hit, and it is unrelated
(`R/fs_mr_inference.R:309`, the word "invariant" in a roxygen sentence about
the per-draw residual). The cut matrix `Z` is built once at
`R/forestsearch_main.R:2769` and the family once at `:3183–3197`; nothing
downstream depends on where the cuts came from. `.fs_family_status()`
(`R/forestsearch_helpers.R:1965–1977`) is purely descriptive and raises
nothing.

**What the source does not settle**, stated plainly: whether the *validity* of
the MR approximation — as a leading-order stand-in for the full bootstrap,
which is what `R/fs_mr_inference.R:1–29` says it is — requires cut invariance
in addition to family fixedness. Under a genuine bootstrap resample,
quantile-derived cuts would move, so the full-bootstrap selection event ranges
over a family MR's fixed `B` does not contain. The code does not need cut
invariance to *run*; whether the correction it computes is the right one
without it is a question the source does not answer. Folded into (d) below.

### 2.d Structural ⇒ follow-on recorded

Condition 2 is structural, so "MR + adaptive screening" is recorded as a
follow-on question under **Follow-on (not to solve)** below. Not solved, not
designed toward, no NOTE file opened.

---

## 3 — 0/50 on one replicate (blocker 3) — PARTIAL, halted by the Task 1 stop

**Execution status.** Step 0 (bundle inspection) and the reading of the harness
source are read-only and were completed before the Task 1 stop was known. Steps
1 and 2 require `devtools::install()` (Task 1 step 6), which is conditional on
a commit that the stop forbids, and the constraints state that a Task 1 stop
condition halts Tasks 2–3. **Steps 1 and 2 were therefore not executed, and
readouts (a), (b), (c), (e) and the executed half of (d) are ABSENT.** They are
not estimated, not inferred, and not reported as if measured.

The diagnostic script is written and ready to run:
`dev/glm-continuous-sims/verification/diagnose_md_harm_pilot_zero_detection.R`.
It reproduces replicate 1 under the harness's own pre-generated seed and prints
Step 0 plus readouts (a)–(e). It has **not been run**. It is untracked and, per
the Task 1 stop instruction to commit and push the report only, is **not
committed** by this session; it is available on request.

### 3.0 Step 0 — the object being diagnosed

Absolute path:

```
/home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/actg175/continuous/mr_sweep_md_harm/md_harm_s50_pilot/fs_md_harm_n1000_res.rds
```

Confirmed as the path named in the handoff. Size **1707 bytes**; mtime
**2026-08-07 19:55:47.259 -0700**. It is the only file in
`mr_sweep_md_harm/` (no grid, no HTML), matching the handoff.

| quantity | value | expectation | met |
|---|---|---|---|
| replicates | 50 | 50 | yes |
| `sum(detected)` | 0 | 0 | yes |
| `sum(mr_ok)` | 0 | — | — |

Per-replicate status tabulation, verbatim:

```
detected:          0
                  50

mr_ok:             0
                  50

betaHhat_status:  ok
                  50

nH_eval:           0
                  50

nHc_eval:       5000
                  50

label all NA  : TRUE
sg_def all NA : TRUE
```

`betaHhat_counts` attribute, verbatim:

```
$n_rules_total       [1] 0
$n_rules_resolved    [1] 0
$n_rules_unresolved  [1] 0
$n_reps_total        [1] 50
$n_reps_resolved     [1] 0
$n_reps_unresolved   [1] 0
$n_reps_undetected   [1] 50
```

`nH_eval + nHc_eval == 0 + 5000 == N` on all 50 rows — the ITT complement of
`ab53239` behaving as designed on undetected rows. `betaHhat_Hc` is a single
value across all 50 rows, `-30.991681593154`, as it must be when every row takes the
ITT complement on the same `df_super`. `n_true` (the true harm count per
replicate) is non-degenerate — min 304.0, median 346.0, mean 343.8, max 382.0
— so the DGM did produce a harm subgroup on every replicate; the failure is
downstream of the data.

Per-replicate wall clock inside `forestsearch()` (`t2_secs`):

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
0.00000 0.00000 0.00100 0.00056 0.00100 0.00100
```

**A median of 1 millisecond inside `forestsearch()` is the diagnosis in one
number.** No search, no fitting, and no MR ran on any replicate. Total cell
elapsed was `61.9 s` for 50 replicates, essentially all of it DGM setup.

Meta signature, verbatim:

```
  n_sample               : 1000
  n_sims                 : 50
  nb_boots               : NULL
  mr_draws               : 5000
  subgroup_method        : consistency
  sim_id_start           : 1
  sim_id_end             : 50
  seed_base              : 8316951
  parallel_mode          : sims
  n_workers              : 115
  elapsed_sec            : 61.9
  run_tag                : md_harm_s50_pilot
  direction              : harm
  focus                  : harm
  outcome_type           : continuous
  effect_measure         : MD
  cal_target_md          : -40
  adverse_outcome        : FALSE
  effect_threshold       : 30
  consistency_threshold  : 10
  threshold_orientation  : positive on the adverse_outcome=FALSE oriented scale (NOTE_threshold_sign_md.md)
  n_super                : 5000
  seed_scheme            : pre-generated table indexed by global sim_id
  pkg_version            : 0.2.0
  pkg_commit             : 79bce8d
  r_version              : 4.6.1
  built_at               : 2026-08-07 19:55:47
```

`truth` block, verbatim:

```
  marg_H         : -40.000000000000
  marg_Hc        : -26.255235876036
  effect_Q       : -40.000000000000
  effect_Qc      : -26.255235876036
  delta          : -26.255235876036
  beta_inter     : -13.744764123964
  prevalence_Q   :   0.344600000000
```

`marg_H == effect_Q` and `marg_Hc == effect_Qc` exactly, so the harness's
closed-form check B (θ† reproduces the DGM's own effects) held on this run.
`marg_H` is exactly `-40`, the calibration target.

The meta signature is internally consistent with the closed decisions:
`cal_target_md = -40`, `focus = "harm"`, `adverse_outcome = FALSE`,
`effect_threshold = +30` positive on the oriented scale per
`NOTE_threshold_sign_md.md`. **The threshold and orientation the bundle records
are correct.** The recorded `pkg_commit` is `79bce8d`, one commit behind the
current tip `5943469` (`5943469` is docs-only, so no `R/` difference).

### 3.1 Source-level cause, from the harness (read-only; not execution)

`quarto/simulations/actg175/continuous/mr_coverage_sweep_md_harm.qmd:125`,
verbatim:

```r
use_lasso <- TRUE; use_grf <- TRUE; use_twostage <- TRUE; is_rct <- TRUE
```

and those flags reach the search at `:282–287`, verbatim:

```r
    vi.grf.min = vi_grf_min, use_lasso = use_lasso, use_grf = use_grf,
    use_twostage = use_twostage, is.RCT = is_rct,
    adverse_outcome = adverse_outcome,
    details = FALSE, quiet = TRUE, seedit = sd_i,
    parallel_args = inner_parallel,
    mr_inference = TRUE,
```

`subgroup_method` is `"consistency"` (`meta$subgroup_method`, and `:81`,
`methods_run <- "consistency"`). That is **exactly** Condition 2 of the guard
traced in §2: `method == "consistency"` with `use_lasso` and `use_grf` both
TRUE and `mr_inference = TRUE`. The gate at `R/forestsearch_main.R:1292–1295`
fires **before any fitting**, and produces the two-front-end message quoted
verbatim in §2.a.

The harness then swallows it. `:275` and `:290–292`, verbatim:

```r
  fs.est <- tryCatch(suppressWarnings(forestsearch(
...
    error = function(e) NULL)
  rec$t2_secs <- proc.time()[3] - t0
  if (is.null(fs.est)) return(rec)
```

`stop()` → `tryCatch` returns `NULL` → `record_replicate()` returns the
`.na_record()` skeleton, whose `detected` field is initialised to `0L` at
`:237`. Every field downstream stays `NA`.

This accounts for every observation in §3.0 with nothing left over: 0/50
detected, 0/50 `mr_ok`, all `label`/`sg_def` `NA`, `n_reps_undetected = 50`
with `n_rules_total = 0`, and a median `t2_secs` of 1 ms — the cost of argument
capture and a `stop()`, with no search.

**This is not a threshold, cut-grid, or orientation failure.** The bundle's own
meta records the correct threshold (`+30`) and the correct orientation
(`adverse_outcome = FALSE`), and neither was ever consulted, because the guard
fires at `R/forestsearch_main.R:1292` — long before the resolution lines
`:1228`, `:1230`, `:1642`, `:1734`, `:1739` that `NOTE_threshold_sign_md.md`
traced, and before `get_fsdata()`/`get_FSdata()` builds any cut grid.

The finding is precisely what the handoff's supersession note anticipated: the
harness still carries the −61.87 sweep's `use_lasso`/`use_grf` front ends, and
those are incompatible with `mr_inference = TRUE` — structurally, per §2.b.

### 3.2 Readouts (a)–(e) — status

| readout | status |
|---|---|
| (a) candidate family size; `age`/`preanti` cut grids; the harness's actual `cut_type`/`conf.cont` arguments | **ABSENT** — not executed. Source-level partial: the qmd passes neither `cut_type` nor `conf.cont`, so both take their `forestsearch()` formal defaults; `confounders.name` is the 12 variables at `:127–128`; `maxk = 2`, `n.min = 60`, `d0.min = d1.min = 12` at `:123`. No grid was measured. |
| (b) bracketing of `age > 34` and `preanti <= 744.5` | **ABSENT** — not executed. No cut grid exists on the failing path, since the guard fires before `get_FSdata()`. |
| (c) truth-adjacent candidates and their oriented effects vs the threshold | **ABSENT** — not executed. No candidate was ever enumerated or fit on the failing path. |
| (d) threshold value, sign, orientation on THIS executing path | **PARTIAL.** From the bundle meta (recorded, not executed): `effect_threshold = 30`, `consistency_threshold = 10`, `adverse_outcome = FALSE`, orientation "positive on the adverse_outcome=FALSE oriented scale". From source: `effect.threshold = md_threshold = 30` is supplied and `hr.threshold` is not, so `user_set_threshold` is `TRUE` at `R/forestsearch_main.R:1228`. **The executed resolution was not measured** — `user_set_threshold` and the resolved `effect_threshold` reaching the engines were not observed at runtime, because the guard at `:1292` aborts before `:1642`. Whether the MR-compatible call resolves through the same lines as the sweep path is therefore **not established by execution**. |
| (e) failure funnel with the first gate that empties the qualifying set | **The funnel is empty at gate 0.** No candidate is ever enumerated, so the search's own `filter_counts` funnel (`R/subgroup_search.R:276–299`) never runs. The first gate that empties the qualifying set is not a search gate at all: it is the **MR configuration guard**, `R/forestsearch_main.R:1292–1295` → `.validate_mr_configuration()` Condition 2, `R/forestsearch_helpers.R:1912–1927`, `stop()` at `:1873–1876`. Everything downstream — enumeration, threshold floor, consistency screen, size gates — is unreached. |

### 3.3 Step 3 — the minimal config fix, STATED, NOT APPLIED

**Provisional**, and flagged as such: it rests on Step 0 plus the source
reading, not on the Step 1–2 execution the stop prevented. Readouts (a)–(c)
and (e) are still owed and could still surface a second, independent problem
behind this one.

**The single smallest change** is in the harness only, at
`quarto/simulations/actg175/continuous/mr_coverage_sweep_md_harm.qmd:125`: set
`use_lasso <- FALSE` and `use_grf <- FALSE`. It follows from the §3.1 readout —
the guard's Condition-2 message names exactly these two flags, and the 1 ms
median `t2_secs` in §3.0 confirms nothing past the guard executed. This is
precisely the MR-compatible configuration the handoff's supersession note
already adopts in principle ("no lasso/GRF screening, fixed candidate family"),
so it enacts a decision already taken rather than introducing a new one.

Deliberately left untouched: `cal_target_md = -40` stands; `focus = "harm"`
stands; no recalibration; `effect_threshold = +30` and
`consistency_threshold = 10` stand; `adverse_outcome = FALSE` stands;
`pconsistency = 0.90`, `fs.splits = 400`, `maxk = 2`, `n.min = 60`,
`d0.min = d1.min = 12` stand; the seed scheme, `mr_draws = 5000` and the
`mr_inference_args` stand; `use_twostage` and `is.RCT` stand. No `R/` change.
No DGM change.

**A second change is likely needed, and it is ranked second, not bundled.** The
harness's `tryCatch(..., error = function(e) NULL)` at `:275`/`:290` inside
`suppressWarnings()` converted a hard `stop()` into a silent `detected = 0`
across all 50 replicates. That is why a configuration error read as a
scientific null result. Ranked second because it does not affect what the next
pilot computes, only whether a future failure is visible; and because
diagnosing it properly needs the (a)–(c)/(e) readouts first, to confirm nothing
else is also being swallowed. Not proposed in detail here, and not applied.

**Nothing was applied. No re-pilot was started. STOP A did not begin.**

---

## Follow-on (not to solve)

- **MR + adaptive screening.** Condition 2 of `.validate_mr_configuration()` is
  structural (§2.b): the multiplier draws reuse one fixed influence matrix `B`
  built from one enumeration, so a data-estimated candidate family has no
  representation in the correction. Making MR valid under lasso/GRF screening
  is a different problem — it needs the family itself to be resampled, not just
  the influence contributions weighted. Recorded, not solved, not designed
  toward. No NOTE file opened.
- **Whether MR's validity needs cut invariance, not merely family fixedness**
  (§2.c). The `family_status` roxygen at `R/forestsearch_helpers.R:1939–1945`
  is explicit that `"no-front-end"` is deliberately weaker than the
  manuscript's §2.1 fixed family, because quantile-derived cuts are not
  resample-invariant. The MR code as written does not consult cuts at all, so
  it does not depend on invariance to run; whether the correction it computes
  is the right one without invariance is not settled by the source. Recorded.
- **`devtools::check()` vs `rcmdcheck::rcmdcheck()` are not interchangeable**
  for `--as-cran` (§1.3, observation 2). The task and handoff documents treat
  them as equivalent. Whichever the project standardises on should be named
  explicitly, so a future "clean check" claim means one thing. Recorded, not
  acted on.
- **Two tracked `man/` pages carry non-LaTeX-safe Unicode** —
  `man/fpr_calibration.Rd` and `man/sg_tables.Rd`, characters `─` (U+2500) and
  `ᶜ` (U+1D9C), last touched at `8a7096c`. Flagged as a side issue only; the
  CRAN-hygiene convention in `CLAUDE.md` is ASCII-only in R source. Not fixed,
  as it is outside this task's scope and the closeout is stopped.

---

## What this session did and did not do

Did: Task 0 preconditions; Task 1 tests; Task 1 check (stopped on
disagreement); Task 2 in full; Task 3 Step 0 and the harness source reading.

Did not: commit anything from the closeout; modify any `R/` file; modify any
test; modify the harness; run `devtools::install()`; execute Task 3 Steps 1–2;
apply any configuration fix; start blocker 4 / STOP A.

Awaiting review of the §1.2 check disagreement against this fetched file.
Committing this report is transport, not approval.
