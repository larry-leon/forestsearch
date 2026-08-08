# ADDENDUM to REPORT_blockers_1-3.md — closeout authorization and execution

Appends to `dev/glm-continuous-sims/REPORT_blockers_1-3.md`, committed at
`1697b5a` and reviewed from the remote. Nothing in the base report is revised;
this addendum records the authorization received, what was executed under it,
and where execution stopped.

---

## A. The authorizing message, recorded verbatim

> Closeout authorized after review of REPORT_blockers_1-3.md at 1697b5a,
> fetched from the remote — from Larry. Record this message verbatim in a
> report addendum.
>
> 1. Amended Task 1 acceptance, authorized: the closeout ONE-COMMIT proceeds
> now. Criterion: 0 errors, test stage OK inside the check, and every
> warning/note classified pre-existing or environmental with file-level
> evidence (report §1.2–1.3). Commit exactly R/betaHhat_truth.R,
> tests/testthat/test-betaHhat-contract.R, NAMESPACE, man/fs_build_eval_frame.Rd,
> man/fs_betaHhat_theta_dagger_check.Rd. The body states the amended criterion
> and cites §1.2–1.3: the 1 warning / 2 notes trace to man/fpr_calibration.Rd
> and man/sg_tables.Rd (U+1D9C, U+2500; pre-existing at 8a7096c), a missing
> tidy binary, and the manual-build byproduct — none attributable to the
> committed files, whose new man pages are pure ASCII. Push. Do NOT touch the
> two offending man pages tonight — the ᶜ is the Ĥᶜ notation and its fix is a
> queued morning decision.
>
> 2. Optional but preferred, in parallel or after: ONE
> devtools::check(args = "--as-cran") run on the same tree, appended verbatim
> to the addendum — this tests the hypothesis that check5.log's 0/0/0 line came
> from the --no-manual surface. One run, nothing modified to influence it,
> result recorded either way.
>
> 3. Then devtools::install() and proceed per the standing amendments: confirm
> the guard error by execution under the harness's own configuration on
> replicate 1 (quote the swallowed condition verbatim and the tryCatch at
> qmd:275/:290–292); then the flags-off replicate-1 run — use_lasso = FALSE,
> use_grf = FALSE, use_dina = FALSE, all else byte-identical — capturing
> readouts (a)–(e) and the status-code verification at
> R/subgroup_search.R:626–630; then overnight O1 under its amended
> MR-compatible configuration and preconditions, then O2, O3. The provenance
> commit (harness qmd as-is pre-fix; bundle if under ~5 MB, else
> size/mtime/sha256) still applies.
>
> 4. Unchanged gates: no harness edits, no fix application, no re-pilot, no
> STOP A, no DGM changes. Finish by committing and pushing the report addendum
> and, if reached, REPORT_overnight_funnel50.md — transport, not approval —
> then stop for the night.

An earlier message in the same exchange, also recorded verbatim, authorized the
diagnostic-script commit:

> On the script question — from Larry. Your deviation was correct: the Task 1
> stop condition's "report only" governs over the general wrap-up line, and
> resolving the conflict toward the harder stop was right.
>
> 1. Authorized: commit
> dev/glm-continuous-sims/verification/diagnose_md_harm_pilot_zero_detection.R
> now, as written — its own small commit, body stating it is committed
> as-designed and NOT yet executed, so any runtime corrections land as diffs
> against this version. Push.
>
> 2. Then proceed per the separate message beginning "Closeout authorized after
> review of REPORT_blockers_1-3.md at 1697b5a" — the amended Task 1 acceptance
> and closeout one-commit, the optional single devtools::check(args =
> "--as-cran") confirmatory run, devtools::install(), the guard-error
> confirmation under the harness config, the flags-off replicate-1 run with
> readouts (a)–(e) and the status-code verification, then amended O1, O2, O3,
> under the unchanged gates. If that message has not reached you, stop after
> the script commit and push, and wait for it.

---

## B. Commits made under this authorization

| commit | contents |
|---|---|
| `b86664b5` | `dev/glm-continuous-sims/verification/diagnose_md_harm_pilot_zero_detection.R`, committed as designed and **not yet executed**, per item 1 of the script message. |
| `ef504707` | The closeout ONE-COMMIT: exactly `R/betaHhat_truth.R`, `tests/testthat/test-betaHhat-contract.R`, `NAMESPACE`, `man/fs_build_eval_frame.Rd`, `man/fs_betaHhat_theta_dagger_check.Rd`. Nothing else. |

`man/fpr_calibration.Rd` and `man/sg_tables.Rd` were **not touched**, per the
explicit instruction; `git status --porcelain` on both is empty. The `ᶜ` is
noted as the Ĥᶜ notation, its disposition a queued morning decision.

### B.1 A defect in the script, committed rather than silently fixed

Before committing `b86664b5` a parse check (`Rscript -e 'parse(<file>)'`)
failed:

```
dev/glm-continuous-sims/verification/diagnose_md_harm_pilot_zero_detection.R:345:9: unexpected 'else'
344: .tbl <- if (is.list(fg) && is.list(fg$out.found)) fg$out.found$hr.subgroups
345:         else
             ^
```

A top-level `if`/`else` split across lines, which R cannot parse. The
instruction was to commit the script **as written** so that runtime corrections
land as diffs. The defect was therefore left in place and named in the commit
body rather than folded in invisibly. It was corrected as a diff before
execution, together with the two additions item 3 requires (`use_dina = FALSE`
on the flags-off run; the status-code verification at
`R/subgroup_search.R:626-630`).

---

## C. Item 2 — the single `devtools::check(args = "--as-cran")` run

One run, on the same tree, nothing modified to influence it. **The hypothesis
is confirmed: `check5.log`'s `0/0/0` line came from the `--no-manual`
surface.**

The invocation line, verbatim — `devtools::check()` injected `--no-manual`
exactly as `formals(devtools::check)$manual == FALSE` predicts:

```
* using options ‘--no-manual --as-cran’
```

Final status, verbatim:

```
── R CMD check results ───────────────────────────────── forestsearch 0.2.0 ────
Duration: 15m 30.3s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
```

Test stage, verbatim:

```
* checking tests ...
  Running ‘testthat.R’ [657s/487s]
 [657s/487s] OK
```

And the note that fired under `rcmdcheck` does not fire here, because the
manual build that creates its byproduct never runs:

```
* checking for non-standard things in the check directory ... OK
```

### The three surfaces side by side

| | `rcmdcheck::rcmdcheck(args="--as-cran")` | `devtools::check(args="--as-cran")` | `check5.log` (handoff) |
|---|---|---|---|
| options | `--as-cran` | `--no-manual --as-cran` | not recorded |
| result | `0 errors ✔ \| 1 warning ✖ \| 2 notes ✖` | `0 errors ✔ \| 0 warnings ✔ \| 0 notes ✔` | `0 errors ✔ \| 0 warnings ✔ \| 0 notes ✔` |
| duration | 11m 15.6s | 15m 30.3s | 14m 46.8s |
| tests | `[237s/187s] OK` | `[657s/487s] OK` | `[654s/484s] OK` |

`devtools::check()` reproduces the certified line exactly. Its test-stage
timing, `[657s/487s]`, sits within 0.5% of `check5.log`'s `[654s/484s]` —
consistent with the same surface on the same machine, though timing alone is
corroboration, not proof. The decisive evidence is the `--no-manual` option
line together with the identical `0/0/0` result.

**Conclusion.** The base report's §1.2 disagreement was a difference in check
surface, not a regression in the package or in the closeout files. Under the
amended acceptance criterion the closeout was correct to proceed. The
underlying condition — two tracked Rd pages carrying LaTeX-unsafe Unicode —
is real and remains unfixed by design (the `ᶜ` is the Ĥᶜ notation; disposition
is the queued morning decision). It is invisible to `devtools::check()` and
visible to a bare `R CMD check --as-cran`, which is the surface CRAN uses.

---

## D. Item 3 — execution

### D.1 Guard error confirmed by execution

`devtools::install()` completed (`* DONE (forestsearch)`) against the committed
state `ef504707`, so everything below runs against the installed package, as the
harness does.

Replicate 1 under the harness's configuration **verbatim** — `use_lasso = TRUE`,
`use_grf = TRUE`, `mr_inference = TRUE`, `subgroup_method` resolving to
`"consistency"` — with the harness's own pre-generated seed
`seed_for(1) = 1530735852`:

```
elapsed (s)                  : 0.0020
RESULT                       : ERROR (forestsearch() stopped)
--- verbatim error message ---
forestsearch(mr_inference = TRUE) requires an identifier aligned with the reported effect.
  The model-based front ends use_lasso and use_grf are set TRUE, so the candidate family is estimated from the same data the search runs on and is not fixed.
  Set use_lasso = FALSE, use_grf = FALSE and re-run.
--- end error message ---
```

**2 milliseconds.** The `stop()` fires at `R/forestsearch_main.R:1292–1295`
before any fitting, which is why the pilot's 50 replicates each cost ~1 ms
inside `forestsearch()` (base report §3.0).

The swallowing site, verbatim from the harness:

```r
# quarto/simulations/actg175/continuous/mr_coverage_sweep_md_harm.qmd:275
  fs.est <- tryCatch(suppressWarnings(forestsearch(

# :290-292
    error = function(e) NULL)
  rec$t2_secs <- proc.time()[3] - t0
  if (is.null(fs.est)) return(rec)
```

`stop()` → `tryCatch` yields `NULL` → early return at `:292` → the
`.na_record()` skeleton, whose `detected` is initialised `0L` at `:237`. The
configuration error is recorded as a scientific null result. Blocker 3's
proximate cause is **confirmed by execution**, not merely by source reading.

### D.2 Reproduction fidelity — a finding that affects the re-pilot

The seed lookup reproduces exactly, and so does the super-population:

```
seed_for(1)                  : 1530735852
bundle results$seed[1]       : 1530735852
seed reproduction matches    : TRUE
df_super n                   : 5000
P(Q) on df_super             : 0.344600      (== meta prevalence_Q 0.3446)
```

The **realized replicate does not**:

```
replicate 1 n_true (flag_harm): 334
bundle results$n_true[1]      : 341
n_true reproduction matches   : FALSE
```

**Cause, confirmed by direct test.** `simulate_from_glm_dgm()` calls
`set.seed(seed)` at `R/simulate_from_glm_dgm.R:78` with `kind = NULL`, which
per `?set.seed` *keeps the current RNG kind* rather than resetting it. The
harness runs replicates through `%dofuture%` with
`.options.future = list(seed = TRUE)` (qmd:348), and future sets workers to
L'Ecuyer-CMRG. A plain `Rscript` runs under Mersenne-Twister. The same integer
seed therefore names two different streams:

```
seed_for(1): 1530735852
bundle n_true[1] (target to reproduce): 341

RNGkind = Mersenne-Twister    n_true = 334
RNGkind = L'Ecuyer-CMRG       n_true = 341
```

L'Ecuyer-CMRG reproduces the bundle exactly. The hypothesis is **confirmed**,
not merely plausible.

**Why this matters beyond bookkeeping.** The harness's stated requirement
(qmd:87–89) is that "results must be invariant to worker count and to batch
position, so the seed for replicate i is looked up, never derived from loop
position." The lookup does deliver that. But the realized *data* additionally
depends on the parallel backend being active and on which RNG kind it installs
— so a replicate rerun sequentially for debugging, or under a different
`future` configuration, is not the same replicate. Recorded as a finding for
the re-pilot; **not fixed, and no harness edit made.**

**Consequence for this session's readouts.** The first flags-off run was
executed under Mersenne-Twister and therefore describes the 334-case draw, not
the pilot's replicate 1. It was re-run under L'Ecuyer-CMRG for fidelity; §D.3
reports the faithful run, and flags any place the two differ.

### D.3 Flags-off replicate-1 run: readouts (a)–(e)

Configuration: `use_lasso = FALSE`, `use_grf = FALSE`, `use_dina = FALSE`,
everything else byte-identical to the harness, seed `1530735852`, executed under
`RNGkind("L'Ecuyer-CMRG")` so the draw matches the pilot's worker context. The
draw is confirmed faithful:

```
n_true (flag_harm) : 341   (bundle results$n_true[1]: 341)
```

Returned: `family_status = "no-front-end"`, `sg.harm = !{gender} & !{race}`,
`|sg.harm.id == 1| = 85`.

**(a) Candidate family and cut grids.** The harness supplies neither `cut_type`
nor `conf.cont`, so both take their `forestsearch()` formal defaults
(`cut_type = "default"`). 27 cut labels over the 12 confounders:

```
 [1] "age <= 35"      "age <= 34"      "age <= 30"      "age <= 40"
 [5] "preanti <= 391" "preanti <= 139" "preanti <= 0"   "preanti <= 777"
 [9] "wtkg <= 74"     "wtkg <= 66"     "wtkg <= 82"     "karnof <= 95"
[13] "karnof <= 90"   "cd40 <= 348"    "cd40 <= 334"    "cd40 <= 257"
[17] "cd40 <= 421"    "cd80 <= 986"    "cd80 <= 903"    "cd80 <= 648"
[21] "cd80 <= 1184"   "hemo"           "homo"           "drugs"
[25] "race"           "gender"         "symptom"
```

`dummy()` expands each into both directions, giving 52 indicator columns, and
`maxk = 2` yields **1485 combinations evaluated**.

- age grid: `age <= 35 | age <= 34 | age <= 30 | age <= 40`
- preanti grid: `preanti <= 391 | preanti <= 139 | preanti <= 0 | preanti <= 777`

**(b) Bracketing.** Both true boundaries are bracketed.

| boundary | nearest at/below | nearest above | brackets |
|---|---|---|---|
| `age > 34` | `34` (exact) | `35` | yes |
| `preanti <= 744.5` | `391` | `777` | yes |

`age <= 34` is on the grid exactly, so `age > 34` is exactly representable as
its negation. `preanti <= 744.5` is not directly representable; the nearest
expression is `preanti <= 777` (355 subjects vs the true 341).

**(c) Truth-adjacent candidates and their oriented effects — the substantive
readout.** Oriented MD, i.e. on the `adverse_outcome = FALSE` negated scale
where positive means harm, computed directly from the replicate:

```
region                                    n orientedMD       SE vs floor 30
TRUE Q (flag_harm == 1)                 341    25.1970  13.6786 REJECT (status 6)
!{age <= 34} & {preanti <= 777}         355    23.7416  13.4793 REJECT (status 6)
!{age <= 34}  (age > 34) alone          488    17.6522  11.7319 REJECT (status 6)
{preanti <= 777} alone                  752    31.5521   9.4987 PASS
!{age <= 35} & {preanti <= 777}         326    26.2101  14.0353 REJECT (status 6)
!{age <= 30} & {preanti <= 777}         498    31.2177  11.4512 PASS
!{age <= 34} & {preanti <= 391}         268    30.6495  15.4415 PASS
SELECTED: !{gender} & !{race}            85   105.7582  25.0732 PASS
```

**The true region Q itself scores 25.20 against a floor of 30, and is rejected
at status 6.** Its population oriented value is 40 (the `cal_target_md = -40`
calibration), so this replicate's estimate sits about 1.08 SE below the truth —
ordinary sampling variation, not a defect. The grid analogue
`!{age <= 34} & {preanti <= 777}` scores 23.74 and is rejected with it.
Confirming this is not a grid artefact: the truth-adjacent conjunction is
absent from the 609 qualifying candidates.

Meanwhile the selected winner is an 85-subject region scoring **105.76**, and
609 candidates clear the floor with oriented effects ranging **30.03 to
105.76**. That is textbook selection optimism over 1485 candidates at
n = 1000 — and precisely the quantity MR exists to correct, which is why §D.4
(MR never running) is the more consequential of the two configuration defects.

**A forward-looking caveat, from one replicate and labelled as such.** The
true region's oriented MD carries SE ≈ 13.7 at n = 1000 while the floor sits at
30 and the truth at 40. The floor is therefore only ~0.73 SE below the true
effect, so rejecting the truth is an ordinary outcome rather than a rare one,
and the search will preferentially select small high-variance regions instead.
This bears on what the re-pilot's identification readouts (sensitivity, PPV)
will show. It is **not** the cause of the 0/50 — that is entirely the guard —
and nothing here proposes changing the threshold, `cal_target_md`, or `n`. One
replicate cannot establish a rate; recorded so the re-pilot's numbers are not a
surprise.

**(d) Threshold, sign, orientation on this executing path.** Verbatim from
`threshold_config` on the returned object:

```
  outcome_type             : continuous
  effect_measure           : MD
  screening                : 30
  consistency              : 10
  screening_natural        : 30
  consistency_natural      : 10
  scale                    : identity
  pconsistency             : 0.9
  screening_description    : MD >= 30.0000
  consistency_description  : MD >= 10.0000 per split
  adverse_outcome executed : FALSE
  admission$effect_floor   : 30
  sg_focus                 : hr
```

`effect.threshold = 30` is supplied and `hr.threshold` is not, so
`user_set_threshold` is TRUE at `R/forestsearch_main.R:1228` and `hr.threshold`
takes 30 at `:1230`. `scale: identity` confirms the MD path takes no log at
`:1734` and is caught by `is_identity` at `:1739` — **the MR-compatible call
resolves through the same lines the sign NOTE traced for the sweep path.** The
value reaching the engines is `+30`, sign included, on the oriented scale where
positive means harm. Blocker 3's orientation question is answered: the
orientation is correct and was never the problem.

**(e) The failure funnel.** Nothing empties the qualifying set:

```
  n_evaluated            :   1485
  n_passed_variance      :   1428   (-57)
  n_passed_prevalence    :   1428   (-0)
  n_passed_redundancy    :   1394   (-34)
  n_passed_events        :   1394   (-0)
  n_passed_sample_size   :   1293   (-101)
  n_passed_cox           :   1293   (-0)
  n_passed_hr            :    609   (-684)
```

**Status-code verification**, `R/subgroup_search.R:626–630`, quoted from source:

```r
    # Status 6: Check effect threshold.  Skipped when the effect floor is
    # disabled (sg_focus = "maxeff" retains the full estimable family so the
    # argmax is unconditional -- see forestsearch() search_overrides).
    if (!disable_effect_floor && glm_result$hr <= hr.threshold) {
      return(list(status = 6L, result = NULL))
    }
```

The comparison is `<=`, so a candidate exactly at the floor is rejected;
`n_passed_hr` counts candidates strictly above 30. The floor is in force on
this path (`sg_focus = "hr"`, `admission$effect_floor = 30`, not `"maxeff"`),
and it is the largest single gate, removing 684 of 1293. It is also the gate
that removed the true region.

**So on the flags-off path there is no funnel failure at all**: 609 candidates
qualify and a subgroup is selected. The pilot's 0/50 is not a funnel outcome —
it is the guard, upstream of the funnel, exactly as base report §3.2 recorded.

### D.4 A second, independent configuration defect: MR never runs

RUN B returned with **`mr_inference` NULL and `mr_harm_confirmed = NA`** even
though the guard had passed. This is independent of the guard and would survive
the §3.3 fix.

Eligibility, `R/forestsearch_main.R:3142–3146`:

```r
  .mr_glm_ok <- consistency_method == "resample" && !is.null(estimator_fn)
  .mr_cox_ok <- outcome_type == "survival" && is.null(estimator_fn)
  .mr_eligible <- isTRUE(mr_inference) && !is.null(sg.harm) &&
    (.mr_glm_ok || .mr_cox_ok) &&
    !is.null(grp.consistency) && !is.null(grp.consistency$sg.harm.id)
```

Measured on the run:

```
mr_inference on object      : FALSE
mr_harm_confirmed           : NA
consistency_method resolved : split
consistency_method FORMAL default (match.arg takes the first) : c("split", "resample")
```

The harness never passes `consistency_method`, so `match.arg()` resolves it to
`"split"`, `.mr_glm_ok` is FALSE, and the continuous/GLM path has no MR route.
The skip is announced by `.mr_skip()` (`R/forestsearch_main.R:3152–3168`) —
captured verbatim by re-running with `quiet = FALSE`:

```
mr_inference = TRUE but multiplier resampling was NOT performed: MR on a GLM outcome requires consistency_method = "resample"; this analysis used consistency_method = "split". Re-run with consistency_method = "resample" to obtain MR.
```

The harness passes `quiet = TRUE` (qmd:285), so this message never reached the
pilot log. The announcement machinery worked exactly as designed; the harness
silenced it.

**Consequence.** Fixing only `use_lasso`/`use_grf` yields a pilot that detects
subgroups on every replicate and still reports `mr_ok = 0` on every replicate —
a differently-shaped null result, with the MR columns (`t2_H_est`, `t2_H_lo`,
`t2_H_hi`, `t2_H_se_ij`, `nv_H_est`, the complement block) all `NA`, since
qmd:317 gates them on `!is.null(g)`. Coverage and unbiasedness — the harness's
entire purpose — would be uncomputable.

**Second fix, STATED, NOT APPLIED.** Add `consistency_method = "resample"` to
the `forestsearch()` call in the harness. It follows from the verbatim skip
message above, which names the required value. It is not a substitute for the
first fix: the guard fires earlier, so both are needed, and neither was applied.

The two fixes are independent and both necessary:

| # | change | evidence | gate it clears |
|---|---|---|---|
| 1 | `use_lasso <- FALSE; use_grf <- FALSE` (qmd:125) | §D.1 guard error, verbatim | `.validate_mr_configuration()` Condition 2, `R/forestsearch_main.R:1292` |
| 2 | add `consistency_method = "resample"` to the call | §D.4 skip message, verbatim | `.mr_glm_ok`, `R/forestsearch_main.R:3142` |

Ranked third, unchanged from base report §3.3: the harness's
`tryCatch(error = function(e) NULL)` and `quiet = TRUE` between them converted
two hard, well-worded diagnostics into a silent `detected = 0` / `mr_ok = 0`.
Both messages existed and both were suppressed.

---

## E. Provenance

Committed as-is, **pre-fix**, per the standing provenance requirement. Both
artifacts are far under the ~5 MB threshold, so the files themselves are
committed rather than only their digests.

| artifact | size | mtime | sha256 |
|---|---|---|---|
| `quarto/simulations/actg175/continuous/mr_coverage_sweep_md_harm.qmd` | 28,563 B | 2026-08-07 19:54:19.590 -0700 | `1083c0b8d9292335c5c507b80caa40b600ea6a920f4cda70444c4ffd65df1098` |
| `.../mr_sweep_md_harm/md_harm_s50_pilot/fs_md_harm_n1000_res.rds` | 1,707 B | 2026-08-07 19:55:47.259 -0700 | `6b9834869f60b5e9e6b63a958e65764d295094469f3f36ac7252428974adc6d4` |

The qmd is committed with `use_lasso <- TRUE; use_grf <- TRUE` at `:125` and no
`consistency_method` argument — i.e. exactly the configuration that produced
the 0/50 bundle, **unmodified**. Neither stated fix is applied.

---

## F. Where execution stops, and why

### F.1 Completed under this authorization

- Item 1 — closeout ONE-COMMIT (`ef504707`), pushed. Exactly the five files.
- Item 2 — the single `devtools::check(args = "--as-cran")` run (§C).
  Hypothesis confirmed.
- Item 3, first half — `devtools::install()`; guard error confirmed by
  execution (§D.1); flags-off replicate-1 run with readouts (a)–(e) and the
  status-code verification at `R/subgroup_search.R:626–630` (§D.3).
- Provenance commit (§E).

### F.2 Not started: O1, O2, O3

**These are not defined in any document reachable from this repository.** The
authorizing message refers to "overnight O1 under its amended MR-compatible
configuration and preconditions, then O2, O3", to "the standing amendments",
and to `REPORT_overnight_funnel50.md`. Searches performed:

```
$ grep -rln "overnight\|funnel50\|standing amendment" dev/ quarto/
  (no match in any .md; hits only in quarto/.quarto/ index JSON and unrelated
   legacy shell/qmd files)

$ grep -rn "O1\b.*O2\b\|O1:" dev/*.md dev/*/*.md
  (no match)

$ ls dev/glm-continuous-sims/
  ADDENDUM_glm_pathway_context.md  CC_TASK_md_mr_harness.md  design-checks/
  HANDOFF_glm_continuous_simulations.md  HANDOFF_md_mr_harness_session.md
  NOTE_survival_betaHhat_na.md  NOTE_target_is_collapsibility.md
  NOTE_threshold_sign_md.md  README.md  REPORT_blockers_1-3.md
  SPEC_betaHhat_md.md  verification/
```

Neither `HANDOFF_md_mr_harness_session.md` (v2) nor
`CC_TASK_md_mr_harness.md` defines an O1/O2/O3 scheme; the handoff's numbered
items are blockers 1–4 and its phases are STOP A / STOP B.

Guessing at the content of three overnight runs — their configurations,
preconditions, replicate counts, and acceptance criteria — would be exactly the
kind of fabrication the handoff's "incidents this handoff exists to prevent
repeating" section exists to stop. **They are therefore not started.** Send the
definitions and they can be picked up directly; everything they depend on
(installed package at `ef504707`, confirmed diagnosis, provenance committed) is
in place.

`REPORT_overnight_funnel50.md` is consequently **not** written, per the "if
reached" qualifier in the authorizing message.

### F.3 Gates honoured

No harness edit. No fix applied — both stated fixes (§D.4 table) remain
proposals. No re-pilot. No STOP A. No DGM change. No `R/` change. The two
offending man pages untouched.

### F.4 What a reviewer should decide next

1. Whether both stated fixes go in together, given that fix 1 alone yields a
   pilot with `detected` high and `mr_ok = 0` throughout (§D.4).
2. Whether the RNG-kind finding (§D.2) warrants any change before the
   re-pilot, given that replicate data currently depends on the parallel
   backend being active.
3. Whether the floor-vs-truth margin at n = 1000 (§D.3, one replicate:
   truth 40, floor 30, SE ≈ 13.7) is accepted as the phenomenon — consistent
   with the closed decision that the selection rate is "on the record, not
   tuned away" — or warrants comment in the re-pilot's readouts.
4. O1/O2/O3 definitions.

Committing this addendum is transport, not approval.
