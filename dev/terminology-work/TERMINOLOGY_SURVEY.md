# TERMINOLOGY SURVEY — "gate" and "tier"

**Status:** survey only. Nothing was renamed. No file outside `dev/terminology-work/`
was created or modified.

**Date:** 2026-07-29 · **Package:** forestsearch 0.2.0 · **R:** 4.6.1 ·
**Branch:** `feature/sg-focus-maxeffcons`

---

## 0. Guard verdict

| item | value |
|------|-------|
| hashing function | `fs_hash_sources()` from `dev/efficiency-eval/R/00_guard.R` (reused, unmodified) |
| files hashed | 133 (installed package tree + worktree `R/`) |
| verdict | **`ok = TRUE` — 133 source files verified unchanged** |
| `git status` | one new untracked tree, `?? dev/terminology-work/`. **No tracked file modified, deleted, or renamed.** |
| `sed` / bulk substitution | none used, not even as a dry run |

`git status --porcelain` also lists three pre-existing untracked entries under
`dev/sg-focus-work/` (`verify_batchsize_fix.R` ×2 and the brief itself). These
are not this survey's doing — the first two predate the session and the third is
the brief. Everything this survey wrote is under `dev/terminology-work/`.

Environment capture: `dev/terminology-work/out/env_capture.csv`.
Before/after digests: `out/hash_before.rds`, `out/hash_after.rds`.

**One process note.** The brief's launch line names
`dev/terminology-work/CC_BRIEF_terminology_survey.md`; at the time of the first
invocation the brief existed only at
`dev/sg-focus-work/CC_BRIEF_terminology_survey.md`. Both copies are now present
and are byte-identical (`md5 03893a218bf7fef8542c59debffefa6f`). The survey was
executed against that content.

---

## 1. Headline findings

Ordered by consequence. Items 1–4 are behaviour findings that surfaced while
answering §5A; they are reported because the brief asks for them, and they matter
to the rename only insofar as a name should not describe behaviour the code does
not have.

**F1 — `debias_gate = TRUE` is a silent no-op for every GLM outcome under the
default `consistency_method`.**
The consistency branch runs the gate only when
`.dg_glm_ok || .dg_cox_ok` (`R/forestsearch_main.R:2817-2820`):

```r
.dg_glm_ok <- consistency_method == "resample" && !is.null(estimator_fn)
.dg_cox_ok <- outcome_type == "survival"       &&  is.null(estimator_fn)
```

`consistency_method` defaults to `"split"`. For a GLM outcome (`estimator_fn`
non-`NULL`, `outcome_type != "survival"`) both tests are `FALSE`, so the block is
skipped — with **no warning**. Executed evidence
(`dev/terminology-work/R/12_spec_probe3.R`):

| outcome | `consistency_method` | subgroup found | gate ran | `harm_flag_debiased` | warned? |
|---------|---------------------|----------------|----------|----------------------|---------|
| binary | `split` (default) | TRUE | **FALSE** | `NA` | **FALSE** |
| binary | `resample` | TRUE | TRUE | TRUE | FALSE |
| count | `split` (default) | TRUE | **FALSE** | `NA` | **FALSE** |
| count | `resample` | TRUE | TRUE | TRUE | FALSE |
| continuous | `split` / `resample` | FALSE | FALSE | `NA` | FALSE |

This is exactly the failure mode §5A.2.5 asks about, in its worst form: the user
opted in, nothing ran, nothing said so, and `harm_flag_debiased` is `NA` — which
`isTRUE()` renders as "harm not confirmed". The distinction between *not
confirmed* and *not computed* is unrecoverable from the return value alone.
(The `continuous` rows are uninformative for the gate: this DGM produced no
subgroup, so the branch was unreachable either way.)

**F2 — the gate is recomputed inside every bootstrap replicate and every CV
fold, and the result is discarded.**
`debias_gate` is captured into `args_call_all` (via
`mget(args_names, ...)` at `R/forestsearch_main.R:1087` — it is *not* in any
`.sync_args_call_all()` list, it arrives with the bulk capture). Both
consumers re-run `forestsearch()` from that list and neither overrides
`debias_gate`:

* `R/bootstrap_analysis_dofuture.R:397` `args_FS_template <- fs.est$args_call_all`,
  mutated at `:549-578` (`df.analysis`, `details`, `plot.*`, `quiet`, `grf_*`,
  `dina_*`, `parallel_args` — not `debias_gate`), then `:585 do.call(forestsearch, args_FS_boot)`.
* `R/forestsearch_cross_validation.R:315` `fs_args <- fs.est$args_call_all`,
  then `:401 do.call(forestsearch, cv_args)`.

Measured by tracing `fs_debias_gate()` (`dev/terminology-work/R/12_spec_probe3.R`,
`13_spec_probe4.R`):

| driver | replicates/folds | `fs_debias_gate()` calls |
|--------|------------------|--------------------------|
| `forestsearch_bootstrap_dofuture(nb_boots = 2)` | 2 | **2** |
| `forestsearch_Kfold(Kfolds = 2)` | 2 | **2** |

At the documented default `draws = 2000L`, a 500-replicate full bootstrap
performs 500 extra multiplier-resampling runs whose output is never read. This
contradicts §5A.2.6's expected answer ("only once on the original analysis") and
is a pure cost. It is a one-line fix (`args_FS_boot$debias_gate <- FALSE`), but
it is a **behaviour** change and therefore out of this brief's scope.

**F3 — the re-selection family does not replay the identifier's admissibility,
for any focus.** Three independent mismatches, only one of which was previously
recorded:

| criterion | identifier applies | gate family applies | consequence |
|-----------|--------------------|---------------------|-------------|
| `maxk`, `n.min` | yes | yes (`main.R:2840-2843`) | — |
| `d0.min` / `d1.min` (per-arm event minima) | yes (`subgroup_search.R:591-601, 639`, `meets_event_criteria()`) | **no** | gate family is a superset; bias term explores candidates the search rejected |
| `max_subgroups_search` | yes (truncates the ranked pool) | **no** | largest mismatch in practice |
| consistency floor under `sg_focus = "maxeff"` | **no** (ungated by contract) | yes (via `t_g`) | the recorded mismatch, confirmed |

Measured on one 400-subject survival run (`sg_focus = "maxSG"`,
`max_subgroups_search = 5`):

```
gate n_family                        : 112
search filter_counts:
  n_evaluated 136 → n_passed_events 116 → n_passed_sample_size 106
              → n_passed_cox 106 → n_passed_hr 23
then max_subgroups_search = 5 truncated the ranked pool 22 → 5
```

So the gate re-selected over **112** candidates while the identifier chose among
**5**. The `112 > 106` gap is the missing event-minima filter; the `112` vs `5`
gap is the un-replayed `max_subgroups_search`. Both inflate `selection_bias`,
i.e. over-correct the de-biased estimate. This is a behaviour finding, not a
naming one, and is out of scope — but it means the phrase "faithful to the
Tier-1 search" (`main.R:2828`) is not currently accurate.

**F4 — the `tryCatch` does fire in practice.** Requested by §5A.2.5. On the GRF
branch with factor covariates the family reconstruction throws and the gate
returns `NULL` with a warning:

```
consistency / maxSG      found=TRUE  gate_ran=TRUE   reselection=maxSG
dina        / maxSG      found=TRUE  gate_ran=TRUE   reselection=maxSG
grf         / maxSG      found=TRUE  gate_ran=FALSE  →  "'>' not meaningful for
                                        factors | debias_gate (survival) failed"
```

`.fs_dg_family_from_table()` builds comparisons with `op_right = ">"`
(`main.R:2038`) and `.fs_dg_members_from_conj()` evaluates them against the raw
frame, where `stage`/`sex` are factors. Unlike F1 this *does* warn, but it still
lands on `harm_flag_debiased = NA`, so a script reading `isTRUE(...)` sees
"not confirmed".

**F5 — two roxygen statements are contradicted by the code.** The brief predicted
one; there are at least three.

| location | claim | actual |
|----------|-------|--------|
| `fs_debias_gate_methods.R:162` | `reselection_default` is `"maxcons"` for consistency | derived from `sg_focus` at `main.R:2877-2879`; observed `"maxSG"` for `sg_focus = "maxSG"`, `"maxeff"` for `maxeffCons` |
| `fs_debias_gate.R:265` | `include_complement` default `FALSE` | every internal caller passes `TRUE` (`methods.R:187`, `main.R:2884`); the complement is always computed in practice |
| `gate_estimates_table.R:28-31` | "de-biased CI uses the subgroup robust SE" | uses the IJ SE under the default `ci_method = "ij"`; robust SE only under `"wald"` |

**F6 — `fs_debias_gate()` performs MR only. It does not span FB.** No split is
needed; the §5.3 structural concern does not materialise. See §5A.5.

**F7 — no deprecation shim is required: the `debias_gate` family never shipped
on CRAN. CONFIRMED (2026-07-29) from two independent directions.**

The brief's premise is correct — **forestsearch 0.1.0 is on CRAN, published
2026-03-23.** What does not follow is the inference that a hard break to
`debias_gate` invalidates existing user scripts, because `debias_gate` was not
in that release.

*Evidence 1 — the CRAN 0.1.0 reference manual.* Fetched
`https://cran.r-project.org/web/packages/forestsearch/forestsearch.pdf` and
extracted with `pdftotext` (134 383 characters). Header confirms
`Package 'forestsearch' / Version 0.1.0`; positive control
`forestsearch_bootstrap_dofuture` appears 12 times, so the extraction is sound.
**The string `gate` occurs zero times in the entire manual** — not
`debias_gate`, not `gate_estimates_table`, not the word itself.

*Evidence 2 — the git record.* `DESCRIPTION` went `0.0.0.9000` → `0.1.0`
(`c4bd485`, 2026-02-13) → `0.2.0` (`9409d2c`, 2026-03-27), with **no `0.1.x`
patch release** in between. Across the whole 0.1.0 window, the number of files
containing each name is zero:

| name | files at 0.1.0 (`git grep -c … 9409d2c^`) | first introduced |
|------|------------------------------------------|------------------|
| `debias_gate` | 0 | `a63900f` **2026-06-08** |
| `fs_debias_gate` | 0 | 2026-06-08/09 |
| `gate_estimates_table` | 0 | `44385df` 2026-06-09 |
| `harm_flag_debiased` | 0 | `6af7e99` 2026-06-09 |
| `gate_draws` | 0 | 2026-06-08 |

The family first appears **2026-06-08, roughly eleven weeks after CRAN published
0.1.0** and ten weeks after the 0.2.0 bump. A March tarball cannot contain a
June identifier.

**Conclusion: `debias_gate`, `debias_gate_args`, `fs_debias_gate()`,
`gate_estimates_table()`, `harm_flag_debiased` and `out$debias_gate` are all
0.2.0-only and unreleased. No deprecation shim is required for any of them.**
The blast radius is real but is entirely the author's own `quarto/` notebooks,
not third-party scripts. This materially simplifies §6.

---

## 2. The vocabulary, corrected

The brief's four-way table is right about the de-biasing code. Surveying the
whole package turned up **two further, independent uses** of "gate" that a purge
would damage:

| # | term | what it is | genuinely a "gate"? | in scope? |
|---|------|-----------|--------------------|-----------|
| 1 | Admissibility criteria | `hr.threshold`, `pconsistency.threshold`, `n.min`, `d0.min`, `d1.min`, `maxk` | no | rename |
| 2 | Selection rule | `sg_focus` | no | rename |
| 3 | Post-selection inference | FB/MR de-biasing + interval estimation | **no** — the bulk of `debias_gate` | rename |
| 4 | Harm confirmation | the decision that a de-biased estimate still indicates harm | **yes** | **keep the word** |
| 5 | **Calibration tolerance gate** | hard pass/fail on whether a calibrator hit its target, then `stop()` | **yes** | **out of scope — do not touch** |
| 6 | **Generic control-flow guard** | informal prose: "the gate keys on `sg.harm`", "threshold-gated trimming" | loosely | leave; reword only opportunistically |

Sense 5 occupies 8 occurrences across `calibrate_helpers.R:15,79`,
`calibrate_glm_interaction.R:144`, `calibrate_k_treat.R:147`,
`find_k_inter_main.R:129`, `simulate_from_dgm.R:790`, `sim_aft_gbsg.R:1039`. It
is a genuine pass/fail applied to a result and has nothing to do with subgroup
identification. **A blanket `gate` → `mr` substitution would corrupt it.** This
is the concrete instance of the risk §3 warns about.

Sense 6 covers `forestsearch_cross_validation.R:75,144`,
`subgroup_consistency_main.R:641`, `subgroup_consistency_helpers.R:566`,
`simulation_tables.R:756,1183`, `guohe_algorithm3.R:177,272`,
`guohe_adaptive_r.R:43`, `dina_subgroup.R:18,969`.

Confirming the brief's own two hypotheses:

* **`debias_gate` names sense 4 while chiefly performing sense 3** — confirmed.
  Of its returned object, one field (`harm_flag`) is sense 4; the other 15 are
  sense 3.
* **`gate_draws` counts MR draws, not anything belonging to a gate** — confirmed,
  and it is not a package name at all (§5.1).

---

## 3. Occurrence inventory

### 3.1 Scope actually searched

All of `R/`, `man/`, `NAMESPACE`, `tests/`, `vignettes/`, `quarto/`, `NEWS.md`,
`README.md`, `DESCRIPTION`, `_pkgdown.yml`, `dev/`. Within `quarto/`, generated
artefacts were excluded (`.quarto/_freeze`, `*_files/`, `_site`, rendered
`.html`) — **681** source `.qmd`/`.md`/`.R` files remain. An unfiltered sweep
returns ~30 800 `gate` lines in `quarto/`; that number is rendered-output noise
and is not used anywhere below.

### 3.2 `gate` in `R/` — 208 occurrences after removing substring false positives

Excluded as unrelated substrings: `aggregate`, `delegate`, `surrogate`,
`propagate`, `investigate`, `negate`, `mitigate`, `gatekeeper`.

| file | n | dominant sense |
|------|---|----------------|
| `forestsearch_main.R` | 65 | 3 (+ arg surface) |
| `fs_fdr_report.R` | 38 | 4 (`c_gate`) |
| `fs_debias_gate_methods.R` | 31 | 3 |
| `fs_debias_gate.R` | 21 | 3, with 4 at the flag |
| `forestsearch_helpers.R` | 18 | 3 (comments) |
| `gate_estimates_table.R` | 13 | 3 |
| `grf_main.R`, `grf_subg_harm_glm.R` | 2 + 2 | 3 (comments) |
| `calibrate_helpers.R`, `calibrate_glm_interaction.R`, `calibrate_k_treat.R`, `find_k_inter_main.R`, `simulate_from_dgm.R`, `sim_aft_gbsg.R` | 8 | **5 — out of scope** |
| `forestsearch_cross_validation.R`, `subgroup_consistency_main.R`, `subgroup_consistency_helpers.R`, `simulation_tables.R`, `guohe_algorithm3.R`, `guohe_adaptive_r.R`, `dina_subgroup.R` | 12 | **6 — leave** |

Split by class: **118 Docs** (comment or roxygen) / **90 code**.
Raw listing: `out/inv_gate_R_filtered.txt`.

### 3.3 `tier` — the FB/MR rename surface is prose-only inside the package

There is **no `tier1`/`tier2`/`t1`/`t2` identifier anywhere in `R/` that refers
to FB or MR.** "Tier-1"/"Tier-2" appears in the package only in comments,
roxygen, and one user-visible string:

| class | count | locations |
|-------|-------|-----------|
| Docs (comments/roxygen) | 34 in `R/`, 13 in `man/` | `fs_debias_gate.R:4,10,28,243,251`; `fs_debias_gate_methods.R:2,35,154`; `forestsearch_main.R:294,340,614,624,745,749,983,1868,1999,2798,2828,2918`; `gate_estimates_table.R:1,7,8,28,67`; `grf_main.R:259,410`; `grf_subg_harm_glm.R:490,702`; `forestsearch_helpers.R:1033,1073,1091,1359,1439,1448,1722` |
| **User-visible string** | 1 | `gate_estimates_table.R:106` — `gt` subtitle `"Tier-2 de-biased estimate (multiplier bootstrap, %s draws)%s"` |
| History | 1 | `NEWS.md:42` — **do not rewrite** |

**`tier` FALSE POSITIVE — the most important one in this survey.**
`R/glm_effect_estimators.R:342-441` uses `tier1`, `tier2`, and the comments
"Tier 1 / Tier 2 / Tier 3" for the **risk-difference estimation fallback
ladder** (identity-link binomial GLM → G-computation from logistic → raw
proportions). This has nothing to do with FB/MR. It is the only place in `R/`
where `tier1`/`tier2` exist as actual identifiers, so a name-based rename would
hit precisely the wrong target. **Exclude `glm_effect_estimators.R` from every
tier rename.**

Other `tier` substring false positives, all unrelated: `frontier` (2911),
`frontier_params` (1616), `pareto_frontier_table` (424), `compute_frontier_cis`,
`plot_pareto_frontier`, `dina_frontier`, `grf_frontier_dmin`.

### 3.4 `dg` — false positives reported rather than dropped

`dg` as a debias-gate marker appears as `.fs_dg_*` (13 helpers, §5.6), `.dg_*`
locals (12 distinct), and `_dg` suffixes (`L_dg`, `combo_dg`, `tot_dg`,
`adj_dg`, `c_screen_dg`, `c_consistency_dg`, `.g_dg`).

Unrelated `dg` substring matches, judged and excluded: **`dgm`** (849
occurrences — data-generating mechanism, the single largest false-positive
class), `edge`, `position_dodge`, `dodge_width`, `nudge_y`, `budget`, `bridge`,
`Bridge`, `ledger`.

Bare `t1`/`t2` in `R/` are all unrelated: `test_hte_crump.R:252-299` (HTE test
result objects) and `oc_analyses.R:377` (`theta_1`).

### 3.5 `tests/`, `vignettes/`, metadata

**Verified: zero occurrences of `debias_gate`, `debias_gate_args`,
`gate_draws`, or `fs_debias_gate` in `tests/` or `vignettes/`.**

`tests/` mentions "gate" only as ordinary English in comments
(`test-cv-mixed-and-fatal.R:241,244`, `test-cv-no-subgroup-dina-grf.R:11,14`,
`test-no-subgroup-*.R`) and in `test-sg-focus-selection.R:91-142`, which tests
the *mapping* `.fs_dg_reselection_from_focus()` without naming the API. Class:
Docs/Internal. `test-return-shape-contract.R:212` "this test tier" is unrelated.

`_pkgdown.yml`, `DESCRIPTION`, `README.md`: **no** `gate`/`tier` occurrences —
nothing to change.

`NAMESPACE`: one line, `export(gate_estimates_table)` (`:138`). `fs_debias_gate`
is **not exported**.

---

## 4. Answers to the six questions of §5

### 5.1 Full API surface

"API" below means *a name a user writes or reads*. Note the split: three of
these are package API; four are simulation-harness names that only look like it.

**Package API (user sets):**

| name | kind | where |
|------|------|-------|
| `debias_gate` | formal of `forestsearch()`, default `FALSE` | `main.R:984` |
| `debias_gate_args` | formal, `list()` | `main.R:985` |
| `debias_gate_args$t_gate` | list element, effect scale, `NULL` → near-null | `main.R:2875` |
| `debias_gate_args$gate` | `"point"` / `"ci"` | `main.R:2876` |
| `debias_gate_args$reselection` | overrides the `sg_focus`-derived rule | `main.R:2877` |
| `debias_gate_args$selection_rule` | `"neighborhood"`/`"pareto"`/`"both"` | `main.R:2881` |
| `debias_gate_args$draws` | `2000L` | `main.R:2882` |
| `debias_gate_args$multiplier` | `"poisson"`/`"gaussian"`/`"rademacher"` | `main.R:2883` |
| `debias_gate_args$include_complement` | `TRUE` as passed | `main.R:2884` |
| `debias_gate_args$ci_method` | `"ij"`/`"wald"` | `main.R:2885` |
| `debias_gate_args$seed` | | `main.R:2886` |
| `debias_gate_args$family_native_neighborhood` | DINA/GRF only | `main.R:1896`, `:2027` |
| `gate_estimates_table()` | **exported** function | `NAMESPACE:138` |
| `c_gate` | formal of `fs_fdr_report()`, `c(0.9, 1.0, 1.25)` | `fs_fdr_report.R:177` |

**Package API (user reads):**

| name | kind |
|------|------|
| `out$debias_gate` | the whole returned list (16 fields, §5A.4.1) |
| `out$harm_flag_debiased` | logical / `NA` |
| `rep$fdr$c_gate` | column of the `fs_fdr_report()` result |

**Not package API — simulation-harness locals** (they never appear in `R/`, so a
package rename does not reach them, and renaming them does not touch the
package): `gate_draws` (2910 occurrences, 211 files), `gate_ok` (501, 177),
`run_debias_gate` (109), `gate_sec` (91), `gate_naive_est` (6), `run_tier1`
(3162), `nb_boots` (7033 — but this *is* also a package formal of
`forestsearch_bootstrap_dofuture()`), and the `t1_*`/`t2_*` record columns.

**Internal:** `fs_debias_gate()` (not exported), the 13 `.fs_dg_*` helpers,
`.fs_apply_debias_gate()`, all `.dg_*` locals.

### 5.2 Which occurrences are FB and which are MR?

Every "Tier-1" occurrence in `R/` refers to the full bootstrap
(`forestsearch_bootstrap_dofuture()`); every "Tier-2" refers to the multiplier
gate (`fs_debias_gate()`). **No mislabelled site was found.** Two are worth
flagging as *imprecise* rather than wrong:

* `main.R:2828` "faithful to the Tier-1 search" — the family is *not* faithful
  (F3). The label is right; the claim is not.
* `fs_debias_gate.R:28`, `:243`, `:251` call the IJ variance "the leading-order
  analogue of the Tier-1 interval" — accurate, but note the fallback chain
  (`ij` → `ij_raw` → `wald_fallback`) means the realised interval is not always
  that analogue. `ij_source` records which was used.

### 5.3 Does `fs_debias_gate()` do FB, MR, or both?

**MR only.** Confirmed by reading it end to end (`R/fs_debias_gate.R:256-422`):
its only stochastic input is `.fs_dg_multipliers()` (`:307`), it never refits a
bootstrap replicate, and it never calls `forestsearch()`. FB lives entirely in
`forestsearch_bootstrap_dofuture()` / `bootstrap_analysis_dofuture.R`, driven by
`nb_boots`, which re-runs the whole search per replicate
(`bootstrap_analysis_dofuture.R:585`).

**The structural concern of §5.3 does not arise; no split is needed.**
`debias_gate_args$draws` and `ci_method = "ij"` do indicate MR, as suspected.
Corroborating evidence from the author's own harness: the simulation cells
already label the two paths `FB` and `MR` —
`est_keys <- c(oracle = "or", naive = "nv", FB = "t1", MR = "t2")`
(`quarto/simulations/gbsg_redux/fs_t1_t2_*_batch_1_20.qmd:880`) — so the FB/MR
vocabulary is already adopted at the display layer, over `t1`/`t2` column
prefixes.

### 5.4 Is the `gate` argument accurately named?

**Yes. It keeps the word.** `gate` selects between two harm-confirmation
decision rules and nothing else (`fs_debias_gate.R:337-339`):

```r
gate_cmp <- if (log_scale) log(t_gate) else t_gate
ci_lo_1s <- beta_deb - stats::qnorm(0.95) * se
flag <- if (gate == "point") (beta_deb >= gate_cmp) else (ci_lo_1s >= gate_cmp)
```

* `"point"` — the de-biased point estimate must clear `t_gate`.
* `"ci"` — the one-sided 95% selection-adjusted lower bound must clear `t_gate`.

That is sense 4 exactly: a pass/fail applied to a completed result. `t_gate`
(the threshold) and `c_gate` (the vector of thresholds swept by
`fs_fdr_report()`) are the same sense and are likewise accurately named.
Verified rather than assumed, as the brief asks.

### 5.5 Blast radius per API name

Occurrence counts (not file counts) over source files only:

| name | `R/` | `man/` | `tests/` | `vignettes/` | `quarto/` | `dev/` | quarto files |
|------|-----:|-------:|---------:|-------------:|----------:|-------:|-------------:|
| `debias_gate` | 32 | 9 | **0** | **0** | 815 | 55 | 218 |
| `debias_gate_args` | 36 | 4 | **0** | **0** | 383 | 14 | 207 |
| `fs_debias_gate` | 7 | 14 | **0** | **0** | 152 | 26 | — |
| `t_gate` | 16 | 5 | 0 | 0 | 59 | 4 | — |
| `harm_flag_debiased` | 4 | 1 | 0 | 0 | 43 | 9 | — |
| `family_native_neighborhood` | 5 | 1 | 0 | 0 | 66 | 0 | — |
| `gate_estimates_table` | 1 | 4 | 0 | 0 | 43 | 0 | — |
| `c_gate` | 12 | 4 | 0 | 0 | 7 | 0 | 2 |
| `gate_draws` | **0** | 0 | 0 | 0 | 2910 | 4 | 211 |
| `gate_ok` | **0** | 0 | 0 | 0 | 501 | 8 | 177 |
| `run_tier1` | **0** | 0 | 0 | 0 | 3162 | 0 | — |

Three facts dominate the deprecation question:

1. **`tests/` and `vignettes/` are entirely clean.** No test or vignette will
   break. There is also therefore **no regression test covering `debias_gate`** —
   which is why F1 (a silent no-op on the default GLM path) went unnoticed.
2. **`quarto/` carries the whole load** — ~218 files for `debias_gate` alone.
   Many are near-duplicate variants of the same analysis (8 near-identical
   ACTG175 files, 6 GBSG, 201 files in `gbsg_redux/`), so the *distinct edit*
   count is far smaller than the file count. These are the author's own
   notebooks, not third-party code.
3. **Combined with F7 (nothing shipped in 0.1.0), no rename here breaks a CRAN
   user.** The choice is a mechanical-cost question, not a compatibility one.

### 5.6 The `.fs_dg_*` prefix

`dg` = "debias gate" at every site. Thirteen helpers, all internal, all
`@keywords internal`, all with generated `man/dot-fs_dg_*.Rd` pages:

| helper | file:line | what it does | sense |
|--------|-----------|--------------|-------|
| `.fs_dg_gate_null` | `fs_debias_gate.R:39` | near-null default for `t_gate` | **4** |
| `.fs_dg_pieces` | `:58` | Cox/GLM influence dispatch | 3 |
| `.fs_dg_assemble` | `:85` | build the `N x S` influence matrix | 3 |
| `.fs_dg_multipliers` | `:107` | draw MR multipliers | 3 (MR) |
| `.fs_dg_select` | `:129` | apply a re-selection rule on one draw | 2/3 |
| `.fs_dg_ij_var` | `:171` | IJ variance | 3 |
| `.fs_dg_se_from_ij` | `:195` | SE with fallback chain | 3 |
| `.fs_dg_members_from_conj` | `methods.R:26` | row indices for a conjunction | 3 |
| `.fs_dg_restrict_native` | `:61` | band the family on the native statistic | 3 |
| `.fs_dg_family_from_table` | `:87` | family from `(v1,d1,c1,v2,d2,c2)` | 3 |
| `.fs_dg_reselection_from_focus` | `:130` | `sg_focus` → re-selection rule | 2 |
| `.fs_dg_spec` | `:203` | build the Cox/GLM `spec` | 3 |
| `.fs_apply_debias_gate` | `:164` | wrapper + `tryCatch` | 3 |

Twelve of thirteen are sense 3 (MR machinery); only `.fs_dg_gate_null` is sense
4. **Proposed prefix: `.fs_mr_*`**, with `.fs_dg_gate_null` → `.fs_mr_gate_null`
(it computes the *gate*'s null default, so it keeps `gate` in the stem while
moving to the `mr` family). All internal ⇒ safe, no deprecation, one commit.

---

## 5A. Functional specification of the current "gate"

*Written to stand alone as a description of behaviour at 0.2.0, independent of
the rename question. Every claim cites file:line; claims marked **[executed]**
were confirmed by running the code — see `dev/terminology-work/R/1*.R`.*

### 5A.1 What is computed

`fs_debias_gate(df, candidates, spec, selected_members, ...)`
(`R/fs_debias_gate.R:256-422`), in order:

1. **Normalise the family** (`:272-283`). Ensure the observed selected subgroup
   is a member: match by `setequal()` against each candidate's row indices; if
   absent, append it as `".selected_H"`.
2. **Fit each candidate once** (`.fs_dg_assemble`, `:85-102`). For candidate
   *g*, fit the treatment effect on its members and scatter the per-subject
   treatment `dfbeta` into column *g* of an `N x S` matrix `B` (zero off-subgroup).
   Retain `beta_hat`, `sigma_D` (robust SE), size. Drop candidates with `< 6`
   members, a failed fit, or a `dfbeta`/index length mismatch (`:92-94`).
   Dispatch is Cox vs GLM via `.fs_dg_pieces()` → `.consistency_cox_pieces()` /
   `.consistency_glm_pieces()` — the same influence pieces the resample
   consistency engine uses.
3. **Admissibility threshold** (`:302-304`):
   `z <- qnorm((1 + p_star)/2)`; `t_g <- pmax(c_screen, c_consistency + z * sigma_D)`.
   A per-candidate floor combining the screening threshold with a
   normal-approximation stand-in for "consistency rate ≥ `p_star`".
4. **Draw multipliers** (`:307`). `Xi` is `N x draws`; `poisson` is
   `rpois(n, 1) - 1` (`:111`) — mean-zero, unit-variance, mimicking
   nonparametric-bootstrap multiplicities.
5. **Perturb and re-select** (`:308-319`). `P <- crossprod(B, Xi)` gives
   `D_g(b)`; `beta_star <- beta_hat + P`. Per draw: `pass <- which(bs >= t_g)`;
   if non-empty, `.fs_dg_select()` picks the winner under `reselection`; record
   `sel_bias[b] <- P[winner, b]`.
6. **Assemble the correction** (`:322-327`):
   `selection_bias = mean_b D_{H*_b}(b)`, `fixed_bias = mean_b D_{Hhat}(b)`,
   `beta_deb = beta_naive - selection_bias - fixed_bias`.
7. **IJ variance** (`:330-335`).
   `r_b = (selection_bias + fixed_bias) - D_{H*_b}(b) - D_{Hhat}(b)`;
   `.fs_dg_ij_var()` returns `tilde_V = sum_i cov_i^2` with
   `cov_i = (1/B) sum_b (Xi_ib - mean_b Xi_ib) r_b`, and the Wager-corrected
   `hat_V = tilde_V - (N/B) * mean(r_b^2)` (`:171-182`).
8. **Resolve the SE** (`.fs_dg_se_from_ij`, `:195-201`): `hat_V` if positive,
   else `tilde_V`, else `sigma_D`; recorded in `ij_source` as
   `"ij"` / `"ij_raw"` / `"wald_fallback"`.
9. **Harm confirmation** (`:337-339`), per §5.4.
10. **Complement** (`:348-400`), optional — see below.

Answers to the specific sub-questions:

* **What is de-biased, on what scale?** The subgroup's treatment effect
  coefficient, on the **working (comparison) scale** — log for ratio measures
  (HR/OR/RR/IRR), identity for differences (RD/MD). `log_scale` comes from the
  influence pieces (`:97`). All arithmetic — bias subtraction, IJ variance, CI
  construction, the gate comparison — happens there; `to_eff()` (`:301`)
  exponentiates only on output. **[executed]** `log_scale = TRUE`,
  `measure = "HR"` on the Cox path.
* **What interval, by which method?** A symmetric Wald interval on the working
  scale, exponentiated: `beta_deb ± qnorm(0.975) * se` (`:409-411`), plus a
  one-sided `lower_1s` at `qnorm(0.95)` (`:412`) which is what `gate = "ci"`
  uses. `ci_method = "ij"` (default) takes `se` from step 8; `"wald"` uses
  `sigma_D`. The `naive` interval always uses `sigma_D` (`:406-408`).
  **[executed]** `est 1.15 (0.811, 1.62)`, `se_ij 0.177` vs `se_wald 0.138`,
  `ij_source "ij"`, `ij_draws 298` of 300 — the IJ interval is materially wider
  than the robust-SE one.
* **What does `gate` select among?** `"point"` and `"ci"` — §5.4.
* **What is `t_gate`?** The harm-confirmation threshold on the **effect** scale
  (not the working scale — it is converted at `:337`). `NULL` → `.fs_dg_gate_null()`:
  `1` for `OR`/`RR`/`IRR`/`HR`, `0` otherwise (`:39-42`), applied at `:286`. The
  roxygen (`:224-226`) explains the choice: near the null rather than at the
  screen, because the correction over-shrinks true effects. **[executed]**
  `gate$t_gate = 1` on the Cox path.
* **What is `include_complement` for, and does it change the returned
  structure?** It additionally de-biases the complement (everyone outside the
  selected subgroup). Because the complement is *induced* by the selection
  rather than chosen, its bias is the perturbation of the complement *of the
  re-selected winner* on each draw (`:364-371`); complements are fit only for
  candidates that actually win, plus the selected one (`:352-363`). **Yes, it
  changes the structure**: `complement` is `NULL` when `FALSE`, otherwise a
  5-element list (`naive`, `debiased`, `selection_bias`, `fixed_bias`, `n`), or
  `list(note = "complement subgroup could not be fit")` on failure (`:398`).
  **[executed]** present, with its own `ij_source`/`ij_draws`.
* **How many draws; what is `multiplier = "poisson"` doing?** `draws = 2000L`
  by default at every caller (`main.R:2882`, `methods.R:185`). `"poisson"` is
  `rpois(n, 1) - 1`: centred Poisson(1) weights are the multiplier-bootstrap
  analogue of multinomial resampling multiplicities, which is why it is the
  default over `gaussian`/`rademacher`.

### 5A.2 Is it always calculated?

Established by execution, not by reading defaults.

1. **Default.** `debias_gate = FALSE` (`main.R:984`); `debias_gate_args = list()`.
   **Off out of the box.** **[executed]** `formals(forestsearch)$debias_gate` →
   `FALSE`.
2. **All three engines have a branch, but by different routes.**
   DINA (`main.R:1868-1916`) and GRF (`:1999-2048`) build a candidate table,
   optionally band it with `.fs_dg_restrict_native()`, convert it with
   `.fs_dg_family_from_table()`, and call the shared wrapper
   `.fs_apply_debias_gate()`. Consistency (`:2798-2900`) instead enumerates the
   full `≤ maxk` combination space from the cut matrix `Z` with the search's own
   helpers (`generate_combination_indices`, `get_covs_in`,
   `get_subgroup_membership`) and calls `fs_debias_gate()` **directly**.
   **Why the paths differ:** the wrapper's contract starts from a pre-built
   `candidates` list, which suits DINA/GRF (they already have candidate tables
   with `(v1,d1,c1,v2,d2,c2)` columns). The consistency branch has no such
   table — it must reconstruct the family from `Z` — and it also needs
   path-specific `spec`/threshold logic (internal `Y`/`Event`/`Treat` names and
   `log(hr.threshold)` on the Cox path vs user names and `effect_threshold` on
   the GLM path, `:2850-2867`). Two consequences of not sharing the wrapper:
   the consistency branch passes the search's own `selection_rule` (`:2881`)
   where the wrapper defaults to `"neighborhood"` (`methods.R:167`), and it
   carries the extra eligibility test in item 5 that DINA/GRF do not have.
   **[executed]** consistency ✓, DINA ✓ (`n_family = 37`), GRF ✗ (failed — F4).
3. **No subgroup identified → clean skip.** The guard requires
   `!is.null(sg.harm) && !is.null(grp.consistency) && !is.null(grp.consistency$sg.harm.id)`
   (`:2819-2821`). **[executed]** with `hr.threshold = 3.0`,
   `pconsistency.threshold = 0.99`: `sg.harm` `NULL`, `debias_gate` `NULL`,
   `harm_flag_debiased` `NA`, **no error and no warning**. DINA gates on
   `isTRUE(dsel$found)` (`:1872`); GRF omits a `found` test and relies on
   `sg.harm.id` (`:2005`).
4. **Under `sg_focus = "maxeff"`** the gate still runs, and
   `.fs_dg_reselection_from_focus()` maps it to the gate's `"maxeff"`.
   **[executed]** `gate_ran = TRUE`, `reselection = "maxeff"`, `n_family = 112`,
   `n_selected = 41`. But the gate's `"maxeff"` is
   `passers[which.max(beta[passers])]` — argmax **among passers** — while the
   `sg_focus = "maxeff"` contract is "argmax effect over ALL candidates,
   ungated" (`tests/testthat/test-sg-focus-selection.R:5`). The recorded
   mismatch, confirmed. Note the same gate rule `"maxeff"` is *correct* for
   `maxeffCons`, whose contract is argmax-among-qualifiers — which is why one
   gate rule serves two foci (`methods.R:143-148`).
5. **When does the `tryCatch` fire?** Two distinct `NULL`-producing mechanisms,
   which a user cannot tell apart from the return value:
   * **Silent skip (no `tryCatch` involved).** The consistency branch's
     eligibility test `.dg_glm_ok || .dg_cox_ok` (`:2817-2820`). **[executed]**
     — this is **F1**: for binary and count outcomes under the default
     `consistency_method = "split"`, the block never executes and nothing warns.
   * **Caught failure.** `.fs_apply_debias_gate()`'s `tryCatch`
     (`methods.R:190-194`) and the consistency branch's own
     (`main.R:2887-2889`). **[executed]** — this is **F4**: the GRF branch on
     factor covariates raises `'>' not meaningful for factors` and warns
     `debias_gate (survival) failed: ...`.

   In both cases `debias_gate` is `NULL` and `harm_flag_debiased` is `NA`. The
   brief's concern is well founded and applies to *both* paths: `NA` read through
   `isTRUE()` becomes `FALSE`, i.e. "harm not confirmed", when the truth is
   "not computed". The one mitigation is that `harm_flag_debiased` is `NA` rather
   than `FALSE`, so a caller who tests `is.na()` can distinguish them — but
   nothing in the package or its docs directs them to.
6. **Bootstrap replicates and CV folds: yes, it is recomputed.** **[executed]**
   2 calls in a 2-replicate bootstrap, 2 in a 2-fold CV — see **F2**. The
   intended answer ("only once on the original analysis") is what the *reported*
   result reflects, but not what the code executes.

### 5A.3 How is it utilised

1. **Does anything downstream branch on the result?** **No — the maintainer's
   position is honoured inside `forestsearch()`.** `harm_flag_debiased` is
   written at three sites (`main.R:1918`, `:2050`, `:2920`) and **read nowhere
   in `R/`** (only `man/forestsearch.Rd:729` documents it). The gate block sits
   at `:2798-2900`, after `sg.harm` is assigned at `:2751`, so it cannot
   influence selection even in principle. **[executed]** with the gate ON vs
   OFF on identical data: `sg.harm`, `grp.consistency$sg.harm.id`, `max_sg_est`,
   and `dim(df.est)` all identical; `names(out)` identical; the only difference
   is `debias_gate` `NULL` → populated and `harm_flag_debiased` `NA` → `TRUE`.
   **No defect here.**
2. **Purely reported, or feeding a computation?** Reported, with one genuine
   consumer: `fs_fdr_report()` sets `fs_params$debias_gate <- TRUE`
   (`fs_fdr_report.R:232`, comment "B needs the de-biased HR"), reads
   `fs$debias_gate$debiased$est` (`:375-377`), and scores declaration rates
   against each `c_gate` threshold (`:438-456`). That is *post-hoc analysis of a
   completed identification*, not part of identification — consistent with the
   maintainer's framing. `gate_estimates_table()` is presentation only.
3. **Does the re-selection replay the identifier faithfully?** **No.** Two
   dimensions:

   *Selection rule* — faithfully mapped. **[executed]**
   `.fs_dg_reselection_from_focus()` over all 11 accepted foci:

   | `sg_focus` | engine `"consistency"` | engine `"effect"` |
   |------------|------------------------|-------------------|
   | `hr`, `eff`, `maxcons` | `maxcons` | `maxeff` |
   | `maxeff` | `maxeff` | `maxeff` |
   | `maxeffCons` | `maxeff` | `maxeff` |
   | `maxSG` / `minSG` | `maxSG` / `minSG` | same |
   | `hrMaxSG` / `hrMinSG` | `effMaxSG` / `effMinSG` | same |
   | `effMaxSG` / `effMinSG` | `effMaxSG` / `effMinSG` | same |

   Every focus reaches a named rule; the fall-through (`methods.R:151`) is never
   the only path for an accepted value. Verified end-to-end: `sg_focus = "maxSG"`
   → `gate$reselection = "maxSG"`; `"maxeffCons"` → `"maxeff"`; `"hr"` →
   `"maxcons"`.

   *Admissibility* — **not** faithfully replayed, for any focus. Per-focus
   tabulation as §5A.3.3 requests:

   | `sg_focus` | consistency floor: identifier | consistency floor: gate | `d0.min`/`d1.min` | `max_subgroups_search` |
   |------------|------------------------------|-------------------------|-------------------|------------------------|
   | `hr`/`eff`/`maxcons` | yes (qualifiers) | yes (`t_g`) — match | **gate omits** | **gate omits** |
   | `maxeffCons` | yes (qualifiers) | yes — match | **gate omits** | **gate omits** |
   | `maxSG`/`minSG` | yes (qualifiers) | yes — match | **gate omits** | **gate omits** |
   | `hrMaxSG`/`hrMinSG`/`effMaxSG`/`effMinSG` | yes (qualifiers) | yes — match | **gate omits** | **gate omits** |
   | `maxeff` | **no** (ungated) | **yes** — **mismatch** | **gate omits** | **gate omits** |

   The consistency-floor mismatch is unique to `maxeff` (the previously recorded
   one). The `d0.min`/`d1.min` and `max_subgroups_search` omissions affect
   **every** focus and were not previously recorded — see **F3**. All three
   inflate `selection_bias`. Out of scope to fix here.

   A fourth, milder divergence: DINA and GRF pass `c_consistency = 0`
   (`main.R:1910`, `:2042`) where the consistency branch passes the resolved
   `consistency_threshold` (`:2866`). With `c_consistency = 0` the floor becomes
   `pmax(c_screen, z * sigma_D)` — still a significance-flavoured floor, not an
   absent one.

### 5A.4 How is it output

**1. Full field inventory** of `out$debias_gate` — 16 top-level fields,
enumerated by `str()` on a live result **[executed]**:

| field | type | meaning | scale |
|-------|------|---------|-------|
| `selected_index` | int | column of the selected candidate in the assembled family | — |
| `selected_label` | chr | that candidate's label | — |
| `measure` | chr | `spec$effect_measure` (e.g. `"HR"`) | — |
| `log_scale` | lgl | working scale is log | — |
| `ci_method` | chr | `"ij"` / `"wald"` as resolved | — |
| `naive$est/lower/upper` | num | uncorrected effect, robust-SE 95% CI | **effect** |
| `debiased$est/lower/upper` | num | de-biased effect + 95% CI | **effect** |
| `debiased$lower_1s` | num | one-sided 95% lower bound (used by `gate = "ci"`) | **effect** |
| `debiased$se` | num | SE actually used for the CI | **working** |
| `debiased$se_ij` | num | IJ SE | **working** |
| `debiased$se_wald` | num | robust `sigma_D` | **working** |
| `debiased$var_ij` | num | IJ variance | **working²** |
| `debiased$ij_source` | chr | `"ij"` / `"ij_raw"` / `"wald_fallback"` | — |
| `debiased$ij_draws` | int | draws with a re-selected winner | — |
| `selection_bias` | num | `mean_b D_{H*_b}(b)` | **working** |
| `fixed_bias` | num | `mean_b D_{Hhat}(b)` | **working** |
| `selection_rate` | num | fraction of draws yielding a winner | — |
| `complement` | list/NULL | same shape, plus `n`; or `note` | as above |
| `gate$t_gate` | num | harm-confirmation threshold | **effect** |
| `gate$type` | chr | `"point"` / `"ci"` | — |
| `gate$reselection` | chr | rule used | — |
| `gate$selection_rule` | chr | band construction | — |
| `gate$multiplier` | chr | multiplier family | — |
| `gate$draws` | int | draws requested | — |
| `harm_flag` | lgl | **the gate decision** (`isTRUE`-coerced) | — |
| `n_family` | int | candidates successfully fit | — |
| `n_selected` | int | size of the selected subgroup | — |
| `timing_seconds` | num | re-selection loop only | — |

Early-return variant (`:292-297`, selected subgroup unfittable): only
`selected_index = NA`, `selected_label`, `harm_flag = NA`, `gate`, `note`,
`n_family`. **Consumers must tolerate this shorter shape** — `gate_estimates_table()`
does not (it would index `g$debiased$est` as `NULL`).

Note a **mixed-scale trap**: `est`/`lower`/`upper`/`t_gate` are on the effect
scale while `se`/`se_ij`/`se_wald`/`selection_bias`/`fixed_bias` are on the
working (log) scale, in the same flat list, with no suffix marking which.

**2. Which fields does user code read?** From the `gbsg_redux` simulation cells
(recorder at `fs_t1_t2_*_batch_1_20.qmd:438-470`): `selected_label`,
`n_selected`, `harm_flag`, `debiased$ij_source`, `naive$est/lower/upper`,
`debiased$se_wald`, `debiased$est/lower/upper`, `debiased$se_ij`, and the whole
`complement$naive`/`complement$debiased` block. `gate$draws` is read by
`gate_estimates_table.R:96`. `debiased$est` by `fs_fdr_report.R:377`.
**Every one of these is therefore API** by the brief's definition and carries the
same rename question as a formal argument. `tests/` reads none of them.

**Tracing `gate_ok`** (as §5A.4.2 asks): it is produced at
`fs_t1_t2_*_batch_1_20.qmd:454` —

```r
rec$gate_ok <- as.integer(!is.null(g))
```

so **`gate_ok` records that the gate *computed*, not that it *passed*.** The
harm decision is stored separately as `rec$t2_harm_flag <- as.integer(isTRUE(g$harm_flag))`
(`:456`). The harness is correct — the comment at `:444-446` says so explicitly
("gate availability is tracked separately in `gate_ok` so a gate failure never
masks a true detection") — but the *name* reads as "the gate was passed". Given
F1 and F4, a `gate_ok = 0` row means "no MR result", which is precisely the
distinction that must not be blurred. Recommend renaming on clarity grounds
(§6), independent of the FB/MR work.

**3. `print()` / `summary()` / plotting.** **None.** `R/forestsearch_methods.R`
contains zero `gate` occurrences, and no `print`/`summary`/`plot` method
displays gate output. The only user-visible label produced by the package is the
`gt` subtitle at `gate_estimates_table.R:106`:
`"Tier-2 de-biased estimate (multiplier bootstrap, %s draws)%s"` — plus the
`cat()` diagnostic at `main.R:2893-2898` under `quiet = FALSE`:
`"De-biased gate: HR = 1.150 (gate point 1.00) -> consistent with harm"`.
Both are strings needing the same treatment as an argument name.

### 5A.5 FB and MR: which code path is which

1. **FB (full bootstrap)** — `forestsearch_bootstrap_dofuture()`
   (`R/bootstrap_dofuture_main.R`, analysis in `R/bootstrap_analysis_dofuture.R`),
   driven by `nb_boots`, re-running the entire search per replicate
   (`bootstrap_analysis_dofuture.R:585`). Supporting files:
   `bootstrap_calculations_helpers.R`, `bootstrap_summaries_helpers.R`,
   `summarize_bootstrap_results.R`. Yes — `nb_boots` in the simulation `.qmd`
   is exactly this path.
2. **MR (multiplier resampling)** — `fs_debias_gate()` (`R/fs_debias_gate.R`)
   plus its application layer (`R/fs_debias_gate_methods.R`), driven by
   `debias_gate_args$draws` (`gate_draws` in the sim).
3. **Does `fs_debias_gate()` span both?** **No — MR only.** See §5.3. No split
   needed; no structural finding.
4. **Are "tier-1"/"tier-2" used in the package at all?** Only as **prose**: 34
   comment/roxygen occurrences in `R/`, 13 in `man/`, one user-visible string
   (`gate_estimates_table.R:106`), one `NEWS.md` history line. **There is no
   `tier1`/`tier2` identifier in `R/` that means FB or MR** — the only such
   identifiers belong to the unrelated RD fallback ladder in
   `glm_effect_estimators.R` (§3.3). The `run_tier1` / `t1_*` / `t2_*`
   identifiers live exclusively in the simulation `.qmd` files, which already
   *display* them as `FB`/`MR`. **Scope of the tier rename inside the package is
   therefore documentation-only, plus one string.**

---

## 6. Proposed rename map

**Proposals, not decisions.** Where several defensible options exist they are
listed with the tradeoff and no choice is made.

Deprecation column reflects **F7** (confirmed, not inferred): nothing in the
`debias_gate` family shipped in CRAN 0.1.0, so "N" means *no CRAN obligation*,
not *no work* — the author's `quarto/` notebooks still need mechanical updating.

### 6.1 Package API

| current | proposed | class | call sites (R/man/tests/vign/quarto) | deprecate? | rationale |
|---------|----------|-------|--------------------------------------|-----------|-----------|
| `debias_gate` | `mr_inference` · or `mr_debias` · or `run_mr` | API | 32/9/0/0/815 | **N** (F7) | names sense 3, which is what it does. **Options:** `mr_inference` is accurate and reads as a noun-phrase toggle; `mr_debias` foregrounds the correction but under-sells the interval; `run_mr` reads best as a logical flag but is opaque without the FB/MR glossary. |
| `debias_gate_args` | `mr_args` (match whichever stem is chosen) | API | 36/4/0/0/383 | **N** | must track the flag |
| `harm_flag_debiased` | `harm_confirmed` · or keep | API (read) | 4/1/0/0/43 | **N** | already sense-4 and accurate; only "debiased" is loose. **Weakest case for change** in the table. |
| `out$debias_gate` | `out$mr` (track the flag) | API (read) | — | **N** | see F7 |
| `gate_estimates_table()` | `mr_estimates_table()` | **API, exported** | 1/4/0/0/43 | **N** | renders MR estimates, not a gate decision. Only exported name in the family. |
| `debias_gate_args$gate` | **keep `gate`** | API | — | — | accurately names sense 4 (§5.4) |
| `debias_gate_args$t_gate` | **keep `t_gate`** | API | 16/5/0/0/59 | — | sense 4 threshold |
| `c_gate` (`fs_fdr_report()`) | **keep** | API | 12/4/0/0/7 | — | sense 4 threshold sweep |
| `debias_gate_args$draws` etc. | **keep** | API | — | — | already sense-3-neutral |
| `family_native_neighborhood` | **keep** | API | 5/1/0/0/66 | — | no gate/tier content |

### 6.2 Internal

| current | proposed | class | deprecate? | rationale |
|---------|----------|-------|-----------|-----------|
| `fs_debias_gate()` | `fs_mr_inference()` · or `fs_mr_debias()` | Internal (**not exported**) | N | not exported ⇒ free. Tradeoff mirrors the flag: `fs_mr_inference()` covers both the correction and the interval (which the function does); `fs_mr_debias()` is shorter but names only half the job. Keeping a short form (`fs_mr()`) is a third option — least typing, least self-describing. |
| `.fs_dg_*` (12 helpers) | `.fs_mr_*` | Internal | N | `dg` = "debias gate" everywhere; 12 of 13 are MR machinery (§5.6) |
| `.fs_dg_gate_null` | `.fs_mr_gate_null` | Internal | N | the one sense-4 helper — keeps `gate` in the stem |
| `.fs_apply_debias_gate()` | `.fs_apply_mr()` | Internal | N | wrapper |
| `.dg_*` locals, `_dg` suffixes | `.mr_*`, `_mr` | Internal | N | cosmetic, zero risk |
| `.dg_glm_ok` / `.dg_cox_ok` | `.mr_glm_ok` / `.mr_cox_ok` | Internal | N | see F1 — worth a comment either way |

### 6.3 Docs and strings

| current | proposed | class | note |
|---------|----------|-------|------|
| "Tier-2" in comments/roxygen (34 in `R/`) | "multiplier resampling (MR)" | Docs | mechanical; regenerate `man/` with `devtools::document()` |
| "Tier-1" in comments/roxygen | "full bootstrap (FB)" | Docs | same |
| `man/*.Rd` (13 tier, 108 gate lines) | — | Docs | **generated — never hand-edit** |
| `gate_estimates_table.R:106` subtitle | "MR de-biased estimate (multiplier bootstrap, %s draws)%s" | **User-visible string** | changes rendered output |
| `main.R:2893` `cat()` "De-biased gate: ..." | e.g. "MR harm confirmation: ..." | **User-visible string** | `quiet = FALSE` only |
| `NEWS.md:42` | **leave** | **History** | do not rewrite existing entries |
| Senses 5 & 6 (20 occurrences) | **leave** | Docs | tolerance gates and control-flow guards — §2 |
| `glm_effect_estimators.R:342-441` `tier1`/`tier2` | **leave** | Internal | **false positive** — RD fallback ladder |

### 6.4 Simulation harness (not package API)

Reachable only by editing `quarto/`; no package change involved. Listed because
the brief asks for the returned-object indexers.

| current | proposed | files | rationale |
|---------|----------|-------|-----------|
| `gate_draws` | `mr_draws` | 211 | counts MR draws, not a gate quantity — the brief's own diagnosis, confirmed |
| `gate_ok` | `mr_ok` · or `mr_computed` | 177 | means "MR produced a result"; `mr_computed` is unambiguous, `mr_ok` keeps the diff small but retains the pass/fail overtone |
| `run_tier1` | `run_fb` | — | FB toggle |
| `t1_*` / `t2_*` record columns | `fb_*` / `mr_*` | — | `est_keys` already displays `FB`/`MR`, so this closes an existing gap |
| `gate_sec`, `gate_naive_est`, `run_debias_gate` | `mr_sec`, `mr_naive_est`, `run_mr` | — | consistency |
| `*_t1_t2_*.qmd` filenames (201 in `gbsg_redux/`) | `*_fb_mr_*` | 201 | **recommend not renaming**: filenames are keyed to `.rds` result artefacts; churn for no readability gain |

### 6.5 Collisions

Checked against the current namespace and `R/` identifiers:

* **No collision** for `mr_inference`, `mr_args`, `mr_debias`, `fs_mr_inference`,
  `fs_mr_debias`, `.fs_mr_*`, `mr_estimates_table`, `harm_confirmed`.
* **Caution, not collision:** the prefix `mr_` is already used in the *simulation*
  filenames `mr_coverage_sweep_*.qmd` and `_sim_mr_*.qmd`, and `mrct_simulation.R`
  concerns *multi-regional clinical trials* — an unrelated "MR". Nothing breaks,
  but "MR" is now overloaded in the repo; a one-line glossary in the package
  docs would earn its keep.
* `run_mr` as the flag name would sit oddly beside `run_tier1`→`run_fb` in the
  sim, where `run_*` means "execute this stage" rather than "compute this
  quantity" — a mild argument for `mr_inference` over `run_mr`.

---

## 7. Suggested execution order

Each step is independently committable and leaves the package loadable. Steps
1–3 carry no user-visible change at all.

| step | scope | independent? | shim? | notes |
|------|-------|--------------|-------|-------|
| **1** | `.fs_dg_*` → `.fs_mr_*`, `.fs_apply_debias_gate` → `.fs_apply_mr`, `.dg_*` locals | yes | no | pure-internal; `devtools::document()` regenerates the 13 `man/dot-fs_dg_*.Rd` (**delete the stale ones**) |
| **2** | "Tier-1"/"Tier-2" → FB/MR in comments and roxygen | yes | no | **exclude `glm_effect_estimators.R`**; leave `NEWS.md:42` |
| **3** | `fs_debias_gate()` → chosen internal name | yes | no | not exported; 7 `R/` sites + `man/` regeneration |
| **4** | `debias_gate` + `debias_gate_args` + `out$debias_gate` + `harm_flag_debiased` | **must land together** | no (F7) | the API step. All four are one contract; a partial rename leaves an incoherent surface |
| **5** | `gate_estimates_table()` → `mr_estimates_table()` | after 4 | no (F7) | only exported name; update `NAMESPACE` via `document()` |
| **6** | the two user-visible strings | with 5 | — | rendered output changes |
| **7** | `quarto/` mechanical update | after 4–6 | — | ~218 files for `debias_gate`; many are near-duplicates, so do it per-family (ACTG175, GBSG, `gbsg_redux`) and re-render one file per family as a check |
| **8** | sim-harness locals (`gate_draws`, `gate_ok`, `run_tier1`, `t1_*`/`t2_*`) | yes, any time | — | independent of the package entirely; `gate_ok` is worth doing on its own merits |

**Must land together:** step 4's four names. **Everything else is independent.**
**No step needs a deprecation shim** (F7) — but see §8 for the one caveat that
would change that.

**Recommended sequencing note.** Steps 1–3 are free and remove most of the
confusion the brief identifies. If the API rename (4–7) is deferred, the package
is still left in a strictly better state, with the FB/MR vocabulary consistent
in all prose and internals.

**Strong recommendation, outside this brief's scope:** add a regression test for
`debias_gate` before step 4. `tests/` currently has **zero** coverage of it
(§5.5), which is why F1 survives. Renaming an untested surface across 218 files
is the riskiest part of this plan, and the risk is removable.

---

## 8. What could not be determined, and why

1. ~~**Whether `debias_gate` truly never shipped.**~~ **RESOLVED 2026-07-29 —
   see F7.** Originally left open because the preceding task's brief forbade
   inspecting git history. On the maintainer's explicit request this was then
   confirmed two ways: the CRAN 0.1.0 reference manual contains the string
   `gate` zero times, and the 0.1.0-era git tree contains none of the names
   (first introduced 2026-06-08, ~11 weeks after the 2026-03-23 CRAN publish;
   no `0.1.x` patch release exists). **F7 stands; step 4 needs no shim.**
2. **Whether the GRF gate failure (F4) is universal or specific to factor
   covariates.** Reproduced once, on the synthetic survival DGM whose `stage`
   and `sex` are factors. I did not test an all-numeric GRF case, so I cannot
   say whether the GRF branch works at all in practice or only fails here.
3. **`continuous` outcomes and the gate.** My continuous DGM produced no
   subgroup under either `consistency_method`, so the F1 table's `continuous`
   rows are uninformative. The code path is identical to `binary`/`count`, so
   F1 should apply, but that is inference, not measurement.
4. **The true distinct-edit count in `quarto/`.** I counted occurrences and
   files (218 for `debias_gate`), but many files are near-identical variants
   (8 ACTG175, 6 GBSG, 201 in `gbsg_redux/`). I did not cluster them, so the
   real editing effort for step 7 is materially lower than 218 and I cannot put
   a number on it.
5. **Whether any `quarto/` notebook reads a field the early-return shape omits.**
   I confirmed the sim recorder guards on `!is.null(g)` before indexing, but did
   not audit all 218 files for unguarded `$debiased$est` access. The
   `gate_estimates_table()` case (§5A.4.1) shows the pattern exists in `R/`.
6. **`selected_label` vs `sg.harm` vocabularies.** **[executed]** the gate
   reports `selected_label = "q7.0 & q3.1"` (internal cut names from
   `colnames(Z)`) where `sg.harm` is `"!{stage} & {age <= 63}"`. These are the
   same subgroup in two notations. Whether that is intended or an oversight is a
   presentation question I could not resolve from the code; flagging it as a
   loose end rather than a finding.

---

## 9. Artefacts

All under `dev/terminology-work/`:

```
TERMINOLOGY_SURVEY.md            this document
R/00_guard_before.R              before-hash snapshot
R/00_guard_after.R               after-hash + fs_guard_verify()
R/10_spec_probe.R                §5A.1/5A.4 field inventory (Cox/consistency)
R/11_spec_probe2.R               §5A.2 engines x foci, no-subgroup, reselection map
R/12_spec_probe3.R               F1 (GLM split vs resample), F3 (family sizes), F2 (bootstrap)
R/13_spec_probe4.R               F2 (CV folds), §5A.3.1 no-branching check
out/hash_before.rds              133 digests
out/hash_after.rds               133 digests
out/env_capture.csv              R/BLAS/package versions
out/inv_gate_R_filtered.txt      208 classified gate occurrences in R/
out/engines_x_focus.rds          §5A.2.2 results
out/glm_split_vs_resample.rds    F1 results
out/fs3_cox_consistency.rds      the live result used for §5A.4.1
out/d_survival.rds               synthetic data for reproduction
```
