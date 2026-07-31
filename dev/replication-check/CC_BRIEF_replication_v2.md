# CC BRIEF — GBSG re-run, discrepancy catalogue, and code assessment

```
claude "Read dev/replication-check/CC_BRIEF_replication_v2.md and execute it."
```

**Supersedes `CC_BRIEF_replication.md`.** That brief carried three parts beyond
the goal (an RNG-neutrality gate, a simulation characterisation, and a repo-wide
pinning pass). They are dropped. This brief does two things: re-run the analysis on the current branch and
catalogue every difference against the legacy renders, then assess whether any
difference — or anything found along the way — indicates something in the code
that needs re-evaluating.

**Exact numerical agreement is not the objective.** Differences are observations
to be explained, not failures. For each one the question is: *is this intended
behaviour, or does it reveal a defect?* A clean replication that hides a latent
code problem is the worse outcome.

Run from a worktree on `feature/mr-in-replicates`. **Change no package code, no
`.qmd`, and no test.** Write only under `dev/replication-check/`. Hash package
sources before and after with `fs_hash_sources()` / `fs_guard_verify()` from
`dev/efficiency-eval/R/00_guard.R` and report the verdict.

---

## 1. Inputs

| role | file |
|---|---|
| notebook to render | `quarto/applications/gbsg/analysis_gbsg_cox_multimethod_psi_v2_2new.qmd` |
| **baseline B** (primary) | `quarto/applications/gbsg/analysis_gbsg_cox_multimethod_psi_v2_2new.html` |
| **baseline A** (legacy, as provided) | `quarto/applications/gbsg/analysis_gbsg_cox_multimethod_psi_v2_2A_linux.html` |
| baseline A source, faithful | `dev/replication-check/legacy_v2_2A_reconstructed.qmd` (see §2) |
| baseline B source, faithful | `dev/replication-check/v2_2new_rendered_source_prerename.qmd` (see §2) |

**Two comparisons are required.** See §5 for why baseline B is the primary one.

Confirmed current by the maintainer: `psi_v2_2new`. Do not render
`psi_v2_2`, `psi_v2_2A`, `psi_v2_2A_linux`, or `psi_v3a`.

Install the package before rendering — `devtools::install()`, not `load_all()`;
`doFuture` workers spawn separate R processes and need the installed package.
Record the commit SHA and `packageVersion("forestsearch")`.

---

## 2. A caveat about the legacy source

`analysis_gbsg_cox_multimethod_psi_v2_2A_linux.qmd` in the repo was rewritten in
place by the terminology purge — it now reads `run_mr` / `mr_inference`, and
differs by 81 lines from the source that actually produced
`..._v2_2A_linux.html`. Do not treat it as the legacy source.

The same applies to `..._v2_2new.qmd`: the on-disk copy is post-rename, but
`..._v2_2new.html` was rendered from pre-rename source.

Faithful reconstructions of both, recovered from the embedded Quarto source
modals in the respective HTML files, are provided at:

* `dev/replication-check/legacy_v2_2A_reconstructed.qmd`
* `dev/replication-check/v2_2new_rendered_source_prerename.qmd`

Use those if you need to inspect what produced either baseline. The **results**
baselines remain the `.html` files, which are unaffected.

---

## 3. Configuration — leave as-is

The legacy render and `v2_2new` share identical identification settings.
Confirm before rendering, change nothing:

```r
fs_sg_focus  <- "effMaxSG"   ;  fs_hr_threshold <- 1.0
fs_hr_consistency <- 1.0     ;  fs_max_subgroups_search <- 50
fs_conf_force <- c("er <= 0", "pgr <= 0")
fs_consistency_method <- "resample"
fs_subgroup_method <- "consistency"
NB <- 1000                   ;  mr_draws <- 5000L
run_cv <- FALSE              ;  run_loo <- FALSE
```

Note `run_cv` and `run_loo` are `FALSE` and every CV/LOO chunk is gated
`eval = run_cv` / `eval = run_loo`. Those sections do not execute — there are no
CV or LOO metrics to compare. Do not switch them on.

---

## 4. Render

```bash
quarto render quarto/applications/gbsg/analysis_gbsg_cox_multimethod_psi_v2_2new.qmd
```

Writes `gbsg_table2new_payload.rds`. Keep it.

---

## 5. The two baselines, and what each isolates

The two existing renders differ by **three lines of source**: `consistency_method`
added to the (A) call, and the payload filename. **Both are pre-rename** — their
embedded sources use `run_debias_gate` / `debias_gate`, with zero occurrences of
`mr_inference`. This has been verified, not assumed.

That makes the two comparisons cleanly separable:

| comparison | isolates | expectation |
|---|---|---|
| **new run vs baseline B** (`v2_2new.html`) | the gate→MR rename + `mr_in_replicates` | exact match on every quantity (confirm precondition below) |
| **new run vs baseline A** (`v2_2A_linux.html`) | the above + `consistency_method` | as B, except (A) differs as tabulated |

**Baseline B is the primary comparison**, because `consistency_method` is held
constant across it.

### Precondition: confirm `mr_in_replicates` is RNG-neutral (cheap, do first)

Commit `5bad3df7` changed what reaches a bootstrap replicate. Replicate
arguments are a copy of the parent fit's call (`args_FS_template <-
fs.est$args_call_all`), and `mr_inference` was previously *inherited* from that
call rather than overridden — so a fit with `mr_inference = TRUE`, as this
notebook has, ran MR inside every replicate. The commit adds
`if (!isTRUE(mr_in_replicates)) args_FS_boot$mr_inference <- FALSE`, default
`FALSE`. (A fit with `mr_inference = FALSE` was never affected.)

Reading the code, this **should** be RNG-neutral:

* the `%dofuture%` loop sets `.options.future = list(seed = TRUE)`, so each
  replicate has its own L'Ecuyer-CMRG substream — draws in replicate *i* cannot
  reach replicate *j*;
* `boot_index_mat` and `Ystar_mat` are built up front in the driver from
  `seed` and passed in, so resample selection is deterministic and independent of
  anything a replicate does.

`fs_mr_inference()` nonetheless takes `seed = NULL` and only calls `set.seed()`
when given one, and the notebook supplies none — so MR does draw from whatever
stream it lands in. Confirm neutrality rather than assuming it:

> With the outer seed fixed, run a short bootstrap twice (`nb_boots = 5` is
> enough) — once `mr_in_replicates = TRUE`, once `FALSE` — and compare
> `Ystar_mat`, the per-replicate selected subgroups, and the bias-corrected
> estimates.

If identical, MR is RNG-neutral and **exact agreement is the expectation for
every quantity below**. If not, say so in the findings and treat FB quantities as
expected-to-move; that is a result the maintainer needs before merging. Either
way, proceed — this is a precondition, not a gate.

### Established baseline results

| analysis | baseline A (`v2_2A_linux`) | baseline B (`v2_2new`) |
|---|---|---|
| FS main (`fs`) | `{er <= 0} & {size <= 35}`, 61 (8.9%) | identical |
| (A) DINA-screened (`fs_dina_screen`) | `{er <= 0} & {nodes <= 7}`, 61 (8.9%) | **`{er <= 0} & {size <= 35}`, 61 (8.9%)** |
| (C) DINA-selected (`fs_dina`) | `{grade3 >= 1} & {pgr <= 10}`, 89 (13.0%) | identical |
| (D) GRF-selected (`fs_grf`) | `{er <= 0}`, 82 (12.0%) | identical |
| (B) standalone DINA frontier | `grade3 >= 1.0`, `pgr >= 32.5` | identical |

Timings — A then B: FS search 25.3 s / 25.2 s; bootstrap 10.9 min / 10.9 min;
DINA-mode bootstrap 2.7 min / 2.7 min; GRF-mode bootstrap 3.3 min / 3.2 min.
Close enough that both were very likely run on the same host.

**The `consistency_method` effect is therefore already characterised**: switching
(A) from the `"split"` default to `"resample"` moved it off
`{er <= 0} & {nodes <= 7}` and onto the main FS subgroup
`{er <= 0} & {size <= 35}`, at unchanged size. That is what the source comment
(`ALIGN: match main run`) intended. Do not re-derive this — confirm the new run
agrees with baseline B on (A), and flag it if it does not.

## 6. Compare and catalogue

Run **both** comparisons. For each of the four analyses, report new vs
baseline B and new vs baseline A side by side, per quantity. Mark each as
`same` / `differs` — these are descriptive labels, not pass/fail:

* selected subgroup label — the primary criterion;
* subgroup size `n` and % of ITT;
* naive effect estimate and interval, H and Hc;
* bootstrap bias-corrected (FB) estimate and interval, H and Hc;
* MR de-biased estimate and interval, and `mr_harm_confirmed`;
* wall-clock: search, bootstrap, MR, and the MR-vs-bootstrap speed-up.

Source the new side from the payload (`P$table`, `P$meta$timings`), not by
scraping the HTML.

### What each comparison should show

**Against baseline B — expect exact agreement on every quantity**, subject to the
§5 precondition. Selected subgroup is the primary criterion and is settled before
any resampling, so a difference there is the single most important thing this
exercise can surface: report it at the top. If the precondition showed MR is not
RNG-neutral, FB quantities move for that reason and should be reported as
magnitudes rather than pass/fail.

**Against baseline A — expect (A) to differ, exactly as tabulated in §5, and
nothing else.** If a second analysis differs against A but not against B, that
is `consistency_method` reaching further than expected and should be called out.

The rename surface, for attribution: `run_debias_gate`/`gate_draws` →
`run_mr`/`mr_draws`, `debias_gate=` → `mr_inference=`,
`fit$debias_gate` → `fit$mr_inference`, `g$gate$t_gate` → `g$settings$t_confirm`,
`fs$harm_flag_debiased` → `fs$mr_harm_confirmed`,
`gate_estimates_table()` → `mr_estimates_table()`.

For each difference found, assign one of: *consistency_method (A)*,
*mr_in_replicates RNG shift* (only if the §5 precondition established one),
*rename* (spurious — the rename should be inert, so investigate),
*parallel scheduling*, or **unexplained**. Anything landing in
*unexplained* goes at the top of the report.

**Do not adjust the notebook to make a difference go away.** Report it.

---

## 6b. Code assessment — inherited replicate arguments

The `mr_inference` problem was structural, not a one-off: bootstrap replicate
arguments are a bulk copy of the parent fit (`args_FS_template <-
fs.est$args_call_all`), and only 12 of `forestsearch()`'s 85 formals are
overridden per replicate — `df.analysis`, `df.predict`, `grf_res`, `grf_cuts`,
`dina_res`, `dina_cuts`, `details`, `quiet`, `show_candidate_summary`,
`plot.sg`, `plot.grf`, and now `mr_inference`. The other 73 are inherited.

Most *should* be inherited — the bootstrap exists to re-run the same algorithm on
a resample, so thresholds, `sg_focus`, `maxk`, screening flags and variable names
all belong. The risky class is anything carrying a **pre-computed, full-data
object** or a **post-selection layer**. That class is exactly what the overrides
already cover, with one apparent gap.

### Confirmed gap: `ps_hat`

`forestsearch_cross_validation.R` nulls it in both CV paths (lines 396, 859)
— *"user-supplied ps_hat (wrong length for training fold)"* — and the CV summary
reports `ps_hat: re-estimated on each fold`.

`bootstrap_analysis_dofuture.R` never touches `ps_hat`, `ps_method`,
`ps_adjust_method` or propensity scores at all. A user-supplied `ps_hat` is
therefore inherited whole into every replicate. Because a bootstrap resample has
the same row count as the original, this does **not** error the way a CV fold
would — it silently applies full-data propensity scores, indexed by original row,
to resampled rows.

This is the same class of defect as the cached `grf_res` / `dina_res` that the
code deliberately nulls, and the same class as `mr_inference`.

**It does not affect this analysis** — GBSG is an RCT, `is.RCT = TRUE`, no PS
adjustment, `ps_method` and `ps_hat` both `NULL`. Report it as a latent issue for
observational / IPTW use, not as a finding about this run. Do not fix it here.

### Open question, unverified: `seedit`

`seedit` is inherited, so every replicate runs `forestsearch()` with the same
internal seed (8316951 here). Whether that is intended — holding the internal
split structure fixed across replicates — or whether it suppresses variability
the bootstrap is meant to capture, has **not** been checked. Read what `seedit`
governs and state which it is. Do not change it.

### Not an issue

`parallel_args` looks inherited but is overridden by sub-field
(`$plan`, `$workers`, `$show_message`). No action.

## 7. Deliverable

`dev/replication-check/REPLICATION_FINDINGS.md` — a standalone summary document,
readable without this brief in hand:

1. **Verdict up front.** One paragraph: is `mr_in_replicates` RNG-neutral, do
   the selected subgroups replicate, and is anything blocking the merge into
   `feature/glm-extension`?
2. **Provenance.** Commit SHA, `packageVersion("forestsearch")`, render date,
   host, worker count, guard verdict, and `git status` confined to
   `dev/replication-check/`.
3. **What was compared.** The two baselines, what each isolates, and the
   three-line source difference between them.
4. **Comparison 1 — new vs baseline B** (`v2_2new.html`, rename isolated).
   Full table: per analysis, per quantity, old / new / verdict.
5. **Comparison 2 — new vs baseline A** (`v2_2A_linux.html`, legacy as provided).
   Same table shape.
6. **Timing table.** Search, bootstrap, MR, and MR-vs-bootstrap speed-up, all
   three runs.
7. **Attribution.** Every difference assigned a cause, unexplained items first.
8. **Code assessment.** The `ps_hat` gap, the `seedit` question resolved, and
   any other inherited-argument concern found. For each: what it is, whether it
   affects this analysis, and whether it needs action before or after the merge.
9. **Merge recommendation.** Blocking or not, and anything the maintainer must
   decide before merging.

## 8. Out of scope

* Fixing anything — including the `ps_hat` gap. Report, do not patch.
* Pinning or annotating arguments.
* Simulations, `gbsg_redux`, RNG-neutrality investigation.
* Any other GBSG notebook.
* Any judgement about which `sg_focus` the analysis should use.
