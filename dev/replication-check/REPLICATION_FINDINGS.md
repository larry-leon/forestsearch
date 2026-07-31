# GBSG replication check — findings

Re-run of `analysis_gbsg_cox_multimethod_psi_v2_2new.qmd` on
`feature/mr-in-replicates`, compared against the two existing rendered
baselines, plus an assessment of what the comparison exposes about the code.

> **Post-hoc path update (filename migration).** The harness paths in `R/01_precondition_rng.R` and `R/04_compare.R` were repointed to `analysis_gbsg_survival_multimethod.qmd` and `gbsg_survival_multimethod_payload.rds` after these findings were recorded; that notebook is *not* the file that produced them — it additionally carries parameter-provenance tables, fit-integrity guards, and a corrected `stop_threshold` comment, all intended to be inert.

---

## 1. Verdict

**`mr_in_replicates` is RNG-neutral, every selected subgroup replicates
exactly, and nothing blocks the merge into `feature/glm-extension`.**

All 66 cells of the manuscript table (3 engines x {H, Hc} x 11 quantities —
sizes, naive effects and intervals, bootstrap bias-corrected (FB) effects and
intervals, MR de-biased effects and intervals) are **identical to full double
precision against both baselines**. All five selected subgroups match baseline
B; against baseline A only analysis (A) differs, exactly as the brief
predicted. The MR harm-confirmation flag is `TRUE` in all three engines in all
three runs, with identical de-biased estimates and selection-bias values. No
difference anywhere is unexplained.

The one substantive change is **speed**: the main bootstrap fell from 653 s to
289 s (**2.26x faster**) because MR no longer runs inside the 1000 replicates.
That is the entire point of commit `5bad3df7`, and it is achieved with zero
movement in any estimate.

Two things the maintainer should know, neither blocking:

1. The brief's stated reason for expecting RNG-neutrality is not the actual
   mechanism, and the actual mechanism is more fragile than the substream
   argument suggests (§9.1).
2. The `ps_hat` gap is real and confirmed (§9.2). It cannot affect this
   analysis, but it will silently corrupt any observational/IPTW bootstrap.

---

## 2. Provenance

| item | value |
|---|---|
| commit | `5bad3df710a92a68d663691b81ef8c71b92ea496` (`feat: add mr_in_replicates to the three resampling consumers`) |
| branch | `feature/mr-in-replicates`, run from worktree `.claude/worktrees/replication-check` |
| `packageVersion("forestsearch")` | 0.2.0 (installed via `devtools::install()`, not `load_all()`) |
| render date | 2026-07-30 |
| host | `pop-os`, AMD Ryzen Threadripper PRO 5995WX (64C/128T), 251 GB RAM, Pop!_OS 24.04 |
| R / BLAS | R 4.6.1 (2026-06-24); `libblas.so.3.12.0` / `liblapack.so.3.12.0` |
| workers | 102 (`floor(0.80 * 128)`) — same as both baselines |
| Quarto | 1.9.38 |
| guard verdict | *(see §2.1)* |
| `git status` | only `?? dev/replication-check/` — no tracked file added, modified, or deleted |

**Both baselines ran on this same physical host**, confirmed from
`payload$meta$machine`: identical hostname, CPU model, RAM, kernel, R version,
BLAS and LAPACK, and `cores_used = 102`. The brief inferred this from timing
similarity; it is now established from the recorded fingerprints, which makes
the wall-clock comparison in §7 a like-for-like measurement rather than a
suggestive one.

### 2.1 Read-only guard

Package sources were hashed with `fs_hash_sources()` before and after the run
and compared with `fs_guard_verify()`.

```
$ok      TRUE
$n_files 134
$note    "134 source files verified unchanged"
```

**Verdict: PASS.** No file in the installed package directory or in `R/`
changed across the entire exercise.

The guard baseline was taken **after** `devtools::install()`, because
`fs_hash_sources()` hashes the installed package directory as well as `R/`;
installing first and hashing second is what makes the snapshot describe the
state the run actually used.

No package code, `.qmd`, or test was modified. Everything written by this
exercise lives under `dev/replication-check/`.

---

## 3. What was compared

Two rendered baselines exist. They differ by exactly two hunks of source —
verified by diffing the two faithful reconstructions supplied with the brief:

```
990a991
>   consistency_method = fs_consistency_method,  # ALIGN: match main run + ACTG175
2141c2142
< .payload_file <- ... "gbsg_table2Alinux_payload.rds"
---
> .payload_file <- ... "gbsg_table2new_payload.rds"
```

The added line sits inside the **(A) `fs_dina_screen` call only**, so
`consistency_method` cannot reach any other analysis by construction.

| baseline | what it isolates when compared to the new run |
|---|---|
| **B** — `..._v2_2new.html` (primary) | the gate→MR rename + `mr_in_replicates`, with `consistency_method` held constant |
| **A** — `..._v2_2A_linux.html` | the above **plus** `consistency_method` on (A) |

**Both baselines confirmed pre-rename**, not assumed. Token counts in the
rendered HTML: `debias_gate` 35, `run_debias_gate` 11, `gate_draws` 11,
`harm_flag_debiased` 3, `t_gate` 3; and `mr_inference` / `mr_draws` /
`mr_harm_confirmed` / `t_confirm` **all zero**. The new render inverts this
exactly: all legacy tokens 0, `mr_estimates_table` 3, `gate_estimates_table` 0.

**Sources of truth.** The new side is read from the payload
(`P$table`, `P$meta$timings`), never scraped. Three quantities the payload does
not carry were taken from the rendered HTML console output, in all three runs
identically: analysis (A) (`fs_dina_screen` is never written to the payload),
the harm-confirmation flag, and the search / bootstrap wall-clock lines. The
baseline payloads are git-tracked and contemporaneous with their HTML
(`built_at` 2026-07-08 13:58:55 and 08:18:06, matching the HTML mtimes), so
they are genuine baseline artifacts rather than reconstructions.

### 3.1 The notebook was rendered twice

Two independent renders of the same commit were produced, which turns out to be
useful for separating a real effect from scheduling noise:

| run | where rendered | why |
|---|---|---|
| **run 1** | `dev/replication-check/render/` (copy of the `.qmd`) | first attempt, kept the baseline HTML out of harm's way |
| **run 2** | `quarto/applications/gbsg/` (the brief's §4 command, in place) | canonical location, inherits the project config |

Run 1 sits **outside** the `quarto/` project root, so it does not inherit
`quarto/_quarto.yml` — a project-level format block that adds `code-tools`,
TOC/banner chrome and fonts. That file contains **no `execute:` options**, so it
cannot affect computation, and the results bear that out: run 1 and run 2 agree
on every statistical quantity, and both agree with baseline B on all 66 table
cells. The only visible consequence is that run 1's HTML lacks the "View
Source" button.

Run 2 is the faithful replicate and is what §5–§7 report. Where the two runs
differ at all it is wall-clock only (main bootstrap 290.5 s vs 289.0 s; FS
search 26.5 s vs 26.3 s — under 1%), and that agreement between independent
renders is itself evidence that the timing change in §7 is a real effect rather
than one machine-load artifact.

**Configuration confirmed unchanged before rendering** (§3 of the brief):
`sg_focus = "effMaxSG"`, `hr.threshold = 1.0`, `hr.consistency = 1.0`,
`max_subgroups_search = 50`, `conf_force = c("er <= 0", "pgr <= 0")`,
`consistency_method = "resample"`, `subgroup_method = "consistency"`,
`NB = 1000`, `mr_draws = 5000L`, `run_cv = FALSE`, `run_loo = FALSE`.
CV and LOO reported as skipped in all three runs — there are no CV/LOO metrics
to compare, as expected.

---

## 4. Precondition — is `mr_in_replicates` RNG-neutral?

**Yes.** With the outer seed fixed, a 5-replicate bootstrap was run twice from
the same parent fit, once with `mr_in_replicates = TRUE` and once `FALSE`:

| quantity | result |
|---|---|
| `Ystar_mat` | **bit-identical** (`identical()`) |
| per-replicate results table, all 32 statistical columns | **bit-identical** |
| — including `g_sg`, `hr_sg`, `N_sg`, `E_sg`, `K_sg`, `m_sg` (per-replicate selected subgroup) | identical |
| — including `H_biasadj_1/2`, `Hc_biasadj_1/2`, `max_sg_est`, `Pcons`, `M.1`–`M.7` | identical |
| `H_estimates`, `Hc_estimates`, `SG_CIs`, `FSsg_tab` | identical |
| `tmins_search`, `tmins_iteration` | **differ** (wall-clock only: 8.6–9.3 s vs 23.2–26.7 s per replicate) |
| `mr_replicates` | `NULL` under FALSE; 5 retained objects under TRUE |

The only columns that move are the two timing columns. **Exact agreement was
therefore the expectation for every quantity in §5–§6, and that is what was
observed.**

The mechanism differs from the brief's, in a way worth recording — see §9.1.

---

## 5. Comparison 1 — new run vs baseline B (rename isolated)

### 5.1 Selected subgroups — the primary criterion

| analysis | new | baseline B | verdict | size new / B |
|---|---|---|---|---|
| FS main (`fs`) | `{er <= 0} & {size <= 35}` | same | **same** | 61 (8.9%) / 61 (8.9%) |
| (A) DINA-screened (`fs_dina_screen`) | `{er <= 0} & {size <= 35}` | same | **same** | 61 (8.9%) / 61 (8.9%) |
| (C) DINA-selected (`fs_dina`) | `{grade3 >= 1} & {pgr <= 10}` | same | **same** | 89 (13.0%) / 89 (13.0%) |
| (D) GRF-selected (`fs_grf`) | `{er <= 0}` | same | **same** | 82 (12.0%) / 82 (12.0%) |
| (B) standalone DINA frontier | `grade3 >= 1.0`, `pgr >= 32.5` | same | **same** | 127 / 127 |

Percentages are `n / 686`. The payload stores `pct` rounded to whole numbers
(`9`, `13`, `12`); the values above are the unrounded quotients, so the FS row
reads 8.9% here and `9` in the payload. The stored `pct` values are identical
across all three runs — the rounding is a display choice, not a discrepancy.

The full candidate-evaluation summaries also match line for line, including
the ranking, `Pcons`, Pareto-frontier and in-band flags for every candidate:

* FS main: `Evaluated: 50  Passed: 4  On frontier: 3  In band: 2  Selected: m=1`
* (A): `Evaluated: 30  Passed: 11  On frontier: 9  In band: 2  Selected: m=1`

both identical to baseline B.

### 5.2 Manuscript table — all 66 cells

`n`, `pct`, `naive_est/lo/hi`, `fb_est/lo/hi`, `mr_est/lo/hi` for each of
FS/H, FS/Hc, DINA/H, DINA/Hc, GRF/H, GRF/Hc.

**Differing cells vs baseline B: 0 of 66.** Values agree to full stored
precision. Representative rows:

| row | quantity | new | baseline B |
|---|---|---|---|
| FS/H | naive | 2.5369176 (1.2453976, 5.1677879) | identical |
| FS/H | FB (bootstrap) | 1.9353002 (0.8895315, 4.2105164) | identical |
| FS/H | MR de-biased | 1.4485071 (0.6220092, 3.3732186) | identical |
| FS/Hc | FB | 0.6304172 (0.4271978, 0.9303087) | identical |
| DINA/H | FB | 0.9154569 (0.3659104, 2.2903460) | identical |
| GRF/H | FB | 1.0439733 (0.4829750, 2.2565975) | identical |

That the **FB** column is bit-identical is the strongest single piece of
evidence: it is the output of all 1000 bootstrap replicates plus the two-source
bias correction and the IJ variance. It could not match if `mr_in_replicates`
had perturbed the replicate RNG.

### 5.3 MR de-biasing and harm confirmation

| engine | de-biased HR | selection bias (log-HR) | candidate family | harm confirmed | verdict |
|---|---|---|---|---|---|
| FS | 1.449 | 0.5556458 | 1744 | TRUE | **same** |
| DINA | 1.225 | 0.2632054 | 1 | TRUE | **same** |
| GRF | 1.624 | 0.1770723 | 1 | TRUE | **same** |

Identical in all three runs. The console wording changed with the rename
(`De-biased gate: HR = 1.449 ... -> consistent with harm` became
`MR harm confirmation: HR = 1.449 (point rule vs 1.00) -> consistent with harm`;
`Gate harm_flag_debiased = TRUE` became `Gate mr_harm_confirmed = TRUE`), with
the same numbers.

---

## 6. Comparison 2 — new run vs baseline A (legacy as provided)

### 6.1 Selected subgroups

| analysis | new | baseline A | verdict |
|---|---|---|---|
| FS main | `{er <= 0} & {size <= 35}` | `{er <= 0} & {size <= 35}` | **same** |
| **(A) DINA-screened** | `{er <= 0} & {size <= 35}`, 61 (8.9%) | **`{er <= 0} & {nodes <= 7}`, 61 (8.9%)** | **differs** |
| (C) DINA-selected | `{grade3 >= 1} & {pgr <= 10}` | same | **same** |
| (D) GRF-selected | `{er <= 0}` | same | **same** |
| (B) standalone DINA | `grade3 >= 1.0`, `pgr >= 32.5` | same | **same** |

**Exactly one analysis differs, and it is the predicted one.** No second
analysis moved — `consistency_method` did not reach further than the (A) call
it was added to.

### 6.2 Manuscript table

**Differing cells vs baseline A: 0 of 66.** (A) is not represented in the
payload, so the table is unaffected by the one difference above.

### 6.3 Mechanism of the (A) difference

The (A) candidate tables show precisely why it moves. Two candidates lead, tied
at `N = 61` and both inside the effect-size band; the selection turns on
`Pcons`, which is exactly what `consistency_method` estimates:

**Baseline A — `consistency_method = "split"` (the default; line absent):**

| rank | effect | N | E | Pcons | on frontier | in band | selected | subgroup |
|---|---|---|---|---|---|---|---|---|
| 1 | 2.335 | 61 | 31 | **1.00** | - | * | **\*** | `{er <= 0} & {nodes <= 7}` |
| 2 | 2.537 | 61 | 34 | 0.99 | * | * | - | `{er <= 0} & {size <= 35}` |

`Evaluated: 30  Passed: 9  On frontier: 8  In band: 2  Selected: m=2`

**Baseline B and the new run — `consistency_method = "resample"`:**

| rank | effect | N | E | Pcons | on frontier | in band | selected | subgroup |
|---|---|---|---|---|---|---|---|---|
| 1 | 2.537 | 61 | 34 | **1.00** | * | * | **\*** | `{er <= 0} & {size <= 35}` |
| 2 | 2.335 | 61 | 31 | 0.99 | - | * | - | `{er <= 0} & {nodes <= 7}` |

`Evaluated: 30  Passed: 11  On frontier: 9  In band: 2  Selected: m=1`

The two `Pcons` values swap (1.00/0.99 → 0.99/1.00) when the consistency
estimator changes, flipping the ranking between two otherwise tied candidates
and with it the selection. The candidate pool is the same size (30 evaluated)
in both; the pass count differs (9 vs 11) because `Pcons` is what the screen
tests. This is the `ALIGN: match main run` comment doing what it says: under
`"resample"`, (A) lands on the same subgroup as the main FS run.

---

## 7. Timing

All three runs on the same host with 102 workers.

| quantity | new (run 2, canonical) | new (run 1) | baseline B | baseline A |
|---|---|---|---|---|
| FS search (console) | 26.3 s | 26.5 s | 25.2 s | 25.3 s |
| **Main bootstrap (1000 boots)** | **289.0 s (4.8 min)** | **290.5 s (4.8 min)** | **653.4 s (10.9 min)** | **656.0 s (10.9 min)** |
| MR, FS engine | 4.51 s | 4.60 s | 4.44 s | 4.48 s |
| MR-vs-bootstrap speed-up, FS | **64x** | 63x | 147x | 147x |
| DINA-mode bootstrap | 154.9 s | 156.9 s | 159.0 s | 159.1 s |
| MR, DINA engine | 0.097 s | 0.107 s | 0.100 s | 0.093 s |
| MR-vs-bootstrap speed-up, DINA | 1597x | 1467x | 1590x | 1711x |
| GRF-mode bootstrap | 192.8 s | 194.2 s | 194.7 s | 195.1 s |
| MR, GRF engine | 0.089 s | 0.096 s | 0.096 s | 0.088 s |
| MR-vs-bootstrap speed-up, GRF | 2166x | 2023x | 2029x | 2217x |

Reading this table:

* **The main bootstrap is 2.26x faster** (653 s → 289 s). This is the
  `mr_in_replicates` change and nothing else: with `mr_inference = TRUE` on the
  parent fit, the pre-commit code ran a 5000-draw multiplier bootstrap inside
  every one of the 1000 replicates. The precondition measured the per-replicate
  cost directly — search time per replicate went 8.6–9.3 s → 23.2–26.7 s when
  MR was switched on, a ~2.7x factor consistent with the 2.26x observed at
  scale.
* **The DINA- and GRF-mode bootstraps barely move** (159→157 s, 195→194 s),
  which is the same story seen from the other side. For those engines the
  candidate family collapses to `n_family = 1`, so MR costs ~0.1 s rather than
  the FS engine's 4.6 s over a 1744-member family. There was almost nothing to
  remove.
* **The MR-vs-bootstrap speed-up ratios fall** (147x → 64x for FS). This is
  arithmetic, not a regression: the ratio is `boot_sec / mr_sec`, and the
  numerator shrank while MR itself is unchanged (4.44 s → 4.51 s, noise). The
  headline claim the ratio supports — MR approximates the bootstrap at a small
  fraction of the cost — is undamaged; the honest comparison is now against a
  bootstrap that is no longer doing redundant work.
* FS search time (26.3 s vs 25.2 s) and all MR times differ by a few percent —
  ordinary run-to-run wall-clock variation on a shared 128-thread host.

---

## 8. Attribution of every difference

Nothing landed in **unexplained**.

| difference | cause | note |
|---|---|---|
| (A) selects `{er <= 0} & {size <= 35}` vs A's `{er <= 0} & {nodes <= 7}` | **consistency_method (A)** | Predicted by the brief; mechanism established in §6.3. Not present vs baseline B. |
| Main bootstrap 653 s → 289 s; FS speed-up 147x → 64x | **`mr_in_replicates` (intended effect)** | Not an RNG shift — see §4. Pure removal of redundant computation. |
| Console vocabulary: `Gate harm_flag_debiased` → `Gate mr_harm_confirmed`; `De-biased gate: HR =` → `MR harm confirmation: HR =`; `gate_estimates_table` → `mr_estimates_table` | **rename** | Inert: every number attached to these labels is identical. |
| Payload field `meta$gate_draws` → `meta$mr_draws`; `timings$*$gate_sec` → `mr_sec` | **rename** | Same quantity, same values (`draws` 5000 in all three runs). No in-repo consumer — see §10. |
| Small wall-clock jitter in search and MR times (±5%) | **parallel scheduling** | Same host, 102 workers, shared machine. |

There is no residual category. In particular, the rename is **inert on
results**, which is what the brief asked to be checked rather than assumed: had
the rename silently changed a code path, the FB column could not have matched
to full precision.

---

## 9. Code assessment — inherited replicate arguments

`forestsearch()` has **85 formals**. `bootstrap_analysis_dofuture.R` overrides
**12** per replicate (`df.analysis`, `df.predict`, `grf_res`, `grf_cuts`,
`dina_res`, `dina_cuts`, `details`, `quiet`, `show_candidate_summary`,
`plot.sg`, `plot.grf`, and now `mr_inference`); the other **73** are inherited
wholesale from `fs.est$args_call_all`. Counts confirmed directly.

### 9.1 The RNG-neutrality mechanism is not the one the brief assumed

The brief expected neutrality from the substream design, noting that
`fs_mr_inference()` "takes `seed = NULL` and only calls `set.seed()` when given
one, and the notebook supplies none — so MR does draw from whatever stream it
lands in."

That is true of `fs_mr_inference()` in isolation but not of how it is reached.
`forestsearch()` passes

```r
seed = .g_mr(mr_inference_args$seed, seedit)      # forestsearch_main.R:3052
```

so with no `seed` in `mr_inference_args` — which is this notebook's case — MR
receives `seedit` and **does** call `set.seed(8316951)`. MR never draws from
the ambient stream at all.

The consequence is a latent fragility rather than a defect. `set.seed()` fired
inside a `%dofuture%` worker **discards that worker's L'Ecuyer-CMRG substream**
and replaces it with a fixed Mersenne-Twister state. It is harmless today only
because MR is the last stochastic step in `forestsearch()` — nothing but output
assembly (SECTION 10) follows it. The empirical check in §4 confirms this
holds now. But the protection is ordering, not the substream architecture the
brief credited: any future stochastic step added after the MR block, inside a
replicate, would silently draw from a **constant** seed in every replicate.

The two structural guarantees the brief cited *are* real and were verified:
`boot_index_mat` is built up front by `.make_boot_index_matrix(seed = boot_seed)`
in the driver and `Ystar_mat` derived from it deterministically, so resample
selection cannot be perturbed by anything a replicate does; and the loop does
set `.options.future = list(seed = TRUE)`.

*Action: none required for the merge.* Worth a comment at the MR call site
noting that it re-seeds and must remain last.

### 9.2 Confirmed gap: `ps_hat` (latent, does not affect this analysis)

Confirmed exactly as described. `bootstrap_analysis_dofuture.R` and
`bootstrap_dofuture_main.R` contain **zero** references to `ps_hat`,
`ps_method`, `ps_adjust_method`, or propensity scores. Both CV paths null it
explicitly (`forestsearch_cross_validation.R:396` and `:859`), with the comment
*"PS must be re-estimated on each training fold — NULL out any user-supplied
ps_hat (wrong length for training fold)"*, and the CV summary reports
`ps_hat: re-estimated on each fold`.

The failure is silent because the only guard is a length check
(`forestsearch_main.R:1828`):

```r
if (length(ps_hat) != nrow(df.analysis))
  stop("ps_hat must have length equal to nrow(df.analysis)")
...
df.analysis$ps_hat <- ps_hat
df.analysis$sw   <- ifelse(W == 1, p_treat / ps_hat, (1 - p_treat) / (1 - ps_hat))
df.analysis$ips_covar <- ifelse(W == 1, 1 / ps_hat, 1 / (1 - ps_hat))
```

A CV fold has fewer rows, so CV would error — which is why it was noticed
there. A bootstrap resample has **the same row count**, so the check passes and
the assignment binds **positionally**: replicate row *i* receives
`ps_hat[i]`, but replicate row *i* is original row `in_boot[i]`. The score is
therefore attached to the wrong subject for every row where
`in_boot[i] != i` — essentially all of them — and the IPTW weights `sw` and
`ips_covar` are computed from the mismatched values. Note also that
`ps_method` is mirrored into `args_call_all` (`forestsearch_main.R:1801`), so
the replicate takes the user-supplied-PS branch rather than re-estimating.

This is the same defect class as the cached `grf_res` / `dina_res` the code
already nulls deliberately, and the same class as `mr_inference`.

**Does not affect this analysis.** GBSG is an RCT: `is.RCT = TRUE`, and
`ps_method` and `ps_hat` are both `NULL` (confirmed in the notebook call and by
the formal defaults), so the whole branch is inert. Reported as a latent issue
for observational / IPTW use. **Not fixed here**, per the brief.

*Action: fix before any observational/IPTW use of the bootstrap; not blocking
this merge.* The one-line fix mirrors CV — null `ps_hat` in the replicate args
and let `ps_method` re-estimate per resample.

### 9.3 Open question resolved: `seedit`

**`seedit` does not govern the consistency search's split structure at all.**
The brief's framing — "every replicate runs `forestsearch()` with the same
internal seed … holding the internal split structure fixed" — attributes to
`seedit` an influence it does not have.

`forestsearch()` has **no `seed` formal** (it has `seedit`).
`subgroup.consistency()` **does**, with its own hard default `seed = 8316951`.
The consistency stage is invoked as

```r
sc_filtered_args <- filter_call_args(args_call_all, subgroup.consistency,
                                     consistency_overrides)
```

and `filter_call_args()` matches **by exact name**
(`summary_utility_functions.R:679-687`). `seedit` and `seed` are different
names, `consistency_overrides` contains no `seed`, and `.sync_args_call_all()`
only ever mirrors names that are already `forestsearch()` formals
(`stop_threshold`, `dmin.grf`, `ps_method`, `adjust_covariates`,
`event.name`) — so no `seed` key is ever injected. Verified against the
installed package: `"seed" %in% names(formals(forestsearch))` is `FALSE`, and
`seed` is absent from the formals intersection.

So `subgroup.consistency()` runs on its **own default seed, 8316951, in the
parent fit and in every replicate, whatever `seedit` is set to**. That the two
constants coincide (`seedit`'s default is also 8316951) makes this easy to miss
and is presumably why it read as inheritance.

**Which is it — intended, or suppressed variability?** It is the benign case,
and it would remain benign even if routed through `seedit`. The consistency
search's RNG is applied to *the replicate's own resampled data*, so the
bootstrap's genuine sampling variability is fully captured; what is held fixed
across replicates is only the Monte-Carlo noise of the internal splitting
procedure. That is a common-random-numbers variance reduction: it makes the
bias estimate less noisy rather than optimistically stable. Holding it fixed is
the defensible choice.

The caveat is discoverability, not correctness: because `seedit` does not reach
the consistency search, a user who changes `seedit` hoping to vary the internal
split structure will observe no change, and there is no message saying so.

*Action: none required. Worth documenting that `seedit` does not control the
consistency-search seed.*

### 9.4 Other inherited arguments — surveyed, no further defects

The risky class is anything carrying a **pre-computed full-data object** or a
**post-selection layer**. All 73 inherited formals were classified:

* **`df.test`** — the only other inherited full-data object. Unlike `ps_hat` it
  is **not row-aligned** to `df.analysis` (it is an independent cohort scored
  through `get_dfpred()`), so there is no index-misalignment failure. Each
  replicate re-scores it and the result lands in `df.test_out`, which the
  bootstrap never reads. Cost is wasted compute when `df.test` is supplied, not
  a wrong number. `NULL` here.
* **`conf.cont_medians` / `conf.cont_medians_force`** — these look like the
  risky class but are **character vectors of confounder names**, not
  pre-computed median values (`get_fsdata.R:289-314`); the cut points are
  recomputed on each resample. Not a risk.
* **`conf.cont_jcuts`** (`list(er = 10, pgr = 10)` here) — a count of cuts, not
  values. Not a risk.
* **`parallel_args`** — as the brief notes, overridden by sub-field
  (`$plan`, `$workers`, `$show_message`). Confirmed, no action.
* **`max.minutes`** (3 here) — inherited as a per-replicate budget, so a slow
  replicate can return a truncated search. Behaviour is by design, not an
  inheritance defect; flagged only because it is a soft failure mode that would
  not announce itself.
* The remaining ~66 are algorithm configuration — thresholds, `sg_focus`,
  `maxk`, screening flags, variable names, GLM/GRF/DINA tuning — which *should*
  be inherited, since the bootstrap's purpose is to re-run the same algorithm
  on a resample.

No further defect found.

---

## 10. Merge recommendation

**Not blocking. Recommend merging `feature/mr-in-replicates` into
`feature/glm-extension`.**

The evidence: `mr_in_replicates` is RNG-neutral by direct measurement; all five
selected subgroups and all 66 manuscript-table cells reproduce the primary
baseline exactly; the sole difference against the legacy baseline is the
already-characterised `consistency_method` effect on (A), confined to (A); and
the change delivers a 2.26x faster main bootstrap with no movement in any
estimate.

Decisions for the maintainer, none of which gate the merge:

1. **Payload field rename — check any reader living outside this repo.**
   `meta$gate_draws` → `meta$mr_draws` and `timings$*$gate_sec` → `mr_sec`.
   **No in-repo consumer exists**: the only `readRDS(...payload...)` calls
   anywhere in the tree are this exercise's own comparison scripts, so nothing
   inside the repo breaks. The payload is written as an export (the notebook
   comment anticipates *"your gbsg reader/manuscript doc"*), so if such a
   document is maintained elsewhere it needs the same field update — and it
   would fail silently as `NULL` rather than erroring. The writers are already
   internally consistent: the gbsg and ACTG175 notebooks all emit the new
   names.
2. **`ps_hat` (§9.2)** — schedule the one-line fix before any
   observational/IPTW bootstrap use. Harmless for RCT analyses today.
3. **MR re-seeds inside the worker (§9.1)** — consider a comment pinning the
   invariant that MR must remain the last stochastic step inside
   `forestsearch()`.
4. **`seedit` does not reach the consistency search (§9.3)** — a documentation
   fix, not a code change, unless routing `seedit` through to
   `subgroup.consistency(seed=)` is actually wanted. If it is, note that doing
   so would change results for any user whose `seedit` differs from 8316951.

Also worth deciding separately: baseline A's `{er <= 0} & {nodes <= 7}` for (A)
is now unreachable with `consistency_method = "resample"` in the source. If the
legacy (A) result matters for any manuscript comparison, it needs to be
regenerated deliberately rather than recovered later.

---

## 11. Files

Everything under `dev/replication-check/`:

| path | content |
|---|---|
| `REPLICATION_FINDINGS.md` | this document |
| `R/00_guard_capture.R` | guard hashing + environment capture |
| `R/01_precondition_rng.R` | the §4 RNG-neutrality test |
| `R/02_extract_html.R`, `R/03_extract2.R` | HTML → comparable quantities |
| `R/04_compare.R` | the two comparisons |
| `out/precond_rng.rds`, `out/precond_rng.log` | precondition results |
| `out/comparison_run2.rds`, `out/comparison_run2.log` | **canonical comparison (run 2) — the one §5–§7 report** |
| `out/comparison.rds`, `out/comparison.log` | run-1 comparison |
| `out/guard_before.rds`, `out/guard_after.rds`, `out/guard_verdict.rds` | read-only guard |
| `out/env_capture.rds` | environment fingerprint |
| `render_canonical/` | run 2 HTML + payload (rendered in place, then the worktree's tracked copies restored with `git checkout`) |
| `render/` | run 1 HTML + payload + log |
| `out/render2.log` | run 2 render log |

`R/04_compare.R` takes the new-run directory and an output tag as arguments, so
either run can be re-compared:

```bash
Rscript dev/replication-check/R/04_compare.R quarto/applications/gbsg _run2
Rscript dev/replication-check/R/04_compare.R dev/replication-check/render
```

Nothing outside `dev/replication-check/` was left modified: the canonical render
necessarily overwrote the worktree's own copies of
`analysis_gbsg_cox_multimethod_psi_v2_2new.html` and
`gbsg_table2new_payload.rds`, and both were restored from git after the outputs
were preserved above. The baselines compared against are the **main repository's**
copies, which were never touched.
