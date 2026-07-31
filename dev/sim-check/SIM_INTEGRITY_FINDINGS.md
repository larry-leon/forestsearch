# SIM INTEGRITY FINDINGS

Response to `dev/sim-check/CC_BRIEF_sim_integrity.md`.
Target: `quarto/simulations/gbsg_redux/sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd`
Comparator: `quarto/simulations/gbsg_redux/fs_t1_t2_m1_h10_knoise0_n500_batch_1_20.qmd`

Nothing outside `dev/sim-check/` was modified. Both `.qmd` files under
`quarto/simulations/gbsg_redux/` are byte-for-byte as found.

Run environment: R 4.6.1, `forestsearch` 0.2.0 (installed library, not `load_all`),
128 physical cores → `n_workers = 115`.

---

## 1. Verdict

**The edits are inert.** All five replicates produced **bit-identical**
per-replicate records across every non-wall-clock column (50 of 52 columns
identical at `tolerance = 0`; the two that differ are `fb_secs`/`fit_mr_secs`,
which the brief exempts). The column sets match exactly under the stated
mapping — no dangling name in either direction. Every pinned parameter reaches
the engine at its expected value, no reset warning is emitted on either route,
and `mr_in_replicates = FALSE` vs `TRUE` leaves `H_estimates`/`Hc_estimates`
identical.

**Nothing is functionally wrong.** Four things are wrong *as documented*, and
one consequence of the stem change is larger than the comment implies:

| # | Issue | Severity |
|---|---|---|
| A | The `mr_in_replicates = FALSE` comment's rationale is factually false — the package default is already `FALSE`, so the pin bought nothing | comment only |
| B | `focus_tag` normalizes `"eff"` but not `"hr"`, so the *same* rule under its canonical spelling writes a different stem | latent |
| C | `meta` still records no `sg_focus`; combine-mode's poolability check does not consult it, so the stem is the sole guard and `combine_files` bypasses it | latent |
| D | The stem change orphans a **500-replicate** archived pool for this exact cell, not just "old batches" | scope |
| E | `fit_mr_secs` is fit **+** MR but every label reads "MR", so the headline MR-vs-FB cost comparison overstates MR's cost | pre-existing, now more visible |

Details in §6. None of these changes any number produced by the current file.

---

## 2. §1 — The decisive test: bit-identical records

### Setup

Both files copied to `dev/sim-check/run_pre/pre.qmd` and
`dev/sim-check/run_post/post.qmd`, with `betaHhat_truth.R` alongside each.
Per the brief, in each copy: `n_sims <- 5L`, `sim_id_start <- 1L`,
`nb_boots <- 20L`, `mr_draws <- 500L`. The pre-edit file already carried all
four of those values, so **only** `n_sims` (50L → 5L) was touched in the
post-edit copy. Nothing else was changed in either copy.

Both rendered to completion. Both report `3/5 detections`; wall-clocks 125.9 s
(pre) and 123.0 s (post).

### Column-set mapping

Applied in this order (specific before generic, so `t2_secs` does not become
`mr_secs`): `t2_secs → fit_mr_secs`, then `t1_ → fb_`, then `t2_ → mr_`.

| | result |
|---|---|
| pre columns | 52 |
| post columns | 52 |
| in mapped-PRE but not POST | **(none)** |
| in POST but not mapped-PRE | **(none)** |
| column *order* identical | **TRUE** |

No dangling name — the rename is complete.

### Per-column comparison, `tolerance = 0`, matched on `sim_id`

`sim_id` 1–5 present in both, same order.

| | count |
|---|---|
| columns compared | 52 |
| identical | **50** |
| differing — wall-clock (exempt) | 2 |
| differing — **substantive** | **0** |
| class mismatches | 0 |

The two exempt columns:

```
fb_secs      pre:  17.34, 16.47, 16.49, NA, NA
             post: 17.11, 15.92, 16.46, NA, NA
fit_mr_secs  pre:  18.35, 15.80, 18.05, 11.21, 12.01
             post: 17.96, 15.19, 17.92, 10.78, 11.48
```

Everything else — `sg_def`, `covs`, `label`, `detected`, `mr_ok`, `n_sel`,
`n_harm`, `n_true`, `betaHhat_H`/`betaHhat_Hc`, `fb_err`, all four estimator
blocks (`or_`/`nv_`/`fb_`/`mr_` × `H`/`Hc` × `est`/`lo`/`hi`/`se`),
`mr_harm_flag`, `ij_source`, and `sens`/`spec`/`ppv`/`npv` — is identical.
`identical(pre$truth, post$truth)` is `TRUE`; the `meta` bundles agree on every
key.

Reproduce with:

```
Rscript dev/sim-check/compare_rds.R \
  dev/sim-check/run_pre/fs_t1_t2_m1_h10_knoise0_n500_res_1_5.rds \
  dev/sim-check/run_post/fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_5.rds
```

Full output: `dev/sim-check/compare_output.txt`.

---

## 3. §2 — Do the pinned parameters reach the engine?

One replicate (`sim_id = 1`) fitted with the post-edit `base_args`, reading
`fs.est$args_call_all`. Warnings captured with `withCallingHandlers`, not
suppressed by chunk options.

| parameter | expected | **observed** | |
|---|---|---|---|
| `sg_focus` | `"hr"` | `"hr"` | ✅ |
| `subgroup_method` | `"consistency"` | `"consistency"` | ✅ |
| `consistency_method` | `"resample"` | `"resample"` | ✅ |
| `stop_threshold` | `NULL` | `NULL` (name present in the list, not dropped) | ✅ |
| `selection_rule` | `"neighborhood"` | `"neighborhood"` | ✅ |
| `effect_neighborhood` | `0.10` | `0.1` | ✅ |
| `mr_inference` | `TRUE` | `TRUE` | ✅ |

`"maxcons"` → `"hr"` is confirmed in the installed engine
(`forestsearch:::.normalize_sg_focus("maxcons") == "hr"`, identical to
`"eff"`), and `effect_neighborhood`'s pinned `0.10` is the formal default
(`formals(forestsearch)$effect_neighborhood == 0.1`).

### `stop_threshold` — both routes, side by side

The same fit was re-run with `stop_threshold` and `effect_neighborhood`
**removed from the argument list entirely** and `sg_focus = "eff"` — i.e. the
pre-edit call.

| | post-edit (`stop_threshold = NULL` pinned) | pre-edit (omitted) |
|---|---|---|
| `args_call_all$stop_threshold` | `NULL` | `NULL` |
| `args_call_all$sg_focus` | `"hr"` | `"hr"` |
| `args_call_all$effect_neighborhood` | `0.1` | `0.1` |
| warnings emitted | **0** | **0** |
| any mentioning `stop_threshold` | FALSE | FALSE |
| `sg.harm` identical to the other route | — | **TRUE** |
| `mr_inference$debiased` identical | — | **TRUE** |

**The engine saw `NULL` either way, and the render emits no reset warning.**
The two routes reach it differently, which is why both had to be checked:

* *Pre-edit*: `stop_threshold` defaults to `pconsistency.threshold` (0.90).
  `forestsearch_main.R:1206` (`mget(args_names, ...)`) forces that promise, so
  `args_call_all` briefly holds 0.90. The reset block at line 1402 fires
  (`!is.null(stop_threshold)` and `sg_focus == "hr"`), takes the
  `user_explicit <- !missing(stop_threshold)` branch as **FALSE** — hence no
  warning — sets the local to `NULL`, and re-syncs `args_call_all` at line 1435.
* *Post-edit*: `NULL` is supplied, so `!is.null(stop_threshold)` is FALSE and
  the reset block is **never entered**. `args_call_all` holds `NULL` from the
  `mget` onward.

Both land on `stop_threshold = NULL` in `args_call_all` **and** in the value
passed to the engine (`consistency_overrides$stop_threshold`, line 2776, reads
the local). The brief's expectation is confirmed on both counts.

Full output: `dev/sim-check/run_post/diag_probe.log` (script:
`dev/sim-check/diag_probe.R`).

---

## 4. §3 — Does `mr_in_replicates = FALSE` do what the comment claims?

Same fitted `fs.est`, same seed (`seed_base + sim_id`), `nb_boots = 20`, run
once each way.

| | result |
|---|---|
| `H_estimates` identical | **TRUE** |
| `Hc_estimates` identical | **TRUE** |
| `SG_CIs`, `FSsg_tab`, `Ystar_mat`, `original_sg`, `est.scale` identical | **TRUE** |
| per-replicate table (`boot$results`, 20 × 34) | see below |
| names present only under `TRUE` | (none) |
| `$mr_replicates` | `NULL` (FALSE) vs list of 20 (TRUE) |
| wall-clock | **16.8 s** (FALSE) vs **21.7 s** (TRUE) — **1.30×** |

**Per-replicate table (`boot$results`, 20 × 34):** identical on **32 of 34**
columns, including everything substantive — `H_biasadj_1`, `H_biasadj_2`,
`Hc_biasadj_1`, `Hc_biasadj_2`, `max_sg_est`, `Pcons`, `hr_sg`, `N_sg`, `E_sg`,
`K_sg`, `any_found`, all eight `events_*`, `M.1`–`M.7`. The only two that
differ are the wall-clock columns `tmins_search` and `tmins_iteration`
(max |Δ| ≈ 0.094 min ≈ 5.7 s across 20 rows), plus the object's `timing`
attribute. Under `TRUE` the table additionally carries an `mr_replicates`
attribute (this is how `bootstrap_dofuture_main.R:449` returns the MR payload);
that is the sole structural difference.

**Identical estimates plus a real timing gap — the expected outcome.** The pin
is a pure cost saving *if the alternative were `TRUE`*. It is not: see finding
A below — `FALSE` is already the formal default, so the pre-edit call (which
omitted the argument) ran on exactly this path. The 1.30× is the cost you
avoid by not opting in, not a cost the pin removed.

At production settings the gap is much larger than 1.30×: it scales with
`mr_draws` and `nb_boots`, and this probe ran at `draws = 500`, `B = 20`
against a production configuration of `draws = 5000`, `B = 300`.

---

## 5. §4 — Structural checks

| # | check | verdict |
|---|---|---|
| 1 | Setup chunk sources standalone in a clean session | ✅ **PASS** |
| 2 | `rds_stem` matches the filename; `rds_path`/`combine_glob` consistent | ✅ **PASS** |
| 3 | Old batches invisible to the new `combine_glob` | ✅ **PASS** (with a caveat — see D) |
| 4 | `fit_mr_secs` measures fit **+** MR; both labels agree | ⚠️ **PASS with a caveat** |
| 5 | `fb_*` / `mr_*` field paths are real | ✅ **PASS** |
| 6 | No dangling `t1`/`t2`; no write-only or read-only `fb_*`/`mr_*` | ✅ **PASS** |

### 5.1 Symbol resolution

The setup chunk (lines 67–276) was extracted verbatim and sourced with
`Rscript --vanilla` in `dev/sim-check/setupcheck/` with `betaHhat_truth.R`
present. **Ran to completion, no `object ... not found`.** Resolved values (at
the file's own `n_sims = 50L`):

```
rds_stem     = fs_maxcons_fb_mr_m1_h10_knoise0_n500
rds_path     = fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_50.rds
combine_glob = fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_*.rds
focus_tag    = maxcons
method_tag   = fs
n_workers    = 115
```

This is the check that the already-known `focus_tag`/`sg_focus` ordering defect
would have failed. It now passes: `sg_focus` is defined at line 127, ahead of
`focus_tag` (155) and `rds_stem` (157).

### 5.2 Stem matches filename

`rds_stem` = `fs_maxcons_fb_mr_m1_h10_knoise0_n500` — exactly the brief's
expectation, and exactly the stem of the on-disk
`fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_100.rds`. `rds_path` is what the
run writes: the smoke render wrote
`fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_5.rds` at `n_sims = 5`, matching
`sprintf("%s_res_%d_%d.rds", rds_stem, 1, 5)`. `Sys.glob(combine_glob)` in
`quarto/simulations/gbsg_redux/` returns that stem's files and nothing else.

### 5.3 Old batches are invisible

Evaluated in the real directory:

```
glob    : fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_*.rds
matches : fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_100.rds   (1 file)
```

All **42** `fs_t1_t2_m1_*_res_*.rds` files present in that directory are
excluded. Intended, and stated in the file's own comment at lines 141–143.
The scope of what is excluded is understated — see finding D.

### 5.4 `fit_mr_secs` measures fit **plus** MR

`mr_t0 <- proc.time()[3]` is at line 477; the `forestsearch()` call is at 478;
`rec$fit_mr_secs <- proc.time()[3] - mr_t0` is at 485. **No statement
intervenes** — the recorded value is the whole `forestsearch()` call, fit *and*
MR, as the inline comment says ("forestsearch fit INCLUDING MR").

The two printed summaries **agree with each other** — both call it "MR":

* run-loop (631): `Mean wall-clock: MR %.3f s vs FB bootstrap %.2f s (B=%d).`
* timing projection (754 → 822): `mr_mean <- mean(results$fit_mr_secs)`,
  printed as `Per rep: MR %.1f s`, and labelled `mr_rep_s = "MR/rep (s)"` (838).

Rendered: `Mean wall-clock: MR 14.665 s vs FB bootstrap 16.50 s (B=20)`.

The check as posed passes. The caveat is that "MR" is the wrong name for the
quantity in *all three* places — see finding E.

### 5.5 Field paths

Verified against `str()`-equivalent introspection on a real fit and a real
bootstrap object, not against the source:

```
names(boot)              : results, SG_CIs, FSsg_tab, Ystar_mat, H_estimates,
                           Hc_estimates, nb_boots, original_sg, outcome_type,
                           effect_measure, est.scale, mr_replicates
names(boot$H_estimates)  : H0, sdH0, H0_lower, H0_upper, H1, sdH1, H1_lower,
                           H1_upper, H2, sdH2, H2_lower, H2_upper
  H_estimates$H2         = 1.002853
  H_estimates$H2_lower   = 0.35904019
  H_estimates$H2_upper   = 2.8011187
  H_estimates$sdH2       = 0.52557111
  Hc_estimates$H2        = 0.66065951      (all four present on Hc too)

names(fs.est$mr_inference)           : selected_index, selected_label, measure,
  log_scale, ci_method, naive, debiased, selection_bias, fixed_bias,
  selection_rate, complement, settings, harm_flag, n_family, n_selected,
  timing_seconds
names(...$debiased)                  : est, lower, upper, lower_1s, se, se_ij,
  se_wald, var_ij, ij_source, ij_draws
  mr_inference$debiased$est          = 0.76717884
  mr_inference$debiased$lower        = 0.39999624
  mr_inference$debiased$upper        = 1.4714223
  mr_inference$debiased$se_ij        = 0.33228407
  complement$debiased present        : TRUE
```

Every path the file reads exists and is populated. No silent `NA` path.
Consistent with §1, where all `fb_*`/`mr_*` estimate columns are non-`NA` on
the three detected replicates.

### 5.6 No dangling references

**`t1`/`t2` symbols:** one hit in the whole file, at line 141 — inside a
comment documenting the stem change (`"<engine>_t1_t2_m1_..."`). **No code
symbol survived the rename.** (This is what the already-known `t2m` defect
would have shown; it is fixed — the symbol is now `mr_m`, assigned at 1152 and
read at 1153.)

**`fb_*`/`mr_*` assigned-but-never-read or read-but-never-assigned:** none. All
20 `fb_`/`mr_` record columns declared in `.na_record()` (382–398) are both
assigned and read. Most are read through `paste0(k, "_", suffix, "_est")`-style
construction (`build_block`, `build_block_med`, `agree_block`, the
`betaHhat` callout) rather than literally, which is why a naive grep shows a
count of 2 for them; each constructed name resolves to a declared column.
`.get_se()`'s `k == "mr"` branch correctly selects the `_se_ij` suffix for MR
and `_se` for the rest — the rename kept that branch in step with the column
names.

Symbols renamed consistently across every use: `run_tier1 → run_fb`
(including the combine-mode reassignment at 710), `t1_secs → fb_secs`,
`t2_secs → fit_mr_secs`, `t1_err → fb_err`, `t1_tried/t1_fail →
fb_tried/fb_fail`, `t2_mean → mr_mean`, `t1_mean → fb_mean`,
`t2_rep_s/t1_rep_s → mr_rep_s/fb_rep_s`, `t2_harm_flag → mr_harm_flag`,
`t2_t0 → mr_t0`, `t1_t0 → fb_t0`, and the ggplot columns `t1/t2 → fb/mr`
(1214, 1234).

---

## 6. Beyond the brief

### A. The `mr_in_replicates = FALSE` comment's rationale is false

Lines 542–547 say:

> *"Pinned to the package default. `fs.est` is fitted with `mr_inference = TRUE`,
> and the replicate call is a copy of that fit's argument list, so **without this
> the flag would be inherited and MR would run inside all `nb_boots` replicates
> of every sim**."*

The first sentence is right; the bolded claim is not.
`forestsearch_bootstrap_dofuture()`'s formal default **is** `FALSE`
(`R/bootstrap_dofuture_main.R:193`), and `R/bootstrap_analysis_dofuture.R:597`
strips the inherited flag whenever the argument is not `TRUE`:

```r
if (!isTRUE(mr_in_replicates)) args_FS_boot$mr_inference <- FALSE
```

So the pre-edit call, which omitted the argument, already took this path. The
comment describes a hazard the package closes by default; as written it reads
as a cost saving that this edit won. It did not — hence §1's identical `fb_*`
columns. The final clause ("No-op as written") is correct and is the whole
story; the sentence before it should say the default already does this rather
than that the flag "would be inherited".

### B. `focus_tag` normalizes `"eff"` but not `"hr"`

```r
focus_tag <- switch(sg_focus, eff = "maxcons", sg_focus)
```

`"eff"` and `"maxcons"` both fold to `maxcons`, so those two spellings share a
stem — the stated goal. But `"hr"` is a third legal spelling of the *same*
rule (it is the canonical form both aliases normalize to inside
`forestsearch()`, confirmed in §3), and it falls through the `switch` default
to `focus_tag = "hr"`, giving stem `fs_hr_fb_mr_m1_...`.

The comment's claim that the change means "runs under a different engine or
focus never overwrite each other" holds. The converse does not: **runs under
the same focus can now fail to pool**, silently, if the spelling differs.
Setting `sg_focus <- "hr"` and re-rendering produces a second, disjoint stem
whose batches `combine_glob` will never see. Adding `hr = "maxcons"` to the
`switch` closes it — or, more robustly, deriving the tag from
`.normalize_sg_focus(sg_focus)`, which the comment at 152–154 explicitly
rejects for a naming reason that a small forward map would satisfy anyway
(`hr`/`eff`/`maxcons` → `maxcons`, `hrMaxSG`/`effMaxSG` → `effMaxSG`, …).

### C. `meta` records no `sg_focus`, and combine-mode does not check it

The saved bundle's `meta` (655–660) is `n_sample`, `n_sims`, `nb_boots`,
`mr_draws`, `subgroup_method`, `sim_id_start`, `sim_id_end`, `seed_base`,
`parallel_mode` — **no `sg_focus`, no `focus_tag`**. Combine-mode's poolability
loop (681) validates `n_sample`, `nb_boots`, `mr_draws`, `subgroup_method`,
`seed_base` — again not the focus.

So the *only* thing preventing two different foci from being pooled is the
filename stem, and `combine_files` (164) is an explicit escape hatch that
bypasses the glob entirely. Recording `sg_focus` in `meta` and adding it to the
agreement loop would make the guard hold regardless of filenames — and would
also make B harmless.

### D. The stem change orphans 500 replicates for *this* cell, not just "old batches"

Check 5.3 passes, but the comment at 141–143 ("Batch `.rds` files written under
the old stem are invisible to `combine_glob` — re-run them or rename them to
match") understates what is being set aside. For **this exact cell**
(`m1_h10_knoise0_n500`) the archived pre-rename set is:

```
fs_t1_t2_m1_h10_knoise0_n500_res_1_20.rds        20 reps
fs_t1_t2_m1_h10_knoise0_n500_res_21_100.rds      80
fs_t1_t2_m1_h10_knoise0_n500_res_101_200.rds    100
fs_t1_t2_m1_h10_knoise0_n500_res_201_500.rds    300
                                            ─────────
                                                500 reps  (B = 300, draws = 5000)
plus fs_t1_t2_m1_h10_knoise0_n500_combined_1_500.rds (the pooled bundle)
```

After the stem change, `combine_glob` finds exactly one file for this cell:
`fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_100.rds`, **100 replicates**. A
combine run today silently reports on 100 reps where 500 were available before.
Nothing is corrupted and nothing is lost from disk — but "re-run them or rename
them" is a 500-replicate, B = 300 re-run, and that is worth stating explicitly
before someone renders combine mode and reads the smaller pool as the whole
study.

A straight rename would not suffice either: the archived bundles carry
`meta$gate_draws`, not `meta$mr_draws` (verified on
`fs_t1_t2_m1_h10_knoise0_n500_res_1_20.rds` and on the pooled
`..._combined_1_500.rds`). Combine-mode's `.meta_get(b, "mr_draws")` returns
`NA` for them, so pooling old with new **errors out** at the agreement check
rather than silently mixing — the failure is loud, which is the right
behaviour, but it means renaming alone will not rescue the archive.

### E. `fit_mr_secs` is fit **+** MR, but all three labels read "MR"

Structural check 4 asks only that the labels agree, and they do. They agree on
the wrong name. `fit_mr_secs` times the entire `forestsearch()` call — candidate
construction, enumeration, screening, split-sample consistency, *and* MR — yet
it is printed as `Mean wall-clock: MR 14.665 s`, as `Per rep: MR 14.7 s`, and
as the column `MR/rep (s)`.

The document's stated purpose is an MR-vs-FB head-to-head, and its
`Notes for scale-up` section says "FB is the wall-clock bottleneck (the
per-replicate timing line above quantifies the gap)". That line currently reads
14.665 s vs 16.50 s — a near-tie — because the MR side carries the search cost
that FB's 16.50 s is measured *on top of*. MR's own cost is
`fs.est$mr_inference$timing_seconds`, which the fit already exposes and the
recorder does not capture. Whatever the true ratio is, the reported one is not
it.

This is **pre-existing**, not introduced here: the pre-edit file printed the
same quantity as `Tier-2 gate 15.084 s`. The rename makes it more misleading
only in that "MR" is now a precise claim where "Tier-2 gate" was a vaguer one.
Two options, either of which is a one-line change: record
`g$timing_seconds` into a separate column and report that as MR, or rename the
labels to "fit+MR".

### F. Smaller notes (all pre-existing, none introduced by these edits)

1. **`n_workers` has a vestigial `min()`.** Line 262:
   `ceiling(0.90 * max(1L, min(parallel::detectCores(logical = FALSE) - 1L)))`
   — `min()` of a single value is a no-op, so on this 128-core box it yields
   **115 workers**, spawned per `forestsearch()` call *and* per bootstrap to run
   `nb_boots = 20` tasks. Worker spin-up dominates at smoke-test counts (both
   renders spent ~125 s on 5 replicates). Harmless for correctness; the `min()`
   looks like it lost a cap argument.
2. **Orphaned comment fragment.** Line 179, `# the MR -- see run_fb below.`,
   sits between `proj_scenarios` and `mr_draws`, detached from the `nb_boots`
   sentence it continues (102–103). The pre-edit file has the same fragment
   reading `the Tier-2 gate -- see run_tier1 below`.
3. **The comparator is not the file that produced the archived runs.**
   `fs_t1_t2_m1_h10_knoise0_n500_batch_1_20.qmd` has mtime 2026-07-30 13:17 and
   writes `meta$mr_draws`; the `.rds` files it nominally produced (June) carry
   `meta$gate_draws`. An earlier `gate_draws → mr_draws` pass had already been
   applied to the "pre-edit" file. **This does not weaken §1** — both files were
   rendered today against the same installed engine, so the comparison is
   clean — but it means §1 certifies the most recent edit round, not the whole
   vocabulary migration back to the archived runs.
4. **Two `forestsearch()` warnings surface at fit time**, on both routes
   equally: `coxph.fit: Loglik converged before variable 1; coefficient may be
   infinite` (separation on small subgroups, which `.cox_hr_ci()` explicitly
   guards) and, from the bootstrap, `nb_boots < 100 may produce unreliable
   confidence intervals` — expected at the smoke-test `B = 20`. Neither is
   related to the edits. The run-loop chunk carries `warning=FALSE`, so neither
   appears in the rendered HTML.

---

## 7. Artifacts

Everything below is under `dev/sim-check/`. Nothing outside it was written.

| path | what |
|---|---|
| `compare_rds.R` | §1 comparator (mapping + tolerance-0 per-column diff) |
| `compare_output.txt` | its full output |
| `diag_probe.R` | §2 / §3 / §5 probe |
| `run_post/diag_probe.log` | its full output |
| `run_post/attr_check.R`, `attr_check*.log` | §3 per-replicate-table drill-down |
| `run_pre/pre.qmd`, `run_pre/pre.html`, `run_pre/fs_t1_t2_m1_h10_knoise0_n500_res_1_5.rds` | pre-edit render |
| `run_post/post.qmd`, `run_post/post.html`, `run_post/fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_5.rds` | post-edit render |
| `setupcheck/setup_chunk.R` | extracted setup chunk for the clean-session source test |
