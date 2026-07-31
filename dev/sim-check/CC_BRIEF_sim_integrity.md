# CC BRIEF — integrity check on `sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd`

```
claude "Read dev/sim-check/CC_BRIEF_sim_integrity.md and execute it."
```

The simulation was revised through a sequence of edits that were all asserted to
be inert: a Tier→FB / gate→MR vocabulary rename, `sg_focus "eff"` → `"maxcons"`,
several parameters pinned to values they already resolved to, a definition
reordering, and an output-stem change. Two defects have already been found after
the fact (`t2m` missed by a word-boundary regex; `focus_tag` referencing
`sg_focus` sixty lines before it was defined). This brief establishes whether
anything else was missed.

**Change nothing.** Write only under `dev/sim-check/`. Report and stop.

---

## 1. The decisive test: the edits should be bit-identical

Every edit was either a rename, a reorder, or a pin to a value the call already
resolved to. So the pre-edit and post-edit files, run on the same codebase with
the same seeds, must produce **identical per-replicate records**.

| | file |
|---|---|
| pre-edit | `fs_t1_t2_m1_h10_knoise0_n500_batch_1_20.qmd` (the original) |
| post-edit | `sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd` |

Copy both to `dev/sim-check/`, and in **each copy** set:

```r
n_sims       <- 5L
sim_id_start <- 1L
nb_boots     <- 20L      # enough to exercise the FB path; not a timing test
mr_draws     <- 500L
```

Change nothing else. Render both. Then compare the two result `.rds` files
column by column at `tolerance = 0`, matching records on `sim_id`.

Note the column names differ by construction — pre-edit has `t1_*`/`t2_*`, post
has `fb_*`/`mr_*`, and `t2_secs` is now `fit_mr_secs`. Map them:
`t1_` → `fb_`, `t2_` → `mr_`, `t2_secs` → `fit_mr_secs`. After that mapping the
column **sets** should match exactly; report any name present in one and not the
other, since a dangling name is how a missed rename would surface.

Wall-clock columns will differ and are exempt. Everything else — `sg_def`,
`detected`, `mr_ok`, all estimates, intervals, standard errors, and the
classification metrics — must be identical. **Any difference here means one of
the edits was not inert. Report it with the column and the replicate.**

---

## 2. Do the pinned parameters actually reach the engine?

Run one replicate with `diag_mode = TRUE` and capture `fs.est`. From
`fs.est$args_call_all`, report the resolved value of each:

| parameter | expected |
|---|---|
| `sg_focus` | `"hr"` — `"maxcons"` normalises to the canonical name |
| `subgroup_method` | `"consistency"` |
| `consistency_method` | `"resample"` |
| `stop_threshold` | `NULL` |
| `selection_rule` | `"neighborhood"` |
| `effect_neighborhood` | `0.10` |
| `mr_inference` | `TRUE` |

`stop_threshold` is the one to look at closely. It is now pinned to `NULL`
explicitly, where before it was inherited and reset internally. Both routes
should end at `NULL`; confirm the engine saw `NULL` either way, and report
whether the render emits the reset warning (it should **not**, because passing
`NULL` means the reset branch is never entered).

---

## 3. Does `mr_in_replicates = FALSE` do what the comment claims?

The `forestsearch_bootstrap_dofuture()` call now pins it. The fit carries
`mr_inference = TRUE`, so without the pin the flag would be inherited into every
replicate.

Confirm by running the FB bootstrap twice on one fitted `fs.est`, same seed,
once each way:

* `mr_in_replicates = FALSE` (as pinned)
* `mr_in_replicates = TRUE`

Report whether `H_estimates`, `Hc_estimates` and the per-replicate table are
identical, and the wall-clock difference. Identical results plus a large timing
gap is the expected outcome — it confirms the pin is a pure cost saving.

---

## 4. Structural checks

1. **Symbol resolution.** Source the setup chunk alone in a clean session with
   `betaHhat_truth.R` present. It must run to completion with no
   `object ... not found`. Report the resolved values of `rds_stem`, `rds_path`
   and `combine_glob`.
2. **Stem matches filename.** `rds_stem` should be
   `fs_maxcons_fb_mr_m1_h10_knoise0_n500`. Confirm `rds_path` is what the run
   actually writes, and that `combine_glob` matches it.
3. **Old batches are invisible.** Confirm any pre-existing
   `fs_t1_t2_m1_*_res_*.rds` files are NOT picked up by the new `combine_glob`.
   This is intended — the stem changed — but it must be stated rather than
   discovered later by a short pooled run.
4. **`fit_mr_secs` measures what it says.** Confirm `mr_t0` is set immediately
   before the `forestsearch()` call, so the recorded value is fit **plus** MR,
   not MR alone. It is used in two printed summaries; confirm both labels agree.
5. **Field paths.** Verify against `str()` on a real fit that `fb_*` reads from
   `boot$H_estimates$H2` / `$H2_lower` / `$H2_upper` / `$sdH2`, and `mr_*` from
   `fs.est$mr_inference$debiased$est` / `$lower` / `$upper` / `$se_ij`. A wrong
   path yields silent `NA`, not an error.
6. **No dangling references.** Scan the chunks for any `t1`/`t2` symbol that
   survived the rename, and for any `fb_*`/`mr_*` name that is assigned but
   never read, or read but never assigned.

---

## 5. Deliverable

`dev/sim-check/SIM_INTEGRITY_FINDINGS.md`:

1. **Verdict up front.** Are the edits inert? Is anything wrong?
2. §1 comparison table: columns compared, columns differing, with the mapping applied.
3. §2 resolved-parameter table.
4. §3 result identity plus timings.
5. §4 structural checks, pass/fail each.
6. Anything found that this brief did not anticipate.

Do not fix anything. If a check fails, report it and continue with the rest.
