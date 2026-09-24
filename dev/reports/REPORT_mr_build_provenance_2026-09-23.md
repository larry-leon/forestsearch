# REPORT — Build provenance: the MR alignment fits, and the cert20 before arm

- **Task:** `dev/tasks/TASK_mr_build_provenance_2026-09-23.md` (`51699cd5`), with Larry's amendment reversing its
  "do not run `devtools::install()`" line (verbatim below).
- **§2: the before/after table of `REPORT_mr_admission_alignment_2026-09-23.md` stands. No re-run.** MR's
  admission floor is evaluated in the calling process, never on a worker, and the calling process ran the
  modified source tree.
- **§3: `cert20` is a valid before arm.** 153 / 153 columns identical, maximum difference 0. The same holds for
  `e1stud`, the source of the other two 31% cells.

## 0. Record

| item | value |
|---|---|
| HEAD at start | `51699cd5` (task doc); `R/` last changed at `7713942e` |
| R, platform | R 4.6.1, pop-os, Linux 7.1.5 |
| Installed forestsearch (before the amendment's install) | `/home/larryleon/R/x86_64-pc-linux-gnu-library/4.6/forestsearch`, 0.3.5.9000, packaged 2026-09-23 02:32:06 UTC, built 02:32:08 UTC; `pconsistency.digits` in `formals(fs_mr_inference)`: **FALSE** (stale) |
| Untracked before the task (never staged) | the two `actg175/.../_d5000/` directories, `actg175/binary_020/smoke_redes.html`, `smoke_relaunch.html`, `gbsg_020/scripts_dinamr/logs/nullmr_findings.err` |

## Amendment (Larry, 2026-09-23; recorded verbatim)

> its "do not run devtools::install()" line is reversed. After completing its Steps 1-4, run devtools::install()
> to put HEAD into the main library, and confirm afterwards that the installed fs_mr_inference has
> pconsistency.digits. Run the install LAST, once every other run in this session has finished, so nothing is
> loading from the library while it is being written. Report the installed version and build timestamp.
>
> Once HEAD is installed, note in the report that the scratch-library route via R_LIBS is no longer needed for
> future runs.

## 1. Where MR's admission runs — from source

**MR re-selection, admission floor and field draws: the calling process.** `forestsearch()` calls MR directly
in its own body, not inside any `future`:

```r
# R/forestsearch_main.R:2703 (consistency path; DINA/GRF likewise via the same wrapper)
      out$mr_inference <- .fs_apply_mr(
        df = .mr_df, candidates = .mr_fam,
# R/forestsearch_main.R:3779
      fs_mr_inference(
        df = df.fs, candidates = fam, spec = gspec,
```

`.fs_apply_mr()` (`R/fs_mr_inference_methods.R:141`) is a thin wrapper that calls `fs_mr_inference(` at `:157`.
The file contains no `future`, `foreach`, `%do%`, `mclapply` or `parLapply`. `R/fs_mr_inference.R` contains no
parallel dispatch either: its four matches for "future" are doc comments naming
`forestsearch_bootstrap_dofuture()` (`:10, :238, :576, :585`). The admission floor is built and applied in-line:

```r
# R/fs_mr_inference.R:675-681
    digits <- pconsistency.digits
    if (is.null(digits))
      digits <- eval(formals(subgroup.consistency)$pconsistency.digits)
    z   <- stats::qnorm((1 + .fs_pcons_eff(admission$consistency$p_star,
                                           as.integer(digits))) / 2)
    t_g <- pmax(admission$effect_floor, c_cons + z * sdv)
    .admit <- function(bs) which(bs >= t_g)
```

It is applied over the field draws by `.admit(bs)` at `:713` (re-selection) and `:930`, with the multipliers
drawn at `:706` (`.fs_mr_multipliers()`), all sequential R in the same frame.

**Candidate search: multisession workers.** `R/subgroup_search.R:164`
`future::plan(future::multisession, workers = parallel_workers)`, then `:275`
`future.apply::future_lapply(seq_len(tot_counts), function(kk) {`. **Consistency splits: workers too**,
`R/subgroup_consistency_main.R:987` `future.apply::future_lapply(batch_indices, eval_fun, ...)`. Namespace
functions called on those workers resolve from the workers' library, so under `load_all()` they run the
**installed** build. The alignment report already records this (`:194`).

## 2. The six MR alignment fits

**The hypothesis holds, from §1.** The stale build's `fs_mr_inference()` has no `pconsistency.digits` formal
(confirmed on the installed build, §0), so stale code in the admitting process would apply the exact cutoff for
every `digits`. The `digits = 2` after fit and the `digits = 6` fit would then both equal the baseline, and
there would be no before/after difference at all. What was observed is the aligned code's prediction on both
counts: `digits = 2` differs (26 / 5,000 winners), and `digits = 6` is `identical()` to the baseline. That second
point holds because `.fs_pcons_eff(0.90, 6)` = 0.8999995, a z within rounding of the exact cutoff.

**A mixture is ruled out by §1, not inferred from the change.** A main-aligned / worker-stale mixture would need
the admission floor to be evaluated on workers, and it never is. The workers ran only the candidate search and
the consistency splits, which the alignment does not touch. Both the before and after fits used the same
workers, and the alignment report verified the FS fit's captured family and analysis frame `identical()`
before and after (its `:24`).

**The calling process ran the modified code.** From the saved fit objects (`~/Downloads/mr_admission_align_{before_fs,after_fs,after_fs_d6}.rds`):

- `session$otherPkgs$forestsearch` has `attr(, "file")` = `/home/larryleon/Documents/GitHub/forestsearch/DESCRIPTION`,
  with no `Built` field and `pkgload`/`devtools` loaded. That is a `load_all()` of the source tree, not the
  library.
- `args_call_all$pconsistency.digits` is `2`, `2` and `6`. That argument exists on `forestsearch()` only from
  `06ac5391`, which the installed build predates (task premise; the installed `fs_mr_inference` also lacks it).
- `pkg_head` is `ba595f4b` in all three: the fix was an uncommitted working-tree change when the fits ran,
  committed afterwards as `7713942e`.
- The alignment report's replay (`:139-141`): MR's re-selection loop re-run on the traced inputs reproduces the
  before winners at `z_exact` and the after winners at `z_eff`, exactly.

**Disposition: the first of the three.** The admission floor is evaluated in the calling process, which ran the
modified code, so the before/after table stands. No re-run, and no amendment to the alignment report.

## 3. The cert20 before arm

Run before the Section 5 sweep, recorded in `REPORT_mr_alignment_section5_sweep_2026-09-23.md` §1, and
repeated here as this task's deliverable.

- **Export:** `git archive ba595f4b` to a scratch directory outside the repo.
- **Deviation from the task text:** the task says "run it with `devtools::load_all()`". §1 of this same task
  establishes that `load_all()` does not reach multisession workers, and the template's replicates run on
  workers. So a `load_all()` run would put the installed stale build on the workers, not `ba595f4b`. Instead the
  export was **installed** to its own library
  (`/tmp/claude-1000/-home-larryleon-Documents-GitHub-forestsearch/ab0a4121-64a0-4c61-972b-998577a64b16/scratchpad/Rlib_pre`)
  and put first via `R_LIBS`. Asserted inside all 64 workers: that library path, and `pconsistency.digits`
  absent from `fs_mr_inference` (`mrs5pre/logs/mrs5pre_C31_h100_n500_gateP_batch_1.buildcheck`).
- **Cell:** `cert20` HR 1.00, n 500, prevalence 31% (`z1q60`); 10 replicates, sim_id 1–10, seeds
  `8316951 + sim_id`; the committed template via `render.sh`, `cert20`'s knob set.

| committed bundle | columns compared | identical | max abs diff |
|---|---|---|---|
| `results/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n500_z1q60_nb20_cert20_combined_1_2000.rds` | 153 (every shared column except the five `*_secs` timings) | **153** | **0** |
| `results/fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_e1stud_combined_1_2000.rds` (the sweep's other 31% source) | 153 | **153** | **0** |

**Identical: `cert20`'s committed payload is a valid before arm, and the sweep's 31% cells stand.** Payloads:
`results/..._n500_z1q60_nb20_mrs5pre_res_1_10.rds` (h100, h150), committed in `7dd22c9e`.

## Post-conditions

1. §1's dispatch quoted from source, with the candidate search, the consistency splits, and MR (re-selection,
   admission, field draws) reported separately: met.
2. §2 reaches disposition 1, stated, on evidence from source and from the fit objects, not from "the results
   changed": met.
3. No re-run was needed: n/a.
4. §3 reports 153 / 153 and max diff 0: met.
5. No `R/` file modified: met (`git status --short -- R/` empty).
6. `devtools::install()`: **reversed by the amendment**; see §5, run after the commits above.
