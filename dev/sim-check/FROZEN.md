# FROZEN — concluded investigation, do not modernise

This directory is a **completed before/after simulation-integrity
investigation**. Its findings are written up in `SIM_INTEGRITY_FINDINGS.md`
(with the rendered `SIM_INTEGRITY_FINDINGS.html`), and the brief that
commissioned it is `CC_BRIEF_sim_integrity.md`.

The rendered logs (`run_pre/render_pre.log`, `run_post/render_post.log`,
`run_post/attr_check*.log`), the `.rds` results in `run_pre/` and `run_post/`,
and `compare_rds.R` / `compare_output.txt` are **historical record**. They are
the evidence the findings rest on, and they are only meaningful as the exact
artefacts that were produced at the time.

## The four `betaHhat_truth.R` copies here are deliberate

```
run_pre/betaHhat_truth.R      run_post/betaHhat_truth.R
fixcheck/betaHhat_truth.R     setupcheck/betaHhat_truth.R
```

All four are byte-identical to the **PRE-consolidation** survival module — the
version before `R/betaHhat_truth.R` existed and before the shims landed.

They are kept for two reasons:

1. **Each document sources its own copy by bare relative path.** `pre.qmd:82`,
   `post.qmd:85`, `setup_chunk.R:19` (both copies), `attr_check.R:3`,
   `diag_probe.R:10` all call `source("betaHhat_truth.R")`, which resolves
   against the document's own working directory. That is why the copies exist:
   they are a per-directory render dependency, not accidental duplication.
2. **The pre/post comparison depends on them being what they were.** The whole
   point of `run_pre` versus `run_post` is that both sides scored targets with
   the same implementation. Replacing either with a shim would silently change
   what "before" means and invalidate the comparison the findings rest on.

## They carry D1 and D2 by design

Same status as
`dev/betaHhat-consolidation/archive/betaHhat_truth_legacy_survival.R`:

- **D1** — the `" & "` split runs ahead of disjunction dispatch, so a GRF rule
  like `"(a & b) | (c & d)"` shreds and the target returns `NA`.
- **D2** — partial-`NA` membership is consumed rather than rejected, so a rule
  the frame cannot express can yield a finite value for the wrong region.

Both are fixed in `R/betaHhat_truth.R`. Neither is fixed here, and neither
should be.

## Do not

- Do not shim these files.
- Do not update them to call the package.
- Do not re-point any document in this directory.
- Do not "fix" the defects above.

If the harness is ever re-run against the consolidated package, **that is a new
investigation and belongs in a new directory.** Re-running it in place would
overwrite the record that makes the original findings checkable.

## Known-broken at freeze time: the three top-level sourcers

These three call `source("betaHhat_truth.R")` and it does **not** resolve,
because `dev/sim-check/betaHhat_truth.R` has never existed in any commit:

```
fs_t1_t2_m1_h10_knoise0_n500_batch_1_20.qmd:82
sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd:85
diag_probe.R:10
```

**This is not damage and must not be repaired.** These are the *inputs* the
investigation copied in to compare, named as such in its own documents:

- `CC_BRIEF_sim_integrity.md:27-28` designates
  `fs_t1_t2_..._batch_1_20.qmd` as **pre-edit (the original)** and
  `sim_fs_maxcons_..._batch_1_100.qmd` as **post-edit**.
- `SIM_INTEGRITY_FINDINGS.md:4-5` names them Target and Comparator.
- `SIM_INTEGRITY_FINDINGS.md:164-165` cites `diag_probe.R` as the script behind
  a logged result.

The renders that produced the findings were run from `run_pre/` and `run_post/`,
which have their own `pre.qmd` / `post.qmd` **and** their own
`betaHhat_truth.R`. That is why those resolve and these never needed to.

Restoring a copy here would place a **current shim beside four deliberately
pre-consolidation ones**, in a directory frozen precisely so that "before"
stays what it was. Recorded as known-broken at freeze; left alone.
