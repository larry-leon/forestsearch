# `legacy/` — archived `t1_t2`-generation combine documents

These 6 documents (each with its `.html` render, 12 files) belong to the
`t1_t2` generation, superseded by the `sim_fs_maxcons_fb_mr_*` template one
level up. They are archived for reference, not deleted, and nothing here is
regenerated or maintained.

All six are `run_mode <- "combine"`. No batch document is archived here.

## Why these six

They lack the κ / β(H) treatment: the `beta-true-region` chunk, suppression of
the oracle's `b_betaHhat` / `Cov_betaHhat` cells, and the identification-factor
footnote `κ = β(H)/β(Ĥ)`. The split across the sixteen documents named
`*_combine_*` is exactly binary — 1231 lines with no κ versus 1291 lines with
κ. Rather than backfill a generation being retired, the κ treatment now lives
in the `sim_fs_maxcons_*` template.

The four κ-carrying combine documents that stay at the parent level are
`fs_t1_t2_m1_h10_knoise0_n500_combine_1_500`,
`fs_t1_t2_m1_h15_knoise0_n1000_combine_1_500`,
`dina_t1_t2_m1_h10_knoise0_n500_combine_1_500` and
`grf_t1_t2_m1_h10_knoise0_n500_combine_1_500`.

## Not runnable in place

Every document here uses paths relative to its own directory:

* `source("betaHhat_truth.R")` — that file is one level up.
* `combine_glob <- sprintf("%s_res_*.rds", rds_stem)` — the `.rds` they read
  are one level up.

Re-rendering from `legacy/` would need those paths adjusted. Nothing in this
directory is on a render path today.

## Their outputs stay at the parent level

No `.rds` was moved. In particular the nine `*_combined_*.rds` are the outputs
of these documents and remain beside their siblings, where the cross-directory
consumers in `dev/sg-focus-work/` expect them.

## Four cells now have no combine document at the parent level

Intended. They are `h10_knoise0_n1000`, `h15_knoise0_n500`, `h15_knoise3_n1500`
and `h20_knoise3_n1500`. Their pooled `*_combined_*.rds` are unaffected and
still readable.

## Two of these cannot produce the ranges they are named for

`fs_t1_t2_m1_h20_knoise3_n1500_combine_1_250` and
`fs_t1_t2_m1_h20_knoise3_n1500_combine_1_500` both have `combine_files <- NULL`,
so the glob takes every `*_res_*.rds` for the stem — 500 replicates, not 250.
The same applies to `fs_t1_t2_m1_h10_knoise0_n1000_combine_1_400`, which pools
500 rather than 400.

Consequently `fs_t1_t2_m1_h20_knoise3_n1500_combined_1_250.rds` at the parent
level has no source that reproduces it. Preserve it; do not attempt to
regenerate it.

Several documents here also carry stale `sim_id_start` / `n_sims` values left
over from a previous batch run. They are inert in combine mode, which derives
its range from the pooled files, but they make the filename unverifiable from
the setup chunk alone. **Treat the filenames here as labels of unknown
accuracy** — the `run_mode` / `sim_id_start` / `n_sims` triple inside each file
is the only reliable description of what it does.
