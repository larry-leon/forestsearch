# `legacy/` — archived `t1_t2`-generation documents

These 12 documents (each with its `.html` render, 24 files) belong to the
`t1_t2` generation, superseded by the `sim_fs_maxcons_fb_mr_*` template one
level up. They are archived for reference, not deleted, and nothing here is
regenerated or maintained.

## Why these twelve

They lack the κ / β(H) treatment: the `beta-true-region` chunk, suppression of
the oracle's `b_betaHhat` / `Cov_betaHhat` cells, and the identification-factor
footnote `κ = β(H)/β(Ĥ)`. The split across the sixteen documents named
`*_combine_*` is exactly binary — 1231 lines with no κ (these twelve) versus
1291 lines with κ (the four that remain at the parent level). Rather than
backfill a generation being retired, the κ treatment now lives in the
`sim_fs_maxcons_*` template.

The four κ-carrying documents that stay at the parent level are
`fs_t1_t2_m1_h10_knoise0_n500_combine_1_500`,
`fs_t1_t2_m1_h15_knoise0_n1000_combine_1_500`,
`dina_t1_t2_m1_h10_knoise0_n500_combine_1_500` and
`grf_t1_t2_m1_h10_knoise0_n500_combine_1_500`.

## Not runnable in place

Every document here uses paths relative to its own directory:

* `source("betaHhat_truth.R")` — that file is one level up.
* `combine_glob <- sprintf("%s_res_*.rds", rds_stem)` and `rds_path` — the
  `.rds` they read and write are one level up.

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

## Six of these are batch documents, despite their names

Worth knowing before anything here is revived. Six files named `..._combine_*`
have `run_mode <- "batch"` and are batch runners, not combine runners. They are
the sole sources for six of the seven `res_*.rds` in the `h20_knoise3_n1500`
cell:

| file (all `fs_t1_t2_m1_h20_knoise3_n1500_…`) | `run_mode` | `sim_id_start` / `n_sims` | writes |
|---|---|---|---|
| `combine_1_100` | **batch** | 1 / 100 | `…_res_1_100.rds` |
| `combine_101_150` | **batch** | 101 / 50 | `…_res_101_150.rds` |
| `combine_151_200` | **batch** | 151 / 50 | `…_res_151_200.rds` |
| `combine_201_250` | **batch** | 201 / 50 | `…_res_201_250.rds` |
| `combine_251_350` | **batch** | 251 / 100 | `…_res_251_350.rds` |
| `combine_351_400` | **batch** | 351 / 50 | `…_res_351_400.rds` |
| `combine_1_250` | combine | — | pools |
| `combine_1_500` | combine | — | pools |

The seventh, `…_res_401_500.rds`, comes from
`fs_t1_t2_m1_h20_knoise3_n1500_batch_401_500.qmd`, which stays at the parent
level.

So archiving these moves six batch sources out of the top level even though the
`t1_t2` batch family is otherwise being kept in place. Nothing is lost — the
files are here and the `.rds` they produced are untouched — but if that cell is
ever re-run, its batch sources are in this directory, not beside the others.

Two further name/content mismatches in the same set: `combine_1_400` and
`combine_1_250` cannot produce the ranges they are named for. Both have
`combine_files <- NULL`, so the glob takes every `*_res_*.rds` for the stem —
500 replicates, not 400 or 250. `fs_t1_t2_m1_h20_knoise3_n1500_combined_1_250.rds`
at the parent level therefore has no source that reproduces it; preserve it, do
not attempt to regenerate it.

**Treat the filenames here as labels of unknown accuracy.** The
`run_mode` / `sim_id_start` / `n_sims` triple inside each file is the only
reliable description of what it does.
