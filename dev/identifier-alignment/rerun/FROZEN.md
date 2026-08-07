# FROZEN — concluded audit re-run, frozen against churn

This directory is the **audit-era re-run** that regenerated the ACTG175 and
GBSG multimethod analyses. Its context is
`dev/identifier-alignment/NOTE_actg175_naive_interval_rootcause.md`, which
records the re-run as already done and living here.

The rendered `analysis_*.html` files and `actg175_table_payload.rds` are the
outputs that audit produced. They are historical record.

## `betaHhat_truth.R` here is a CURRENT SHIM, not a defective copy

This is the important difference from `dev/sim-check/`, which is frozen for the
opposite reason.

The copy here is byte-identical to the live shim at
`quarto/simulations/gbsg_redux/betaHhat_truth.R`. It was re-synced when the
shims landed, so it **delegates to `R/betaHhat_truth.R` and carries no
defect** — not D1, not D2. Sourcing it is safe.

It is kept because **three sibling documents source it by bare relative
path**:

```
sim_dina_maxcons_fb_mr_m1_h10_knoise0_n500.qmd:97   source("betaHhat_truth.R")
sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500.qmd:97     source("betaHhat_truth.R")
sim_grf_maxcons_fb_mr_m1_h10_knoise0_n500.qmd:97    source("betaHhat_truth.R")
```

Those calls resolve against this directory. Deleting the file breaks all three
renders. A path-based grep for `identifier-alignment/rerun` does **not** find
them — the reference is a bare filename — which is exactly how this copy was
nearly deleted as unreferenced.

## Frozen against churn, not against defects

- Do not delete `betaHhat_truth.R` — three documents depend on it.
- Do not re-point the documents at `library(forestsearch)`; editing a completed
  audit's inputs is what freezing is meant to prevent.
- Do not add logic here. If the package shim changes, re-sync this copy from
  `quarto/simulations/gbsg_redux/betaHhat_truth.R` rather than editing in place.

**If the audit is ever re-run, the shim delegates correctly as-is** — no
migration work is needed first. That is the whole benefit of it being a shim
rather than a frozen pre-consolidation copy.

The same reasoning and the same policy apply to the sibling copy at
`dev/identifier-alignment/sim_analyses/betaHhat_truth.R`, restored 2026-08-07
after being found absent while its document sourced it.
