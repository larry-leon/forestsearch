# Extreme-subgroup simulation documents

Research documents for the behaviour of subgroup estimates under a uniform
treatment effect. They use `run_subgroup_sims()`, `subgroup_cox()` /
`subgroup_glm()` and `compare_subgroup_sims()`. The published pkgdown article
is a separate, adapted copy: `vignettes/articles/extreme_subgroups.qmd`.

| Document | Data / endpoint | Baseline design |
|---|---|---|
| `gbsg_survival_resample.qmd` | GBSG, recurrence-free survival | random-X (`baseline = "resample"`) |
| `gbsg_survival_fixed.qmd` | GBSG, recurrence-free survival | fixed-X (`baseline = "fixed"`) |
| `gbsg_design_comparison.qmd` | reads the two GBSG payloads | fixed vs random comparison (no simulation) |
| `actg175_continuous_fixed.qmd` | ACTG 175 (arms 1 vs 3), CD4 change | fixed-X |

## Rendering

The documents default to 10,000 replicates. Rendering is a compute step.

    mkdir -p results
    quarto render gbsg_survival_resample.qmd            # writes results/extreme_sims_resample_10000_payload.rds
    quarto render gbsg_survival_fixed.qmd               # writes results/extreme_sims_fixed_10000_payload.rds
    quarto render gbsg_design_comparison.qmd            # reads both payloads
    quarto render gbsg_survival_fixed.qmd -P n_sims:50  # quick check

`results/`, `*.rds` and rendered `*.html` are git-ignored. Earlier versions,
renders, slides and the August archive are in git history at tag
`extreme-subgroups-pre-reorg`. Slides and the v1 feasibility document now
live in the private `hta-subpopulations` repository.
