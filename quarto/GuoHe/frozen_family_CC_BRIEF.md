---
title: "Frozen-family FB vs MR check — Claude Code brief"
bibliography: []
---

# Claude Code brief: frozen-family FB vs MR correspondence

## Purpose

Test the estimand claim that FB and MR diverge in the GBSG maxeff analysis
BECAUSE FB regenerates the candidate family on each bootstrap resample
(recomputing quantile cutpoints), while MR holds the fitted family fixed. The
check builds a family in which every cut is a FIXED threshold, so FB has nothing
to regenerate, and asks whether FB and MR then correspond. If they do, it
corroborates that the headline FB-MR gap is the family-regeneration component.

This is corroboration of an estimand distinction, not a new method. Application-
style single run; no simulation.

## File

```
quarto/GuoHe/analysis_gbsg_frozen_family.qmd   # NEW, self-contained
```

Requires the forestsearch package installed. Nothing else is added or modified;
`guohe_algorithm3.R`, `fs_to_guohe.R`, and the package are untouched (GH is not
even used here -- this is FB vs MR only).

## How the family is frozen

- Continuous factors (age, size, nodes, pgr, er) are DROPPED from
  `confounders.name` so forestsearch never auto-cuts them.
- Each continuous factor is instead represented by a FIXED median cut computed
  once on the observed data (`age <= 53`, `size <= 25`, `nodes <= 3`,
  `pgr <= 32.5`, `er <= 36`), plus the clinical receptor-negative cuts
  `er <= 0` and `pgr <= 0`, all passed via `conf_force`.
- `conf.cont_jcuts` is OMITTED -- no quantile cuts anywhere.
- Categorical factors meno and grade3 remain in `confounders.name`.

Fixed-threshold `conf_force` expressions are replayed verbatim under the
bootstrap (they are exempt from cut collapse; verified in get_fsdata.R), so the
family is identical on every resample.

## Gates

### Gate 0 -- the family is actually frozen (HARD)

Built into the qmd (`gate0` chunk). After the fit, it inspects the realized
candidate family and errors if any continuous-factor cut is present that is NOT
one of the fixed thresholds -- i.e. if anything quantile-derived slipped in. If
Gate 0 fails, STOP: the comparison is not interpretable, and the freezing
construction needs revisiting (do not work around it). Report the offending cut.

### Gate 1 -- render

```bash
quarto render quarto/GuoHe/analysis_gbsg_frozen_family.qmd
```

## Prediction (stated before running, so it is not rationalized after)

- Frozen family: FB debiased HR ~ MR debiased HR, up to Monte Carlo error and
  MR's leading-order approximation. The `gap` chunk prints FB-MR on the log-HR
  scale and the HR ratio; near-zero / near-one is the corroboration.
- Contrast: the headline maxeff run gave FB 1.95 vs MR 1.45 (a real gap) in a
  regenerating family. Frozen should show the gap substantially reduced.
- A small residual gap may remain from MR's leading-order approximation -- a
  distinct term. A LARGE residual would implicate that approximation instead and
  would be its own finding. Report the number either way; do not force a
  narrative.

## Known-fragile fields

The table reads the same FB/MR schema paths validated in the headline run:
`fs_bc$H_estimates$H2 / H2_lower / H2_upper` (FB) and
`g$debiased$est / lower / upper` (MR). These are confirmed correct from the
headline analysis, so they should not need adjustment; if the schema somehow
differs, fixing the read is in scope (display only, estimates unchanged).

## Reporting

Report: Gate 0 result, the selected subgroup (it will differ from the headline
maxeff subgroup -- expected, coarser family), the three-row Naive/FB/MR table,
and the FB-MR gap (log-HR difference and HR ratio). Do not compare point
estimates to the headline run -- only the within-family FB-MR gap is comparable
across the two. Do not commit; send the rendered result back for review.
```
