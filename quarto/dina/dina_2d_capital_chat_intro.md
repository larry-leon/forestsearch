# New chat introduction — DINA 1-D → 2-D extension; evaluate CAPITAL

> Start this in a **new chat inside the forestsearch Project** so the
> architecture context (memory + project knowledge) carries over while the
> bootstrap-debugging history is shed.

**Goal:** Extend DINA subgroup identification in the `forestsearch` package
from 1-D to 2-D, and first evaluate whether algorithms from the CAPITAL
methodology can be reused for it. I'll attach the CAPITAL codebase in this chat.

## Background

DINA is the model-based subgroup engine in forestsearch (`dina.R`,
`dina_subgroup.R`, `dina_df.R`, plus the Pareto-frontier machinery). It fits a
cross-fitted DINA model and selects a subgroup via `dina_subgroup()` under
forestsearch's own criteria — `sg_focus`, `selection_rule`,
`effect_neighborhood`, `n.min`, and a harm floor `m_diff`
(= `log(hr.threshold)` for ratio families, identity for gaussian). It's
reachable two ways:

- **screening contributor** (`use_dina = TRUE`, which now defaults to
  contributing the single `dina_subgroup()`-selected cut via
  `dina_args$selected_only = TRUE`), and
- **selection engine** (`subgroup_method = "dina"`).

Both are currently **depth-1**: a single covariate threshold, e.g.
`{nodes >= 10}`. The candidate search lives in `subgroup_search.R` /
the `.dina_collect_candidates()` + `dina_frontier()` path; membership/labels
use the `{var op q}` convention via `get_dfpred`.

## What I want to build

Depth-2 DINA: AND-conjunctions of two covariate thresholds (e.g.
`{nodes >= q1 & age >= q2}`), selected under the *same*
`sg_focus`/`selection_rule`/`m_diff`/`n.min` criteria, producing the same
interpretable `{...}` label form so it composes into existing membership,
frontier display, and the bootstrap bias-correction without conceptual change.

## First task — evaluate CAPITAL (before any implementation)

Once I attach the CAPITAL codebase, assess whether its algorithms are reusable
for the 1-D → 2-D extension, against four questions:

1. **Search strategy** over the 2-D candidate space — exhaustive grid,
   greedy/sequential refinement, policy-tree, pruned — and whether complexity
   stays tractable at forestsearch scale (K covariates, N candidates).
2. **Objective alignment** — does CAPITAL optimize something compatible with
   DINA's cross-fitted tau-hat estimand and harm floor `m_diff`, or a
   different target (RMST, constrained policy value) needing translation?
3. **Output form** — interpretable axis-aligned conjunction (composes into
   `{var op q1 & var2 op q2}` + `get_dfpred`), or a more general rule that
   wouldn't.
4. **Inference interaction** — does CAPITAL's selection-adjustment conflict
   with or complement the existing forestsearch bootstrap bias-correction?

Please start by reviewing the CAPITAL code I attach and giving a structured
**"reusable / adapt / not a fit"** verdict against those four questions, with
the relevant CAPITAL functions identified — before proposing any design. Work
from the attached code, not priors about CAPITAL.

## Conventions reminder

- Most recent in-chat file version is the source of truth.
- Modify only what's requested; ask if scope is unclear.
- End any file-producing turn with the files in the panel as downloadable
  artifacts.
- Quarto living-doc smoke tests over testthat.
- `devtools::install()` (not `load_all()`) so parallel workers see changes.
