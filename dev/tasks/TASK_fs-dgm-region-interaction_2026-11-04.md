# TASK_fs-dgm-region-interaction_2026-11-04 — region-dependent treatment effect in the DGM builder

Repo: **forestsearch** (Pop!_OS box). `export GIT_TERMINAL_PROMPT=0`; CC never pushes. Standing
first action: copy this task document from `~/Downloads/` into `dev/tasks/` and commit it alone.
Short-lived branch off `main` (check `git branch -r` first), fast-forward, delete. NOTES to
`~/Downloads`.

## R/ callout (required by convention)

This task TOUCHES R/. Classification: **adds capability to the simulation DGM builder via new
optional arguments whose defaults reproduce current behaviour exactly** — it moves no existing
code and changes no inference method. Approved by Larry 2026-11-04 (the §2 extension discussion:
region-dependent ψ written as treatment-by-region interactions).

## What is added

The extended patient-level log hazard ratio, approved wording in the manuscript task:
ψ(X) = log b + log(h/b)·1{X ∈ H*} + γ_R·1{X ∈ R} + γ_{H*R}·1{X ∈ H*}·1{X ∈ R},
with γ_R = log(b_R/b) and γ_{H*R} = log(h_R/h) − log(b_R/b). Semantics: treatment hazard ratio is
b outside / b_R inside the region for patients outside H*, and h outside / h_R inside the region
for patients in H*. Base model = b_R = b, h_R = h (both γ's zero).

## Steps

1. Read the CURRENT builder source (`inject_mrct_structure()` and whatever helpers construct the
   θ/ψ columns) and the analytic-layer entry points (`ahr_population()`, `ahr_limit()`,
   `ahr_design()`) — verify from source, never from this description. Record in NOTES whether the
   analytic layer consumes per-row (θ0, θ1) (expected, per the finite-sum design) or parametric
   (b, h)-style inputs anywhere on its computation path.
2. Extend the builder with two new optional arguments for the region-dependent effect — named per
   the repo's argument style; semantics exactly as above (accept either (b_R, h_R) or (γ_R,
   γ_{H*R}) parameterization, one of them, matching house style); **defaults reproduce the current
   builder exactly** (γ's = 0).
3. If step 1 found any parametric (b, h) assumption on the analytic path (e.g., a closed-form AHR
   or CDE hard-coded as b^{1−p}h^{p} / (1−ω⁰)b + ω⁰h), extend it minimally to the four-cell form
   (shares over H* × region); if the path is per-row throughout, record "no analytic-layer change
   needed" and change nothing there.
4. Roxygen for the new arguments; NEWS entry; version bump per repo practice (dev version; release
   packaging into the pending 0.4.2 is Larry's separate call — record the resulting version string
   in NOTES for the manuscript's fs_min_version).
5. Tests (testthat), all new:
   - T1 Backward compatibility: with defaults, the builder's ψ/θ columns on the shipped
     super-population are IDENTICAL (exact equality) to the current builder's output.
   - T2 Regional null: with b_R = h_R = 1, every region row has ψ = 0, and `ahr_limit()` on the
     region set returns |β*_R| < 1e-6 (margHR = 1).
   - T3 Four-cell AHR identity: on a set straddling the boundary, exp(mean ψ) equals
     b^{p00} h^{p10} b_R^{p01} h_R^{p11} with the four H*×region shares (exact).
   - T4 Interaction algebra: per-row ψ equals the γ-form and the (b_R, h_R)-form for a grid of
     parameter values (exact).
6. Run the repo's standard package check target and the existing test suite.

## Gates (stop-on-failure)

- G1 All existing tests pass unchanged; new tests T1–T4 pass; standard package check green by the
  repo's usual route.
- G2 T1 is exact: default invocation output byte-identical to current behaviour.
- G3 Scope: `git diff --name-only` confined to the builder/analytic R files touched per steps 2–3,
  `man/`, `tests/`, `NEWS`/`DESCRIPTION` (+ `dev/tasks/`). Nothing else.
- G4 NOTES records: per-row vs parametric finding (step 1), files touched, new version string.
- G5 Clean close: fast-forward `main`, delete branch, `git status --porcelain` empty.

No compute: this task runs no simulation campaign. Downstream: TASK_mrct-sim-dgm-section-predictions
(mrct-consistency repo) depends on this task being green.
