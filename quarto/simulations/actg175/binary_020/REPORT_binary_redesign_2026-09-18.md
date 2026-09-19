# REPORT — Binary study redesign: a feasible planted prevalence, the oracle four-cell condition, and a Stage 0 gate

Date: 2026-09-18 (UTC). Machine: `pop-os` (64 physical cores; R 4.6.1, forestsearch 0.3.5).
Branch `feature/glm-extension`. Task: `dev/tasks/TASK_binary_study_redesign_2026-09-18.md`
(committed as received, `625b10d0`). HEAD at the start of the work: `625b10d0`.

**Scope, as directed.** No `R/` edit. No campaign launch — the task stops before Stage 2. The two
committed binary cells are marked superseded, not deleted and not re-run.

**Prerequisites at HEAD — both present.**

- `fs_dgm_feasibility()` — `R/fs_dgm_feasibility.R`, committed `f9b794f6`, exported.
- The estimator boundary — `dev/reports/REPORT_estimability_boundary_2026-09-18.md` and
  `.fs_existence_reason()` / `.fs_nonestimable()` in `R/glm_effect_estimators.R`, committed
  `1719056f`; documented in `e55c5da6`.

**One thing had to be done that the task does not name: the installed package was stale.** The
template runs against the *installed* forestsearch (its `doFuture` workers are separate R
processes), and the installed build was `2026-09-17 04:47 UTC` — older than both prerequisite
commits, so it had no `fs_dgm_feasibility()`. `devtools::install(dependencies = FALSE, quick =
TRUE)` was run (22 s) before Step 5. The version string is unchanged at 0.3.5, so the smoke's
`pkg_version` gate is unaffected; the installed estimator code is now HEAD's, which includes the
estimability boundary.

**Stage naming.** The task says "Stage 0" for what this directory's own records call the Stage 1
pre-flight (smoke + calibration artefacts), and "Stage 2" for the production campaign, which is
this directory's Stage 2. The gate added below is the pre-flight gate and refuses the production
campaign; both names are used where each belongs.

---

## 1. Step 1 — the design chosen from the table

### 1.1 How prevalence is moved

The planted region is H = {`wtkg` > q} ∩ {`cd40` > q} at a **common quantile level** `sg_quantile`
on both cut variables, so the prevalence knob is `sg_quantile` and nothing else in the recipe
changes. The cuts are taken on the 1083-row ACTG175 analysis sample, which carries ties, so
prevalence is a **decreasing step function of q**, not a continuous one. `sg_quantile` is therefore
fixed at the **midpoint of a plateau**, never at a bisection root — a root sits on a discontinuity
and a sixth-decimal rounding moves it to the next step. The plateau table is
`prevalence_plateaus_binary_redesign_2026-09-18.csv` (35 plateaus over q ∈ [0.60, 0.70], with the
achieved super-population prevalence at each midpoint).

Each nominal prevalence in Larry's list is mapped to the plateau whose achieved prevalence is
nearest it. **What is recorded everywhere below is the achieved prevalence, never the nominal
target.**

| nominal | `sg_quantile` (plateau midpoint) | plateau | achieved prevalence(H) |
|---|---|---|---|
| 12.0% | 0.66775 | [0.6673, 0.6682] | 0.11923 |
| 12.4% | 0.65850 | [0.6581, 0.6589] | 0.12404 |
| 13.0% | 0.65160 | [0.6498, 0.6534] | 0.13215 |
| 14.0% | 0.63955 | [0.6368, 0.6423] | 0.13686 |
| 15.0% | 0.62850 | [0.6267, 0.6303] | 0.14917 |

(The 14% row is the honest one: the step structure has no plateau nearer 14% than 13.686% on one
side and 14.635% on the other.)

### 1.2 How the table was produced

`calibrate_glm_interaction()` with the study's literals — `n_super` 100000, `seed` 8316951,
`k_inter_range` (0.3, 1.5), `grid_step` 0.025, `adverse_outcome = FALSE` on `y_neg`,
`k_treat` 1 — at each (prevalence, OR) pair, then `fs_dgm_feasibility(dgm, n = c(500, 750, 1000,
2000), n.min = 60, d0.min = 10, d1.min = 10, n_rep = 200, tolerance = 0.05, effect_measure = "OR",
seed = 20260918)`.

**The DGM is built and calibrated before any RNG-kind switch.** Nothing in this task switches the
generator: every run above is under R's default Mersenne-Twister, which is the state the template's
seed table is drawn in and is drawn *before* `record_replicate()` switches to L'Ecuyer-CMRG
per replicate. `RNGkind()` was confirmed unchanged at the end of the sweep.

**The n grid.** The task names 500 / 750 / 1000 / 2000 and asks it be confirmed. The template's `n`
is a single knob (`FS_OR_N`); the *campaign* runs two sizes, n = 500 and n = 2000, and the
committed study swept 500…2000 by 250 (seven). The task's four-point grid was used: it contains the
campaign's two, and n = 500 is the binding point at every prevalence, so the selection is the same
under either grid.

### 1.3 The full table

`feasibility_binary_redesign_2026-09-18.csv` carries every column of `fs_dgm_feasibility()$table`
for all 60 rows (5 prevalences × 3 OR points × 4 sample sizes). The shares:

| prevalence(H) | OR | n=500 undecl / under-events | n=750 | n=1000 | n=2000 | feasible |
|---|---|---|---|---|---|---|
| 0.11923 | 0.75 | **0.565** / 0.135 | 0.000 / 0.010 | 0.000 / 0.000 | 0.000 / 0.000 | FALSE |
| 0.11923 | 1.00 | **0.565** / 0.095 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 | FALSE |
| 0.11923 | 1.50 | **0.565** / 0.040 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 | FALSE |
| 0.12404 | 0.75 | **0.430** / 0.095 | 0.000 / 0.005 | 0.000 / 0.000 | 0.000 / 0.000 | FALSE |
| 0.12404 | 1.00 | **0.430** / 0.065 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 | FALSE |
| 0.12404 | 1.50 | **0.430** / 0.025 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 | FALSE |
| 0.13215 | 0.75 | **0.220** / 0.055 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 | FALSE |
| 0.13215 | 1.00 | **0.220** / 0.035 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 | FALSE |
| 0.13215 | 1.50 | **0.220** / 0.015 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 | FALSE |
| 0.13686 | 0.75 | **0.100** / 0.045 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 | FALSE |
| 0.13686 | 1.00 | **0.100** / 0.025 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 | FALSE |
| 0.13686 | 1.50 | **0.100** / 0.010 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 | FALSE |
| **0.14917** | 0.75 | **0.015** / 0.010 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 | **TRUE** |
| **0.14917** | 1.00 | **0.015** / 0.005 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 | **TRUE** |
| **0.14917** | 1.50 | **0.015** / 0.005 | 0.000 / 0.000 | 0.000 / 0.000 | 0.000 / 0.000 | **TRUE** |

`share_nonestimable` is **0.000 in all 60 rows**: on this design the OR exists on the true region in
every replicate drawn, at every prevalence and every sample size. Region sizes at the chosen design:
mean 75.3 (min 56) at n = 500, 112.4 (91) at 750, 150.3 (120) at 1000, 299.3 (261) at 2000, against
`n.min = 60`.

Two structural facts the table makes plain:

1. **`feasible` is invariant to the OR design point**, because the undeclarable share depends only
   on the size of the planted region, which the covariate distribution fixes. The OR point moves
   the events floor and nothing else.
2. **n = 500 is the whole question.** At n ≥ 750 the undeclarable share is exactly 0 at every
   prevalence in the range.

### 1.4 The selection, applied mechanically

Feasible at **every** n at tolerance 0.05: **only 14.917%** (nominal 15%).

- 12.4% — the survival three-identifier grid's prevalence, the preferred point — is **not**
  feasible: 0.430 undeclarable at n = 500, more than eight times the tolerance. It is not chosen.
- The rule's fallback is the smallest feasible prevalence in Larry's range.

**Design of record: `sg_quantile = 0.62850`, prevalence(H) = 14.917%.** It is inside Larry's
12–15% range, so the STOP branch of the rule is not reached and the table is not being handed back
for a decision.

### 1.5 The DGM of record, recalibrated at the chosen prevalence

`calibrate_glm_interaction()`, same seed discipline (`seed = 8316951`, `n_super = 100000`,
`k_inter_range = (0.3, 1.5)`), at the three design points. All three solve inside the range.

| target OR(H) | k_inter | prevalence(H) | OR causal (ITT) | θ† H | θ† Hᶜ | θ‡ H | θ‡ Hᶜ |
|---|---|---|---|---|---|---|---|
| 0.75 | 0.1513020 | 0.149170 | 0.6722232577 | 0.7499999940 | 0.6564077470 | 0.7345268290 | 0.6313905111 |
| 1.00 | 0.4598307 | 0.149170 | 0.7043635080 | 0.9999999999 | 0.6564077470 | 0.9999999999 | 0.6313905111 |
| 1.50 | 0.8925310 | 0.149170 | 0.7501835133 | 1.5000000002 | 0.6564077470 | 1.5414142152 | 0.6313905111 |

The complement's references are unchanged across design points, as before: only the planted region
carries the calibrated interaction.

The DGM of record is committed as the template's knob, not as a serialized object — the template
rebuilds it deterministically from `sg_quantile` and `seed_base`, and every bundle's `meta` carries
`sg_quantile` and the five truths as poolability keys. The feasibility table sits beside it in the
study directory.

---

## 2. Step 2 — the two committed cells are superseded, not deleted

`orfs_or075_n500` and `orfs_or075_n2000` were produced at the **old** planted prevalence
**9.632%** (`sg_quantile` 0.70); the design of record is **14.917%** (`sg_quantile` 0.62850). They
are marked **superseded by design change** in `status_curated.md` §2, with both prevalences named.
No bundle is deleted, nothing is re-run, and the note says explicitly that they pool with nothing
under the new design and may not be placed beside a new cell without saying so.

The same note records that the committed study's own `mr_sweep/` payloads are at 9.632% and are
therefore not a comparator for anything run under the new design.

---

## 3. Step 3 — the oracle helper

> **Updated 2026-09-18, pre-launch (`TASK_binary_launch_stage1`, Part 1 Step 3): one in-scope
> copy, not two.** Larry's disposition is that the historical sweep driver
> `maxeffCons_mr_coverage_sweep_or075.qmd` stays **byte-identical to what produced its committed
> cells**, so the four-cell condition described in this section was reverted out of it and the
> driver is back to its pre-`625b10d0` content (blob `e9eceabe`; `git diff --stat` against
> `e55c5da6` empty). **The template of record keeps the new helper**, and it is now the single
> in-scope copy: the Stage 0 assertion compares nothing and instead requires that copy to be
> present exactly once, printing `found in 1 of 1 in-scope copies`. The reason for the split is
> that the driver's committed payloads are a historical record — an edited driver would no longer
> be the code that produced them — while the template is what every new cell runs. The paragraphs
> below record the redesign render as it stood; read “both in-scope copies” and “2 of 2” as the
> state on 2026-09-18 before this revert. `status_curated.md` §1's exception no longer applies:
> the driver is once again unedited by any campaign in this directory.


**Every copy, found by search.** Within the study directory
(`quarto/simulations/actg175/binary_020/`, the task's `Where`) there are exactly **two** copies:
`sim_fs_mr_field_or_template.qmd:471` and `maxeffCons_mr_coverage_sweep_or075.qmd:363`. A
repository-wide search finds 15 further copies in `quarto/simulations/actg175/binary/` (the legacy
sweep drivers) and more in `.claude/worktrees/`; those are outside the task's scope and were **not**
touched. Both in-scope copies were rewritten to one identical text.

**What changed.** The legacy pooled guard is kept exactly as it was — at least 5 events and 5
non-events pooled over arms — and the four-cell condition is added beside it. Both conditions
apply; either failing returns the NA quadruple, which the study's existing non-convergence
convention already counts per estimator.

```r
  # LEGACY guard, unchanged: at least 5 events and 5 non-events POOLED over arms.
  if (sum(y == 1L) < 5L || sum(y == 0L) < 5L) return(na4)
  # FOUR-CELL condition: the odds ratio does not exist unless all four
  # arm x outcome cells are non-empty -- the condition the package's estimator
  # boundary applies for OR (.fs_existence_reason(), R/glm_effect_estimators.R).
  if (sum(tt == 0L & y == 1L, na.rm = TRUE) == 0L ||
      sum(tt == 0L & y == 0L, na.rm = TRUE) == 0L ||
      sum(tt == 1L & y == 1L, na.rm = TRUE) == 0L ||
      sum(tt == 1L & y == 0L, na.rm = TRUE) == 0L) return(na4)
```

The rest of the sweep copy's diff is the template copy's existing `na4` refactor (a named NA
quadruple in place of five repeated literals, plus an explicit `!nrow(df)` test that a 0-row frame
already failed on the `length(unique(...)) < 2L` line). It is behaviour-preserving; adopting it is
what makes the two copies one text.

**The assertion.** The template's Stage 0 chunk reads `.logit_or_ci()` out of both files on disk,
compares the definitions verbatim, and **stops** if they differ. Observed in the Stage 0 render:

```
HELPER ASSERTION: .logit_or_ci() found in 2 of 2 copies; identical: TRUE
```

**Flagged:** this edits `maxeffCons_mr_coverage_sweep_or075.qmd`, which `status_curated.md` §1
recorded as *not edited by any campaign in this directory*. The task directs one helper identical
in every copy, so the edit was made and §1 now records the exception explicitly: the driver's DGM,
seeds, thresholds and rule are untouched, its committed payloads are not rewritten, and it is not
re-run. A re-run would differ from the committed payloads only on replicates whose true region has
an empty arm × outcome cell.

---

## 4. Step 4 — the Stage 0 gate

A new chunk `stage0-gate` in the template, placed immediately after the DGM build and the
evaluation frame, before any per-replicate machinery:

- calls `forestsearch::fs_dgm_feasibility()` on the DGM of record at the study's n grid, with
  **`tolerance` passed explicitly** (`FS_OR_FEAS_TOL`, default 0.05; `FS_OR_FEAS_N`, default
  `500,750,1000,2000`; `FS_OR_FEAS_NREP`, default 200; seed 20260918), and the search's own floors
  as configured in the render (`n_min`, `d0_min`, `d1_min`);
- prints the whole table and the verdict into the Stage 0 record;
- **stops the render** — Stage 2 refused — unless `feasible` is TRUE;
- the only way past a FALSE is `FS_OR_FEAS_OVERRIDE=TRUE`, which prints an explicit OVERRIDE line
  naming the shares it is overriding and is carried in the bundle as `meta$feas_override`. Same
  pattern as the DINA caps: possible deliberately, never silently.

`meta` gains `feas_feasible`, `feas_tolerance`, `feas_n_grid`, `feas_n_rep`, `feas_seed`,
`feas_override`, `feas_share_undeclarable` and `helper_identical`. `gate2.R` and `smoke_identity.R`
check named meta keys, not set equality, so the additions are compatible.

Observed in the Stage 0 render (target 0.75, n 500):

```
STAGE 0 FEASIBILITY [prevalence(H) 0.14917, sg_quantile 0.62850, target_or_h 0.75]:
  tolerance 0.050, n_rep 200, seed 20260918
     n n_rep size_mean size_q05 size_q95 size_min share_undeclarable share_under_events share_nonestimable
1  500   200    75.275    64.00    88.05       56              0.015               0.01                  0
2  750   200   112.445    97.00   128.05       91              0.000               0.00                  0
3 1000   200   150.345   131.95   169.00      120              0.000               0.00                  0
4 2000   200   299.340   272.95   325.00      261              0.000               0.00                  0
STAGE 0 FEASIBILITY: feasible = TRUE (every undeclarable share <= 0.050)
```

The in-render table reproduces the offline Step 1 table exactly, as it must — same DGM, same seed,
same `n_rep`.

---

## 5. Step 5 — the Stage 0 smoke on the new design

The study's existing machinery, unchanged in shape: one 20-replicate render of the template
(`FS_OR_TARGET=0.75 FS_OR_N=500 FS_OR_CAMPAIGN=redes FS_OR_NSIMS=20 FS_OR_WORKERS=20`, the campaign
rule `effMaxSG` / eps 0.20 / `neighborhood`, `ci_method` field), then
`scripts_or/smoke_identity.R 0.75 500 redes 20 fs`.

**Render: rc = 0, 106 s.** **Checker: PASS** — every gate, including all of §1.5(d)'s construction
checks, the meta set, zero factor-comparison and zero NA-membership warnings, 0 MR failures on 14
declared replicates, and every field identity at ≤ 2.3e-16.

**One check had to be changed, and this is the one thing in this report beyond the task's four
steps.** `smoke_identity.R` gates the DGM's truths against the committed study bundle within 1e-8
whenever the design token is `or075`. Under the redesign that comparison is between two *different*
super-populations, so it failed — correctly, and by construction, not because anything is wrong:

```
[FAIL] truth targets match the committed study within 1e-08 relative (max 0.00851)
```

Rather than leave a gate that can only fail, the comparison is now conditioned on the two bundles
planting at the same `sg_quantile` (both carry it in `meta`). When they do, the 1e-8 gate applies
unchanged. When they do not, the difference is **reported** with both values named and is not
gated — which is Step 2's "superseded by design change" applied at the point where a checker
consults the superseded cells:

```
SUPERSEDED BY DESIGN CHANGE: the committed study plants at sg_quantile 0.7; this DGM of
record plants at 0.6285, prevalence(H) 0.149170. The committed truths are not a comparator
for it, and are not gated.
SMOKE fs target_or_h=0.75 n=500: PASS
```

**Gate, as the task states it:**

| gate | result |
|---|---|
| smoke green | **PASS** (`SMOKE fs target_or_h=0.75 n=500: PASS`) |
| the assertion check passes | **PASS** (`.logit_or_ci()` found in 2 of 2 copies; identical: TRUE) — superseded, see §3: now 1 of 1 in-scope copies |
| feasibility table shows `feasible = TRUE` at every n | **PASS** (0.015 / 0 / 0 / 0 vs tolerance 0.05) |

Old vs new smoke at the same cell and the same 20 seeds, for the record (per-replicate cost is what
the launch estimate depends on):

| | declared | `n_true` mean | `n_sel` mean | family mean | `fit_mr_secs` mean | `id_secs` mean | sens/ppv |
|---|---|---|---|---|---|---|---|
| `orsmoke` (9.632%) | 15/20 | 50.6 | 106.1 | 2200 | 33.0 | 5.86 | 0.220 / 0.124 |
| `redes` (14.917%) | 14/20 | 77.5 | 107.0 | 2197 | 31.0 | 5.89 | 0.214 / 0.172 |

The planted region is half again as large and PPV rises accordingly; **per-replicate cost is
unchanged**, so the prior campaign's wall clock transfers.

Artefacts are Stage 0 and deliberately untracked, as `orsmoke`'s were:
`smoke_redes.html`, `logs_or/smoke_redes.log`,
`mr_or_harm/fs_effMaxSG_mr_field_or075_n500_nb20_redes_d5000/`.

---

## 6. For the launch go/no-go — projected campaign size

**Not a launch request.** These are the numbers the go/no-go needs.

- **Cells: 18** — 3 identifiers (`orfs` consistency, `orgrf`, `ordina`) × 6 (3 design points OR
  0.75 / 1.00 / 1.50 × 2 sizes n = 500 / 2000). Nine cells at n = 500, nine at n = 2000.
- **Replicates per cell: 2,000**, run as two seed-disjoint batches of 1,000 plus a combine, with
  `gate2.R` per cell. Total 36,000 replicates.
- **Wall clock, from this study's own prior runs** (`LOG_or_progress.txt`, 63 workers on `pop-os`,
  the only cells ever measured on this design — both `orfs`):

  | cell | batch 1 | batch 2 | combine | cell wall |
  |---|---|---|---|---|
  | `orfs_or075_n500` | 1093 s | 1129 s | 26 s | **2248 s** (37 min) |
  | `orfs_or075_n2000` | 5021 s | 4939 s | 25 s | **9986 s** (2 h 46 min) |

  Projection at 63 workers: 9 × 2248 + 9 × 9986 = **110,106 s ≈ 30.6 hours** of render wall clock,
  plus per-cell gating. The campaign's previous ceiling was `ORSG_CEILING=105165` s, i.e. **the
  projection exceeds the ceiling that was set for it** — the ceiling needs raising to about
  125,000 s or the cell list narrowing. Peak memory was 78 GB at n = 500 and 123 GB at n = 2000 on
  63 workers.

  **Three caveats on that number, all honest.** (i) Both measured cells are `orfs`; no `orgrf` or
  `ordina` cell has ever completed on this design, and the projection assumes identifier parity —
  defensible because the MR gate dominates (31 s vs 5.9 s identification in the smoke), but
  unmeasured. (ii) `orfs_or150_n500` aborted twice, once on a timeout at 9915 s and once
  unexplained; that instability is not in the projection. (iii) The measured cells ran at the old
  prevalence; §5 shows per-replicate cost essentially unchanged, which is the evidence for carrying
  them over.

---

## 7. Compute

Hard cap 20 minutes, as directed. **Used ≈ 8.0 minutes**: feasibility and calibration runs ≈ 4.7 min
(the plateau map, the 15-cell sweep at `n_rep` 200, and the three recalibrations), the package
install 22 s, the Stage 0 smoke render 106 s, the checker runs ≈ 15 s. No fits beyond the Stage 0
smoke. Nothing pushed; nothing launched.

---

## 8. Findings — no task attached

1. **The feasibility boundary inside Larry's range was not resolved.** At n = 500 and tolerance
   0.05 the undeclarable share is 0.100 at prevalence 13.686% and 0.015 at 14.917%. The crossing is
   somewhere between, and the plateau table has exactly two untested steps in that interval:
   14.635% (`sg_quantile` 0.63630) and 14.732% (0.63310). The selection rule ran over the five nominal prevalences only, as
   directed, so this was not searched. If 15% is felt to be too far from the survival grid's 12.4%,
   the cheapest next question is whether 14.635% clears the tolerance.
2. **Nothing in the 12–15% range makes n = 500 comfortable.** Even at the chosen design the mean
   region is 75 subjects against a floor of 60, with a 5th percentile of 64 — the search has about
   four subjects of headroom in the lower tail. n = 750 is the first size with real margin. Whether
   n = 500 belongs in the grid at all is a design question the feasibility table can inform but
   does not answer.
3. **`recipe` mode of `smoke_identity.R` is now inapplicable** and was not run. It compares
   rule-independent data-level columns (`n_true`, the oracle quadruples) against the committed
   study bundles, which are at the old prevalence, so those columns must differ. It is not part of
   the Stage 0 gate. Retiring it, or re-pointing it at a new reference bundle produced under the
   design of record, is an open decision.
4. **`share_nonestimable` is 0 in all 60 rows**, so the four-cell condition added to the oracle in
   Step 3 will, on this design, essentially never fire on the *true* region. Its value is on the
   *selected* region and on the small-n tail, and as a statement that the study's helper and the
   package's estimator boundary now apply the same existence condition — not as a change that will
   move any number in the campaign.
5. **The installed package was two commits behind the tree** and nobody would have noticed until a
   render failed on a missing function. The campaign runner's preflight logs `packageDescription()$Built`
   but does not compare it with `git HEAD`; it could.
