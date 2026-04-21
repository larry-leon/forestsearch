# Legacy vs Current ForestSearch: Code Audit

**Scope.**  Six legacy files uploaded by Larry (paper-vintage codebase) vs
the current forestsearch KB.

| Legacy file                            | Current counterpart (KB)              |
|----------------------------------------|----------------------------------------|
| `forest_search_v0.R`                   | `R/forestsearch_main.R`                |
| `subgroup_search_v0.R`                 | `R/subgroup_search.R`                  |
| `subgroup_consistency_v0.R`            | `R/subgroup_consistency_main.R` + helpers |
| `forest_search_cross-validation_v0.R`  | `R/forestsearch_cross_validation.R`    |
| `get_FSdata_v0.R`                      | `R/get_fsdata.R` + `get_FSdata_helpers.R` |
| `forestsearch_functions_v0.R`          | `R/forestsearch_helpers.R` (partially) |

**What the user has pre-flagged as an intentional difference:**
the legacy "stopping for futility" (`pstop_futile`) mechanism is gone
in the current code.  That is not treated as a bug in this audit.

**Method.**  Legacy files read in full; current files re-read via
`project_knowledge_search` on each functional area.  For each area I
note (i) *what changed*, (ii) *whether the change is algorithmic or
architectural*, (iii) *whether it plausibly affects paper
reproducibility*.

---

## 1. Findings ordered by expected reproducibility impact

### §1.1 — HIGH: Sample-splitting mechanism (50/50 vs independent Bernoulli)

This is the single most important difference I found, and it is not
limited by `use_twostage`.

**Legacy** (`subgroup_consistency_v0.R`, the loop at line 157-169):

```r
n.split   <- round(nrow(df.x) / 2, 0)
datax.id  <- c(1:nrow(df.x))
idx.split1 <- sample(datax.id, size = n.split, replace = FALSE)   # <<<
df.x.split1 <- df.x[idx.split1, ]
idx.split2  <- setdiff(datax.id, idx.split1)
df.x.split2 <- df.x[idx.split2, ]
```

Every split is **exactly 50/50** (within rounding).  `sample(...,
replace = FALSE)` draws `n/2` unique row indices; the other half is
the set complement.

**Current** (`run_single_consistency_split()` in
`subgroup_consistency_helpers.R`):

```r
in.split1 <- sample(c(TRUE, FALSE), N.x, replace = TRUE,
                    prob = c(0.5, 0.5))        # <<<
df.x.split1 <- df.x[in.split1 == TRUE,  ]
df.x.split2 <- df.x[in.split1 == FALSE, ]
```

Every observation is flipped **independently** with probability 0.5.
Split sizes are `Binomial(N_x, 0.5)` — typically differ by a few
subjects, occasionally by tens.  Standard deviation of the size
imbalance is `√(N_x/4)` ≈ 5 for a subgroup of N_x = 100.

**Why this matters for CV reproducibility.**  Each consistency
split's HR estimate has larger variance under the Bernoulli scheme,
particularly for small subgroups where the extra noise from unequal
splits is substantial relative to the HR signal.  Pcons therefore
fluctuates more across splits in the current code than in the legacy
code.  For the GBSG "er <= 0" subgroup (n ≈ 82 in the harm arm,
n ≈ 604 in the complement), the harm subgroup is the noise-vulnerable
one.  Across 200 10-fold CV simulations, this extra per-split variance
translates directly into lower fold-to-fold stability of the
identified subgroup — exactly the pattern Larry observed (sens_H ≈
51% vs the paper's 73%).

**Status.**  Algorithmic change, not documented as paper-breaking.
Worth either (a) restoring the exact 50/50 split as an option, or
(b) recalibrating the paper's reported CV metrics to the new splitter.

### §1.2 — HIGH: Pcons denominator differs when a split is non-estimable

**Legacy** (same loop as §1.1):

```r
flag.consistency <- rep(0, n.splits)             # initialize to 0
for (bb in 1:n.splits) {
  ...
  if (!inherits(hr.split1, "try-error") &
      !inherits(hr.split2, "try-error")) {
    flag.consistency[bb] <- (hr.split1[1, 1] > hr.consistency &
                             hr.split2[1, 1] > hr.consistency)
  }
}
p.consistency <- mean(flag.consistency, na.rm = TRUE)
```

Non-estimable splits **stay at 0** and are counted in the denominator
via `mean(...)`.  I.e., Pcons = (# consistent splits) / n.splits.

**Current** (`run_single_consistency_split()` returns `NA_real_` for
non-estimable splits; `evaluate_consistency_*` accumulates
`n_total_valid` and `n_total_success` and computes):

```r
p.consistency <- n_total_success / n_total_valid
```

Non-estimable splits are **excluded** from the denominator.

**Why this matters.**  For a subgroup where some splits are
non-estimable (common for small harm subgroups with few events per
arm), these two formulas differ.  Current systematically produces
**higher** Pcons than legacy on the same subgroup + same splits.

That sounds like it would help the current code find subgroups more
often — but the real effect is subtler.  In the GBSG example with
`sg_focus = "maxSG"`, the algorithm picks the LARGEST candidate
passing Pcons ≥ 0.9.  Under the current inflated-Pcons convention,
marginally small candidates can sneak past the 0.9 bar when they
would have failed in the legacy code.  For LOO folds (N ≈ 685) a
different, smaller-than-"er <= 0" subgroup can win on some folds,
driving up the number of distinct subgroups found across folds and
suppressing the exact-match rate.

**Status.**  Algorithmic change.  Not an obvious bug — both
conventions are defensible — but they are not equivalent.

### §1.3 — HIGH: Seeding within the split loop

**Legacy** (line 156):

```r
# If not estimable (both splits) then consider NOT consistent
set.seed(8316951)
for (bb in 1:n.splits) {
  ...
}
```

**The seed is fixed before each candidate's split loop.**  Across
different candidate subgroups, the sequence of splits is identical
(same seed → same `sample()` draws).  Two candidates of the same size
therefore see exactly the same partition sequence.

**Current.**  No fixed-seed-per-candidate.  The L'Ecuyer-CMRG streams
installed by `future::plan()` advance across all the random splits
done in `evaluate_consistency_fixed()` / `evaluate_consistency_twostage()`,
so each candidate gets a different split sequence.

**Why this matters.**  Legacy's paired-split design reduces the
variance in Pcons **differences** between candidates — a candidate
that edges out another in legacy will also edge it out on a slightly
different draw, because they share the random splits.  Current's
independent streams introduce extra variance into the ranking.  For
`sg_focus = "maxSG"` where Pcons merely has to clear a threshold and
size is the tiebreaker, this matters less; for `sg_focus = "hr"`
where HR and Pcons jointly determine the ranking, it matters more.

**Status.**  Deliberate architectural choice driven by
parallelisation.  Reintroducing a per-candidate seed would be cheap
and would tighten reproducibility against the legacy result.

### §1.4 — MEDIUM: `sg_focus` naming and preview-sort semantics

**Legacy** had five values: `"hr"`, `"Nsg"`, `"Msg"`, `"Nsg_only"`,
`"Msg_only"`.  The `_only` suffix controlled the **evaluation order**:

```r
# line 107-111 of subgroup_consistency_v0.R
if (!(sg_focus %in% c("Nsg_only","Msg_only")))
  found.hrs <- found.hrs[order(found.hrs$HR, decreasing = TRUE), ]
if (sg_focus == "Nsg_only")
  found.hrs <- found.hrs[order(found.hrs$n,  decreasing = TRUE), ]
if (sg_focus == "Msg_only")
  found.hrs <- found.hrs[order(found.hrs$n,  decreasing = FALSE), ]
```

So legacy `"Nsg"` = evaluate in **HR-desc** order (with early-stop
at `stop.threshold`), then select largest among those passing.
Legacy `"Nsg_only"` = evaluate in **size-desc** order (with the same
early stop), select largest passing.

**Current** collapses the two into `"maxSG"`, and `sort_subgroups_preview()`
orders evaluation as `setorder(result_new, -n, K)` — i.e.,
**size-desc**.  This matches legacy `"Nsg_only"`, not legacy `"Nsg"`.

**Why this matters.**  With `stop.threshold = 0.95` + size-desc
evaluation, the current `"maxSG"` stops at the first large-enough
candidate passing 0.95 — exactly what legacy `"Nsg_only"` did.  But
if the paper's GBSG run used legacy `"Nsg"` (which is plausible —
it's the non-`_only` default), the evaluation order would have been
HR-desc, and whichever HR-high candidate first cleared 0.95 would
have terminated the loop early.  Different early-stop point → could
pick a different subgroup than current `"maxSG"`.

*(Not fully verifiable by inspection.)*  Which of `"Nsg"` vs
`"Nsg_only"` the paper's GBSG analysis used is not clear from the
paper text, and it is not encoded in the uploaded files.  If the
paper used `"Nsg"`, current `"maxSG"` is a semantics change and a
`"maxSG_hrorder"` option might be needed to reproduce.

### §1.5 — MEDIUM: `subgroup.consistency()` runs all three foci in legacy

**Legacy** (`subgroup_consistency_v0.R` lines 218-227): after the
main candidate loop, runs `sg_consistency_out()` **three times**
(`sg_focus = "hr"`, `"Nsg"`, `"Msg"`) and populates `out_hr`,
`out_Nsg`, `out_Msg` on the return object.  Only one is selected as
`sg.harm`; the other two are there for inspection / comparison.

**Current.**  Runs `sg_consistency_out()` once with the focus_type
resolved via `sg_focus_map`.  No parallel comparison across foci.

**Why this matters for reproducibility.**  It does not, directly —
the selected `sg.harm` is determined by the user's `sg_focus` in
both.  But the legacy output had richer diagnostic content; if Larry
or any reader produced the paper tables by poking at `out_hr`,
`out_Nsg`, or `out_Msg`, they may have been looking at a focus
different from the one that `sg.harm` was populated from.  Worth
double-checking how the paper numbers were produced.

### §1.6 — MEDIUM: CV argument passing (hand-picked vs `args_call_all`)

**Legacy** (`forest_search_cross-validation_v0.R` lines 68-80, 462-479):
explicitly lists ~25 named arguments when calling `forestsearch()`
inside the per-fold foreach body.  Arguments NOT passed include
`m1.threshold`, `vi.grf.min`, `minp`, `seedit`, `grf_res`, and
`grf_cuts` — these inherit the function's own defaults.

**Current** (`forestsearch_cross_validation.R`): sets
`cv_args <- fs_args` where `fs_args <- fs.est$args_call_all`, then
applies four forced overrides (sequential plan, `details = FALSE`,
`quiet = TRUE`, `plot.sg = FALSE`, `ps_hat = NULL`).  **Every other
parameter propagates.**

**Why this matters.**  Mostly the current approach is *more* correct —
it preserves exactly what the user configured.  But two subtle
divergences:

1. *Post-paper parameters propagate automatically*.  If the user sets
   `use_twostage = TRUE` (default!) or any of the GLM / PS-adjustment
   args in the main call, they flow into CV folds too.  Legacy had
   no such options so it couldn't propagate them.  This is the lever
   behind the CV fixes you already applied (setting
   `use_twostage = FALSE` to match the paper).
2. *`grf_res` / `grf_cuts` hazard*.  The audit I wrote in
   `cv_code_audit.md` §2.1 flagged this.  In legacy, those arguments
   are never passed, so each CV fold re-runs GRF from scratch.  In
   current, they would propagate if the user had set them — a silent
   correctness hazard that is dormant for Larry's GBSG run (he
   doesn't pass them) but real for any user who does.

### §1.7 — LOW (user-flagged): `pstop_futile` gone

Legacy had `pstop_futile = 0.7` by default.  The loop stopped
evaluating further candidates the first time `Pcons <= 0.7` was hit,
relying on the HR-desc evaluation order to ensure remaining
candidates would not beat the threshold.  Current has no futility
stop; it evaluates up to `max_subgroups_search` candidates (default
10) instead.

**Impact.**  In the GBSG analysis with ~263 screened candidates,
legacy futility stop would typically kick in early (most candidates
Pcons < 0.7), saving ~80-90% of the compute but not changing the
*selected* subgroup (because selected subgroup passes Pcons ≥ 0.9
and is encountered before futility kicks in).  Current compensates
by capping at `max_subgroups_search = 10` evaluations.  Net: should
produce the same selected subgroup.

User has already pre-flagged this as OK.  Agreed — unlikely to
explain the CV discrepancy by itself.

### §1.8 — LOW: Parallelism layout (folds vs sims)

**Legacy `forestsearch_tenfold()`:** outer `for (ksim in 1:sims)`
loop is **sequential**.  The inner fold loop is parallel
(`%dofuture%`).

**Current `forestsearch_tenfold()`:** outer sim loop is **parallel**
(`%dofuture%`).  The inner fold loop is sequential.

**Why this matters.**  It does not — either layout produces the same
results.  The only consequence is scheduling: legacy uses workers =
`min(Kfolds, requested)` = typically 10; current uses workers =
`min(sims, requested)` = up to the full core count.  Current scales
better for many sims + few cores.

### §1.9 — LOW: `forestsearch_Kfold` validation checks

Both versions do the same two post-loop validations (every subject
appears exactly once, fold size ranges match).  The current version
also calls `forestsearch_KfoldOut()` internally to pre-compute
`sens_summary` / `find_summary` and attach them to the returned
object.  Purely additive, no behaviour change.

### §1.10 — INFORMATIONAL: post-paper feature additions

Tracked here for completeness; none of these are relevant to the
paper's GBSG CV numbers except insofar as they change defaults.

- **GLM outcome support** (`outcome_type`, `effect_measure`,
  `offset.name`, `adverse_outcome`, `overdispersion`,
  `grf_count_transform`).  Default `outcome_type = "survival"`; no
  effect on Larry's GBSG run.
- **Propensity-score adjustment** (`ps_method`, `ps_adjust_method`,
  `ps_hat`).  Default off; no effect unless used.
- **Two-stage consistency** (`use_twostage`, `twostage_args`).
  Default **TRUE** in current — this is the flip you've already
  addressed by setting it to FALSE in the QC document.
- **Early-stopping knob** (`stop_threshold`).  Legacy called the
  same thing `stop.threshold` (dot, not underscore).  Semantics
  equivalent.
- **`max_subgroups_search`**.  Hard cap on candidates evaluated;
  legacy had no such cap (only the futility stop + timing cap).
- **`vi.grf.min` variable-importance filter**.  Already present in
  legacy with the same default (-0.2).  No change.
- **Bootstrap bias correction / L_eff / FPR calibration machinery**.
  Entirely new file(s) (`bootstrap_*`, `fpr_calibration.R`,
  `calibrate_*`), not in the paper-era codebase.
- **`quiet` argument on `forestsearch()`**.  New in current.  Used
  by CV to suppress startup banners per fold.  No semantic effect.

---

## 2. Summary: which differences most plausibly explain Larry's CV gap

Ranking against the observed signal (current sens_H ≈ 51% vs paper's
73%; 60% fold-recovery rate vs paper's ≈ 90% at 10-fold):

| Rank | Difference | Direction of expected effect | Mechanism |
|------|------------|------------------------------|-----------|
| 1 | **§1.1 Bernoulli vs exact 50/50 split** | Reduces CV stability | Extra per-split variance on small subgroups destabilises "er <= 0" recovery across folds |
| 2 | **§1.3 Seed not fixed per candidate** | Reduces CV stability | Independent streams per candidate + random splits → less consistent ranking across sims |
| 3 | **§1.2 Pcons denominator convention** | Indirectly inflates non-intended candidates | Current's higher Pcons lets smaller candidates pass the 0.9 bar, reducing exact-match rate |
| 4 | §1.4 maxSG = Nsg_only vs Nsg semantics | Could pick different subgroup | Only if the paper used "Nsg" rather than "Nsg_only" — not verifiable from the paper text alone |
| 5 | §1.7 `pstop_futile` gone | Likely neutral | Current caps at `max_subgroups_search = 10` which typically matches the point where legacy's futility stop would trigger |

Items 1-3 together are, in my assessment, the most likely source of
the residual sens_H gap after the configuration fixes already
applied.  None is a bug.  They are deliberate post-paper design
choices that have changed the numerical characteristics of the
consistency evaluation.

---

## 3. Suggested path forward

None of these require a package code change right now.  Listed in
order of cost / benefit:

1. **Run the QC document with Larry's current config** (`use_twostage
   = FALSE`, `sg_focus = "maxSG"`, `Ksims_cv = 200`, `.oob_kfolds =
   nrow(df.analysis)`) and observe the actual gap against the paper.
   If the gap is small (≤ 5 pp on sens_H), the differences above are
   tolerable Monte Carlo noise and no action is needed.
2. **If the gap remains material (>5-10 pp)**, the cheapest next step
   is a debug chunk in the QC QMD that swaps out
   `run_single_consistency_split()`'s splitter for the exact-50/50
   variant via a function override, re-runs 10-fold CV, and reports
   the sens_H.  If the gap closes, §1.1 is confirmed as the dominant
   cause.
3. **For a proper v0.3.0 reproducibility option**, add
   `split_scheme = c("bernoulli", "exact5050")` and
   `per_candidate_seed = TRUE/FALSE` flags to
   `subgroup.consistency()` / `forestsearch()`, defaulting to the
   current behaviour but letting users opt into the legacy
   convention for paper-match runs.  The CV functions would
   propagate these automatically via `args_call_all`.
4. **Documentation**: a brief "reproducing León et al. (2024)"
   vignette that sets these flags and reports the expected paper
   numbers.  Useful for CRAN users who want to verify the published
   result.

---

## 4. Out of scope for this audit

- `get_fsdata.R` current vs `get_FSdata_v0.R`.  Structurally similar
  (LASSO + median/default cuts + GRF cuts integration + q-name
  dummification).  No differences I traced through that would affect
  GBSG paper reproducibility, but I did not read the current helper
  in full.
- `subgroup_search.R` vs `subgroup_search_v0.R`.  Same exhaustive
  combinatorial search up to `maxk`; current has minor refactoring
  (`remove_near_duplicate_subgroups()` and `sort_subgroups_preview()`
  added) but the candidate-enumeration logic is preserved.
- `forestsearch_functions_v0.R`.  `dummy`, `dummy2`, `FS_labels`, and
  helper cutters.  Current split across
  `forestsearch_helpers.R` and `get_FSdata_helpers.R`.  No material
  behaviour change observed in the areas I checked, but not
  exhaustively read.
