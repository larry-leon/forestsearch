# CC TASK — the structural-null counterpart of the GBSG survival simulations (identification only)

**Opened:** 2026-09-21 · **Repository:** forestsearch, branch `feature/glm-extension` ·
**Authorized by:** Larry, 2026-09-21 · **Machine:** pop-os (the survival campaigns' machine).
**Study directory:** the survival campaign behind the manuscript's eighteen cells — the manuscript
names it `gbsg_020` (main text p. 7); confirm the path from source in Step 0 and use it throughout.

**What this is.** The manuscript says twice that the design has no strict-null cell and that no
false-positive rate against the complete absence of an effect is measured (§5.1 p. 27, §5.4 p. 34).
This task builds that cell: the same GBSG-based survival design with **no planted region and a
uniform treatment benefit**, run for **identification and classification only**.

**Scope fence.** **No MR anywhere**: no multiplier resampling, no field, field-s, IJ or Bonferroni
products, no bootstrap, no cross-validation. **No edit to `R/`** — if the structural null cannot be
expressed through the campaign's existing parameters, that is a proposal for Larry, not part of this
task (Step 0 says what to do). No change to any of the eighteen committed cells, and nothing is
re-run or re-rendered there. No R CMD check, no vignette build, no test suite. **CC never fetches,
pulls or pushes.** Every `git add` names its paths; pre-existing untracked files are never staged.

**Fixed by Larry:** uniform hazard ratio **0.657** and **0.721** (the complement effects of the
narrow and broad attenuated-benefit cells, main text p. 27) · **n = 500, 1000, 1500** ·
**2,000 replicates** per cell · identifiers **FS, DINA, GRF** on identical draws within a cell ·
rule **`max_A N(ε)` at ε = 0.20** · all other search settings exactly the campaign's.

**Gates are stop-on-failure, never stop-to-ask.** A failing cell writes its halt record and the run
continues to the next cell.

---

## Step 0 — read the campaign from source, print it, decide feasibility

Read (do not modify) the study directory's template, runner and recorder, and print:

1. the template and runner paths, and every knob they expose (the `FS_*` environment knobs and their
   current defaults);
2. **where the planted region is constructed** — the covariate rule, its quantile knob, and the
   parameter that sets the interaction — and **where the marginal / complement treatment effect is
   set**; quote the lines with file:line;
3. the search settings the campaign actually uses: screening threshold, per-split threshold, `p⋆`,
   the selection rule and ε, the eligibility minimum, the per-arm event minimum, and **DINA's and
   GRF's proposal floors as set here** (the manuscript's p. 7 TODO says Table 1 gives DINA log 1.0
   while this campaign ran it at log 0.90 — record what the template sets, do not change it);
4. how inference is switched off for an identification-only run — the mode the selection-rule
   comparison of S3.1 used — and confirm that mode disables MR entirely;
5. what the recorder stores per replicate, so Step 3's list can be satisfied without new code;
6. the worker setting and seed scheme.

**Feasibility.** If a uniform-effect DGM with no planted region can be expressed through existing
parameters (interaction zero, treatment effect uniform), proceed. If it cannot without editing `R/`,
write a HALT record naming the file and line that would have to change, and stop — no compute.

## Step 1 — build the structural null, verified from the generator, not from description

Add the null design points at the template/campaign level only: **no planted region (interaction
zero), treatment effect uniform**, at HR 0.657 and 0.721. Before any replicate runs, print these
three checks from the generator's own truth object:

- the planted-region prevalence is zero or absent;
- over the enumerated candidate family on the super-population, `max_g |β(g) − log HR| ≤ 1e-8` — every
  candidate carries the uniform effect, which is what "no subgroup exists" means here;
- the super-population marginal hazard ratio equals the target to the same tolerance.

A failure of any of the three is a halt before compute.

## Step 2 — smoke, then continue without stopping

Run one cell (n = 500, HR 0.657, all three identifiers) at 20 replicates with MR off. Print the
per-replicate cost and the resulting projection for the full grid, then **continue straight into
Step 3** — the projection is recorded, not a decision point.

## Step 3 — the grid: 6 cells × 3 identifiers × 2,000 replicates

Cells in this order: `null0657_n500`, `null0721_n500`, `null0657_n1000`, `null0721_n1000`,
`null0657_n1500`, `null0721_n1500`. Identifiers run on identical draws within a cell, as the campaign
does. Commit each cell as it completes. Record per replicate, using fields the recorder already
provides (a field it does not provide is a template-level addition; it is never an `R/` change, and
if it would be, leave it out and say so in the report):

- whether a region was declared, and the identifier and rule that declared it;
- `|Ĥ|` and `|Ĥ|/n`, and the selected region's definition — its covariates, cuts and directions;
- the family counts: candidates enumerated, candidates clearing the floor, candidates passing the
  consistency screen;
- the consistency rate of the selected candidate and the maximum over the family;
- **`max_g T_g`**, the standardized screen statistic over the family — the quantity Section 4's
  calibrated cutoff is built from — if the recorder exposes it;
- the unadjusted within-region estimate, its robust standard error and its one-sided 95% Wald lower
  bound (these come from the search itself; they are not MR products);
- `sim_id`, seed, and the cell's wall clock.

If the campaign's machinery already records the other selection rules' picks from the same search, as
the S3.1 sweep does, record them too — they cost no extra search. Do not add a second search to get
them.

## Step 4 — the record

`REPORT_null_gbsg_identification_2026-09-21.md` in the study directory. Per-cell items as one table
(cells × identifiers), findings as bullets:

- **declaration rate** with a Wilson 95% interval — under this DGM it is the false-declaration rate,
  and it is the headline;
- `|Ĥ|` and `|Ĥ|/n`: mean with Monte Carlo standard error, and the quartiles;
- **specificity (benefit retained)** under **both conventions** the paper uses — unconditional, where
  a replicate declaring nothing scores 1, and conditional, over declaring replicates only (S3.4);
- the family counts, and the share of replicates where the consistency screen declined a region that
  had cleared the floor;
- the quantiles of `max_g T_g` against the conventional cutoff `z_0.95 = 1.645`, if recorded;
- the unadjusted within-region estimate: median and quartiles, and the share of declaring replicates
  whose unadjusted one-sided lower bound reaches **HR 1.00** and **HR 1.25**;
- the composition of Ĥ: the covariates and cut directions that appear, with their frequencies;
- measured wall per cell, and the machine, R version and package build.

State plainly in the report that **sensitivity and PPV are undefined with an empty planted region**
and are therefore not reported, that NPV is 1 by construction, and that the identified-to-planted
size ratio has no denominator here. Rates carry Wilson intervals and means carry Monte Carlo standard
errors computed across replicates, as S3.7 requires. Findings only beyond that.

Then the closeout: regenerate the study directory's `current_status.md` as the last action if that
directory carries one, under its existing pin convention. Render logs follow the directory's own
convention; any index-only git operation is a bare commit (no pathspec) with a staged-set assertion
before it and an `ls-tree HEAD` check after.

## Not in this task

MR and every interval product. The evaluation of Section 4's calibrated cutoff — the recorded
`max_g T_g` is what makes that possible later without a re-run. A boundary null at uniform HR 1.00 —
three more cells on this machinery, Larry's call, not this task. Any change to `R/`, to the eighteen
committed cells, or to the manuscript. Fetching, pulling, pushing.
