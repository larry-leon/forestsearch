# CC TASK — campaign (B): the structural null at c1 0.90 / c2 0.80, with MR post-selection inference

**Opened:** 2026-09-21 · **Repository:** forestsearch, branch `feature/glm-extension` ·
**Authorized by:** Larry, 2026-09-21 · **Machine:** pop-os · **Study directory:**
`quarto/simulations/gbsg_020` · **Campaign tag:** `nullmr`.

**What this is.** The six `nullid` cells at `nullid`'s screen, now with MR inference. Identification
is `nullid`'s, since MR cannot change which subgroup is identified. What this campaign adds is where
the selection-adjusted bounds on a false region sit, and how they cover their targets.

**Fixed by Larry.**
- **Screen:** (c1, c2, p⋆) = (0.90, 0.80, 0.90), i.e. `FS_S7_C1` / `FS_S7_C2` unset, as `nullid`.
  DINA's floor is log 0.90; GRF's effect-scale floor is 0.90; `dmin.grf` is 0.0.
- **As `nullid`:** uniform super-population marginal HR 0.657 and 0.721 × n 500 / 1000 / 1500;
  2,000 replicates; FS, DINA and GRF on identical draws; `effMaxSG` at ε = 0.20; seeds
  `8316951 + sim_id`.
- **Targets:** β(Ĥ) and β(Ĥᶜ), i.e. `betaHhat_H` / `betaHhat_Hc`, the super-population marginal Cox HRs
  of the selected region and its complement. With no planted region there is no θ†(H) and no oracle
  refit on H. DINA and GRF results are labelled **conditional on the proposed family**.
- **Products evaluated:**
  - naive;
  - IJ two-term;
  - the field lower bound on β(Ĥ);
  - the field-s upper bound on β(Ĥᶜ);
  - the Bonferroni pair of the last two.

  **Not evaluated or reported:** the oracle, any winner-only or winner-floor variant, and the full
  bootstrap (`FS_S7_FB=none`).
- **MR settings:** per engine, exactly those the committed survival campaigns ran (Step 1).

**Scope fence.** **No edit to `R/`. No install or reinstall.** The template is untouched except for
the one guard Step 3 allows. No change to `nullid`, `nullc125` or any committed cell; nothing
committed is re-run. No R CMD check, no vignette build, no test suite. CC never fetches, pulls or
pushes. Every `git add` names its paths; pre-existing untracked files are never staged.

**Gates are stop-on-failure, never stop-to-ask, and assert invariants only** — never a guessed rate,
count or share. Coverage and `mr_ok` rates are printed, never gated. A failing cell writes its halt
record and the campaign continues at the next cell.

---

## Step 0 — preconditions (any failure stops the task before compute)

1. `hostname` is pop-os; branch `feature/glm-extension`; HEAD contains `c4baf796`; no tracked
   modifications; no R, Rscript or quarto process running.
2. Copy this document from `~/Downloads` to `dev/tasks/TASK_null_gbsg_mr_2026-09-21.md` and commit it
   alone.
3. The installed forestsearch build is not older than the last commit touching `R/`, as in
   `nullc125` §0. If it is older, stop — no install inside this task.

## Step 1 — MR settings and coverage definitions, from source

1. **MR knobs.** Per engine, read the knobs the template reads only on the MR path, as the committed
   campaigns set them:
   - FS: `scripts_p12x20/campaign_p12x20.sh`, cross-checked against cert20's recorded environment
     (`REPORT_cert20_2026-09-08.md:9-17`);
   - DINA: `scripts_dinamr/campaign.sh`;
   - GRF: `scripts_dinamr/grfmr.sh` and `grfmrC.sh`;
   - any MR literals in the template itself (draw count, for instance).

   Record each with path:line. A disagreement within an engine is a stop. Any knob that only adds
   winner-only or winner-floor products stays off; record where a committed campaign had it on.
   `FS_S7_FB=none`.
2. **Coverage.** Compute it exactly as the committed survival campaigns' reports compute it: name
   the script and the lines reused, and use them unchanged.

## Step 2 — driver and gates

1. **Driver:** `scripts_dinamr/nullmr.sh`, a transplant of `nullthr.sh`. Change only these:
   - `FS_S7_MR=TRUE`, plus the per-engine MR knobs from Step 1;
   - campaign tag `nullmr`;
   - threshold knobs unset;
   - timeouts: 6 h per render, 24 h for the campaign.

   Everything else is `nullthr.sh`'s, including the halt-and-continue behaviour and the
   commit-per-cell by named paths.
2. **Gate A, per run** (`nullmr_gateA.R`): `nullthr_gateA.R`'s invariants, except that its
   MR-absence checks invert:
   - every MR product is NA wherever nothing was declared;
   - wherever `mr_ok == 1`, the five evaluated products are finite;
   - every full-bootstrap product is NA;
   - the meta records c1 0.90, c2 0.80 and p⋆ 0.90, and the resolved thresholds equal them on every
     replicate;
   - `betaHhat_H` and `betaHhat_Hc` are finite on every declaring replicate.

   The `mr_ok` rate among declaring replicates is printed, not gated.
3. **Gate C, per cell:** `nullthr_gateC.R` as is: `itt_est` / `itt_se` identical across the three
   identifiers on all 2,000 rows.

## Step 3 — smoke, identity gate, then the grid without stopping

1. **Smoke:** `null0657_n500` and `null0721_n1500`, `sim_id` 1–20, three identifiers, MR on,
   through the real driver.
2. **Identity gate (same machine).** The `null0657_n500` smoke must be identical to `nullc125`'s
   committed MR-off render of the same replicates (`results/*_nullc125inertunset_quickrun_res_1_20.rds`)
   on every shared column except timing and the MR / field / IJ products, and on the truth object.
   This is the machine-checkable form of "MR cannot change identification".
3. **One allowed template change.** If the MR path fails under `FS_S7_DGM=null` only because a
   recorder line reads a planted-region truth that is NULL under the null:
   - guard that line add-only (NA under the null, unchanged under `alt`);
   - prove it inert on the `alt` path by rendering an `alt` cell at `sim_id` 1–20 before and after
     the edit on this machine (identical on every shared non-timing column);
   - commit it and record it.

   Anything beyond that is a stop.
4. **Projection.** Project the grid from both smoke cells, with the n-scaling measured rather than
   assumed. If every gate so far is green and the projection is under 16 h, run the grid straight
   away; otherwise record the projection and stop.
5. **Grid:** 6 cells × 3 identifiers × 2,000 replicates = 18 renders, in `nullid`'s cell order,
   committing each cell (its three bundles, renders and logs) as it completes.

## Step 4 — the record

`REPORT_null_gbsg_mr_2026-09-21.md` in the study directory. Tables first, each followed by a short
plain-language reading in bullets.

- **Coverage table**, cell × identifier rows, product columns:
  - the field lower bound on β(Ĥ);
  - the field-s upper bound on β(Ĥᶜ);
  - the Bonferroni pair, jointly;
  - naive and IJ two-term, each as the committed reports define it.

  Each entry carries a Wilson 95% interval and the count of declaring replicates it rests on.
- **Bound-location table**, same rows. For each lower-bound product (naive, IJ two-term, field, and
  the Bonferroni pair's lower member), give:
  - the median lower bound on the HR scale;
  - the share of declaring replicates whose bound reaches HR 1.00, and HR 1.25;
  - the same two shares as unconditional rates over all 2,000.

  Print `nullid`'s unadjusted within-region share (its `nv_H` refit) beside the naive row for
  reference. These are locations against 1.00 and 1.25, not tests.
- **Context per cell:**
  - median `itt_est` against the target, and median true β(Ĥ) and β(Ĥᶜ). The truths are marginal Cox
    HRs on uncensored potential outcomes, while the trial fits censored data (`nullc125` §5), so
    coverage is read beside that gap;
  - the `mr_ok` rate among declaring replicates.
- **Identification against `nullid`, recorded, not gated:** per cell and identifier, the number of
  replicates whose declaration, label or |Ĥ| differs from `nullid`'s bundle.
- **Wall per cell**, and the machine, R version and build.
- **Not reported, and why:**
  - the oracle, which is undefined with no planted region;
  - the winner variants;
  - the full bootstrap;
  - sensitivity and PPV, which are undefined here;
  - the identified-to-planted size ratio, which has no denominator.
- **Findings as bullets.**

## Step 5 — closeout, then the bundle for chat

1. Add §2.10 to `status_curated.md` for `nullmr`, and regenerate `current_status.md` as the last repo
   action, its pin equal to HEAD at commit time. `check_current_status.sh` needs zsh: check its
   conditions by hand. The local branch name and the gitignored-file lines are accepted, as in
   `cabbd3f3`.
2. Write `~/Downloads/bundle_nullmr_2026-09-21.zip` from HEAD, containing:
   - this document and the report;
   - the Step 1 record;
   - the gate, smoke and render logs;
   - the driver, cells files and scripts;
   - any template diff;
   - `git log --stat` for this task's commits;
   - `current_status.md`.

   No `.rds` or `.html`. Print its path and listing.

## Not in this task

Section 4's calibrated cutoff. A boundary null. MR at the 1.25 / 1.00 screen. Any change to `R/`, to
committed cells, or to the manuscript.
