# CC TASK — campaign (A): the structural null at c1 = 1.25, c2 = 1.00 (identification only)

**Revision:** v2, superseding v1 of the same date. One threshold pair, (1.25, 1.00), not two.

**Opened:** 2026-09-21 · **Repository:** forestsearch, branch `feature/glm-extension` ·
**Authorized by:** Larry, 2026-09-21 · **Machine:** pop-os · **Study directory:**
`quarto/simulations/gbsg_020` · **Campaign tag:** `nullc125`.

**What this is.** The `nullid` grid (`dev/tasks/TASK_null_gbsg_identification_2026-09-21.md`,
`REPORT_null_gbsg_identification_2026-09-21.md`) ran the structural null at the template's
c1 0.90 / c2 0.80. This task runs the same six cells at a stricter screen. Everything except the
thresholds is exactly `nullid`.

**Fixed by Larry.**
- **(c1, c2, p⋆) = (1.25, 1.00, 0.90)**, c1 and c2 passed explicitly. DINA's proposal floor follows
  c1 (`m_diff = log(1.25)`), as does GRF's effect-scale admission floor; `dmin.grf` stays 0.0.
- **As `nullid`:** uniform super-population marginal HR 0.657 and 0.721 × n 500 / 1000 / 1500;
  2,000 replicates; FS, DINA and GRF on identical draws; `effMaxSG` at ε = 0.20; seeds
  `8316951 + sim_id`; the design-point gate as `nullid` ran it, including its check-(3) fallback;
  every other setting the template's.

**Scope fence.** No MR anywhere: no multiplier resampling, no field, field-s, IJ or Bonferroni
products, no bootstrap, no cross-validation. **No edit to `R/`. No install or reinstall.** No change
to `nullid` or to any committed cell; nothing committed is re-run. No R CMD check, no vignette
build, no test suite. CC never fetches, pulls or pushes. Every `git add` names its paths;
pre-existing untracked files are never staged.

**Gates are stop-on-failure, never stop-to-ask, and they assert invariants only** — never a guessed
rate, count or share. Both of `nullid`'s gate bugs were guessed thresholds. A strict screen that
rarely declares is an expected outcome, not a failure. A failing cell writes its halt record and the
campaign continues at the next cell.

---

## Step 0 — preconditions (any failure stops the task before compute)

1. `hostname` is pop-os; branch `feature/glm-extension`; HEAD contains `cabbd3f3`; no tracked
   modifications (`git status --porcelain --untracked-files=no` empty); no R, Rscript or quarto
   process running.
2. Copy this document from `~/Downloads` to
   `dev/tasks/TASK_null_gbsg_thresholds_2026-09-21_v2.md` and commit it alone.
3. The installed forestsearch build is not older than the last commit touching `R/`
   (`git log -1 --format=%cI -- R/` against `packageDescription("forestsearch")$Built`). Record the
   version, build time and R version. If the install is older, stop — no install inside this task.

## Step 1 — threshold knobs in the template (add-only, inert at the default)

1. **Knobs.** `FS_S7_C1` and `FS_S7_C2`, defaulting to the literals at `:601-602` (0.90, 0.80), fed
   to the same `hr.threshold` / `hr.consistency` the call at `:1225` passes. Both or neither: setting
   one alone stops the render, and so does C2 > C1. p⋆ stays the literal.
2. **Stem.** A threshold token in the stem only when the knobs differ from the defaults, empty at the
   default, so no committed stem changes and no stricter-screen bundle can pool with `nullid`.
3. **Meta (add-only).** c1, c2, p⋆, `dmin.grf` and DINA's `m_diff`, closing the gap `nullid` found:
   no committed bundle records its thresholds.
4. **Recorder (add-only).**
   - The thresholds the fit actually resolved, per replicate, from the result's `args_call_all`; if the
     result does not carry them on an engine, leave NA and say so in the report.
   - `itt_est` / `itt_se`: the Cox fit of treatment alone on the whole trial (`.cox_hr_ci()`, as the
     oracle block uses), computed before the no-detection return so it exists on every replicate. It
     is Gate C's draw fingerprint at any declaration rate, and it gives each cell's realized ITT HR.
5. **Inertness, same machine, before and after.** Before editing, render `null0657_n500` at `sim_id`
   1–20 on all three identifiers from HEAD. After editing, render the same with the knobs unset, and
   again with them set explicitly to 0.90 / 0.80. All three must be identical on every shared
   non-timing results column and on the truth object. Give each render its own smoke campaign tag with
   `FS_S7_QUICKRUN=TRUE`; nothing is written to a committed stem.

Commit the template, the inertness bundles and a short inertness record in `scripts_dinamr/logs/`.

## Step 2 — driver and gates

1. **Driver.** `nullid.sh` is zsh with macOS calls (`sysctl`, `memory_pressure`, BSD `date -r`), and
   pop-os has no zsh. Port it to bash as `scripts_dinamr/nullthr.sh`, reading the threshold pair and
   the cells in run order from a cells file. Take pop-os's `render.sh` invocation and worker count from
   the `p12x20` driver that ran on this machine (`campaign_p12x20.sh`). Keep `nullid.sh`'s pinned knob
   set, its unsets (`FS_S7_Z1Q`, `FS_S7_ER_JCUTS`), `FS_S7_DGM=null`, `FS_S7_MR=FALSE`,
   `FS_S7_FB=none`, and its halt-and-continue behaviour.
2. **Gate A, per run:** `nullid_gateA.R`'s invariants, plus: the meta thresholds equal the run's
   knobs; the resolved thresholds equal them on every replicate where they are recorded; `itt_est` is
   finite on every replicate.
3. **Gate C, per cell:** `itt_est` and `itt_se` identical across the three identifiers on all 2,000
   rows.
4. **Timeouts:** 90 minutes per render; 6 hours for the campaign from the driver's start.

## Step 3 — smoke, then the grid without stopping

1. **Smoke:** `null0657_n500` at (1.25, 1.00), 20 replicates, three identifiers. Record the
   per-replicate cost and the grid projection.
2. If every gate so far is green and the projection for the grid is under 4 hours, run the grid
   straight away; otherwise record the projection and stop.
3. **Grid:** 6 cells × 3 identifiers × 2,000 replicates = 18 renders, in `nullid`'s cell order.
   Commit each cell — its three bundles and its render and gate logs — as it completes.

## Step 4 — the record

`REPORT_null_gbsg_thresholds_2026-09-21.md` in the study directory.

- **One table, screen × cell × identifier**, with `nullid`'s 0.90 / 0.80 rows read from its committed
  bundles (read only) so the two screens sit side by side:
  - declaration rate with a Wilson 95% interval — every declaration is false here;
  - mean |Ĥ|/n with Monte Carlo standard error;
  - specificity under both conventions (unconditional, where a replicate declaring nothing scores 1,
    and conditional);
  - median unadjusted within-region HR, and median true β(Ĥ) (`betaHhat_H`) where populated;
  - the share of declaring replicates whose unadjusted one-sided 95% lower bound reaches HR 1.00 and
    HR 1.25, and the same as an unconditional rate over all 2,000;
  - every conditional summary carries the count of declaring replicates it rests on.
- **FS:** the family counts, the share of replicates where the consistency screen declined a
  floor-clearing family, and the quantiles of `max_g T_g` against 1.645.
- **Per cell:** the median and quartiles of `itt_est`.
- **The composition of Ĥ** per identifier and screen.
- **Wall per cell**, and the machine, R version and build of each campaign, `nullid` included.
- **The statements of `nullid`'s §6:** sensitivity and PPV undefined with an empty planted region,
  NPV 1 by construction, no identified-to-planted size ratio.
- **Findings as bullets.**

## Step 5 — closeout, then the bundle for chat

1. Add §2.9 to `status_curated.md` for `nullc125`, and regenerate `current_status.md` as the last
   repo action, its pin equal to HEAD at commit time. `check_current_status.sh` needs zsh: check its
   conditions by hand, as in `cabbd3f3`. The local branch name and the gitignored-file lines of the
   inventory are accepted, as there.
2. Write `~/Downloads/bundle_nullc125_2026-09-21.zip` from HEAD: this document, the report, the
   inertness record, the gate and render logs, the driver, cells file and scripts, this task's template
   diff, `git log --stat` for this task's commits, and `current_status.md`. No `.rds` or `.html`.
   Print its path and listing.

## Not in this task

MR and every interval product — that is campaign (B), a separate task. Section 4's calibrated cutoff.
A boundary null. Any change to `R/`, to `nullid`, to any committed cell, or to the manuscript.
