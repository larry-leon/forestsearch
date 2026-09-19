# CC TASK — binary campaign stage 2: the ten GRF/DINA cells (v3 — Mac Studio; run to completion)

**Opened:** 2026-09-19 · **Repository:** forestsearch · **Machine:** Mac Studio · **Authorized by:** Larry, 2026-09-19.
**Supersedes `TASK_binary_stage2_v2_2026-09-19.md` and `TASK_binary_stage2_2026-09-19.md` — delete any earlier
copy from `dev/tasks/` and `~/Downloads` and commit this one in its place.**
**Study:** `quarto/simulations/actg175/binary_020/`.
**Predecessors (run on pop-os; must be at HEAD in this clone):** the six `orfs` cells (Gate 2 74/74) and the two
stage-2a timing cells (`orgrf_or150_n500` 937 s, 87/87; `ordina_or150_n500` 1,359 s, 88/88).
**v3 is v2 moved from pop-os to this Mac. Only the machine-specific elements changed; the design, cells,
estimators and order did not.**

**This task has no gates, no thresholds to clear, no checkpoints and no decision points. It runs all ten
cells to completion in one unattended pass and reports at the end.** The only way it stops before compute is a
Step 0 failure — this clone not carrying the campaign state (0a), the install failing (0c), or a Linux-only
construct that cannot be adapted (0d) — which writes a HALT record naming the cause. There is no R CMD check,
no vignette build, no test suite, and no re-running or re-verifying of any committed cell.
**No edit to `R/`.** `run_or.sh` and `gate2.R`: Step 0d's machine adaptations only. **CC never fetches, pulls
or pushes.** Every `git add` names its paths; pre-existing untracked files are never staged.

**Fixed by Larry:** replicates **1,000** per cell · **MR only** — unadjusted, oracle, IJ two-term, field,
field-s, Bonferroni; no bootstrap, no cross-validation anywhere · workers **W = this Mac's physical cores − 1**
(`sysctl -n hw.physicalcpu` − 1; 13 on this 14-core Mac Studio, where 13 beat 10 in mdf1) — the same rule that
gave 63 on pop-os's 64.

**The cumulative render ceiling must not be able to stop this run: set it to 90,000 × 63 / W s, rounded up to
the next 10,000 s (W = 13 → 440,000 s)** — pop-os's 90,000 s scaled by the worker ratio. Every other time cap
the runner applies (per-render timeouts, heartbeat staleness limits) scales by the same 63 / W. The caps exist
only so a pathological run is not unbounded; Larry stops the run himself if he wants it stopped.

---

## Step 0 — prepare this Mac, discover state, print it, launch

**0a. Campaign state.** This clone must carry pop-os's stage-1 state: the eight predecessor cells committed with
Gate 2 PASS, and the design of record in the runner configuration and the stage-1 bundle meta — `sg_quantile`
0.62850 (prevalence 0.149170), 1,000 replicates, MR only (`fb_mode` none), `seed_base` 8316951. Anything
missing or different means the clone is not synced (Larry pushes from pop-os and pulls here; CC does not pull):
write a HALT record naming what is missing and stop, before any install or compute.

**0b. Discover.** Which of the ten cells are already committed with Gate 2 PASS — skip those. A cell whose bundle
is complete on disk but uncommitted: render its Gate 2 and commit it; do not re-run it. Any active run: monitor
it to completion, then continue with whatever of the ten remains. Any `HALT_*`: clear it following this
directory's precedent and proceed.

**0c. Install forestsearch — always, except under an active run.** The Mac's installed build is not the
campaign's by construction, and a dev version string (0.3.5.9000) does not identify a build. Let S1 = the newest
of the eight stage-1 cell commits. If the package source is unchanged from S1 to HEAD
(`git diff --quiet S1 HEAD -- R DESCRIPTION NAMESPACE src inst data`), the install commit C = HEAD; otherwise
C = S1, so all eighteen cells run one package source, and the diff goes in the report as a finding. Install from
a temporary `git worktree` at C with `devtools::install(<worktree>, upgrade = "never")`, then remove the worktree
— no uncommitted change in this clone can enter the build, and a non-interactive install cannot upgrade
dependencies. Print the installed version, build time and C. Never install while a run is active (0b): the pass
keeps the build it started with.

**0d. Machine adaptations.** Read `run_or.sh`, the heartbeat, `gate2.R`, the cell templates and anything they
call, and print every pop-os-specific element:

- the worker count (63 / 64) and every time cap — set from W as above, through the runner's existing knobs;
- Linux-only constructs, which fail on macOS: GNU `timeout` (macOS has none), `nproc`, `sed -i` without a suffix
  argument, `date +%N`, `stat -c`, `readlink -f`, `setsid`, `flock`, `/proc`, `free`, `grep -P`,
  `find -printf`, absolute Linux paths, and bash ≥ 4 features under macOS's `/bin/bash` 3.2.

Resolve each through the environment first: a shim directory outside the repo, prepended to `PATH` for this pass
— `timeout` via `gtimeout` if Homebrew coreutils is installed, otherwise a Perl `alarm` wrapper (the approach
this Mac's earlier campaigns used), exiting 124 on expiry as GNU `timeout` does; an installed Homebrew `bash` if
bash ≥ 4 is needed. Where the environment cannot supply it, or a value has no knob, a minimal portability edit to
`run_or.sh` or `gate2.R` is permitted: identical behaviour on pop-os, `bash -n` / `parse()` clean, its own commit
before any compute, diff in the report. Nothing else in either file changes, and `gate2.R`'s render-line glob
stays untouched. A construct that can be neither supplied nor edited that way is a Step 0 HALT.

For the whole pass: `VECLIB_MAXIMUM_THREADS=1` in the runner's environment (Apple Accelerate is not fork-safe
after dense work — recorded on this Mac), and a `caffeinate -is` assertion held for the life of the runner
(e.g. `caffeinate -is -w <runner pid>`) so idle sleep cannot suspend it.

Seeds: confirm from the runner source that a replicate's seed and data depend only on `seed_base` and the
replicate index, not on the worker count, and print the lines. If they depend on W, run anyway and record it —
the Mac cells then do not share DGM draws with the pop-os cells.

**0e. Print, then launch.** Machine (model, P/E cores, memory), macOS and R versions, BLAS/LAPACK, the
installed build and C, MR only, W, the ceiling and the scaled caps, and the adaptations made in 0d.

## Step 1 — run the ten cells in this order, without interruption

The table gives the order and the **pop-os reference**: the measured `orfs` wall at the same design point and
n × the identifier ratio measured at `or150_n500` (GRF 0.758, DINA 1.099), all at 63 workers on pop-os. **It is
not a Mac prediction.** At W = 13 the same work is 63 / 13 = 4.85× the load per worker — ≈ 163,500 s (≈ 45 h)
if a Mac core matched a pop-os core; no Mac wall has been measured for this workload. Print the table, then run.
**A cell that fails writes its halt record and the run continues to the next cell** — per-cell stop-on-failure
only, exactly as the existing runner behaves; no new retry logic, and one cell's failure does not end the pass.

| order | cell | orfs wall (s) | × ratio | pop-os projection (s) |
|---|---|---|---|---|
| 1 | `orgrf_or075_n500` | 1,105 | 0.758 | 837 |
| 2 | `orgrf_or100_n500` | 1,135 | 0.758 | 860 |
| 3 | `orgrf_or075_n2000` | 5,054 | 0.758 | 3,831 |
| 4 | `orgrf_or100_n2000` | 5,321 | 0.758 | 4,033 |
| 5 | `orgrf_or150_n2000` | 5,555 | 0.758 | 4,211 |
| 6 | `ordina_or075_n500` | 1,105 | 1.099 | 1,214 |
| 7 | `ordina_or100_n500` | 1,135 | 1.099 | 1,247 |
| 8 | `ordina_or075_n2000` | 5,054 | 1.099 | 5,554 |
| 9 | `ordina_or100_n2000` | 5,321 | 1.099 | 5,848 |
| 10 | `ordina_or150_n2000` | 5,555 | 1.099 | 6,105 |

**pop-os total 33,740 s** — GRF's five first (13,772 s), so a complete second identifier is banked early, then
DINA's five (19,968 s). GRF's kept family is tight at n = 500 (K 1,036–1,298); DINA's spans K = 1 to 6,075 with
median 873 and p90 3,881, and MR cost scales with K, so DINA's n = 2000 cells are where a wall may run long.
**That is a number to report, not a reason to stop.**

Existing runner, heartbeat and Gate 2 conventions unchanged apart from 0d. Every bundle's meta carries version,
build time, prevalence, replicates, workers, wall clock, and the non-estimable / NA-oracle counts. Commit each
cell as it completes, as stage 1 did. Record swap in use (`sysctl vm.swapusage`) at the end of each cell — this
Mac has 36 GB and has run out of application memory before; swap growth is a finding, not a stop.

**The pass outlives this CC session.** Launch the runner detached (`nohup`, background) under the `caffeinate`
assertion and poll sparsely — about every 30 minutes, printing only the heartbeat's latest line — so a
multi-day monitor does not exhaust the session. If the session ends mid-pass, the same kickoff in a fresh
session resumes it through 0b.

## Step 2 — report, then four housekeeping commits

`REPORT_binary_stage2_2026-09-19.md` in the study directory — per-cell items as one table (cells × items),
findings as bullets:

- **Provenance, once.** The eighteen cells now span two machines: the eight predecessors on pop-os at 63
  workers, these ten on this Mac at W. Record machine (model, P/E cores, memory), macOS, R version, BLAS/LAPACK,
  `sessionInfo()` with forestsearch loaded, the installed build and C (0c), and each 0d adaptation with its
  diff; pop-os's values where the committed stage-1 record carries them.
- **Per cell:** Gate 2 result; measured wall beside its pop-os projection (reference only — their ratio mixes
  machine and identifier); declaration rate stated against the cell's design point; the kept-family K
  distribution; the non-estimable / NA-oracle / MR-failure / CONFIG-ERROR counts; swap at cell end.
- **Identifier cost at n = 2000, stated machine-free:** for GRF and DINA at `or075` and `or100`, the n-scaling
  (n = 2000 wall ÷ n = 500 wall, both measured on this Mac) against FS's on pop-os (4.57 at `or075`, 4.69 at
  `or100`). Equal n-scaling means the n = 500 identifier ratio held at n = 2000. This replaces v2's direct ratio
  against the `orfs` walls, which would now mix machine and identifier.
- Any coverage figure for GRF or DINA is coverage of the conditional-on-proposed-family estimand and says so.
  Findings only beyond that.

Then, after the compute and the report, four record-keeping commits:

1. Move every 2026-09-17 `logs_or/` file on this machine into `logs_or/superseded_prev09632/` (`git mv` if
   tracked, `mv` if not), so the Gate 2 render-line glob can only match the current run. **Do not edit the
   glob** — editing `gate2.R` is what halted cell 1 of stage 1. Untracked 2026-09-17 logs that exist only on
   pop-os are out of reach here: record that; do not reconstruct them.
2. `git rm --cached` the tracked `logs_or/` files per §3 of the study convention, plus the `.gitignore` line —
   a bare commit (no pathspec), with a staged-set assertion before it and an `ls-tree HEAD` check after. Files
   stay on disk.
3. Banner `REPORT_actg175_or_gate2_2026-09-17.md`: it interleaves both designs under identical headings —
   `orfs_or075_n500` and `orfs_or075_n2000` each appear twice. Name the superseded sections and the two
   discriminating fields (`pkg_version`: 0.3.5 superseded, 0.3.5.9000 current; and `prev`). Record text only;
   delete or rewrite nothing.
4. Regenerate the directory's `current_status.md` as the last action, per the standing convention.

**If any of these fails, record it and stop attempting it. None of them may hold up or reverse the compute,
which is already committed by then.**

## Not in this task

The 18-cell bias-and-coverage synthesis — a separate task once these ten land. Any edit to `R/`. Any edit to
`run_or.sh` or `gate2.R` beyond 0d. Any change to floors, thresholds, the DGM of record or replicates, or to the
worker count beyond the W rule above. Fetching, pulling, pushing.
