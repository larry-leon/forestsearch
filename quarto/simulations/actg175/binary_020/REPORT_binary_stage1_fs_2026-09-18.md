# REPORT — Binary campaign, stage 1: the six `orfs` cells at 1,000 replicates

Task: `dev/tasks/TASK_binary_launch_v2_2026-09-18.md` (committed as received, `04c5aaf1`), Steps 0–3.
Machine `pop-os` (128 physical cores), R 4.6.1, **forestsearch 0.3.5.9000, built 2026-09-19
03:54:02 UTC**. Branch `feature/glm-extension`. HEAD at launch `236ef82a`; HEAD at completion
`439c6363`. **Nothing pushed.**

**Design of record.** `sg_quantile = 0.62850`, prevalence(H) = 0.149170, from
`REPORT_binary_redesign_2026-09-18.md`. Every cell's `meta` carries both, and Gate 2 asserts them.

**Fixed by Larry, and honoured in every cell:** 1,000 replicates · 63 workers · MR only — the
evaluated set is unadjusted, oracle, IJ two-term, field, field-s, Bonferroni; no bootstrap, no
cross-validation anywhere (`nb_boots = NULL`, `fb_mode = "none"`, zero `forestsearch_bootstrap` /
`_Kfold` / `_tenfold` / `cv_*` calls in the template).

---

## 1. Gate 2, per cell

All six cells **PASS**, `run=74 passed=74 failed=0` in every one. The checker is
`scripts_or/gate2.R`; the per-cell records are `logs_or/GATE2_orfs_<cell>.txt` and are reproduced
verbatim in `REPORT_actg175_or_gate2_2026-09-17.md`.

| # | cell | target OR(H) | n | Gate 2 | declaration | sens | ppv | mean \|Ĥ\| | commit |
|---|---|---|---|---|---|---|---|---|---|
| 1 | `orfs_or150_n500`  | 1.50 | 500  | **74/74** | 0.9310 (931/1000) | 0.4068 | 0.3256 | 99.9  | `13d5b6c4` |
| 2 | `orfs_or075_n500`  | 0.75 | 500  | **74/74** | 0.7810 (781/1000) | 0.2401 | 0.1861 | 101.7 | `86f543b1` |
| 3 | `orfs_or100_n500`  | 1.00 | 500  | **74/74** | 0.8520 (852/1000) | 0.3027 | 0.2374 | 101.3 | `3ee0db52` |
| 4 | `orfs_or075_n2000` | 0.75 | 2000 | **74/74** | 0.9100 (910/1000) | 0.0888 | 0.1943 | 138.9 | `27170d87` |
| 5 | `orfs_or150_n2000` | 1.50 | 2000 | **74/74** | 0.9980 (998/1000) | 0.2805 | 0.5875 | 137.5 | `ed607c82` |
| 6 | `orfs_or100_n2000` | 1.00 | 2000 | **74/74** | 0.9610 (961/1000) | 0.1462 | 0.3188 | 139.9 | `6da4df13` |

Run in the task's order: `orfs_or150_n500` first — the configuration that halted twice on
2026-09-17 — then the other n = 500 cells, then the n = 2000 cells.

Gate 2 grew from 69 checks to 74. The five additions are the launch record Step 3 requires:
`n_sims == 1000`, the Stage 0 feasibility verdict carried in `meta` (`feasible` TRUE and **not**
overridden), `helper_identical`, the wall clock, and the estimability counts.

Two facts that hold in all six cells and are worth stating once: **0 CONFIG-ERROR replicates** and
**0 MR failures on the declared replicates** (the gate's ceiling is 40).

### 1.1 Bound location, for reading

Read against the ladder, never as significance at OR = 1 (`status_curated.md` §4).

| cell | mean field lower 1s on Ĥ | share ≥ 1.0 | mean field-s upper 1s on Ĥᶜ | share ≤ 1.0 |
|---|---|---|---|---|
| `orfs_or150_n500`  | 0.4108 | 0.020 | 0.9821 | 0.556 |
| `orfs_or075_n500`  | 0.3434 | 0.004 | 0.9213 | 0.703 |
| `orfs_or100_n500`  | 0.3634 | 0.007 | 0.9504 | 0.641 |
| `orfs_or075_n2000` | 0.3395 | 0.002 | 0.7748 | 0.997 |
| `orfs_or150_n2000` | 0.4917 | 0.044 | 0.8445 | 0.961 |
| `orfs_or100_n2000` | 0.3737 | 0.010 | 0.8045 | 0.991 |

---

## 2. Measured wall clock against the estimate

The task's estimate is linear scaling from the two measured 2,000-replicate cells
(`orfs_or075_n500` 2,248 s and `orfs_or075_n2000` 9,986 s, `REPORT_binary_redesign_2026-09-18` §6):
**1,124 s** per n = 500 cell and **4,993 s** per n = 2000 cell, so 3 × 1,124 + 3 × 4,993 =
**18,351 s ≈ 5.10 h** — the task's 5.1 h.

| cell | estimate | measured cell wall | ratio | batch render | combine | peak MB |
|---|---|---|---|---|---|---|
| `orfs_or150_n500`  | 1,124 s | **1,237 s** | 1.10 | 1,216 s | 21 s | 79,041 |
| `orfs_or075_n500`  | 1,124 s | **1,105 s** | 0.98 | 1,079 s | 25 s | 78,865 |
| `orfs_or100_n500`  | 1,124 s | **1,135 s** | 1.01 | 1,114 s | 20 s | 79,664 |
| `orfs_or075_n2000` | 4,993 s | **5,054 s** | 1.01 | 5,028 s | 25 s | 123,664 |
| `orfs_or150_n2000` | 4,993 s | **5,555 s** | 1.11 | 5,529 s | 25 s | 126,339 |
| `orfs_or100_n2000` | 4,993 s | **5,321 s** | 1.07 | 5,300 s | 20 s | 125,426 |
| **total** | **18,351 s** | **19,407 s** | **1.058** | — | — | — |

Cumulative render wall as the runner counts it (render time only, excluding gating and commits):
**19,402 s = 5 h 23 min**, against the `ORSG_CEILING=30000` s set for this stage. Elapsed
wall-clock from launch to `CAMPAIGN COMPLETE`: 04:33:39Z → 09:57:06Z = 5 h 23 min.

**The estimate holds to 5.8% in aggregate**, and the per-cell scatter is small and one-sided in a
readable way: the two OR 1.5 cells are the slowest relative to their estimate (1.10 and 1.11) and
both OR 0.75 cells are the fastest (0.98 and 1.01). The estimate was built entirely from OR 0.75
cells, so this is the design point's own cost, not noise — the harm point declares far more often
(0.931 and 0.998 vs 0.781 and 0.910), and each declaration buys an MR gate. Any projection onto the
OR 1.5 cells should carry the ~10% surcharge.

Peak memory is unchanged from the prior campaign's profile: ~79 GB at n = 500, ~123–126 GB at
n = 2000, on 63 workers.

---

## 3. The non-estimable / NA-oracle counts

Step 3 requires these per cell. They are computed in the template and carried in every bundle's
`meta`; Gate 2 asserts they are present and prints them.

- **NA-oracle** (`n_na_oracle_H`, `n_na_oracle_Hc`): replicates whose oracle refit on the **true**
  region returned the NA quadruple — the legacy pooled 5-events/5-non-events guard, or the
  four-cell existence condition added in `TASK_binary_study_redesign_2026-09-18` Step 3. Counted
  over all 1,000 rows.
- **Non-estimable** (`n_nonestimable_H`, `n_nonestimable_Hc`): **declared** replicates whose
  unadjusted plug-in on the **selected** region is not finite — the package's own estimator
  boundary (`.fs_existence_reason()`, `R/glm_effect_estimators.R`) where it can actually bind.

| cell | NA-oracle Ĥ | NA-oracle Ĥᶜ | non-estimable Ĥ | non-estimable Ĥᶜ | MR failures | CONFIG-ERROR |
|---|---|---|---|---|---|---|
| `orfs_or150_n500`  | 0 / 1000 | 0 / 1000 | 0 | 0 | 0 | 0 |
| `orfs_or075_n500`  | 0 / 1000 | 0 / 1000 | 0 | 0 | 0 | 0 |
| `orfs_or100_n500`  | 0 / 1000 | 0 / 1000 | 0 | 0 | 0 | 0 |
| `orfs_or075_n2000` | 0 / 1000 | 0 / 1000 | 0 | 0 | 0 | 0 |
| `orfs_or150_n2000` | 0 / 1000 | 0 / 1000 | 0 | 0 | 0 | 0 |
| `orfs_or100_n2000` | 0 / 1000 | 0 / 1000 | 0 | 0 | 0 | 0 |

**Every count is zero, in all 6,000 replicates.** This is the expected result, not a surprise, and
it should be read as confirmation rather than as evidence the machinery is idle: the redesign's
feasibility sweep found `share_nonestimable` **0.000 in all 60 rows** at every prevalence and every
n, and predicted that on this design the four-cell condition would essentially never fire on the
true region (`REPORT_binary_redesign_2026-09-18` §8, finding 4). The campaign confirms that at
1,000 replicates per cell. It does **not** demonstrate that the boundary works — the design never
puts it under load. The mean true-region size is 74.3 at n = 500 and 298.1 at n = 2000, against a
floor of 60, and the realized selected regions average 99.9–101.7 and 137.5–139.9; there is no tail
here thin enough to empty an arm × outcome cell.

The same is true of the declaration side: the Stage 0 gate recorded `share_undeclarable` 0.015 at
n = 500 and 0.000 above it, and no cell produced a CONFIG-ERROR.

---

## 4. What happened during the run, in full

### 4.1 The two `or075` cells could not run where their superseded bundles sat — resolved by a move

At 1,000 replicates a cell writes `<stem>_res_1_1000.rds`. For `orfs_or075_n500` and
`orfs_or075_n2000` **that exact path was git-tracked**: it is the first batch of the superseded
2,000-replicate run at the old 9.632% prevalence. The template's `.refuse_if_tracked()` guard —
correctly — refuses to overwrite a git-tracked bundle, so two of the six cells were hard-blocked
before the launch could reach them.

`status_curated.md` §2 records that those bundles are **not deleted**. They were therefore moved,
not removed: both directories were `git mv`'d whole from `mr_or_harm/` to
`mr_or_harm_superseded_prev09632/` (`81e671a9`). Git records all six files as renames; the payloads
are byte-identical and still in the worktree; one `git mv` back restores the previous layout
exactly. Independently of the write guard, the move takes them out of the template's combine glob,
where a 9.632% batch could otherwise pool with a 14.917% one — which §2 already forbids in prose
and which now holds mechanically. Their combine renders (`fs_…_orfs_combine_1_2000.html`) stay in
the directory root and are distinguishable by their `_1_2000` suffix.

**This is the one change to the repository that the task does not name, and it is Larry's to
reverse if he prefers the alternative** — running the design-of-record campaign under a fresh
campaign tag, which would leave the superseded bundles untouched but rename every cell the task
names, `orgrf_or150_n500` and `ordina_or150_n500` included.

### 4.2 A halt on cell 1 — a bug in the Gate 2 edit, not in the run

The first launch (04:11:45Z) halted at 04:32:27Z on `orfs_or150_n500` with
`GATE_COUNTS run=73 passed=72 failed=1`. The failing check was `checker ran without error` with
`object 'f1' not found`: the Step 1.5 edit (`56936b85`) replaced `gate2.R`'s hard-coded two-batch
file names with a vector derived from `ORSG_NSIMS`, but the payload-size loop at the tail of the
cell checker still named the removed `f1` / `f2`. That loop is the last block in the function, so
all 72 substantive checks had already run and passed.

Fixed in `557e173f`. Re-running the fixed checker against the bundles left on disk returned
**74/74** — the render itself was sound, 1,000 replicates in 1,242 s. The halt file was cleared in
`236ef82a`, following this directory's precedent (`6f80f292`), and the campaign relaunched at
04:33:39Z. The cell re-rendered rather than being hand-committed: its bundles were untracked, so
the runner's skip rule does not cover them, and a resume path would be exactly the retry logic the
task rules out. Cost ≈ 21 min. **No cell was committed on a failed or partial gate.**

### 4.3 Three Gate 2 entries carry render lines from the superseded 2026-09-17 run

The runner composes each cell's Gate 2 entry with
`grep -h '^WALL_SECONDS=' "$LOGD/${CELL}_"*.log`. Three cells — `orfs_or150_n500`,
`orfs_or075_n500`, `orfs_or075_n2000` — also ran on 2026-09-17 under the 2 × 1,000 layout, which
left `<cell>_batch_1001_2000.log` and `<cell>_combine_1_2000.log` in `logs_or/`. The glob picks
those up, so those three entries in `REPORT_actg175_or_gate2_2026-09-17.md` list **extra `Render:`
lines belonging to the old run** (e.g. `orfs_or075_n2000` shows `batch_1001_2000` at 4,939 s and
`combine_1_2000` at 25 s beside today's `batch_1_1000` at 5,028 s).

This is additive noise in one bullet of a report, nothing more: it does not touch any bundle, any
gate result, or any number in §2 above. The `Cell wall:` line and the heartbeat's `cell_wall_s` are
computed from timestamps, not from the glob, and are authoritative — they are what §2 reports. The
glob is the original runner's and predates this task's edits; the 2 × 1,000 → 1 × 1,000 change is
what exposed it. It was **not** fixed mid-run, because editing a running bash script risks
corrupting its execution, and the stage-2a cells carry different cell names so they cannot be
affected. See §6.

### 4.4 Two stale pins that would have failed Gate 2 on every cell

Found in Step 0 and fixed in `56936b85`, both env-overridable so the superseded cells stay
checkable: `gate2.R` asserted `meta$sg_quantile == 0.70` — the **old** design — and both `gate2.R`
and `smoke_identity.R` pinned `pkg_version 0.3.5`, which the Step 1.2 bump to 0.3.5.9000 had
already invalidated. The second of these surfaced as the only `[FAIL]` in the Stage 0 re-render.

---

## 5. Step 0 / Step 1 / Step 2, for the record

**Step 0** found no active run (no `R`/`Rscript`/`quarto` process, no `HALT_*`; the last runner
entry was the 2026-09-18T16:38:45Z start of `orfs_or150_n500`, which ended in the `OR positivity
violated -- non-positive finite values in: or_H_lo (1)` render halt). Already done: the version
bump (`66c4469b`, to 0.3.5.9000 — the repo's convention, not the task's provisional `.9001`), the
install (0.3.5.9000 built 2026-09-19 03:54:02 UTC; no `R/` commit postdates it), and the driver
revert (`23b9714d`, byte-identical to `625b10d0^`, sha256 `408b864b…8285ad`). Missing: the
replicates change, and a Stage 0 record on the current build.

**Step 1** completed only what was missing. Replicates → 1,000 as **one** batch over global
`sim_id` 1–1000 plus the combine (the seed table has `MAX_SIMS = 5000` and 1,000 distinct seeds in
`1:1000`, so coverage is exact and the draws stay identical across identifiers cell for cell);
`gate2.R` derives its batch layout from the same knob, so a superseded 2 × 1,000 bundle is still
checkable with `ORSG_NSIMS=2000`. The build and its provenance were recorded in `status_curated.md`
beside the superseded-cells note, `smoke_identity.R`'s `recipe` mode was recorded as **deferred**
with no code change, and §1's "one recorded exception" for the sweep driver was withdrawn — after
`23b9714d` the template is the single in-scope copy of `.logit_or_ci()`.

**Step 1.7 — Stage 0 re-render**, 105 s against a 300 s cap, campaign tag `relaunch`:
`feasible = TRUE` (undeclarable shares 0.015 / 0 / 0 / 0 against tolerance 0.05), `feas_override
FALSE`, `helper_identical TRUE`, `sg_quantile 0.62850`, prevalence 0.149170, pkg 0.3.5.9000, and
`SMOKE fs target_or_h=0.75 n=500: PASS` on every check including all field and field-s identities
(≤ 2.2e-16) and both γ intervals. 14 of 20 declared, matching the pre-install `redes` smoke exactly.

**Step 2 — the gate: all six items green**, and the campaign launched.

---

## 6. Open, no task attached

1. **The Gate 2 render-line glob** (§4.3). A one-line fix — glob only the labels the cell rendered,
   or clear stale per-cell logs at the start of a cell — would stop this recurring on any future
   resume across a layout change. Not done here: out of the task's scope, and it was unsafe to
   touch mid-run.
2. **The superseded-bundle move** (§4.1) is the one repository change the task does not name.
3. **The estimability counts are all zero and the design cannot make them otherwise** (§3). If the
   boundary is meant to be exercised rather than merely confirmed, that needs a design where the
   true or selected region can empty an arm × outcome cell — a smaller prevalence, a smaller n, or
   a more extreme OR — not this one.
4. **The OR 1.5 design point costs ~10% more wall clock** than the OR 0.75 cells the projection was
   built from (§2). The stage-2 projection should carry that, since 4 of the 12 remaining
   GRF/DINA cells are at OR 1.5.
5. **`ORSG_MRFAIL` is unchanged at 40** — the convention was written for 2,000-replicate cells, so
   at 1,000 it is twice as permissive in rate. It never bound: every cell recorded 0.
