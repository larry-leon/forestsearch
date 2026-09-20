# REPORT — ACTG175 binary / OR stage 2: the ten GRF/DINA cells (Mac Studio)

Task: `dev/tasks/TASK_binary_stage2_v3_2026-09-19.md`. Opened 2026-09-19, run to completion in one
unattended pass on this Mac Studio; **all ten cells Gate 2 PASS**, no halt, no retry, no cell re-run.
Compute 2026-09-19T20:26:20Z – 2026-09-20T00:39:35Z, driver cumulative wall **15,195 s (4 h 13 m)**.

The eighteen cells of the binary/OR campaign now span two machines: the eight predecessors on `pop-os`
at 63 workers, these ten on this Mac at W = 13. The 18-cell bias-and-coverage synthesis is a separate
task.

**Every coverage figure for `orgrf` and `ordina` is coverage of the estimand CONDITIONAL ON THE
PROPOSED FAMILY.** GRF's and DINA's candidate families are generated from fitted surfaces, so the
fixed-family condition does not hold; FS's family is the prespecified cut grid. Nothing below is a
coverage figure — the synthesis task carries those — but the condition travels with the bundles.

---

## 1. Provenance (once, for all ten cells)

| | this Mac | `pop-os` (the eight predecessors) |
|---|---|---|
| machine | Mac16,9 Mac Studio, Apple M4 Max, 14 physical cores (10 P + 4 E), 14 logical, 36 GiB | 128 physical cores (as the committed stage-1 and stage-2a records state) |
| OS | macOS 27.0 (26A428), arm64 | — (not carried by the stage-1 record) |
| R | 4.5.2 (2025-10-31), platform `aarch64-apple-darwin20` | 4.6.1 |
| BLAS | Accelerate `vecLib .../libBLAS.dylib` | — (not carried by the stage-1 record) |
| LAPACK | `R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib` | — (not carried by the stage-1 record) |
| forestsearch | 0.3.5.9000, Built `R 4.5.2; ; 2026-09-19 20:23:40 UTC; unix` | 0.3.5.9000, built 2026-09-19 |
| workers | 13 | 63 |
| nodename | `Mac-Studio-3.local` | `pop-os` |

`sessionInfo()` with forestsearch loaded: R version 4.5.2 (2025-10-31); Platform `aarch64-apple-darwin20`;
Running `macOS 27.0`; BLAS `/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib`;
LAPACK `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib`; forestsearch 0.3.5.9000.

**Installed build and C (0c).** S1 = `adfd783f` (`ordina_or150_n500`, the newest of the eight stage-1 cell
commits). `git diff --quiet S1 HEAD -- R DESCRIPTION NAMESPACE src inst data` was **empty**, so **C = HEAD =
`b3d4921d`** and no package-source diff is reported as a finding. Installed from a temporary `git worktree`
at C with `devtools::install(upgrade = FALSE, dependencies = FALSE)`; the worktree was removed afterwards.
The build it replaced was **0.3.5, Built 2026-09-11 14:53:59 UTC** — not the campaign's, exactly as the task
anticipated. `devtools` 2.5.2 here rejects `upgrade = "never"` (it wants a single TRUE/FALSE/NA); `FALSE` is
the same instruction and was used.

**Design of record, verified in this clone before any compute (0a).** `sg_quantile` 0.62850
(template `sim_fs_mr_field_or_template.qmd:307`; bundle meta 0.6285), prevalence(H) 0.149170, 1,000
replicates, `fb_mode` `none` (MR only), `seed_base` 8316951, scheme "pre-generated table indexed by global
sim_id". Eight predecessors present and committed with Gate 2 PASS.

**Seeds do not depend on W.** `sim_fs_mr_field_or_template.qmd:155` states it, and `:164–172` build
`SEED_TABLE <- sample.int(...)` from `set.seed(seed_base)` alone, read by `seed_for(sim_id)` on the global
replicate id. The Mac cells therefore share DGM draws with the pop-os cells; §4 verifies this per replicate
rather than resting on the construction.

### 1.1 Machine adaptations (0d)

**No edit was made to `run_or.sh` or `gate2.R`** — every pop-os-specific element was resolved in the
environment, so there is no pre-compute portability commit and no diff to show. `gate2.R`'s render-line glob
was not touched. The adaptations:

| element | pop-os | resolution on this Mac |
|---|---|---|
| `timeout` (`run_or.sh:116`) | GNU coreutils | **absent on macOS, and no Homebrew coreutils/`gtimeout` installed** → Perl `alarm` wrapper on a shim `PATH` (see below) |
| `pgrep -fc` (`run_or.sh:133`) | procps-ng `-c` | **BSD `pgrep` has no `-c`** → shim counts lines and mirrors GNU's exit status; without it the preflight fails outright |
| `PATH` prepend (`run_or.sh:68`) | `/usr/lib/rstudio/.../quarto/bin` | directory does not exist here; inert. quarto resolves at `/usr/local/bin/quarto` (1.10.18). The shim dir stays ahead of `/usr/bin` |
| `ps -eo rss,comm` (`mem_sampler.sh:8`) | Linux `ps` | accepted by macOS `ps` (`-e` is `-A`); verified, no change |
| bash | ≥ 4 | `/bin/bash` 3.2.57 runs `run_or.sh` as written; `bash -n` clean. Empty-array iteration under `set -u` would fail on 3.2, but the guard uses `${#a[@]}` (fine on 3.2) and neither array is empty in this pass. No Homebrew bash needed |
| `nproc`, `sed -i`, `date +%N`, `stat -c`, `readlink -f`, `flock`, `/proc`, `free`, `grep -P`, `find -printf` | — | **none executed**; `setsid` appears only in a header comment |
| worker count | 63 | **W = 13** = `sysctl -n hw.physicalcpu` (14) − 1. The template's own cap `.phys_minus1` (`:142–143`) is also 13, so `FS_OR_WORKERS=13` is taken verbatim |
| `hostname` meta assertion | default `pop-os` | `ORSG_HOST=Mac-Studio-3.local` exported into the runner, inherited by `gate2.R` (`gate2.R:29`, asserted at `:163`) |
| Accelerate | n/a | `VECLIB_MAXIMUM_THREADS=1` in the runner environment (and again at `run_or.sh:69`); the future plan is `multisession` (socket, not forked) |
| idle sleep | n/a | `caffeinate -is -w <driver pid>` held for the life of the run |

The shim directory is `~/.orsg_shim_2026-09-19`, prepended to `PATH` for the pass; its contents are
committed under `scripts_or_mac/` (§6). The `timeout` shim is a stricter relative of this Mac's earlier
`../../continuous/scripts_mdf1/tmo.sh`: the child runs in its own process group, expiry sends TERM to the
group, waits a 10 s grace, then KILLs and exits **124** as GNU `timeout` does; otherwise the child's status
is propagated (128+N for a signal death). Tested before launch for rc 0, rc 3, rc 124 on expiry, and for
leaving no straggler behind a nested child.

### 1.2 Time caps

The cumulative render ceiling is **440,000 s** = 90,000 × 63/13 = 436,154, rounded up to the next 10,000 s.
The per-render timeout is **52,339 s** = stage 2a's 10,800 × 63/13. `run_or.sh` applies no other time cap —
it has no heartbeat staleness limit. Nothing came near either: the longest single render was 2,054 s.

**Supplement, applied mid-pass:** DINA's five cells (6–10) were given a per-render timeout of **157,017 s**
(3 × 52,339); GRF's five kept 52,339 s. `run_or.sh` reads `ORSG_TIMEOUT` once per invocation and may not be
edited, so the override lives in the `timeout` shim, which is exec'd afresh for every render and so applied
without restarting anything. It identifies the campaign from the render's own environment (`FS_OR_METHOD`,
set by `run_or.sh:112–117` as a prefix assignment on the same simple command as the wrapper) or the
`_ordina_` output stem, **and** requires the wrapped command to be `quarto render`, so no other use of the
wrapper is re-capped. Every render appends the passed and applied caps to `timeout_applied.log`, committed
in §6: **10 GRF renders at 52,339 s and 10 DINA renders at 157,017 s**, plus two pre-launch test lines. One
campaign render is missing from that log — cell 1's batch render ran at 20:26, before the audit line was
added at 20:47; it ran under the original shim at 52,339 s, which is what GRF gets in any case.

---

## 2. Per cell

Walls are `cell_wall_s` from the heartbeat. The **pop-os projection** is the task's table (that cell's
measured `orfs` wall × the `or150_n500` identifier ratio, GRF 0.758 / DINA 1.099, at 63 workers) and is a
reference only: its ratio to the measured wall mixes machine and identifier and is not an identifier result.
K is MR's kept family. "non-est." is non-estimable on the selected region among declared replicates.

| # | cell | Gate 2 | wall (s) | pop-os proj. (s) | ratio | declaration | K: min / med / p90 / max | NA-oracle H/Hc | non-est. H/Hc | MR fail | CONFIG-ERROR | swap at cell end |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `orgrf_or075_n500` | **PASS 87/87** | 1,428 | 837 | 1.71 | 0.9990 (999/1000) | 1,036 / 1,066 / 1,146 / 1,298 | 0 / 0 | 0 / 0 | 0 | 0 | 0.00M |
| 2 | `orgrf_or100_n500` | **PASS 87/87** | 1,444 | 860 | 1.68 | 0.9990 (999/1000) | 1,036 / 1,066 / 1,146 / 1,298 | 0 / 0 | 0 / 0 | 0 | 0 | 0.00M |
| 3 | `orgrf_or075_n2000` | **PASS 87/87** | 1,954 | 3,831 | 0.51 | 0.9990 (999/1000) | 1,339 / 1,351 / 1,447 / 1,552 | 0 / 0 | 0 / 0 | 0 | 0 | 0.00M |
| 4 | `orgrf_or100_n2000` | **PASS 87/87** | 1,963 | 4,033 | 0.49 | 1.0000 (1000/1000) | 1,339 / 1,351 / 1,447 / 1,552 | 0 / 0 | 0 / 0 | 0 | 0 | 0.00M |
| 5 | `orgrf_or150_n2000` | **PASS 87/87** | 2,070 | 4,211 | 0.49 | 1.0000 (1000/1000) | 1,339 / 1,351 / 1,447 / 1,552 | 0 / 0 | 0 / 0 | 0 | 0 | 0.00M |
| 6 | `ordina_or075_n500` | **PASS 88/88** | 1,330 | 1,214 | 1.10 | 0.9440 (944/1000) | 1 / 427 / 2,292.8 / 6,112 | 0 / 0 | 0 / 0 | 0 | 0 | 0.00M |
| 7 | `ordina_or100_n500` | **PASS 88/88** | 1,539 | 1,247 | 1.23 | 0.9680 (968/1000) | 1 / 540 / 2,842 / 6,082 | 0 / 0 | 0 / 0 | 0 | 0 | 0.00M |
| 8 | `ordina_or075_n2000` | **PASS 88/88** | 733 | 5,554 | 0.13 | 0.8700 (870/1000) | 1 / 130 / 614.2 / 2,836 | 0 / 0 | 0 / 0 | 0 | 0 | 0.00M |
| 9 | `ordina_or100_n2000` | **PASS 88/88** | 991 | 5,848 | 0.17 | 0.9340 (934/1000) | 1 / 260.5 / 979 / 4,439 | 0 / 0 | 0 / 0 | 0 | 0 | 0.00M |
| 10 | `ordina_or150_n2000` | **PASS 88/88** | 1,690 | 6,105 | 0.28 | 0.9870 (987/1000) | 3 / 662 / 2,034 / 6,987 | 0 / 0 | 0 / 0 | 0 | 0 | 0.00M |

Declaration is stated against each cell's design point. All three run with `adverse_outcome = TRUE` and a
0.90 harm-screen threshold on the OR scale, so `or075` is a true OR of 0.75 in H — a benefit, below the
screen; `or100` is the null at 1.0, already above the screen; `or150` is planted harm at 1.5 (present here
only at n = 2000, its n = 500 cell being a stage-2a predecessor). GRF declares at 0.999–1.000 in every cell,
including the benefit point `or075`. DINA declares at 0.870–0.987, falling with n and rising with the design
point at n = 2000 (0.870 → 0.934 → 0.987 across `or075` → `or100` → `or150`).

Peak memory (summed R/quarto RSS, `mem_sampler.sh`): **16.8–18.1 GB** across all ten batch renders on a
36 GiB machine, against ~79 GB (n = 500) and ~123–126 GB (n = 2000) at 63 workers on pop-os. **Swap stayed
at 0.00M for the whole pass**, measured at every cell end.

---

## 3. Identifier cost at n = 2000, stated machine-free

n-scaling = (n = 2000 cell wall) ÷ (n = 500 cell wall), **both measured on this Mac**, against FS's on
pop-os. Equal n-scaling would mean the n = 500 identifier ratio held at n = 2000.

| identifier | at `or075` | at `or100` |
|---|---|---|
| FS (pop-os, 63 workers) | 4.57 | 4.69 |
| GRF (this Mac, 13 workers) | **1.37** (1,954 / 1,428) | **1.36** (1,963 / 1,444) |
| DINA (this Mac, 13 workers) | **0.55** (733 / 1,330) | **0.64** (991 / 1,539) |

Neither identifier's n = 500 ratio held at n = 2000, and both moved the same way: GRF's n-scaling is about
**one third** of FS's, DINA's about **one eighth**, and DINA's is **below 1** — its n = 2000 cells ran
*faster* than its n = 500 cells. Quotients of two walls on one machine, so the machine largely cancels; what
does not cancel is any difference in how the two machines scale from n = 500 to n = 2000 (memory bandwidth,
worker count), so read these as identifier comparisons on this Mac rather than as pop-os numbers.

The K column in §2 accounts for the direction. MR's cost scales with the kept family K, and **DINA's K falls
as n rises** — median 427 → 130 at `or075` and 540 → 260.5 at `or100`, with p90 falling 2,292.8 → 614.2 and
2,842 → 979 — while GRF's K rises only modestly (median 1,066 → 1,351). The task anticipated DINA's n = 2000
cells as "where a wall may run long" on the strength of DINA's wide K span at n = 500; the opposite occurred,
and the shrinking family is why.

---

## 4. Cross-machine pairing, verified per replicate

The bundles carry data-only columns, so the shared seed table is checked rather than assumed. For each of
the ten cells, the Mac bundle was compared **per replicate** against the committed **pop-os** `orfs` bundle
at the same design point and n, over `seed`, `n_true` (the true-region size) and the 8 oracle columns
(the oracle refits on the TRUE region) — all rule-independent by construction. Tolerance-based at 1e-8
rather than `identical()`, because the machines differ in BLAS and R version.

**All ten cells: 10,000/10,000 column-replicate comparisons matched (1.0000)**, max relative difference
1.69e-14 at n = 500 and 4.5e-14 at n = 2000. `seed` matches under `identical()` on all 1,000 rows of every
cell; `or_H_est` does not (`identical()` FALSE) while every replicate matches at 1e-8 — the expected
signature of the same draws through a different BLAS and R version. `gate2.R`'s own `same_draws()` block
(`gate2.R:71–93`) makes the same comparison in both directions inside every cell's Gate 2, and passed there too.

---

## 5. Findings

- **All ten cells passed Gate 2 first time** — 87/87 for each GRF cell, 88/88 for each DINA cell. No halt
  record was written, no cell was re-run, and no Gate 2 was re-rendered. The supplement's contingencies for a
  halt commit and for a hostname-only Gate 2 failure were therefore never exercised.
- **The pass cost 15,195 s (4 h 13 m), not the 163,500 s (≈ 45 h) that equal per-core speed would have
  implied.** Against the pop-os projections the ten cells ran at ratios 1.71 down to 0.13; the M4 Max more
  than absorbs the 63 → 13 worker reduction on this workload. No cap was approached.
- **Memory was never a constraint**, contrary to the extrapolation from pop-os's 123–126 GB at n = 2000:
  peak 16.8–18.1 GB at 13 workers, swap 0.00M throughout, on 36 GiB.
- **DINA gets cheaper as n grows** (n-scaling 0.55 and 0.64), because its kept family K shrinks with n. This
  is the pass's main cost finding and it reverses the task's expectation.
- **DINA's declaration rate falls with n and rises with the design point**: 0.944/0.968 at n = 500 against
  0.870/0.934/0.987 at n = 2000. GRF declares at 0.999–1.000 in every cell including the benign design
  points. Both are descriptive; the synthesis task carries the inferential reading.
- **GRF's K is far tighter than DINA's** at both sizes (1,036–1,552 across all five GRF cells; 1–6,987
  across DINA's), as the task's stage-2a note anticipated.
- **The Mac's installed build was 0.3.5 from 2026-09-11** before this pass — eight days stale and a different
  version from the campaign's. The 0c install rule caught it, which is the point of installing unconditionally.
- **Procedural, for the next Mac pass:** the two constructs that would have stopped the run outright are
  `pgrep -fc` (preflight, fails immediately) and `timeout` (every render). Both are one-file shims; nothing
  in `run_or.sh` or `gate2.R` needed changing, and nothing else in the runner was Linux-only.
- **Cell order required one departure from a single runner invocation.** `run_or.sh` iterates its own
  `CELLS_ALL` order, which puts `or150_n2000` before `or100_n2000`, and `ORSG_CELLS` subsets without
  reordering, so the task's table order cannot be produced by one invocation. The driver (§6) invokes the
  unchanged runner **once per cell** in the task's order. Consequence for the record: `run_or.sh`'s
  cumulative `CUM` resets at each invocation, so the ceiling was effectively per cell, and each invocation
  emitted its own `actg175 or: campaigns <tag> complete` commit whose "cumulative render wall" is that one
  cell. The driver kept the true cumulative — **15,195 s** — and that is the figure quoted throughout.
