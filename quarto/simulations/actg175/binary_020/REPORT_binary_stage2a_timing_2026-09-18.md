# REPORT — Binary campaign, stage 2a: the two timing cells, and the projection for the remaining 10

Task: `dev/tasks/TASK_binary_launch_v2_2026-09-18.md` Step 4. Machine `pop-os` (128 physical cores,
251 GB RAM), R 4.6.1, **forestsearch 0.3.5.9000, built 2026-09-19 03:54:02 UTC**. Branch
`feature/glm-extension`, HEAD at launch `29689584`. 1,000 replicates, 63 workers, MR only.
**This task issued no push, and nothing was launched beyond these two cells** — the stage-2
go/no-go is Larry's. An external push to `origin` occurred mid-run and is recorded in
`REPORT_binary_stage1_fs_2026-09-18.md` §7; neither of these two cells nor this report is on
`origin`.

---

## 1. The two cells, measured

Both completed far inside the 3 h cap, and both passed their gate outright.

| cell | Gate 2 | cell wall | **vs the 3 h cap** | batch render | combine | peak MB | commit |
|---|---|---|---|---|---|---|---|
| `orgrf_or150_n500`  | **87/87** | **937 s** (15 min 37 s) | **8.7%** | 916 s | 20 s | 61,489 | `5ef6b748` |
| `ordina_or150_n500` | **88/88** | **1,359 s** (22 min 39 s) | **12.6%** | 1,338 s | 21 s | 83,300 | `adfd783f` |

Cumulative render wall 2,295 s against `ORSG_CEILING=22000`. **Neither cell was aborted; the cap
was never approached**, so there is nothing to record under the task's abort clause.

Gate 2 runs more checks here than on `orfs` (87 and 88 against 74) because both identifiers are
checked for **same draws against `orfs` in both directions** — the per-replicate seed, the
true-region size and the eight oracle columns identical within 1e-8 — plus each identifier's own
accounting. All passed.

### 1.1 The three identifiers at the same cell, side by side

Descriptive only. GRF's and DINA's candidate families are generated from fitted surfaces, so every
coverage figure of those two campaigns is coverage of the estimand **conditional on the proposed
family**; FS's family is the prespecified cut grid. These rows are read side by side and never
ranked (`status_curated.md` §5).

| | `orfs` | `orgrf` | `ordina` |
|---|---|---|---|
| cell wall | 1,237 s | 937 s | 1,359 s |
| **cost relative to `orfs`** | 1.000 | **0.758** | **1.099** |
| declaration | 0.9310 | 1.0000 | 0.9790 |
| `fit_mr_secs` mean / median / p90 / max | 57.0 / 60.8 / 69.6 / 83.0 | 40.3 / 41.2 / 44.6 / 49.1 | 52.1 / **37.4** / **119.4** / **211.6** |
| `id_secs` mean / median / max | 9.95 / 9.72 / 19.23 | 5.49 / 5.62 / 7.14 | 5.27 / 3.08 / 28.50 |
| MR family K: min / med / p90 / max | 1850 / 2226 / 2260 / 2293 | 1036 / 1066 / 1146 / 1298 | **1 / 873 / 3881 / 6075** |
| peak MB | 79,041 | 61,489 | 83,300 |
| NA-oracle, non-estimable, MR failures | 0, 0, 0 | 0, 0, 0 | 0, 0, 0 |

Two things in that table drive everything below.

**The MR gate dominates and it scales with the kept family K.** `fit_mr_secs` is 4–8× `id_secs` for
every identifier, so identification cost is nearly irrelevant to the total; what matters is how
large a family each identifier hands the gate. GRF keeps the smallest (median 1,066) and is
correspondingly the cheapest. That is why `orgrf` runs *faster* than `orfs` despite declaring on
every replicate.

**DINA's cost is not summarised by its mean.** Its kept family runs from **1 to 6,075** and its
`fit_mr_secs` has a long right tail — median 37.4 s but p90 119.4 s and max 211.6 s, against a
median above its mean for both other identifiers. The cell-level wall still averages out over 1,000
replicates, but the variance is a real feature of the projection below, not a rounding concern.

---

## 2. Projecting the remaining 10 cells

**Method.** The identifier factor is measured at one cell only (`or150_n500`): 0.758 for `orgrf`,
1.099 for `ordina`, against the `orfs` cell at the same design point and size. Each remaining cell
is projected as *the measured `orfs` cell at that (design point, n)* × *that identifier's factor* —
so the design-point and sample-size structure comes from six measured `orfs` cells rather than from
a single scaling law, and only the identifier ratio is extrapolated.

| cell | `orfs` measured | × factor | **projected** |
|---|---|---|---|
| `orgrf_or075_n500`  | 1,105 s | 0.758 | 838 s |
| `orgrf_or100_n500`  | 1,135 s | 0.758 | 860 s |
| `orgrf_or075_n2000` | 5,054 s | 0.758 | 3,831 s |
| `orgrf_or150_n2000` | 5,555 s | 0.758 | 4,211 s |
| `orgrf_or100_n2000` | 5,321 s | 0.758 | 4,033 s |
| **`orgrf` subtotal (5 cells)** | | | **13,773 s ≈ 3 h 50 min** |
| `ordina_or075_n500`  | 1,105 s | 1.099 | 1,214 s |
| `ordina_or100_n500`  | 1,135 s | 1.099 | 1,247 s |
| `ordina_or075_n2000` | 5,054 s | 1.099 | 5,554 s |
| `ordina_or150_n2000` | 5,555 s | 1.099 | 6,105 s |
| `ordina_or100_n2000` | 5,321 s | 1.099 | 5,848 s |
| **`ordina` subtotal (5 cells)** | | | **19,968 s ≈ 5 h 33 min** |
| **remaining 10 cells** | | | **33,741 s ≈ 9 h 22 min** |

Adding the two cells already run, the full 12-cell GRF/DINA half of the campaign is **36,036 s
≈ 10 h 01 min** of render wall, of which 2,295 s is done.

**A ceiling for it** should sit near **45,000 s** for the remaining 10 (≈ 33% headroom), on the same
reasoning that made 30,000 s right for stage 1's 18,351 s estimate.

**Peak memory is not a constraint.** The identifier memory factors at n = 500 are 0.78 (`orgrf`)
and 1.05 (`ordina`) against `orfs`. Applied to the measured `orfs` n = 2000 peak of 126 GB, the
worst projected cell is `ordina_or150_n2000` at ≈ 133 GB, against 251 GB of RAM — the same ~50%
headroom stage 1 ran at, at 63 workers.

### 2.1 Four caveats, all honest

1. **The identifier factor rests on one cell.** It is measured at `or150_n500` only. It is applied
   to two other design points and, more consequentially, to n = 2000 — which is 15 of the 10 cells'
   19,968 s for `ordina`. The mechanism is at least the right one (the MR gate dominates, and it
   scales with the kept family K, which the identifier fixes), but it is one measurement.
2. **DINA's family may not scale like FS's.** From n = 500 to n = 2000, FS's median K rose from
   2,226 to 2,985 (×1.34). DINA's `admitted_n` at n = 500 already spans 1 to 5,596 with a median of
   543 and a p90 of 2,627; whether its median K grows by ×1.34 or faster at n = 2000 is
   **unmeasured**, and the projection assumes the former by construction. If DINA's family grows
   faster than FS's, its three n = 2000 cells — 17,507 s of the 33,741 s — are underestimated. This
   is the single largest uncertainty in the number, and it is one cheap `ordina_or075_n2000` run
   away from being resolved.
3. **The OR 1.5 surcharge is inside the base, not on top of it.** The `orfs` cells this projection
   scales are the measured ones, so the ~10% that OR 1.5 cost over the linear estimate in stage 1
   is already carried. No further adjustment is needed.
4. **Stage 1's own overrun was +5.8%.** Applying the same slippage to 33,741 s gives ≈ 35,700 s
   ≈ 9 h 55 min as a realistic upper-middle figure rather than a floor.

---

## 3. Where this stops

Per the task, **stage 2 is not launched**. The two timing cells are run, gated, committed and
reported; the remaining 10 GRF/DINA cells wait on Larry's go/no-go against the numbers above.

Recommended launch line if it is a go — the runner sequences both campaigns over the five remaining
cells each and skips the two already committed:

```
cd quarto/simulations/actg175/binary_020
ORSG_WORKERS=63 ORSG_TIMEOUT=10800 ORSG_CEILING=45000 ORSG_NSIMS=1000 \
  ORSG_CAMPAIGNS="orgrf ordina" \
  setsid nohup bash scripts_or/run_or.sh > logs_or/runner_stage2.log 2>&1 &
```

The cheapest way to retire caveat 2 before committing to all ten is to run
`ORSG_CELLS="or075_n2000"` first (projected 5,554 s ≈ 1 h 33 min) and check the measured wall
against it.
