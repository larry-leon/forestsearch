# REPORT — Gate 1b pilot (T2 Adaptive column), Mac Studio, 2026-09-09

**PILOT ONLY. No production was launched, and none will be from this record — Larry chooses the
configuration.** Two cells, `t7_beta2_00` and `t7_beta2_03`, 10 replicates each, `B = 2000`,
`orient = +1`, `v = 5`, per v1 §5 as amended by v2 A5.

Driver: `quarto/GuoHe/guohe_sec52_adaptive_run.R`. No `git fetch`, `pull` or `push`.

---

## 1. Machine and usable cores

| quantity | value |
|---|---|
| `parallel::detectCores(logical = FALSE)` | **14** |
| `parallel::detectCores(logical = TRUE)` | **14** |
| `hw.physicalcpu` / `hw.logicalcpu` | 14 / 14 |
| `hw.memsize` | 38,654,705,664 (36 GiB) |
| **used here** | **12 workers** |

This Mac Studio (M4 Max) reports **no hyperthreading** — logical equals physical at 14 — so 14 is
the hard ceiling, not 28. T1 ran at 10; **12 is usable and was used here**, leaving 2 cores for the
OS. 14 would be available at the cost of running the machine fully saturated. Memory is not a
constraint: peak resident set stayed far below the 36 GiB ceiling (each worker holds a 400-row
design and a `400 × 2000` bootstrap matrix).

`devtools::install()` was run first, so the workers see the installed package.

## 2. What was measured, and why both

Larry has fixed **r = 1/12**; Guo & He's Table 7 varies little across 1/12, 1/21 and 1/30.

**The driver's Algorithm-2 path selects r by cross-validation, so a single r leaves nothing to
select.** Running the CV machinery over a one-element grid would fit 5 folds × 1 candidate and then
"choose" the only candidate — a degenerate measurement that reports CV overhead without any CV.
Both configurations were therefore measured separately:

- **(a) FIXED r = 1/12** — `guohe_algorithm3()` at one r, **no CV at all**. The CV path is skipped
  entirely rather than run degenerately.
- **(b) ADAPTIVE, `r_grid = c(1/3, 1/12)`** — the genuine CV path with 1/12 in the candidate set,
  so the adaptive marginal cost is **measured, not inferred**. 1/3 is the natural partner: it is
  the far end of the published grid, so the CV has a real choice to make.

**(b) required only driver-level argument changes**, so both are reported. Two optional flags were
added, `--fixed-r=` and `--r-grid=`, **both defaulting to the committed behaviour** (the full
published grid, CV path, `B = 2000`). Omit them and the driver runs exactly as committed; no
default was changed permanently. Output filenames carry a mode suffix so the pilot bundles cannot
be mistaken for production ones.

## 3. Pairing proof — all 20 replicates, under the v4 N1 standard

Per-replicate cost never approached the 20-minute abort threshold ((a) 0.6 min, (b) 6.0 min), so
all 20 replicates ran in both configurations.

| check | `t7_beta2_00` | `t7_beta2_03` |
|---|---|---|
| **Selection keys, plain `identical()`** — `c_hat_naive`, `c_hat_gh`, `n_sel`, `naive_cover`, `seed_data` | **ALL PASS** | **ALL PASS** |
| **Floats, `all.equal` 1e-8** — `naive_point`, `naive_lower`, `gamma_s` | **ALL PASS**, worst \|diff\| 0 | **ALL PASS**, worst \|diff\| 1.11e-16 |
| `ad_sel_ok` (adaptive refit lands on the stored selection) | 10/10 | 10/10 |
| `sel_mm` | 0 | 0 |
| `r_hat_offgrid` | 0 | 0 |

**The driver's own `pair_ok` counter reported `pair_mm 2` on `t7_beta2_03` in both runs. That is
not a pairing failure.** `pair_ok` predates the v4 N1 standard and applies plain `identical()` to
float columns. Broken down:

| replicate | component | `identical()` | \|diff\| | `all.equal` 1e-8 |
|---|---|---|---|---|
| m = 3 | `naive_point` | FALSE | 1.11e-16 | TRUE |
| m = 3 | `naive_lower` | FALSE | 1.11e-16 | TRUE |
| m = 4 | `naive_lower` | FALSE | 5.55e-17 | TRUE |

Every selection key and every truth lookup at those replicates is bit-identical; only the two
naive float columns differ, at the 1e-16 scale, exactly as in T1. Under N1 both replicates pass.

## 4. Measured cost

Serial core-seconds per replicate, pooled over the two cells (n = 20 each):

| configuration | mean | median | max | per-cell wall (10 reps, 12 cores) |
|---|---|---|---|---|
| **(a) fixed r = 1/12** | **34.67 core-s** | 34.3 s | 38.6 s | 0.63–0.65 min |
| **(b) adaptive CV, `r_grid = c(1/3, 1/12)`** | **356.54 core-s** | 353 s | 397 s | 6.55–6.63 min |

**Adaptive marginal cost: ×10.28 over fixed-r** (356.5 vs 34.7 core-s). That is close to the
structural expectation — the CV path runs `v × |r_grid| = 5 × 2 = 10` inner Algorithm-3 fits plus
one final refit, all at `B = 2000`, against the single fit of (a).

**A note on the cost the reproduction ledger would have implied.** The committed t7 reproduction
bundles record ~688 s per replicate for four fixed-r fits (`guohe_repro_t7_beta2_00.rds`,
median 688.2 s; `_03`, median 698.0 s), implying ~172 s per fit. This Mac does the same fit in
**34.7 s — roughly 5× faster**. The projections below are this machine's measurements, not the
ledger's.

### Projections, 12 workers

| scenario | core-h | wall @ 12 cores |
|---|---|---|
| **(a) fixed r = 1/12, 2 cells × 2000** | **38.5** | **3.21 h** |
| **(b) adaptive CV, 2 cells × 2000** | **396.2** | **33.01 h** |
| (a) fixed r = 1/12, 2 cells × 500 | 9.6 | 0.80 h |
| (b) adaptive CV, 2 cells × 500 | 99.0 | 8.25 h |

### `B = 200` sensitivity — **labelled, not measured**

Scaled from the measurements above assuming bootstrap cost is linear in `B`. That assumption is
untested here; it ignores the fixed per-fit overhead (the ~151-candidate Cox loop), which does not
shrink with `B`, so these are **optimistic lower bounds** and the true figures will be somewhat
higher.

| scenario | core-h | wall @ 12 cores |
|---|---|---|
| (a) fixed r = 1/12, `B = 200`, 2 cells × 2000 | 3.9 | 0.32 h |
| (b) adaptive CV, `B = 200`, 2 cells × 2000 | 39.6 | 3.30 h |
| (a) fixed r = 1/12, `B = 200`, 2 cells × 500 | 1.0 | 0.08 h |
| (b) adaptive CV, `B = 200`, 2 cells × 500 | 9.9 | 0.83 h |

## 5. (b) — r̂ and the per-candidate CV objective

Ten replicates per cell, `r_grid = c(1/3, 1/12)`, `v = 5`, `B = 2000`.

### `t7_beta2_00`

- r̂ per replicate: 0.0833, 0.3333, 0.3333, 0.0833, 0.3333, 0.0833, 0.3333, 0.0833, 0.0833, 0.0833
- **r̂ distribution: 1/12 six times, 1/3 four times**
- CV objective at r = 1/3: mean −0.05987, range [−0.19818, +0.68328]
- CV objective at r = 1/12: mean −0.06460, range [−0.20704, +0.65642]
- objective difference (1/12 − 1/3), negative favours 1/12: **mean −0.00473; 1/12 wins 6/10**
- coverage: primary 10/10, secondary 9/10; primary-flag agreement with the stored column 10/10

### `t7_beta2_03`

- r̂ per replicate: 0.3333, 0.0833, 0.3333, 0.0833, 0.3333, 0.0833, 0.3333, 0.0833, 0.3333, 0.0833
- **r̂ distribution: 1/12 five times, 1/3 five times**
- CV objective at r = 1/3: mean −0.08315, range [−0.17842, +0.10462]
- CV objective at r = 1/12: mean −0.08300, range [−0.18037, +0.10807]
- objective difference (1/12 − 1/3): **mean +0.00015; 1/12 wins 5/10**
- coverage: primary 9/10, secondary 9/10; primary-flag agreement with the stored column 10/10

### Reading — the finding of this pilot

**The CV objective barely separates the two candidates, and r̂ is close to a coin flip.** The mean
objective gap is −0.0047 at β₂ = 0 and **+0.00015** at β₂ = 0.3 — the latter is three orders
smaller than the objective's own spread across replicates (range ≈ 0.28) — and r̂ splits 6/4 and
5/5. On 20 replicates that is indistinguishable from selecting r at random between 1/3 and 1/12.

This is the mechanism `guohe_reproduction_RUN.md:131-136` warns about: `guohe_adaptive_r()` draws
**independently across r**, with no common random numbers on the grid, so the between-candidate
Monte Carlo noise enters Var(r̂) directly — and the Adaptive column is precisely a measurement of
r̂. It is also consistent with Guo & He's own Table 6 caution on the adaptive procedure.

The practical consequence for the decision: **paying ×10.28 buys a selection that, on this design
and this grid, is close to arbitrary.** Coverage is not the discriminator either — (a) and (b)
land within one replicate of each other on both cells (primary 10/10 and 9/10 in both).

### (a) coverage, fixed r = 1/12

| cell | primary | secondary | `ad_sel_ok` |
|---|---|---|---|
| `t7_beta2_00` | 10/10 | 9/10 | 10/10 |
| `t7_beta2_03` | 9/10 | 9/10 | 10/10 |

Ten replicates per cell carry no useful coverage precision (MCSE ≈ 0.07 at 0.95); these are
reported as pilot sanity checks, not as coverage estimates.

## 6. Caveats carried forward

1. **Guo & He's own Table 6 caution** on the adaptive procedure stands.
2. **Independent draws across r.** `guohe_adaptive_r()` uses no common random numbers across the
   grid, inflating Var(r̂) — §5 shows this is not hypothetical on this design.
3. **The `B = 200` figures are scaled, not measured**, and are optimistic for the reason given.
4. **Two cells, 10 replicates.** Everything here is a cost and mechanism measurement. No coverage
   claim is made.

## 7. Status

- Gate 1b pilot complete for two cells at r = 1/12, both configurations.
- Pairing proof passes on all 20 replicates under N1 in both configurations.
- **No production launched. No default changed permanently.** The two new driver flags default to
  the committed behaviour.
- **Waiting on Larry to choose the configuration.**

Pilot bundles written beside this record:
`guohe_adaptive_t7_beta2_00_fixedr00833.rds`, `guohe_adaptive_t7_beta2_03_fixedr00833.rds`,
`guohe_adaptive_t7_beta2_00_grid2.rds`, `guohe_adaptive_t7_beta2_03_grid2.rds`. They carry mode
suffixes and sit outside the production `guohe_adaptive_t7_beta2_NN.rds` naming, so they cannot be
picked up by the T3 qmd's Adaptive-column glob, which still reads "pending Phase B".

One cosmetic defect noted, not fixed: the driver's closing `=== T2 GATE TALLY ===` block looks for
the un-suffixed production filename and therefore prints `MISSING` for both cells after a
mode-suffixed run. The per-cell `[done]` lines above it carry the real tallies.

---

# APPENDED — OPEN ITEM: r̂ tracks the parity of the replicate index

Found while checking the r̂ sequences before closing the record. **Recorded as an observation with
a plausible mechanism; the mechanism is not proven, and no further compute was spent on it.**

## The observation

`t7_beta2_03`, r̂ against replicate index m:

| m | r̂ | obj(1/3) | obj(1/12) | diff (1/12 − 1/3) |
|---|---|---|---|---|
| 1 | 0.3333 | +0.052039 | +0.054375 | +2.337e-03 |
| 2 | **0.0833** | −0.178419 | −0.180367 | −1.948e-03 |
| 3 | 0.3333 | −0.136588 | −0.129679 | +6.909e-03 |
| 4 | **0.0833** | −0.032957 | −0.041704 | −8.747e-03 |
| 5 | 0.3333 | +0.104624 | +0.108071 | +3.447e-03 |
| 6 | **0.0833** | −0.159534 | −0.161324 | −1.790e-03 |
| 7 | 0.3333 | −0.128499 | −0.124991 | +3.509e-03 |
| 8 | **0.0833** | −0.130256 | −0.130703 | −4.473e-04 |
| 9 | 0.3333 | −0.053566 | −0.053285 | +2.813e-04 |
| 10 | **0.0833** | −0.168341 | −0.170355 | −2.014e-03 |

**r̂ = 1/12 on 5/5 even m and 0/5 odd m — perfect alternation.** Under an even coin that split has
probability ≈ 0.2% (two-sided). `t7_beta2_00` leans the same way without being perfect: 1/12 on
4/5 even and 2/5 odd.

The objective differences driving these flips are tiny — median |diff| 2.18e-03, max 8.75e-03,
against an objective whose spread across replicates is ≈ 0.28. Every selection is a near-tie, and
**the sign of the near-tie tracks the parity of m.**

## Plausible mechanism, not verified

The adaptive call is seeded `seed_ad = base + m + 600000L`, so consecutive replicates differ by 1
in the seed. `guohe_adaptive_r()` draws its v-fold assignment from that seed. Consecutive
Mersenne-Twister seeds are a well-known source of correlated early draws, so a parity-linked fold
assignment feeding a near-tied CV comparison is a credible explanation. **This was not tested** —
doing so means more runs, and this is a pilot.

## Why it matters for the decision

§5 read r̂ as "close to a coin flip". This is worse than a coin flip: if the association holds, r̂
is **partly a deterministic function of the replicate index rather than of the data.** An Adaptive
column at 2000 replicates would then return ≈ 50/50 r̂ by construction, and its coverage would be a
blend of the two fixed-r columns in a ratio set by the seed grid — not a measurement of adaptive
selection.

That strengthens, rather than changes, §5's conclusion: the ×10.28 buys a selection that on this
design and grid is close to arbitrary.

## Cheapest test, if Larry wants it

Re-run the same 10 replicates of `t7_beta2_03` with a different seed offset (e.g.
`GHA_SEED_OFFSET` 600000 → 600001, or a stride of 7 instead of 1). If the parity alignment moves
with the offset, it is the seed grid; if r̂ is unchanged, it is the data. Cost: ~6 min wall at 12
workers, one cell, no production implications. **Not run — awaiting instruction.**

---

# APPENDED — seed-offset test: the parity IS a seed-grid artefact (CONFIRMED)

The discriminating test proposed above was authorized and run. **Cell `t7_beta2_03`, the same 10
replicates, everything identical except the adaptive seed offset: `+600000L` → `+700000L`.**

The committed driver was **not modified**. The test replicates its per-replicate logic in a
scratch script — same data regeneration from `base + m`, same `guohe_adaptive_r()` call at
`orient = +1`, `r_grid = c(1/3, 1/12)`, `v = 5`, `B = 2000`, `min_events = 5`, `refit = TRUE` —
changing only the offset. Output: `guohe_adaptive_t7_beta2_03_grid2_seed700k.rds`.

## Result

| m | seed (+600000) | r̂ @ 600k | seed (+700000) | r̂ @ 700k | same? |
|---|---|---|---|---|---|
| 1 | 93802767 | 0.3333 | 93902767 | **0.0833** | no |
| 2 | 93802768 | 0.0833 | 93902768 | **0.3333** | no |
| 3 | 93802769 | 0.3333 | 93902769 | **0.0833** | no |
| 4 | 93802770 | 0.0833 | 93902770 | 0.0833 | YES |
| 5 | 93802771 | 0.3333 | 93902771 | 0.3333 | YES |
| 6 | 93802772 | 0.0833 | 93902772 | 0.0833 | YES |
| 7 | 93802773 | 0.3333 | 93902773 | **0.0833** | no |
| 8 | 93802774 | 0.0833 | 93902774 | 0.0833 | YES |
| 9 | 93802775 | 0.3333 | 93902775 | **0.0833** | no |
| 10 | 93802776 | 0.0833 | 93902776 | 0.0833 | YES |

| | r̂ = 1/12 on even m | r̂ = 1/12 on odd m |
|---|---|---|
| offset **+600000L** | **5/5** | **0/5** |
| offset **+700000L** | 4/5 | 4/5 |

- **The parity alignment moved with the offset.** Perfect at +600000L, gone at +700000L.
- **r̂ changed on 5 of 10 replicates from a seed change alone** — the data, the candidate family
  and the selection were untouched.
- Selection keys reproduced 10/10 (`c_hat_naive`, `naive_cover`), so the data regeneration is
  unaffected; only the CV's internal draws moved.
- Cost unchanged: mean 367.9 s per replicate (vs 357.3 s at the committed offset).

## Verdict

**Seed-grid artefact, confirmed.** r̂ on this design and grid is substantially determined by the
adaptive seed rather than by the data. The parity pattern was the visible symptom; the
5-of-10 flip is the direct measurement.

This does not change the Gate 1b recommendation — it sharpens it. §5 read r̂ as "close to a coin
flip"; the test shows the coin is weighted by the seed. Recorded in full at
`dev/notes/NOTE_adaptive_seed_parity_2026-09-09.md`. **No repair made; the seed derivation in the
committed driver is untouched; disposition is Larry's.**

## Bundle naming

The four pilot bundles were renamed to carry an explicit `_pilot` suffix before the fixed-r
production run, so production output at the canonical `..._fixedr00833.rds` names cannot
overwrite the pilot evidence:

- `guohe_adaptive_t7_beta2_00_fixedr00833_pilot.rds`
- `guohe_adaptive_t7_beta2_03_fixedr00833_pilot.rds`
- `guohe_adaptive_t7_beta2_00_grid2_pilot.rds`
- `guohe_adaptive_t7_beta2_03_grid2_pilot.rds`
- `guohe_adaptive_t7_beta2_03_grid2_seed700k.rds` (this test)

All five sit outside the T3 Adaptive-column glob.
