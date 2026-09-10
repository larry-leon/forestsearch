# REPORT — Guo & He fixed r = 1/12 on t7, production (Phase B configuration), 2026-09-09

Larry's configuration: **fixed r = 1/12, no adaptive CV**, cells `t7_beta2_00` and `t7_beta2_03`,
2000 replicates each, `B = 2000`, `orient = +1`, `v` unused (no CV). Driver
`quarto/GuoHe/guohe_sec52_adaptive_run.R` in `--fixed-r=1/12` mode, 12 workers.

**RESULT: both cells complete and clean. 2000/2000 replicates each, 0 errors, 0 selection
mismatches, no float outside `all.equal` 1e-8.** Coverage of γ_ĉ at r = 1/12 is
**0.9515** (β₂ = 0.0) and **0.9610** (β₂ = 0.3).

---

## OPEN ITEMS

Per the standing authorization: recorded, none blocks a number, none is a STOP.

1. **`pair_mm` is non-zero (170 and 210) and is not a pairing failure.** The driver's `pair_ok`
   predates the v4 N1 standard and applies plain `identical()` to float columns. Verified
   explicitly: among those replicates the differing components are **only** `naive_point` and
   `naive_lower`; `naive_dist`, `naive_bias`, `gamma_s_naive` and `cens_rate` differ on **zero**
   of them, and **every selection key remains bit-identical**. Under N1 both cells pass — see §4.
2. **Cost ran 14% over projection**: 43.9 core-h against the 38.5 core-h Gate 1b projection
   (per-replicate 39.1 core-s realized vs 34.67 measured on the 10-replicate pilot). Well inside
   the 60 core-h STOP threshold. The pilot's 20 replicates simply under-sampled the cost.
3. **Stale metadata label in fixed-r mode.** The bundles record
   `secondary_bound = "guohe_adaptive_r() own final refit bound"`. In `--fixed-r` mode the value
   actually comes from the direct `guohe_algorithm3()` call — the **computation is correct**, only
   the label is inherited from the CV path. Not corrected, to avoid touching a committed driver
   mid-run.
4. **The driver's closing `=== T2 GATE TALLY ===` prints `MISSING` for both cells.** It looks for
   the un-suffixed production filename while `--fixed-r` writes a mode-suffixed one. Cosmetic; the
   per-cell `[done]` lines carry the real tallies.
5. **The Adaptive column was not produced and is not pending.** Per the Gate 1b pilot and the
   seed-offset test (`dev/notes/NOTE_adaptive_seed_parity_2026-09-09.md`), r̂ on this design is
   substantially seed-determined; the fixed-r column is reported instead. B2 states this.

---

## 1. Provenance

```
driver : quarto/GuoHe/guohe_sec52_adaptive_run.R  --cells=t7_beta2_00,t7_beta2_03
         --reps=2000 --fixed-r=1/12 --cores=12 --force
mode   : fixed-r      fixed_r 0.0833333333333333      B 2000      orient +1
machine: Mac Studio M4 Max, 14 physical = 14 logical cores, 36 GiB; arm64, R 4.5.2, Accelerate
env    : VECLIB_MAXIMUM_THREADS=1  OMP_NUM_THREADS=1
```

No `git fetch`, `pull` or `push`. Bundles: `guohe_adaptive_t7_beta2_00_fixedr00833.rds`,
`guohe_adaptive_t7_beta2_03_fixedr00833.rds`.

## 2. Coverage of γ_ĉ at r = 1/12, with Wilson intervals

Two bounds are recorded per replicate and both are reported, because they are **not** the same
quantity:

- **Primary** — the stored `B = 2000` Algorithm-3 bound at r = 1/12, recovered from the committed
  reproduction bundle as `γ_ĉ − r2_dist`.
- **Secondary** — this run's **own** Algorithm-3 fit at r = 1/12, under its own bootstrap stream
  (`seed_data + 600000L`) rather than the reproduction's (`seed_data + 500000L`).

| cell | β₂ | coverage, primary (Wilson) | coverage, secondary (Wilson) | mean margin | bias |
|---|---|---|---|---|---|
| `t7_beta2_00` | 0.0 | **0.9515 (0.9412, 0.9601)** | 0.9515 (0.9412, 0.9601) | 0.2922 / 0.2925 | +0.01376 / +0.01379 |
| `t7_beta2_03` | 0.3 | **0.9610 (0.9516, 0.9686)** | 0.9620 (0.9527, 0.9695) | 0.3121 / 0.3121 | −0.00734 / −0.00726 |

Margin and bias are given primary / secondary. MCSE at 2000 replicates near 0.95 is ≈ 0.0049.

Both cells sit at or slightly above nominal, consistent with the published Table 7 columns running
conservative (0.947–0.973): adjacent nested candidates differ by a single subject, so the effective
number of independent comparisons is far below the ~151 family size.

## 3. Agreement with the committed fixed-r column

Same replicates, same `m`, compared against `r2_*` in `guohe_repro_t7_beta2_0{0,3}.rds`
(`r2` is r = 1/12 in `GH52_R_GRID`).

| cell | stored `r2_cover` | primary flag vs stored | secondary flag vs stored |
|---|---|---|---|
| `t7_beta2_00` | 0.9515 (0.9412, 0.9601) | **2000/2000, `identical()` TRUE** | **1980/2000 (99.00%)** |
| `t7_beta2_03` | 0.9610 (0.9516, 0.9686) | **2000/2000, `identical()` TRUE** | **1994/2000 (99.70%)** |

**They agree — but read the two rows differently.**

- **The primary agreement is tautological and is reported as a integrity check, not as evidence.**
  The primary bound *is* the stored bound, recovered by arithmetic from the stored `r2_dist`.
  `identical()` holding at 2000/2000 (and `primary dist identical to stored r2_dist` TRUE in both
  cells) confirms the lookup and the `γ_ĉ − dist` reconstruction are exact. It says nothing about
  reproducibility of the method.
- **The secondary agreement is the real comparison.** This run refits Algorithm 3 from scratch at
  r = 1/12 under an independent bootstrap stream at the same `B = 2000`, and reaches the **same
  coverage decision on 99.0% and 99.7% of replicates**. Bound locations track closely: mean
  difference in margin +0.00028 and −0.00009, with worst-case single-replicate differences of
  0.058 and 0.051 log-HR.

The residual 1.0% / 0.3% disagreement is Monte Carlo, not method: at `B = 2000` two independent
bootstrap streams place the bound a little differently, and a replicate whose bound sits within
that jitter of γ_ĉ can flip its coverage indicator. **Conclusion: the committed fixed-r column at
r = 1/12 reproduces.**

## 4. Selection tally and N1 deviation summary

| check | `t7_beta2_00` | `t7_beta2_03` |
|---|---|---|
| replicates | 2000/2000 | 2000/2000 |
| errored | 0 | 0 |
| `ad_err` | 0 | 0 |
| **selection mismatch** (`ad_sel_ok == 0`) | **0** | **0** |
| `r_hat_offgrid` | 0 | 0 |
| **N1 selection keys, plain `identical()`** — `c_hat_naive`, `c_hat_gh`, `n_sel`, `naive_cover`, `seed_data` | **ALL PASS** | **ALL PASS** |
| **N1 floats, `all.equal` 1e-8** — `naive_point`, `naive_lower`, `gamma_s`, `cens_rate` | **ALL PASS** | **ALL PASS** |
| worst absolute float deviation | 2.60e-15 (`naive_point`, m = 298) | 3.55e-15 (`naive_point`, m = 724) |
| float values bit-identical | 7768/8000 (97.1%) | 7707/8000 (96.3%) |
| `naive_point` bit-identical | 1903/2000 | 1906/2000 |
| driver `pair_ok == 0` (pre-N1 flag) | 170 | 210 |

**Zero selection mismatches in 4000 replicates.** The 170 and 210 `pair_ok` failures decompose as:

| cell | `naive_point` differs | `naive_lower` differs | `naive_dist` | `naive_bias` | `gamma_s_naive` | `cens_rate` | selection keys still identical |
|---|---|---|---|---|---|---|---|
| `t7_beta2_00` | 97 | 135 | 0 | 0 | 0 | 0 | **TRUE** |
| `t7_beta2_03` | 94 | 199 | 0 | 0 | 0 | 0 | **TRUE** |

Confined entirely to the two naive float columns, at ≤ 3.55e-15 — the same arm64-vs-x86_64
signature measured throughout T1.

## 5. Cost against projection

| quantity | Gate 1b projection | realized | ratio |
|---|---|---|---|
| per-replicate (serial) | 34.67 core-s | **39.1 core-s** | ×1.13 |
| total, 2 cells × 2000 | 38.5 core-h | **43.9 core-h** | ×1.14 |
| wall at 12 workers | 3.21 h | **3.66 h** | ×1.14 |
| STOP threshold | — | 60 core-h | **not reached** |

Per cell: 110.21 min (`t7_beta2_00`, 39.26 core-s/rep) and 109.22 min (`t7_beta2_03`,
38.89 core-s/rep). The overrun is a pilot sampling artefact — 20 replicates estimated the mean
per-replicate cost 13% low — not a change in configuration.

## 6. Reading

- **The fixed-r bound at r = 1/12 is reproducible and behaves as published.** Coverage 0.9515 and
  0.9610 against a nominal 0.95, on a design whose published columns are themselves conservative.
- **An independent refit at the same `B` reaches the same coverage decision on 99.0–99.7% of
  replicates**, with mean bound locations within 0.0003 log-HR. The committed reproduction stands.
- **Nothing here required the adaptive path.** The Gate 1b pilot measured that path at ×10.28 the
  cost, and the seed-offset test showed its selection is substantially seed-determined on this
  design. Reporting r = 1/12 costs 1/10 as much and is a stated choice rather than an unstable one.
