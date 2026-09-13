# REVIEW — `grfmr` completion: the two deferred harm cells and the six HR 1.00 cells

- **Date:** 2026-09-12
- **Reviews:** `REPORT_grfmr_completion_2026-09-12.md`, commit range `6385b42a..9d4f4405`, not pushed. No `R/` or template change.
- **State:** `grfmr` is complete at **18 of 18 cells**, matching the DINA grid cell for cell. Eight cells at 2,000 replicates, nothing deferred or dropped, 6.425 h against 6.744 h projected (0.953).
- Both model-generated identifiers now have full grids on the same DGM draws, with FS beside them at every cell.

---

## 1. Verification — accepted

- Stage 1 smoke 30/30. Gate 3 7 passes on every one of 16 batches. Gate 2 37 checks per cell, all passing; corrected identity ≤ 3.33e-16.
- The designated comparator resolved for every cell — `cert20` at 31%, `tier2` or `p12ext` at 12.4%, never `map1` or `s7`. The check CC added for this was worth adding: `e1stud` has no HR 1.00 bundle, and the substitution it guards against would have been silent.
- Amendment 3 holds on all eight: `n_true` identical on every row, truth within 1e-8.
- The 10 h watchdog never fired. Five of the six null cells were costed from the more expensive HR 1.50 probes, which is why they came in 6.6% under.
- Closeout clean: 1435/1435 tracked, largest 7.86 MB, `current_status.md` regenerated last with the pin check passing either way.

## 2. The null-cell result

### 2.1 Selection rate, with Wilson intervals at 2,000

| Cell | GRF | FS (same draws) | DINA (same draws) |
|---|---|---|---|
| 12.4% n 500 | 0.9840 [0.978, 0.989] | 0.681 | 0.714 |
| 12.4% n 1000 | 0.9575 [0.948, 0.966] | 0.660 | 0.524 |
| 12.4% n 1500 | 0.8675 [0.852, 0.882] | 0.624 | 0.344 |
| 31% n 500 | 0.9970 [0.994, 0.999] | 0.921 | 0.930 |
| 31% n 1000 | 0.9980 [0.995, 0.999] | 0.955 | 0.926 |
| 31% n 1500 | 0.9935 [0.989, 0.996] | 0.959 | 0.888 |

- **At 31% GRF's selection rate does not separate a null region from a harm region at all** — it is ≥ 0.993 at every null cell and 1.0000 at every harm cell. Reported alone, it would say nothing.
- At 12.4% it falls with n but far less steeply than the others: 0.984 → 0.868, against FS 0.681 → 0.624 and DINA 0.714 → 0.344.
- **The three identifiers sit at three points on the same axis.** DINA declines to return a region as n grows; FS declines moderately and stably; GRF nearly always returns one. That is an identifier-level property, and it is the cleanest three-way contrast either campaign produced.

### 2.2 Bound location — the diagnostic that does separate

| Cell | GRF share ≥ 1.00 [Wilson] | share ≥ 1.25 |
|---|---|---|
| 12.4% n 500 | 0.0061 [0.0035, 0.0106] | ≤ 0.0065 |
| 12.4% n 1000 | 0.0042 [0.0021, 0.0083] | across the six |
| 12.4% n 1500 | 0.0046 [0.0023, 0.0091] | null cells |
| 31% n 500 | 0.0196 [0.0144, 0.0267] | |
| 31% n 1000 | 0.0110 [0.0073, 0.0166] | |
| 31% n 1500 | 0.0091 [0.0058, 0.0143] | |

- At matched prevalence and n, the harm cells' ≥ 1.00 share is **3.8 to 90 times** the null share. Location separates harm from null where the selection rate does not.
- At 12.4% the share is flat within its intervals; at 31% it falls with n.
- **This is the reporting standard demonstrated on the identifier that most needed it.** On GRF the selection rate would read as a 87–100% "error rate" at cells where the bound reaches HR 1.00 on 0.4–2.0% of replicates.

### 2.3 Criterion-matched comparison at 31%

- Share ≥ 1.00: GRF 0.0196 / 0.0110 / 0.0091 against FS 0.0060 / 0.0105 / 0.0094 and DINA 0.0301 / 0.0157 / 0.0096.
- At n 500 GRF sits above FS with non-overlapping intervals; at n 1000 and 1500 the three are indistinguishable. DINA is highest at n 500 and converges by n 1500.
- The confound travels: identifier, family construction, detection set. At 12.4% the criterion differs too (`maxeffCons` ε 0.10 against `effMaxSG` ε 0.20).

### 2.4 One n-trend caveat, and only one

- The 12.4% null cells are the only GRF cells whose rate moves with n by more than 0.01, so a trend read there conditions on a shrinking detected set. Every other GRF n-trend in the campaign is free of that qualification — unlike DINA, where it attaches to all of them.

## 3. The empty admitted set fires for the first time

- 405 non-detections across the six null cells, all `NO-DETECTION`, none an error.
- 392 carry `admitted_n` NA; **13 carry `admitted_n` = 0**.
- Across the ten harm cells and the 180 probe replicates, `admitted_n` was never 0. The path exists and is reachable — at a null, at low prevalence, where nothing clears the floor.
- **This bears on the §8 frontier-filter question.** The earlier finding was that the frontier band cannot empty at `dmin.grf = 0.0` and the protected branch is unreachable. That remains true of the *band*; what these 13 replicates show is the *floor* emptying, which is the distinct mechanism the probe analysis identified. Consistent, not contradictory — and worth stating precisely wherever the asymmetry is discussed, since "unreachable" without the qualifier would now be wrong.
- 13 of 12,000 replicates is too few to characterise; recorded, not analysed.

## 4. `admitted_n` at the null cells

| Cell | median (q25–q75) | Enumerated pool |
|---|---|---|
| 12.4% n 500 / 1000 / 1500 | 61.5 (32–113) / 38 (17–70) / 23 (11–42) | 776–830, outcome-independent |
| 31% n 500 / 1000 / 1500 | 166 (95.5–281) / 141 (82–225) / 115 (66–176) | same |
| the two late harm cells | 487.5 and 493.5 | same |

- The qualified count falls with n at both prevalences and is roughly 3× larger at 31% than 12.4%, against an enumerated pool that does not move with either. The stratifier substitution continues to be the right one.
- At the null cells the admitted share runs 3–21% of the pool, against roughly 50% at the 31% harm cells.

## 5. Dispositions

- Report accepted. No re-run. Nothing proposed for CC beyond section 6.
- The `gate2G.R` extension and its comparator-resolution check: accepted, and regression-tested against the ten committed harm cells before use.
- The two stale citations CC corrected in the null section: accepted.

## 6. Outstanding

- **Leftover DINA wording in `summary_grfmr.qmd`**, flagged by CC and not fixed: line 26 "the same one GRF carries", line 41 and the strata section "GRF's family-size stratifier", and the strat-miss caption "proposed-family-size". The first three should say DINA; the caption's stratifier is `admitted_n`. No numbers affected, but under the tracking rule this document is read by other repositories, and "GRF's family-size stratifier" is exactly the error the substitution exists to prevent.
- A criterion-matched FS comparator at 12.4% is still the only compute gap in the grid.
