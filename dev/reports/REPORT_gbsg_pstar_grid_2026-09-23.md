# REPORT — p\* calibration grid under a uniform marginal Cox HR of 0.75, GBSG application design

Task: `dev/tasks/TASK_gbsg_pstar_grid_2026-09-23.md` (commit `13573a74`). FS only, no MR, no `R/` change.
Driver `quarto/simulations/gbsg_app_null/run_gbsg_pstar_grid.R` (`8803326c`). Payloads, logs and the read-out script
`pstar_grid_findings.R` (output `logs/pstar_grid_findings.txt`) are in `quarto/simulations/gbsg_app_null/` (`fea641d7`).

## Result

Declaration rate = `mean(max Pcons >= p*)` over 5,000 replicates per cell, with a 90% Wilson interval (alpha = 0.10) and the count.

| p\* | Cell 1: arm re-randomized, P(treat) = 0.359 | Cell 2: observed arm fixed |
|---|---|---|
| 0.90 | 0.448 [0.437, 0.460] 2,241 | 0.432 [0.421, 0.444] 2,161 |
| 0.95 | 0.270 [0.260, 0.280] 1,350 | 0.255 [0.245, 0.265] 1,273 |
| 0.96 | 0.227 [0.218, 0.237] 1,137 | 0.215 [0.206, 0.225] 1,075 |
| 0.97 | 0.185 [0.176, 0.194] 924 | 0.171 [0.163, 0.180] 856 |
| 0.98 | 0.137 [0.130, 0.146] 687 | 0.122 [0.115, 0.130] 610 |
| 0.99 | **0.081 [0.075, 0.087] 404** | **0.074 [0.068, 0.080] 370** |

- **The rate first falls to or below 0.10 at p\* = 0.99 in both cells.** The bracketing points are 0.98 (0.137 / 0.122,
  both intervals above 0.10) and 0.99 (0.081 / 0.074, both below).
- **That is the finest resolution the design supports.** `Pcons` is reported to two decimals
  (`pconsistency.digits = 2`), so any p\* in (0.98, 0.99] declares exactly as 0.99 does, and no finer p\* is defined.
  The 0.925 grid point is dropped for the same reason. The full CDF at 0.01 steps over [0.50, 1.00] is in
  `logs/pstar_grid_findings.txt`. At p\* = 1.00 the rate is still 0.026 / 0.022, because a `Pcons` that rounds to 1.00
  clears it.
- **At p\* = 0.90:** 0.448 (Cell 1) and 0.432 (Cell 2), beside FŴ₀.₁₀(0.75) = 0.651 and the 2026-09-23
  super-population figure of 0.329. This is the first matched-family comparison: same covariates, same allocation, same
  n as the family FŴ was computed on.

**Reading.**
- On the application's own covariates the executed p\* = 0.90 screen declares in about 44% of null trials. That is higher
  than the 33% the super-population design gave and still below FŴ's 0.651.
- Holding the rate to 10% under a uniform marginal HR of 0.75 takes p\* = 0.99, the top of the two-decimal scale short
  of 1.00. The rate falls steeply through that range, from 0.137 to 0.081 between 0.98 and 0.99 in Cell 1.
- The two cells agree to about 1.5 points at every grid point. Cell 2 is slightly lower throughout, so re-randomizing
  the arm contributes little. Which subjects are treated barely moves the rate; the null outcome noise and the family
  drive it.

## Cost (a result, not only a gate)

- **One run per cell gives the whole p\* grid.** Each replicate runs once at the floor p\* = 0.50. The rate at every
  p\* ≥ 0.50 is then read off the per-replicate `max Pcons`, so no run is needed per grid point. Gate B shows this is
  exact (below).
- Machine: pop-os, 128 logical cores, 48 multisession workers, `forestsearch()` sequential inside each worker.

| run | replicates | workers | wall clock |
|---|---|---|---|
| Cell 1 (floor) | 5,000 | 48 | 718.5 s = 12.0 min |
| Cell 2 (floor) | 5,000 | 48 | 735.2 s = 12.3 min |
| **both cells** | **10,000** | 48 | **1,453.7 s = 24.2 min** |
| Gate B (p\* 0.90 and floor, 200 each) | 400 | 48 | 67.0 s |
| Gate C (floor) | 10 | 48 | 7.2 s |

- The mean time per fit is 5.5 s. Throughput was 6.9–7.1 worker-seconds per replicate, against the 14 min per cell
  projected from Gate B's loaded batches.

## Gates

- **Gate A: PASS, no recalibration.** `k_treat` = 1.048469, as calibrated on `df_super` (marginal HR 0.750000 there).
  - On `df_source` (the observed 686 rows, both potential outcomes under a common extreme-value error, 20 draws
    stacked) the marginal Cox HR is **0.7463**, 0.50% below 0.750. Across five seeds for the draws it ranges
    0.7443–0.7495.
  - `flag_harm` is identically 0.
  - The patient-level conditional HR is uniform at 0.676, recorded for vocabulary only; it is not the target.
- **Gate B: PASS.** Two runs of 200 replicates on Cell 1 with the same seeds, one at p\* = 0.90 and one at 0.50.
  - `I(max Pcons >= 0.90)` at the floor equals the p\* = 0.90 declaration indicator on **200 / 200** replicates; no
    replicate disagrees. Both runs have 87 declarations.
  - `n_candidates_total` matches on every replicate. There were 0 errors.
- **Gate C: PASS.** 10 floor replicates took 7.2 s on 48 workers.
  - Ten replicates cannot keep 48 workers busy, so the go was sized on Gate B's loaded throughput: about 14 min per
    cell, against a measured 12.0–12.3 min.

## Candidate family and truth on the application's covariates

- **Family.** The prior driver's route on `df_source` gives 66 cut columns and 2,211 enumerated combinations. It drops
  142 empty, 0 minp, 126 rmin and 323 below size; 1,620 are kept, 113 duplicate an earlier membership, and
  **M = 1,507** remain.
  - The application's `family_size_prereduction` = 1,744 comes from a different route: the MR block of `forestsearch()`,
    which keeps every combination of up to two factors with at least 60 members, with no other floor and no merging of
    identical memberships.
  - That route on `df_source` gives **1,744 exactly**, so the fixed baseline reproduces the application's covariates.
  - The 13.6% gap between 1,507 and 1,744 comes entirely from the empty / rmin floors and the membership merge.
- **Truth: the uncensored population marginal Cox HR of each candidate on `df_source`.** It uses the potential-outcome
  construction above.
  - Over the 1,507 candidates: **min 0.674, median 0.733, max 0.798**.
  - **364 of 1,507 (24.2%) lie above HR 0.75.**

## Post-conditions

1. PASS. The marginal Cox HR on `df_source` is 0.7463, within 1% of 0.750; `flag_harm` is identically 0.
2. PASS. `args_call_all` shows `stop_threshold = NULL` and every setting in §1 on **10,410 / 10,410 fits** (gates and
   both cells). `sg_focus` is stored under its standard name `hrMaxSG`; `effect_measure` is NULL.
3. PASS. Gate B agrees exactly on 200 / 200 replicates, and `n_candidates_total` matches on each.
4. PASS. The Gate C wall clock is recorded; the cells took 719 s and 735 s, far inside the 2 h abort.
5. PASS. Errors are 0 of 5,000 in each cell. Each denominator is 5,000, the replicate count.
6. PASS. Cell 2's treated count is 246 in 5,000 / 5,000 replicates. Cell 1's mean treated fraction is 0.3586.
   - `rand_ratio = 246/440` gives n·p = 246 exactly, so Cell 1's count is also 246 in every replicate. Only which
     subjects are treated varies.
7. PASS. Replicate 1 re-run sequentially under `RNGkind("L'Ecuyer-CMRG")` reproduces the parallel result, and equals
   Cell 1 replicate 1.
   - Cell 1's rows 1–200 are also identical to the Gate B floor run, and rows 1–10 to Gate C.
8. PASS. The mean event rate is 0.4376 (Cell 1) and 0.4377 (Cell 2), against GBSG's 0.436.
9. The catalogue pin is checked at commit time with `scripts_dinamr/check_current_status.sh` (see the closeout commits).
10. `R/` is unmodified. Run output is written only to `quarto/simulations/gbsg_app_null/` and `dev/tasks/`. §8 itself
    requires three writes outside those directories: this report, in `dev/reports/` as on 2026-09-23; the catalogue line
    in `quarto/simulations/gbsg_020/status_curated.md`, together with the "Updated" line of `current_status_regen.R`; and
    the regenerated `current_status.md`.

## Pins

- HEAD at the cells: `8803326c`. Gates A–C ran at `13573a74` on the identical driver; the driver was committed unchanged
  after them, and Cell 1 reproduces Gates B and C row for row.
- forestsearch 0.3.5.9000, installed build 2026-09-23 02:32:08 UTC. That postdates the last `R/` commit (`845fea56`), so
  the installed package is HEAD's `R/`.
- R 4.6.1, x86_64-pc-linux-gnu.
- Seeds: DGM 8316951; replicate b uses `simulate_from_dgm` seed = `seedit` = 8316951 + b; the truth draws use 20260923.
- `simulate_from_dgm(baseline = "fixed", n = NULL, analysis_time = Inf, cens_adjust = 0)`.
- Pre-existing untracked files (the `actg175/binary_020` d5000 directories and smoke html files, and
  `gbsg_020/scripts_dinamr/logs/nullmr_findings.err`) were never staged.
