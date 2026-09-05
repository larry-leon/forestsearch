# TASK — MR de-biased gate vs Guo & He (2021) on the replicated Guo–He simulations

Date: 2026-09-04. Author: chat (spec). Executor: Claude Code. Approver: Larry.

## Protocol

- First action: copy this file to `dev/tasks/` and commit it. Do not push.
- Gates are stop-on-failure. A green gate needs no chat round trip; a red gate stops with a report naming what failed.
- No compute beyond the Stage 1 smoke test without Larry's explicit go on Stage 1's projection.
- Verify everything from the current source. Nothing in this document about function names, storage layouts, or settings is authoritative; read the files and quote them.
- Reports and verification records go in the directory of the Guo–He replication document, beside any existing `REPORT_*` files — not in `dev/tasks/`.
- Do not modify the Guo–He replication code, its bundles, or its committed numbers. Do not re-run anything already committed unless Stage 0 shows it is the only way to obtain per-replicate values (D5).

## Objective

A paired comparison, on identical simulated replicates, of three estimates of the treatment effect in the selected subgroup: (i) naive, (ii) the Guo–He (2021) proposed estimator and confidence interval as replicated, (iii) the multiplier-resampling (MR) two-term de-biased estimate with its infinitesimal-jackknife interval. The target for all three is θ_{k̂}, the true effect of the subgroup actually selected — the same conditional estimand in both papers. Report bias, empirical SD, mean model SE, the SE/SD ratio, coverage with Monte Carlo uncertainty, and interval length, per scenario.

Non-goals: no change to the MR method; no Andrews–Kitagawa–McCloskey comparator; no full bootstrap unless D4 says it is free; no new DGMs.

## Stage 0 — Discovery (no compute, no edits under R/)

0a. Locate the Guo–He replication: `grep -ril "guo" quarto/ dev/ R/ vignettes/ inst/ tests/` and any bundle directories. Identify the replication document(s), the DGM function(s), the scenario grid, the saved bundles, and the record of how the replication's numbers were matched to the paper's tables. Quote paths and the current commit hash. If more than one replication document exists and it is not evident which is current, STOP.

0b. From source, record: outcome type; K and the subgroup definitions (disjoint or overlapping); the selection rule (pure argmax, or argmax subject to a threshold); the within-subgroup estimator (difference in means, regression coefficient, other); orientation (benefit = larger is selected); the Guo–He bootstrap size, tuning-parameter rule and calibration; n per scenario; replicates per scenario; seed scheme; whether per-replicate values are stored (selected index, truth θ_{k̂}, naive estimate, G&H estimate and CI, seed) or only summaries; wall-clock per replicate and per cell.

0c. Locate the MR gate entry point and quote its signature: the function(s) that build the per-candidate effect vector and influence matrix, draw the multipliers, re-select per draw, form bias_sel and bias_fix, form the IJ variance with its finite-B correction, and return the interval. Determine whether the gate can be called with an externally supplied family (a list of membership indicators or logical columns), a pure argmax rule, benefit orientation, and the outcome type found in 0b. If the continuous/MD path is required, confirm it is present on the branch in use.

0d. Determine whether the replication's per-replicate data are reproducible from stored seeds (or stored outright), so MR can be evaluated on the identical replicates.

Gate 0 (stop-on-failure): 0a–0d resolved with quoted paths and signatures. Any of: ambiguous replication source; gate not callable on an external family without changing its arithmetic; continuous path absent when needed; per-replicate reproduction impossible → STOP and report which. Larry decides.

Output: `REPORT_mr_vs_guohe_stage0_<date>.md`.

## Stage 1 — Adapter, identities, smoke test

1a. Adapter. If the gate accepts an external family, write the adapter in the comparison `.qmd` and touch nothing under R/. If not, add one function under R/ (working name `mr_gate_fixed_family(data, family, outcome, treat, rule = "argmax", orientation, B)`) that calls the existing internals unchanged. R/ classification: adds code; moves nothing; changes no existing behaviour (the existing test suite must pass unchanged); changes no method. Anything requiring a change to the gate's arithmetic is a method proposal — STOP and write it up with justification; do not implement.

1b. Identities (machine-checked, tolerances stated in the record):
- Within every candidate, the influence values sum to zero and their sum of squares equals the robust (sandwich) variance of the within-subgroup estimator to relative tolerance 1e-6.
- With K = 1 the de-biased estimate equals the naive estimate up to Monte Carlo error at B draws, and V̂/σ²_D lies in [3.6, 4.4] at B = 5,000 (the separated-regime limit of the gate's variance).
- With all K candidates identical, the selection frequencies are exchangeable (chi-square test on the counts, p > 0.01).

1c. Smoke test: 5 replicates of one scenario, paired by seed, all three estimators; assert the naive column reproduces the replication's stored naive values exactly, and that MR ran on every replicate. Record wall-clock for MR and, if D5 applies, for G&H. Project the full grid: cells × replicates × seconds per replicate → wall-clock, stated with the worker count.

Gate 1: 1b identities pass; 1c runs and the naive identity holds; projection reported. STOP for Larry's compute go on the projection (this gate is stop-to-ask by design: it is the compute go/no-go).

Output: `REPORT_mr_vs_guohe_stage1_<date>.md`.

## Stage 2 — Full run (after the compute go)

Per scenario, paired replicates. Per-replicate record: seed; selected index under each method (assert identical when both use argmax on the same estimator; otherwise record both and the agreement rate); θ_{k̂}; naive est/CI; G&H est/CI; MR est/CI, se_ij, bias_sel, bias_fix, ij_source; re-selection frequencies p̂ for the top three candidates; M_eff; the (A6) diagnostic value; timings. Bundle naming per the repository's existing convention.

Gate 2: replicate counts complete per cell; MR NA count zero or documented per replicate; naive column equals the replication's stored naive values across all cells; if G&H was re-run, its cell summaries match the committed replication numbers within stated Monte Carlo tolerance.

## Stage 3 — Report

`REPORT_mr_vs_guohe_<date>.md` plus the rendered `.qmd`. Per scenario, rows naive / G&H / MR; columns mean estimate, bias against θ_{k̂}, empirical SD, mean SE, SE/SD, coverage with Wilson interval, mean and median length. Add: the fraction of the naive optimism retained by MR, against the tie constant 1 − 2^(−1/2) = 0.293 in tie scenarios; a per-replicate scatter of MR against G&H estimates; a cost table. Interpretation limited to what the numbers show. Findings go in the record; propose a task only if something blocks.

## Decisions required from Larry (defaults in brackets)

- D1 Scenarios: all replicated Guo–He scenarios [default] / a named subset.
- D2 MR rule: pure argmax on the replication's within-subgroup estimator, benefit orientation, no screen or consistency [default]; a second MR row using the package's default gate [default: no].
- D3 Multiplier draws: B = 5,000, centred Poisson [default].
- D4 Full bootstrap row: only if obtainable at zero extra cost from the replication's existing bootstrap [default: no].
- D5 If per-replicate G&H values are not stored: re-run G&H at the replication's settings — separate compute go on Stage 1's projection.

## Done means

Stage 3 report committed, Gate 2 record beside it, every identity in 1b listed with its measured value, and the branch left unpushed for Larry.
