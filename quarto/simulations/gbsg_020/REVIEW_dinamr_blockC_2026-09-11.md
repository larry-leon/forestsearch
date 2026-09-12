# REVIEW — `dinamr` Block C, the deferred Block B cell, and the GRF cost probes

- **Date:** 2026-09-11
- **Reviews:** `REPORT_dinamr_blockC_2026-09-11.md` (commit range `ab7d7dbb..322daebd`, 12 commits, committed not pushed, no `R/` change).
- **Campaign state:** `dinamr` is complete — 18 of 18 cells at 2,000 replicates, none deferred, none dropped. The handoff's expectation of a partial campaign is closed.
- **Governing statement, unchanged:** every coverage number is coverage of β(Ĥ) conditional on the proposed DINA family, over selected replicates. DINA fails the fixed-family condition, so nothing here certifies DINA. Block C is differentially null against a benefiting complement, not a global null.

---

## 1. Verification — accepted

- Gate 0 resolved from source, not assumed: `FS_S7_HR` calibrates `k_inter` to a target Cox HR inside the planted region (template:340); the region rule depends on `FS_S7_Z1Q` alone (template:344–346); `dgm_model <- "alt"` is a literal (template:488), so the harness's global-null path is unreachable. Calibration confirms HR(H) = 1.00 against HR(Hᶜ) = 0.657 / 0.721.
- Gate 2: 204 passes / 0 failures on Block C and again on Block B; corrected identity 1.67e-16 to 3.33e-16; γ inside [0.025, 0.028] on both joints.
- Amendment 3 holds at all six Block C cells: `n_true` `identical()` on every row, truth `all.equal()` at 1e-8, maximum absolute discrepancy 1.443e-15.
- Structural columns reported as structural on every DINA cell, not on account of HR 1.00.
- Walls taken from the driver's own per-render figures; Part A 4.557 h realized against 4.800 h projected (0.949) inside a 9 h ceiling.
- Gate 1's compute model was accurate to 0.972 all along; the Block A overrun was render overhead, and the overhead's shape is now corrected (flat 9 s combine, the rest in the two batch renders, scaling with prevalence and inversely with n).
- Part B: Gate 3 PASS, 5.7 minutes against a 1.5 h cap.

## 2. Block C coverage — the reading

Absolute levels, conditional on the proposed family and on selection, with Wilson intervals at the realized selected counts. Nominal 0.95.

| Cell | Selected | Field lower on β(Ĥ) | Field-s upper on β(Ĥᶜ) | Bonferroni joint | IJ two-sided |
|---|---|---|---|---|---|
| 12.4% · n 500 | 1,427 | 0.875 [0.856, 0.891] | 0.930 [0.915, 0.942] | 0.908 [0.891, 0.921] | 0.992 |
| 12.4% · n 1000 | 1,049 | 0.849 [0.826, 0.870] | 0.947 [0.931, 0.959] | 0.888 [0.867, 0.905] | 0.989 |
| 12.4% · n 1500 | 687 | 0.825 [0.795, 0.852] | 0.953 [0.935, 0.967] | 0.882 [0.856, 0.904] | 0.990 |
| 31% · n 500 | 1,860 | 0.906 [0.892, 0.918] | 0.903 [0.888, 0.915] | 0.903 [0.888, 0.915] | 0.992 |
| 31% · n 1000 | 1,852 | 0.933 [0.921, 0.944] | 0.932 [0.920, 0.943] | 0.933 [0.920, 0.943] | 0.994 |
| 31% · n 1500 | 1,776 | 0.947 [0.936, 0.957] | 0.936 [0.923, 0.946] | 0.941 [0.929, 0.951] | 0.997 |

- **The region bound degrades with n at 12.4% and improves at 31%.** 0.875 → 0.849 → 0.825 against 0.906 → 0.933 → 0.947. The 12.4% trend reverses Block A's at the same prevalence (0.902 → 0.915 → 0.929) and its Wilson intervals are disjoint end to end.
- **The complement bound improves with n at both prevalences** (0.930 → 0.953 and 0.903 → 0.936). The degradation is specific to the region bound at 12.4%.
- **Block C's region bound is the worst in the campaign.** 0.825–0.947 against Block A's 0.902–0.944, and below nominal at every cell.
- **IJ covers by width everywhere**, 0.989–0.997 on the region and ≈ 1.000 on the complement. Its behaviour carries no information about the family component here.
- **Bonferroni joint tracks the weaker of its two arms**: 0.882–0.941, below the FS reference at every cell.

### 2.1 What drives the 12.4% degradation — reading, not established

- Selection collapses with n: 0.714 → 0.524 → 0.344. By n 1500 barely a third of replicates select, so the summary conditions on an increasingly noise-selected subset.
- The bound becomes less conservative as n grows: bound / θ(Ĥ) rises 0.721 → 0.806 → 0.842. A bound far below the target covers often but uninformatively; as it tightens it crosses more.
- So the fall in coverage and the fall in the location shares are the same movement, not opposing ones.
- At 31% selection stays near 0.9, the conditioning is weak, and coverage improves with n in the ordinary way.
- **Distinguishing detection-conditioning from the family component needs the strata, which the report does not carry** (section 4.1).

### 2.2 The between-cell family ordering

- Block C has the tightest families in the campaign: 12.4% median 137 → 66 → 36, CV 1.422–1.458; 31% median 485.5 → 349 → 278.5, CV 1.008–1.083. FS on identical draws is median 1,223–1,300 at CV 0.036–0.054.
- The cell with the tightest family (12.4% n 1500, median 36) carries the worst region coverage (0.825); the cells with the largest families carry the best. That is the manuscript's ordering — the gap grows when the family is a tight surface-selected handful.
- It does not hold within a prevalence block: at 31% the family shrinks with n while coverage improves. Family size alone does not order the six cells, and prevalence is confounded with it.
- This continues to sit against Block A's within-cell finding, where the largest-family tertile covered worst. The two stratifications measure different things — between cells, family size indexes closeness to the enumerated space; within a cell, it indexes winner dominance.

### 2.3 The location split — the naming decision vindicated

| | 12.4% (n 500 → 1500) | 31% (n 500 → 1500) |
|---|---|---|
| Selection rate | 0.714 → 0.524 → 0.344 | 0.930 → 0.926 → 0.888 |
| Share of lower bounds ≥ 1.00 | 1.9% → 1.6% → 0.9% | 3.0% → 1.6% → 1.0% |
| Share ≥ 1.25 | 0.35% → 0.10% → 0.15% | 0.86% → 0.32% → 0.17% |

- Selection is frequent and admissible; the bound reaches a level that would assert harm on 0.9–3.0% of replicates, and 0.10–0.86% at HR 1.25, falling with n at both prevalences.
- Had these been labelled false selections, the record would have reported a 34–93% error rate for what is 1–3% by location.

### 2.4 FS beside DINA, criterion-matched at 31%

- Region bound: DINA 0.906 / 0.933 / 0.947 against FS 0.973 / 0.956 / 0.967. FS is above nominal at all three, DINA below at all three; the gaps are 6.8, 2.3 and 1.9 points and the Wilson intervals are disjoint at n 500 and n 1000.
- Complement bound: DINA 0.903 / 0.932 / 0.936 against FS 0.912 / 0.928 / 0.927 — overlapping, no separation.
- Selection: DINA 0.930 / 0.926 / 0.888 against FS 0.921 / 0.955 / 0.959.
- At 12.4% the comparison carries the criterion confound (`maxeffCons` ε 0.10) and is not read here.

## 3. Part B — what the probes establish

- Cost: 13.5–19.8 s median per replicate, rising with n and with prevalence; wall independent of family size (|ρ| ≤ 0.171). This is the opposite of DINA's family-size-driven profile.
- Peak 15.8–19.7 GB summed over the process tree at 12 workers (1.7–1.9 GB parent alone; the summed figure over-counts shared pages).
- Selection 1.0000 at all four harm corners and 0.9722 at the null corner.
- Recovery columns, p̂ and ρᶜ present and populated on the GRF path at every probe.
- The frontier band was empty on 0 of 180 replicates, and the recomputed eligible count matched the code's own `admitted_n` on 180 of 180.

### 3.1 The two-floor finding

- `dmin.grf = 0.0` is a DR-score pre-filter on the eligible set, consumed only by the native frontier select (`grf_subgroup_labels.R:358`), where the effect is a mean DR contrast in RMST units. The decision's premise is exactly right for that filter.
- The binding effect-scale floor on the re-selection path is `hr.threshold = 0.90` — the same floor DINA carries, because it is a property of the resolved admission set rather than of the engine (`forestsearch_helpers.R:1632–1635`, `:2345`; `forestsearch_main.R:2020–2030`, `:2442`; template:531).
- The decision stands as made. What changes is what it is a decision about: `dmin.grf = 0.0` does not leave GRF unfloored relative to DINA, it makes the DR pre-filter maximally permissive and leaves the binding decision to the shared 0.90.
- **The "GRF must not be called FS-analogous" line in the task document needs narrowing.** The criterion-scale difference is real at the DR pre-filter and absent at admission. The GRF campaign's framing should say which of the two it means, rather than carrying the unqualified claim.

### 3.2 The frontier-filter asymmetry (§8, open)

- The band cannot empty under this configuration: the documented mechanism needs a negative maximum over the eligible set, and the eligible set is non-negative on the DR path (`dmin` 0.0) and strictly positive on the re-selection path (hazard ratios).
- What empties is the floor, recorded as `admitted_n = 0L`.
- So the asymmetry has no reachable consequence at `dmin.grf = 0.0`. It is unreachable rather than acceptable — a different basis for the decision than the one §8 assumed, and it is contingent on `dmin.grf` staying at 0.0.
- CC left the decision open; correct.

### 3.3 GRF's `n_family` is a different object from DINA's — the item that most affects a GRF campaign

- GRF's family is nearly constant (776 / 784 / 853 at n 500; 828 / 834 / 838 at n 1500) and **identical across prevalences at matched n**, which differ only in `k_inter`, i.e. only in the outcome.
- CC's source reading: `.grf_dr_candidates(X, dr_scores, n_min)` enumerates on the covariate matrix subject to `n.min`, so the candidate set is a function of (X, `n_min`) and not of the outcome.
- The outcome-dependent narrowing is the admission step: the DR pool runs 760–955 while the admitted set runs 3–792 across corners.
- **Consequences, each to be verified from source before a GRF campaign:**
  - `n_family` on GRF appears to be the outcome-independent enumerated pool, whereas on DINA it is the surface-proposed set. If so, family-size comparisons across engines compare different objects, and the DINA-vs-GRF family contrast in any record needs restating in those terms.
  - Stratifying GRF on `n_family` would stratify on something nearly constant and outcome-independent. `admitted_n` is the analogue of DINA's family-size stratifier.
  - Handoff §3 calls GRF "the sharpest case of the family caveat" because "the forest itself determines which conjunctions qualify". That wording is about qualification, not enumeration, so the finding refines it rather than contradicting it — but it is the enumerated pool, not the qualifying set, that a record calling GRF's family "forest-generated" would be describing.
- **Which set MR resampling re-evaluates is the open question**: the enumerated pool or the admitted subset. It decides what "the family" is for the conditional estimand on GRF, and it is not answered in the report.

### 3.4 The unreproducible fit

- GRF fits are identical run-to-run within a session (3 of 3) but not across contexts: 0 of 179 detected rows matched between the diagnostic and the pipeline, differing by up to 75 and 131 candidates, with nearly identical distributions.
- `seedit` does reach the forest (`grf_helpers.R:74-85`), so the cause lies elsewhere — the `future` worker context and how `n.min = NULL` resolves are the untested candidates.
- Not chased, correctly: chasing it risks an `R/` change the task forbids.
- **For a GRF campaign this is more than a loose end.** DINA's Amendment 3 held bit-exactly, which is what let the same-draws pairing serve as the proof of comparability. If GRF fits are not reproducible across contexts, a GRF Gate 2 cannot assert the same pairing and needs a different verification design. This should be settled before a GRF campaign is specified, not during one.
- The band finding is unaffected: it is structural and held on all 180 replicates regardless of which context produced the fit.

## 4. Gaps in the report

- **4.1 No Block C stratified tables.** The report carries no coverage or bias by `n_family` tertile or by p̂ bin for the Block C cells, and no joint count table. Block C has the campaign's most volatile families (CV 1.008–1.458) and its most extreme selection gradient, so this is where the stratification bears most. Without it, the 12.4% degradation cannot be separated into detection-conditioning and family-generation components. They may be in the rendered `summary_dinamr.qmd` and simply not quoted.
- **4.2 No retained bias, error SD, marginal SD or mean SE for Block C.** The Block A decomposition (bias-driven at small n, variance-deficit at large n) cannot be repeated here without them.
- **4.3 The `SP` export in `campaign.sh`** is recorded as not fixed; the Part A run exported it explicitly. Fix pending.
- **4.4 The candidate-list identity** across prevalences is argued from source and from size identity, not measured as a symmetric difference on the stored lists.

## 5. Dispositions

- Report accepted. No re-run, no re-verification of committed work.
- No `R/` change proposed.
- The campaign's substantive open items go to the GRF go/no-go (3.1, 3.3, 3.4), not to a follow-up DINA task.
- §8 decisions after this report: `dmin.grf` closed (Larry, 2026-09-11); DINA/GRF field defaults still open and now informed by the full 18-cell grid; the frontier-filter asymmetry unreachable at `dmin.grf = 0.0` and left open; the commit pin still open.
