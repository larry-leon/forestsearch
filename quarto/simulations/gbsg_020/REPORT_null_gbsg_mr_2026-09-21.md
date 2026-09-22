# REPORT — campaign (B): the structural null at c1 0.90 / c2 0.80, with MR post-selection inference

**Task:** `dev/tasks/TASK_null_gbsg_mr_2026-09-21.md` (committed `25798342`) · **Opened:** 2026-09-21,
run overnight into 2026-09-22 · **Branch:** `feature/glm-extension` · **Study directory:**
`quarto/simulations/gbsg_020` · **Campaign tag:** `nullmr` · **Machine:** pop-os.

**What this is.** The six `nullid` cells — no planted region, uniform treatment benefit at
super-population marginal Cox HR 0.657 and 0.721, n 500 / 1000 / 1500, 2,000 replicates, FS / DINA /
GRF on identical draws, `effMaxSG` at ε = 0.20, seeds `8316951 + sim_id` — at `nullid`'s screen
**(c1, c2, p⋆) = (0.90, 0.80, 0.90)** (threshold knobs unset), now **with MR post-selection
inference on**. Identification is `nullid`'s; what this adds is where the selection-adjusted bounds on
a false region sit, and how they cover their targets **β(Ĥ)** and **β(Ĥᶜ)** (`betaHhat_H` /
`betaHhat_Hc`: the super-population marginal Cox HRs of the selected region and its complement).
With no planted region there is no θ†(H) and no oracle. **DINA and GRF results are conditional on
the proposed family.** **Every declaration is a false declaration.**

**Scope held.** No edit to `R/`; no install. **The template is untouched** — the Step 3.3 guard was
not needed (the MR path runs under `FS_S7_DGM=null` as it stands). `nullid`, `nullc125` and every
committed cell untouched; nothing committed was re-run. No R CMD check, vignette build or test
suite. Nothing fetched, pulled or pushed. Every `git add` named its paths. `FS_S7_FB=none`.

**Outcome.** 18 of 18 runs; **Gate A PASS on every run, Gate C PASS on every cell; no halt record.**
Identity gate PASS. Grid wall 19,858 s of render (driver 19,865 s end to end, **5.52 h**) against a
**5.55 h** projection.

---

## 0. Preconditions and build

| check | result |
|---|---|
| host / branch | `pop-os` / `feature/glm-extension` |
| HEAD contains `c4baf796` | yes |
| tracked modifications / R, Rscript, quarto processes | none / none |
| task committed alone | `25798342` |
| last commit touching `R/` | `f9b794f6`, 2026-09-18T20:02:11-07:00 (= 2026-09-19 03:02 UTC) |
| installed forestsearch | **0.3.5.9000**, `Built: R 4.6.1; ; 2026-09-19 03:54:02 UTC; unix` — newer than the last `R/` commit; no install |
| R | 4.6.1 (2026-06-24) |
| machine | AMD Ryzen Threadripper PRO 5995WX, 64 physical cores, 251 GB; 64 workers |

## 1. MR settings and coverage definitions (Step 1)

Full record with path:line: `scripts_dinamr/logs/nullmr_step1_record.md` (commit `40db0402`).

- **Per engine, as committed.** Common to all three: `FS_S7_MR=TRUE`, `FIELD_COMPLEMENT=TRUE`,
  `FIELD_SCALEC=selected`, `FIELD_DECOMP=TRUE`, `IJ_RESIDUAL=two_term`, `FB=none`. **FS**
  (`scripts_p12x20/campaign_p12x20.sh:24-27`, cross-checked against cert20's recorded environment,
  `REPORT_cert20_2026-09-08.md:11-15`): adds `RETURN_RESEL=TRUE`, leaves `FIELD_RECOV` unset
  (FALSE). **DINA** (`scripts_dinamr/campaign.sh:9-11`) and **GRF** (`grfmr.sh:19-21`,
  `grfmrC.sh:17-19`, identical): add `FIELD_RECOV=TRUE`, leave `RETURN_RESEL` unset (TRUE).
  `UNIFORM`, `WINNER_ROWS` unset everywhere.
- **No disagreement within an engine.** `FIELD_RECOV` postdates cert20, so its absence there is not a
  disagreement; cert20's and the dinamr / grfmr bundles' own meta confirm the rest.
- **Template MR literals:** 5,000 draws (T:545), `ci_method = "field"` (T:673), confirm rule
  `"point"`, `t_confirm` NULL, `include_complement = TRUE`; field R_out / R_in = 1000 / 500 (package
  defaults, not forwarded).
- **Winner-only / winner-floor:** no knob turns them on — they are recorded whatever
  `IJ_RESIDUAL` selects — and **no committed survival driver set `FS_S7_WINNER_ROWS`**. Not evaluated
  or reported.
- **Coverage** is `summary_grfmr.qmd`'s `cov-fns` and `wilson-fn` chunks, **extracted from the
  committed .qmd and eval'd verbatim** by `scripts_dinamr/nullmr_findings.R` (the mechanism of
  `grfmr_tables.R:12-21`); the same definitions are in `summary_cert20.qmd:74-89, :131-161` and
  `summary_dinamr.qmd`. Field lower `betaHhat_H >= fld_H_lo1s`; field-s upper `betaHhat_Hc <=
  fld_Hc_up1s_s`; Bonferroni pair `joint_s_bonf`; IJ two-term two-sided on `mr_*_lo / mr_*_hi`; all
  on `covs_cell`'s row set. Naive by `fs_sim_bias_coverage()` (one-sided on the exposed side, and
  two-sided).

## 2. Driver and gates (Step 2)

- **`scripts_dinamr/nullmr.sh`** — transplant of `nullthr.sh`. Changes: `FS_S7_MR=TRUE` plus the
  per-engine knobs above; tag `nullmr`; threshold knobs never passed; timeouts 6 h per render, 24 h
  for the campaign. Mechanical consequences: the stem has no `_nomr` token and no threshold tag; the
  gates are the nullmr ones; the halt record names this task. Everything else is `nullthr.sh`'s
  (64 workers, inherited `FS_S7_*` unset first, halt-and-continue, commit-per-cell by named paths).
- **Gate A** (`nullmr_gateA.R`; 53 checks FS / 49 DINA / 50 GRF): `nullthr_gateA.R` with the MR
  checks inverted — MR meta (`mr_inference` TRUE, `ci_method` field, 5,000 draws, the knob values,
  `field_recovery` per engine); `mr_ok ∈ {0,1}` and 0 on every non-declaring row; **every
  `mr_*` / `fld_*` / `n_family` / `p_hat_*` column NA on every non-declaring row**; **wherever
  `mr_ok == 1`, naive (`nv_H/Hc_est, lo, hi, se`), IJ two-term (`mr_H/Hc_est, lo, hi, se_ij`),
  `fld_H_lo1s`, `fld_Hc_up1s_s` and the Bonferroni pair all finite**; every full-bootstrap product
  NA; meta c1 0.90 / c2 0.80 / p⋆ 0.90 and resolved c1 / c2 equal to them **on every replicate**
  (non-NA required); `betaHhat_H` and `betaHhat_Hc` finite on every declaring replicate. Before use,
  the finite-products invariant was checked against the 43 committed MR bundles (dinamr, grfmr,
  cert20): it holds on every one.
- **Three `nullthr_gateA.R` checks read columns the template fills only on its MR-off branch** —
  `nv_H_lo1s` (T:1408-1409), `p_sel` and `p_max_qual` (T:1391-1395). With MR on they are NA by
  construction, so those checks became "NA on every replicate"; the unadjusted-estimate pairing keeps
  estimate and SE. Recorded here for review.
- **Gate C** (`nullmr_gateC.R`, 16 checks): `nullthr_gateC.R` as is, with the bundle path's `_nomr`
  token dropped (the only change). `itt_est` / `itt_se` `identical()` across identifiers on all
  2,000 rows, same `k_treat`, same truth.
- Tooling commit `40db0402`.

## 3. Smoke, identity gate, projection, grid (Step 3)

**Smoke** (`scripts_dinamr/logs/nullmrsmoke*`, commit `e13e3484`): `null0657_n500` and
`null0721_n1500`, `sim_id` 1–20, three identifiers, MR on, through the real driver. Gate A PASS × 6,
Gate C PASS × 2. The MR path ran under the null without failure, so **the one allowed template change
(Step 3.3) was not needed and not made.**

**Identity gate** (`scripts_dinamr/nullmr_identity.R`, log `logs/nullmr_identity.txt`): the MR-on
smoke of `null0657_n500` against `nullc125`'s committed MR-off render of the same replicates
(`*_nomr_nullc125inertunset_quickrun_res_1_20.rds`, same machine, same screen). The excluded classes
were fixed in the script before either output was read: timing; MR / field / IJ products (`mr_*`,
`fld_*`, `ij_source`, `n_family`, `p_hat_*`); the naive block `nv_*`, which with MR on is MR's own
`g$naive` and with MR off a template Cox refit (printed, not gated); and the MR-off-only fields
`p_sel`, `p_max_qual`, `nv_H_lo1s`. **Identical on 46 of 46 gated columns and on the truth object for
all three identifiers** — declarations 6 / 9 / 20 of 20, the same labels, sizes and rules, `betaHhat`,
classification, band and family counts, `maxT`, resolved thresholds, `itt_*`. **PASS.**

- **Printed, not gated:** `nv_H_est` / `nv_Hc_est` are identical. The naive **SE and interval
  differ**: MR's `se_wald` against the MR-off `coxph` model SE, up to 0.012 on the log scale for
  `nv_H_se` and 0.0024 for `nv_Hc_se`. So "naive" below (MR's) and `nullid`'s unadjusted bound share
  the point estimate but not the SE.

**Projection** (`scripts_dinamr/nullmr_project.R`, log `logs/nullmr_smoke_projection.txt`). Per
identifier, the per-replicate cost is split by declaration, because MR runs only on declaring
replicates. The declaring and non-declaring means were measured at n 500 and at n 1500 (linear
between) and weighted by each cell's `nullid` declaration rate. Measured n-scaling of the
all-replicate mean, n 1500 / n 500: FS 2.73, DINA 1.63, GRF 1.49 — measured, not assumed. The load
factor at 64 workers is the one `nullc125` used (63w 1.17 at n 500 → 3.27 at n 2000, linear).
**Projection 5.55 h** (3.08 h with no load factor, 9.42 h at 3.27 throughout). It was under 16 h with
every gate green, so the grid ran straight away.

**Grid:** 18 renders in `nullid`'s cell order; driver log `scripts_dinamr/logs/nullmr.driver.log`.
Cell commits `dbf210dc`, `0c87e7f3`, `fb9bc481`, `469dd838`, `1c232574`, `5213222b`. No failure, no
halt, no retry. Measured against the projection: the n 500 cells ran at 1.56× and 1.64× (the load
factor was light there), the n 1000 cells at 0.96× and 0.94×, and the n 1500 cells at 0.75× and
0.75×. The total came to 5.52 h against 5.55 h.

---

## 4. Results

All tables are produced by `scripts_dinamr/nullmr_findings.R`; the captured output is
`scripts_dinamr/logs/nullmr_findings.md`, reproduced here.

### Table 1 — coverage of β(Ĥ) and β(Ĥᶜ), over declaring replicates

Each entry is a rate with its Wilson 95% interval, and k is the number of declaring replicates it
rests on. Field, field-s, Bonferroni and IJ come from `covs_cell`, on its row set (declared, both
targets finite, both one-sided field bounds finite). Naive comes from `fs_sim_bias_coverage()`: the
one-sided 95% rate on the exposed side, and the two-sided rate. IJ two-term is two-sided, as the
committed reports define it. The nominal level is 0.95 throughout (the Bonferroni pair is jointly
0.95).

| cell | identifier | declared | field lower β(Ĥ) | field-s upper β(Ĥᶜ) | Bonferroni pair (joint) | IJ two-term β(Ĥ) 2-sided | IJ two-term β(Ĥᶜ) 2-sided | naive β(Ĥ) 1-sided lower | naive β(Ĥᶜ) 1-sided upper | naive β(Ĥ) 2-sided | naive β(Ĥᶜ) 2-sided |
|---|---|---|---|---|---|---|---|---|---|---|---|
| null0657_n500 | FS | 894 | 0.9564 [0.9409, 0.9679] (k=894) | 0.9329 [0.9146, 0.9475] (k=894) | 0.9508 [0.9346, 0.9631] (k=894) | 0.9922 [0.9839, 0.9962] (k=894) | 1.0000 [0.9957, 1.0000] (k=894) | 0.0000 [0.0000, 0.0043] (k=894) | 0.8400 [0.8146, 0.8626] (k=894) | 0.0000 [0.0000, 0.0043] (k=894) | 0.9072 [0.8864, 0.9245] (k=894) |
| null0657_n500 | DINA | 1116 | 0.8629 [0.8415, 0.8818] (k=1116) | 0.9382 [0.9225, 0.9509] (k=1116) | 0.8996 [0.8806, 0.9159] (k=1116) | 0.9928 [0.9859, 0.9964] (k=1116) | 1.0000 [0.9966, 1.0000] (k=1116) | 0.2115 [0.1885, 0.2364] (k=1116) | 0.8970 [0.8777, 0.9134] (k=1116) | 0.3728 [0.3449, 0.4015] (k=1116) | 0.9462 [0.9314, 0.9580] (k=1116) |
| null0657_n500 | GRF | 1883 | 0.9426 [0.9312, 0.9523] (k=1883) | 0.9368 [0.9249, 0.9469] (k=1883) | 0.9379 [0.9260, 0.9479] (k=1883) | 0.9920 [0.9869, 0.9952] (k=1883) | 1.0000 [0.9980, 1.0000] (k=1883) | 0.1354 [0.1207, 0.1516] (k=1883) | 0.7488 [0.7287, 0.7679] (k=1883) | 0.3218 [0.3011, 0.3433] (k=1883) | 0.8423 [0.8251, 0.8580] (k=1883) |
| null0721_n500 | FS | 1356 | 0.9735 [0.9635, 0.9808] (k=1356) | 0.9226 [0.9071, 0.9356] (k=1356) | 0.9528 [0.9402, 0.9629] (k=1356) | 0.9926 [0.9865, 0.9960] (k=1356) | 1.0000 [0.9972, 1.0000] (k=1356) | 0.0000 [0.0000, 0.0028] (k=1356) | 0.8038 [0.7819, 0.8241] (k=1356) | 0.0007 [0.0001, 0.0042] (k=1356) | 0.8768 [0.8583, 0.8933] (k=1356) |
| null0721_n500 | DINA | 1480 | 0.8703 [0.8522, 0.8864] (k=1480) | 0.9216 [0.9068, 0.9342] (k=1480) | 0.9014 [0.8851, 0.9155] (k=1480) | 0.9926 [0.9867, 0.9958] (k=1480) | 1.0000 [0.9974, 1.0000] (k=1480) | 0.2277 [0.2071, 0.2498] (k=1480) | 0.8480 [0.8288, 0.8654] (k=1480) | 0.3554 [0.3314, 0.3801] (k=1480) | 0.9068 [0.8909, 0.9205] (k=1480) |
| null0721_n500 | GRF | 1969 | 0.9426 [0.9314, 0.9520] (k=1969) | 0.9223 [0.9096, 0.9333] (k=1969) | 0.9335 [0.9216, 0.9437] (k=1969) | 0.9929 [0.9881, 0.9958] (k=1969) | 0.9995 [0.9971, 0.9999] (k=1969) | 0.1762 [0.1600, 0.1937] (k=1969) | 0.7232 [0.7030, 0.7425] (k=1969) | 0.3530 [0.3322, 0.3743] (k=1969) | 0.8207 [0.8032, 0.8370] (k=1969) |
| null0657_n1000 | FS | 637 | 0.8791 [0.8515, 0.9022] (k=637) | 0.9655 [0.9483, 0.9771] (k=637) | 0.9294 [0.9068, 0.9468] (k=637) | 0.9780 [0.9634, 0.9869] (k=637) | 1.0000 [0.9940, 1.0000] (k=637) | 0.0000 [0.0000, 0.0060] (k=637) | 0.9089 [0.8841, 0.9289] (k=637) | 0.0000 [0.0000, 0.0060] (k=637) | 0.9576 [0.9390, 0.9707] (k=637) |
| null0657_n1000 | DINA | 593 | 0.8061 [0.7723, 0.8359] (k=593) | 0.9646 [0.9465, 0.9767] (k=593) | 0.8685 [0.8389, 0.8933] (k=593) | 0.9865 [0.9736, 0.9931] (k=593) | 1.0000 [0.9936, 1.0000] (k=593) | 0.1906 [0.1610, 0.2241] (k=593) | 0.9410 [0.9190, 0.9573] (k=593) | 0.3474 [0.3102, 0.3866] (k=593) | 0.9646 [0.9465, 0.9767] (k=593) |
| null0657_n1000 | GRF | 1663 | 0.9405 [0.9280, 0.9509] (k=1663) | 0.9525 [0.9412, 0.9617] (k=1663) | 0.9465 [0.9346, 0.9563] (k=1663) | 0.9874 [0.9808, 0.9917] (k=1663) | 1.0000 [0.9977, 1.0000] (k=1663) | 0.0794 [0.0673, 0.0934] (k=1663) | 0.8172 [0.7979, 0.8350] (k=1663) | 0.2387 [0.2189, 0.2598] (k=1663) | 0.8827 [0.8664, 0.8973] (k=1663) |
| null0721_n1000 | FS | 1234 | 0.9408 [0.9263, 0.9527] (k=1234) | 0.9611 [0.9488, 0.9705] (k=1234) | 0.9473 [0.9334, 0.9585] (k=1234) | 0.9903 [0.9831, 0.9944] (k=1234) | 1.0000 [0.9969, 1.0000] (k=1234) | 0.0000 [0.0000, 0.0031] (k=1234) | 0.8395 [0.8180, 0.8590] (k=1234) | 0.0000 [0.0000, 0.0031] (k=1234) | 0.9206 [0.9042, 0.9344] (k=1234) |
| null0721_n1000 | DINA | 1152 | 0.8785 [0.8583, 0.8961] (k=1152) | 0.9401 [0.9249, 0.9524] (k=1152) | 0.9019 [0.8834, 0.9178] (k=1152) | 0.9939 [0.9875, 0.9971] (k=1152) | 0.9991 [0.9951, 0.9998] (k=1152) | 0.2943 [0.2687, 0.3212] (k=1152) | 0.9054 [0.8871, 0.9210] (k=1152) | 0.4375 [0.4091, 0.4663] (k=1152) | 0.9505 [0.9364, 0.9616] (k=1152) |
| null0721_n1000 | GRF | 1886 | 0.9438 [0.9325, 0.9533] (k=1886) | 0.9401 [0.9285, 0.9499] (k=1886) | 0.9390 [0.9273, 0.9490] (k=1886) | 0.9899 [0.9843, 0.9935] (k=1886) | 1.0000 [0.9980, 1.0000] (k=1886) | 0.1458 [0.1306, 0.1625] (k=1886) | 0.7471 [0.7270, 0.7662] (k=1886) | 0.3388 [0.3178, 0.3605] (k=1886) | 0.8547 [0.8381, 0.8699] (k=1886) |
| null0657_n1500 | FS | 337 | 0.7923 [0.7458, 0.8322] (k=337) | 0.9585 [0.9315, 0.9751] (k=337) | 0.8991 [0.8623, 0.9269] (k=337) | 0.9288 [0.8962, 0.9517] (k=337) | 0.9970 [0.9834, 0.9995] (k=337) | 0.0000 [0.0000, 0.0113] (k=337) | 0.9080 [0.8724, 0.9344] (k=337) | 0.0000 [0.0000, 0.0113] (k=337) | 0.9496 [0.9207, 0.9683] (k=337) |
| null0657_n1500 | DINA | 278 | 0.7626 [0.7092, 0.8088] (k=278) | 0.9568 [0.9261, 0.9751] (k=278) | 0.8273 [0.7785, 0.8672] (k=278) | 0.9928 [0.9742, 0.9980] (k=278) | 1.0000 [0.9864, 1.0000] (k=278) | 0.0827 [0.0558, 0.1211] (k=278) | 0.9424 [0.9086, 0.9643] (k=278) | 0.2914 [0.2411, 0.3473] (k=278) | 0.9604 [0.9305, 0.9778] (k=278) |
| null0657_n1500 | GRF | 1129 | 0.9389 [0.9234, 0.9514] (k=1129) | 0.9610 [0.9481, 0.9708] (k=1129) | 0.9486 [0.9342, 0.9601] (k=1129) | 0.9885 [0.9804, 0.9933] (k=1129) | 1.0000 [0.9966, 1.0000] (k=1129) | 0.0230 [0.0158, 0.0335] (k=1129) | 0.8450 [0.8227, 0.8649] (k=1129) | 0.1541 [0.1342, 0.1763] (k=1129) | 0.9221 [0.9049, 0.9363] (k=1129) |
| null0721_n1500 | FS | 969 | 0.9340 [0.9165, 0.9479] (k=969) | 0.9370 [0.9200, 0.9507] (k=969) | 0.9474 [0.9315, 0.9597] (k=969) | 0.9938 [0.9866, 0.9972] (k=969) | 1.0000 [0.9961, 1.0000] (k=969) | 0.0000 [0.0000, 0.0039] (k=969) | 0.8338 [0.8091, 0.8560] (k=969) | 0.0000 [0.0000, 0.0039] (k=969) | 0.9071 [0.8872, 0.9238] (k=969) |
| null0721_n1500 | DINA | 748 | 0.8489 [0.8215, 0.8728] (k=748) | 0.9505 [0.9326, 0.9639] (k=748) | 0.8917 [0.8674, 0.9120] (k=748) | 0.9947 [0.9863, 0.9979] (k=748) | 1.0000 [0.9949, 1.0000] (k=748) | 0.2968 [0.2652, 0.3305] (k=748) | 0.9131 [0.8907, 0.9312] (k=748) | 0.4706 [0.4351, 0.5064] (k=748) | 0.9639 [0.9480, 0.9751] (k=748) |
| null0721_n1500 | GRF | 1690 | 0.9574 [0.9467, 0.9660] (k=1690) | 0.9544 [0.9434, 0.9634] (k=1690) | 0.9550 [0.9441, 0.9639] (k=1690) | 0.9923 [0.9869, 0.9955] (k=1690) | 1.0000 [0.9977, 1.0000] (k=1690) | 0.1201 [0.1055, 0.1365] (k=1690) | 0.7615 [0.7406, 0.7812] (k=1690) | 0.3148 [0.2931, 0.3373] (k=1690) | 0.8592 [0.8418, 0.8749] (k=1690) |

(For FS at `null0721_n1000`, the naive rows' Wilson lower limit prints as −0.0000; that is the
function's floating-point output at p = 0.)

- **The naive bound almost never covers the truth in a false region.** For FS the one-sided naive
  lower bound covers β(Ĥ) on 0 of 5,427 declarations across the six cells. On DINA the rate is
  0.08–0.30 and on GRF 0.02–0.18. Selection puts the unadjusted estimate far above the region's true
  (benefit-side) HR.
- **The field lower bound restores much of the coverage, but not uniformly.** GRF: 0.939–0.957 in
  every cell. FS: 0.934–0.974 in four cells, but **0.879 at `null0657_n1000` and 0.792 at
  `null0657_n1500`**. DINA: **0.763–0.879, below nominal in every cell** — conditional on DINA's
  proposed family.
- **The field-s upper bound on β(Ĥᶜ) is 0.92–0.97 on all three identifiers**, close to nominal.
- **The Bonferroni pair tracks the weaker of its two members.** FS 0.899–0.953, GRF 0.934–0.955,
  DINA 0.827–0.902.
- **IJ two-term (two-sided) over-covers.** β(Ĥ) is covered at 0.929–0.995 and β(Ĥᶜ) at 0.997–1.000;
  Table 2 shows the IJ interval is wide.

### Table 2 — where the lower bounds on β(Ĥ) sit (HR scale)

These are locations against 1.00 and 1.25, not tests.

- **Bounds.** Naive and IJ two-term use the one-sided 95% Gaussian bound
  `exp(log est − z_0.95 · SE)` that `fs_sim_bias_coverage()` scores, with SE `nv_H_se` (MR's
  `se_wald`) and `mr_H_se_ij` respectively. Field uses `fld_H_lo1s`. For the Bonferroni pair, the
  table uses its lower member `fld_joint_s_bonf_loH` (97.5%).
- **Shares.** *cond.* is over declaring replicates with the bound finite (k). *uncond.* is over all
  2,000, where a replicate that declares nothing scores 0.
- **Reference row.** The `nullid` row is `nullid`'s committed unadjusted within-region bound
  (`nv_H_lo1s`, from its MR-off `nv_H` refit), printed beside naive for reference.

| cell | identifier | product | k | median LB | LB ≥ 1.00 cond. | LB ≥ 1.25 cond. | LB ≥ 1.00 uncond. | LB ≥ 1.25 uncond. |
|---|---|---|---|---|---|---|---|---|
| null0657_n500 | FS | naive | 894 | 0.848 | 101/894 = 0.1130 | 7/894 = 0.0078 | 101/2000 = 0.0505 [0.0417, 0.0610] | 7/2000 = 0.0035 [0.0017, 0.0072] |
| null0657_n500 | FS | nullid nv_H (MR off, ref.) | 894 | 0.843 | 93/894 = 0.1040 | 7/894 = 0.0078 | 93/2000 = 0.0465 [0.0381, 0.0566] | 7/2000 = 0.0035 [0.0017, 0.0072] |
| null0657_n500 | FS | IJ two-term | 894 | 0.434 | 0/894 = 0.0000 | 0/894 = 0.0000 | 0/2000 = 0.0000 [0.0000, 0.0019] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n500 | FS | field | 894 | 0.402 | 1/894 = 0.0011 | 0/894 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n500 | FS | Bonferroni lower | 894 | 0.347 | 0/894 = 0.0000 | 0/894 = 0.0000 | 0/2000 = 0.0000 [0.0000, 0.0019] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n500 | DINA | naive | 1116 | 0.724 | 97/1116 = 0.0869 | 17/1116 = 0.0152 | 97/2000 = 0.0485 [0.0399, 0.0588] | 17/2000 = 0.0085 [0.0053, 0.0136] |
| null0657_n500 | DINA | nullid nv_H (MR off, ref.) | 1116 | 0.719 | 93/1116 = 0.0833 | 16/1116 = 0.0143 | 93/2000 = 0.0465 [0.0381, 0.0566] | 16/2000 = 0.0080 [0.0049, 0.0130] |
| null0657_n500 | DINA | IJ two-term | 1116 | 0.408 | 1/1116 = 0.0009 | 0/1116 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n500 | DINA | field | 1116 | 0.450 | 7/1116 = 0.0063 | 4/1116 = 0.0036 | 7/2000 = 0.0035 [0.0017, 0.0072] | 4/2000 = 0.0020 [0.0008, 0.0051] |
| null0657_n500 | DINA | Bonferroni lower | 1116 | 0.399 | 4/1116 = 0.0036 | 1/1116 = 0.0009 | 4/2000 = 0.0020 [0.0008, 0.0051] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0657_n500 | GRF | naive | 1883 | 0.726 | 121/1883 = 0.0643 | 22/1883 = 0.0117 | 121/2000 = 0.0605 [0.0509, 0.0718] | 22/2000 = 0.0110 [0.0073, 0.0166] |
| null0657_n500 | GRF | nullid nv_H (MR off, ref.) | 1884 | 0.722 | 113/1884 = 0.0600 | 20/1884 = 0.0106 | 113/2000 = 0.0565 [0.0472, 0.0675] | 20/2000 = 0.0100 [0.0065, 0.0154] |
| null0657_n500 | GRF | IJ two-term | 1883 | 0.349 | 0/1883 = 0.0000 | 0/1883 = 0.0000 | 0/2000 = 0.0000 [0.0000, 0.0019] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n500 | GRF | field | 1883 | 0.351 | 11/1883 = 0.0058 | 3/1883 = 0.0016 | 11/2000 = 0.0055 [0.0031, 0.0098] | 3/2000 = 0.0015 [0.0005, 0.0044] |
| null0657_n500 | GRF | Bonferroni lower | 1883 | 0.304 | 3/1883 = 0.0016 | 0/1883 = 0.0000 | 3/2000 = 0.0015 [0.0005, 0.0044] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n500 | FS | naive | 1356 | 0.871 | 237/1356 = 0.1748 | 24/1356 = 0.0177 | 237/2000 = 0.1185 [0.1051, 0.1334] | 24/2000 = 0.0120 [0.0081, 0.0178] |
| null0721_n500 | FS | nullid nv_H (MR off, ref.) | 1356 | 0.865 | 225/1356 = 0.1659 | 21/1356 = 0.0155 | 225/2000 = 0.1125 [0.0994, 0.1271] | 21/2000 = 0.0105 [0.0069, 0.0160] |
| null0721_n500 | FS | IJ two-term | 1356 | 0.443 | 0/1356 = 0.0000 | 0/1356 = 0.0000 | 0/2000 = 0.0000 [0.0000, 0.0019] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n500 | FS | field | 1356 | 0.422 | 2/1356 = 0.0015 | 0/1356 = 0.0000 | 2/2000 = 0.0010 [0.0003, 0.0036] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n500 | FS | Bonferroni lower | 1356 | 0.366 | 0/1356 = 0.0000 | 0/1356 = 0.0000 | 0/2000 = 0.0000 [0.0000, 0.0019] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n500 | DINA | naive | 1480 | 0.815 | 278/1480 = 0.1878 | 53/1480 = 0.0358 | 278/2000 = 0.1390 [0.1245, 0.1549] | 53/2000 = 0.0265 [0.0203, 0.0345] |
| null0721_n500 | DINA | nullid nv_H (MR off, ref.) | 1480 | 0.809 | 269/1480 = 0.1818 | 46/1480 = 0.0311 | 269/2000 = 0.1345 [0.1202, 0.1502] | 46/2000 = 0.0230 [0.0173, 0.0305] |
| null0721_n500 | DINA | IJ two-term | 1480 | 0.438 | 4/1480 = 0.0027 | 0/1480 = 0.0000 | 4/2000 = 0.0020 [0.0008, 0.0051] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n500 | DINA | field | 1480 | 0.490 | 18/1480 = 0.0122 | 8/1480 = 0.0054 | 18/2000 = 0.0090 [0.0057, 0.0142] | 8/2000 = 0.0040 [0.0020, 0.0079] |
| null0721_n500 | DINA | Bonferroni lower | 1480 | 0.433 | 11/1480 = 0.0074 | 3/1480 = 0.0020 | 11/2000 = 0.0055 [0.0031, 0.0098] | 3/2000 = 0.0015 [0.0005, 0.0044] |
| null0721_n500 | GRF | naive | 1969 | 0.803 | 284/1969 = 0.1442 | 53/1969 = 0.0269 | 284/2000 = 0.1420 [0.1274, 0.1580] | 53/2000 = 0.0265 [0.0203, 0.0345] |
| null0721_n500 | GRF | nullid nv_H (MR off, ref.) | 1969 | 0.799 | 268/1969 = 0.1361 | 49/1969 = 0.0249 | 268/2000 = 0.1340 [0.1198, 0.1496] | 49/2000 = 0.0245 [0.0186, 0.0322] |
| null0721_n500 | GRF | IJ two-term | 1969 | 0.395 | 0/1969 = 0.0000 | 0/1969 = 0.0000 | 0/2000 = 0.0000 [0.0000, 0.0019] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n500 | GRF | field | 1969 | 0.398 | 21/1969 = 0.0107 | 5/1969 = 0.0025 | 21/2000 = 0.0105 [0.0069, 0.0160] | 5/2000 = 0.0025 [0.0011, 0.0058] |
| null0721_n500 | GRF | Bonferroni lower | 1969 | 0.346 | 8/1969 = 0.0041 | 1/1969 = 0.0005 | 8/2000 = 0.0040 [0.0020, 0.0079] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0657_n1000 | FS | naive | 637 | 0.832 | 28/637 = 0.0440 | 5/637 = 0.0078 | 28/2000 = 0.0140 [0.0097, 0.0202] | 5/2000 = 0.0025 [0.0011, 0.0058] |
| null0657_n1000 | FS | nullid nv_H (MR off, ref.) | 637 | 0.829 | 26/637 = 0.0408 | 4/637 = 0.0063 | 26/2000 = 0.0130 [0.0089, 0.0190] | 4/2000 = 0.0020 [0.0008, 0.0051] |
| null0657_n1000 | FS | IJ two-term | 637 | 0.521 | 1/637 = 0.0016 | 0/637 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1000 | FS | field | 637 | 0.506 | 3/637 = 0.0047 | 1/637 = 0.0016 | 3/2000 = 0.0015 [0.0005, 0.0044] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0657_n1000 | FS | Bonferroni lower | 637 | 0.456 | 1/637 = 0.0016 | 0/637 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1000 | DINA | naive | 593 | 0.705 | 15/593 = 0.0253 | 5/593 = 0.0084 | 15/2000 = 0.0075 [0.0046, 0.0123] | 5/2000 = 0.0025 [0.0011, 0.0058] |
| null0657_n1000 | DINA | nullid nv_H (MR off, ref.) | 593 | 0.704 | 12/593 = 0.0202 | 4/593 = 0.0067 | 12/2000 = 0.0060 [0.0034, 0.0105] | 4/2000 = 0.0020 [0.0008, 0.0051] |
| null0657_n1000 | DINA | IJ two-term | 593 | 0.472 | 0/593 = 0.0000 | 0/593 = 0.0000 | 0/2000 = 0.0000 [0.0000, 0.0019] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1000 | DINA | field | 593 | 0.520 | 5/593 = 0.0084 | 1/593 = 0.0017 | 5/2000 = 0.0025 [0.0011, 0.0058] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0657_n1000 | DINA | Bonferroni lower | 593 | 0.475 | 4/593 = 0.0067 | 0/593 = 0.0000 | 4/2000 = 0.0020 [0.0008, 0.0051] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1000 | GRF | naive | 1663 | 0.709 | 40/1663 = 0.0241 | 7/1663 = 0.0042 | 40/2000 = 0.0200 [0.0147, 0.0271] | 7/2000 = 0.0035 [0.0017, 0.0072] |
| null0657_n1000 | GRF | nullid nv_H (MR off, ref.) | 1662 | 0.707 | 38/1662 = 0.0229 | 6/1662 = 0.0036 | 38/2000 = 0.0190 [0.0139, 0.0260] | 6/2000 = 0.0030 [0.0014, 0.0065] |
| null0657_n1000 | GRF | IJ two-term | 1663 | 0.418 | 1/1663 = 0.0006 | 0/1663 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1000 | GRF | field | 1663 | 0.411 | 5/1663 = 0.0030 | 1/1663 = 0.0006 | 5/2000 = 0.0025 [0.0011, 0.0058] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0657_n1000 | GRF | Bonferroni lower | 1663 | 0.369 | 2/1663 = 0.0012 | 0/1663 = 0.0000 | 2/2000 = 0.0010 [0.0003, 0.0036] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n1000 | FS | naive | 1234 | 0.841 | 90/1234 = 0.0729 | 9/1234 = 0.0073 | 90/2000 = 0.0450 [0.0368, 0.0550] | 9/2000 = 0.0045 [0.0024, 0.0085] |
| null0721_n1000 | FS | nullid nv_H (MR off, ref.) | 1234 | 0.838 | 85/1234 = 0.0689 | 9/1234 = 0.0073 | 85/2000 = 0.0425 [0.0345, 0.0523] | 9/2000 = 0.0045 [0.0024, 0.0085] |
| null0721_n1000 | FS | IJ two-term | 1234 | 0.523 | 1/1234 = 0.0008 | 0/1234 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n1000 | FS | field | 1234 | 0.514 | 6/1234 = 0.0049 | 0/1234 = 0.0000 | 6/2000 = 0.0030 [0.0014, 0.0065] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n1000 | FS | Bonferroni lower | 1234 | 0.463 | 2/1234 = 0.0016 | 0/1234 = 0.0000 | 2/2000 = 0.0010 [0.0003, 0.0036] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n1000 | DINA | naive | 1152 | 0.754 | 76/1152 = 0.0660 | 7/1152 = 0.0061 | 76/2000 = 0.0380 [0.0305, 0.0473] | 7/2000 = 0.0035 [0.0017, 0.0072] |
| null0721_n1000 | DINA | nullid nv_H (MR off, ref.) | 1152 | 0.753 | 75/1152 = 0.0651 | 7/1152 = 0.0061 | 75/2000 = 0.0375 [0.0300, 0.0468] | 7/2000 = 0.0035 [0.0017, 0.0072] |
| null0721_n1000 | DINA | IJ two-term | 1152 | 0.501 | 1/1152 = 0.0009 | 0/1152 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n1000 | DINA | field | 1152 | 0.540 | 10/1152 = 0.0087 | 2/1152 = 0.0017 | 10/2000 = 0.0050 [0.0027, 0.0092] | 2/2000 = 0.0010 [0.0003, 0.0036] |
| null0721_n1000 | DINA | Bonferroni lower | 1152 | 0.493 | 6/1152 = 0.0052 | 0/1152 = 0.0000 | 6/2000 = 0.0030 [0.0014, 0.0065] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n1000 | GRF | naive | 1886 | 0.767 | 101/1886 = 0.0536 | 14/1886 = 0.0074 | 101/2000 = 0.0505 [0.0417, 0.0610] | 14/2000 = 0.0070 [0.0042, 0.0117] |
| null0721_n1000 | GRF | nullid nv_H (MR off, ref.) | 1886 | 0.765 | 97/1886 = 0.0514 | 11/1886 = 0.0058 | 97/2000 = 0.0485 [0.0399, 0.0588] | 11/2000 = 0.0055 [0.0031, 0.0098] |
| null0721_n1000 | GRF | IJ two-term | 1886 | 0.455 | 2/1886 = 0.0011 | 0/1886 = 0.0000 | 2/2000 = 0.0010 [0.0003, 0.0036] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n1000 | GRF | field | 1886 | 0.446 | 10/1886 = 0.0053 | 2/1886 = 0.0011 | 10/2000 = 0.0050 [0.0027, 0.0092] | 2/2000 = 0.0010 [0.0003, 0.0036] |
| null0721_n1000 | GRF | Bonferroni lower | 1886 | 0.400 | 4/1886 = 0.0021 | 1/1886 = 0.0005 | 4/2000 = 0.0020 [0.0008, 0.0051] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0657_n1500 | FS | naive | 337 | 0.818 | 6/337 = 0.0178 | 0/337 = 0.0000 | 6/2000 = 0.0030 [0.0014, 0.0065] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | FS | nullid nv_H (MR off, ref.) | 337 | 0.816 | 6/337 = 0.0178 | 0/337 = 0.0000 | 6/2000 = 0.0030 [0.0014, 0.0065] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | FS | IJ two-term | 337 | 0.577 | 0/337 = 0.0000 | 0/337 = 0.0000 | 0/2000 = 0.0000 [0.0000, 0.0019] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | FS | field | 337 | 0.562 | 1/337 = 0.0030 | 0/337 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | FS | Bonferroni lower | 337 | 0.518 | 0/337 = 0.0000 | 0/337 = 0.0000 | 0/2000 = 0.0000 [0.0000, 0.0019] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | DINA | naive | 278 | 0.697 | 1/278 = 0.0036 | 0/278 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | DINA | nullid nv_H (MR off, ref.) | 278 | 0.694 | 1/278 = 0.0036 | 0/278 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | DINA | IJ two-term | 278 | 0.507 | 0/278 = 0.0000 | 0/278 = 0.0000 | 0/2000 = 0.0000 [0.0000, 0.0019] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | DINA | field | 278 | 0.541 | 0/278 = 0.0000 | 0/278 = 0.0000 | 0/2000 = 0.0000 [0.0000, 0.0019] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | DINA | Bonferroni lower | 278 | 0.498 | 0/278 = 0.0000 | 0/278 = 0.0000 | 0/2000 = 0.0000 [0.0000, 0.0019] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | GRF | naive | 1129 | 0.703 | 6/1129 = 0.0053 | 1/1129 = 0.0009 | 6/2000 = 0.0030 [0.0014, 0.0065] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0657_n1500 | GRF | nullid nv_H (MR off, ref.) | 1129 | 0.701 | 6/1129 = 0.0053 | 1/1129 = 0.0009 | 6/2000 = 0.0030 [0.0014, 0.0065] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0657_n1500 | GRF | IJ two-term | 1129 | 0.471 | 0/1129 = 0.0000 | 0/1129 = 0.0000 | 0/2000 = 0.0000 [0.0000, 0.0019] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | GRF | field | 1129 | 0.475 | 3/1129 = 0.0027 | 0/1129 = 0.0000 | 3/2000 = 0.0015 [0.0005, 0.0044] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | GRF | Bonferroni lower | 1129 | 0.435 | 1/1129 = 0.0009 | 0/1129 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n1500 | FS | naive | 969 | 0.823 | 19/969 = 0.0196 | 1/969 = 0.0010 | 19/2000 = 0.0095 [0.0061, 0.0148] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0721_n1500 | FS | nullid nv_H (MR off, ref.) | 969 | 0.822 | 17/969 = 0.0175 | 1/969 = 0.0010 | 17/2000 = 0.0085 [0.0053, 0.0136] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0721_n1500 | FS | IJ two-term | 969 | 0.576 | 1/969 = 0.0010 | 0/969 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n1500 | FS | field | 969 | 0.573 | 1/969 = 0.0010 | 0/969 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n1500 | FS | Bonferroni lower | 969 | 0.527 | 1/969 = 0.0010 | 0/969 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n1500 | DINA | naive | 748 | 0.735 | 14/748 = 0.0187 | 1/748 = 0.0013 | 14/2000 = 0.0070 [0.0042, 0.0117] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0721_n1500 | DINA | nullid nv_H (MR off, ref.) | 748 | 0.733 | 14/748 = 0.0187 | 1/748 = 0.0013 | 14/2000 = 0.0070 [0.0042, 0.0117] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0721_n1500 | DINA | IJ two-term | 748 | 0.541 | 0/748 = 0.0000 | 0/748 = 0.0000 | 0/2000 = 0.0000 [0.0000, 0.0019] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n1500 | DINA | field | 748 | 0.581 | 3/748 = 0.0040 | 1/748 = 0.0013 | 3/2000 = 0.0015 [0.0005, 0.0044] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0721_n1500 | DINA | Bonferroni lower | 748 | 0.538 | 1/748 = 0.0013 | 0/748 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n1500 | GRF | naive | 1690 | 0.746 | 29/1690 = 0.0172 | 2/1690 = 0.0012 | 29/2000 = 0.0145 [0.0101, 0.0207] | 2/2000 = 0.0010 [0.0003, 0.0036] |
| null0721_n1500 | GRF | nullid nv_H (MR off, ref.) | 1689 | 0.746 | 29/1689 = 0.0172 | 2/1689 = 0.0012 | 29/2000 = 0.0145 [0.0101, 0.0207] | 2/2000 = 0.0010 [0.0003, 0.0036] |
| null0721_n1500 | GRF | IJ two-term | 1690 | 0.503 | 2/1690 = 0.0012 | 0/1690 = 0.0000 | 2/2000 = 0.0010 [0.0003, 0.0036] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n1500 | GRF | field | 1690 | 0.499 | 2/1690 = 0.0012 | 1/1690 = 0.0006 | 2/2000 = 0.0010 [0.0003, 0.0036] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0721_n1500 | GRF | Bonferroni lower | 1690 | 0.458 | 2/1690 = 0.0012 | 0/1690 = 0.0000 | 2/2000 = 0.0010 [0.0003, 0.0036] | 0/2000 = 0.0000 [0.0000, 0.0019] |

- **The adjusted bounds on a false region almost never reach HR 1.00.** The field lower bound does so
  on 0–21 of 2,000 replicates per cell (unconditional 0.0000–0.0105). The Bonferroni lower member
  does so on 0–11 (0.0000–0.0055), and the IJ bound on 0–4 (0.0000–0.0020). At HR 1.25 the counts
  are 0–8 (field), 0–3 (Bonferroni) and 0 (IJ).
- **The naive bound reaches 1.00 far more often:** unconditionally 0.0005–0.1420, highest at n 500
  and at HR 0.721 (FS 0.1185, DINA 0.1390, GRF 0.1420 at `null0721_n500`).
- **Median adjusted lower bounds sit far below 1.00.** Field medians are 0.35–0.58 on the HR scale,
  against naive medians of 0.70–0.87.
- **The naive row runs slightly above `nullid`'s reference row everywhere** (e.g. 101 against 93
  claims at FS `null0657_n500`). The point estimates are identical (identity gate). MR's naive SE
  (`se_wald`) is somewhat smaller than the MR-off `coxph` model SE.

### Table 3 — context per cell: the realized trial against the target, the truths in Ĥ and Ĥᶜ, and mr_ok

`itt_est` is the whole-trial Cox HR, which Gate C shows is engine-independent. The truths β(Ĥ) and
β(Ĥᶜ) are super-population marginal Cox HRs on uncensored potential outcomes; the trial fits censored
data (`nullc125` §5). Medians are over declaring replicates.

| cell | identifier | target HR | itt_est median [Q1, Q3] | declared | median β(Ĥ) | median β(Ĥᶜ) | median β(Ĥᶜ) − target | mr_ok among declaring |
|---|---|---|---|---|---|---|---|---|
| null0657_n500 | FS | 0.657 | 0.6287 [0.5713, 0.6877] | 894 | 0.6162 | 0.6284 | -0.0286 | 894/894 = 1.0000 [0.9957, 1.0000] |
| null0657_n500 | DINA | 0.657 | 0.6287 [0.5713, 0.6877] | 1116 | 0.6179 | 0.6276 | -0.0294 | 1116/1116 = 1.0000 [0.9966, 1.0000] |
| null0657_n500 | GRF | 0.657 | 0.6287 [0.5713, 0.6877] | 1883 | 0.6125 | 0.6287 | -0.0283 | 1883/1883 = 1.0000 [0.9980, 1.0000] |
| null0721_n500 | FS | 0.721 | 0.6987 [0.6365, 0.7602] | 1356 | 0.6867 | 0.6982 | -0.0228 | 1356/1356 = 1.0000 [0.9972, 1.0000] |
| null0721_n500 | DINA | 0.721 | 0.6987 [0.6365, 0.7602] | 1480 | 0.6881 | 0.6977 | -0.0233 | 1480/1480 = 1.0000 [0.9974, 1.0000] |
| null0721_n500 | GRF | 0.721 | 0.6987 [0.6365, 0.7602] | 1969 | 0.6844 | 0.6983 | -0.0227 | 1969/1969 = 1.0000 [0.9981, 1.0000] |
| null0657_n1000 | FS | 0.657 | 0.6282 [0.5922, 0.6687] | 637 | 0.6158 | 0.6278 | -0.0292 | 637/637 = 1.0000 [0.9940, 1.0000] |
| null0657_n1000 | DINA | 0.657 | 0.6282 [0.5922, 0.6687] | 593 | 0.6150 | 0.6276 | -0.0294 | 593/593 = 1.0000 [0.9936, 1.0000] |
| null0657_n1000 | GRF | 0.657 | 0.6282 [0.5922, 0.6687] | 1663 | 0.6133 | 0.6285 | -0.0285 | 1663/1663 = 1.0000 [0.9977, 1.0000] |
| null0721_n1000 | FS | 0.721 | 0.6966 [0.6566, 0.7387] | 1234 | 0.6883 | 0.6978 | -0.0232 | 1234/1234 = 1.0000 [0.9969, 1.0000] |
| null0721_n1000 | DINA | 0.721 | 0.6966 [0.6566, 0.7387] | 1152 | 0.6880 | 0.6976 | -0.0234 | 1152/1152 = 1.0000 [0.9967, 1.0000] |
| null0721_n1000 | GRF | 0.721 | 0.6966 [0.6566, 0.7387] | 1886 | 0.6848 | 0.6983 | -0.0227 | 1886/1886 = 1.0000 [0.9980, 1.0000] |
| null0657_n1500 | FS | 0.657 | 0.6266 [0.5967, 0.6593] | 337 | 0.6162 | 0.6277 | -0.0293 | 337/337 = 1.0000 [0.9887, 1.0000] |
| null0657_n1500 | DINA | 0.657 | 0.6266 [0.5967, 0.6593] | 278 | 0.6121 | 0.6274 | -0.0296 | 278/278 = 1.0000 [0.9864, 1.0000] |
| null0657_n1500 | GRF | 0.657 | 0.6266 [0.5967, 0.6593] | 1129 | 0.6127 | 0.6284 | -0.0286 | 1129/1129 = 1.0000 [0.9966, 1.0000] |
| null0721_n1500 | FS | 0.721 | 0.6960 [0.6611, 0.7305] | 969 | 0.6877 | 0.6978 | -0.0232 | 969/969 = 1.0000 [0.9961, 1.0000] |
| null0721_n1500 | DINA | 0.721 | 0.6960 [0.6611, 0.7305] | 748 | 0.6867 | 0.6970 | -0.0240 | 748/748 = 1.0000 [0.9949, 1.0000] |
| null0721_n1500 | GRF | 0.721 | 0.6960 [0.6611, 0.7305] | 1690 | 0.6866 | 0.6983 | -0.0227 | 1690/1690 = 1.0000 [0.9977, 1.0000] |

- **MR returned a result on every declaring replicate**, 21,014 of 21,014.
- **Both truths sit below the design target.** The trial's ITT HR has median 0.627–0.629 (target
  0.657) and 0.696–0.699 (target 0.721). The truths behave the same way: β(Ĥᶜ) medians are
  0.627–0.629 / 0.697–0.698, 0.023–0.030 below target, and β(Ĥ) medians 0.612–0.618 / 0.684–0.688.
  β(·) is evaluated on the super-population under the trials' censored analysis (template :856-861),
  so it tracks the censored trial HR rather than the uncensored target. Coverage is scored against
  β(·) and is unaffected by this gap; it is recorded as `nullc125` §5 recorded it.
- **Every Ĥ carries a benefit-side truth.** In a false "harm" region the target of the lower bound
  lies at HR about 0.61–0.69, the uniform effect.

### Table 4 — identification against `nullid`, recorded, not gated

Each replicate is compared by `sim_id` with the committed `nullid` bundle (Mac Studio, R 4.5.2,
12 workers, MR off). Label, |Ĥ| and rule are compared where both runs declared.

| cell | identifier | declared nullmr / nullid | declaration differs | label differs (both declared) | \|Ĥ\| differs (both declared) | sg_def differs (both declared) |
|---|---|---|---|---|---|---|
| null0657_n500 | FS | 894 / 894 | 0 | 0 / 894 | 0 / 894 | 0 / 894 |
| null0657_n500 | DINA | 1116 / 1116 | 0 | 0 / 1116 | 0 / 1116 | 0 / 1116 |
| null0657_n500 | GRF | 1883 / 1884 | 1 | 0 / 1883 | 0 / 1883 | 0 / 1883 |
| null0721_n500 | FS | 1356 / 1356 | 0 | 0 / 1356 | 0 / 1356 | 0 / 1356 |
| null0721_n500 | DINA | 1480 / 1480 | 0 | 0 / 1480 | 0 / 1480 | 0 / 1480 |
| null0721_n500 | GRF | 1969 / 1969 | 0 | 0 / 1969 | 0 / 1969 | 0 / 1969 |
| null0657_n1000 | FS | 637 / 637 | 0 | 0 / 637 | 0 / 637 | 0 / 637 |
| null0657_n1000 | DINA | 593 / 593 | 0 | 0 / 593 | 0 / 593 | 0 / 593 |
| null0657_n1000 | GRF | 1663 / 1662 | 1 | 0 / 1662 | 0 / 1662 | 0 / 1662 |
| null0721_n1000 | FS | 1234 / 1234 | 0 | 0 / 1234 | 0 / 1234 | 0 / 1234 |
| null0721_n1000 | DINA | 1152 / 1152 | 0 | 0 / 1152 | 0 / 1152 | 0 / 1152 |
| null0721_n1000 | GRF | 1886 / 1886 | 0 | 0 / 1886 | 0 / 1886 | 0 / 1886 |
| null0657_n1500 | FS | 337 / 337 | 0 | 0 / 337 | 0 / 337 | 0 / 337 |
| null0657_n1500 | DINA | 278 / 278 | 0 | 0 / 278 | 0 / 278 | 0 / 278 |
| null0657_n1500 | GRF | 1129 / 1129 | 0 | 0 / 1129 | 0 / 1129 | 0 / 1129 |
| null0721_n1500 | FS | 969 / 969 | 0 | 0 / 969 | 0 / 969 | 0 / 969 |
| null0721_n1500 | DINA | 748 / 748 | 0 | 0 / 748 | 0 / 748 | 0 / 748 |
| null0721_n1500 | GRF | 1690 / 1689 | 1 | 0 / 1689 | 0 / 1689 | 0 / 1689 |

- **FS and DINA reproduce `nullid` exactly across machines** — every declaration, label, size and
  rule in all six cells.
- **GRF differs on 3 of 36,000 replicate-runs, each a declaration flip:**
  - `null0657_n500` `sim_id` 1022: `nullid` declared `size > 19 & grade > 2` (`admitted_n` 4);
    nullmr declared nothing.
  - `null0657_n1000` `sim_id` 216: nullmr declared `age > 62 & size > 20` (`admitted_n` 4); `nullid`
    did not.
  - `null0721_n1500` `sim_id` 1405: nullmr declared `age > 44 & size <= 24` (`admitted_n` 5);
    `nullid` did not.
  On the non-declaring side `admitted_n` is NA, not 0. That points to the GRF step returning no
  admitted count on that machine, not to an empty admitted set. MR cannot change identification (the
  identity gate shows it on the same machine), so these are cross-machine differences (Mac R 4.5.2
  vs pop-os R 4.6.1). This is **not verified per replicate**: no same-machine MR-off render of those
  `sim_id`s exists, and making one would be a re-run, out of scope.

### Table 5 — wall per cell, and the machine, R version and build

| cell | wall s FS / DINA / GRF | cell s |
|---|---|---|
| null0657_n500 | 900 / 598 / 1372 | 2870 |
| null0721_n500 | 1219 / 910 / 1442 | 3571 |
| null0657_n1000 | 1002 / 337 / 1501 | 2840 |
| null0721_n1000 | 1559 / 628 / 1645 | 3832 |
| null0657_n1500 | 796 / 233 / 1431 | 2460 |
| null0721_n1500 | 1925 / 442 / 1918 | 4285 |
| **total (18 renders)** | | **19858 s = 5.516 h** |

| hostnames | R | forestsearch | workers | bundles |
|---|---|---|---|---|
| pop-os | 4.6.1 | 0.3.5.9000 (built 2026-09-19 03:54:02 UTC) | 64 | 18 |

- **Wall tracks the declaration count**, because MR runs only where something was declared. The n
  1500 cell at 0.657, with the fewest declarations, is the cheapest cell, not the dearest.

---

## 5. What is not reported, and why

- **The oracle** is undefined with no planted region: `or_H_*` is a fit on an empty set. Its
  columns are in the bundles and are not read.
- **Winner-only / winner-floor variants** (`mr_*_se_w`, `_wf` and their bounds) are recorded by the
  template whatever the knob, but were rejected on theoretical grounds (template :407-410) and are
  not evaluated.
- **The full bootstrap** was not run (`FS_S7_FB=none`); Gate A asserts every `fb_*` product NA.
- **Sensitivity and PPV are undefined** with an empty planted region: sensitivity is NA everywhere
  (Gate A) and PPV is 0 by construction.
- **The identified-to-planted size ratio** has no denominator.
- The unscaled complement field (`fld_Hc_up1s`) and the unscaled Bonferroni pair (`joint_bonf`) are
  computed by `covs_cell` but are not the products this task evaluates. They are in
  `scripts_dinamr/nullmr_findings.rds`.

## 6. Findings

- **In a false region, the selection-adjusted lower bounds almost never support a harm claim.** The
  field lower bound on β(Ĥ) reaches HR 1.00 on at most 21 of 2,000 replicates in any cell
  (unconditional 0.0000–0.0105, highest GRF `null0721_n500`). The Bonferroni lower member reaches it
  on at most 11 (0.0055), the IJ bound on at most 4 (0.0020). At HR 1.25: at most 8, 3 and 0. Median
  field bounds are 0.35–0.58. **The naive bound reaches 1.00 on up to 14.2% of all replicates**
  (GRF `null0721_n500`; FS 11.9%, DINA 13.9% there). The ratio of naive to field counts runs from about 2× (GRF `null0657_n1500`, 6 against 3) to over
  100× (FS `null0721_n500`, 237 against 2), and is 10× or more in 10 of the 17 cell-identifiers where
  the field count is non-zero. That gap is what MR buys against false claims.
- **Coverage of the false region's own target is not uniform.**
  - **GRF**: field lower 0.939–0.957, field-s upper 0.922–0.961, Bonferroni pair 0.934–0.955 — near
    nominal in all six cells.
  - **FS**: field lower 0.934–0.974 in four cells, but **0.879 [0.852, 0.902] at `null0657_n1000`
    and 0.792 [0.746, 0.832] at `null0657_n1500`**. The two cells below nominal are FS's two lowest
    declaration rates (0.1685 → 0.792, 0.3185 → 0.879), and the highest rate gives the highest
    coverage (0.678 → 0.974). The order is not exact in the middle (0.447 → 0.956 against
    0.4845 → 0.934). A reading consistent with this: the rarer the declaration, the more extreme the
    selection it conditions on, and the more the field under-corrects. It is an observed association,
    not a tested mechanism.
  - **DINA**: field lower **0.763–0.879, below nominal in every cell** (conditional on the proposed
    family). It rises with the declaration rate (0.139 → 0.763 … 0.576 → 0.879) except for the top
    pair (0.740 → 0.870). The Bonferroni pair
    inherits it (0.827–0.902).
- **The undercoverage is on the safe side for a harm claim.** A field lower bound that misses β(Ĥ)
  misses it by sitting above the truth, and the truth here is a benefit (HR about 0.61–0.69). Yet
  Table 2 shows those bounds still sit below 1.00 on at least 98.8% of declaring replicates in every cell. Undercoverage of
  β(Ĥ) here is not a rate of false harm claims.
- **The complement's upper bound behaves on every engine** (field-s 0.92–0.97). A benefit claim on
  the complement ("HR in Ĥᶜ at most U") is close to nominal at the null.
- **IJ two-term (two-sided) is conservative here**: 0.929–0.995 on β(Ĥ), 0.997–1.000 on β(Ĥᶜ).
- **Naive inference fails as expected.** The one-sided naive lower bound covers β(Ĥ) on 0 of 5,427
  FS declarations, and on 8–30% (DINA) and 2–18% (GRF).
- **MR changes nothing about identification.** On the same machine, 46 of 46 identification columns
  and the truth object are identical to `nullc125`'s MR-off render (identity gate). Across machines,
  FS and DINA reproduce `nullid` exactly on all 12,000 replicate-runs each. GRF differs on 3 of
  12,000, which cannot be attributed per replicate within this scope.
- **MR returned a result on every one of the 21,014 declarations**; every evaluated product was finite
  wherever it did (Gate A).
- **The truths sit 0.02–0.03 below the design target**, consistent with β(·) being evaluated under the
  censored trial analysis. Coverage is scored against β(·) itself, so this gap is context, not a bias
  in the coverage.

## 7. Files

| what | path |
|---|---|
| bundles (18) | `results/{fs,dina,grf}_effMaxSG_fb_mr_field_m1_h0{66,72}_knoise0_n{500,1000,1500}_null{657,721}_nb20_nullmr_res_1_2000.rds` |
| smoke bundles (6) | `results/*_nullmrsmoke_quickrun_res_1_20.rds` |
| driver, cells, gates, identity gate, projection, findings | `scripts_dinamr/nullmr.sh`, `nullmr.cells`, `nullmr_smoke.cells`, `nullmr_gateA.R`, `nullmr_gateC.R`, `nullmr_identity.R`, `nullmr_project.R`, `nullmr_findings.R` (+ `nullmr_findings.rds`) |
| Step 1 record | `scripts_dinamr/logs/nullmr_step1_record.md` |
| logs | `scripts_dinamr/logs/nullmr*` (driver, per-render, Gate A, Gate C, identity, smoke/projection, findings) |
| renders (24; committed at closeout, following `nullc125`) | `nullmr_*.html` (18), `nullmrsmoke_*.html` (6) in the study directory |
| template diff | none — the template was not changed |
