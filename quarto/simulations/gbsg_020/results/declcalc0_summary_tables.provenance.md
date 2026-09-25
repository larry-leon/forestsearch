# Provenance — `declcalc0_summary_tables.csv`

Produced by `dev/analysis/declcal_c0_summary/declcal_c0_summary.R` (run from the repository root, no arguments) for `TASK_declcal_c0_summary_2026-09-25` (`dev/tasks/`). The script only reads the committed payloads, and every value comes from a payload column or `meta` field. The only constant is the fixed cutoff `k_fix = 2.0`, which the task specifies. Column semantics follow `scripts_dinamr/declcalc0_run.R` at `a46bf9b7` and `dev/reports/REPORT_declcal_c0_campaign_reads_2026-09-24.md`, Q6.

## Columns of the CSV

`table`, `row`, `col` label the quantity. `alpha` is 0.05 or 0.10, or blank where it does not apply. `c0` is the protected level (`c070`, `c075`, `c080`, `c085`), or blank for the comparators. `cell` is the cell, or the `;`-joined cell set a summary ranges over. `design_hr` and `n` come from each payload's `meta$target_hr` and `meta$n`, not from its label. `value` is at full precision. `count` is the number of ones, for rates only. `denominator` is the number of replicates (rates), non-missing replicates (pooled medians) or cells (medians of per-cell medians). `source_column` is the payload column or expression, and `source_payload` the file(s) read.

## Cell sets

The cells are read from `meta`. The `declcalc0` payloads carry the protected-level rule and the `declcal` payloads carry the claim-threshold rule. Both sets share cell ids, designs, `n` and seeds.

| payload | cell | dgm | design HR | n |
|---|---|---|---|---|
| `declcalc0_inull_B1_res_1_2000.rds` | B1 | null0657 | 0.657 | 500 |
| `declcalc0_inull_B2_res_1_2000.rds` | B2 | null0657 | 0.657 | 1000 |
| `declcalc0_inull_B3_res_1_2000.rds` | B3 | null0657 | 0.657 | 1500 |
| `declcalc0_inull_B4_res_1_2000.rds` | B4 | null0721 | 0.721 | 500 |
| `declcalc0_inull_B5_res_1_2000.rds` | B5 | null0721 | 0.721 | 1000 |
| `declcalc0_inull_B6_res_1_2000.rds` | B6 | null0721 | 0.721 | 1500 |
| `declcalc0_power_C1_res_1_2000.rds` | C1 | alt_h150 | 1.5 | 1000 |
| `declcalc0_power_C2_res_1_2000.rds` | C2 | alt_h150 | 1.5 | 1500 |
| `declcalc0_power_C3_res_1_2000.rds` | C3 | alt_h200 | 2 | 1000 |
| `declcalc0_power_C4_res_1_2000.rds` | C4 | alt_h200 | 2 | 1500 |
| `declcal_inull_B1_res_1_2000.rds` | B1 | null0657 | 0.657 | 500 |
| `declcal_inull_B2_res_1_2000.rds` | B2 | null0657 | 0.657 | 1000 |
| `declcal_inull_B3_res_1_2000.rds` | B3 | null0657 | 0.657 | 1500 |
| `declcal_inull_B4_res_1_2000.rds` | B4 | null0721 | 0.721 | 500 |
| `declcal_inull_B5_res_1_2000.rds` | B5 | null0721 | 0.721 | 1000 |
| `declcal_inull_B6_res_1_2000.rds` | B6 | null0721 | 0.721 | 1500 |
| `declcal_power_C1_res_1_2000.rds` | C1 | alt_h150 | 1.5 | 1000 |
| `declcal_power_C2_res_1_2000.rds` | C2 | alt_h150 | 1.5 | 1500 |
| `declcal_power_C3_res_1_2000.rds` | C3 | alt_h200 | 2 | 1000 |
| `declcal_power_C4_res_1_2000.rds` | C4 | alt_h200 | 2 | 1500 |

- **B (uniform benefit):** B1–B6. **C (planted harm):** C1–C4.
- The `declcal` Block A payloads (`declcal_bnull_A*`) are checksummed below but not read, because no S1.8 table uses them.

## Definitions

Every rate is the mean of a 0/1 indicator over a cell's 2,000 replicates. A missing value (`NA`) counts as not declared.

- **Table S3 (`row = conventional`).** The rate of `declared_conv` (the conventional screen as executed: rounded rate, post-reduction family) in each B cell, from `declcalc0`.
- **Tables S4 (α 0.10) and S5 (α 0.05).**
  - `calibrated_<c0>`: the rate of `declared_cal<α>_<c0>` from `declcalc0`, in every B cell and every C cell. `col = max_B` is the maximum over B1–B6, and `cell` names every B cell attaining it (ties joined by `;`).
  - `conventional`: `declared_conv`, as `max_B` and in each C cell.
  - `fixed_k2.0`: `max_T_post >= 2.0` (`NA` counts as not declared), in each B cell, as `max_B` and in each C cell.
  - `claim_threshold`: the unshifted `declared_cal05` / `declared_cal10` from the `declcal` payloads, over the same B and C cells.
  - `implied_pstar` rows:
    - conventional: `meta$p_star`, checked to be identical in all ten cells;
    - fixed: `2*pnorm(2.0)-1`;
    - claim threshold, over every cell and every n: `2*pnorm(k)-1`, where `k` is the pooled median of `kappa_hat_<α>` (`implied_pstar_pooled_median_allcells_alln`) or the median of the ten per-cell medians (`implied_pstar_median_of_cell_medians_allcells_alln`).
  - The comparator rows appear under both S4 and S5. Only the claim threshold's values depend on α.
- **Table S6 (`kappa_*`, `level_*`).** At each n (500, 1000, 1500), for each α and `<c0>` (column `kappa_hat_<α>_<c0>`, `declcalc0`) and for the claim threshold (`kappa_hat_<α>`, `declcal`), four candidate summaries:
  - `pooled_median_allcells`: the median of the per-replicate `kappa_hat` pooled over every cell at that n (B and C);
  - `median_of_cell_medians_allcells`: the median of the per-cell medians over those cells;
  - `pooled_median_Bcells` and `median_of_cell_medians_Bcells`: the same two, over the B cells only.
  - `level_<def>` is `2*pnorm(kappa)-1` of the corresponding summary.
  - At n 500 the all-cells set equals the B set (B1, B4), because there is no C cell at n 500.
  - The definition that reproduces the typed table is `median_of_cell_medians_allcells` (see the report).
- **Footprint (item d).** Over C1–C4, among replicates with `declared_cal05_<c0> == 1`: the share with `n_admitted_cal05_<c0> > 1`, per `<c0>` and pooled over the four `<c0>` (the counts are summed).

## Payload checksums (SHA-256)

These were taken before and after the run and are identical. All 23 payloads are tracked and show no `git diff`.

| payload | sha256 |
|---|---|
| `declcal_bnull_A1_res_1_2000.rds` | `9190508210880832787ec090f0579f9430071bfbaa2363636c65ac18cf6c2cb3` |
| `declcal_bnull_A2_res_1_2000.rds` | `8b787b76e9289466c32ccce08f1262064ca2e430f1c18bd652ee3107dba0157f` |
| `declcal_bnull_A3_res_1_2000.rds` | `30c0b946f3c650da2aabfdd75527753f3a45c813042ec898430ec91f9fbd3699` |
| `declcal_inull_B1_res_1_2000.rds` | `e3e0067c11a234aca43454127d3b345b6b54272fb47cb09684f538e65117d981` |
| `declcal_inull_B2_res_1_2000.rds` | `e99a1e9d5774851e8a5d86ebd495e947ec729ebc1e5b3ebc9412d55ee1a7e00f` |
| `declcal_inull_B3_res_1_2000.rds` | `a1faf9f422de2b1347d9820253444fcac7a7cb60ac2e0d31291f07ff0a5dc6fb` |
| `declcal_inull_B4_res_1_2000.rds` | `2fb4f731d8a8e8706bf51cf7bfca3ab10bcab09a3a1529b2afd82cc937a78d6f` |
| `declcal_inull_B5_res_1_2000.rds` | `354bb03a5397ea8eb7e350b8a9c2809aba0f32f53240799ddf7595ec96bbde87` |
| `declcal_inull_B6_res_1_2000.rds` | `8d2dcd85f20b978773e6dc71b21f73084a7e5228edbb1f0ca0bb0007e3a1af7f` |
| `declcal_power_C1_res_1_2000.rds` | `6a139f61bdbb38cf0cd895feee296080bcac65777df70793fe331b193b9da35b` |
| `declcal_power_C2_res_1_2000.rds` | `8d87bee590ae784b114ff0f95a700cf3404db8d2ee2a83ee53e8a04209e07dc4` |
| `declcal_power_C3_res_1_2000.rds` | `774b40450da5c15a183f436a62d6fdbf04b88a22bf32fa3a1375972af7cc48c6` |
| `declcal_power_C4_res_1_2000.rds` | `29e0bd24ccf88cbe9c967c9acb4f515264c2bf6c4876483cc7506200016449ae` |
| `declcalc0_inull_B1_res_1_2000.rds` | `eca8347440ec871130a66bb02c48e433a436d9961d12a4c69b15a980b4975048` |
| `declcalc0_inull_B2_res_1_2000.rds` | `ed99001f40f0812207e6133ec560761ce5760268edda89851ec1dce88fe636a9` |
| `declcalc0_inull_B3_res_1_2000.rds` | `f018e646972c57f30d6e4d7f7e621b6e13c11ed87b196ec7ddbffed1a8bdd4b3` |
| `declcalc0_inull_B4_res_1_2000.rds` | `aa0ecdb962f5ccf78f2cc36e817179434add9d525b374e3bde4af11c598b9f78` |
| `declcalc0_inull_B5_res_1_2000.rds` | `983b7ab87b520a96963527b8f3ef995511f83e1cea812fd2171a22b7220f5e1a` |
| `declcalc0_inull_B6_res_1_2000.rds` | `7162308da0eef522d2ca20ec85fcb277c2c9104f254dcd1b23f007f049c4624b` |
| `declcalc0_power_C1_res_1_2000.rds` | `bb9b3de1992095f181f4f6053da6bc6ce3020cc67ee6330150f409136dc93f9b` |
| `declcalc0_power_C2_res_1_2000.rds` | `614a4b3d569bbbe3171cf481539787647afc35e12adb5131a03e543256465771` |
| `declcalc0_power_C3_res_1_2000.rds` | `3518118ad00a18ef16b5ee0635876b8da1e0b64e3ba086baf178d16be5bf9872` |
| `declcalc0_power_C4_res_1_2000.rds` | `68357c1b8dd3ab8d9ee6d018b6c216b0c978631a0baafe17865ef043c5e339f0` |
