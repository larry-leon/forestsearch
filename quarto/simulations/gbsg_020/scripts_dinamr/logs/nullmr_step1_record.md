# nullmr Step 1 record — MR settings and coverage definitions, from source

Task `dev/tasks/TASK_null_gbsg_mr_2026-09-21.md` (committed `25798342`). Read on pop-os, 2026-09-21,
HEAD `25798342`. Paths are relative to `quarto/simulations/gbsg_020/`; the template is
`sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` (T).

## 1.1 The MR knobs

**Knobs the template reads only on the MR path.** Each feeds `mr_inference_args` (T:741-749), which
reaches `forestsearch()` only through `mr_inference_args =` (T:1275) and is consulted only when
`mr_inference = TRUE` (T:1274). With `FS_S7_MR=FALSE` every one is inert (that is what `nullid` /
`nullc125` ran).

| knob | template read | default |
|---|---|---|
| `FS_S7_FIELD_COMPLEMENT` | T:690 | TRUE |
| `FS_S7_FIELD_SCALEC` | T:716 | selected |
| `FS_S7_FIELD_DECOMP` | T:706 | FALSE |
| `FS_S7_FIELD_RECOV` | T:726 | FALSE |
| `FS_S7_IJ_RESIDUAL` | T:697 | two_term |
| `FS_S7_UNIFORM` | T:680 | FALSE |
| `FS_S7_RETURN_RESEL` | T:740 | TRUE |
| `FS_S7_WINNER_ROWS` | T:410 (display only: hides winner rows from the render's tables; recorder columns unaffected) | FALSE |

**MR literals in the template (not knobs):** `mr_draws <- 5000L` (T:545); `mr_ci_method <- "field"`
(T:673); `mr_t_confirm <- NULL` (T:674, package near-null default 1 for HR); `mr_confirm_rule <-
"point"` (T:675); `include_complement = TRUE` (T:742); field Monte Carlo sizes R_out 1000 / R_in 500
are package defaults, not forwarded (T:1914, recorded in meta `field_R`).

**Per engine, as the committed survival campaigns set them.**

| knob | FS: `scripts_p12x20/campaign_p12x20.sh` | FS cross-check: cert20 env (`REPORT_cert20_2026-09-08.md:11-15`) | DINA: `scripts_dinamr/campaign.sh` | GRF: `scripts_dinamr/grfmr.sh` | GRF: `scripts_dinamr/grfmrC.sh` |
|---|---|---|---|---|---|
| `FS_S7_MR` | unset = TRUE (:13-15 comment; not in KN :24-27) | not set (knob postdates cert20; template default TRUE) | unset = TRUE (not in KN :9-11) | unset = TRUE (:19-21) | unset = TRUE (:17-19) |
| `FS_S7_FIELD_COMPLEMENT` | TRUE (:24) | TRUE (:12) | TRUE (:9) | TRUE (:19) | TRUE (:17) |
| `FS_S7_FIELD_SCALEC` | selected (:25) | selected (:13) | selected (:10) | selected (:20) | selected (:18) |
| `FS_S7_FIELD_DECOMP` | TRUE (:25) | TRUE (:13) | TRUE (:10) | TRUE (:20) | TRUE (:18) |
| `FS_S7_FIELD_RECOV` | **unset = FALSE** (:13-15, "stay UNSET") | not set (knob postdates cert20) | **TRUE** (:10) | **TRUE** (:20) | **TRUE** (:18) |
| `FS_S7_IJ_RESIDUAL` | two_term (:26) | two_term (:13) | two_term (:11) | two_term (:21) | two_term (:19) |
| `FS_S7_RETURN_RESEL` | TRUE (:27) | TRUE (:14) | unset = TRUE | unset = TRUE | unset = TRUE |
| `FS_S7_UNIFORM` | unset = FALSE | not set = FALSE | unset = FALSE | unset = FALSE | unset = FALSE |
| `FS_S7_WINNER_ROWS` | unset = FALSE | not set = FALSE | unset = FALSE | unset = FALSE | unset = FALSE |
| `FS_S7_FB` | none (:26) | none (:14) | none (:11) | none (:21) | none (:19) |

**No disagreement within an engine.** FS: p12x20 and cert20 agree on every knob both carry;
`FIELD_RECOV` is absent from cert20 because it postdates it, and p12x20 leaves it unset (FALSE) — not a
disagreement. GRF: `grfmr.sh` and `grfmrC.sh` carry the identical KN line. Cross-check against the
committed bundles' own meta (batch `res_1_1000` of one cell per engine): cert20 FS
`field_complement TRUE, field_decompose TRUE, field_scale_complement selected, ij_residual two_term,
field_uniform FALSE, mr_draws 5000, ci_method field, fb_mode none` (no `field_recovery` field — it
predates the knob); dinamr DINA and grfmr GRF the same plus `field_recovery TRUE`.

**Between engines** the only difference is `FS_S7_FIELD_RECOV` (FS FALSE; DINA and GRF TRUE). It adds
descriptive recovery columns (`fld_recov_*`, T:718-726) that no construction reads; it is carried per
engine as committed.

**Winner-only / winner-floor.** No knob turns them on or off: `mr_*_se_w` / `_wf` and their bounds are
recorded whatever `FS_S7_IJ_RESIDUAL` selects (T:691-697, T:1133-1140); `FS_S7_IJ_RESIDUAL=two_term`
keeps them out of the reported columns and `FS_S7_WINNER_ROWS` (display) stays off. **No committed
survival campaign driver set `FS_S7_WINNER_ROWS`** (grep over `scripts_dinamr/*.sh`,
`scripts_p12x20/*.sh`: no hit). They are not evaluated or reported here. `FS_S7_FB=none`.

**Resolved for nullmr** (in `scripts_dinamr/nullmr.sh`):
- common: `FS_S7_MR=TRUE FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE
  FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none`;
- FS adds `FS_S7_RETURN_RESEL=TRUE`, leaves `FS_S7_FIELD_RECOV` unset (FALSE);
- DINA and GRF add `FS_S7_FIELD_RECOV=TRUE`, leave `FS_S7_RETURN_RESEL` unset (TRUE);
- `FS_S7_UNIFORM`, `FS_S7_WINNER_ROWS` unset (FALSE) on all three.

## 1.2 Coverage, exactly as the committed survival reports compute it

The definitions are the committed summary documents' own chunks. They are identical in form in
`summary_cert20.qmd` (FS: `cov` :74-89, `wilson-fn` :131-161) and `summary_grfmr.qmd` (GRF: `cov-fns`
:212-242, `wilson-fn` :326-359), and the same functions are in `summary_dinamr.qmd` (`cov-fns`
:158-189, `wilson-fn` :271-305) and `scripts_dinamr/blockC_numbers.R:9-32, :62-76`.
`scripts_dinamr/nullmr_findings.R` evaluates **`summary_grfmr.qmd`'s `cov-fns` and the `wil` / `covs_cell`
body of `wilson-fn` verbatim** by extracting the chunk text from the committed .qmd (the mechanism
`scripts_dinamr/grfmr_tables.R:12-21` uses), so the arithmetic is the committed arithmetic, not a
re-implementation.

- **Field lower bound on β(Ĥ):** `mean(betaHhat_H >= fld_H_lo1s)` over rows with `detected == 1` and
  `betaHhat_H`, `betaHhat_Hc`, `fld_H_lo1s`, `fld_Hc_up1s_s` all finite (`covs_cell`,
  summary_grfmr.qmd:332-337, row filter :334-335, product :337).
- **Field-s upper bound on β(Ĥᶜ):** `mean(betaHhat_Hc <= fld_Hc_up1s_s)`, same rows (:338).
- **Bonferroni pair, jointly:** `mean(betaHhat_H >= fld_joint_s_bonf_loH & betaHhat_Hc <=
  fld_joint_s_bonf_upHc)`, same rows (:346, `joint_s_bonf`) — the field-s joint's Bonferroni pair,
  i.e. the pair of the two products above at 97.5% each.
- **IJ two-term:** two-sided, `mean(betaHhat_H >= mr_H_lo & betaHhat_H <= mr_H_hi)` and the Hc
  twin, same rows (:340-341, `ij2_H` / `ij2_Hc`).
- **Naive:** `fs_sim_bias_coverage()` (R/fs_bias_coverage.R:87-181) via `cov_block()`
  (summary_grfmr.qmd:221-227): `cov1` = one-sided 95% on the exposed side (lower for β(Ĥ), upper
  for β(Ĥᶜ)) with bound `exp(log nv_est ∓ z_0.95 · nv_se)` (fs_bias_coverage.R:139-148), over
  detected rows with the bound and target finite; `cov2` = two-sided on `nv_lo / nv_hi`. The FS
  report reads the naive block as the one-sided rate (REPORT_cert20_2026-09-08.md:119, "the naive
  one-sided rate on the harm block"); both are reported.
- **Wilson 95%:** `wil()` (summary_grfmr.qmd:327-330) for the `covs_cell` products; the
  `fs_sim_bias_coverage` internal Wilson (R/fs_bias_coverage.R:116-121, same formula) for naive.
- **Bound location** (Step 4) follows `summary_grfmr.qmd` `null-location` (:1176-1210) and
  `grfmr_tables.R:116-132`: medians on the HR scale and shares `>= 1.00`, `>= 1.25` over the
  declaring rows the product is finite on, with Wilson limits. For naive and IJ two-term the lower
  bound is the one-sided 95% Gaussian bound `b1` that `fs_sim_bias_coverage()` scores `cov1` with
  (fs_bias_coverage.R:144-146; SE `nv_H_se` / `mr_H_se_ij`); for field `fld_H_lo1s`; for the
  Bonferroni pair `fld_joint_s_bonf_loH`.
