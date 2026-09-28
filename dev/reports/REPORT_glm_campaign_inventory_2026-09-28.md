# REPORT — Inventory of the GLM simulation campaigns (ACTG 175 binary and continuous), read-only (2026-09-28)

Task: `dev/tasks/TASK_glm_campaign_inventory_2026-09-28.md` (committed `65d5cea2`). Repo `forestsearch`,
branch `feature/glm-extension`, HEAD at the start `5b7b09f3`; machine `pop-os`.

**Read-only.** Nothing was run, rendered or computed. No R session was opened. Every fact below is read from
a committed record, a template or driver (cited `file:line`), a bundle directory listing, or `git log` /
`git ls-files`. Where a record does not settle a point the entry says *not established from source*. No
recommendations.

**Untracked files present at the start, read in place and left untouched** (`git status --porcelain`):
`quarto/simulations/actg175/binary_020/mr_or_harm/fs_effMaxSG_mr_field_or075_n500_nb20_redes_d5000/`
(one file, `…_redes_res_1_20.rds`, 2026-09-18 20:37),
`…/fs_effMaxSG_mr_field_or075_n500_nb20_relaunch_d5000/` (one file, `…_relaunch_res_1_20.rds`, 2026-09-18 21:09),
`quarto/simulations/actg175/binary_020/smoke_redes.html`, `…/smoke_relaunch.html`, and
`quarto/simulations/gbsg_020/scripts_dinamr/logs/nullmr_findings.err` (survival; out of scope). The three
`logs_*` directories under `continuous/` and `logs_or/` under `binary_020/` are untracked by design
(the catalogs say so) and are not listed by `git status` because they are ignored.

Paths below are relative to `quarto/simulations/actg175/` unless they start with `dev/` or `R/`.
"Catalog" means the directory's `current_status.md`.

---

## 0. Summary table

| # | campaign (directory) | outcome / measure | status | committed? | cells × replicates | results path |
|---|---|---|---|---|---|---|
| 1 | **Committed binary study** `binary_020/mr_sweep/maxeffCons_actg175_or075_seedtab_s1000/` | binary / OR | complete; the supplement's Figures S9, S10 | yes (22 files) | 21 × 1,000 (FS, DINA, GRF × n 500…2000 by 250) | `mr_coverage_grid_maxeffCons_actg175_or075_seedtab_s1000.rds` (378 rows, 21 cells); driver render |
| 2 | **`orfs`** `binary_020/mr_or_harm/fs_effMaxSG_mr_field_or*_nb20_orfs_d5000/` | binary / OR | compute complete, Gate 2 PASS 6/6; **Stage 3 (summary, extract, record, closeout) not done** | yes (6 batch + 6 combined) | 6 × 1,000 (OR 0.75 / 1.0 / 1.5 × n 500 / 2000) | 6 combined bundles; no `or_metrics.csv`; no summary render |
| 3 | **`orgrf`** `binary_020/mr_or_harm/grf_…_orgrf_d5000/` | binary / OR | compute complete, Gate 2 PASS 6/6; Stage 3 not done | yes (12) | 6 × 1,000 | 6 combined bundles; no extract |
| 4 | **`ordina`** `binary_020/mr_or_harm/dina_…_ordina_d5000/` | binary / OR | compute complete, Gate 2 PASS 6/6; Stage 3 not done | yes (12) | 6 × 1,000 | 6 combined bundles; no extract |
| 5 | **`orfs` superseded cells** `binary_020/mr_or_harm_superseded_prev09632/` | binary / OR | superseded by design change (prevalence 9.632 % → 14.917 %) | yes (6, moved aside by `git mv`) | 2 × 2,000 (OR 0.75, n 500 / 2000) | 2 combined bundles; combine renders `fs_…_orfs_combine_1_2000.html` |
| 6 | **Binary Stage 0/1 artefacts** (`orsmoke` ×4, `orcal*` ×5, `admsmoke`, `orassert`, `fs_maxeffCons_…_orsmoke`; `redes`, `relaunch`) | binary / OR | smoke / calibration; catalog says "deleted at closeout", still present | 12 tracked; `redes`, `relaunch` untracked | 20-replicate smokes; calibration batches | not results |
| 7 | **Older binary sweeps** `binary/mr_sweep/` (9 run tags), `binary/mr_sweep_legacy/` (6) | binary / OR | complete as committed 2026-06-30 … 2026-08-17; predate the task/record protocol | yes | see §5 | per-tag `*_mr_n*_res.rds` + `mr_coverage_grid_*.rds` |
| 8 | **Document-only binary directories** `binary_legacy/`, `binary_methods/`, `_archive_pre_fsparams/` | binary / OR | render-only; no payload in the repo | qmd/html tracked; no `.rds` | n/a | renders only |
| 9 | **`mdf1`** `continuous/mr_md_harm/fs_maxeffCons_mr_field_*_mdf1_d5000/` | continuous / MD | complete (Stage 3 record) | yes | 4 × 2,000 | `summary_continuous_field_mdf1.html`; no CSV extract |
| 10 | **`mdsgnb20`** `continuous/mr_md_harm/fs_effMaxSG_…_mdsgnb20_d5000/` | continuous / MD | complete | yes | 4 × 2,000 | `md_field_metrics.csv` + `COLUMNS_md_field.md`; summary html |
| 11 | **`mdgrf`** `continuous/mr_md_harm/grf_…_mdgrf_d5000/` | continuous / MD | complete | yes | 4 × 2,000 | `md_grf_metrics.csv` + `COLUMNS_md_grf.md`; summary html |
| 12 | **`mddina`** `continuous/mr_md_harm/dina_…_mddina_d5000/` | continuous / MD | complete | yes | 4 × 2,000 | `md_dina_metrics.csv` + `COLUMNS_md_dina.md`; summary html |
| 13 | **Continuous twin (IJ-only)** `continuous/mr_md_harm/*_s1000_d5000/`, `*_s100_d5000/` | continuous / MD | complete as committed 2026-08-10 … 08-29 | yes | 5 × 1,000 (md40 n500, md40 n700, md120 n500, md120 c1*, null n500) + 1 × 100 | batch renders `sim_fs_maxeffCons_mr_md*_batch_1_1000.html` |
| 14 | **FB bundles** `continuous/fb_mr_md_harm/` | continuous / MD | complete; FB joined by `mdf1` on md40 n500 sim_id 1–100; never re-run | yes (4) | 1 × 100 (`nb_boots` 300) + 3 quick-runs (20 / 100 / 1,000 at `nb_boots` 20) | bundles only |
| 15 | **Early continuous grid** `continuous/mr_coverage_sweep_md_harm.qmd` | continuous / MD | render only; its payload directory `mr_sweep_md_harm/` is absent from the repo | qmd + html tracked; no payload | FS × n 500 / 1000 / 2000 / 4000 × 1,000 (driver) | `mr_coverage_sweep_md_harm.html` |
| 16 | **Drivers with no payload in the repo**: `continuous/actg175_continuous_simulations.qmd` (MD), `survival/gbsg_poisson_simulations.qmd` (count / IRR) | MD; count / IRR | no campaign record; write to `_data` / `FORESTSEARCH_SIM_DIR`, neither present | qmd only | n/a | none |

No campaign on the binary RR or RD measures exists (`binary_020/current_status.md:275`: "The other effect
measures on the binary path (RR, RD) are untouched by this directory"). No count / IRR campaign has a payload
or record in the repo (§13).

---

## 1. The committed binary study (`maxeffCons_actg175_or075_seedtab_s1000`)

1. **Identity.** `binary_020/mr_sweep/maxeffCons_actg175_or075_seedtab_s1000/` (21 `<method>_mr_n<n>_res.rds` + the
   grid). Driver `binary_020/maxeffCons_mr_coverage_sweep_or075.qmd` (`1d42f6da`, 2026-08-17), figure fragment
   `binary_020/_sim_mr_coverage_or075.qmd`, render `maxeffCons_mr_coverage_sweep_or075.html` (grid "built
   2026-08-18 02:00:21"). No task document launched it (it predates the protocol; `binary/REPORT_actg175_binary_stage0_2026-09-17.md:77`:
   "No record (`REPORT_*`) for this study existed before this one"). Records: that Stage 0 report (read-only, §3 = the
   design as committed) and `binary_020/current_status.md` §1, §3.6. Pin: catalog `d356bd8a`.
2. **Status.** Complete. 1,000 replicates per cell in every one of the 21 bundles (`meta$n_sims = 1000`;
   Stage 0 §3, "Replicates per cell actually run: 1,000"). It is the payload behind the supplement's Figures S9 and S10
   (Stage 0 §2.3). Stage 0 F3: the supplement text's "500 simulations per cell" and the captions' "NULL MR draws" do not
   match the payload (1,000; 5,000 draws).
3. **Git.** All 22 files tracked (`current_status.md` §3.6: "22 files, 4.77 MB; 22 tracked"); "read, never rewritten".
4. **Design** (driver, per Stage 0 §3): ACTG175 arms 1 (ZDV+ddI) vs 3 (ddI), missing `cd420` dropped (`:235–238`);
   outcome `y_neg = 1 − 1{cd420 > cd40}` (`:239`); OR, `adverse_outcome = TRUE` in the analysis so OR > 1 is harm
   (`:176`); planted H = {wtkg > q70} ∩ {cd40 > q70} (`:164–168`), prevalence 9.632 %, super-population 100,000
   (`:152`); calibrated to a marginal OR of 0.75 in H (`calibrate_glm_interaction(target_effect = 0.75, …)`,
   `:249–266`), `k_inter` 0.1480184. Truths: θ†(H) 0.7499999955, θ†(Hᶜ) 0.6560116720, θ‡(H) 0.7321189340,
   θ‡(Hᶜ) 0.6313905111, **overall marginal OR 0.6665537395** (`current_status.md` §3.6; Stage 0 §3). n grid
   `seq(500, 2000, by = 250)` (`:71`); 3 identifiers × 7 sizes = 21 cells; `n_sims <- 1000L` (`:83`).
   Seeds: pre-generated table from `seed_base 8316951` indexed by global `sim_id` (`:124–133`); DGM built before any
   RNG-kind switch.
5. **Analysis settings** (Stage 0 §3, driver `:174–207`, `:484–522`): identifiers `consistency`, `dina`, `grf`;
   `sg_focus = "maxeffCons"`, `selection_rule = "neighborhood"`, `effect_neighborhood = 0.10` (inert under
   `maxeffCons`, template note `binary_020/sim_fs_mr_field_or_template.qmd:192–193`); thresholds c1 = 0.90,
   c2 = 0.80, p⋆ = 0.90 (`:192–194`); `n.min 60`, `d0/d1.min 10`, `maxk 2`, `fs.splits 500`, `use_twostage TRUE`
   (`:195–196`); GRF `frontier` / `effect` / depth 2 / `dmin.grf 0`; DINA `dina_args = list()`, `effect`.
   **MR on**, `ci_method = "ij"`, 5,000 draws, `include_complement = TRUE` (`:510–511`); FB wired but dormant
   (`nb_boots <- NULL`, `:84`). Bounds computed: **unadjusted (`nv_*`), oracle (`ora_*`), IJ two-term (`t2_*`)
   two-sided on H and Hᶜ**; no field, field-s, Bonferroni or p̂ columns (Stage 0 §3, "Estimators recorded";
   `current_status.md` §1: "It carries **no field columns**").
6. **Constructions current?** No, by its own record: `current_status.md` §1 "the unadjusted plug-in, the per-replicate
   oracle, and the two-term de-biased multiplier-resampling estimate with its IJ interval"; §2 "the campaigns … supersede
   nothing in the committed study, which ran a different rule and different constructions" (§6). The Stage 0 report
   §8 states the rule and settings "previously used vs the current campaign standard". Built on forestsearch 0.2.0 at
   `1d42f6da`, Mac M1 Ultra, 12 workers (Stage 0 §3).
7. **Results available.** The coverage grid (378 rows over 21 cells: coverage and bias of naive and MR against
   `C_betaHhat`, `C_dagger`, `C_ddagger`, `C_oracle`, over detected replicates; Stage 0 §3 "Detection and targets") and
   its render. Detection rates per identifier and n are in Stage 0 §3. No CSV extract. Summary: the driver's render and
   the manuscript fragments (Stage 0 §2.1).
8. **Flags.** Stage 0 F1 (driver lives in `binary_020/`, not `binary/`), F2 (manuscript `payload_manifest.qmd` names
   the older grid `mr_coverage_grid_actg175_or075_s1000.rds`), F3 (caption/text mismatches). `current_status.md` §2:
   its payloads "are likewise at 9.632 % and are not a comparator for anything run under the new design".

## 2. The three binary campaigns under the current constructions: `orfs`, `orgrf`, `ordina`

Shared design and settings first; per-campaign identity and status in §2.1–§2.3.

**Task and record documents.** Launched by `dev/tasks/TASK_actg175_binary_campaign_2026-09-17.md` (Stages 1–3),
redesigned by `dev/tasks/TASK_binary_study_redesign_2026-09-18.md`, relaunched by
`dev/tasks/TASK_binary_launch_v2_2026-09-18.md` (six `orfs` cells + stage-2a timing cells), completed by
`dev/tasks/TASK_binary_stage2_v3_2026-09-19.md` (ten GRF/DINA cells on the Mac Studio). Preceded by
`dev/tasks/TASK_actg175_binary_stage0_2026-09-17_v2.md` and `…_stage0_resume_2026-09-17.md` (read-only Stage 0 in
`binary/`). Side tasks that touched this directory: `TASK_binary_admission_check_2026-09-18.md`,
`TASK_audit_criterion_and_defaults_2026-09-18.md`, `TASK_threshold_naming_docs_2026-09-18.md`,
`TASK_binary_default_or_2026-09-18.md`. Records, all in `binary_020/`: `REPORT_actg175_or_stage1_2026-09-17.md`
(Gate 1), `REPORT_actg175_or_gate2_2026-09-17.md` (append-only Gate 2, 218 KB, both designs interleaved — its
own index at `:19–51`), `REPORT_binary_redesign_2026-09-18.md`, `REPORT_binary_stage1_fs_2026-09-18.md`,
`REPORT_binary_stage2a_timing_2026-09-18.md`, `REPORT_binary_stage2_2026-09-19.md`, plus
`REPORT_binary_admission_check_2026-09-18.md`, `REPORT_criterion_and_defaults_audit_2026-09-18.md`,
`REPORT_threshold_docs_2026-09-18.md`, `REPORT_binary_default_or_2026-09-18.md`. Catalog
`binary_020/current_status.md` (pin `d356bd8a`; last regenerated `e1c58844`), curated half `status_curated.md`.
Heartbeat `LOG_or_progress.txt` (last line: `2026-09-20T00:39:35Z campaign complete`).

**Design** (template of record `binary_020/sim_fs_mr_field_or_template.qmd`; catalog §2): the committed study's data
recipe, DGM builder, truths and evaluation frame verbatim ("[S] :233–341 verbatim", template `:98`), with two changes:
- `target_effect = FS_OR_TARGET` (`:207`, `:395`), giving three design points: **protective OR 0.75, borderline
  null OR 1.0, harm OR 1.5** in H against a complement that inherits the fitted ACTG175 effect (catalog §2
  "Design points"; the driver's homogeneous `dgm_model = "null"` branch is not used).
- **Planted prevalence changed 2026-09-18**: `sg_quantile <- 0.62850` (`:307`), prevalence(H) = 14.917 %,
  H = {wtkg > q} ∩ {cd40 > q} on the same two variables (`:311–313`); reason and selection in
  `REPORT_binary_redesign_2026-09-18.md` (feasibility table `feasibility_binary_redesign_2026-09-18.csv`, plateaus
  `prevalence_plateaus_binary_redesign_2026-09-18.csv`). Catalog §2: the original 9.632 % region "is undeclarable …
  in 96 % of replicates at n = 500".
- Sizes n = 500 and n = 2000 (`:210`; catalog §2 "the ends of the study's sweep"); 6 cells per campaign, 18 in all.
- `n_super 100000` (`:296`), `eval_seed 20260628` (`:315`), `outcome_type "binary"`, `effect_measure "OR"`
  (`:316–317`), `adverse_outcome = TRUE` for the analysis (`:293`), DGM constructed with `adverse_outcome = FALSE`
  (`:401`, intentional, `:383–388`).
- Replicates: **1,000 per cell, one batch over sim_id 1–1000** (`n_sims` `:146`; catalog §2 "the campaigns' 2 × 1,000
  layout was superseded on 2026-09-18"); seeds: table from `seed_base 8316951` (`:164–166`) indexed by global
  `sim_id`, so identical across identifiers and machines (stage-2 report §4 verifies 10,000/10,000 column-replicate
  matches at 1e-8).
- Trial-wide effect: not recorded per design point in the catalog. For OR 0.75 the study's overall marginal OR
  0.6665537395 applies to the 9.632 % design only; for the 14.917 % design and for OR 1.0 / 1.5: *not established
  from source* (every bundle's `meta` carries "the five truths", template `:107`, but no record prints them).
- Stage 0 feasibility gate in the template (`:500–510`): `fs_dgm_feasibility()` at n 500/750/1000/2000, tolerance
  0.05, `n_rep 200`; refuses to proceed unless `feasible` is TRUE (catalog §2).

**Analysis settings** (template; catalog §2 "Rule, in every cell"): identifier per campaign via `FS_OR_METHOD`
(`:177`); `sg_focus = "effMaxSG"` (`:179`), `effect_neighborhood = 0.20` (`:186`), `selection_rule = "neighborhood"`
(`:268`); thresholds **c1 = 0.90, c2 = 0.80, p⋆ = 0.90** on the OR scale (`:277–279`); `fs.splits 500`, `maxk 2`,
`n.min 60`, `d0/d1.min 10` (`:280`); `consistency_method = "resample"` (`:267`); `use_lasso/use_dina/use_grf FALSE`,
`use_twostage TRUE` (`:281`); GRF `frontier` / `effect` / depth 2 / `dmin.grf 0.0` (`:287–290`); DINA
`dina_args = list()`, `dina_select_statistic = "effect"` (`:291–292`). **MR on**: `ci_method = "field"` (`:331`),
5,000 draws (`:149`), `include_complement`, field complement on (`:336`), `field_scale_complement = "selected"`
(`:339`), `ij_residual = "two_term"` (`:341`), `return_reselection = TRUE` (`:343`), field R_out/R_in 1000/500 package
defaults (`meta$field_R`, `:1210`). **FB never run** (`nb_boots <- NULL` `:148`, `fb_mode <- "none"` `:354`); no
bootstrap, no cross-validation (stage-1 record header). **Bounds computed** (catalog §2 "Constructions"): the field
one-sided lower bound on Ĥ (`fld_H_*`, `:887`), the field-s one-sided upper bound on Ĥᶜ (`fld_Hc_*_s`, `:646`),
their Bonferroni pair (`joint_s`), with unadjusted, oracle and IJ two-term two-sided as references. Excluded: IJ
winner-only / winner-floor, κ / uniform, covariate adjustment, tuned inflation. The oracle helper `.logit_or_ci()` in
the template (`:552–580`) carries the legacy pooled 5/5 guard and the four-cell existence condition
(`REPORT_code_review_revalidation_2026-09-19.md` A14); the study driver's copy was reverted to byte-identical
(`23b9714d`).

**Constructions current?** Per the record, yes as of the run: catalog §2 "the three campaigns that re-run that design
under the current constructions"; every cell's `meta` records `pkg_version 0.3.5.9000`, built 2026-09-19 03:54:02 UTC
(pop-os) or 2026-09-19 20:23:40 UTC (Mac), and Gate 2 asserts `pkg_version` (catalog §2 "Which build produced
what"). That build "carries the estimability boundary (`1719056f`) and `fs_dgm_feasibility()` (`f9b794f6`)".
The template is the single in-scope copy of record (`REPORT_binary_stage1_fs_2026-09-18.md` §5, Step 1). **See §14
for the package change of 2026-09-23 that postdates these builds.**

**Results available.** Per-replicate rows in the 18 combined bundles (`*_combined_1_1000.rds`, 0.62–0.75 MB each,
catalog §3.2). Per-cell numbers in the records: declaration, sensitivity, PPV, mean |Ĥ| for `orfs`
(`REPORT_binary_stage1_fs_2026-09-18.md` §1), mean field lower bound on Ĥ / share ≥ 1.0 and mean field-s upper bound
on Ĥᶜ / share ≤ 1.0 for `orfs` (§1.1), NA-oracle / non-estimable / MR-failure counts (all zero; §3); declaration
and family-size K for the ten GRF/DINA cells (`REPORT_binary_stage2_2026-09-19.md` §2); declaration for the two
stage-2a cells (`REPORT_binary_stage2a_timing_2026-09-18.md`). **No coverage or bias table exists for any of the 18
cells**: `REPORT_binary_stage2_2026-09-19.md:8` "The 18-cell bias-and-coverage synthesis is a separate task";
`:28` "Nothing below is a coverage figure — the synthesis task carries those"; `TASK_binary_stage2_v3_2026-09-19.md:158`
"Not in this task: The 18-cell bias-and-coverage synthesis — a separate task once these ten land". The Stage 3
deliverables named in `TASK_actg175_binary_campaign_2026-09-17.md` §3.1–3.5 — the render of `summary_actg175_or.qmd`,
`or_metrics.csv`, `COLUMNS_or.md`, `REPORT_actg175_or_2026-09-17.md`, the closeout deletion of the smoke and calibration
outputs — **do not exist** (no `summary_actg175_or*.html`, no `or_metrics.csv`, no `COLUMNS_or.md` anywhere in the
tree; `git log` on `summary_actg175_or.qmd`: last commit `adafc72c`, its first commit message "Stage 3 material
prepared during Stage 1's waits — DRY-RUN VERIFIED ONLY, not run"). `summary_actg175_or.qmd` as committed reads
`*_combined_1_2000.rds` (`:55`, `:64` `glob_new <- … "combined_1_2000"`), the superseded 2 × 1,000 layout, and its
subtitle says "the eighteen committed combined bundles … writes or_metrics.csv"; the committed cells are
`*_combined_1_1000.rds`.

**Flags common to the three.** (i) Stage 3 not done (above). (ii) Catalog §3.5 says the `orsmoke` and `orcal*` rows
"are Stage 1 artefacts and are deleted at closeout"; they are still present and tracked (§4). (iii) Catalog §7 open
work: only the two ends of the sweep are run; one design family; the 12–15 % feasibility boundary unresolved;
`smoke_identity.R` `recipe` mode deferred; "No genuine global null"; RR and RD untouched; `fs_identification_summary()`
not used. (iv) `REPORT_binary_stage1_fs_2026-09-18.md` §6: Gate 2 render-line glob picks up stale logs; estimability
counts all zero because the design never loads the boundary; OR 1.5 cells cost ~10 % more; §7: an external push mid-run.
(v) `REPORT_actg175_or_gate2_2026-09-17.md` `:19–51`: the record interleaves superseded and current sections under
identical headings; discriminate by `pkg_version` and `prev`. (vi) Every `orgrf` / `ordina` figure is "conditional on
the proposed family" (catalog §5; stage-2 report header).

### 2.1 `orfs` (forest search, consistency engine)

1. **Identity.** `mr_or_harm/fs_effMaxSG_mr_field_or{075,100,150}_n{500,2000}_nb20_orfs_d5000/`; combine renders
   `fs_…_orfs_combine_1_1000.html` (6, directory root; the two `_combine_1_2000.html` belong to the superseded cells).
   Task `TASK_binary_launch_v2_2026-09-18.md` Steps 0–3; record `REPORT_binary_stage1_fs_2026-09-18.md`.
2. **Status.** Complete: 6/6 cells Gate 2 PASS `74/74`, 1,000 of 1,000 planned replicates per cell; run on `pop-os`
   at 63 workers, 2026-09-19 (catalog §3.2). Declared: 781 / 852 / 931 (n 500, OR 0.75 / 1.0 / 1.5) and 910 / 961 / 998
   (n 2000). Cell commits `13d5b6c4`, `86f543b1`, `3ee0db52`, `27170d87`, `ed607c82`, `6da4df13`.
3. **Git.** All 12 bundles tracked (catalog §3.5 `6/6`, `6/6`); 8 commits `59b144ea` … `6da4df13`. Nothing untracked
   under these six directories.
4–6. As §2.
7. **Results.** Declaration / sens / PPV / |Ĥ|, bound-location means and shares, and zero-count tables in
   `REPORT_binary_stage1_fs_2026-09-18.md` §1, §1.1, §3. No coverage, no bias, no extract.
8. **Flags.** One halt on cell 1 from a `gate2.R` edit, re-rendered (§4.2); the two `or075` cells could not run where
   their superseded bundles sat (§4.1, the move to `mr_or_harm_superseded_prev09632/`).

### 2.2 `orgrf` (GRF identifier)

1. **Identity.** `mr_or_harm/grf_effMaxSG_mr_field_or*_n*_nb20_orgrf_d5000/`; combine renders `grf_…_orgrf_combine_1_1000.html`
   (6). `orgrf_or150_n500` by `TASK_binary_launch_v2_2026-09-18.md` Step 4 (`REPORT_binary_stage2a_timing_2026-09-18.md`,
   pop-os, 63 workers, commit `5ef6b748`); the other five by `TASK_binary_stage2_v3_2026-09-19.md`
   (`REPORT_binary_stage2_2026-09-19.md`, Mac Studio `Mac-Studio-3.local`, 13 workers, 2026-09-19/20).
2. **Status.** Complete: 6/6 Gate 2 PASS `87/87`; 1,000 of 1,000 per cell. Declared 999 / 999 / 1000 (n 500) and
   999 / 1000 / 1000 (n 2000) (catalog §3.2). Stage-2 report §5: "All ten cells passed Gate 2 first time".
3. **Git.** 12 bundles tracked; 6 commits `5ef6b748` … `4353aab9`.
4–6. As §2, with `subgroup_method = "grf"`. Stage-2 report §1: "The Mac's installed build was 0.3.5 from 2026-09-11
   before this pass"; replaced by 0.3.5.9000 built 2026-09-19 20:23:40 UTC from a worktree at `b3d4921d` = HEAD,
   package source identical to the pop-os build's commit.
7. **Results.** Declaration and K in stage-2 §2; cross-machine pairing §4. No coverage / bias / extract.
8. **Flags.** Conditional-on-family reading; GRF's family is a covariate-quantile grid, its selection is surface-driven
   (stage-2 header, `R/grf_subgroup_labels.R:258–310`). GRF "declares at 0.999–1.000 in every cell including the
   benign design points" (stage-2 §5, "descriptive; the synthesis task carries the inferential reading").

### 2.3 `ordina` (DINA identifier)

1. **Identity.** `mr_or_harm/dina_effMaxSG_mr_field_or*_n*_nb20_ordina_d5000/`; renders `dina_…_ordina_combine_1_1000.html`
   (6). `ordina_or150_n500` from the stage-2a task (pop-os, `adfd783f`); the other five from
   `TASK_binary_stage2_v3_2026-09-19.md` (Mac).
2. **Status.** Complete: 6/6 Gate 2 PASS `88/88`; 1,000 of 1,000 per cell. Declared 944 / 968 / 979 (n 500) and
   870 / 934 / 987 (n 2000).
3. **Git.** 12 bundles tracked; 6 commits `adfd783f` … `50fdc6ec`.
4–6. As §2, with `subgroup_method = "dina"`; DINA's proposal floor is `log(effect.threshold)` on the link scale
   (catalog §4). Ran after the P2 orientation fix (`064fce91`; `binary/REPORT_actg175_binary_stage0_2026-09-17.md` §4).
7. **Results.** Declaration and K per cell in stage-2 §2 (K min 1 … max 6,987). No coverage / bias / extract.
8. **Flags.** Conditional-on-family reading. "DINA gets cheaper as n grows" because K shrinks (stage-2 §5);
   declaration "falls with n and rises with the design point". The `use_dina` screening path under FS "still applies the
   unoriented floor" (continuous catalog §7, thresholds STATUS C2) — not used here.

## 3. The two superseded `orfs` cells (planted prevalence 9.632 %)

1. **Identity.** `binary_020/mr_or_harm_superseded_prev09632/fs_effMaxSG_mr_field_or075_n{500,2000}_nb20_orfs_d5000/`
   (each: `res_1_1000`, `res_1001_2000`, `combined_1_2000`). Task `TASK_actg175_binary_campaign_2026-09-17.md` Stage 2;
   records `REPORT_actg175_or_gate2_2026-09-17.md` `:106`, `:208` and `REPORT_actg175_or_stage1_2026-09-17.md`.
   Combine renders `fs_…_orfs_combine_1_2000.html` (2, directory root).
2. **Status.** Complete at 2 × 1,000 = 2,000 replicates each, Gate 2 PASS; **superseded by design change**
   (catalog §2 "Two committed cells are SUPERSEDED BY DESIGN CHANGE, and they have been MOVED ASIDE").
3. **Git.** All 6 files tracked; moved by `git mv` whole (committed `50c25059` and earlier; move recorded in
   `REPORT_binary_stage1_fs_2026-09-18.md` §4.1).
4. **Design.** As §2 but `sg_quantile 0.70`, prevalence 9.632 %; OR 0.75 only; n 500 and 2000.
5. **Settings.** As §2 (same rule and MR settings); forestsearch **0.3.5 built 2026-09-17 04:47:31 UTC**, which
   "carries neither" the estimability boundary nor `fs_dgm_feasibility()` (catalog §2).
6. **Constructions.** Field / field-s / Bonferroni as §2; the build predates the launch build.
7. **Results.** Per-replicate rows in the two combined bundles; the two combine renders. No extract.
8. **Flags.** "No table may place them beside a new cell without saying so"; "they pool with nothing under the new
   design"; "One `git mv` back restores the old layout exactly" (catalog §2).

## 4. Binary Stage 0 / Stage 1 artefacts (smoke, calibration, assertion; tracked and untracked)

1. **Identity.** Under `binary_020/mr_or_harm/`: `{fs,grf,dina}_effMaxSG_mr_field_or075_n500_nb20_orsmoke_d5000/`,
   `fs_maxeffCons_mr_field_or075_n500_orsmoke_d5000/` (Stage 1 smokes, `REPORT_actg175_or_stage1_2026-09-17.md`);
   `fs_effMaxSG_mr_field_or150_n2000_nb20_orcal{16,32,32b,63}_d5000/`, `…_or150_n500_nb20_orcal63n500_d5000/`
   (calibration batches; renders `cal_orcal*.html`); `…_or075_n500_nb20_admsmoke_d5000/` (`REPORT_binary_admission_check_2026-09-18.md`,
   render `admsmoke.html`); `…_or075_n500_nb20_orassert_d5000/` (render `assert_check.html`); **untracked**
   `…_or075_n500_nb20_redes_d5000/` (the redesign's Stage 0 smoke, `REPORT_binary_redesign_2026-09-18.md` §5, render
   `smoke_redes.html`) and `…_or075_n500_nb20_relaunch_d5000/` (the launch's Stage 0 re-render,
   `REPORT_binary_stage1_fs_2026-09-18.md` §5 Step 1.7, render `smoke_relaunch.html`). Other smoke renders in the root:
   `smoke_a_recipe.html`, `smoke_consistency.html`, `smoke_dina.html`, `smoke_grf.html`.
2. **Status.** Working outputs, not results. `redes`: 14 of 20 declared; `relaunch`: 14 of 20, "matching the pre-install
   `redes` smoke exactly" (stage-1 record §5).
3. **Git.** 12 bundle directories tracked (one file each; committed with the `f9b794f6` and adjacent commits, which also
   committed the `cal_*.html`, `admsmoke.html`, `assert_check.html` renders); `redes` and `relaunch` bundles and their two
   renders untracked. `REPORT_binary_redesign_2026-09-18.md` §5: "Artefacts are Stage 0 and deliberately untracked, as
   `orsmoke`'s were" — yet the `orsmoke` bundles are tracked at HEAD (`git ls-files`).
4–6. Same template and settings as §2 (the `redes` / `relaunch` smokes at the 14.917 % design, the `orsmoke` / `orcal` /
   `admsmoke` / `orassert` at 9.632 %); `orcal*` at OR 1.5.
7. **Results.** None intended.
8. **Flags.** Catalog §3.5: "`orsmoke` and `orcal*` rows are Stage 1 artefacts and are deleted at closeout, so they are
   expected to be absent in a closed-out state" — the closeout (task §3.5) has not run. The untracked `redes` / `relaunch`
   items are listed as "pre-existing untracked (never staged)" in every later report from 2026-09-23 on.

## 5. The older binary sweeps in `binary/` (2026-06-30 … 2026-08-17)

1. **Identity.** Run tags and drivers (all `binary/`): `mr_sweep/actg175_or075_s100`, `_s500` (`mr_coverage_sweep_or075_s100.qmd`,
   `_s500.qmd`; committed `4cc13ec2` 2026-07-02); `mr_sweep/legacy-actg175_or075_s1000` (`mr_coverage_sweep_or075_s1k.qmd`,
   run_tag `actg175_or075_s1000` `:97`; `c59f8eb7` 2026-08-04 — this is the grid the manuscript's `payload_manifest.qmd`
   names, Stage 0 F2); `mr_sweep/effMaxSG_actg175_or075_s1000` (`effMaxSG_mr_coverage_sweep_or075.qmd`; `c59f8eb7`);
   `mr_sweep/maxcons_actg175_or075_s1000` (`maxcons_mr_coverage_sweep_or075.qmd`; `c59f8eb7`…`ac1cc50b`);
   `mr_sweep/maxeffCons_actg175_or075_s1000` (`binary/maxeffCons_mr_coverage_sweep_or075.qmd`, an earlier 100-replicate
   version of the study driver, `:79` `n_sims 100`, `:91` run_tag `_s100`; the tracked payload directory is named `_s1000`
   and holds 6 FS files, n 500–1750, no grid; `47eb3f6a` 2026-08-04); `mr_sweep/actg175_or20_s1000`
   (`mr_coverage_sweep_or20.qmd`; 3 FS files n 500/750/1000, no grid; `6b5db6a5` 2026-07-14);
   `mr_sweep/actg175_or35_effMaxSG_s100`, `_s500` (`mr_coverage_sweep_or35_effMaxSG.qmd`; `5748a2d2` 2026-07-08);
   `mr_sweep/actg175_or35_s1000` (`mr_coverage_sweep_or35_s1k.qmd`; `7da89d3a` 2026-07-09);
   `mr_sweep_legacy/actg175_or075_s1000`, `_s20`, `actg175_or10_s1000`, `actg175_or15_effMaxSG_s100`, `_s1000`,
   `actg175_or15_s500` (`mr_coverage_sweep_or10.qmd`, `_or15.qmd`, `_or15_effMaxSG.qmd`; all `a4bd2393` 2026-06-30).
   **No task document and no report** for any of them; the only record is the table in
   `binary/REPORT_actg175_binary_stage0_2026-09-17.md` §2.2 and its §5 "Harm and null designs, as facts".
2. **Status.** Complete as committed; replicates per cell per driver: `s100` 100, `s500` 500, `s1000` 1,000, `s20` 20
   (each driver's `n_sims`, e.g. `mr_coverage_sweep_or075_s1k.qmd:85`, `_or15.qmd:79` = 500, `_or35_effMaxSG.qmd:79` = 500).
   Rows per bundle: *not established from source* (not opened).
3. **Git.** Every `.rds` under `binary/mr_sweep/` and `binary/mr_sweep_legacy/` is tracked; nothing untracked under
   `binary/`.
4. **Design.** Same builder as §1 with `target_or_h` 0.75 (`_or075_s1k.qmd:131`), 1.0 (`_or10.qmd:125`), 1.5
   (`_or15.qmd:125`, `_or15_effMaxSG.qmd:125`), 2.0 (`_or20.qmd:131`), 3.5 (`_or35_s1k.qmd:131`, `_or35_effMaxSG.qmd:125`);
   `k_inter_range = c(0.3, 1.5)`; H = {wtkg > q} ∩ {cd40 > q} with `sg_quantile <- 0.70` (`_or075_s1k.qmd:132`, `_or15.qmd:126`, `_or10.qmd:126`, `_or35_s1k.qmd:132`, `effMaxSG_…:132`; `subgroup_cuts` at `:143` / `:137`); `n_super` 100,000
   (`_or075_s1k.qmd:130`) but **25,000 in `_or15.qmd:124`**; n grid 500…2000 by 250 (`:77` / `:71`), except
   `effMaxSG_…` and `maxcons_…` at n 500 and 2000 only (`:77`, "by = 1500"); identifiers `consistency, dina, grf` (`:76` /
   `:70`); `seed_base 8316951` (`:112` / `:106`); `dgm_model "alt"` (`:129` / `:123`). Trial-wide effect: *not established
   from source* for the non-0.75 points.
5. **Settings.** `sg_focus = "eff"` (`_or075_s1k.qmd:158`, `_or10.qmd:152`, `_or15.qmd:152`, `_or20.qmd:158`,
   `_or35_s1k.qmd:158`; "GLM-natural alias of hrMaxSG"), `"effMaxSG"` (`effMaxSG_…:158`, `_or15_effMaxSG.qmd:152`,
   `_or35_effMaxSG.qmd:152`), `"maxcons"` (`maxcons_…:158`), `"maxeffCons"` (`binary/maxeffCons_…:152`);
   `effect_neighborhood 0.10` (`:160` / `:154`); c1 0.90, c2 0.80, p⋆ 0.90 (`:161–163` / `:155–157`); **MR
   `ci_method = "ij"`** (`:436` / `:409`, `:430`), 5,000 draws; FB dormant. Bounds: unadjusted, oracle, IJ two-sided only.
6. **Constructions current?** No (IJ-only, pre-field; Stage 0 §8 "Rule and settings previously used vs the current campaign standard").
7. **Results.** Per-tag coverage grids `mr_coverage_grid_<tag>.rds` (absent for `maxeffCons_actg175_or075_s1000` and
   `actg175_or20_s1000`); renders `mr_coverage_sweep_or*.html`, `effMaxSG_…html`, `maxcons_…html`, `_sim_mr_coverage_or075.html`,
   `_sim_mr_coverage_or15.html`. No CSV extracts. `mr_coverage_sweep_or15.qmd` was the manuscript's inactive (escaped)
   include (Stage 0 §2.1).
8. **Flags.** Stage 0 F1/F2 (which driver / grid the manuscript names); the `_s1000`-named `maxeffCons` payload holds a
   100-replicate driver's output — per the driver, not verified from the bundle.

## 6. Document-only binary directories (`binary_legacy/`, `binary_methods/`, `_archive_pre_fsparams/`)

1. **Identity.** `binary_legacy/actg175_binary_m1_harm_<focus>_fs{1..8}.qmd` (48 per-config articles, focus × FS-bundle grid,
   the `_fsparams/` bundle system described in `README.qmd:43–102`), `binary_legacy/_summaries/`; `binary_methods/` (the
   `m1b`–`m1f` "methods" documents, `adjusted` / `resample` variants, `fdr_*.R`, `fs_profile_harness.R`, verification scripts);
   `_archive_pre_fsparams/` (3 qmd). No task or report documents.
2. **Status.** *Not established from source* beyond "rendered" (html present for some); README says bundles are "built in
   full mode" as `.rds` beside the qmd — none exist in the tree.
3. **Git.** qmd / html / R / sh tracked; **no `.rds` tracked or on disk** under these three directories.
4–6. Pre-field, pre-`effMaxSG` era; the `binary_methods/fdr_mr_inference.R` and `fdr_family_multiplier.R` scripts carry
   the 2026-09-27 header "Verifies pre-alignment (exact-cutoff) behaviour, before the 2026-09-23 alignment (`96f84ad8`,
   `7713942e`); not valid against the aligned package" (`b8970c47`; `dev/reports/REPORT_declcal_consumers_2026-09-24_v3.md` §9.1 class (a)).
7. **Results.** Renders only.
8. **Flags.** As 4–6.

## 7. Continuous campaign `mdf1` (maxeffCons, unstudentized complement field)

1. **Identity.** `continuous/mr_md_harm/fs_maxeffCons_mr_field_{md40_knoise0_n500, md120_knoise0_n500, mdnull_knoise0_n500, md40_knoise0_n700}_mdf1_d5000/`
   (each: `res_1_1000`, `res_1001_2000`, `combined_1_2000`, `gate2_flips.txt`); combine renders
   `fs_maxeffCons_mr_field_<cell>_mdf1_combine_1_2000.html`; summary `summary_continuous_field_mdf1.{qmd,html}`; scripts
   `scripts_mdf1/` (15). Task `dev/tasks/TASK_continuous_field_mac_2026-09-07.md` (handoff
   `dev/tasks/HANDOFF_continuous_field_mac_2026-09-07.md`). Records `REPORT_continuous_field_{stage0,stage1,gate2}_2026-09-07.md`,
   `REPORT_continuous_field_2026-09-07.md` (Stage 3). Catalog `continuous/current_status.md` (pin `c767027d`).
2. **Status.** Complete: 4 cells × 2,000 replicates (2 batches of 1,000), declared 1998 / 2000 / 1993 / 1999 (catalog §3.2);
   Mac Studio M4 Max, 13 workers, 2026-09-07, branch `feature/glm-extension-mac`, forestsearch 0.3.5.
3. **Git.** 12 bundles + 4 `gate2_flips.txt` tracked (catalog §3.5); 4 commits `c13d0882` … `c80a1291`.
4. **Design** (catalog §1; Stage 3 record "Design and conventions"; template `continuous/sim_fs_maxeffCons_mr_field_md_template.qmd`):
   ACTG175 arms 1 and 3, `cd4_change = cd420 − cd40`, MD, ddI = 1, `adverse_outcome = FALSE` (`:284`; harm = negative raw
   MD, every recorded estimate harm-oriented); true region `age > 34 & preanti <= 744.5` (`:289–290`, `:387–388`),
   prevalence 0.345; `n_super 5000` (`:291`); planted harm-region MD 40 (`calibrate_glm_interaction(target_effect = −40)`,
   `:430–433`; `beta_inter −13.7`) at n 500 and 700, MD 120 built directly at the locked `k_inter −93.7447641240`
   (`:417–424`), and the null cell (`model = "null"`, `:408–410`; homogeneous −26.26 raw); complement −26.26 raw in every
   cell; **ITT −30.99 for md40** (`REPORT_continuous_field_stage0_2026-09-07.md`, "DGM and harm orientation"); ITT for
   md120 / null: *not established from source*. Seeds `8316951 + sim_id` (`:131`, `:663`), L'Ecuyer-CMRG; sim_id 1–1,000
   are the twin's committed seeds. `k_random_noise 0` (`:123`).
5. **Settings** (catalog §1; template): consistency engine, `consistency_method = "resample"` (`:242`), `sg_focus =
   "maxeffCons"` (`:145`), `effect_neighborhood 0.10` inert (`:154`, `:253–255`), `selection_rule "neighborhood"` (`:245`);
   thresholds **c1 = 30, c2 = 10** on the harm-oriented MD (`:257`, `:259`; `:259` carries "OPEN QUESTION D2"), **p⋆ = 0.90**
   (`:261`); `fs.splits 400`, `maxk 2`, `n.min 60`, `d0/d1.min 12` (`:262`); J = 10 cuts on `age`, `preanti` (`:268`),
   `str2` in the pool (`:302`); MR `ci_method "field"` (`:317`), 5,000 draws (`:122`), complement field on (`:325`),
   `ij_residual two_term` (`:334`), `return_reselection TRUE` (`:338`), field R_out/R_in 1000/500; **FB joined**
   (`FS_MD_FB=join`, `:355–359`) on md40 n500 sim_id 1–100 from the committed FB bundle. Bounds computed: unadjusted,
   oracle, IJ two-term, **field lower on Ĥ, unstudentized field upper on Ĥᶜ (`fld_Hc_up1s`)**, joint pair separate /
   Bonferroni / calibrated (Stage 3 "Design and conventions"). No field-s columns (catalog §5: "its recorder has no `_s`
   column").
6. **Constructions current?** **No, per the catalog** (§5 "Superseded — do not quote as the current product"): its
   unstudentized complement bound is "Superseded by the studentized complement field (field-s, `fld_Hc_up1s_s`), the
   gate's default construction since `fb62705c`"; its rule "is not the survival grid's". Its Ĥ-block field figures are
   "read beside `mdsgnb20`'s in the paired rule contrast".
7. **Results.** Full tables in `REPORT_continuous_field_2026-09-07.md` (Table-2 layout both blocks; bound location; joint
   pair; regime; display; findings 1–9) and `summary_continuous_field_mdf1.html`; **no CSV extract** (the `md_field_metrics.csv`
   schema was introduced by `mdsgnb20`, which carries 28 paired `mdf1` comparator rows). Headline (record finding 1): field
   one-sided lower coverage on Ĥ 0.950 / 0.948 / 0.950 / 0.947; complement upper 0.924 / 0.934 / 0.926 / 0.936.
8. **Flags.** Superseded as the complement product (6). `scripts_mdf1/reselection_check.R` carries the 2026-09-27
   "pre-alignment (exact-cutoff) … not valid against the aligned package" header (`b8970c47`). Its render logs live in the
   Mac scratchpad, not the repo (catalog §6).

## 8. Continuous campaign `mdsgnb20` (effMaxSG ε 0.20, field-s)

1. **Identity.** `continuous/mr_md_harm/fs_effMaxSG_mr_field_<cell>_nb20_mdsgnb20_d5000/` (4 cells as §7); renders
   `fs_effMaxSG_mr_field_<cell>_nb20_mdsgnb20_combine_1_2000.html`; summary `summary_continuous_field_mdsgnb20.{qmd,html}`;
   extract `md_field_metrics.csv` + `COLUMNS_md_field.md`; heartbeat `LOG_mdsgnb20_progress.txt`; scripts `scripts_mdsgnb20/`
   (also the catalog generator / checker). Tasks `dev/tasks/TASK_md_field_rerun_stage0_2026-09-15.md`,
   `TASK_md_field_rerun_2026-09-15.md`. Records `REPORT_md_field_rerun_{stage0,stage1,gate2}_2026-09-15.md`,
   `REPORT_md_field_rerun_2026-09-15.md`.
2. **Status.** Complete: 4 × 2,000; Gate 2 58/58 per cell; declared 1998 / 2000 / 1993 / 1999; `pop-os`, 63 workers,
   2026-09-15/16; forestsearch 0.3.5 built 2026-09-16 05:57:14 UTC; 9,783 s total. Cell commits `72e6f129`, `b12f5983`,
   `dbeee517`, `f528d188`; campaign complete `2e4540c6`.
3. **Git.** 12 bundles tracked; 4 commits `72e6f129` … `f528d188`. `logs_mdsgnb20/` (38 files) untracked by design.
4. **Design.** As §7, same seeds and draws (Gate 2 checks `n_true` identical and oracle columns ≤ 2.8e-11 vs `mdf1`, both
   directions).
5. **Settings.** As §7 except `sg_focus = "effMaxSG"`, `effect_neighborhood = 0.20`, `selection_rule = "neighborhood"`
   (`FS_MD_FOCUS` / `FS_MD_NBHD`; record "Gate 0 dispositions" D1); `field_scale_complement = "selected"` (`:330`);
   **FB off** (`fb_mode "none"`). Bounds computed (record header): "the field one-sided lower bound on β(Ĥ), the field-s
   one-sided upper bound on β(Ĥᶜ), and their Bonferroni pair, with unadjusted, oracle and IJ two-term as references";
   the unstudentized field recorded beside field-s (catalog §2).
6. **Constructions current?** Yes, per its record and the catalog ("re-runs `mdf1`'s four cells … and records the current
   constructions", catalog §2; "the gate's default construction since `fb62705c`", §5). Build 0.3.5 (2026-09-16). See §14.
7. **Results.** `md_field_metrics.csv` (1,114 rows: 1,086 `mdsgnb20` + 28 paired `mdf1` comparators; one row per cell ×
   block × estimator × metric with Monte Carlo SEs; definitions `COLUMNS_md_field.md`); tables pasted in
   `REPORT_md_field_rerun_2026-09-15.md` (declaration, Ĥ, Ĥᶜ, ladder, joint, identification, rule contrast);
   `summary_continuous_field_mdsgnb20.html`; figures `fig_mdsgnb20_bias_coverage_display_{H,Hc}.png`.
8. **Flags** (record "Findings"): checker halt on cell 1 (fixed `34580fb0`); field-s equals the unstudentized field on this
   design to within 0.003; complement one-sided coverage below nominal at n = 500 (0.914–0.917); md120 field two-sided 0.916
   vs one-sided lower 0.967. Scope: "They do not verify condition (A3) on the GLM paths, and no construction is promoted".

## 9. Continuous campaign `mdgrf` (GRF identifier)

1. **Identity.** `continuous/mr_md_harm/grf_effMaxSG_mr_field_<cell>_nb20_mdgrf_d5000/`; renders `grf_…_mdgrf_combine_1_2000.html`;
   summary `summary_continuous_field_mdgrf.{qmd,html}`; extract `md_grf_metrics.csv` + `COLUMNS_md_grf.md`; heartbeat
   `LOG_mdgrf_progress.txt`; scripts `scripts_mdgrf/`. Tasks `dev/tasks/TASK_md_dina_grf_stage0_2026-09-16.md`,
   `TASK_md_grf_2026-09-16.md`, `TASK_md_grf_resume_2026-09-16.md`, and the fix task `TASK_grf_dina_fixes_2026-09-16.md`.
   Records `REPORT_md_dina_grf_stage0_2026-09-16.md`, `REPORT_md_grf_stage1_2026-09-16.md` (stopped),
   `REPORT_grf_dina_fixes_2026-09-16.md`, `REPORT_md_grf_stage1_resume_2026-09-16.md`, `REPORT_md_grf_gate2_2026-09-16.md`,
   `REPORT_md_grf_2026-09-16.md`.
2. **Status.** Complete: 4 × 2,000; Gate 2 66/66 per cell; declared 2,000 in every cell; `pop-os`, 63 workers, 2026-09-17;
   forestsearch 0.3.5 built 2026-09-17 04:47:31 UTC (after fix P1 `0cd33f7b`); 7,423 s. Cell commits `c4c49572`,
   `f84261e3`, `15bee7cd`, `432f69dc`; complete `1a65b28a`.
3. **Git.** 12 bundles tracked; `logs_mdgrf/` (39) untracked by design.
4. **Design.** As §7/§8, same draws (Gate 2: `n_true` identical, oracle max relative difference 0).
5. **Settings.** `subgroup_method = "grf"` (`FS_MD_METHOD`, `:139`), `grf_selection "frontier"`, `grf_select_statistic
   "effect"`, `grf_depth 2` (`:273–275`), **`dmin.grf = 30`** (`:280`; floors the DR-score pre-filter only); `effMaxSG` /
   0.20 / `neighborhood`; MR as §8; FB off. Bounds: field on Ĥ, field-s on Ĥᶜ, Bonferroni; unadjusted, oracle, IJ references.
6. **Constructions current?** Yes per its record ("the same constructions"); build 0.3.5 (2026-09-17). See §14.
7. **Results.** `md_grf_metrics.csv` (`md_field_metrics.csv` schema plus `identifier`; FS comparator rows copied);
   `REPORT_md_grf_2026-09-16.md` tables; `summary_continuous_field_mdgrf.html`; `fig_mdgrf_…png`.
8. **Flags.** "Conditional on the proposed family" on every figure; complement one-sided coverage below nominal in every
   cell (Wilson upper 0.931–0.942); family size identical across the three n = 500 cells; `status_inventory.R` has no
   `mdgrf` rows (record finding 9).

## 10. Continuous campaign `mddina` (DINA identifier)

1. **Identity.** `continuous/mr_md_harm/dina_effMaxSG_mr_field_<cell>_nb20_mddina_d5000/`; renders `dina_…_mddina_combine_1_2000.html`;
   summary `summary_continuous_field_mddina.{qmd,html}`; extract `md_dina_metrics.csv` + `COLUMNS_md_dina.md`; heartbeat
   `LOG_mddina_progress.txt`; scripts `scripts_mddina/`. Task `dev/tasks/TASK_md_dina_campaign_2026-09-17.md`. Records
   `REPORT_md_dina_stage1_2026-09-17.md`, `REPORT_md_dina_gate2_2026-09-17.md`, `REPORT_md_dina_2026-09-17.md`.
2. **Status.** Complete: 4 × 2,000; Gate 2 65/65 per cell; declared 1996 / 2000 / 1987 / 1997; `pop-os`, 63 workers,
   2026-09-17; forestsearch 0.3.5 built 2026-09-17 04:47:31 UTC (after fix P2 `064fce91`); 27,727 s. Cell commits
   `1856c107`, `c9f1b948`, `9051b702`, `21e48ee5`; complete `30b282eb`; catalog closeout `1d122b21`.
3. **Git.** 12 bundles tracked; `logs_mddina/` (41) untracked by design.
4. **Design.** As §7/§8, same draws.
5. **Settings.** `subgroup_method = "dina"`, `dina_args = list()`, `dina_select_statistic = "effect"` (`:282–283`);
   proposal and admission floors 30 on the harm-oriented MD after P2; `effMaxSG` / 0.20 / `neighborhood`; MR as §8; FB off.
   Bounds as §8.
6. **Constructions current?** Yes per its record; build 0.3.5 (2026-09-17). See §14.
7. **Results.** `md_dina_metrics.csv` (`md_grf_metrics.csv` schema; FS and GRF comparator rows copied);
   `REPORT_md_dina_2026-09-17.md` tables incl. the three-identifier table; `summary_continuous_field_mddina.html`;
   `fig_mddina_…png`.
8. **Flags.** Conditional-on-family; MR's family is DINA's whole proposed family (single-member families occur);
   complement one-sided coverage below nominal in every cell (Wilson upper 0.906–0.947); the `use_dina` screening path under
   FS still applies the unoriented floor (out of scope); cost 100–305 s per replicate, peak 123 GB.

## 11. The continuous twin (IJ-only) and the FB bundles

1. **Identity.** Twin bundles `continuous/mr_md_harm/fs_maxeffCons_mr_{md40_knoise0_n500, md40_knoise0_n700, md120_knoise0_n500, mdnull_knoise0_n500}_s1000_d5000/`
   (`…_res_1_1000.rds`), `…_md120_knoise0_n500_c1star_s1000_d5000/` (`…_c1star_res_1_1000.rds`), `…_md40_knoise0_n500_s100_d5000/`;
   drivers `sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd` (the twin of record, `REPORT_continuous_field_stage0_2026-09-07.md:7`),
   `…_n700_batch_1_1000.qmd` (differs in `n_sample <- 700L`, `:133`), `…_md120_knoise0_n500_batch_1_1000.qmd` (adds the
   Stage-2 direct build and the `_c1star` stem, `:130`, `:172`), `…_batch_1_100.qmd`; renders `sim_fs_…_batch_1_1000.html`
   (md40 n500, md40 n700; the md120 and null renders are not in the tree). Tasks: the `cc_task_oc_*` series
   (2026-08-28 … 09-01) reads them; `cc_task_oc_breadth_stage2_2026-08-31.md` produced the md120 and c1* cells
   (`6504e0ea`, 2026-08-29); md40 cells `2315b8b4` (2026-08-10, "re-run both production cells with computed scale and
   corrected oracle bounds" `bb75cca6`); null `2b180813` (2026-08-26). No `REPORT_*` of their own in `continuous/`; the
   catalog §2 lists them as "Earlier material". FB bundles `continuous/fb_mr_md_harm/fs_maxeffCons_fb_mr_md40_knoise0_n500_s100_d5000/`
   (`nb_boots 300`, pkg 0.2.0) and three `_quickrun_s{20,100,1000}_d5000/` (`nb_boots 20`, "illustrative only",
   driver `:19–20`); driver `o2_fb_mr_batches.R`.
2. **Status.** Complete as committed: twin 1,000 rows each, detected 1000 / 999 / 1000 / 786 (c1*) / 998
   (`REPORT_continuous_field_stage0_2026-09-07.md`, "Committed cells" table); packages 0.2.2 / 0.3.1. FB: 100 replicates
   at `nb_boots 300` on md40 n500 sim_id 1–100; "never re-run" (catalog §2).
3. **Git.** All tracked (catalog §3.5: `_s1000_d5000/**` 5/5, `_s100_d5000/**` 1/1, `fb_mr_md_harm/**` 4/4).
   `…_mr_field_md40_knoise0_n500_s0dg_d5000/` is an **empty tracked-nothing directory** (0 files; the Stage 0 stem used by
   `REPORT_md_dina_grf_stage0_2026-09-16.md:52` with "nothing saved").
4. **Design.** As §7 (the same DGM; the template copies the twin). `c1*`: `effect.threshold = 135.741` (breadth forecast
   scoring, "a non-standard threshold", Stage 0 report).
5. **Settings.** `sg_focus "maxeffCons"` (`:127`), `effect_neighborhood 0.10` pinned (`:169`), c1 30, c2 10 (`:171`, `:173`),
   p⋆ 0.90, `n.min 60`, `maxk 2`; **MR `ci_method = "ij"`** (`:590`), 5,000 draws (`:115`); `nb_boots NULL` (`:112`).
   Bounds: unadjusted, oracle, IJ two-sided (and FB in the FB bundles).
6. **Constructions current?** No (IJ-only; the template's `ci_method = "ij"` "reproduces the committed bundles", template `:42`).
7. **Results.** The batch renders' summary layers (md40 n500, n700); the OC-wrapper documents (`oc_wrapper_verification.qmd`,
   `analytic_verification_and_prediction_md_harm.qmd`) read them. No extract.
8. **Flags.** `oc_wrapper_verification.qmd` carries the 2026-09-27 "Values computed at the exact cutoff, before the
   2026-09-23 alignment" header (`b8970c47`).

## 12. Early continuous grid `mr_coverage_sweep_md_harm.qmd` (render only)

1. **Identity.** `continuous/mr_coverage_sweep_md_harm.qmd` (`run_tag "md_grid_s1000"`, `:76`; results to
   `mr_sweep_md_harm/<run_tag>/`, `:77`) and its render `mr_coverage_sweep_md_harm.html` (5.24 MB; commit `c489c574`
   2026-08-08 "docs(stopB): md/harm grid — MR residual bias GROWS with n, coverage falls to 0.731"); summary
   `md_harm_mr_simulation_summary.{qmd,html}` reads `mr_sweep_md_harm/md_grid_s1000` and a `md_harm_s50_pilot` (`:51–52`).
2. **Status.** Rendered; **the payload directory `mr_sweep_md_harm/` does not exist in the repo.**
3. **Git.** qmd and html tracked; no payload.
4. **Design.** MD harm −40 (`:112`), FS only (`:81`), n 500 / 1000 / 2000 / 4000 (`:84`), 1,000 per cell (`:85`),
   `seed_base 8316951` (`:92`); thresholds 30 / 10 (`:121–122`).
5. **Settings.** MR `ci_method "ij"` (`:323`); FB dormant (`:87`).
6. Not current (IJ-only, `maxeffCons`).
7. **Results.** In the two renders only.
8. **Flags.** The commit message's finding ("coverage falls to 0.731") is the only recorded reading.

## 13. Drivers with no payload in the repo (MD, and the count / IRR path)

- `continuous/actg175_continuous_simulations.qmd` (`nsims_alt 2500`, `nsims_null 5000`, `:27–28`; writes to `_data`,
  `:185`; last commit `53def0de` 2026-05-20). No `_data` directory, no render, no record.
- `survival/gbsg_poisson_simulations.qmd` (Poisson / IRR on GBSG recurrence with a person-time offset; `nsims_alt 1000`,
  `nsims_null 2000`; writes to `FORESTSEARCH_SIM_DIR/quarto/_data`, `:106`; last commit `cb1fd19d` "Redirect heavy outputs
  to FORESTSEARCH_SIM_DIR (machine-local)"). No payload under `quarto/_data`, no render in `survival/`, no task or record.
  **This is the only count-outcome simulation driver found; no count campaign exists in the repo.** (`quarto/applications/count_data_demo.qmd`
  and `quarto/guides/count_data_hte_summary.qmd` are single-dataset demonstrations, not campaigns.)
- Status, cells run, results: *not established from source* (nothing in the repo records a run).

## 14. Cross-cutting: what the record flags after the campaigns closed

- **Package change after every GLM campaign build.** The MR admission rule was aligned with the screen's rounded rule on
  2026-09-23 (`dev/reports/REPORT_mr_admission_alignment_2026-09-23.md`; commits `7713942e`, `96f84ad8`, `06ac5391`
  per `REPORT_section5_full_rerun_2026-09-24.md` §0), with `pconsistency.digits` passed through `forestsearch()` to
  `fs_mr_inference()`. Every GLM campaign bundle above was produced on an earlier build (binary: 0.3.5.9000 built
  2026-09-19; continuous: 0.3.5 built 2026-09-07 / 09-16 / 09-17; twin and older sweeps earlier still). The alignment
  report's own scope was the GBSG survival application; it states "The Section 5 simulation re-run is a separate
  decision" and "Whether the Section 5 re-run is a refresh or a finding depends on how often the simulated fits carry
  candidates near the band". The Section 5 re-run that followed (`REPORT_section5_full_rerun_2026-09-24.md`) covers the
  18 GBSG survival cells only: "declaration and selected region are identical on every replicate of every cell, but the
  re-run does move three published range endpoints by one replicate each (0.0005)". **No record re-runs, re-reads or
  classifies any ACTG175 binary or continuous campaign bundle under the aligned package**; the declcal-consumers work
  (`REPORT_declcal_consumers_2026-09-24_v3.md` §9.1) touched only four ACTG175 files, all verification / theory scripts
  (`binary_methods/fdr_mr_inference.R`, `fdr_family_multiplier.R`, `continuous/scripts_mdf1/reselection_check.R`,
  `continuous/oc_wrapper_verification.qmd`), adding the header "before the 2026-09-23 alignment … not valid against the
  aligned package". Whether the aligned rule changes any number in §§1–12: *not established from source*.
- **Thresholds workstream (2026-09-18)** (`dev/reports/STATUS_thresholds_workstream_2026-09-18.md`): binary default
  `effect_measure` resolved to `"OR"`; replicates resolve the parent fit's thresholds; `c2 > c1` stops and a silent `c2`
  derives `0.80 × c1` on the consistency path for HR / OR; DINA refuses identity-scale estimands. The binary campaigns
  (2026-09-19) ran on a build containing these; the continuous campaigns (≤ 2026-09-17) predate them. The continuous
  campaigns pass both thresholds explicitly (template `:694`), and the binary template passes both (`:754–755`).
- **Binary Stage 3 undone** (§2): no coverage / bias / extract / record / closeout for the 18 cells; the summary source
  is written against the superseded `_1_2000` layout.
- **Untracked working files** (§4) unchanged since 2026-09-18, listed as pre-existing in every report since 2026-09-23.
- **Catalog pins.** `binary_020/current_status.md` pin `d356bd8a`, regenerated `e1c58844`; no bundle commit after it.
  `continuous/current_status.md` pin `c767027d` (`1d122b21`); the only later commits touching `continuous/` are
  `f9b794f6` (no bundle change) and `b8970c47` (header lines).

---

*Written 2026-09-28. Nothing run, nothing computed, no file outside `dev/tasks/` and `dev/reports/` written; the
untracked files listed at the top are unchanged.*
