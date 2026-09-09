# HANDOFF — Linux work stream: MR field development (survival), state at 2026-09-08

Date: 2026-09-08. From: the Linux MR-field chat (2026-09-08, took over from `HANDOFF_mr_field_linux_2026-09-07.md`). For: a fresh chat that will review CC reports and write specs for the Linux (Pop!_OS, ~127 cores, 100-worker) stream. Self-contained; every record cited lives in the repository (`larry-leon/forestsearch`, branch `feature/glm-extension`, HEAD 6d8d36d9 once pushed; the Mac branch `feature/glm-extension-mac` is merged in at c0f48a7c). Records are under `quarto/simulations/gbsg_020/` unless noted; task documents under `dev/tasks/`; the decision NOTE under `dev/notes/`.

## 1. Roles and protocol (unchanged, with additions)

- Larry decides and approves; the chat writes task documents, proposals, and reviews; Claude Code (CC) does every repository operation. CC never pushes; Larry pushes via GitHub Desktop.
- **Summaries in bullet form** (standing procedure, Larry 2026-09-08): reviews, campaign summaries, status reports — one item per bullet, short; detail lives in the accompanying file.
- Task documents travel as `.md` via `~/Downloads`; CC copies to `dev/tasks/` and commits first. **The Downloads transport has failed four times** — every kickoff must state the full Stage-0/gate delta so CC can reconstruct the spec exactly; if reconstruction happens, the versioned file is committed afterwards as the governing spec (precedent: `TASK_complement_location_2026-09-08_v2.md`, 6d8d36d9). Amended task documents get new filenames (`_v2`).
- Kickoffs short and pasteable, the absolute final content of the message. Compute is a separate go: cost stated (cells, replicates, wall), pre-authorized only by Larry pasting the kickoff, verification renders (≤ 5 replicates) distinguished from campaign compute.
- Gates stop on failure without asking; Gate 1 = compute go/no-go against a ceiling; Gate 2 per cell; identity to committed bundles is the pairing proof (all pre-existing non-timing columns `identical()`, `truth` identical; timing columns excluded). Never re-run committed work.
- Any `R/` change classified: adds code / changes behaviour / changes the method; add-only with byte-identical defaults is the standard; method constructions enter **add-beside** (new outputs on an enabled path, defaults byte-identical), the pattern by which `field` and `field-s` entered. `devtools::install()` before parallel runs, never `load_all()`.
- Verify from source; stop-point reviews require the record's verbatim numbers; no vacuous checks; transplant-first for documents.
- Interpretation: bounds by location against clinically meaningful effect sizes, never significance at the null; Wilson intervals; bias in SD units — **marginal SD (across replicates, carries target spread) and error SD (of the de-biased error) shown side by side** (Part A, A0).
- **Closed lines (do not reopen):** the submitted JRSS-B paper; κ variants and the hybrid κ; IJ winner-only and winner-floor (excluded from every table, figure, and report); covariate adjustment or alternative working models; attribution of cross-version selection flips; interpolation between dial points (Larry: not trustworthy); ε > 0.25 as an adoptable band; pursuit of Xu & Guo (2026) or Bargagli-Stoffi & Melikechi (2026) beyond citation (Larry 2026-09-08; the only future use of Xu–Guo is a cast into an interpretable identifier with our inference, no comparison); tuned inflation factors of any kind for the complement's residual.

## 2. Constructions and standing dispositions

- **Harm subgroup Ĥ, lower bound on β(Ĥ): the field** (`ci_method = "field"`; est₂; one-sided lower). Coverage 0.936–0.980 across every survival cell (s7, map1, p30, p30sg, nb20, banddial, e1stud); inside ε ≤ 0.25: 0.964–0.976 (n = 500), 0.956–0.959 (n = 1000). Field SE/error-SD 0.96–1.12; the re-selection mixture supplies the selection-inflated scale the naive SE lacks (naive SE/error-SD 0.83–0.89; naive one-sided 0.38–0.68). Retained bias grows with the band (−0.09/−0.10 at ε 0.30). **Field-s is not applied to Ĥ** — the harm error is not naive-SE-scaled. IJ two-term: conservative one-sided option (0.97–0.999, 1.2–1.6× wider) and the **two-sided reporting default** (`ci_method = "ij"`; the field's two-sided under-covers in harm cells, 0.91–0.94).
- **Complement Ĥᶜ, upper bound on β(Ĥᶜ): field-s** (`field_scale_complement = "selected"`, R1 studentized complement field; `dev/notes/NOTE_complement_product_2026-09-08.md`, decided by Larry 2026-09-08). Coverage 0.912–0.921 at n = 500 in the band cells (field 0.897–0.912; naive 0.774–0.826; IJ two-term 0.995–0.996 at 1.6–1.8× width, rules out nothing); two-sided 0.933–0.941; never worse than field in any cell; ends unchanged (0.926 maxSG, 0.929 minSG). **Stated shortfall ≈ 3 points below 0.95**: ≈ 0.6 from residual scale (SE/error-SD 0.965), ≈ 3 from location (bias −0.019 to −0.032 log-HR), concentrated in the stable-pick third (0.877–0.887 there), where 0.41–0.44 of the complement's naive optimism goes uncorrected, stable across ε and HR (location record). Location + scale explain the whole shape (Gaussian-implied within 0.016 in 48/48 strata). Field retained as the unscaled comparator; naive never reported as a product. Campaign convention from here: `FS_S7_FIELD_SCALEC=selected`. **Package default stays `"none"`** so committed bundles stay byte-reproducible; a default flip is a separate decision.
- **Joint two-subgroup claim:** Bonferroni, γ = 0.025 each side (calibrated = Bonferroni; corr(Λ*, Λ*ᶜ) +0.01 to +0.05); `joint_s` with field-s 0.939–0.952 against 0.95 (`joint` 0.932–0.935). Rule for two-subgroup claims unchanged.
- **Diagnostics (add-beside, defaults off):** `field_decompose = TRUE` → `scale_sel` (= naive robust SE, identity to 4 digits), `scale_win_mean`, `scale_win_cv`, `scale_ratio_c` (ρᶜ), the field's own Var(ζᶜ_G)/Var(m̂ᶜ)/Cov; `return_reselection = TRUE` → p̂(Ĥ). Recorder columns `fld_Hc_scale_*`, `fld_Hc_*_s`, `fld_joint_s_*`.
- **Analysis-time flags:** for field-s, a **high p̂(Ĥ)** (stable pick) is the caution regime on the complement — the reverse of the unscaled field's flag; `se_field_s` should track the naive SE (corr 0.94–0.96). For the unscaled field, λ-SDᶜ/naive SE < 1 flagged the scale deficit.
- **Caveats on record:** field-s inherits the naive SE's calibration to the error — holds in the band regimes (0.99–1.00), fails under a pure-size pick (maxSG: 0.876); the maxSG within-cell redistribution is that inheritance made visible. ε > 0.25 not adoptable. The two-term correction and the field's correction are both driven by re-selection variability while the complement's optimism is driven by the realized winner's extremity; they decouple when the pick is stable — the conditional-vs-unconditional distinction as a regime.

## 3. Results map (survival; records beside the results)

| Campaign / record | Design | Record | Headline |
|---|---|---|---|
| Guo–He | 16 prespecified-family cells, pure argmax | `REPORT_mr_vs_guohe_2026-09-04.md`, `REPORT_mr_field_vs_guohe_2026-09-05.md` | tie constant 0.293; field at G&H's width, within 1–2 points |
| s7 / map1 / s7c / map1c / s7w / map1w | 12.5% prevalence; complement field; refinements | as in the 09-07 handoff §3 | field one-sided 0.92–0.98; complement dominated at 12.5% |
| p30 / p30sg | 31% prevalence; maxeffCons; effMaxSG ε 0.10 J 10 | `REPORT_p30_2026-09-06.md`, `REPORT_p30sg_2026-09-07.md` | identifier returns half the region; complement regime moved (SD/naive 1.08–1.12) |
| nb20 | ε 0.20, J 10 (arm A) and J 20 (arm B) | `REPORT_nb20_2026-09-07.md`, `REPORT_nb20_gate2_2026-09-07.md` | band is the lever, grid is not; harm field 0.956–0.976; complement field 0.897–0.912 (n 500), 0.938 (n 1000) |
| Part A | decomposition on 17 committed cells (no compute) | `REPORT_complement_variance_2026-09-07.md`, `summary_complement_variance.qmd` | error ≈ naive SE (1.00–1.02); field's λ-SDᶜ = family-average scale; marginal SD ≠ error SD |
| Part C | template hygiene | `REPORT_template_hygiene_2026-09-07.md` | pooled meta carries ε and J; tag guard accepts `_` |
| banddial | ε 0.30, maxSG, minSG on the J 10 harm cells | `REPORT_banddial_2026-09-07.md`, `REPORT_banddial_gate2_2026-09-07.md`, `summary_banddial.qmd` | frontier convex (exchange rate 5.8–9.1 → 2.6–2.7 → 0.26–0.36); ε 0.30 first overshoot and first β(Ĥ) dilution; maxSG returns 73–80% of the sample with β(Ĥ) ≈ 1.0; complement coverage recovers to 0.92–0.93 wherever the pick stops varying |
| E0 | instrumented smoke, 200 replicates | `REPORT_field_studentize_stage1_2026-09-08.md`, `REPORT_field_studentize_e0_2026-09-08.md` | ρᶜ predicts the deficit (corr −0.93); corrected ratio 0.996–1.003 in every tertile; s_sel ≡ naive SE |
| E1 | field-s on six cells, identical replicates | `REPORT_field_studentize_e1_stage1_2026-09-08.md`, `REPORT_field_studentize_e1_2026-09-08.md`, `summary_e1stud.qmd`, `e1stud_findings.R` | +0.7 to +1.8 points to 0.912–0.921; ≥ 0.92-with-Wilson bar not met; shape flat on |Ĥ|/|H|, still falling on p̂ |
| Location | means by stratum on e1stud (no compute) | `REPORT_complement_location_2026-09-08.md`, `summary_complement_variance.qmd` §A5 | uncorrected fraction ≈ 0 / 0.13–0.22 / 0.41–0.44 across p̂ tertiles; location + scale explain the shape |
| Chat-side reviews (not in the repo unless Larry routes them) | | `REVIEW_partB_banddial_2026-09-08.md`, `REVIEW_partsAC_close_2026-09-08.md`, `REVIEW_E1_fields_2026-09-08.md`, `MEMO_xuguo2026_reading_2026-09-08.md`, `PROPOSAL_complement_field_scale_2026-09-08_v2.md` | the standard tables and decision trails |

## 4. Decisions on record (2026-09-08)

- effMaxSG band cap: **ε ≤ 0.25** for adoption; working set **ε 0.20 (working band) and ε 0.30 (mapped stress comparator only)**; no ε 0.25 mapping for now; the formal 0.10-vs-0.20 choice remains open with 0.20 as the working setting.
- Field-s adopted as the complement product with its shortfall stated (D-1); variant R1 (canonical studentized form; R1 = global rescale to 1% wherever winner scales are homogeneous, separates only at maxSG).
- Xu & Guo (2026) and Bargagli-Stoffi & Melikechi (2026): not pursued; cite only (positioning sentence in the memo §7).
- Three `R/` additions this week, all add-only with byte-identical defaults: `field_decompose`, `field_scale_complement` (with `joint_s`), and the two one-line pass-throughs in `forestsearch_main.R`; plus the Mac merge's `scale = "identity"` in `R/fs_bias_coverage.R` (reporting path; post-merge regression check = exact reproduction of E1 Finding 2).

## 5. Open items (not proposed; Larry's call)

- **Field-s at n = 1000 and at the null (HR 1.00):** not yet run (E1 was n = 500 harm cells only); the unscaled field's n = 1000 complement record is 0.938.
- **The residual:** a *derived* conditional correction for the stable-pick regime (from what the field already knows about ζᶜ at the observed winner versus the re-selected ones) is the research thread; never a tuned factor. Adjacent: the level dimension for the one-sided guarantee under the screen (null-n1000 0.920).
- Rolling field / field-s into the headline GBSG (maxeff, regenerating family) and ACTG175 analyses.
- Table conventions: whether the standard tables adopt error-scale SD columns permanently (Part A side observation 1).
- Template: an `FS_S7_COMBINED_PATH` knob (Part C side issue); the A3 chunk renders empty HR 1.00 / n 1000 tables on the e1stud set (cosmetic).
- Housekeeping: stale memory-monitor loop pid 2086198 (harmless); the Downloads transport.
- The GLM paths (continuous/MD primary, binary/OR) are the Mac stream's / the GLM program's, not this stream's.

## 6. Cost facts

- Search + gate ≈ 28–40 s/replicate at n = 500 under 100-worker load; field ≈ 15 s; complement field ≈ 2–3 s; **field-s and the decompose diagnostics cost nothing measurable** (fit+MR 36.1 vs 36.3 s light load; complement block 2.0 vs 1.9 s).
- Campaign walls at 100 workers, 2,000 replicates: band cells 31–32 min (n = 500), size-rule cells 23–24 min; n = 1000 cells 77–85 min (J = 20); J = 20 costs ×1.5–2.2 for no gain. Projections from a 5-replicate timing have been accurate to within a few minutes every time.
- Identity/verification renders: 5 replicates ≈ 67–69 s at 5 workers.
- The batch save refuses git-tracked paths (`.refuse_if_tracked()`); pooled metas now carry `effect_neighborhood`, `er_jcuts`, `field_decompose`, `field_scale_complement`.

## 7. Repository state and inventory

- Branch `feature/glm-extension`, HEAD 6d8d36d9 (five commits ahead of origin at the time of writing: task, NOTE, G0 gate bundle, qmd/HTML/report, record fix); Mac merge c0f48a7c in HEAD; forestsearch 0.3.5 installed and matching HEAD by `deparse()` on 661/661 functions.
- Template `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`: knobs `FS_S7_FOCUS` ∈ {maxeffCons, effMaxSG, maxSG, minSG}, `FS_S7_NBHD` (ε), `FS_S7_ER_JCUTS` (J), `FS_S7_Z1Q`, `FS_S7_FIELD_COMPLEMENT`, `FS_S7_IJ_RESIDUAL`, `FS_S7_FB`, `FS_S7_FIELD_DECOMP`, `FS_S7_FIELD_SCALEC`, `FS_S7_CAMPAIGN` (alphanumerics and `_`), `FS_S7_NSIMS`, `FS_S7_START`, `FS_S7_WORKERS`, `FS_S7_MODE=combine`, `FS_S7_SAVE_COMBINED`.
- Documents: `summary_complement_variance.qmd` (A0–A5; any bundle set via `FS_SUMCV_GLOBS`), `summary_banddial.qmd`, `summary_e1stud.qmd`, `e1stud_findings.R`, `summary_bias_coverage.qmd`.
- Gate bundles committed beside the campaigns: `stud1inert`, `stud1decomp`, `e1inert`, `e1scale`, `e1regr`, `postmerge`, `e0stud` (200 replicates).
