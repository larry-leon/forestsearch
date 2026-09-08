# HANDOFF — Mac continuous/MD workstream: the field method on the continuous path

Date: 2026-09-07. From: the CI-construction chat (2026-09-04 → 09-07). For: a fresh chat that will write specs and review CC reports for the Mac Studio work stream. This document is self-contained; every record it cites is in the repository on `origin/feature/glm-extension` (pushed 2026-09-07 at 22c5f854) and readable by CC.

## 1. Roles and protocol (unchanged)

- Larry decides and approves; the chat writes specs (task documents) and reviews; Claude Code (CC) does every repository operation. CC never pushes; Larry pushes via GitHub Desktop.
- Task documents travel as `.md` files via ~/Downloads; CC's first action is to copy the document into `dev/tasks/` and commit it. Kickoff messages are short, pasteable, and point CC at the file; always give the kickoff inline as cut-and-paste text.
- Gates are stop-on-failure; Gate 1 is the compute go/no-go with cell count, replicates and wall-clock stated; compute can be pre-authorized in the kickoff for unattended runs with a wall ceiling and a hard timeout.
- Verify from source, never from description or memory; CC quotes paths, signatures and line numbers in Stage 0 records. Records (`REPORT_*`) go beside the results, not in `dev/tasks/`.
- Any R/ change is called out and classified: adds code / changes behaviour / changes the method. Method changes are proposals; add-only with byte-identical defaults is the standard. Do not re-run committed work; identity checks against committed bundles are the pairing proof.
- Findings go in the record; a task is proposed only when something blocks. Do not attach a follow-up task to every finding.
- Analysis documents must be self-contained (no reads from `dev/` or simulation directories), with the interpretation inside the document.
- Interpretation convention: read bounds by location against clinically meaningful effect sizes, never as significance at the null; report the spread across adjusted bounds as the price of selection; coverage tables carry Wilson intervals; bias in SD units beside the natural scale.
- The paper is not a topic in this work stream (it is under review); the work is CI construction only.

## 2. What the constructions are (as they stand)

Selection: forestsearch identifies a harm subgroup Ĥ from an enumerated family of covariate-cut conjunctions; the reported effect is the standard within-subgroup analysis (Cox coefficient on survival; mean difference on the continuous path). The target throughout is β(Ĥ), the true effect in the region actually found (and β(Ĥᶜ) for its complement), with the family held fixed.

- **MR (IJ two-term)** — the paper's method: multiplier resampling over the fixed family, two-term de-biased estimate β̃ = β̂ − bias_sel − bias_fix, infinitesimal-jackknife variance from the two-term residual. Conservative by construction (SE/SD 1.3–2.0 in most regimes; ≈ 2 when one candidate dominates); the reported two-sided interval by default (`ci_method = "ij"`), with one named failure regime (moderate harm at n = 1,500, two-sided 0.90).
- **MR (field)** — `ci_method = "field"`, add-only since 2026-09-05: simulate the de-biased estimator's error on a "shrunk field" (the observed candidate effects with the winner's entry replaced by β̃), Λ* = ζ*_{G(w+ζ*)} − m̂(w+ζ*), with Gaussian multipliers through the influence matrix; returns est₂ = β̃ − mean(Λ*), the one-sided bound β̃ − q₀.₉₅(Λ*), the plain two-sided quantile interval, λ-SD. Its **one-sided bound is the standard directional product** (harm: lower bound; complement: upper bound). Its two-sided interval under-covers in harm regimes and is not a product.
- **Complement field** — `field_complement = TRUE` (2026-09-06): the complement's Λ*ᶜ from the same multipliers and winners; primary product the one-sided **upper** bound (benefit claim). Costs ~3 s/replicate.
- **Joint pair** — `field$joint`: simultaneous (Ĥ lower, Ĥᶜ upper) from the aligned draws; the calibrated pair equals Bonferroni (γ = 0.025 each) because the two errors are independent to within ±0.05. Rule for a claim on both subgroups: Bonferroni.
- **Diagnostics** — `return_reselection = TRUE` gives p̂(Ĥ) (the winner's re-selection frequency: near 1 = settled selection, below ~0.5 = tie regime); for the complement the regime diagnostics SD(β̃ᶜ)/naive SE and λ-SDᶜ/naive SE.
- **Display** — `fs_sim_bias_coverage()` / `fs_plot_bias_coverage()`: coverage against residual bias (SD units) with a Gaussian reference Φ(1.645·r ∓ b) / Φ(1.96·r − b) − Φ(−1.96·r − b) from each cell's own bias, SD and mean SE; `side = "lower"/"upper"`. Needs a `scale = "identity"` option for the MD path (the Mac task's one R/ change).
- **Closed lines** (do not re-open or re-propose): the κ-calibrated two-sided interval (`field_uniform`; valid but wider than IJ and full-bootstrap-class cost — a documented research option); the hybrid κ; the IJ winner-only variance (rejected 2026-09-07 — under-covers wherever selection isn't settled; this under-performance is what motivated the two-term approach) and the winner-floor variant (rejected on theoretical grounds — an unadjusted naive SE cannot be promoted from simulation performance); covariate adjustment or alternative working models (not on the table: the method must apply to standard drug-development analyses); the vintage-flip attribution task (dropped; not to be overstated from current experience).

## 3. Survival results to compare against (all cells 2,000 replicates; records in `quarto/simulations/gbsg_020/`)

| Block / product | Range across the 12.5%-prevalence cells (s7, map1, s7c) | At 31% prevalence, effMaxSG (p30sg) |
|---|---|---|
| Ĥ field one-sided lower coverage | 0.920–0.981 (one soft cell: tie regime at n = 1,000, 0.920 — the screen-conditioning effect) | 0.946–0.971 |
| Ĥ IJ two-sided coverage / SE-SD | 0.90–0.99 / 1.08–1.96 | 0.977–0.991 / 1.34–1.67 |
| Ĥ field two-sided | 0.85–0.92 in harm cells (over-correction tail) | 0.932–0.985 |
| Ĥᶜ field one-sided upper coverage / λ-SDᶜ-to-naive | 0.930–0.937 at n = 500, 0.953–0.956 at n ≥ 1,000 / 0.97–0.99 | 0.911–0.932 (regime moved: SD(β̃ᶜ)/naive 1.08–1.12) |
| Ĥᶜ IJ two-term SE/SD | 1.80–1.85 (upper bound ≈ 1.0 at n = 500, uninformative) | 1.65–1.85 |
| Joint Bonferroni | 0.940–0.960 | 0.938–0.942 |
| Retained bias, Ĥ: naive → IJ → field (SD units) | +2–6 → +0.1–1.6 → −0.4–0.9 | +2.1–3.7 → +0.1–0.7 → −0.1–0.3 |
| Complement display points | all on the Gaussian reference (no tail) | — |

Guo–He comparison (prespecified families, pure argmax, `quarto/GuoHe/`): field at G&H's width and within 1–2 points of its calibration; retains 8–13% of the tie optimism vs G&H's 0–3% (IJ 26–31%); the tie constant 1 − 2^(−1/2) = 0.293 confirmed to three decimals on disjoint ties.

Key records: `REPORT_mr_field_s7_2026-09-05.md`, `REPORT_mr_field_ocmap_2026-09-05.md`, `REPORT_mr_field_complement_2026-09-06.md`, `REPORT_complement_refinements_2026-09-06.md`, `REPORT_p30_2026-09-06.md`, `REPORT_p30sg_2026-09-07.md`, `REPORT_bias_coverage_display_2026-09-06.md`, `SUMMARY_ij_vs_field_2026-09-05.md`, and `quarto/applications/gbsg/analysis_gbsg_survival_frozen_family.qmd` (the real-data illustration: adjusted lower bounds 0.79–0.96 against naive 1.30 on the harm region; complement field upper 0.79 vs IJ 0.94; joint pair (0.76, 0.82)).

## 4. What the continuous/MD path has and lacks

Has (committed, pushed): the continuous twin template and bundles from the MD workstream (`HANDOFF_continuous_2026-08-27_v5`), the OC wrapper and the applied OC evaluation on ACTG175 continuous (`analysis_actg175_continuous_oc_evaluation.qmd`), MR (IJ) intervals. Lacks: everything in §2 after the first bullet — the field block, complement field, joint pair, p̂ recording, display, campaign-tag stems, `FS_*` knobs, the `.refuse_if_tracked()` save guard. The gate is outcome-agnostic (effect vector + influence matrix), so the field should run unchanged; the continuous path is where the theory is cleanest (linearization closest to exact for a mean difference, Λ* closest to Gaussian), which is why it is the right validation. Binary/OR port status: opened 08-30/31 in a separate chat; not in this record — check the repo before assuming.

## 5. The Mac task

`TASK_continuous_field_mac_2026-09-07.md` (in ~/Downloads on the Mac; copy to `dev/tasks/` first). Branch `feature/glm-extension-mac` from `origin/feature/glm-extension` at 22c5f854; push that branch; never push to `feature/glm-extension`. The Linux box is concurrently running `TASK_p30sg_nb20_2026-09-07.md` on `feature/glm-extension` (effMaxSG band ε = 0.20; arm A two harm cells on the J = 10 er grid, arm B five cells on J = 20 incl. harm 1.5 at n = 1,000; ~6 h). The Mac task must not touch: `quarto/simulations/gbsg_020/*`, `R/fs_mr_inference.R`, `R/forestsearch_main.R`. Its one permitted R/ change: `scale = c("log","identity")` on `fs_sim_bias_coverage()` in `R/fs_bias_coverage.R`, add-only, the 14-point survival fixture still passing. Attended through Gate 1 (Mac worker count and projection unknown until measured); pre-authorization from Gate 1 with a ceiling set from the measurement. Merge back: `git merge feature/glm-extension-mac` on Linux when both sides are quiet — trivial because the task adds files only.

Kickoff (after `git checkout -b feature/glm-extension-mac origin/feature/glm-extension`, package installed, file in ~/Downloads):

```
Task (Mac Studio): continuous/MD path — port the field additions to the continuous twin and evaluate the current constructions.
Task document: ~/Downloads/TASK_continuous_field_mac_2026-09-07.md
Branch: feature/glm-extension-mac (already checked out from origin/feature/glm-extension at 22c5f854); push this branch as commits land; never push to feature/glm-extension.
First action: copy the task document to dev/tasks/ and commit. The "Standing conventions" section of the document governs this session.
R/ change acknowledged: a scale = c("log","identity") argument on fs_sim_bias_coverage() only, add-only, default preserving the 14-point fixture; nothing else under R/, and no file listed as off-limits in the document.
Decisions M-1–M-4 at defaults. Execute Stage 0 → Stage 1; stop at Gate 1 with the identities and the measured Mac projection. M-5 decided there.
```

## 6. Things the reviewing chat should do

- At Gate 0: check the cells CC lists (M-1) and the harm orientation on the MD scale; confirm the identity anchors (which pre-existing columns the committed bundles carry).
- At Gate 1: read the worker calibration and the per-cell projection; set the pre-authorization ceiling from them (a Mac Studio is several times slower per cell than the 100-worker Linux box).
- At Stage 3: read the field's one-sided coverage on both blocks against §3's survival ranges; expect the smallest Gaussian-reference departures of any path; report bias on the MD scale and in SD units; bound locations against stated MD thresholds; the complement's regime diagnostics. Correct any "below/above the null" framing to bound-location language.
- Do not propose: κ variants, hybrid κ, winner-only or winner-floor rows, covariate adjustment, attribution of cross-version selection flips.
- Cost facts to reuse: the field pass ≈ 15 s/replicate on survival n = 500 (complement ≈ 3 s); dense-matrix Monte Carlo (κ, M_eff) is memory-bandwidth-bound and scales badly with workers; the search itself scales well.

## 7. Open items elsewhere (for awareness, not for this stream)

The Linux nb20 campaign (identifier band/grid); the level-dimension protection for the one-sided guarantee under the screen (research, unopened); the planted-boundary-on-grid stress design (noted, not proposed); rolling the new bounds into the headline GBSG and ACTG175 analyses; the binary/OR path.
