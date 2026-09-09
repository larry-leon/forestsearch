# TASK — C-4: fixed-p̂ diagnostics across the thirteen committed cells, the harm-block two-sided decay, and detection-conditioning (zero compute)

Date: 2026-09-09. Author: chat (spec). Executor: Claude Code (Linux). Approver: Larry (C-1…C-4 agreed 2026-09-09). Reviewer: the Linux MR-field chat.
Predecessors: `REPORT_cert20_2026-09-08.md` (31% prevalence, sample-size profile), `REPORT_tier2_2026-09-08.md` (12.4%, dominated regime), `REPORT_field_studentize_e1_2026-09-08.md`, `REPORT_harm_location_2026-09-08.md`, `REPORT_complement_location_2026-09-08.md`, `summary_complement_variance.qmd` §A5/§A6 (the means-by-stratum machinery), `dev/notes/NOTE_complement_product_2026-09-08.md`.

**Framing, to be stated in the record's opening paragraph.** p̂(Ĥ) is a **recorded diagnostic**, computed from draws already made after the bounds are formed; **no construction in this package reads it, and nothing in this task proposes that any should.** The fixed-p̂ bands below are a re-cut of committed simulation results for analysis only — the same replicates, grouped by a stratum definition that means the same thing at every n and prevalence, because the per-cell tertiles used so far shift their boundaries with n (12.4% T3 is [0.268, 0.991] at n = 500 and [0.594, 0.996] at n = 1500) and therefore cannot answer an n-comparison. Reading p̂ as an analysis-time flag (high p̂ ⇒ the stable-pick regime, bounds somewhat optimistic) is the existing documented rule and is unchanged. Note in the record that p̂ is estimated from the same draws that build the bound, so it is not an independent instrument — adequate for a flag and for diagnostics, and a reason not to let any interval depend on it.

## Protocol

- First action: archive stale variants in `~/Downloads` to `~/Downloads/cc_archive/`; copy this file to `dev/tasks/` and commit. Do not push. If the file is missing from `~/Downloads`, reconstruct it from the kickoff (which carries the band edges and every section), commit it under this name with a reconstruction note, and proceed.
- **No compute. No `R/` change. No bundle written** (assert at the end: `git status` shows no new `.rds` under `results/`; the seven pre-existing untracked files are not this task's and must not be committed).
- Transplant-first: build the new sections by copying the committed §A5/§A6 chunk pattern and the `cert20` / `tier2` findings scripts — no fresh authorship of loading, filtering or table machinery.
- Standing conventions: winner-only and winner-floor excluded; bounds by location; Wilson intervals; marginal and error SDs side by side; verify from source. **No repair proposal, no recommendation change; report and wait.**

## Data: the thirteen committed cells

`e1stud` ε 0.20 HR 1.50 / 1.75 n500 (31%); `cert20` HR 1.50 / 1.75 at n1000 and n1500, HR 1.00 at n500 / n1000 / n1500 (31%); `tier2` HR 1.75 at n500 / n1000 / n1500 and HR 1.00 n500 (12.4%). State per cell: campaign, prevalence, focus (`effMaxSG` ε 0.20 vs `maxeffCons`), HR, n, detections. Cross-machine note: `tier2` was produced on the Mac — comparisons across campaigns are of coverage rates and summary statistics, never `identical()` on fitted values.

## Part F — Coverage at fixed p̂ bands (both blocks)

Bands (fixed, common to every cell): **[0, 0.05), [0.05, 0.10), [0.10, 0.20), [0.20, 0.35), [0.35, 0.55), [0.55, 1.0]**. Report the band's n per cell; suppress any band with fewer than 30 replicates from the trend reading (state it, do not drop it from the table).

- **F1.** Per cell × band: n, mean p̂, harm field one-sided lower coverage [Wilson], complement field-s one-sided upper [Wilson], IJ two-term two-sided on the harm block [Wilson], and the Gaussian-implied value beside each from the band's own mean and SD of the relevant error and its mean SE.
- **F2.** Pooled across cells **within prevalence and HR**, by band: the same quantities, so the shape of coverage as a function of p̂ is visible with usable n.
- **F3.** The n-comparison at fixed band: for each band, coverage against n (500 / 1000 / 1500) within each prevalence and HR. **This is the question the per-cell tertiles could not answer:** at a fixed p̂ band, does coverage improve, stay flat, or degrade with n?
- **F4.** Composition: per cell, the share of replicates in each band, so the movement of the p̂ distribution with n and prevalence is on record beside the coverage-versus-band shape.
- **F5 — Reading (in the record, not a task).** State which the numbers support: (A) coverage is a **stable function of p̂** — at fixed band it is flat in n, and the cell-average pattern is explained by the composition shift in F4; (B) coverage **degrades at fixed band** as n grows — something is worsening beyond composition; (C) mixed (say where). Do this separately for the harm block, the complement block, and IJ two-sided.

## Part G — The harm-block two-sided decay

The finding: IJ two-term two-sided on the harm block runs 0.9611 → 0.9166 → 0.9129 at 12.4% prevalence (HR 1.75, n 500 → 1000 → 1500) against 0.971–0.981 at every 31% cell; the unscaled field's two-sided moves the same way (0.887 → 0.850 → 0.878). The complement block is 0.9995–1.0000 throughout. The hypothesis on record is a persistent harm-block bias (field −0.075 / −0.091 / −0.075; IJ −0.010 / −0.081 / −0.093 log units) consuming a two-sided interval as the SE shrinks.

- **G1.** Per cell (all thirteen), harm block, for the IJ two-term and the field: the two-sided **miss rate split by side** — share of replicates whose interval lies entirely above β(Ĥ) and share entirely below — beside the bias (log), the marginal and error SDs, the mean SE, and SE/error-SD. If the mechanism holds, the miss is one-sided and grows with n at 12.4%.
- **G2.** The trajectory: bias, SE, bias/SE, and the two-sided miss by side against n, within each prevalence, one table per block. State whether bias/SE grows with n (the signature of a bias that does not shrink while the SE does).
- **G3.** Cross with Part F: the two-sided miss by side, by fixed p̂ band — does the decay live in the stable-pick bands, or is it uniform across p̂? This is the test of whether the two-sided decay and the stable-pick caveat are the same phenomenon.
- **G4 — Reading (in the record, not a task).** Whether G1–G3 support the persistent-bias mechanism, and whether the decay is stratum-localized or uniform. No repair, no recommendation on `ci_method`.

## Part H — Detection conditioning (Larry's question, made testable)

All coverage is computed on detected replicates, and detection rates move with n (12.4% harm: 0.950 / 0.995 / 0.999; 31% null: 0.920 / 0.955 / 0.959; 12.4% null n500: 0.680), so the analysed set's composition changes with n independently of anything the constructions do.

- **H1.** Per cell: detection rate; and for the cells sharing DGM draws across n within a prevalence and HR, the set of sim_ids detected at **every** n in that series ("always-detected"). Report its size.
- **H2.** For each such series: coverage (both blocks, field / field-s / IJ two-sided) and the harm bias computed (i) on all detected replicates, as reported so far, and (ii) restricted to the always-detected set — the like-for-like comparison across n. State whether the n-pattern changes when composition is held fixed.
- **H3.** The undetected replicates' character: for cells with meaningful non-detection (12.4% null n500, 31% nulls, 12.4% HR 1.75 n500), compare `n_true`, the naive harm effect and its SE between detected and undetected replicates, to say what the screen is removing.
- **H4 — Reading (in the record, not a task).** Whether detection conditioning contributes to the harm-block bias trajectory and to the two-sided decay, or is second-order beside the p̂ composition shift of F4.

**Documents.** Extend `summary_complement_variance.qmd` with sections **A7 (fixed-p̂ bands, Part F)**, **A8 (two-sided decay, Part G)** and **A9 (detection conditioning, Part H)**, transplanted from the A5/A6 chunk pattern; the bundle set comes in through `FS_SUMCV_GLOBS` (all thirteen cells) with a cell-metadata table naming campaign, prevalence, focus, HR, n. Render. Output `REPORT_fixedphat_ij2s_2026-09-09.md` beside the `cert20` records with F1–F5, G1–G4, H1–H4, every number verbatim from the rendered document, and the framing paragraph above at the top.

## Done means

The three sections added and rendered; the record committed; no bundle written; the seven pre-existing untracked files left untouched; branch left unpushed; one-paragraph closing summary giving the F5, G4 and H4 readings in one line each and the commit range. Out of scope: any repair, any compute, any `R/` change, any recommendation on `ci_method` or on the documented rules.
