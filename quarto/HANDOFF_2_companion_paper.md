# Hand-off — Companion Identification Paper (Calibrated Declaration + Pareto Frontier)

*Paste into a new chat to start companion-paper work with full context. Self-contained.*

---

## Who / what

Larry F. León, Merck biostatistician, author of `forestsearch` (CRAN v0.1.0; GLM extension → v0.2.0). Co-author Keaven M. Anderson. The **submitted** JASA paper is "Post-Selection Inference for the Standard Trial Analysis of Data-Driven Subgroups" (León & Anderson) — it corrects winner's-curse bias in data-driven-subgroup effect estimates via refit-free multiplier resampling (MR) yielding infinitesimal-jackknife (IJ) intervals. Prior work: León et al. (2024) "Exploratory subgroup identification in the heterogeneous Cox model," *Stat Med* 43(20):3921–3942.

**This chat is the COMPANION paper** — a *new* paper, distinct from the submitted one.

**Notation discipline (hard):** FS / DINA / GRF are the three identifiers. β(Ĥ) = conditional estimand (Theorem-2 target), kept STRICTLY distinct from θ†(H) = marginal truth. MR = multiplier resampling (refit-free); FB = full bootstrap. Never conflate log-HR (β) and HR (ratio) scales. Prefer "consistent with"/"corroborated by" over stronger verbs for simulation findings. FS uses a fixed candidate family; DINA/GRF regenerate families per replicate.

**Working rules:** read files before editing; byte-level edits with count assertions; never fabricate values or citations (six .bib entries carry `% VERIFY:` markers — do not invent bibliographic details); every file-modifying response ends with a downloadable artifact; formal academic register, LaTeX math in `$...$`; recommend before implementing; ground claims in actual source.

## The thesis (verified against the v5_6 supplement)

**The splitting-consistency screen IS stability selection whose stability probability has a CLOSED FORM; therefore the family-wise false-declaration rate can be CALIBRATED EXACTLY, not merely bounded.** This speaks to two literatures — subgroup identification and stability selection. The charter is v5_6's own closing disclaimer that the correction "does not control the rate of false declaration," plus León (2024)'s literal promise to evaluate "the Guo et al bootstrap calibration procedure" and its claim that type-1 error is "quickly approximated by numerical integration" — the companion replaces that heuristic with an exact device.

### The three-line reduction (proved in the v5_6 supplement)
1. (§2.3) fair-coin split = Rademacher multiplier bootstrap; the two halves are a mirror pair β̂ ± D_g, with Var(D_g) = σ²_{D,g}.
2. "both halves pass" ⟺ |D_g| ≤ β̂(g) − c_cons → closed form `p_cons(g) = 2Φ((β̂(g) − c_cons)/σ_{D,g}) − 1`.
3. (S3.2) the screen ⟺ an exact one-sided z-test. This is PER-CANDIDATE; FS declares if ANY candidate survives → need **max-over-candidates** control.
- Device: **Minnier-Tian-Cai (2012, JASA 106(496):1371–1382)** simultaneous region — perturb → standardize → max → (1−α) quantile. The v5_6 supplement already computes `P = B_effᵀ Ξ` (an S×draws matrix), which is exactly this null law.
- CCK (Chernozhukov-Chetverikov-Kato 2015) Thms 2–3 are already invoked for the re-selection law; the non-asymptotic rate is DEFERRED there → the companion picks it up. **Shao (1996, JASA 91(434):655–665)**: n-out-of-n bootstrap model selection is INCONSISTENT → the multiplier route is NECESSARY (already cited in the package's `compute_frontier_cis()`).

## Prototype status — VALIDATED

A max-statistic calibration prototype (`fs_dg_calibrate.R` in outputs) was unit-tested: κ (the family-wise cutoff = (1−α) quantile of max_g D_g/σ_{D,g}) matches the brute-force quantile to machine precision (2.795), and κ > z always. Null Monte Carlo (self-consistent null, α=0.05, 1000 reps, three overlap regimes): the calibrated screen controls the false-declaration rate at 0.05 across **all** overlap regimes — disjoint (cor 0.00) 0.051, random (cor 0.40) 0.049, nested (cor 0.82) 0.048 — versus ~0.78 uncalibrated. **The nested/correlated case (cor 0.82) is the technical crux and it held.** A methodological wrinkle to note in the paper's sim design: the observed statistic must be drawn from the true null law (Rademacher through the same influence matrix), not Gaussian; a mismatch inflates the rate and does NOT shrink with more draws (the tell). In the real sweep this is automatic (observed β̂ and multiplier draws share the same B).

## Literature map (all read from PDFs; roles established)

- **Minnier-Tian-Cai (2012)** — perturbation root + max-statistic device (the core tool).
- **Shao (1996)** — bootstrap model selection inconsistent; multiplier necessity.
- **Guo & He (2021, JASA)** — debiasing/CI anchor, tuning r ∈ (0, ½).
- **Xu & Guo (2026, arXiv 2605.03141)** — Xinzhou Guo 2nd author; keeps the fitted score FIXED; §3.1 declares repeated-split / subsampling / m-out-of-n INVALID when the object regenerates. This is **independent corroboration of Larry's regenerating-family argument, from Guo himself.** Analyzes ACTG175 (continuous CD4). Demolishes sample-splitting for inference ("loses power… unstable, inviting selective reporting") — the contrast to Chiu's splitting approach.
- **Zhao-Ivanova-Fine (2023, Clin Trials 20(4):370–379)** — measured that augmented/re-modeling gives "less bias but poorer coverage," theory "beyond scope" = the anomaly the conditional estimand explains. **LEAD MOTIVATION.** (bib page fixed 394→370.)
- **Zhao-Fine-Ivanova (2024, Stat Med 43:2487–2500)** — multi-outcome best subgroup.
- **Zhao-Tian-Cai-Claggett-Wei (2013, JASA 108(502):527–539)** — perturbation for subgroup selection; explicitly PUNTS on the search step; ancestor of MR.
- **Shah-Samworth (2013, JRSS-B 75(1):55–80)** — complementary pairs stability selection + error control; Thm 1(i) is WORST-CASE-IN-EXPECTATION; the companion's smooth AND-both-pass gives a CLOSED FORM → EXACT. **This is the template to sharpen.**
- **Meinshausen-Bühlmann (2010, JRSS-B 72(4):417–473)** — origin of stability selection.
- **Isotonic Subgroup Selection (Müller-Reeve-Cannings-Samworth, arXiv 2305.04852, `ISS` pkg)** — IS a declaration method with a non-asymptotic uniform Type-I guarantee + minimax power (a STRONGER guarantee than the companion claims, so do not overclaim). BUT: verified **0 occurrences of censor/survival/hazard/Cox/Kaplan** in the paper — they binarize ACTG175's censored composite; `forestsearch` is a SURVIVAL method → **decisive differentiation via censoring**. Also requires a prespecified monotone covariate direction (indefensible in exploratory search). Monotonicity (shape constraint on TRUTH) vs the companion's smoothness (regularity on the ESTIMATOR) are ORTHOGONAL.
- **Müller et al. (Feb 2026, arXiv 2512.15676)** — Samworth+Novartis "Data-driven controlled subgroup selection," direct competition — track it.

### Pareto-frontier lineage (the second contribution)
- **Wang & Rudin (2022, INFORMS J Computing 34(3):1626–1643)** "Causal rule sets…" — ORIGIN of the "efficient treatment-frontier" idea.
- **Chiu (2025, arXiv 2507.09494)** §3.3 Pareto Efficient Frontier over rule sets — gets inference by SAMPLE SPLITTING (the contrast: Xu & Guo demolish splitting).
- CAPITAL / Muysers — display only.
- **NOVEL CONTRIBUTION: post-selection INFERENCE on a Pareto frontier — nobody has done valid in-sample inference on it.** MR-compatibility on the frontier has been VERIFIED against the package source (the bootstrap re-runs the full search per replicate at `bootstrap_analysis_dofuture.R:563`, confirming the Shao regime and that MR is the stable conditional alternative). Two open research items: (a) simultaneous bands over the frontier via the MTC γ_α device; (b) frontier-membership stability.

Reserve/tiered references (roles noted in `references.bib`, 77 entries): Knaus (2022), Liu-Mazumder-Radchenko (2026), Liu-Wang (2018), Sun-Sechidis-Bornkamp (2024, `benchtm`), Lipkovich (2017), Stallard (2024), Andrews-Kitagawa-McCloskey (2024), Wang-He (2023).

## Codebase facts (verified — reuse, don't re-derive)

- `effect_neighborhood`: KEEP IT. `.compute_inclusion_band()` (subgroup_consistency_helpers.R:736) is a pure functional of {β̂(g)}; sort keys all fixed/functional. Satisfies fixed-family + functional-of-effects for all three selection rules. `"pareto"` alone gives a zero-tuning-parameter path.
- Closed-form Pcons ALREADY IMPLEMENTED: `consistency_resample.R:342` `out$rate_closed <- max(0, 2*pnorm(delta/sigma_D)-1)` = exactly S3.2's p_cons, with a validator (`consistency_resample_compare()`) checking vs literal splits.
- `args_call_all` on the forestsearch return (forestsearch_main.R:996) = complete replayable arg-list; c1 = hr.threshold, c2 = hr.consistency.
- `compute_frontier_cis()` (frontier_cis.R) returns 3 POINTWISE CIs (naive / split-citing-Shao / FSBC); the simultaneous-band gap over the frontier is the open contribution.

## Datasets / DGM

GBSG (survival/Cox) and ACTG175 (binary/GLM/odds-ratio). DGM truth: θ†(Hᶜ) = 0.6569 (≈ 0.66). `setup_gbsg_dgm()` is a survreg AFT; harm subgroup H = {z1=1 AND z3=1}; `k_inter=0` gives true homogeneity.

## Deliverables already built (in outputs, reuse)

- `companion_paper_plan.md` (v3; v1/v2 archived) — full plan: thesis/two-audiences, theorem-reachable reduction (3-line + MTC device + Shao necessity + CCK crux), architecture (3 families + Xu&Guo convergence + ISS competitor with censoring gap), Pareto-frontier part (citation debt + MR-compat verified + simultaneous bands + frontier stability), three-way ACTG175 table, assets, deferred package work.
- `companion_reference_list.md` — tiered acquisition list.
- `references.bib` — 77 entries, parse-verified, six with `% VERIFY:` markers.
- `fs_dg_calibrate.R` — the validated prototype (`.fs_dg_calibrate()` + `fs_declare_calibrate()`).
- `SMOKE_TEST_RESULTS.md` — the validation write-up.

## The empirical case just got STRONGER (new, important)

The submitted paper's FDR numbers were biased LOW by a now-fixed bug (intermittent subscript crash in `subgroup.consistency` that silently suppressed declarations). Corrected borderline-null false-discovery rates: **H0-hom 0.138** [0.118, 0.161], **H0-harm 0.278** [0.251, 0.307] — roughly DOUBLE the previously-believed 0.10–0.17. Larry deems ~10–15% tolerable for exploratory identification but **>30% an issue**; H0-harm at 0.278 sits at the edge. **This materially strengthens the case that the calibration is needed, not optional** — the earlier "calibration is polish" read was based on the biased-low numbers. This is the empirical hook for the companion paper's motivation. (The bug is fixed and merged to `feature/glm-extension`; grids verified unaffected because they use the parallel path.)

## Open research work (the paper's contributions to complete)

1. Prove the CCK coupling for the finite-family max-statistic (the non-asymptotic rate the v5_6 supplement deferred).
2. Run the max-statistic calibration prototype at scale on knoise0new on the 127-core box; confirm it reproduces the idealized 0.05.
3. Simultaneous frontier bands via the MTC γ_α device; frontier-membership stability.
4. Sharpen the Shah-Samworth worst-case-in-expectation bound to the exact closed form.
5. Three-way ACTG175 identification table (FS/DINA/GRF).
6. Acquire the reference stragglers (the six `% VERIFY:` entries; Reeve/Cannings/Samworth; Wang & Rudin 2022; arXiv 2507.09494; Muysers/CAPITAL details).

## First action for the new chat

Pick up from `companion_paper_plan.md` v3. The natural next step is either (a) drafting the calibration section around the validated prototype + the strengthened empirical motivation (the corrected ~0.28 FDR), or (b) formalizing the CCK coupling. Ask Larry which. Do NOT re-derive the reduction or re-establish the literature — both are done above.
