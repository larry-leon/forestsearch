# REVIEW — `grfmr` campaign (GRF MR inference), ten harm cells

- **Date:** 2026-09-12
- **Reviews:** CC's `grfmr` report and per-cell extraction (`TABLES_grfmr_percell_2026-09-12.md`), commit range `a451230d..5e3f71a6`.
- **Campaign:** 10 of 12 harm cells at 2,000 replicates — 12.4% and 31% × HR 1.50 × n 500/1000/1500, 12.4% HR 1.75 × n 500/1000/1500, 31% HR 1.75 n 500. Deferred: 31% HR 1.75 at n 1000 and n 1500. Dropped: none. 7.939 h against a 9 h ceiling, ratio 0.896 uniform across cells.
- **Captions as run:** every coverage number is labelled coverage of β(Ĥ) conditional on the proposed family, over detected replicates. Section 2 is why that label may be too weak, and why nothing should be changed until it is settled.

---

## 1. Verification — accepted

- Gate 2 PASS on all ten cells; Gate 3 7/7 on all twenty batches.
- Corrected identity 2.22e-16 to 3.33e-16; Bonferroni identity exactly 0 on the γ-at-floor rows; γ ∈ [0.02500, 0.02800].
- Amendment 3 holds everywhere: `n_true` `identical()` on all 2,000 rows of every cell, truth `all.equal()` at 1e-8 at every cell. The three 12.4% HR 1.75 cells are bit-identical; the rest differ by 6.66e-16 to 8.88e-15. No DGM-path finding.
- Part T2 add-only, Gate T2 7/7 with 163 non-timing columns `identical()`, and the column confirmed populated by smoke rather than declared.
- Detection 1.0000 at all four 31% cells and 0.9975–1.0000 at 12.4%, flat in n. **No n-trend in this campaign carries a detection-conditioning caveat** — the qualification that attached to every DINA n-trend is absent here.
- `admitted_n` never 0 across 20,000 replicates, and never 1. The empty-admitted-set path remains unexercised, consistent with the probes.

### 1.1 A defect CC found and fixed mid-extraction

- The `strat-xtab` chunk built its own tertiles on `n_family` and does not call `strata_of()`, so the transplant's stratifier substitution missed it. Cross-tabulating p̂ against the outcome-independent enumerated pool would have been uninformative.
- Fixed to `admitted_n`, summary re-rendered, committed `5e3f71a6`. All numbers below are post-fix and are produced by re-executing the `.qmd`'s own definition chunks, so they are the rendered document's values.
- Worth noting as a pattern: the substitution was applied to the helper, and one chunk bypassed the helper. Any future stratifier change should grep for direct `tert*()` calls, not only for `strata_of()`.

---

## 2. Gate 0a — the determination that may change what this campaign establishes

Three source findings, taken together:

- The enumerated candidate pool comes from quantiles of X subject to `n_min`, with DR scores entering only the effect column (`R/grf_subgroup_labels.R:255–277`). `n_min` resolves from sample size alone.
- MR re-evaluates **the full enumerated pool**, not the forest-admitted subset: the unfiltered pool is what is attached to the result (`grf_main.R:331`, `:419`); the admitted copy is local and discarded (`forestsearch_helpers.R:1625–1627`).
- Admission is applied **per draw**, as a filter on perturbed effects (`fs_mr_inference.R:661–666`).

**What that is.** A fixed, finite candidate list, fixed before any outcome is seen, re-evaluated in full on every draw, with selection re-run inside the draw. That is the structure the fixed-family condition asks for.

**Why it matters.** Handoff §3 calls GRF "the sharpest case of the family caveat" on the grounds that a literal bootstrap re-fitting the forest would propose a different family each resample. MR does not re-fit the forest. If the reading holds, the reasoning does not reach the MR construction, and GRF's numbers would not be merely conditional on the proposed family.

**The numbers agree with the source reading**, which is what makes it worth taking seriously rather than filing as a curiosity:

| | GRF | DINA (harm cells) | FS (certified) |
|---|---|---|---|
| Field lower on β(Ĥ) | 0.930–0.969 | 0.902–0.944 | 0.941–0.980 |
| Field retained bias, log scale | \|bias\| ≤ 0.049 at every cell | +0.05 to +0.48 error-SD | over-corrects at low p̂ |
| IJ two-sided at 12.4% | mild decay, 0.989 → 0.980 | flat and high, 0.989–0.995 | strong decay, 0.977 → 0.901 |
| IJ miss side where it decays | entirely upper | balanced, ≤ 0.009 | entirely upper |

- GRF overlaps FS's certified band and sits above DINA's at every cell.
- GRF's field estimate is nearly unbiased — \|bias_log\| ≤ 0.049 across all ten cells, against DINA's large positive retained bias.
- **GRF reproduces FS's qualitative signature, not DINA's**: a mild two-sided decay with n that is entirely upper-limit misses (12.4% HR 1.75: above 0.0035 → 0.0145 → 0.0200 while below stays 0.0065 → 0.0080 → 0.0045). DINA showed no decay and balanced misses.
- That is what a fixed family predicts. It is consistent with the source reading and not explained by it — a distinct identifier could behave this way for other reasons.

**Not to be acted on yet.** This is a code trace set against a manuscript description, and the manuscript may describe a different GRF configuration, or the condition may carry a requirement this trace does not address. CC recorded the mechanism, drew no conclusion and changed no caption; that was correct. Three source checks would settle it:

- whether anything in the enumerated pool can move under MR resampling (X and N are fixed, so the cut quantiles and `n_min` should be — confirm rather than infer);
- whether the DR pre-filter is genuinely outside the MR loop, not merely discarded at the one site traced;
- whether the alignment repair makes the within-draw selection rule identical to the rule that produced Ĥ, so the resampling reproduces the original selection.

If it holds, the captions change from conditional to unconditional, GRF's standing differs from DINA's, and handoff §3 needs revising. **Larry's call, against the manuscript.**

---

## 3. The per-cell reading

### 3.1 Field lower bound on β(Ĥ)

| Cell | Detection | Field lower [Wilson] | Field-s upper on β(Ĥᶜ) | Bonferroni joint-s |
|---|---|---|---|---|
| 12.4% HR 1.50 n 500 | 0.9985 | 0.9434 [0.932, 0.953] | 0.9319 | 0.9359 |
| 12.4% HR 1.50 n 1000 | 0.9975 | 0.9404 [0.929, 0.950] | 0.9519 | 0.9398 |
| 12.4% HR 1.50 n 1500 | 0.9975 | 0.9579 [0.948, 0.966] | 0.9529 | 0.9544 |
| 12.4% HR 1.75 n 500 | 0.9995 | 0.9415 [0.930, 0.951] | 0.9335 | 0.9350 |
| 12.4% HR 1.75 n 1000 | 0.9990 | 0.9364 [0.925, 0.946] | 0.9540 | 0.9399 |
| 12.4% HR 1.75 n 1500 | 1.0000 | 0.9535 [0.943, 0.962] | 0.9560 | 0.9575 |
| 31% HR 1.50 n 500 | 1.0000 | 0.9295 [0.917, 0.940] | 0.9185 | 0.9230 |
| 31% HR 1.50 n 1000 | 1.0000 | 0.9525 [0.942, 0.961] | 0.9460 | 0.9505 |
| 31% HR 1.50 n 1500 | 1.0000 | 0.9690 [0.961, 0.976] | 0.9480 | 0.9565 |
| 31% HR 1.75 n 500 | 1.0000 | 0.9300 [0.918, 0.940] | 0.9210 | 0.9260 |

- Rises with n at both prevalences and both hazard ratios — unlike DINA, whose 12.4% Block C bound degraded with n.
- Below nominal at n 500 (0.930–0.942) and reaching or passing it by n 1500 (0.954–0.969).
- The complement bound and the joint follow the same n-trend.

### 3.2 The two stratifications carry independent information

- **By `admitted_n` tertile:** monotone in all ten cells — coverage falls and retained bias rises from T1 to T3. Gaps 0.062–0.097 at 12.4% and 0.048–0.077 at 31%, narrowing with n. Conservative where few candidates clear the floor, anti-conservative where many do.
- **By p̂ tertile:** sharper. Gaps 0.075–0.198, retained bias swinging −0.232 to +0.313. Both attenuate with n.
- **The joint counts are near-uniform** (≈ 222 per cell against 222 expected under independence), with only a weak diagonal at 12.4% and a faint anti-diagonal at 31%. So the p̂ gradient is not a relabelling of the `admitted_n` gradient — this is the opposite of DINA, where the two were strongly anti-diagonal at ρ = −0.47.
- **Consequence:** on GRF the two are separate axes and both must be reported. On DINA they were confounded and the joint table was needed to resolve them.
- `admitted_n` = 1 never occurs; `admitted_n` ≤ 5 occurs on 1–11 replicates per cell, always with coverage 1.000 and strong over-correction (retained bias −0.16 to −0.58). Too few to read.

### 3.3 Location

- The field lower bound sits below the realized θ(Ĥ) in all ten cells, at 0.54–0.67 of it; the ratio and the paired ratio agree within 0.024 everywhere, so the ordering is not a median artefact.
- Shares reaching HR 1.00 run 0.023–0.405 and HR 1.25 run 0.005–0.110, both rising with n and with prevalence.
- Median naive exceeds θ(Ĥ) at every cell and falls toward it with n; est2 sits between the bound and θ(Ĥ) throughout.

### 3.4 The stratifier substitution, measured

- At each n the enumerated pool is identical quantile-for-quantile across both prevalences and both hazard ratios — 712 / 776 / 870 at n 500, 780 / 830 / 913–914 at n 1500 — while `admitted_n`'s median spans 115–449 across those same cells.
- ρ(`admitted_n`, `n_family`) runs +0.018 to +0.042 across all ten cells.
- Admitted share falls with n at 12.4% (0.147 → 0.091) and is flat at 31% (0.493 → 0.482).
- Part T2 was necessary, not merely convenient: stratifying on `n_family` would have stratified on a quantity the outcome never touches.

## 4. FS beside GRF

- Criterion-matched at 31% only (e1stud/cert20, `effMaxSG` ε 0.20). At 12.4% the comparator is `maxeffCons` ε 0.10, so a gap there cannot be read as engine behaviour even in part.
- Matched rows, field lower: GRF 0.9295 / 0.9525 / 0.9690 against FS 0.9745 / 0.9585 / 0.9615 at HR 1.50; GRF 0.9300 against FS 0.9700 at HR 1.75 n 500. FS is higher at n 500 by ~0.045 and the two converge by n 1000–1500, with GRF above FS at n 1500.
- Matched rows, IJ two-sided: indistinguishable (0.9840/0.9670/0.9730 against 0.9810/0.9710/0.9755).
- The confound travels with every row: identifier, family construction, detection set, and at 12.4% the criterion.

## 5. Dispositions

- Report and extraction accepted. No re-run.
- The `strat-xtab` fix accepted.
- Captions unchanged pending §2.
- **Open for Larry:** the Gate 0a determination and its three confirmation checks. Everything else about GRF's standing waits on it.
- Deferred work: 31% HR 1.75 at n 1000 and n 1500, and the six HR 1.00 null cells — roughly 6 h at the realized rate.
