# Hand-off — FDR Reporting Tool for a Given forestsearch Analysis

*Paste this into a new chat to start the FDR-tool work with full context. It is self-contained: it repeats the durable facts so the new chat needs nothing else.*

---

## Who / what

Larry F. León, Merck biostatistician, author of the `forestsearch` R package (CRAN v0.1.0; GLM extension on `feature/glm-extension` → v0.2.0). Works on Pop!_OS Linux (~127 physical cores), RStudio, Quarto, GitHub Desktop. Companion JASA paper "Post-Selection Inference for the Standard Trial Analysis of Data-Driven Subgroups" (León & Anderson) is submitted; this tool is adjacent to it.

**Hard working rules:** read files before editing; byte-level edits with exact-match anchors and occurrence-count assertions; never fabricate values/citations; every file-modifying response ends with a downloadable artifact via present_files; formal register; recommend before implementing (act on "proceed"); provide the files under discussion as downloadable artifacts; don't re-present unchanged files. Ground all claims about code in the actual source — do not reason from structure or memory (this cost real time in prior work: hypotheses about code behavior were repeatedly wrong until verified against a live traceback).

## The goal — state it precisely, because it is a deliberate reframing

**NOT** "find thresholds `{c_screen, c_cons}` that *control* the FDR at a nominal target." (That is the inverse/calibration problem Larry spent weeks on without success — do not pursue it.)

**INSTEAD:** for a **given** analysis with a **fixed** `{c_screen, c_cons}` that Larry actually ran, answer two questions by simulation:

- **(A)** What *is* the false-discovery rate? `A = P(declare any harm subgroup | c_screen, c_cons)` under a null.
- **(B)** Does adding a **de-biased-effect gate threshold `c_gate`** *reduce* it? `B = P(declare AND de-biased HR ≥ c_gate)`. Always `B ≤ A`; the gap `A − B` is what the de-biasing (Larry's submitted-paper method) buys.

So the tool **characterizes** the operating FDR of a specific analysis and **tests whether the gate helps** — it does not solve for thresholds. This dissolves the failure mode of the earlier calibration work.

## The deliverable

One function, roughly:

```r
fs_fdr_report(fs_analysis, df_analysis, c_gate = c(0.9, 1.0, 1.25), nsims = 1000, ...)
```

- Takes a **fitted** `fs_analysis <- forestsearch(df_analysis, ...)` and the dataset.
- **Inherits `{c_screen, c_cons}` and all knobs** from `fs_analysis$args_call_all` (see mechanism below) — so the FDR reported is the FDR *of the analysis actually run*, not a re-specified variant.
- Builds a **DGM null from the analysis** (see null choice below), defaulting sample size to `nrow(df_analysis)`.
- Simulates and reports **A (raw FDR)** and **B (gated FDR at each `c_gate`)** with **Wilson 95% CIs**, plus the **A − B reduction**.
- Output: a small table — one row per `c_gate` — that reads as "your analysis's FDR is A; the de-biased gate at `c_gate` brings it to B."

## What to merge — two existing files are the two halves

1. **`diagnostic_false_discovery.R`** — already contains the **(A)/(B) engine**. In its `one_rep()`: `declared <- !is.null(fs$sg.harm)` gives A; the de-biased HR is read as `fs$debias_gate$debiased$est` (HR scale, via `to_eff()` — NOT `$estimate`), and `gated_<g>` = `declared && deb_hr >= g` gives B at each gate. It reports "(A) declare-any" and "(B) declare AND de-biased HR ≥ g" at gates 0.9/1.0/1.25. It uses `parallel_args = list(plan = "sequential")` inside a `%dofuture%` `multisession` loop with `seed = TRUE`.

2. **`null_comparison_permutation_vs_dgm.qmd`** — contains the **DGM-null construction and the verified harness**: builds the borderline-homogeneous null via `setup_gbsg_dgm()` + `calibrate_k_inter()`, the **visible-error harness** (`declare_once` returns `list(flag, msg)` with flag ∈ {1 declared, 0 clean, −1 failure}; `fpr_over` counts failures and prints distinct messages; ends with a CLEAN/WARNING verdict), Wilson CIs, and a realized-HR table proving the null's structure.

The merge: take the DGM-null-from-analysis + harness from the `.qmd`, and the A/B/gate scoring from the `.R`, into one `fs_fdr_report()`.

## The inheritance mechanism (verified against source)

`forestsearch()` stores all its arguments on the return object: `args_call_all <- mget(args_names, envir = environment())` (forestsearch_main.R:996), returned as `fs_analysis$args_call_all`. Threshold mapping: `c1 = hr.threshold`, `c2 = hr.consistency`. A working inheritance wrapper already exists: **`fs_null_fpr.R`** (in outputs) pulls `c1`/`c2` from `args_call_all`, strips the args the calibrator manages (`df.analysis`, `hr.threshold`, `hr.consistency`, `seedit`, `parallel_args`), and passes the rest as `fs_params`. Reuse this pattern for `fs_fdr_report()`.

Note: `forestsearch::fpr_calibration()` **already exists** (fpr_calibration.R, Larry's v0.2.0 code) and does direct-simulation FPR at a specific `(c1,c2)`, supporting a `null_fn` for a DGM null. The FDR tool can either call it with a DGM `null_fn` or reimplement the loop — but it must add the **(B) gated layer**, which `fpr_calibration()` does not have. Decide which, in-chat.

## The null to simulate — CRITICAL, three distinct nulls (verified from source)

`setup_gbsg_dgm()` is a survreg-based AFT model. HR identity: `HR(H) = exp(−γ_treat/σ − γ_zh/σ)`, `HR(Hc) = exp(−γ_treat/σ)`. Harm subgroup H = {z1=1 AND z3=1}.

- **H0-hom** (`k_inter = 0`): true homogeneity, HR(H) = HR(Hc). Verified realized: HR(H)=0.641, HR(Hc)=0.650. "False positive" = declaring heterogeneity when the effect is uniform. **This is the cleanest FDR-of-heterogeneity null.**
- **H0-harm** (`target_hr_harm = 1.0` → calibrated `k_inter ≈ 0.568`): a REAL z1×z3 interaction tuned so the flagged region has HR=1.0 while treatment stays active elsewhere. Verified realized: HR(H)=1.008, HR(Hc)=0.650. "False positive" = calling the flagged region harmful when it is null there, amid real heterogeneity. The **more adversarial** null.
- **Permutation** (`sample(treat)`): **REMOVED — do not use.** It is the sharp no-effect null (not a subgroup-analysis scenario) AND it is structurally broken: it fails `forestsearch` on ~100% of calls (1000/1000 errors) because shuffling treatment breaks the search. Larry has decided against it in all forms (including within-arm, which still can't set an interpretable `hr_harm`).

For a *given analysis* the tool builds the null from the analysis's own DGM structure at the analysis's `{c_screen, c_cons}`. Default to H0-hom and/or H0-harm; let the user pick.

## The reference baseline the tool must reproduce

Clean end-to-end run (fix installed, nesting bug fixed, 1000 valid reps each, `fs_errors = 0`, verdict CLEAN):

- **H0-hom FDR = 0.138** [0.118, 0.161]
- **H0-harm FDR = 0.278** [0.251, 0.307]

Earlier (A)/(B) evidence at a leaner sequential declare-any config: raw ~0.174, gated at HR≥1.25 ~0.092 (gate helps meaningfully), gated at HR≥1.0 ~0.169 (gate barely helps). So the tool should show: **the gate reduces FDR substantially at a meaningful `c_gate`, marginally at HR≥1.0.**

## Architecture decision already made (respect it)

The pipeline is three steps: (1) **calibrate/characterize** the base identifier's declaration rate, (2) **de-bias** the selected subgroup's estimate (always — Larry's submitted paper), (3) **gate** on whether the de-biased HR clears `c_gate` (a threshold that MAY differ from `c_screen`/`c_cons`). `c_gate` feeds back into step 1 **only for null/operating-characteristic description, never at the analysis stage** (at analysis time you don't know truth, so the gate only annotates). In `fs_fdr_report()`, `c_gate` is a null-characterization parameter, inert on real data.

## Critical technical gotchas (learned the hard way)

- **The consistency bug is FIXED and merged** to `feature/glm-extension` (merge commit 6d0e878; fix commit 4f95ede). The fix: in the sequential loop of `subgroup.consistency`, `results_list[[m]] <- NULL` *deleted* the preallocated slot and shrank the list → out-of-bounds. Fix guards with `if (!is.null(res_m)) results_list[[m]] <- res_m`. Confirm the installed package has it: `deparse(forestsearch:::subgroup.consistency) |> grep("res_m", x=_, value=TRUE)` should show the guard.
- **Reinstall + restart R after any branch switch.** `devtools::install()` then Session → Restart R. The `multisession` workers load the INSTALLED package, not the source — skip this and they run stale code.
- **Nested parallelism breaks forestsearch.** If a `%dofuture%`/`multisession` loop calls `forestsearch` WITHOUT `parallel_args = list(plan = "sequential")`, forestsearch spawns multisession-inside-multisession → mass failures (this caused a 2624-error render). The FDR tool's inner forestsearch call MUST set `parallel_args = list(plan = "sequential")`.
- **Make failures visible, never masked.** The consistency bug hid for a long time behind three layers of `tryCatch → NULL` / `.errorhandling = "remove"`. The harness now codes failures (−1), counts them, prints messages, and emits a CLEAN/WARNING verdict. Keep this — silent error-masking is the villain.
- De-biased HR accessor: `fs$debias_gate$debiased$est` (HR scale). NOT `$estimate`.
- Emitted DGM column names: `y_sim` / `event_sim` / `treat_sim` / `id` / `flag_harm`.
- This is Claude Code work once specced — it needs the live namespaced package to run forestsearch. The sandbox cannot run compiled forestsearch.

## The bigger open question this tool informs

The corrected FDR (~0.14 to ~0.28) is roughly **double** the previously-believed 0.10–0.17, because the consistency bug was systematically suppressing declarations (biasing FDR *low*). Larry deems ~10–15% tolerable for exploratory identification, but **>30% an issue**; H0-harm at 0.278 sits at the edge. This reopens whether the max-statistic **calibration** (the companion-paper prototype) is optional polish or genuinely needed. The FDR tool provides the *measurement* that informs that decision — but the calibration-build decision itself is the companion-paper chat's concern, not this tool's.

## Deferred, SEPARATE package decisions — do NOT fold into this chat

Flag these and move on; they need different context:

1. **`Pcons` default** is `"split"` (subgroup_consistency_main.R:351) while the validated closed form `rate_closed` (consistency_resample.R:342) exists. A "measure-then-maybe-flip-default" decision — needs closed-form-vs-split context, not FDR context. (For the Cox/resample path the closed form is already used, so this is effectively settled there.)
2. **Remove `fpr_calibration()`'s permutation default** in favor of a required DGM `null_fn` — a wrapper-level or API decision. Recommend doing it wrapper-level in `fs_null_fpr()`.

## First action for the new chat

Confirm the merge spec, then build `fs_fdr_report()` by merging `diagnostic_false_discovery.R` (A/B engine) and `null_comparison_permutation_vs_dgm.qmd` (DGM null + harness), with `args_call_all` inheritance, the H0-hom/H0-harm null choice, sequential inner parallelism, the visible-error harness, and A/B/gate output with Wilson CIs. Reference baseline to reproduce: H0-hom 0.138, H0-harm 0.278. Then run in Claude Code against the installed (fixed) package.
