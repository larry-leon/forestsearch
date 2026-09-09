# TASK — Reporting the certified survival products: a `print`/`summary` extension and the survival post-selection vignette

Date: 2026-09-09. Author: chat (spec). Executor: Claude Code (Linux), unattended (Larry offline several hours). Approver: Larry (2026-09-09). Reviewer: the Linux MR-field chat.
Predecessors: `dev/notes/NOTE_survival_products_2026-09-09.md` (the certified products and their caveats — **the single source of truth for every claim printed or written below**), `REPORT_cimethod_flip_2026-09-09.md`, `REPORT_cert20_2026-09-08.md`, `REPORT_tier2_2026-09-08.md`, `REPORT_fixedphat_ij2s_2026-09-09.md`.

**Scope note.** `print.forestsearch()` (`R/forestsearch_methods.R:86`) and `summary.forestsearch()` (`:179`) already exist and neither reports the MR/field results. This task **extends** them; it does not author new methods. Transplant-first throughout: reuse `.fs_get()`, `.fs_sg_labels()`, `.fs_format_admission()` and the existing `cat()` idiom rather than inventing a formatting layer.

## Protocol

- First action: archive stale variants in `~/Downloads` to `~/Downloads/cc_archive/` — **do not archive `PROPOSAL_field_recovery_2026-09-09.md` or `HANDOFF_guohe_comparison_2026-09-09.md`**, both current. Copy this file to `dev/tasks/` and commit. Do not push. If missing from `~/Downloads`, reconstruct from the kickoff, commit under this name with a reconstruction note, and proceed.
- Compute: verification renders and one vignette build only. No campaign, no committed simulation bundle.
- `R/` changes confined to `R/forestsearch_methods.R`, classified **adds code; existing output byte-identical when no MR results are present**. No change to `fs_mr_inference.R`, `forestsearch_main.R`, or any construction.
- Gates stop on failure; **on failure revert the touched files to HEAD, re-install, record the failure, and stop.**
- Standing conventions: winner-only and winner-floor excluded; bounds read **by location**, never as significance at the null; no claim printed or written that the NOTE does not support. Leave the seven pre-existing untracked files alone.

## Stage 0 — Discovery (quote from HEAD; STOP if the object model is not as assumed)

1. **Where do the MR results live?** Establish from source whether a `forestsearch()` call with `mr_inference = TRUE` attaches the `fs_mr_inference()` return to the returned object, and under what name — quote the assignment in `R/forestsearch_main.R` and the returned list's construction. Quote the structure of the `fs_mr_inference()` return: the harm block (`debiased$lower_1s` and its siblings), `field$complement` (`upper_1s`, `upper_1s_s`, `se_field_s`), `field$joint` / `field$joint_s` (`bonf_loH`, `bonf_upHc`, `gamma`), the IJ elements (`se_ij`, the two-sided bounds), `p_hat` and `p_hat_sum`. **If the MR object is not attached to the `forestsearch` object, STOP and report** — the methods cannot print what they cannot reach, and the fix is a decision, not an edit.
2. Quote `print.forestsearch()` and `summary.forestsearch()` in full, plus `.fs_get()`, `.fs_sg_labels()` and `.fs_format_admission()`.
3. **Vignette infrastructure:** confirm there is no `vignettes/` directory and quote DESCRIPTION's `Suggests`/`VignetteBuilder` fields. Note whether `quarto` is a declared dependency. Report what building a vignette would add to DESCRIPTION before adding it.
4. Confirm the effect scale and orientation used by the printed quantities (harm block on the HR scale via `to_eff`; larger = more harm), so no bound is printed on the wrong scale.

## Part P — The reporting extension

**P1. `print.forestsearch()`** gains, after the existing subgroup block and only when MR results are present (absent ⇒ output byte-identical to today):

```
Post-selection inference (certified products):
  Harm subgroup H:        one-sided 95% lower bound on HR   <value>
  Complement Hc:          one-sided 95% upper bound on HR   <value>   [field-s]
  Joint (Bonferroni):     H lower <value>, Hc upper <value>  (gamma <value> each side)
  Two-sided (IJ, secondary): H (<lo>, <hi>)
  Re-selection frequency  p-hat(H) = <value>
```

- Print the **certified one-sided products first**, the joint second, the two-sided IJ third and explicitly labelled secondary, p̂ last as a diagnostic.
- One caveat line, printed only when it applies and worded from the NOTE, not invented: a two-sided caveat pointing to the NOTE, and — when p̂ is high — the directional stable-pick note. **Choose the p̂ threshold from the NOTE's own statement (bias crosses zero near p̂ ≈ 0.5) and say in the roxygen that the threshold is descriptive, not calibrated.**
- Never print a bound as significant or not; location only.

**P2. `summary.forestsearch()`** gains the same block plus: the field-s SE and the naive complement SE beside it; the IJ two-term SE; the top-3 re-selection mass (`p_hat`'s three largest entries) and `p_hat_sum`; the constructions' names as they appear in the NOTE. One short paragraph, printed once, naming what is certified and what is not, sourced from the NOTE.

**P3. Roxygen** for both methods updated to describe the new block, name `dev/notes/NOTE_survival_products_2026-09-09.md`, and state that `ci_method = "field"` (the default) is what produces the field and field-s bounds. `devtools::document()`; NAMESPACE unchanged (both are already `@export`ed).

**Gate P** (no new simulation; use the committed `cim_unset` gate bundle or a 2-replicate render of the standing config):
- **Pa — absent-MR invariance:** on a `forestsearch` object with `mr_inference = FALSE` (or MR results absent), `capture.output(print(x))` and `capture.output(summary(x))` are **identical** to the same capture from the pre-change installed package. Build both captures explicitly and diff them.
- **Pb — present:** on an object carrying MR results, every printed number matches the corresponding element of the MR object to the printed precision (assert programmatically, element by element; quote the table in the record).
- **Pc — no invented claims:** every caveat string printed appears in, or is a faithful compression of, the NOTE; quote each printed caveat beside its NOTE source line in the record.

## Part V — The survival post-selection vignette

- `vignettes/survival-post-selection.qmd` (Quarto; **do not convert any existing `.Rmd`/`.qmd`**), a **smoke-test-style living document exercising the pipeline end-to-end** on a public dataset (GBSG, as the applications documents use). Not a unit-test scaffold; not new methodology.
- **Transplant-first:** copy the structure and the fitting/reporting code from the committed `quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd`, cutting it down to a vignette-sized run — no bootstrap, no cross-validation, no Guo & He, no LOO. Keep the runtime inside CRAN's vignette budget and state the measured build time in the record.
- Sections: the identified subgroup; the certified products with the print/summary output shown; how to read the bounds **by location**; p̂ as a diagnostic and what a high value means (directional, from the NOTE); what is **not** certified (the two-sided interval at low prevalence with n ≥ 1000); and a pointer to the NOTE and the campaign records.
- DESCRIPTION plumbing (`VignetteBuilder`, `Suggests`) only as Stage 0 established is required. **CRAN compliance is a standing priority: run `devtools::check()` (or `R CMD check --as-cran`) after the vignette exists and report the full NOTE/WARNING/ERROR tally**, comparing against the same check on HEAD before this task so any new item is attributable.

## Done means

Stage 0 quotes (or a STOP with the reason); Part P committed with Pa/Pb/Pc concrete values; Part V committed with the built vignette, its measured build time, and the before/after `check()` tallies; `REPORT_print_vignette_2026-09-09.md` beside the other records; full test suite re-run and its tally reported (target FAIL 0); branch left unpushed; one-paragraph closing summary with the gate results, the check tallies and the commit range. **Out of scope:** the field-recovery diagnostics, `.fs_apply_mr()`'s `ci_method` default, any construction change, any campaign.
