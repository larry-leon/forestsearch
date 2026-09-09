# TASK — Part D2: the `ci_method` default becomes `"field"`, and the complement NOTE is updated to the certified wording

Date: 2026-09-09. Author: chat (spec). Executor: Claude Code (Linux). Approver: Larry (2026-09-09). Reviewer: the Linux MR-field chat.
Predecessors: `REPORT_defaults_flip_2026-09-08.md` (Part D, which left `ci_method` as the open decision), `REPORT_cert20_2026-09-08.md`, `REPORT_tier2_2026-09-08.md`, `REPORT_fixedphat_ij2s_2026-09-09.md` (C-4), `dev/notes/NOTE_complement_product_2026-09-08.md`.

**Why now.** Part D flipped `field_complement`, `field_scale_complement` and `return_reselection` to the recommendations but left `ci_method` at `"ij"`. The field block is gated on `if (ci_method == "field")`, so under the shipped defaults those three flips are **inert** and the certified products (field one-sided lower on β(Ĥ), field-s one-sided upper on β(Ĥᶜ), the Bonferroni joint pair) are unreachable without the user setting `ci_method` explicitly. C-4 closed the open question. **Larry's constraint, and the point of Gate D2c below: the IJ two-sided interval must not be lost** — the IJ SE is computed unconditionally, before and independent of the field gate, so `"field"` adds the field block without removing anything.

## Protocol

- First action: archive stale variants in `~/Downloads` to `~/Downloads/cc_archive/` (**do not archive `PROPOSAL_field_recovery_2026-09-09.md` — it is current and belongs to a later task**); copy this file to `dev/tasks/` and commit. Do not push. If missing from `~/Downloads`, reconstruct from the kickoff, commit under this name with a reconstruction note, and proceed.
- **One `R/` change, classified *changes behaviour*, in exactly two files** (`R/fs_mr_inference.R`, `R/forestsearch_main.R`), plus `man/`, `NEWS.md` and the NOTE. Nothing else under `R/`.
- Compute: verification renders only (≤ 5 replicates each). No campaign, no bundle beyond the gate bundles.
- Gates stop on failure; **on Gate D2 failure revert the touched files to HEAD, re-install, record the failure and the differing columns, and stop** (there is no Part T to fall through to).
- Standing conventions: winner-only and winner-floor excluded; verify from source; records beside results. Leave the seven pre-existing untracked files alone.

## Stage 0 — Discovery (quote from HEAD)

Quote with current line numbers: the `ci_method` formal of `fs_mr_inference()` and its roxygen `@param`; the gate `if (ci_method == "field")`; the unconditional `se_ij <- .fs_mr_se_from_ij(...)` line and the `se <- if (ci_method == "wald") ...` line, establishing that the IJ SE is computed before and independent of the gate; the `.g_mr(mr_inference_args$ci_method, "ij")` fallback in `forestsearch_main.R`; whether the template exposes a `ci_method` knob (name it) or hard-codes the value, and how the campaigns set it. Also quote `.fs_apply_mr()`'s `ci_method` handling (`R/fs_mr_inference_methods.R`) and state whether the DINA/GRF branches inherit the new default or fall back independently.

## Part D2 — The flip

1. `fs_mr_inference()`: `ci_method = c("field", "ij", "wald")` so `match.arg` yields `"field"`. Roxygen updated: `"field"` is the recommended default and runs the field block (the certified one-sided products); the IJ elements are returned side by side in every call regardless; `"ij"` restores the prior default and omits the field block.
2. `R/forestsearch_main.R`: `ci_method = .g_mr(mr_inference_args$ci_method, "field")`, comment rewritten to name the new default and the value that restores prior behaviour.
3. Template: if a `ci_method` knob exists, its env default becomes `"field"` and the knob echo is updated; if the value is hard-coded, **leave it** and record that campaigns set it explicitly (do not add a knob in this task).
4. `NEWS.md`, one bullet under the development header: the default `ci_method` is now `"field"`, so the recommended constructions are what a default call produces; the IJ SE and its intervals are unchanged and still returned; prior behaviour is `ci_method = "ij"`; note the added per-fit cost of the field block.
5. `devtools::document()`; `devtools::install(dependencies = FALSE)`; `deparse()` of the installed bodies equals source; read the installed formals and the installed `forestsearch()` body's `.g_mr` fallbacks back and quote them.

**Gate D2** (5 replicates, effMaxSG ε 0.20 HR 1.50 n500 config, seeds `8316951 + sim_id`, sim_id 1–5, `FS_S7_FIELD_DECOMP=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none`, tags `cim_unset` / `cim_field` / `cim_ij`; timing columns excluded throughout):

- **D2a — unset ≡ explicitly `"field"`.** Every non-timing column `identical()`, `truth` `identical()`; resolved `meta$ci_method` under *unset* reads `field`.
- **D2b — the unset render reproduces the committed `e1stud` rows 1–5** (which ran `ci_method = "field"` explicitly): every pre-existing non-timing column `identical()`, `truth` `identical()`.
- **D2c — IJ IS NOT LOST (the criterion Larry named).** Render with `ci_method = "ij"` explicitly and check all three of: (i) every `^mr_H_`, `^mr_Hc_` and IJ-derived column is finite and **`identical()` to the same columns of the unset (`"field"`) render** — the field gate does not touch the IJ path; (ii) the field / field-s / joint columns are all `NA` in the `"ij"` render; (iii) in the **unset** render both families are present simultaneously — the IJ two-sided bounds finite **and** the field, field-s and joint columns finite on every detected row. Report the IJ two-sided bounds from both renders side by side for the five replicates.

On PASS commit with `REPORT_cimethod_flip_2026-09-09.md` (Stage 0 quotes, the diff summary, all three gate results with concrete values, and the D2c side-by-side table).

## Part N2 — The complement NOTE, updated to the certified wording

Rewrite `dev/notes/NOTE_complement_product_2026-09-08.md` as `dev/notes/NOTE_survival_products_2026-09-09.md` (new file; leave the old one in place and add a one-line pointer at its top saying it is superseded). Content, as a rule statement:

> **Survival post-selection products (certified 2026-09-09).** Defaults are the recommendations: `ci_method = "field"`, `field_complement = TRUE`, `field_scale_complement = "selected"`, `return_reselection = TRUE`.
> **Harm subgroup Ĥ — one-sided lower bound on β(Ĥ), the field.** Coverage 0.944–0.974 across ten harm cells at 12.4% and 31% prevalence, n = 500 / 1000 / 1500; flat in n.
> **Complement Ĥᶜ — one-sided upper bound on β(Ĥᶜ), field-s** (the studentized complement field, R1). **The small-sample shortfall closes with n:** 0.912 → 0.942 → 0.947 (31%, HR 1.50) and 0.919 → 0.942 → 0.946 (HR 1.75); ≥ 0.94 at every cell of the 12.4% dominated regime (0.941 / 0.956 / 0.961). field-s is at or above the unscaled field in all thirteen cells.
> **Joint two-subgroup claim:** Bonferroni, γ = 0.025 each side; 0.939–0.963 at every cell.
> **Two-sided intervals are not certified.** The IJ two-term two-sided interval is retained and reported as the secondary, conservative option and is the only two-sided construction offered, but its harm-block coverage falls to 0.913–0.917 at 12.4% prevalence with n ≥ 1000 (0.971–0.981 at 31%). The field's own two-sided is lower still (0.878–0.930). Read two-sided statements at low prevalence and large n with that caveat.
> **Analysis-time diagnostic.** p̂(Ĥ), the field's re-selection frequency, is recorded and reported; **no construction reads it**. Harm-block bias is a monotone increasing function of p̂, crossing zero near p̂ ≈ 0.5: **over-correction at low p̂** (bias −0.11 to −0.28 log units at 12.4%, n = 1500) and **under-correction at high p̂** (the stable-pick regime, +0.02). The one-sided products are each exposed to one pole only, and in the conservative direction, which is why they certify while the two-sided interval — exposed to both — does not. Complement coverage is a stable function of p̂ (flat-to-improving at fixed band in all 18 usable band × series combinations); the harm block's is not common across prevalence, so its p̂ flag is directional, not calibrated.
> **Caveats on record.** ε > 0.25 is not adoptable (Larry, 2026-09-08); under a pure-size pick (`maxSG`) the naive SE mis-calibrates (0.876 of the error SD) and field-s inherits it; fixed p̂ bands are fixed in value but not in meaning, since the p̂ distribution shifts with n (replicates below p̂ = 0.20 are 52% of the 12.4% HR 1.75 cell at n = 500 and 8% at n = 1500).

Commit the NOTE with the flip, or as its own commit — CC's choice, stated in the record.

## Done means

Stage 0 quotes; Part D2 committed on PASS with D2a/D2b/D2c concrete values and the D2c side-by-side table (or the failure recorded and the tree reverted); the NOTE committed with the pointer added to the superseded file; branch left unpushed; one-paragraph closing summary with the gate results and the commit range. **Out of scope:** the print/summary method, the vignette, the field-recovery diagnostics (`PROPOSAL_field_recovery_2026-09-09.md`), any campaign, any other `R/` change.
