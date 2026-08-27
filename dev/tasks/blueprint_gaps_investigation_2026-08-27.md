Save this message as dev/tasks/blueprint_gaps_investigation_2026-08-27.md and commit it, then execute. Base dca84c93. Gates internal: stop only on failure.

Rules: anchor by content, never line number. Verify from source, never from description. Never push. Do not touch: dev/glm-continuous-sims/eval_bracket_direct.R, REPORT_blockers_1-3_ADDENDUM.md, a1_anchor_correction.py, the two frozen parent documents, the labelled provenance comment at driver site 1.

PART 1 — housekeeping, one commit.
NEWS.md has no entry for fs_sym_root. Check which of fs_sym_root, fs_dgm_scale, fs_scale_se, fs_mr_oc_summary are missing entries and add what is missing, matching the file's existing style and the v0.2.0 development section. No other change.
Then, read-only, no fix: grep R/ for roxygen raw Rd macros containing a literal % — the \code{...%...} pattern that silently dropped the scale param from fs_sym_root.Rd. Report the count and the files. Hypothesis: markdown escaping does not apply inside raw Rd macros you write yourself, so backticks are safe where \code{} is not. Confirm or refute against one example. Convert nothing.

PART 2 — read-only investigations, no writes, report only.

(a) ONE QUESTION ONLY. .fs_oc_median is being deleted; the decision is made, do not investigate whether to keep it. Its only real use is medians of confidence-interval lengths, and the replacement is an optional flag in fs_mr_oc_summary() swapping mean for median on CI lengths only. Report: is the mean over CI lengths computed at a single site inside fs_mr_oc_summary(), or at several? Quote the site or sites. That is the whole question — it decides whether the flag is a one-line change or not worth doing. Do not delete anything and do not add the flag.

(b) The str2 / preanti <= 0 alias — this is the priority. beta(Hhat) depends on which of two in-sample-equivalent subgroup strings the search emits, shifting super-population membership in five rows per cell. Report: where candidate strings are constructed and emitted; what makes the two forms in-sample-equivalent but super-population-different; whether the choice between them is deterministic or depends on ordering, worker count, or floating-point comparison; the full signature and body of .maxeff_membership_dedup at R/subgroup_consistency_helpers.R:2036 and any evidence about why it has no call site; and what canonical-form options exist (fewest cuts, lowest names.Z index, something else). Do not propose a fix — I need the options, not a choice.

(c) Search reproducibility across worker counts. Read the parallel RNG scheme in subgroup_consistency_main.R and report it exactly: there is a claim it pre-generates one L'Ecuyer-CMRG stream per candidate indexed by global position, which would make results invariant to worker count by construction. Confirm or refute from source. Also report R/subgroup_search.R:252 — whether that loop draws at all, and so whether the missing future.seed matters. Then state what a reproducibility test would need to run and roughly what it would cost.

Report once: the NEWS.md diff, the raw-Rd grep result, and the three items. Then git status --short and git diff --stat dca84c93..HEAD.
