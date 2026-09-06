# REPORT — GBSG frozen-family interval constructions: STOPPED AT GATE 0

**Task:** `dev/tasks/TASK_gbsg_frozen_intervals_2026-09-05.md` (084ec961). I1–I6 at defaults.
**Date:** 2026-09-06. Executor: Claude Code.

---

## GATE 0: STOP — Table 2's source element cannot exist under the wrapper without an `R/` change, which this task forbids

**The finding (0a, verified from source).** The task's change 1 prescribes the gate call
`mr_inference_args = list(draws = mr_draws, ci_method = "field", return_reselection = TRUE, field_uniform = FALSE)`
and Table 2 reads "re-selection diagnostics from `fs$mr_inference$reselection`: p̂_Ĥ; the top three re-selected candidates …". But `forestsearch()` does **not** forward `return_reselection` — zero occurrences in `R/forestsearch_main.R` (the gate call at `:3368–3390` forwards `t_confirm`, `confirm_rule`, `reselection`, `effect_neighborhood`, `selection_rule`, `draws`, `multiplier`, `include_complement`, `ci_method`, `field_uniform`, `seed` and nothing else) and zero in `R/fs_mr_inference_methods.R`. The argument would be silently dropped, `fs_mr_inference()` would run at its `return_reselection = FALSE` default (`R/fs_mr_inference.R:405`), and `fs$mr_inference$reselection` would be `NULL`. Same pass-through gap class as `field_R_out`/`field_R_in` (s7 Stage 0) and `field_uniform` (uniform Stage 0 — closed there by an authorized add-only line).

**Why no workaround is faithful:**
- A document-side replay of the gate with `return_reselection = TRUE` cannot reproduce the wrapper's call: `df.fs`, the enumerated candidate family and the resolved admission are internal to `forestsearch()`. The one exportable family — `fs_to_guohe()`'s — has **116** candidates versus the gate's own **133** (committed payload `extras$mr$n_family`), so p̂ computed on it would describe a different selection problem.
- Omitting p̂_Ĥ and the top-3 frequencies guts Table 2 (the fit does retain `selection_rate` and the field's λ diagnostics, but the task's regime-reading rests on p̂_Ĥ).

**Candidate fix (NOT applied — `R/` is frozen for this task):** one add-only line in the `forestsearch()` gate call, `return_reselection = .g_mr(mr_inference_args$return_reselection, FALSE),` — identical in kind to the `field_uniform` forwarding Larry classified on 2026-09-05 as "add-only pass-through, no behaviour change" (default `FALSE` reproduces today's output byte-for-byte; the element is additive on the return). With it, the task proceeds as written. Alternatives: (a) authorize the line under this task with the usual guards (suite + a fixed-seed ij/field byte-identity spot-check); (b) amend Table 2 to drop p̂ diagnostics (weakens the illustration's stated purpose); (c) defer.

## What Stage 0 established before the stop (for reuse on resumption)

- **0a.** `fs$mr_inference$field` under the wrapper carries `lambda_mean, lambda_sd, q05, q25, q50, q75, q95, q025, q975, n_out_used, n_in_used_mean, est2, lower_1s, lower_2s, upper_2s, se_field, lower_se, upper_se, R_out, R_in, seed_offset, timing_seconds` (source `R/fs_mr_inference.R:691–709`; degenerate case a `note`). `reselection` (when returned) is `list(winner, p_hat)` with `sum(p_hat) == selection_rate` (`:738–742`). `fs_bc$H_estimates` is a 1-row table `H0, sdH0, H0_lower, H0_upper, H1, sdH1, …, H2, sdH2, H2_lower, H2_upper` — **an SE but no draws**, so decision I2 lands on its middle option: FB one-sided `exp(log H2 − 1.645 · sdH2/H2)`.
- **0b.** Committed payload (identity anchor): `_payloads_2026-09-01_complete/analysis_gbsg_survival_frozen_family/analysis_gbsg_survival_frozen_family_payload.rds` — `extras$mr` (naive est/lower/upper 2.22/1.17/4.22; debiased est/lower/upper/lower_1s 1.75/0.68/4.52/0.79, se_ij 0.484, ij_source "ij", n_family 133), `extras$fb$H_estimates` (H2 1.96 [0.81, 4.74], sdH2 0.883), `extras$gh` (naive 2.22, debiased 1.36, bound 0.962, n_family 116), selected subgroup `{pgr <= 32.5} {er <= 0}`. **FB is seed-reproducible:** `forestsearch_bootstrap_dofuture()` has formal default `seed = 8316951L` (the document omits the argument, so the fixed default applies; plan-invariance is the document's own verified note).
- **0c.** LOO cache present (committed): `_payloads_2026-09-01_complete/analysis_gbsg_survival_frozen_family/cv_out_analysis_gbsg_survival_frozen_family_maxeff_neighborhood.rds`; the document's default path looks in `gh_dir` (`quarto/GuoHe/`, where the file no longer lives), so the re-render sets the `LOO_CACHE` env var to the committed path (read-only) — or recomputes in ~18 s (memory-recorded, seeded).

Nothing was rendered, edited, or recomputed; committed payloads untouched. **Stopped at Gate 0.** Decision needed from Larry on the `return_reselection` pass-through.
