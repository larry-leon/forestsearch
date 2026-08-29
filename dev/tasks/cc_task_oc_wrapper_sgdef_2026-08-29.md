# CC task — the realized selection distribution: `sg_def` tabulated, evaluated on the population, and compared against the analytic

**File:** `dev/tasks/cc_task_oc_wrapper_sgdef_2026-08-29.md` · **Issued:** 2026-08-29 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** Larry places this file in `~/Downloads`; CC copies it into `dev/tasks/` and commits it first (§1).
**Follows:** `HANDOFF_oc_wrapper_2026-08-29_v1.md` §5 and §9, and `REPORT_oc_wrapper_null_2026-08-29.md`. Read the handoff's §5 first — this task exists to resolve its (b) and (c).

**What it is for.** The wrapper lands on the measured record but not uniformly within Monte-Carlo noise: `E|Ĥ|`, sensitivity and PPV are low by a small signed amount, and the naive bias is off by ~2% with a sign that flips across cells. Two explanations are currently tangled:

- **within-rule** — for the rule the search actually picked, its *realized* sample size and purity exceed their population values, because the rule was selected partly for being favourable in that sample. The analytic side reports population functionals. If this is the whole story, the gaps are definitional and already predicted by the inventory's B5 row 4.
- **between-rule** — the analytic selection distribution puts weight on *different rules* than the search does. That would mean the family, the argmax, or both differ, and it is a substantive finding rather than a definitional one.

The payloads carry `sg_def` — the winner's rule verbatim — on every replicate of all three cells, and nothing has looked at it. This separates the two.

---

## ⚠ CATEGORY

**No `R/` change of any kind.** No new exports, no edits to any package file. Read-only apart from one scratch script, its `.rds`, and the report, all under `dev/glm-continuous-sims/`.

**Compute: minutes, DIAGNOSTIC only.** Evaluating stored rule strings on `df_super` and tabulating stored columns. No replicates, no simulation, no renders, no `fs_oc_predict()` re-runs — the analytic side is read from the existing `oc_wrapper_grid_2026-08-29.rds`.

---

## 1. Provenance and first commit — GATE

```bash
hostname; pwd
git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD
git status --porcelain
git log --oneline -4
Rscript -e 'cat("forestsearch", as.character(packageVersion("forestsearch")), "\n")'
```

*GATE:* branch `feature/glm-extension`; HEAD descends from `648c2af8`; installed `0.2.6`.

```bash
ls ~/Downloads/cc_task_oc_wrapper_sgdef*
cp ~/Downloads/cc_task_oc_wrapper_sgdef*.md dev/tasks/cc_task_oc_wrapper_sgdef_2026-08-29.md
git add dev/tasks/cc_task_oc_wrapper_sgdef_2026-08-29.md
git commit -m "task doc — realized selection distribution from sg_def"
```

*GATE:* exactly one match (hyphens may be stripped; copy under the name above, report both).

---

## 2. Tabulate the realized selection distribution

For each of the three tracked payloads — alternative n = 500, alternative n = 700, null n = 500 — read `results$sg_def` over the detected replicates and report:

**Two units of aggregation, because they answer different questions:**

- **(a) verbatim string.** Distinct `sg_def` values with counts. Expect heavy fragmentation: the search cuts at *sample* quantiles, so two replicates choosing "the same rule" carry different numeric thresholds (`cd40 <= 415` vs `cd40 <= 419`). Report the number of distinct strings and the top 15 with counts.
- **(b) covariate signature.** The set of variables used and each one's direction, with the numeric thresholds stripped — e.g. `{cd40 ≤ ·} & {cd80 ≤ ·}` negated or not. This is the stable unit and the one to compare against the analytic labels. Report every signature with a count, and the count of distinct signatures.

Also report `K` (how many factors the winner uses) as a distribution, and the fraction of detected replicates whose winner uses each covariate at all.

**The question this answers directly:** does the realized selection distribution live on the same covariates as the analytic family's high-`p_sel` candidates? The one winner quoted in the inventory (`"!{cd40 <= 415} & !{cd80 <= 1022}"`) uses neither of the prediction document's two axes. State plainly whether that is typical.

---

## 3. Evaluate each realized rule on the population

**Use the package's own cut-expression resolver — do not write a parser.** The `sg_def` strings use the search's `{...}` syntax; locate the helper that resolves them (the same machinery `fs_build_eval_frame()` / the eval-frame path uses) and apply it. If no resolver is callable, report that and fall back to a documented, explicitly-stated transformation — and say exactly what you did.

For each **distinct verbatim** `sg_def`, evaluate its membership on the MD40 (and null) `df_super` — the same deterministic DGM rebuild the earlier scripts use, gated against the payload's committed `truth` as they do — and compute its **population** quantities:

| quantity | definition |
|---|---|
| `Pg_pop` | `mean(memb)` on `df_super` |
| `nPg_pop` | `n * Pg_pop` at that cell's `n` |
| `PQg_pop` | `P(g ∩ Q)/P(g)` — the population PPV (`0` in the null cell) |
| `sens_pop` | `P(g ∩ Q)/P(Q)` (`NA` in the null cell) |
| `spec_pop` | `1 − P(g ∩ Qc)/P(Qc)` |
| `beta_pop` | the mixture mean `tauQc + bint * PQg_pop`, from `fs_dgm_scale()` |

Report how many distinct rules failed to evaluate and why; do not silently drop them.

---

## 4. The within-rule comparison — the discriminating test

For each covariate signature with **at least 20 detected replicates**, compare realized against population, over the replicates that selected it:

| column | realized (from `results`) | population (from §3) |
|---|---|---|
| size | `mean(n_harm)` | `mean(nPg_pop)` over those replicates' own rules |
| PPV | `mean(ppv)` | `mean(PQg_pop)` |
| sensitivity | `mean(sens)` | `mean(sens_pop)` |
| specificity | `mean(spec)` | `mean(spec_pop)` |
| effect | `mean(betaHhat_H)` (oriented as the OC summary orients it) | `mean(beta_pop)` |

Report the paired difference and its SE per signature, and pooled across signatures weighted by count.

*This is the discriminating result.* State it plainly:

- **If realized exceeds population within signatures**, by roughly the magnitude of the handoff's §5b gaps, then those gaps are a within-rule sample-realization effect — definitional, already predicted, and the analytic side is not wrong. Say so and give the magnitude.
- **If realized and population agree within signatures**, the §5b gaps are *between*-rule: the analytic selection distribution weights different rules than the search does. That is a substantive finding about the family or the argmax. Say so, and §5 below localises it.
- Report both directions honestly; a mixed result across quantities is a real possibility and should not be forced into one story.

**Do the same for the naive bias** (handoff §5c), which is the unexplained one: per signature, `mean(naive_estimate − betaHhat_H)` realized, against the analytic `Enaive_bias`. Note whether its sign behaves the same way within signatures as it does in aggregate — the aggregate sign flips between n = 500 and n = 700, and this is the first chance to see whether that is a mix effect.

---

## 5. Measured against analytic selection distribution

From `oc_wrapper_grid_2026-08-29.rds`, take the analytic per-candidate selection probabilities `p_sel` (or `sel_c`) at the driver's `(c1, c2)`, both gates, per cell. Map the analytic candidate labels onto §2(b)'s covariate signatures — **by signature, not by string** — and report:

- the top 15 analytic signatures with their probabilities, beside the measured frequencies for the same signatures;
- total analytic probability mass on signatures the search **never** selected;
- total measured frequency on signatures **absent from the analytic family** (if any — say why: filtered by a floor, not representable by the enumeration's cuts, or something else);
- a single overlap summary — the total variation distance between the two distributions over signatures, or an equivalent, with the measure you used stated.

If the label-to-signature mapping is ambiguous for some candidates, say which and how many; do not force a match.

---

## 6. Report

`dev/glm-continuous-sims/REPORT_oc_wrapper_sgdef_2026-08-29.md`:

1. Provenance (raw).
2. §2 tabulations, all three cells.
3. §3 evaluation, with failures accounted.
4. §4 the within-rule table and **the verdict, in the terms §4 sets out**.
5. §5 the distribution comparison.
6. What this settles in the handoff's §5b and §5c, and what it leaves open — stated as facts, no recommendation.
7. Ten-line verdict.

Save the per-rule table as an `.rds` beside the report so it can be re-read without recomputing.

Close out: `git status --porcelain`; stage by explicit path (the script, its `.rds`, the report); commit; `git log --oneline -3`. **No push.** No `devtools::install()` — nothing in `R/` changed.

---

## 7. Out of scope

- No `R/` change, no new exports, no edits to the prediction document, the verification document, or any driver.
- No `fs_oc_predict()` / `fs_oc_grid()` re-runs — the analytic side comes from the stored `.rds`.
- No simulations, no replicate runs, no renders.
- No recommendation about what to do with the result; §4's verdict is a statement of which explanation the evidence supports, not a proposal.
