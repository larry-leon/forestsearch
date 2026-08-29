# CC task (overnight, unattended) — settle the confounder-list mismatch, correct the analytic family, and test whether it closes the residuals

**File:** `dev/tasks/cc_task_oc_wrapper_confs_2026-08-30.md` · **Issued:** 2026-08-29 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** Larry places this file in `~/Downloads`; CC copies it into `dev/tasks/` and commits it first (§1).

**Follows and builds on `REPORT_oc_wrapper_sgdef_2026-08-29.md` (commits `69ccc84d`, `b8bd0f0a`, `f2ca69e6`) — that report stands and is not superseded.** It established, from the measured side alone:

- **Size (`E|Ĥ|`) is between-rule.** Within every signature, realized `n_harm` matches the rule's population size (pooled −0.56 ± 0.24, −0.02 ± 0.26, −0.50 ± 0.24). The gap lies in the weights: the rules the search picks are larger than the analytic `sel_c` weighting expects by +2.0 / +0.6 / +1.5 subjects, and re-weighting measured rules by analytic `sel_c` reproduces the analytic figure.
- **PPV is within-rule at n = 500** (+0.0047 ± 0.0016, ~1%); at n = 700 within noise. Sensitivity in between.
- **Naive bias tracks the winner's size and purity** (corr −0.6 to −0.95), so a few-percent shift in signature mix moves the aggregate by ≈1 — the size of the discrepancy, and consistent with its sign flip.
- **3–5% of measured selection mass sits on `str2` signatures the analytic family cannot contain.** The driver ran `include_str2 <- TRUE` (13 confounders); the analytic family was built on 12.

**Everything unexplained now points at one thing: the analytic family may have been enumerated over the wrong covariate set.** This task settles that and measures whether correcting it closes the residuals.

**Unattended.** Larry is asleep. Every gate is stop-on-failure: stop, write what you have into the report with the failure at the top, commit, end. **Never ask. Never work around.** A stopped session with an honest report is the correct outcome.

---

## ⚠ CATEGORY

**No `R/` change of any kind.** No new exports, no edits to any package file, driver, or document. Re-running an enumeration with a corrected input is recomputation, not a code change.

Writes: scratch scripts, their `.rds`, and the report, all under `dev/glm-continuous-sims/`; the task document under `dev/tasks/`.

**Compute:** ~2 h if §2 returns (A) — six draw sets at ~5–13 min each plus sweeps and inversions under the order-statistic path. Minutes if it returns (B). No replicates, no simulation study, no renders.

---

## 1. Provenance and first commit — GATE

```bash
hostname; pwd
git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD
git status --porcelain
git log --oneline -6
Rscript -e 'cat("forestsearch", as.character(packageVersion("forestsearch")), "\n")'
```

*GATE:* branch `feature/glm-extension`; HEAD descends from `f2ca69e6`; installed `0.2.6`. A dirty tree is not a failure.

```bash
ls ~/Downloads/cc_task_oc_wrapper_confs_2026*
cp ~/Downloads/cc_task_oc_wrapper_confs_2026*.md dev/tasks/cc_task_oc_wrapper_confs_2026-08-30.md
git add dev/tasks/cc_task_oc_wrapper_confs_2026-08-30.md
git commit -m "task doc — confounder-list mismatch and corrected analytic family"
```

*GATE:* exactly one match (hyphens may be stripped; copy under the name above, report both).

---

## 2. Which covariate set did the search actually run? — GATE

Establish **from the driver source and the tracked payloads**, quoting each piece. Do not infer from a count.

1. **The driver's argument.** The exact `confounders.name` vector passed to `forestsearch()` in each of the three drivers (alt n = 500, alt n = 700, null n = 500), verbatim with its length, including how `include_str2` resolves. Say whether the three agree.
2. **Candidate versus evaluated.** `forestsearch()` returns both `confounders.candidate` and `confounders.evaluated`. From `R/forestsearch_main.R`, what can reduce the first to the second — `cont.cutoff`, `exclude_cuts`, `collapse_cuts`, the GRF pre-screen (`use_grf = FALSE` here), degenerate cuts, anything else. **A 13 → 12 drop could be legitimate; establish whether it is.**
3. **What the payloads say.** Do they record either list (`meta`, or elsewhere)? Quote it, or say they do not.
4. **Evidence from realized rules.** The `sg_def` analysis already shows `str2` in 3–5% of measured selection mass — a variable appearing in a realized winner was unambiguously evaluated. Confirm the count from `sgdef_selection_2026-08-29.rds` rather than recomputing.
5. **What the analytic side used.** The `confounders.name` in `oc_wrapper_grid_2026-08-29.R` and the earlier fixture scripts, verbatim with length; and the cut columns `fs_oc_family_enumerate()` actually built.

**Also check every other enumeration-relevant argument** between the drivers and the analytic reconstruction, one line each, agree or differ: `conf.cont_jcuts`, `cut_type`, `cont.cutoff`, `maxk`, `n.min`, `n.min.frac`, `minp`, `rmin`, `exclude_cuts`, `conf_force`, `collapse_cuts`, `sg_focus`, `effect.threshold`, `consistency.threshold`, `pconsistency.threshold`. A second mismatch found now is far cheaper than one found later.

**Verdict, one of three:**

- **(A) The analytic list was wrong** → §3.
- **(B) The analytic list was right** → skip §3; record the finding, note that the `str2` mass must then have another explanation, and go to §5 using the existing `.rds`.
- **(C) Cannot determine** → stop, report, commit, end.

---

## 3. Correct and re-run — only under (A)

New script `dev/glm-continuous-sims/oc_wrapper_grid_corrected_2026-08-30.R`; leave the 08-29 scripts untouched as the record of what was run before. Rebuild with the **drivers'** arguments, **same seeds as the 08-29 run (20260825)** so the candidate space is the only change.

Into one new `.rds`:

- **Alternative** n ∈ {500, 700} × gate ∈ {resample, split} — full `fs_oc_predict()` output with MC SEs, `M`, enumeration stage counts, floor, settings, seed, wall-clock.
- **Null** n = 500 × both gates — false declaration and MC SE, per-candidate rate range, `L_eff`, `E|Ĥ|`, `E[β(Ĥ)]`, naive bias, `mass_below`, `M`.
- Measured columns recomputed from `payload$oc` / `payload$results` — computed, never typed.
- The `c1` sweep and inversions at the driver's `c2`, both gates, each cell, under the order-statistic path.

**Report every corrected table beside the superseded one**, with the delta per row and the new `M` beside 1601 / 1784. Show both; do not quietly replace.

*GATE:* the fidelity harness and the frozen-reference check must still pass — this changes an input, not code. Run them and say so.

---

## 4. Does the correction close the residuals? — the point of the task

Using the corrected `.rds` and the stored `sgdef_selection_2026-08-29.rds`:

1. **`str2` mass.** Corrected analytic selection mass on `str2` signatures, beside the measured 3–5%. Do they agree?
2. **The between-rule size gap.** Recompute the `sg_def` report's decisive comparison on the corrected family: measured-frequency-weighted population size against corrected-analytic-`sel_c`-weighted population size, per cell. The 08-29 figures were +2.0 / +0.6 / +1.5 subjects. Report the corrected gaps beside them.
3. **`E|Ĥ|` against measurement.** Corrected analytic against measured (72.3 ± 0.4 at n = 500; the corresponding n = 700 and null figures). Closed, partly closed, or unchanged?
4. **Naive bias.** Corrected analytic against measured, all three cells, with the sign. The `sg_def` report showed it tracks the winner's size and purity, so a corrected mix should move it; say by how much and whether the n = 500 / n = 700 sign flip survives.
5. **PPV and sensitivity.** Corrected analytic against measured. The within-rule component (+0.0047 at n = 500) is not expected to move; report what actually happens.
6. **The distribution comparison**, redone on the corrected family: top 15 analytic signatures beside measured frequencies; analytic mass on never-selected signatures; measured mass still absent from the analytic family and why; and the total variation distance beside the 08-29 values (0.17–0.20 against a 0.13–0.14 multinomial noise floor).

**State plainly, per residual, whether the confounder correction accounts for it.** A partial or mixed answer is a real outcome — report it as such, do not force one story.

---

## 5. Report

`dev/glm-continuous-sims/REPORT_oc_wrapper_confs_2026-08-30.md`:

0. Provenance (raw).
1. §2's evidence and verdict (A/B/C), plus the full argument-by-argument comparison.
2. §3's corrected tables beside the superseded ones, with deltas and both `M`; the fidelity and frozen-reference results.
3. §4, item by item, with the per-residual verdict.
4. **What is superseded**, listed explicitly: which numbers in `REPORT_oc_wrapper_{fixture_run,gate_and_n700,grid,null}_*.md`, in `REPORT_oc_wrapper_sgdef_2026-08-29.md`, in `oc_wrapper_verification.qmd`, and in `HANDOFF_oc_wrapper_2026-08-29_v1.md` are replaced, and which stand. Larry needs this list to know what he can still quote.
5. Ten-line verdict.

Close out: `git status --porcelain`; stage by explicit path; commit; `git log --oneline -3`; porcelain again. **No push.** No `devtools::install()` — nothing in `R/` changed.

---

## 6. Out of scope

- No `R/` change, no edits to any package file, driver, or document. **The verification `.qmd` is not re-pointed in this session** — it keeps reading the 08-29 `.rds`; §5 item 4 records that it needs updating later.
- No re-derivation of the `sg_def` tabulation, the population evaluation, or the within-rule comparison — those stand; read `sgdef_selection_2026-08-29.rds`.
- No simulations, no replicate runs, no renders, no push.
- No second family constructor, no binary/OR work, no breadth testing on another DGM.
- No recommendation about what to do with the result.
