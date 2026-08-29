# CC task — the two candidate explanations for the remaining residual

**File:** `dev/tasks/cc_task_residual_2026-08-30.md` · **Issued:** 2026-08-30 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Reads:** `REPORT_oc_wrapper_sgdef_2026-08-29.md` + `sgdef_selection_2026-08-29.rds`; `REPORT_oc_wrapper_confs_2026-08-30.md` + `oc_wrapper_grid_corrected_2026-08-30.rds`.

**Where the question stands.** After the confounder correction the analytic prediction still differs from measurement in one direction: `E|Ĥ|` low, PPV and sensitivity low, naive bias off with a sign that flips between n = 500 and n = 700. The `sgdef` report established the gap is **between-rule** — within every covariate signature the realized values match their population values, so it is not a sample-realization effect; the search simply weights *larger* rules than the analytic `sel_c` does, by +2.11 / +0.61 / +1.65 subjects. Correcting the confounder list fixed the `str2` mass and the representability hole but **did not** move that gap.

Two explanations remain. This task tests both. **Neither is chat's to assert** — an earlier chat hypothesis about `se_g` was stated in a form the `str2` evidence did not support, and is recorded here as a hypothesis only.

---

## ⚠ CATEGORY

**No `R/` change of any kind.** No new export, no default change, no edit to any package file, driver or document. Writes: scratch scripts, their `.rds`/`.log`, and the report, under `dev/glm-continuous-sims/`. Plus this task document.

**Compute:** §2 is arithmetic on stored objects — minutes. §3 re-runs `fs_oc_predict()` on modified families, one draw set per cell — estimate 30–60 minutes. No simulation study, no replicates, no renders.

**Unattended.** Gates stop the task; never ask, never work around.

---

## 1. Provenance and first commit — GATE

Standard block plus the installed version (expect **0.3.0**). *GATE:* branch `feature/glm-extension`, clean tree. Copy this document into `dev/tasks/` and commit alone; the `~/Downloads` stem arrives with hyphens stripped, so report both names.

---

## 2. Hypothesis A — is `se_g`'s prevalence scaling wrong for candidates unlike Q?

The wrapper sets `se_g = seQ1000 * sqrt(2) * sqrt(piQ / Pg)`, rescaled to `n`. That takes the **Q region's** effective variance and scales it by prevalence alone, assuming per-subject effective variance is candidate-invariant. The handoff records this as the load-bearing approximation, justified numerically (the alternative cell's three `V_eff` rows spread 0.86%) rather than structurally — `V_eff` carries the within-region variance of baseline means, so it is not region-invariant by construction.

**`fs_dgm_scale()` takes a public `regions` argument** — established in the null task, where it computed every S-row quantity on a null DGM without any edit. Use it.

For each **distinct realized rule** in `sgdef_selection_2026-08-29.rds` (already evaluated on `df_super`; do not re-evaluate), compute at each cell's `n`:

| quantity | how |
|---|---|
| `se_scaled` | the wrapper's value: `seQ1000 * sqrt(2) * sqrt(piQ / Pg)`, rescaled to `n` |
| `se_direct` | from `fs_dgm_scale()` with `regions` defined by **that rule's own membership** — its own `V_eff` at its own `Pg` |
| `ratio` | `se_scaled / se_direct` |

Report the ratio's distribution, and — this is the test — how it moves with the candidate's **prevalence** `Pg` and with its **composition distance from Q** (use purity `PQg`, or the Jaccard with Q; state which). A ratio near 1 with no trend refutes the hypothesis. A ratio drifting systematically — and in particular **below 1 for large, low-purity candidates**, which would mean the wrapper understates their SEs and so under-selects them — supports it.

Also report the same comparison restricted to Q itself and to the highest-purity candidates, as a control: those should agree by construction.

**State the verdict plainly, either way.**

---

## 3. Hypothesis B — does the dedup asymmetry account for the gap?

The search runs near-duplicate removal **unconditionally** — no off switch — and for every focus except `maxeff` it is the **statistics-keyed** `remove_near_duplicate_subgroups()`. The analytic family deliberately excludes it, collapsing only *identical* population membership. So the two families differ in composition by construction.

Read `remove_near_duplicate_subgroups()` and state the exact key it uses. Then build an analytic analogue: collapse candidates in the corrected family whose **population** counterparts of those same statistics match at the same tolerance (e.g. `Pg`, `beta_g`, `K`), keeping the same representative the package keeps. Report how many candidates it removes per cell.

Re-run `fs_oc_predict()` on the deduped family at each cell — n ∈ {500, 700} alternative and the n = 500 null, both gates, **same seeds as the corrected run (20260825)** — and report:

- `M` before and after the dedup;
- the between-rule size gap, against the corrected run's +2.11 / +0.61 / +1.65;
- `E|Ĥ|`, PPV, sensitivity, naive bias, against measurement;
- whether the direction and magnitude move toward measurement, away, or not at all.

If the analogue cannot be built faithfully — the package's key uses fitted statistics the analytic family has no counterpart for — **say so, build the closest defensible version, state exactly how it differs, and label the result approximate.** Do not present an approximation as a test of the real mechanism.

---

## 4. Verdict

State, separately and without hedging:

1. Hypothesis A — supported, refuted, or indeterminate, with the ratio distribution and its trend.
2. Hypothesis B — how much of the +2.11 / +0.61 / +1.65 gap it accounts for.
3. Whether either, both, or neither explains the residual — **and if neither does, say so.** "Both refuted, the residual is unexplained" is a valid and useful outcome; do not manufacture a cause.

No recommendation.

---

## 5. Report

`REPORT_residual_2026-08-30.md`: provenance · §2's table, distribution and trend · §3's key from source, the analogue's construction and its fidelity caveats, the re-run comparison · the three verdicts · ten-line summary. Commit scripts, outputs and report by explicit path. **No push.** No `devtools::install()`.

---

## 6. Out of scope

- No `R/` change, no edit to any driver, application or document. `oc_wrapper_verification.qmd` is not re-pointed here.
- No second family constructor, no breadth testing on another DGM, no binary/OR work.
- No re-derivation of the `sgdef` tabulation or population evaluation — read the stored `.rds`.
- No simulations, no replicate runs, no renders, no push.
