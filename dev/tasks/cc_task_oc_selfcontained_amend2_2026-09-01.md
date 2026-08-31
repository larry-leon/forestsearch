# CC task — amendment 2 to the self-contained applied OC document: each Q-variant gets its own sub-threshold null

**File:** `dev/tasks/cc_task_oc_selfcontained_amend2_2026-09-01.md` · **Issued:** 2026-09-01 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Amends:** `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd` (commit `5ba546b3`), per chat's review of that render. One defect and three hygiene items. The defect is chat's spec error, not CC's implementation — the previous task's wording asked for exactly what was built.

---

## ⚠ CATEGORY

**No `R/` change.** Edits: the one document above, its rendered HTML, its payload. Plus this task document. The archived `analysis_actg175_continuous_oc_evaluation.qmd` is **not touched**.

**Compute:** two additional (variant, rung) jobs — 14 in all, plus the homogeneous grid — inside a re-render. At the measured ≈ 10 min per job and ≈ 5 GB each, one wave of 14 at `n_workers = 14` projects **≈ 40 min, ≈ 75 GB peak** (the 12-wide run measured 62.8 GB on a 251 GB host). *GATE:* if the projection exceeds **1.5 h**, stop and report.

**Unattended.** Gates stop the task; never ask, never work around. **Provenance gate is dirt-tolerant**; stop only if dirt touches this task's own paths.

---

## 1. Provenance and first commit — GATE

Standard block plus the installed version (expect **0.3.2**). *GATE:* branch `feature/glm-extension`; `59166fc2`, `5ba546b3`, `8b4e31c7` in the log. Copy this document into `dev/tasks/` and commit it alone; report both filenames.

## 2. The defect

The document asserts that the `q = 0.01` sub-threshold null is **shared across Q-variants** and uses the primary's grid at that rung as the left endpoint of every variant's calibration curve (`knob-table` chunk: `d <- if (qv == q_rungs[1]) ocp else oc[oc$variant == v, ]`). **It is not shared.** At `q → 0` a candidate's true mean is

    beta_g^v = beta_treat + (q − beta_treat)·PQg^v  ≈  beta_treat·(1 − PQg^v)

which depends on the variant through `PQg^v`. A candidate largely inside `Q_broad` but outside `Q_primary` has `PQg^b ≈ 1` and `PQg^p ≈ 0`, so its mean is ≈ 0 under broad against ≈ `beta_treat` under primary — a ~27-unit difference, not 0.01.

The `null-shared` chunk's check, `max(abs(f0$beta_g − beta_treat * (1 − f0$PQg)))`, compares each family's `beta_g` to **its own Q-dependent background**. That is an algebraic identity whose deviation equals `q·PQg ≤ 0.01` by construction, so it passes for every variant and cannot detect the between-variant difference. It verifies the wrong thing.

**Scope:** the primary curve, type-I, the ladder and every power column are unaffected. Only the **wider and broad calibration curves' left endpoints**, and hence their `q at 0.05` / `q at 0.5` crossings, are wrong. The bias is conservative — the primary's lower rate understates the supersets at the low end, pushing their crossings to higher `q` and *understating* the breadth–severity trade.

## 3. Fix — each variant runs its own sub-threshold null

- `q_shared <- c(0.01, 20, 40, T_obs)` for the supersets (keep `q_rungs` as it is for the primary), so `jobs` becomes 6 + 4 + 4 = **14**; `n_workers` default → **14**.
- In the `knob-table` chunk, delete the shared-null branch: every variant reads its own rows, `d <- oc[oc$variant == v, ]`, for every `q` including `0.01`.
- **Replace the `null-shared` chunk with a `null-not-shared` chunk** that reports the fact rather than the identity: for each variant's `q = 0.01` family, the max and median of `|beta_g^v − beta_g^primary|` across the M members, and each variant's own declaration rate at `(T̂_obs, c2_ratio·10)`. *GATE:* the between-variant max exceeds 1 CD4 unit for at least one superset (if it does not, the premise of this amendment is wrong — stop and report). Two sentences of prose: the sub-threshold null is a property of the planted Q, so each variant carries its own; the earlier shared-null shortcut was invalid and is removed.
- The `knob-figure` caption and the §10 prose lose the "sharing the sub-threshold null point" claim.
- One sentence in §7.1 noting that each superset's sub-threshold null has its own type-I, that the values are in the payload, and that §7.1's headline is the primary's.

## 4. Hygiene

1. **Payload trap.** `fam_hom` inherits `PQg`, `sens_g`, `spec_g` from the primary family, so the homogeneous rows' `Eppv`, `Esens`, `Espec`, `Enpv` are computed against a Q that does not exist under that truth. Set those four columns to `NA_real_` on the homogeneous rows **before** they enter `oc` and the payload, with a one-line comment saying why. `EnH`, `EbetaH` and `mass_below` depend only on `Pg` and `beta_g` and stay as computed.
2. **Per-rung truth gate.** In `eval_job`, add `abs(abs(sc$m_tau[Q-row]) − q) < 1e-9` (read from `fs_dgm_scale()` on that job's DGM, or from the family's stored `orientation$m_tau_Q`) alongside the existing `M` and `s` assertions, so a mis-set `k_inter` fails loudly rather than being caught only by its sign.
3. **Inert argument.** `adverse_outcome = TRUE` in `dgm_at()` has no effect for a continuous outcome (the negation branch in `generate_glm_dgm()` is binary-only). Drop it, or keep it with a comment saying it is inert here — state which you did.

## 5. Render and gates

Re-render, recording wall-clock and peak memory. *GATES:* §3's between-variant gate; all existing gates (anchor, `M = 4508`, `s = 1`, Q-nesting, null monotonicity pointwise, ladder monotonicity); HTML free of `NA ±` and error text in the body; payload written. Confirm the primary column's crossings are unchanged from the previous render (they use the same inputs) — report them beside the corrected wider/broad values.

## 6. Close-out

Commit by explicit path (qmd, HTML, payload, report). **No push. No install. No `R/` change.** Report `REPORT_oc_selfcontained_amend2_2026-09-01.md`: provenance · §2 restated with the measured between-variant differences · the corrected knob table's four crossing columns beside the previous render's · the primary-unchanged check · §4's three items · render wall-clock and memory · commits · ten-line summary.

## 7. Out of scope

- No `R/` change; no new Q-variants, rungs or draw-count change beyond §3's two added jobs; no change to the archived evaluation; no gbsg; no push.
