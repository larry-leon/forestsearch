# CC task — signed orientation for `fs_oc_family_enumerate()`: opposite-sign region effects supported, same-sign cases bit-identical

**File:** `dev/tasks/cc_task_oc_signed_orientation_2026-08-31.md` · **Issued:** 2026-08-31 by chat, commissioned by Larry (pre-authorized 2026-08-31, conditional on stage 0)
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Follows:** `REPORT_oc_applied_stage0_2026-08-31.md` — read its §3 verdict first.

**Premise gate.** This task runs **only if stage 0's report says the sign gate binds** (opposite-sign `m_tau[Q]` / `m_tau[Qc]` on the applied ACTG175 continuous DGM). If the report says it does not bind, copy and commit this document, write a three-line report saying the task is not needed and why, and **stop** — no edit.

**What and why.** `fs_oc_family_enumerate()`'s alternative branch reads the scale through oriented absolute values — `tauQc = |m_tau[Qc]|`, `bint = |m_tau[Q]| − |m_tau[Qc]|`, `beta_g = tauQc + bint·PQg` — and `stop()`s when the two regions' effects differ in sign, because the absolute-value reading is undefined there. That reading is a fixture-ism: the general mixture is the **signed** form. With `s = sign(m_tau[Q])` — the harm direction, so the planted effect is oriented positive —

    beta_g = s · ( m_tau[Qc] + (m_tau[Q] − m_tau[Qc]) · PQg )

which is **algebraically identical** to the current computation whenever `sign(m_tau[Qc]) == s` (both effects the same direction: `s·m_Qc = |m_Qc|`, `s·(m_Q − m_Qc) = |m_Q| − |m_Qc|`), and simply defined when the signs differ (benefit-region candidates get oriented-negative means, which is correct — they should rarely clear a positive `c1`). Everything downstream consumes `beta_g` as a vector of means; nothing else changes.

---

## ⚠ `R/` CALLOUT — **CHANGES BEHAVIOUR** (extends the domain; the existing domain is bit-identical)

Editable: `R/fs_oc_family.R` (the alternative branch's scale reading and its roxygen), **one new test file** `tests/testthat/test-fs-oc-signed.R`, `man/`, `NAMESPACE` as regenerated, `DESCRIPTION`, `NEWS.md`. **Nothing else** — in particular not `R/fs_oc_predict.R`, `R/fs_oc_grid.R`, `R/fs_dgm_scale.R`, the null branch, or any driver or document.

**Compute:** the identity gates (one MD40 enumeration, small synthetic fixtures), `devtools::test()`, `R CMD check`. **Estimate 30–45 minutes.** No simulation study, no renders.

**Unattended.** Gates stop the task; never ask, never work around.

---

## 1. Provenance and first commit — GATE

Standard block plus the installed version (expect **0.3.1**). *GATE:* branch `feature/glm-extension`, clean tree, stage 0's report commit in the log. Read stage 0's §3 verdict; apply the premise gate above. Copy this document into `dev/tasks/` and commit it alone; report both filenames. Run `devtools::test()` for the parity baseline (expect WARN 31).

## 2. Source verification — GATE

Quote from the current `R/fs_oc_family.R`: the same-sign `stop()` and the `tauQc` / `bint` / `seQ1000` block; the `beta_g` construction at the interface-fields step; and — **grep the downstream consumers** (`fs_oc_predict.R`, `fs_oc_grid.R`) for any `abs()` or positivity assumption applied to `beta_g`, `EbetaH`, or the gate. *GATE:* the downstream must consume `beta_g` sign-agnostically (draws `N(beta_g, S_g)`, thresholds compared directly). If any downstream site assumes positivity, **stop and report** — the edit would then be bigger than authorized.

## 3. The edit

In the alternative branch:

- delete the same-sign `stop()`;
- compute `s <- sign(reg$m_tau[iQ])`; if `s == 0`, `stop()` with a message directing the caller to a nonzero planted effect or the null path (a planted effect of exactly zero has no harm direction to orient by);
- `beta_g <- s * (m_Qc + (m_Q − m_Qc) * PQg)` with the **signed** `m_tau` values; store `s`, the signed `m_tau` pair, and oriented `tauQc = s·m_Qc`, `bint = s·(m_Q − m_Qc)` in the returned object so provenance is explicit (note in roxygen that `tauQc` may now be negative for opposite-sign families);
- `seQ1000` unchanged (it is sign-free, from `V_eff[Q]`);
- roxygen for the family object and the enumerate function: orientation is the harm direction `sign(m_tau[Q])`; opposite-sign families are supported; benefit-direction candidates carry oriented-negative `beta_g`.

Quote the full diff in the report.

## 4. The identity gates — the existing domain must not move

1. **MD40, bit-for-bit.** Rebuild the MD40 DGM exactly as the breadth ladder's §2 did (the direct `generate_glm_dgm(k_inter = k40)` route) and enumerate at n = 500. *GATE:* `identical()` to the corrected run's stored family on every stored field (`lab`, `Pg`, `PQg`, `beta_g`, `se_g`, `sens_g`, `spec_g`, `M`).
2. **The standing guards.** Re-run the frozen-reference and fidelity checks the workstream always runs after an `R/` edit. *GATE:* pass.
3. **A same-sign synthetic** (small DGM, both effects negative): enumerate pre-edit (git stash or temp lib — say which) and post-edit. *GATE:* `identical()`.
4. **The new domain.** Enumerate the applied opposite-sign DGM from stage 0's §3 at one rung (`q = 20`): no stop; report M (must equal stage 0's §4 M — same memberships), the `beta_g` range, and an **independent recomputation** of `beta_g` for five candidates spanning purity 0–1, directly from `df_super` memberships and the signed mixture formula, agreeing to 1e-12. Also confirm the coherence direction: the candidate nearest Ĥ has the largest `beta_g`, and pure-benefit candidates are negative.

## 5. Tests

`tests/testthat/test-fs-oc-signed.R`, small and fast: (a) same-sign fixture — post-edit `identical()` to a stored reference built in the test itself pre-edit-equivalent (construct the expected values from the abs formula); (b) opposite-sign fixture — enumerate succeeds, `beta_g` matches the signed mixture computed directly, monotone in `PQg`; (c) `s == 0` stops with the new message; (d) the null path untouched (one null enumeration unchanged in structure). Per the v5 §9 principle, show one test failing against the defect it guards (revert the edit for (b), confirm the old `stop()` fails it, restore).

## 6. Close-out

- `devtools::document()`; `devtools::test()` — *GATE:* FAIL 0; warning parity 31 = 31 apart from the new file's own additions (report its counts); any decrease or unrelated increase is a stop.
- `R CMD check` (`RSTUDIO_PANDOC` on Rscript's PATH). *GATE:* 0/0/0, status verbatim.
- `DESCRIPTION` → **0.3.2**. `NEWS.md`: opposite-sign region effects supported in `fs_oc_family_enumerate()` via signed orientation `s = sign(m_tau[Q])`; same-sign families bit-identical (the MD40 gate); `tauQc` in the family object may be negative for opposite-sign families; a planted effect of exactly zero is rejected with guidance.
- Commit by explicit path; **no push**; `devtools::install()` (`upgrade = FALSE`); confirm **0.3.2**.
- Report `REPORT_oc_signed_orientation_2026-08-31.md`: the premise verdict · §2's quotations and the downstream grep · the diff · §4's four gates with the independent recomputation shown · test/check output raw · commits · ten-line verdict.

## 7. Out of scope

- No change to the null branch, to `se_g`'s construction, to `fs_oc_predict()`/`fs_oc_grid()`/`fs_oc_invert()`, to `fs_dgm_scale()`, or to any driver, application or document.
- No evaluation run — the applied ladder is the next task, gated on Larry's compute go/no-go from stage 0's §5 projection.
- No push, no renders.
