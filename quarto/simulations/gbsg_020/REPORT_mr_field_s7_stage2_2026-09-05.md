# REPORT — MR (field) on Section 7 cells: Stage 2 (Run) — STOPPED, batch failure

**Task:** `dev/tasks/TASK_mr_field_section7_2026-09-05.md` (bb516569). Prior records: Stage 0 ca2996b3; Gate 1 STOP 8293a3ce; fix a6702fd8; Gate 1 PASS 7f94fdc9.
**Date:** 2026-09-05. Executor: Claude Code, under the F7 pre-authorization ("a failed batch stops and reports; no mid-run changes of scope or settings").

---

## STAGE 2 STOPPED — batch h100 1–1,000 failed its own FB-join identity gate; cause: 2 of 1,000 committed h10 replicates flip their *selection* between package vintages 0.2.0 and 0.3.5

**What happened.** Stage 2 started at 100 workers (per the Stage 1 calibration), campaign `s7`, batch order h100 1–1000 (fb join) → h100 1001–2000 → h175 1–1000 → h175 1001–2000 → combines. The first batch ran 1,000 replicates and then died, by design, in the template's join assertion:

> `fb_mode = 'join': naive-estimate identity FAILED on 1 of 329 joined detected rows (max rel diff 3.957e-02)`

Fail-fast held: no `s7` bundle was written, nothing downstream ran, and no settings were changed.

**Diagnosis (fresh diag-campaign renders of the committed ranges, `fb_mode = none`; bundles listed below).**

| Range | Divergent replicates vs the committed bundle | Detection flips |
|---|---|---|
| h100 sims 1–1,000 (vs `fs_maxeffCons_mr_m1_h10_..._res_1_1000.rds`) | **2** — sims **393** and **780** | 0 |
| h175 sims 1–1,000 (vs `fs_maxeffCons_mr_m1_h18_..._res_1_1000.rds`) | **0** — exact on every column checked | 0 |

On the other 998 h10 replicates every checked column (naive, oracle, memberships, MR (IJ), β(Ĥ)) is identical to ≤ 1e-12 (oracle columns exact everywhere, including sims 393/780 — same data, so these are pure **selection flips**, not estimation changes):

- **sim 393:** committed (0.2.0) selected `{er <= 121} & {pgr <= 0}` (n = 62, naive HR 1.397); fresh (0.3.5) selects the near-tie parent `{pgr <= 0}` (n = 64, naive HR 1.342, winner Pcons 0.92). Under maxeffCons the higher-HR two-term rule must have fallen just below the 0.90 consistency screen in the fresh run — a knife-edge screen flip.
- **sim 780:** committed `{pgr <= 38} & {er <= 0}` vs fresh `{pgr <= 0} & {er <= 17}` — different near-tie rule, same size (n = 62 both).

**Reproduction note (mechanism confirmed):** a single-replicate refit reproduces the fresh selection only under `RNGkind("L'Ecuyer-CMRG")` — the kind in force inside `future` workers (`seed = TRUE`), which `set.seed(seedit)` then seeds (the known current-kind gotcha). Both vintages ran under that same regime, so the flip is a genuine 0.2.0→0.3.5 behavior difference on knife-edge candidates — most plausibly RNG-stream or evaluation-order sensitivity in the consistency screen introduced by the 0.3.4 effect-screen reordering (3921ffdd) or the 0.3.5 subset-vector fits (f3975b99); both were verified bit-identical on the *applications*, which contain no knife-edge replicate like these. The 0.3.x fit object retains only the winner row, so the displaced rule's fresh Pcons is not recoverable from the fit.

**Why this stops Stage 2 rather than being routed around.** Gate 2 requires "on sim_id 1–1,000 naive / oracle / detection / MR (IJ) identical to the committed bundles" — MR (IJ) est differs on sims 393/780, which is Larry's explicit STOP ("a difference beyond tolerance in MR (IJ) est/SE under the current package is a STOP"). Skipping or special-casing the two replicates, or relaxing the join gate, is a mid-run change of settings, which F7 forbids.

## What is established despite the stop

- The fresh 0.3.5 pipeline reproduces the committed Section 7 results on **1,998 of 2,000** committed replicates exactly (≤ 1e-12, mostly exact-0), including the **entire h175 cell**, and on 100% of detection flags.
- The 2 divergences are enumerated, characterized as near-tie selection flips in the borderline-null cell, and are estimation-invariant (oracle identical).
- The template, fix, and FB-join machinery all behaved exactly as designed — the join gate caught a real cross-vintage discrepancy at replicate granularity.

## Options for Larry (not actioned)

1. **Accept the enumerated flips**: amend Gate 2 to "identical except on the enumerated vintage-flip replicates (393, 780), each shown to be a knife-edge selection near-tie with identical data (oracle exact)"; FB join then drops those sims (join count 329 → 327 expected on 1–500 — sim 393 is in the joined range; 780 > 500 affects no join) with the exclusion recorded. Re-issue Stage 2 with an explicit amendment (a one-knob join-tolerance/skip-list change to the template, called out as such).
2. **Bisect the vintage flip** (0.3.3 → 3921ffdd → f3975b99 → …) to attribute it precisely before deciding; costly, and attribution alone does not restore identity.
3. **Re-baseline on 0.3.5**: treat the fresh runs as the reference for this study; FB joins only where the per-row naive identity holds (the template's gate, demoted from stop to filter-and-report). Same amendment mechanics as (1).

## Evidence left untracked beside the documents (diag campaign, fb = none)

- `results/fs_maxeffCons_fb_mr_field_m1_h100_knoise0_n500_diag_res_1_500.rds`, `..._res_501_1000.rds`
- `results/fs_maxeffCons_fb_mr_field_m1_h175_knoise0_n500_diag_res_1_1000.rds`
- Renders `diag_h100.html`, `diag_h100_b.html`, `diag_h175.html`; driver/census logs in the session scratchpad (`stage2_driver.sh`, `s2_h100_b1.log`, `diag_*.log`).

No `s7`-campaign bundle exists (the failed batch stopped before its save). Committed simulation documents and bundles remain untouched; the only `R/` change in this arc remains a6702fd8.

**Stopped per the F7 fail-fast rule.** Amendment decision required from Larry before Stage 2 can be re-issued.
