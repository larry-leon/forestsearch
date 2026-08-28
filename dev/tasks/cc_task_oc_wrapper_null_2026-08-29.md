# CC task — the null path: fast inversion, a null branch for the family, and the null-cell comparison

**File:** `dev/tasks/cc_task_oc_wrapper_null_2026-08-29.md` · **Issued:** 2026-08-29 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** Larry places this file in `~/Downloads`; CC copies it into `dev/tasks/` and commits it first (§1).
**Follows:** `REPORT_oc_wrapper_grid_2026-08-29.md` (read first), and the two 08-28 reports.

**Why this closes the workstream's last gap.** The continuous workstream's one unexplained discrepancy is the null cell: the prediction document reports a false-declaration rate on its stylised M = 16 family against a measured 998/1000. The enumerated family closed the alternative-cell gaps; this task asks the same question of the null. It is also the only cell where the two gates can differ, because the alternative saturates.

---

## ⚠ `R/` CALLOUT

Two changes, both to files this workstream created:

1. **`R/fs_oc_grid.R` — changes behaviour of a new function** (§3): the inversion and threshold sweep are replaced by an exact order-statistic reduction. Results must be **unchanged**; the twelve stored slow-path inversions are the guard.
2. **`R/fs_oc_family.R` — adds a null branch** (§4). Additive; the alternative path must be bit-identical afterwards.

**No existing package file may be edited — in particular NOT `R/fs_dgm_scale.R`.** §2 investigates it read-only and reports what a null branch there would take, as a proposal for Larry.

**Compute:** ~25 min. One enumeration and one draw set per gate at n = 500 under the null; the §3 reduction makes the sweep and inversions essentially free. No replicates, no simulation study.

---

## 1. Provenance and first commit — GATE

```bash
hostname; pwd
git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD
git status --porcelain
git log --oneline -6
Rscript -e 'cat("forestsearch", as.character(packageVersion("forestsearch")), "\n")'
```

*GATE:* branch `feature/glm-extension`; HEAD descends from `399ccabe`; installed `0.2.5`.

```bash
ls ~/Downloads/cc_task_oc_wrapper_null*
cp ~/Downloads/cc_task_oc_wrapper_null*.md dev/tasks/cc_task_oc_wrapper_null_2026-08-29.md
git add dev/tasks/cc_task_oc_wrapper_null_2026-08-29.md
git commit -m "task doc — null path: fast inversion, null family branch, null-cell comparison"
```

*GATE:* exactly one match (hyphens may be stripped; copy under the name above, report both). Run `devtools::test()` for the parity baseline.

---

## 2. Read-only — why `fs_dgm_scale()` cannot serve the null, and what it would take

The driver guards it out entirely: `dgm_scale <- if (isTRUE(null_cell)) NULL else fs_dgm_scale(dgm)`, and the null payload carries no `scale` element. Establish from source, quoting it:

1. **What fails.** Read `R/fs_dgm_scale.R`. Does it error under a null DGM, return a malformed object, or silently produce something wrong? Which step depends on a non-empty Q?
2. **What survives.** Are the whole-population (`region == "S"`) quantities — `sigma`, `V_mu0`, `V_tau`, `v_cond0`, `bracket`, `V_eff`, `n_g`, `P_g` — computable when there is no Q/Qc partition? Answer per field.
3. **`V_eff` across regions.** In the alternative-cell scale table the three `V_eff` values differ by well under 1%. State the actual spread from the tracked payload, and say whether `V_eff` is region-invariant by construction or only numerically close here. This decides whether a single whole-population `V_eff` is a sound basis for `se_g` under the null.
4. **Proposal, not an edit.** What a null branch in `fs_dgm_scale()` would take: which lines, what it would return (presumably an S row only), and whether any existing caller could be affected. Write it as a diff in the report. **Do not apply it.**

*GATE:* if (2) shows the whole-population effective variance is **not** computable under the null without editing `R/fs_dgm_scale.R`, stop after §3 — implement the fast inversion, write the report with the proposal, commit, and end. The null branch then needs Larry's word on the `fs_dgm_scale()` change.

---

## 3. `R/fs_oc_grid.R` — the exact order-statistic reduction

Each inversion currently draws its own set (~2.4 h for twelve). It need not. At fixed `(n, gate, c2)` neither the eligible set nor the winner depends on `c1`:

```
resample:  eligible_g = (Bhat_g - c2 >= z_p * se_g)
split:     eligible_g = (W1_g >= c2) & (W2_g >= c2),   Bhat_g = (W1_g + W2_g)/2

per draw:  w = argmax_{g eligible} Bhat_g      (NA if none eligible)
           T = Bhat_w                          (-Inf if none)
```

Declaration at `c1` is exactly `T >= c1`, and `w` is the winner whenever anything declares: if the argmax over eligible falls below `c1` nothing passes, and if it clears `c1` it is also the argmax over passers.

So per draw, store `T` and the winner's statistics once — `Pg[w]`, `PQg[w]`, `sens_g[w]`, `spec_g[w]`, `npv[w]`, `beta_g[w]`, `noise = Bhat_w - beta_g[w]`. Sort draws by `T` descending; the draws declaring at any `c1` are a **prefix** of that order, so every conditional mean is a `cumsum` read off with `findInterval`. `mass_below` remains `c1`-dependent but only through `1{beta_g[w] < c1}` at the stored winner. Inversion becomes an order statistic — `k = round(target * draws)`, `c1` = the k-th largest `T`; the ceiling is `mean(T > -Inf)`, the fraction of draws with any eligible member, with `NA` and the binding threshold named above it.

Keep the comparison operators exactly as they are (`>=` for declaration, `<` for `mass_below`). One pass per `c2` value; a `c2` vector loops.

*GATE:* re-run the twelve stored inversions from `oc_wrapper_grid_2026-08-29.rds` under the reduction. The inverted `c1` values must reproduce the stored ones **to the resolution the draw count supports** — state the comparison you used and the largest discrepancy. A difference beyond one order-statistic step means the reduction is wrong: stop and report. The existing GATE tests (grid/predict identity, monotonicity, round-trip, ceiling) must all still pass.

Report the new wall-clock for the twelve inversions beside the old 2.4 h.

---

## 4. `R/fs_oc_family.R` — the null branch

**Detection first.** Establish from source how a null DGM is reliably identified — candidates include `prevalence_Q == 0`, an absent or all-zero harm flag, `beta_inter == 0`, `is.na(effect_Q)`. Read the null driver and the null payload's `truth` and pick a test that cannot misfire on an alternative DGM. Report what you chose and why. Do not guess a field name.

**Under the null, every candidate has the same true effect**, so the fields become:

| field | null value | note |
|---|---|---|
| `beta_g` | the common effect, for every `g` | from the DGM's own truth (`effect_Qc` / `effect_ITT`); oriented as the alternative path orients |
| `se_g` | from the whole-population effective variance at `(n, Pg)` | same prevalence scaling as the alternative path; §2(3) decides whether this is sound |
| `PQg` | `0` | `P(g ∩ Q) = 0` |
| `sens_g` | `NA_real_` | `0/0` — undefined, **not** zero |
| `spec_g` | `1 - Pg` | `Qc` is the whole population |
| `npv` | `1` | follows from `PQ = 0` |

Enumeration, floors, `ovl`, `Rho`, `Sg` and the root are unchanged — none of them touches Q.

Downstream: `fs_oc_predict()` and `fs_oc_grid()` must propagate `NA` sensitivity without contaminating other quantities (`E[sens]` is `NA`; `E|Ĥ|`, `E[β(Ĥ)]`, naive bias, `mass_below`, `Espec`, `Eppv`, `Enpv` all remain defined). Check this rather than assuming it.

*GATE:* the alternative path is unchanged — re-run the fidelity harness (bit-identical to the document's chunk) and the fixed-seed comparison against the stored reference. Both must still pass.

Tests: a small synthetic null DGM giving `PQg = 0`, `sens_g` all `NA`, `spec_g == 1 - Pg`; declaration monotone in `c1`; and the alternative path untouched.

---

## 5. The null-cell comparison

Rebuild the null DGM deterministically (the driver that produced `mr_md_harm/fs_maxeffCons_mr_mdnull_knoise0_n500_s1000_d5000/fs_maxeffCons_mr_mdnull_knoise0_n500_res_1_1000.rds`, added in `2b180813`), gated against that payload's committed `truth` as the earlier scripts gate against the MD40 truth.

Run at `n = 500`, both gates, the driver's `c1`/`c2`, 2e5 draws, stated seed. Extend `oc_wrapper_grid_2026-08-29.rds` (or write a sibling `.rds`; say which) with:

- false declaration rate and its MC SE, both gates;
- the per-candidate rate range (`min`/`max` of `colMeans(pass)`) and the implied `L_eff = log(1 - fam)/log(1 - max(p1))`, the document's own definition, at `M ≈ 1601` rather than 16;
- `E|Ĥ|`, `E[β(Ĥ)]`, naive bias, `mass_below`;
- `M` and the enumeration stage counts under the null floor;
- the **measured** column from the null payload's `$oc` / `$results` — computed, never typed. The null payload's third element is an `NA`-named `NULL` (a known pre-`d884adbf` artifact); index by name and do not trip on it.

Then a c1 sweep and the inversions at the null, using §3's reduction — under the null this is the type-I error curve, and it is the quantity that decides whether the gates differ anywhere.

**Report explicitly whether the two gates differ here.** They agreed to three decimals at the alternative because the rate saturated. If they still agree, say so plainly; that is a finding, not a failure.

---

## 6. Document

Extend `quarto/simulations/actg175/continuous/oc_wrapper_verification.qmd` with a null section: the null branch's field definitions (the table above), the comparison against the measured null cell, the type-I error curve in the document's existing base-graphics idiom, and two sentences on what the M ≈ 1601 result says about the document's stylised M = 16 figure. Types no number; re-render and report exit code and time.

---

## 7. Close-out

- `devtools::document()`; `devtools::test()` — counts and warning-count parity against §1.
- `R CMD check` (`RSTUDIO_PANDOC` on Rscript's PATH); Status verbatim.
- Version → `0.2.6`; `NEWS.md` naming the null branch and the inversion reduction, and stating that alternative-cell results are unchanged.
- Commit, **explicit paths only**, never `git add -A`. **No push.** `devtools::install()`; confirm version.
- Report `dev/glm-continuous-sims/REPORT_oc_wrapper_null_2026-08-29.md`: §2's source answers and the `fs_dgm_scale()` proposal diff; §3's guard result and new wall-clock; §4's detection choice and the unchanged-alternative guards; §5's tables; the render; commits; ten-line verdict.

---

## 8. Out of scope

- No edit to `R/fs_dgm_scale.R` or any other pre-existing package file — §2 is read-only and its diff is a proposal.
- No binary/OR path, no second family constructor, no changes to the prediction document or any driver.
- No new inversion targets beyond the declaration rate.
- No simulation study, no replicate runs.
- If §2's gate closes, stop after §3 rather than working around the missing scale.
