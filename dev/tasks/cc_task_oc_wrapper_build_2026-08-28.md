# CC task — build `fs_oc_predict()` and the population-enumeration family constructor

**File:** `dev/tasks/cc_task_oc_wrapper_build_2026-08-28.md` · **Issued:** 2026-08-28 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** Larry places this file in `~/Downloads`; CC copies it into `dev/tasks/` and commits it as its first action (§1).
**Evidence base:** `dev/glm-continuous-sims/REPORT_oc_wrapper_inventory_2026-08-28.md` (Parts A/B/C) and `dev/tasks/cc_task_oc_wrapper_inventory_2026-08-28.md`. Read the report's A3 S1–S6 before writing code.

---

## ⚠ THIS TASK TOUCHES `R/` — category, stated plainly

**It ADDS new files to `R/`. It moves no existing package code, changes no existing function's behaviour, and changes no method.**

- No existing file under `R/` may be edited. Not one line. If you believe an existing function must change, **stop and report** — that is a proposal for Larry, not part of this task.
- The wrapper lives inside the package, so it calls internal helpers directly. **Export nothing new except the two functions named in §3 and §4.** Do not add `@export` to any existing internal.
- `DESCRIPTION` (version), `NAMESPACE` (regenerated), `NEWS.md` and `man/` are the only other files that change.

**Compute: none.** No simulations, no replicates, no renders. The only execution is `devtools::document()`, `devtools::test()`, `R CMD check`, and the fidelity harness of §5, which draws from normal distributions and runs in seconds.

---

## 1. Provenance and first commit — GATE

```bash
hostname; pwd
git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD
git status --porcelain
git log --oneline -4
Rscript -e 'cat("forestsearch", as.character(packageVersion("forestsearch")), "\n")'
```

*GATE:* branch `feature/glm-extension`; HEAD descends from `c016df02` (the inventory report commit). Dirty tree is not a failure — list it, leave it alone, work around it.

Then, before anything else:

```bash
ls ~/Downloads/cc_task_oc_wrapper_build*
cp ~/Downloads/cc_task_oc_wrapper_build*.md dev/tasks/cc_task_oc_wrapper_build_2026-08-28.md
git add dev/tasks/cc_task_oc_wrapper_build_2026-08-28.md
git commit -m "task doc — build fs_oc_predict() and the population-enumeration family constructor"
```

*GATE:* exactly one file matches in `~/Downloads` (the stem may have hyphens stripped; copy under the name above regardless and report both).

Also run `devtools::test()` **now**, before any edit, and record the raw counts — this is the parity baseline for §7.

---

## 2. Read first — the `sigma_D` question — GATE

Read `R/consistency_resample.R` in full, plus its call site in `R/subgroup_consistency_main.R`. Establish **from source**, quoting it:

1. The exact rate formula as implemented, the comparison scale `c` it uses, and how `sigma_D` is computed.
2. Whether `sigma_D` is a fixed function of the subgroup's fitted effect SE — i.e. whether `sigma_D = k * SE(beta_hat_g)` for a constant or a simply-derivable factor `k` — or whether it depends on sample quantities (the dfbeta vector) with no population expression at `(n, Pg)`.
3. Whether `fs.splits` is used at all on the `"resample"` path.
4. Whether two-stage screening (`use_twostage = TRUE`, `twostage_args`, `n.splits.screen`) consumes RNG on the `"resample"` path.

*GATE — this decides the scope of §4:*

- **If (2) yields an exact population expression**, implement both gates in §4 and report the derivation with its source lines.
- **If it does not**, implement the `"split"` gate only, leave `consistency_method = "resample"` **unimplemented with an informative `stop()`**, and record in the report exactly what `sigma_D` would need. Choosing an approximation is a *method* change and is Larry's call, not yours. Do not invent one; do not proceed past this point on the `"resample"` path.

Either way, continue with the rest of the task.

---

## 3. `R/fs_oc_family.R` — the population enumeration

New exported function. Suggested signature (adjust names to match existing package conventions, and say so if you do):

```r
fs_oc_family_enumerate(dgm, forestsearch_args, n, max_M = 2000L, verbose = FALSE)
```

**What it does**, per the decision record's family/gate boundary:

1. Take the DGM's population frame (`dgm$df_super`) and Q's true membership indicator on it. Establish from source how the harm indicator is carried on the DGM object — do not guess the field name.
2. Build the candidate cut columns by calling the package's **own** cut machinery (`get_FSdata()` / `cut_var_jq()`) on `df_super`, driven by the cut-related entries of `forestsearch_args`: `confounders.name`, `conf.cont_jcuts`, `cut_type`, `cont.cutoff`, `conf.cont_medians`, `conf.cont_medians_force`, `conf_force`, `defaultcut_names`, `exclude_cuts`, `collapse_cuts`, `collapse_cuts_args`. Both directions of each cut, as the search does. **Cuts land at population quantiles because the frame is `df_super`** — that is the whole point; do not resample.
3. Enumerate all combinations up to `maxk` using the package's own combination machinery (`get_combinations_info()` / `generate_combination_indices()`), exactly as `forestsearch_main.R:3264–3277` composes it for the MR family.
4. Apply the **structural** floors only, all on the population frame: `minp` per constituent factor; the `rmin` redundancy rule; and subgroup size as `Pg >= n.min/n` (establish from source how `n.min` and `n.min.frac` combine, and use that rule).
5. **Do not** apply the GRF pre-screen, the effect screen, the consistency screen, or statistics-keyed near-duplicate removal. Collapse candidates whose population membership vectors are *identical*; keep the rest.
6. Compute and return the nine interface fields, plus provenance:

| field | definition |
|---|---|
| `lab` | candidate label (the rule as a readable expression) |
| `Pg` | `mean(memb_g)` on `df_super` |
| `PQg` | `P(g ∩ Q) / P(g)` |
| `beta_g` | mixture mean from the scale table: `tauQc + bint * PQg` |
| `se_g` | scale-table se at this `n` and `Pg` — derive it the way the document does (A3 S1(b), `se_g <- seQ1000 * sqrt(2) * sqrt(piQ / Pg)`, with the n-rescaling of `appendix-tier2` L2324); take every input from `fs_dgm_scale(dgm)`, never typed |
| `sens_g` | `P(g ∩ Q) / P(Q)` |
| `spec_g` | `1 − P(g ∩ Qc) / P(Qc)` |
| `ovl` | M×M matrix of `P(g_i ∩ g_j)` |
| `M` | `length(lab)` |

plus `scale` (the `fs_dgm_scale` object used), `n`, and the `forestsearch_args` subset actually consumed. Give the returned object a class (e.g. `fs_oc_family`) and a `print()` method reporting M, the prevalence range and the floor applied.

**Size guard:** report M after filtering. If `M > max_M`, `stop()` with the count and the covariance cost before allocating an M×M matrix — the enumerated family is much larger than the document's 16 and `fs_sym_root()` is O(M³).

---

## 4. `R/fs_oc_predict.R` — the assembly

New exported function. Suggested signature:

```r
fs_oc_predict(dgm, forestsearch_args, n, c1 = NULL, c2 = NULL,
              family = NULL, consistency_method = c("resample", "split"),
              draws = 2e5, seed = NULL, ...)
```

- `c1`/`c2` default to `forestsearch_args$effect.threshold` / `$consistency.threshold`; **passing them explicitly overrides**, which is what makes a grid and the inversion possible. `n` likewise overrides any n implied by the args.
- `family = NULL` calls `fs_oc_family_enumerate()`; a supplied `fs_oc_family` object is used as-is. This is the pluggable seam — do not build any other constructor in this task.
- `consistency_method` mirrors the package's own argument, with `"resample"` (closed form) first, subject to §2's gate.

**The body is `predict_scn()` (`worked-sensitivity` L1972–2003) with its family construction removed.** Port it faithfully — read the report's A3 S2–S6 and the document's `worked-predictions` chunk, and keep the same operations in the same order:

`Rho`/`Sg` from `ovl` and `se_g` → `fs_sym_root(Sg, scale = 2)` → `W1`, `W2` → `Bhat = (W1+W2)/2` → the gate → `det`/`det_rate` → `Bmask`/`max.col` argmax (maxeffCons) → `P1`, `p_sel`, `sel_c` → `EnH`, `Esens`, `Espec`, `Eppv`, `Enpv`, `EbetaH`, `Enaive_bias`, `mass_below`.

Gate forms:
- `"split"`: `(W1 + W2 >= 2*c1) & (W1 >= c2) & (W2 >= c2)` — the document's, verbatim.
- `"resample"`: the closed-form screen established in §2, applied to `Bhat`, combined with the effect screen at `c1` — **only if §2's gate opened.**

Return a classed list: the operating characteristics, their Monte-Carlo SEs where they are proportions (reuse the document's `sqrt(p*(1-p)/R)` form), the family's M, the settings used, and the seed. Add a `print()` method.

**Do not** implement a grid, an inversion, a null-DGM path, or a binary/OR path in this task. Note in the report what each would need.

---

## 5. Fidelity harness — the gate that makes the extraction trustworthy — GATE

An extraction is only correct if it reproduces what it was extracted from.

1. In a scratch script under `dev/glm-continuous-sims/` (not in `R/`, not in `tests/`), rebuild the document's **M=16 scenario family** by running the `worked-scenario` chunk's own code from `quarto/simulations/actg175/continuous/analytic_verification_and_prediction_md_harm.qmd` — the nine fields, unchanged.
2. Capture the document's own `worked-predictions` outputs at `Rdraw = 2e5`, `seed = 20260825` by running that chunk's code in the same scratch session. **These are the reference values; read them from the document's code, never from a handoff, a report, or this task.**
3. Call `fs_oc_predict(family = <that family>, n = 500, c1 = 30, c2 = 10, consistency_method = "split", draws = 2e5, seed = 20260825)`.

*GATE:* every returned quantity is **bit-identical** to the chunk's — `identical()` on the numeric vector, not `all.equal()`. If they differ, the port is not faithful: stop, report the first differing quantity with both values and your diagnosis, commit what you have, and end. Do not "fix" by adjusting tolerances.

If the draws cannot be made bit-identical because the RNG consumption order differs, say so explicitly, state the order in each, and report the largest absolute difference — but treat that as a **failure of this gate**, not a pass.

Keep the scratch script; commit it beside the report as the fidelity record.

---

## 6. Tests

`tests/testthat/test-fs-oc-predict.R`, all small and fast (`skip_on_cran()` where a fixture is slow):

1. **Interface contract:** a hand-built minimal `fs_oc_family` (M = 3, simple overlaps) runs through `fs_oc_predict()` and returns all named quantities, finite, with `det_rate` in [0,1] and `sel_c` summing to 1.
2. **Family construction:** on a small DGM, `fs_oc_family_enumerate()` returns `Pg` equal to direct `mean()` of the membership it reports, `ovl` symmetric with `diag(ovl) == Pg`, `PQg` in [0,1], and every kept candidate satisfying `Pg >= n.min/n`.
3. **`conf.cont_jcuts` is honoured:** J = 10 on one variable produces the expected number of cut columns for that variable and J = 4 produces fewer; the cut values equal the corresponding population quantiles of `df_super`.
4. **Determinism:** two calls with the same `(dgm, args, n)` give `identical()` families; two `fs_oc_predict()` calls at the same seed give identical output; different seeds differ.
5. **Override semantics:** explicit `c1`/`c2` override `forestsearch_args`, and `det_rate` is monotone non-increasing in `c1` across a short ladder.
6. **Size guard:** `max_M` below the realized M raises the guard's error.
7. **Gate parity, only if §2 opened:** on a one-candidate family, the `"resample"` declaration probability matches a direct `pnorm()` evaluation of the threshold derived in §2.

Follow the v5 §9 principle: at least one test must be shown to **fail against the defect it guards**. For test 4, verify that removing the seed handling makes it fail; report that you checked.

---

## 7. Close-out

- `devtools::document()`; `devtools::test()` — report raw counts and **warning-count parity** against the §1 baseline.
- `R CMD check` — full run; report Status verbatim. Zero errors/warnings/notes is the target, per CRAN priority.
- Version → `0.2.3`; `NEWS.md` entry naming both new exports and stating that the family is a population enumeration under the `forestsearch()` cut specification, with `consistency_method` mirroring the package's.
- Commit on `feature/glm-extension`, staging **by explicit path only** (`R/fs_oc_family.R`, `R/fs_oc_predict.R`, `tests/testthat/test-fs-oc-predict.R`, `man/`, `NAMESPACE`, `DESCRIPTION`, `NEWS.md`, the fidelity script). Never `git add -A`. Re-run `git status --porcelain` after and confirm any pre-existing dirt is untouched.
- **No push.** Larry pushes via GitHub Desktop.
- `devtools::install()` — `R/` changed, so this is warranted; confirm `packageVersion("forestsearch")`.

Report: `dev/glm-continuous-sims/REPORT_oc_wrapper_build_2026-08-28.md` — provenance; §2's `sigma_D` finding with source quotes and which branch of its gate you took; realized M for the fixture with the enumeration counts at each stage; the §5 fidelity result; test and check output raw; commit hash. Close with a ten-line verdict.

---

## 8. Out of scope — do not do these

- No edits to any existing `R/` file, to the prediction document, or to any driver.
- No grid, no inversion, no null-DGM path, no binary/OR path.
- No second family constructor; no specified-M path.
- No simulations, no replicate runs, no document renders.
- No new exports beyond the two functions named above.
- If §2's gate closes, no approximate `sigma_D` — report and stop on that path.
