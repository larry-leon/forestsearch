# REPORT — reproduction check of the two application references against forestsearch 0.2.2

**Date:** 2026-08-27
**Verdict up front:** the currently installed package (0.2.2, built from `ae6f4025`'s
`R/`, HEAD `a4ec8c6d`) reproduces both pre-GLM-extension application references.
**No selection, estimate, CI, or bootstrap summary changed.** Every difference in
the rendered output is attributed below; none is unexplained.

---

## 1. The references

| application | files | .qmd last commit | .html last commit | rendered |
|---|---|---|---|---|
| GBSG survival multi-method | `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.{qmd,html}` | `cf4d6432` 2026-08-15 21:12 | `43b051b6` 2026-08-16 12:07 | Aug 15, 2026 |
| ACTG175 continuous compare-all | `quarto/applications/actg175/analysis_actg175_continuous_compare_all.{qmd,html}` | `cf4d6432` 2026-08-15 21:12 | `43b051b6` 2026-08-16 12:07 | Aug 16, 2026 |

Neither `.qmd` was modified after its `.html` was rendered, so the comparison is
source-identical. (The task that requested this check named the gbsg file
`analysis_gbsg_survival_method` and dated both references 17-Aug; neither
matches the tree. The files above are the pair committed together in
`43b051b6` "update gbsg and actg175 analysis runs", and they predate every
R/ commit in §3 — a stronger baseline than the one asked for.)

## 2. What kind of comparison this was

The references were rendered on a **different machine, R, and BLAS**:
R 4.5.2 on `aarch64-apple-darwin20` (macOS, Accelerate/vecLib), forestsearch
**0.2.0**. The re-render used R 4.6.1 on x86_64 Linux (reference BLAS/LAPACK
3.12), forestsearch **0.2.2**, with 102/121 workers where the references used
11/13. So this is a **cross-version, cross-platform, cross-BLAS**
reproduction, which is what makes the agreement below meaningful — and what
explains the only numeric deltas found (last-visible-digit floating dust).

Re-render mechanics: scratch-named `.qmd` copies in their own directories,
sequential renders with the installed package, exits 0/0, wall-clocks
12m53.5s (gbsg; reference ~15.3 min) and 6m00.2s (actg175). Payload exports
resolved to scratch-stem `_payloads/` directories, so no tracked artifact was
touched; all scratch files were deleted afterwards and the tree left clean.

## 3. R/ commits between the references and the installed build

| commit | date | existing R/ files edited? |
|---|---|---|
| `ae6f4025` promote `.sym_root()` to `fs_sym_root()` | Aug 27 | new file only (`fs_sym_root.R`) |
| `7f4174f2` add `fs_mr_oc_summary()` | Aug 26 | new file only (`fs_mr_oc_summary.R`) |
| `d00b5ecf` add `fs_dgm_scale()` | Aug 26 | new file only (`fs_dgm_scale.R`) |
| `83bd19de` builders own population noise | Aug 18 | **edits** `generate_glm_dgm.R`, `setup_gbsg_dgm.R` (default `k_random_noise = 0` output-identical; neither application builds a DGM) |
| `16e6bd96` sg_focus transparency + message stream | Aug 17 | **edits** `dina_subgroup.R`, `dina_subgroup_bootstrap.R`, `forestsearch_helpers.R`, `forestsearch_main.R` (+ new `fs_focus_tag.R`) |

Only `16e6bd96` could change rendered output, and only through messages — which
is exactly and exclusively what it did.

## 4. The six attributed differences

Tag-stripped whole-file diffs: 712 changed lines (gbsg), 122 (actg175).
640 lines were noise (per-render random `gt` CSS ids, render dates, wall-clock
timing echoes and timing-table cells, worker/batch topology echoes, package
attach messages, sessionInfo platform lines). Everything substantive:

1. **Alias announcement added** (both files; ×4 in gbsg):
   `sg_focus 'effMaxSG' resolves to canonical rule 'hrMaxSG' (aliases: effMaxSG).`
   → `16e6bd96`, the announcement that commit introduced. A message, not a result.

2. **`MR harm confirmation: MD = 33.615 (point rule vs 0.00) -> consistent with harm` line absent from the new actg render** (present in the reference).
   → `16e6bd96`, verified in its diff: the line moved from `cat()` (stdout,
   rendered) to `.mr_msg()`/message stream. The underlying value is unchanged —
   the full-precision estimate appears in both renders (item 4).

3. **Version strings** `forestsearch 0.2.0` → `0.2.2` ("read at render time"
   lines and sessionInfo). → the new-file commits plus version bump; expected.

4. **Full-precision printed estimates, actg175:** `33.6147229870207` →
   `33.6147229870208`, `-72.0371931311212` → `…213`. → environment, not a
   commit: 13th-significant-digit shifts from macOS/vecLib vs Linux/reference-
   BLAS arithmetic.

5. **Effect-band threshold in two "Selection criterion" prose lines, actg175:**
   `mean difference >= 79.12` → `>= 79.13`. → same class as item 4: the
   10%-of-frontier-maximum threshold inherits the shifted maximum and crosses a
   2-dp print boundary. Band membership and every selection are unchanged —
   no comparison-table line differs.

6. **GBSG GRF leaf table, leaf.node 5 `control.mean`:** `3.85` → `3.86`.
   → same class: grf's compiled arithmetic across platforms shifts the value
   across a 2-dp print boundary. All other leaves and everything downstream
   identical.

## 5. Incidental observations

- **Worker-count invariance, observed in the wild:** the references ran at 11
  (gbsg) and 13 (actg175) workers with 6 consistency batches; the re-renders
  ran at 102 and 121 workers with 1 batch. **No result moved.** This is the
  same invariance the dedicated test now pins
  (`tests/testthat/test-search-reproducibility.R`), here confirmed on the real
  applications across a 9–10x change in parallel topology.

- **Latent 2-dp boundary sensitivity, recorded not fixed:** item 5 shows the
  band-threshold prose (`>= 79.12` / `79.13`) rounds a derived quantity whose
  input differs across BLAS implementations at the 13th significant digit.
  Today the flip is display-only; a candidate sitting exactly on the band edge
  could in principle have its *membership* flip across platforms. Nothing in
  either application is near that edge, and no change is proposed here — this
  is a note for the display-precision / canonicalization queue.

## 6. Installed-build provenance

`packageVersion("forestsearch")` = 0.2.2, Built 2026-08-27 18:03 UTC — from
`ae6f4025`'s tree, the last commit touching `R/` (none since, through HEAD
`a4ec8c6d`). `deparse(body(forestsearch::fs_sym_root))` matches the `R/`
source, and all post-Aug-16 exports (`fs_dgm_scale`, `fs_scale_se`,
`fs_mr_oc_summary`, `fs_focus_tag`, `fs_sym_root`) are present. The installed
build is HEAD's `R/`.
