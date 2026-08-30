# HANDOFF — applications re-run under 0.3.1 — 2026-08-31 — v1

**Written:** 2026-08-31 by the chat that ran this workstream, for the next chat that touches `quarto/applications/`. **Package:** forestsearch 0.3.1, branch `feature/glm-extension`. **Read this before the brief** (`brief_applications_rerun_20260831.md`): several of the brief's premises were corrected below.

---

## 1. What happened (commit shas are post-rewrite; see §4)

Four CC sessions, in order, all on pop-os, all reports committed under `dev/tasks/`:

1. **Read-only inventory** (`dfc8f3c9` task, `2767197c` report → `REPORT_applications_inventory_2026-08-31.md`). Census: 68 tracked `.qmd`, 19 active. Key finds: the 08-27 check rendered **0.2.2** against tracked **0.2.0 macOS** HTML baselines at `43b051b6` (HTML diff, ≈19 min, two docs only); no gbsg payloads existed on disk; one versionless LOO cache (`quarto/GuoHe/cv_out_analysis_gbsg_survival_effMaxSG_hrMaxSG_neighborhood.rds`, N=686); `HANDOFF_continuous_2026-08-27_v5.md` does not exist anywhere in history.
2. **Render-and-compare** (`728ccf5a` task, `d1f0e88e` run_loo flip, `8ad2c3f6` dated payloads, `2c2406a9` reference renders, `0e9d48fa` report → `REPORT_applications_render_0.3.1_2026-08-31.md`). Four docs rendered under 0.3.1 in ≈21m36s: template (0:12), compare_all (6:13), frozen_family (3:10, LOO off), multimethod (12:02). Triage vs 0.2.0 baselines: **0 SUBSTANTIVE hunks**; 12 CLAUSE-ORDER swaps; 4 cross-BLAS print deltas adjudicated VOLATILE (Larry accepted; pushed). Payloads staged to `_payloads_2026-08-31/`, all read back 0.3.1; template's tracked 0.2.0 payload restored beside them. GuoHe byte-identical; no LOO/CV executed; MR ran (5000 draws).
3. **Push size fix** (`4090d1d4` task, `13d5987e` report → `REPORT_push_size_fix_2026-08-31.md`). GitHub rejected the push: `comparison_continuous.rds` = **105.3 MB** (limit 100). History rewritten locally (nothing had reached origin): the payload commit now excludes it; file kept on disk, untracked, byte-identical (blob `15156d3e…`). **Sha map:** `713cd93e→8ad2c3f6`, `e872763d→2c2406a9`, `4ca38eec→0e9d48fa`. Backup branch `backup/pre-size-fix-2026-08-31` (delete after push confirmed; never push it).
4. **Roxygen use_grf correction** (`26cfce21` task, `ec116ffe` edit, `03127472` report → `REPORT_roxygen_use_grf_2026-08-31.md`). Docs-only (every changed `R/` line is `#'`); `man/forestsearch.Rd` + `man/interpret_search_config.Rd` regenerated; NEWS gained a `# forestsearch (development version)` heading + bullet. Package reinstalled on pop-os afterward so installed help matches.

**Branch tip at close: `03127472`.** Larry pushes via GitHub Desktop only.

## 2. Findings of record

- **0.2.0 → 0.3.1 changed no analysis result** in the four documents: identical subgroups (subject sets and clause sets), identical estimates/CIs/counts at printed precision. Headlines: template `{cd40≤421}ᶜ∧{wtkg≤82}ᶜ` N=73; compare_all same selections in all 8 combos; frozen_family harm n=75 HR 2.22 (1.17,4.22), complement 611 HR 0.61 (0.47,0.79), GH family 116; multimethod FS `{er≤0}∧{size≤35}`, DINA `{grade3≥1}∧{pgr≤10}` n=89 (13.0%, τ̂=0.0563), GRF `{er≤0}`.
- **The 08-30 (B) gate is closed empirically.** `vi.grf.min = NULL`'s entire measurable footprint on the applications is label clause order — never membership, never a number. Three-tier pattern, no residue: unpinned consistency-path calls swap (template 6, frozen_family 5); pinned calls don't (multimethod's consistency fits, compare_all's main call at `-0.2`); `subgroup_method = "dina"/"grf"` bypasses the screening block entirely (byte-identical selections).
- **Mechanism (from `R/` source, 29 Aug snapshot; corrected the inventory's reasoning):** the VI block (Section 5, `forestsearch_main.R` ~L2791) is gated **solely** by `!is.null(vi.grf.min)` — independent of `use_grf`, which gates GRF **candidate-cut generation** (Section 3A, ~L2491). At `-0.2` the block *orders* confounders by importance and filters nothing (importance ratios ≥ 0); `max_n_confounders` caps the ordered list. `use_grf = FALSE` does **not** moot the default.
- **Four cross-BLAS print deltas** (macOS/vecLib baseline vs Linux/reference BLAS), each byte-identical to an 08-27-attributed artifact: 13th-digit `33.6147229870207→208`, `-72.0371931311212→213`; band prose `79.12→79.13`; leaf table `3.85→3.86`. Adjudicated VOLATILE; ~1e-15 relative, far below all SEs.

## 3. Corrections superseding parts of committed documents (reports are not retro-edited; this section is the bridge)

- Inventory report §9.3 and render report §5.x claim `use_grf = FALSE` moots the `vi.grf.min` default ("no other unpinned call"). **Superseded** by §2's mechanism: the observed template/frozen_family swaps *are* the NULL path. The dina/grf **bypass** part of §5.x stands.
- Render report §8 lists **pre-rewrite** shas (`713cd93e`, `e872763d`, `4ca38eec`); the size-fix report §3 holds the old→new map.
- The 29 Aug snapshot's roxygen for `use_grf` ("Use GRF for variable importance") was the misstatement, fixed at `ec116ffe`; `@param vi.grf.min` and `@param max_n_confounders` were already correct and untouched.

## 4. Current state and standing values

- `run_loo <- FALSE` in frozen_family is a **standing committed change** (`d1f0e88e`); multimethod already had `run_cv/run_loo FALSE`, `run_mr TRUE` (5000 draws).
- **Payload tracking policy: files ≤ 50 MB are committed; larger are excluded** and manifested. Currently excluded: `quarto/applications/actg175/_payloads_2026-08-31/analysis_actg175_continuous_compare_all/comparison_continuous.rds` (105.3 MB, pop-os disk only).
- Untracked-by-design in the worktree: that file plus `dev/tasks/_baseline_html_2026-08-31/` and `dev/tasks/_diffs_2026-08-31/` (deletable at leisure; report quotes the essential hunks; baselines recoverable at `51fa758d`/`43b051b6`/`9b8d92ae`).
- The four `.html` at HEAD are now the **0.3.1 reference renders**; 0.2.0 vintages live at their baseline commits.
- Payload internals: definitions are **strings** (`labels$sg_harm` sgfocus-shape; `labels$fs/dina/grf` multimethod-shape), no membership vectors — any payload comparison must evaluate memberships on data, and clause order must be canonicalized first.

## 5. Pending (owner: Larry unless noted)

1. Push the last three commits (`26cfce21…03127472`) if not already pushed.
2. **Mac capture**: copy `_payloads` → `_payloads_mac_2026-08-31`; screen with `find quarto/applications -name '*.rds' -path '*_payloads_mac*' -size +50M`; leave hits unstaged in GitHub Desktop; commit + push. Uncheck-in-Changes prevents; it cannot undo a commit.
3. **Payload field-by-field comparison** (owner: next chat) — promised once the Mac vintages land: 0.2.x imported vs `_payloads_2026-08-31`, membership-based, table/meta/extras field-level. This is the precision layer beneath the HTML text diff.
4. Mac-side `devtools::install()` on the next Mac session (installed help there is stale).
5. Delete `backup/pre-size-fix-2026-08-31` after push confirmed.

## 6. Open decisions (parked; none blocks anything)

- **Target (c)** — new canonical results under recommended specs for the manuscript. Note: today's renders are regression evidence, **not** manuscript payloads; the manuscript's ACTG175 section consumes the **binary** multimethod document (`analysis_actg175_binary_multimethod_fixed_family.qmd`), out of scope today. (c) needs its own design conversation.
- The three gbsg siblings (`effMaxSG`, `sgfocus`, `maxeff_mrconfirm`) un-rendered under 0.3.1; `effMaxSG` would **silently load** the stale versionless 0.2.x LOO cache — clear/version-key it before any such render. LOO/K-fold under 0.3.1 entirely unmeasured, by decision.
- Slim what compare_all saves (the excluded object carries 8 full fits, plots, console text) and/or a narrow gitignore for known-huge names, so future runs fit the 50 MB policy natively.
- Crossanalysis doc disposition (binary inputs have no active writer; half-inert) — in scope vs archive, undecided since the inventory.
- multimethod's stale `# default TRUE` comment on `use_grf` (echoed into HTML — fix in the next qmd-touching pass, not standalone).
- `interpret_search_config` message strings tie "variable importance" prose to `use_grf` — ambiguous (may describe VI inside the cut-generation forest); untouched.
- Optional `FS_RESULTS_DIR` env plumbing (pattern: frozen_family's `GH_DIR`) for render-native dated payload dirs.

## 7. Protocol (unchanged)

Larry decides — chat writes specs/task docs — CC executes all repo operations; task docs via `~/Downloads` → `dev/tasks/`, committed first; kickoffs short, inline, final content of the message; gates stop on failure; compute is a go/no-go with wall-clock stated; verify from source; Larry pushes via GitHub Desktop; CC never pushes; one CC task per fresh session; committed work is not re-verified; the binary/OR OC-wrapper port runs in its own parallel chat — do not touch it.
