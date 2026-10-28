# HANDOFF — the two CRAN vignettes rewritten; the 0.4.1 release (from the chats of 2026-10-27/28)

**For:** the next forestsearch chat: finish the vignette workstream on `feature/vignettes-0.4.1`, run Larry's check, decide the ACTG175 article, then the CRAN 0.4.1 submission. **Written:** 2026-10-28, ~12:45 PT, from the CC transcripts and committed reports of the two days. **Supersedes** `HANDOFF_vignettes_release_0.4.1_2026-10-27.md` (repo root, commit `174778d0`) for everything about the vignettes; that file's §2 (what the package does differently vs CRAN 0.4.0) and §5 (release checklist, carried forward and updated in §8 below) remain the record for the pre-vignette work. Re-verify from source before relying on any number or line reference; the project KB is stale for R code.

---

## 0. Read this first

- **Target: submit forestsearch 0.4.1 to CRAN on Monday 10 November 2026.** Version **0.4.1 confirmed** (Larry, 2026-10-27). Justification for a sub-two-month update: the user-visible `xlog` fix (GitHub #7) and `forestploter (>= 1.2.0)`; say so in `cran-comments.md`.
- **The release carries the two rewritten CRAN vignettes.** Both are written, committed and verified on `feature/vignettes-0.4.1`; the branch has not been checked (`devtools::check()`) or merged. That check is the next action (§7).
- **Protocol unchanged** (strict): Larry decides → chat writes specs, vignette text and task docs → CC executes every repo operation → Larry pushes via GitHub Desktop; CC never pushes/fetches/pulls. Task docs travel `~/Downloads` → `dev/tasks/`, committed first and alone; amended docs get new versioned names. Kickoffs short, inline, in a fenced code block, "clear CC first, then paste". In a pasted CC transcript, **everything after `>` or `❯` is CC, not Larry.** Gates are stop-on-failure, never stop-to-ask — and **a failed gate is a stop even when the reason is clear** (Larry's ruling 2026-10-28 after CC proceeded past a failed `diff -w` gate; that instance was accepted after the fact, not as a precedent). Never point CC at project-KB paths. Verify from source. Any change under `R/` is called out separately (none was made in this workstream). Do not re-run committed work.
- **Primary references** (verbatim): **2026a** León, L. F. and Anderson, K. M. (2026a). *Inference for standard trial estimands after data-driven subgroup discovery.* arXiv:2609.38361 — the two-term (MR) estimate, IJ interval, aligned selection, fixed-family scope, bootstrap comparison; GBSG application (Table S2, p. 46). **2026b** León and Anderson (2026b). *Interpretable subgroup identification across survival and generalized-linear outcomes.* arXiv:2610.28379 — the closed-form stability-selection screen, declaration calibration, **field / field-s one-sided bounds and the two-subgroup (Bonferroni) claim** (attribution corrected 2026-10-28: these are 2026b, not 2026a), selection rule as a modelling choice; GBSG Table 1 (p. 19), Tables S13/S15/S17; ACTG175 in S3. Both PDFs are in the project KB (`LeonAnderson_2026_InferenceStandardTrialEstimands.pdf`, `LeonAnderson_2026_InterpretableSubgroupIdentificationSurvival.pdf`, the second is paper + supplement). Manuscript repos (read-only for CC): `~/Documents/GitHub/fs-glms-interpretable` (2026b; `quarto/gbsg/analysis_gbsg_{mr,dina_mr,grf_mr}.qmd` are Table 1's sources, sibling HEAD `19c6ac41`), `~/Documents/GitHub/fs-post-selection` (2026a; Table S2 reads `manuscript/_payloads/analysis_gbsg_survival_multimethod_payload.rds`, built 2026-08-15 by forestsearch 0.2.0 from forestsearch's `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd` at `cf4d6432`).

---

## 1. State at handoff

| item | state |
|---|---|
| `master` | `174778d0` (handoff commit of 2026-10-27, on top of `30441802`); whether Larry pushed it is his to confirm; the pkgdown site is live from `30441802` |
| Branch | `feature/vignettes-0.4.1`, local only, from `174778d0`; **not pushed, not checked, not merged** |
| HEAD (last confirmed) | `25318e8c` "dev: report — vignette math pass (2026-10-28)". One more task was dispatched after it (§1.2) |
| Version | `0.4.0.9000`; CRAN carries 0.4.0 |
| Installed build (Pop!_OS) | `0.4.0.9000` built 2026-10-27 21:01:28 UTC from `936b6734`'s tree; **valid for every task so far because `git diff 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` is empty** (this replaced the NEWS-identity install gate, which fails since NEWS changed; Larry waived it once, 2026-10-28) |
| Rendered HTML on disk | `vignettes/forestsearch.html`, `vignettes/survival-post-selection.html` — git-ignored (`.gitignore:15 vignettes/*.html`); `vignettes/articles/forestsearch-2024.html` — **tracked, committed by Larry in `9567a568`** (4.4 MB); `vignettes/articles/methodology.html`, `paper_simulations.html` — untracked, not ignored (left alone by every task) |
| Project KB | **full** (a write of this handoff was refused for size); nothing filed in the two days — the task docs and reports are in the repo, this handoff travels as an attachment |

### 1.1 Commit chain on `feature/vignettes-0.4.1` (from CC transcripts and reports; `git log --oneline 174778d0..HEAD` is authoritative)

| task | commits | what |
|---|---|---|
| vignettes stage 0 (10-27) | `60d7ca33` task; `30aee979` record + report | Table 1 provenance, package facts, timing ladder (S0–S3, D1/D3, G1/G3; one core), transplanted calls, saved fits `fits/fit_{S3,D3,G3}.rds` |
| primary vignette (10-27) | `0908a3c1` task; `225cac4d` move; `c810dbd8` vignette; `0de9a8e5` `_pkgdown.yml` + NEWS + cross-ref; report commit | `vignettes/forestsearch.Rmd` → `vignettes/articles/forestsearch-2024.Rmd` (pointer note); new `vignettes/forestsearch.qmd` |
| Larry | `9567a568` | tracked `vignettes/articles/forestsearch-2024.html` |
| vignettes align (10-28) | `75bdac25` task; `946849ea` introduction v4; further commits through `c1fb38d3`; amendment `bf508916`, `44663ac8`; log commit | post-selection vignette aligned (p\* 0.90, 5,000 draws, `stop_threshold = NULL`, 2026a/2026b); `methodology.Rmd`, `paper_simulations.Rmd` pointer notes, deprecated args removed, `references.bib` entries; NEWS |
| rendering polish (10-28) | `7442213f` task; `400f7778`; `5a1af95d`; `99e8a64a` report; `d8756010` deviation note | `code-overflow: wrap`, `fontsize: 15px`, content-width tables, short declaration headers (both vignettes) |
| post-selection stage 0 (10-28) | `51fc55cc` task; `db5056f9` record + report | Table S2 provenance; naive + MR cells reproduced at 2 dp for FS/DINA/GRF; bootstrap per-replicate costs; fits `fits/fit_M-{FS,DINA,GRF}.rds` |
| post-selection vignette (10-28) | `c4ce38cb` task; `0384adc0` vignette v2; `437eeefb` attribution paragraph in `forestsearch.qmd`; `525e3813` NEWS; `e7f4fa6d` report | `survival-post-selection.qmd` rewritten (§2.2) |
| math + two-worker bootstrap (10-28) | `aceb413b` task v2; `2f36d0ce` `forestsearch.qmd` v6; `ce9cfb38` `survival-post-selection.qmd` v5; `25318e8c` report | LaTeX math in both vignettes; bootstrap on `n_workers`; bootstrap column as point estimate; plan-invariance recorded |

### 1.2 Dispatched, outcome not seen by this chat

`dev/tasks/cc_task_postsel_caption_2026-10-28.md` (kickoff given 12:35 PT): places `survival-post-selection_v6.qmd` (md5 `660b2f87bd26edd0a5e1fdd9801fd09e`, 499 lines) — the cost-table caption built from `n_workers` ("Wall-clock timings from this run: the multiplier stage under the fits' sequential plan, the bootstrap on 2 workers…"), the paragraph after it saying the bootstrap column is wall-clock on that many workers, the recovery sentence choosing its wording by `rc$q10` ("at least a tenth of the draws retain none of them" when the 10th percentile is exactly 0), "uses at most two, the limit a CRAN check imposes" in two places. Expected: three commits (task doc; `docs: post-selection vignette — cost caption and recovery sentence made exact for the two-worker bootstrap`; `dev: report — post-selection caption fix (2026-10-28)`), copies in `~/Downloads/postsel_caption_2026-10-28/`. **First thing in the new chat: read CC's transcript for it.** If it has not run, run it before anything else (clear CC first, then paste its kickoff: `Execute ~/Downloads/cc_task_postsel_caption_2026-10-28.md exactly as written. Input: ~/Downloads/survival-post-selection_v6.qmd (md5 660b2f87bd26edd0a5e1fdd9801fd09e). Gates stop on failure — a failed gate is a stop even when the reason is clear; report with a proposed fix. No push.`). Its gates: purl diff confined to the `cost` chunk; render; pandoc `--wrap=none -t plain` text checks (the caption phrase once, `Single-core` and `0% or less` absent, the six multiplier strings ≥ 2 each, bootstrap column `1.88 / 0.83 / 1.29`).

---

## 2. What is on the branch (from source; both vignettes are `.qmd`, `%\VignetteEngine{quarto::html}`, `embed-resources: true`, MathJax from Quarto's default CDN)

### 2.1 `vignettes/forestsearch.qmd` — "Interpretable subgroup identification with forestsearch" (final, Larry 2026-10-28; 662 lines at v6)

- **Runs Table 1 of 2026b end to end** on `survival::gbsg` (686 patients): the three identifiers by one `forestsearch()` call each, at the manuscript's settings, `plan = "sequential"`, ~74 s render (fits ≈ 38 + 10 + 23 s).
  - FS (`fit-fs`): `subgroup_method = "consistency"`, `conf.cont_jcuts = list(er = 10, pgr = 10)`, `conf_force = c("er <= 0", "pgr <= 0")`, `n.min 60`, `d0.min/d1.min 10`, `maxk 2`, `consistency_method = "resample"` (the closed form), `hr.threshold 1.00`, `hr.consistency 1.00`, `pconsistency.threshold 0.90`, `stop_threshold = NULL`, `sg_focus = "effMaxSG"`, `selection_rule = "neighborhood"`, `effect_neighborhood = 0.20`, `mr_inference = TRUE`, `mr_inference_args = list(ci_method = "field", draws = 5000L, include_complement = TRUE, field_scale_complement = "selected", field_decompose = TRUE, keep_declaration_field = TRUE, declaration_c0 = c(0.70, 0.75, 0.80, 0.85))`, `quiet = TRUE`.
  - DINA (`fit-dina`): `subgroup_method = "dina"`, `dina_select_statistic = "effect"`, same grid/thresholds. GRF (`fit-grf`): `subgroup_method = "grf"`, `grf_selection = "frontier"`, `grf_select_statistic = "effect"`, `grf_depth = 2L`, `dmin.grf = 0.0`.
- **Regulatory-framed introduction in bullets** (v4, approved): the goal is NOT to search a trial for a subgroup with a worse effect; the search as closed-form stability selection on splitting consistency; the three thresholds c₁ = `hr.threshold`, c₂ = `hr.consistency`, p\* = `pconsistency.threshold`; the selection rule as a modelling choice; selection-adjusted inference (2026a) and the field/field-s bounds, two-subgroup claim and declaration calibration (2026b).
- States that the consistency evaluation is **closed-form** (`consistency_method = "resample"`, `R/consistency_resample.R:455–457`: `max(0, 2*pnorm((beta_hat - log c2)/sigma_D) - 1)`), with random splitting (`fs.splits`, `"split"`) only as a rare fallback when the closed form is not evaluable; stage 0 traced 120/120 evaluated candidates on the closed form, 0 fallbacks.
- Sections: data; the FS fit and the identified region (`\x60` trick for inline backticks); the funnel (`fit$find.grps$filter_counts`: 2,211 / 1,410 / 121 / 9 / 1,744); the candidates passing the screen (`fit$grp.consistency$out_sg$result`, columns `Pcons hr N E g m K M.1 M.2`; definition joined from `M.1`/`M.2`); the two model-based identifiers; the three-way membership partition (= Table S15); **Table 1** built by `table1_column(fit, data)` (nine quantities: region, n/events Ĥ and Ĥᶜ, naive and corrected HRs, field lower Ĥ, field-s upper Ĥᶜ, Bonferroni pair); **declaration calibration** (`fs_declaration_calibration(fit_fs, alpha = 0.10, c0 = c0_grid)` and at c₂; table c₀ / κ̂ / implied p\* = 2Φ(κ̂) − 1 / fw / admits; `T_H <- cal$T_hat[[fld$family_id[fld$meta$selected_index]]]`); reading by location (never a significance dichotomy at HR 1); references (2024, 2026a, 2026b).
- **LaTeX math** (v6): display equations `p_{\mathrm{cons}}(g) = 2\,\Phi\!\left(\frac{\hat\beta(g) - \log c_2}{\sigma_D(g)}\right) - 1` and `T(g) = \frac{\hat\beta(g) - \log c_2}{\sigma_D(g)}`; 41 inline spans; prose-only change (purl byte-identical).
- YAML: `code-overflow: wrap`, `fontsize: 15px`; an `{=html}` style block making tables content-width (`#quarto-document-content table.table { width: auto; margin-left: 0; }`).
- Small text item still open: "The computation takes about a minute and a half on one core" (line ≈124) against a measured 1:14 render — "about a minute" would be exact.

### 2.2 `vignettes/survival-post-selection.qmd` — the GBSG analysis of 2026a with the full bootstrap illustrated (v5 committed `ce9cfb38`; v6 = caption fix, §1.2)

- Three fits at the Table S2 settings (ε = 0.10; `plan = "sequential"`): FS (`ci_method = "field"`, `draws = 5000L`, `field_recovery = TRUE`), DINA (`dina_args = list(seed = 8316951, n_folds = 2)`, `hr.threshold = 0.90`, `ci_method = "ij"`), GRF (frontier / effect / depth 2L / dmin 0.0 / frac.tau 0.8, `ci_method = "ij"`). Verified identical to stage 0's saved M fits on all 51 gated quantities; `field_recovery = TRUE` moves none of the four field bounds (consumes no draws).
- `mr_inference` in full (5,000 draws), then the **full bootstrap at a small B** for illustration: `B <- c(FS = 25L, DINA = 25L, GRF = 10L)`, `forestsearch_bootstrap_dofuture(fs.est, nb_boots, seed = 8316951L, mr_in_replicates = FALSE, parallel_args = boot_plan)` with `n_workers <- min(2L, future::availableCores())` and `boot_plan` multisession on `n_workers` (sequential if 1) — Larry: "use 2 cores" (2026-10-28). Text: the point estimates are an illustration, the intervals are not (finite-B variance needs B in the hundreds; the package warns below 100).
- Shows `summarize_bootstrap_results(sgharm = fit_fs$sg.harm, boot_results = boots$FS, create_plots = FALSE, est.scale = "hr")$table` (gt; naive + full bootstrap) and **`mr_estimates_table(fit_fs, boot_results = boots$FS, est.scale = "hr", digits = 2)`** (gt; naive + MR) — Larry's request; no package function gives the three-column Table S2 layout, so the vignette builds a side-by-side (naive / full bootstrap **point estimate only** at small B / multiplier) and prints Table S2 as published for comparison.
- Cost table (MR `timing_seconds` vs the bootstrap call; seconds per replicate; projection at B = 1,000, linear in B at that worker count); diagnostics: p̂ (`g$reselection$p_hat`) as a flag not an instrument; membership recovery (`g$field$recovery`: `sens_H`, `ppv_H`, q10/q50/q90, `share_equal_1`) with a live sentence; scope (fixed, enumerated family); references (Efron 2014, 2026a, 2026b, 2024). LaTeX math (v5).
- Render 2:07 on two workers (2:51 sequential): `boot_secs` FS 32.6 / DINA 23.6 / GRF 21.8 s (53 / 30.5 / 36 s sequential). Plan invariance (recorded, §2.5 of the math task): the FS bootstrap at B = 25 gives identical point estimates and an identical per-replicate table (32 of 34 columns; only the two timing columns and one timing attribute differ) under `sequential` and `multisession, 2` — the seed makes it deterministic, so the printed values are reproducible by a reader on any core count.

### 2.3 Articles and site

- `vignettes/articles/forestsearch-2024.Rmd` — the old vignette moved as-is (rename, 97 % similarity) with the title "(the 2024 workflow)" and a pointer blockquote (León et al. 2024, *Statistics in Medicine*, doi:10.1002/sim.10163); renders in 19 s.
- `vignettes/articles/methodology.Rmd`, `paper_simulations.Rmd` — pointer notes placing them as the 2024 workflow, citations to 2026a/2026b, deprecated arguments removed from listings, LASSO default comments fixed (amendment `bf508916`/`44663ac8`); `references.bib` has `Leon_Anderson_2026a/b` `@misc` entries.
- `_pkgdown.yml`: Articles navbar gains a `---` separator and "The 2024 workflow (León et al. 2024)" → `articles/forestsearch-2024.html`; `articles:` gains a `Reference` section; `home: description:` replaced (double-quoted — unquoted it is invalid YAML); `pkgdown::check_pkgdown()` "No problems found".
- `NEWS.md`: `### Documentation` under the development heading: the new getting-started vignette; the 2024 vignette as an article; `vignette("survival-post-selection")` reproduces the GBSG analysis of 2026a with the bootstrap at a small B and `mr_estimates_table()`; the methodology/simulation articles' notes.

---

## 3. Verified reproductions (gated; all PASS)

| what | result | record |
|---|---|---|
| Table 1 of 2026b, FS column | all nine quantities exact to full precision against the fs-glms-interpretable payload (region `{er <= 0} & {pgr <= 26}`, N 75, 41 events; naive 2.222 (1.169, 4.221); corrected Ĥ 1.276, Ĥᶜ 0.645; field lower 0.634; field-s upper 0.804; Bonferroni 0.509 / 0.840; n_family 1,744) | `dev/reports/vignettes_stage0_2026-10-27/comparison_table1.csv`; `REPORT_vignette_primary_2026-10-27.md` §2.4 (54/54 identical) |
| Table 1, DINA / GRF columns | exact (DINA `{grade3 >= 1} & {pgr <= 20}` N 113, corrected 1.076 / 0.658, n_family 30; GRF `{er <= 0}` N 82, 1.111 / 0.663, n_family 858) | same |
| Table S15 (three-way partition) | 35/19, 40/22, 78/42, 7/4, 526/212 | `REPORT_vignette_primary` §2.4 |
| Table S17 (declaration) | κ̂(c₀ 0.75) 2.726, implied 0.994, fw 0.667; T_H 2.438; T_max 2.564; κ̂(c₂) 3.503, fw 0.974; admits 0 everywhere | §2.5 |
| Table S2 of 2026a, naive + multiplier | all six rows equal as strings at 2 dp (FS `1.44 (0.62, 3.36)` / `0.64 (0.40, 1.03)`, DINA `1.11 (0.45, 2.75)` / `0.66 (0.40, 1.09)`, GRF `1.05 (0.47, 2.34)` / `0.66 (0.41, 1.06)`); at full precision DINA/GRF MR cells equal the 0.2.0 payload to machine precision, **FS MR cells differ in the 3rd decimal (1.4439 vs 1.4421) — uninvestigated** | `REPORT_postsel_stage0_2026-10-28.md`; `REPORT_postsel_vignette_2026-10-28.md` |
| Table S2 bootstrap column (not gated; B = 25/25/10 vs published B = 1,000) | FS 1.88 (1.94), DINA 0.83 (0.89), GRF 1.29 (1.04) | same |
| Determinism | S1 vs S1r identical; two renders and two plans give the same bootstrap values | stage-0 report §3; math report §2.5 |

Timings for `cran-comments.md`: primary vignette 1:14 on one core; post-selection 2:07 on two workers (2:51 on one); 2024 article 19 s (not built by check).

---

## 4. Rulings and decisions of 2026-10-27/28 (binding)

- Version 0.4.1 (not 0.5.0). Legacy vignette kept as an article (no `.Rmd` → `.qmd` conversion); primary vignette `.qmd`, all three identifiers, **runs the analyses at the manuscript's settings outright** — no precomputed results, no reduced draws, no `params` switch (stage 0 measured 75 s and exact reproduction).
- The vignette must say the consistency evaluation is closed-form, `fs.splits` a rare fallback.
- Introduction framing: the goal of `forestsearch()` is not to search a trial for a subgroup where the effect is worse than the trial's; comprehensive but in bullets.
- Post-selection vignette: reproduce the 2026a primary analysis (Table S2) analogously; a few bootstraps for illustration, full MR; exact reproduction of the bootstrap column not sought; illustrate `mr_estimates_table()`; bootstrap on two workers (`min(2, availableCores())`), not sequential.
- Mathematical typesetting: LaTeX via Quarto's default MathJax (CDN at viewing time; nothing downloaded at build; same arrangement as `rmarkdown::html_vignette`).
- Interpretation of bounds **by location** against clinically meaningful effect sizes; never a significance dichotomy at HR 1 (the 0.80 / 0.85 reading thresholds for the complement bound are clinical, not statistical).
- Attribution: two-term estimate and IJ interval → 2026a; field / field-s bounds, two-subgroup claim, declaration calibration → 2026b.
- Gate failures: stop even when the reason is clear (CC's "flaw in how the gate is written" call on the `diff -w` Table 1 block gate was accepted once, after the fact). The stage-0 ladder crash from CC's own harness (omitted `timings <- list()` setup chunk) was accepted as a harness bug, rerun from a clean record.
- Installed-build gate for documentation tasks: `git diff 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` empty (NEWS identity no longer usable). A reinstall before the release task resets this.
- Rendered HTML is never staged by a task; Larry tracked `forestsearch-2024.html` himself.

---

## 5. Task-doc lessons from this workstream (apply to the next task docs)

- **Gate on `pandoc --wrap=none -t plain` text, not raw HTML**: the raw HTML's MathJax polyfill contains `NaN`; prose wraps at 72 columns without `--wrap=none`, splitting phrases. Count occurrences with `grep -o … | wc -l`, not `grep -c` (lines). `NA (` matches "DINA (" — use `grep -cP '\bNA \('`.
- **Numbered-list hazard**: a prose line starting `2023)` renders as an ordered list; gate `grep -nE '^\s*[0-9]+[.)]( |$)'` on the source and `grep -c '<ol'` = 0 on the HTML.
- **Inline R cannot contain backticks**: `sprintf("\x60%s\x60", …)`.
- **Prose-only changes** are gated by `knitr::purl()` byte-identity (inline R is not purled); confined code changes by asserting every purl-diff line lies inside the named chunk markers.
- **Do not gate on `diff -w` of a rendered block that includes code** (rewrapped chunk code fails it); gate the table text or the identity values instead. Printed timing lines differ between renders — exclude them.
- Narrow literal greps to exact strings (`pconsistency.threshold = 0.80`, not `0.80`).
- Render logs: new file names per run; never overwrite the previous run's logs. Record `Rscript -e 'cat(future::availableCores())'` when a vignette's plan depends on it.
- Caps and measured times: primary render 74 s; post-selection 2:07 (two workers); 2024 article 19 s; `paper_simulations.Rmd` unmeasured here (cap 30 min).

---

## 6. Open findings for Larry's disposition (recorded in the reports; none queued, nothing blocks)

1. **`mr_estimates_table()` roxygen is stale**: says `"ij"` is the default; the FS path now defaults `ci_method` to `"field"` (`forestsearch_main.R:≈3854`), DINA/GRF keep `"ij"` by design (`fs_mr_inference_methods.R:137–139, 171`). Documentation-only fix under `R/` + `man/` regen — a candidate for the release task. The two-term estimate and IJ interval are identical under both methods; `"ij"` saves 15.9 s on the FS fit.
2. `ci_method` default differs by identifier path (above) — by design; consider one sentence in `?forestsearch`.
3. `field_R_out` / `field_R_in` are forwarded on the DINA/GRF MR path but not on the FS path (fixed 1000 × 500) — asymmetry; also why a `params` switch on `draws` saves ≤ 20 % on FS.
4. `print.forestsearch` shows `hrMaxSG` for `sg_focus = "effMaxSG"` (normalized name); cosmetic.
5. FS multiplier cells vs the 0.2.0 payload differ in the 3rd decimal (DINA/GRF exact) — uninvestigated.
6. MR stage timing: `timing_seconds` 4.6 s (FS) vs the manuscript's 0.38 s under 0.2.0 — the FS fit now runs the field block (R_out × R_in); the vignette says what each number measures.
7. DINA bootstrap per replicate 1.22 s vs the stage-0 projection 0.33 s (3×); FS and GRF within 10 %.
8. `remove_near_duplicate_subgroups()` keys on rounded summary statistics (`subgroup_consistency_helpers.R:1133`): drops exact-membership duplicates (the 121 → 120 candidate) and would collapse two distinct subgroups agreeing to 0.001 on ten statistics; called for every `sg_focus` but `"maxeff"` (`subgroup_consistency_main.R:618`).
9. The Table S2 document's DINA floor is `hr.threshold = 0.90`, and its DINA/GRF calls pass neither `conf_force` nor `conf.cont_jcuts` (DINA family 84 vs 30 in the Table 1 call) — check 2026a's text describes it so.
10. `vignettes/articles/forestsearch-2024.html` (4.4 MB) is tracked; `methodology.html` / `paper_simulations.html` untracked and not ignored — decide `.gitignore` for `vignettes/articles/*.html` before the release.
11. "about a minute and a half on one core" in `forestsearch.qmd` vs 1:14 measured (cosmetic).

---

## 7. Next steps, in order

1. **Confirm the caption task ran** (§1.2); if not, run it.
2. **Larry's `devtools::check()` on `feature/vignettes-0.4.1`** (`--as-cran`; the first R CMD check build of the quarto vignettes with MathJax and the two-worker bootstrap). Expect ≈ 3.5 min of vignette building where the old pair took 45 s, inside a check that took 12 min on `30441802`. Things to watch: `future::availableCores()` honours `_R_CHECK_LIMIT_CORES_` (→ 2) so `n_workers` resolves to 2 under check; no NOTE on vignette build time; confirm `.Rbuildignore` covers the root `HANDOFF_*.md` files and `dev/` (a "non-standard files at top level" NOTE would say otherwise).
3. **Decide task 3 vs release**: the ACTG175 GLM pkgdown article (binary OR and continuous MD, S3 of 2026b; `speff2trial` already in Suggests; a pkgdown article, not a CRAN vignette — Larry 2026-10-27) can follow the release; the release task (§8) is a half day. The window to 10 Nov allows both; the release does not depend on the article.
4. Merge `feature/vignettes-0.4.1` → `master` and push (Larry) — the pkgdown Action publishes the site, including the new vignettes and the 2024 article.

---

## 8. Release checklist for 0.4.1 (one CC task + Larry's checks; carried from the 10-27 handoff §5, updated)

1. Branch `release/0.4.1` from `master` (after the merge in §7.4). Task doc. Reinstall (`devtools::install()`) so the installed-build gate is simply "installed = HEAD".
2. `DESCRIPTION`: `Version: 0.4.1`. `NEWS.md`: `# forestsearch (development version)` → `# forestsearch 0.4.1`; re-read Deprecated / Performance / Internal / Bug fixes / Documentation for a CRAN reader.
3. Disposition of §6 item 1 (the stale `mr_estimates_table()` roxygen) — a one-line `R/` documentation change, called out separately.
4. `cran-comments.md`: changes since 0.4.0 (bug fix #7 + `forestploter >= 1.2.0`; deprecations `vi.grf.min`, `max_n_confounders`, `max.minutes` with warnings; performance; the two rewritten vignettes reproducing the manuscripts' GBSG analyses, build ≈ 3.5 min); why four weeks after 0.4.0; test environments; `0 errors | 0 warnings | 0 notes`; reverse dependencies none (forestploter is a forward dependency).
5. `devtools::document()`; `devtools::check(remote = TRUE, manual = TRUE)`; `check_win_devel()` / `check_win_release()`; optionally `rhub::rhub_check()`.
6. `devtools::spell_check()` (the vignettes add many new words: León, Bonferroni, multiplier, GBSG…), `urlchecker::url_check()` (arXiv DOIs, pkgdown URL).
7. Merge, push, `devtools::submit_cran()` — Larry only. After acceptance: tag `v0.4.1`, bump to `0.4.1.9000`, NEWS "(development version)" heading back. Deprecated formals removed in the release after.

---

## 9. Where everything lives

- **Task docs** (`dev/tasks/`): `cc_task_vignettes_stage0_2026-10-27.md`, `cc_task_vignette_primary_2026-10-27.md`, `cc_task_vignettes_align_2026-10-28.md`, `cc_task_vignette_render_2026-10-28.md`, `cc_task_postsel_stage0_2026-10-28.md`, `cc_task_postsel_vignette_2026-10-28.md`, `cc_task_vignettes_math_2026-10-28.md` (v1, superseded, never executed), `cc_task_vignettes_math_2026-10-28_v2.md`, `cc_task_postsel_caption_2026-10-28.md` (if run).
- **Reports and records** (`dev/reports/`): `REPORT_vignettes_stage0_2026-10-27.md` + `vignettes_stage0_2026-10-27/` (`package_facts.md`, `ladder.csv`, `ladder_quantities.csv`, `comparison_table1.csv`, `table1_calls/{fs,dina,grf}_call.R`, `fits/fit_{S3,D3,G3}.rds`); `REPORT_vignette_primary_2026-10-27.md` + `vignette_primary_2026-10-27/` (`gate_identity.R`, `check_pkgdown.txt`); `REPORT_vignettes_align_2026-10-28.md` + `vignettes_align_2026-10-28/`; `REPORT_vignette_render_2026-10-28.md` + record (`table1_block_diff_w.txt`); `REPORT_postsel_stage0_2026-10-28.md` + `postsel_stage0_2026-10-28/` (`package_facts_s2.md`, `table_s2_calls/`, `manuscript_source/`, `fits/fit_M-{FS,DINA,GRF}.rds`); `REPORT_postsel_vignette_2026-10-28.md` (with the math-pass and caption-fix sections appended) + `postsel_vignette_2026-10-28/`.
- **Vignette drafts** (chat's, in `~/Downloads`): `forestsearch{,_v2…_v6}.qmd` (v6 committed = `2f36d0ce`), `survival-post-selection_v{2…6}.qmd` (v5 committed = `ce9cfb38`; v6 the caption fix).
- **Close-out bundles** in `~/Downloads`: `vignettes_stage0_2026-10-27/`, `vignette_primary_2026-10-27/`, `vignette_render_2026-10-28/`, `postsel_stage0_2026-10-28/`, `postsel_vignette_2026-10-28/`, `vignettes_math_2026-10-28/`, `postsel_caption_2026-10-28/` (each with the rendered HTMLs and the report).
- **Machines:** Pop!_OS/System76 (128 cores), R 4.6.1, quarto CLI 1.10.18 at `/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto` (not on PATH; `RSTUDIO_PANDOC=…/quarto/bin/tools/x86_64`, pandoc 3.10), pkgdown 2.2.1, knitr 1.52, gt 1.4.0, future 1.80.0; CC from the repo root.
- **Prior handoff:** `HANDOFF_vignettes_release_0.4.1_2026-10-27.md` (repo root, `174778d0`) — §2 package changes since 0.4.0, §4 standing decisions (max.minutes mechanics, roxygen 8.1.1, never overwrite `quarto/applications/**`), §6 archive location.
