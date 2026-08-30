# REPORT — applications render-and-compare under 0.3.1 — 2026-08-31

All gates passed (G1–G7); triage **PASSED** — 0 SUBSTANTIVE hunks across the four documents. Total render wall-clock ≈ 21m36s, under the ceiling. Payloads staged and committed; renders committed; no push.

## 0 Environment and version gate

| item | value |
|---|---|
| hostname | `pop-os` (128 physical cores — worker echoes below confirm this is the same host as the 08-27 re-render) |
| branch | `feature/glm-extension` |
| `git log -1 --oneline` (pre-task HEAD) | `2767197c applications: read-only inventory ahead of 0.3.1 re-run (report)`; task doc committed on arrival as `728ccf5a` |
| `git status --porcelain \| wc -l` | `0` at G1 |
| installed forestsearch | `0.3.1` (`Rscript -e 'packageVersion'`) |
| `DESCRIPTION` Version at HEAD | `Version: 0.3.1` (line 3) — **G2 pass, no install needed** |
| R | `R version 4.6.1 (2026-06-24)` |
| quarto | `/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto` (present; used for all renders), `--version` = `1.9.38` |

## 1 Flag audit and the run_loo edit

Verbatim grep hits before the edit (`grep -n -E 'run_(loo|cv|cv_kfold|mr)\s*<-'`; extended to `run_gh|run_pareto|run_mr_confirm`, no hits for those as assignments in the four docs):

- `template_actg175_continuous.qmd` — no hits (as expected).
- `analysis_actg175_continuous_compare_all.qmd` — no hits (as expected).
- `analysis_gbsg_survival_frozen_family.qmd`:
  - `90:run_loo    <- TRUE      # leave-one-out (N-fold) honest out-of-sample subgroup estimate`
  - `98:run_cv_kfold <- FALSE    # repeated 10-fold cross-validation (stability diagnostics)`
- `analysis_gbsg_survival_multimethod.qmd`:
  - `144:run_cv  <- FALSE    # K-Fold cross-validation (kfold-cv + diagnostics)`
  - `145:run_loo <- FALSE    # Leave-one-out / N-fold OOB (oob-cv)`
  - `180:run_mr <- TRUE     # add mr_inference = TRUE to the forestsearch call` (exempt; stays TRUE; `mr_draws <- 5000L` at L181)

Edit applied — exactly one line, alignment and comment preserved (`git diff --stat`: `1 file changed, 1 insertion(+), 1 deletion(-)`):

```
-run_loo    <- TRUE      # leave-one-out (N-fold) honest out-of-sample subgroup estimate
+run_loo    <- FALSE     # leave-one-out (N-fold) honest out-of-sample subgroup estimate
```

Commit: **`d1f0e88e`** `gbsg frozen_family: run_loo -> FALSE (core analyses only; LOO/CV off per decision 2026-08-31)`.

**G3 after the edit** — `grep -n -E '^run_(loo|cv|cv_kfold)\s*<-'` across the four:
```
analysis_gbsg_survival_frozen_family.qmd:90:run_loo    <- FALSE     # leave-one-out (N-fold) honest out-of-sample subgroup estimate
analysis_gbsg_survival_frozen_family.qmd:98:run_cv_kfold <- FALSE    # repeated 10-fold cross-validation (stability diagnostics)
analysis_gbsg_survival_multimethod.qmd:144:run_cv  <- FALSE    # K-Fold cross-validation (kfold-cv + diagnostics)
analysis_gbsg_survival_multimethod.qmd:145:run_loo <- FALSE    # Leave-one-out / N-fold OOB (oob-cv)
```
Every CV-type flag FALSE; `run_mr` TRUE (exempt). **G3 pass.**

## 2 Baselines

All four `.html` siblings are tracked; snapshots taken from `git show HEAD:<path>` into `dev/tasks/_baseline_html_2026-08-31/` (uncommitted scratch). Sentinel: `touch /tmp/fs_render_sentinel_20260831` (2026-08-30 12:05, after snapshots, before renders).

| doc | tracked? | baseline commit | in-file stamp |
|---|---|---|---|
| template_actg175_continuous | yes | `51fa758d 2026-08-17 actg175 continuous template ... need to revisit` | no legible stamp |
| analysis_actg175_continuous_compare_all | yes | `43b051b6 2026-08-16 update gbsg and actg175 analysis runs` | no legible stamp |
| analysis_gbsg_survival_frozen_family | yes | `9b8d92ae 2026-08-13 gbsg analyses update` | `forestsearch 0.2.0` |
| analysis_gbsg_survival_multimethod | yes | `43b051b6 2026-08-16 update gbsg and actg175 analysis runs` | `forestsearch 0.2.0`, `forestsearch_0.2.0` |

Note (matters for §5.3): the frozen_family **baseline HTML predates its `.qmd`'s later commits** — the payload-export section (`cf4d6432`, 2026-08-15), the doc/focus/rule LOO cache key (`48c419eb`, 2026-08-15), and today's `run_loo` flip (`d1f0e88e`) — so its diff necessarily contains source-evolution hunks on top of 0.2.0→0.3.1 behavior. The other three baselines are source-identical with their `.qmd` at render time (multimethod & compare_all per the 08-27 repro report §1; template rendered from `61dd99df`'s qmd, unchanged since).

## 3 Render timings

`/usr/bin/time -f "elapsed %E maxRSS %MkB" $QUARTO render <qmd> --to html`, repo root, scope order:

| # | doc | elapsed | maxRSS | exit | fresh outputs (all `-newer` sentinel) |
|---|---|---|---|---|---|
| 1 | template_actg175_continuous | **0:11.74** | 746 MB | 0 | `.html`; `_payloads/template_actg175_continuous/template_actg175_continuous_payload.rds` |
| 2 | analysis_actg175_continuous_compare_all | **6:13.48** | 1862 MB | 0 | `.html`; `_payloads/analysis_actg175_continuous_compare_all/{selected_subgroups_continuous.rds,comparison_continuous.rds}` |
| 3 | analysis_gbsg_survival_frozen_family | **3:09.55** | 726 MB | 0 | `.html`; `gbsg/_payloads/analysis_gbsg_survival_frozen_family/analysis_gbsg_survival_frozen_family_payload.rds` |
| 4 | analysis_gbsg_survival_multimethod | **12:01.57** | 890 MB | 0 | `.html`; `gbsg/_payloads/analysis_gbsg_survival_multimethod/analysis_gbsg_survival_multimethod_payload.rds` |

**Total: ≈ 21m36s** (0:11.74 + 6:13.48 + 3:09.55 + 12:01.57) — within the ≈35–55 min estimate's floor and well under the ~60 min ceiling. The 08-27 per-doc measurements transfer cleanly: multimethod 12m02 (vs 12m53), compare_all 6m13 (vs 6m00).

**G4 pass** for docs 1–3 (exit 0, `.html` and payload outputs newer than sentinel, verified per render). Template's 11.7s (vs "minutes" expected) is a genuine fresh render: its payload read back `0.3.1 | 2026-08-30 12:06:05 | md` immediately after, and no freeze/cache is configured (`grep -c freeze` = 0 in the qmd and `quarto/_quarto.yml`). frozen_family at 3m09s sits below the 10–30 min guess because LOO (previously ~686 refits) is now off.

## 4 Invariants — G5 pass

After all four renders: `git status --porcelain quarto/GuoHe` → empty (byte-identical, still only the one tracked cache); `find . -name 'cv_out_*.rds' -newer /tmp/fs_render_sentinel_20260831 -not -path './.git/*'` → nothing. No LOO/CV executed anywhere. (Also checked interim after render 3.)

Worktree change-set at that point (`git status --porcelain`), exactly the expected files and nothing else:
```
 M quarto/applications/actg175/_payloads/template_actg175_continuous/template_actg175_continuous_payload.rds
 M quarto/applications/actg175/analysis_actg175_continuous_compare_all.html
 M quarto/applications/actg175/template_actg175_continuous.html
 M quarto/applications/gbsg/analysis_gbsg_survival_frozen_family.html
 M quarto/applications/gbsg/analysis_gbsg_survival_multimethod.html
?? dev/tasks/REPORT_applications_render_0.3.1_2026-08-31.md
?? dev/tasks/_baseline_html_2026-08-31/
?? dev/tasks/_diffs_2026-08-31/
?? quarto/applications/actg175/_payloads/analysis_actg175_continuous_compare_all/
?? quarto/applications/gbsg/_payloads/
```

## 5 Diff triage

Method per task: `sed -E 's/data:image[^"]*//g'` — **widened to `s/data:[a-zA-Z0-9/+-]+[^"]*//g'`** because the pages embed fonts/CSS as `data:text/css`/`data:font` URIs that the narrower pattern left in (first template diff attempt was dominated by base64 font noise); the widened strip removes only data-URI payloads, no rendered text. Diffs live in `dev/tasks/_diffs_2026-08-31/<name>.diff` (uncommitted scratch).

### 5.1 template_actg175_continuous — 54 hunks · SUBSTANTIVE 0 · CLAUSE-ORDER 6 · VOLATILE 48

- **VOLATILE (48)** — timing echoes (`Subgroup search completed in 0.01→0.03 minutes`; `Seconds and minutes forestsearch overall = 1.641 0.0274 → 2.12 0.0353`); the 0.3.x message-stream change replacing the batch echo block (`Parallel config: workers = 1 , batch_size = 2` + 12 `Batch i / 12` lines) with `Parallel config: sequential plan -- running the consistency screen as a plain loop` + 23 per-candidate `Subgroup k: resample Pcons=…` echoes — the three PASSED candidates (Pcons 0.8, 0.93, 0.92) match the baseline's candidate-table values (0.800, 0.930, 0.920) exactly; gt CSS id churn (`yzulaxshcj→lppjflwfsy`, `gyvjgwfihi→xbsappnmlh`, ~44 hunks of id-scoped CSS); payload path `/Users/…` → `/home/…`.
- **CLAUSE-ORDER (6)** — same clause set, order swapped; all adjacent numeric cells (N, E, d1, Pcons) identical:
```
-Subgroup identified: !{cd40 <= 421} & !{wtkg <= 82}
+Subgroup identified: !{wtkg <= 82} & !{cd40 <= 421}
```
```
-<td headers="M.1">!{cd40 <= 421}</td> / -<td headers="M.2">!{wtkg <= 82}</td>   (row: N=73, E=73, d1=36, Pcons=0.800)
+<td headers="M.1">!{wtkg <= 82}</td> / +<td headers="M.2">!{cd40 <= 421}</td>
-<td headers="M.1">!{cd40 <= 421}</td> / -<td headers="M.2">{age <= 34}</td>    (row: N=137, E=137, d1=69, Pcons=0.930)
+<td headers="M.1">{age <= 34}</td>   / +<td headers="M.2">!{cd40 <= 421}</td>
-<td headers="M.1">!{cd40 <= 421}</td> / -<td headers="M.2">{age <= 35}</td>    (row: N=151, E=151, d1=73)
+<td headers="M.1">{age <= 35}</td>   / +<td headers="M.2">!{cd40 <= 421}</td>
```
```
-…<strong>Identified subgroup:</strong> !{cd40 <= 421} & !{wtkg <= 82}   (×2 gt footnotes)
+…<strong>Identified subgroup:</strong> !{wtkg <= 82} & !{cd40 <= 421}
```
- **SUBSTANTIVE: none.** Image blocks: 0 / 0.

### 5.2 analysis_actg175_continuous_compare_all — 46 hunks · SUBSTANTIVE 0 · CLAUSE-ORDER 1 · VOLATILE 45

- **CLAUSE-ORDER (1):**
```
-<pre><code>definition : !{cd40 <= 507} & {age <= 37}</code></pre>
+<pre><code>definition : {age <= 37} & !{cd40 <= 507}</code></pre>
```
- **VOLATILE (45)**, examples of each family:
  - Timing: `maxeffCons MR anchor elapsed: 8.8 s → 61.9 s`; per-combo `Seconds and minutes forestsearch overall = 12.x → 18–19.x` (×8); timing-table cells `143.3/2.4 → 275.4/4.6`, `175.9/2.9 → 362.4/6.0`; `Comparison completed in 2.4 → 4.6 minutes`.
  - Worker/batch topology: `Parallel config: workers = 13 , batch_size = 26 → workers = 121 , batch_size = 141`; `Batch 1..6 / 6 → Batch 1 / 1` (×7 each) — 121 = 0.95·128, confirming same-host class as the 08-27 re-render (102/121 workers).
  - Message-stream relocation (repro-report item 2, baseline predates `16e6bd96`): `MR harm confirmation: MD = 33.615 (point rule vs 0.00) -> consistent with harm` present in baseline stdout, absent from the new render; the code-listing split around it heals, cascading `cb37→cb35`-style renumbering through every later listing (the bulk of the 1,887 structural residue lines).
  - Cross-BLAS last-digit dust, **matching repro-report item 4 verbatim** (macOS/vecLib baseline vs Linux/reference-BLAS):
```
- "33.6147229870207" "-72.0371931311212"  "139.266639105163" "-55.0511600557639"
+ "33.6147229870208" "-72.0371931311213"  "139.266639105163" "-55.0511600557639"
```
  - Band-threshold 2-dp print boundary, **matching repro-report item 5 verbatim** (both "Selection criterion" prose lines, hrMaxSG and hrMinSG): `mean difference >= 79.12 → >= 79.13`. Band membership unchanged — no comparison-table line differs.
  - Paths `/Users/…` → `/home/…`.
- **SUBSTANTIVE: none.** Classification note: the two numeric items above are printed values that changed, but each reproduces byte-for-byte the difference the 08-27 repro check attributed to BLAS arithmetic against this **same 0.2.0 macOS baseline** (inventory report §6.3, repro §4 items 4–5); they are classified VOLATILE as pre-attributed environment echoes, and flagged here so Larry can overrule. Image blocks: 4 / 4, not compared.

### 5.3 analysis_gbsg_survival_frozen_family — 27 hunks · SUBSTANTIVE 0 · CLAUSE-ORDER 5 · VOLATILE 12 · attributed-to-source/flag 10

Baseline (`9b8d92ae`, 0.2.0) predates three committed source changes (§2 note), so hunks fall in four groups:

- **CLAUSE-ORDER (5)** — same clause sets, all counts/HRs identical:
```
-   1   2.222     75    41    2  *   {er <= 0} & {pgr <= 32.5}      →  {pgr <= 32.5} & {er <= 0}
-   6   1.133    138    73    2  *   {pgr <= 32.5} & !{meno}        →  !{meno} & {pgr <= 32.5}
-   8   1.087     92    47    2  -   !{size <= 25} & {grade3}       →  {grade3} & !{size <= 25}
    (each appears twice: pre- and post-filter candidate tables — same values both times)
-*** Subgroup found: {er <= 0} {pgr <= 32.5}                        →  {pgr <= 32.5} {er <= 0}
-Subgroup identified: {er <= 0} & {pgr <= 32.5}                     →  {pgr <= 32.5} & {er <= 0}
-Selected subgroup: {er <= 0} & {pgr <= 32.5}   |  fit 0 min        →  {pgr <= 32.5} & {er <= 0}   |  fit 0 min
```
- **VOLATILE (12)** — workers/host (`Using 11 workers (80% of 14 physical cores) → Using 102 workers (80% of 128 physical cores)`, `System max cores available: 14 → 128`, bootstrap/plan echoes, timing-table cells `1.082/0.02 → 2.622/0.04`, `92.111/1.54 → 147.142/2.45`, `0.147/0.00 → 0.613/0.01`, `12.575/0.21 → 30.356/0.51`; `FB completed in 1.5 → 2.5 min`; `MR gate: 0.15 s → 0.61 s over 5,000 draws`; `GH family size: 116 | selected: sg_0001 | GH 0.2 → 0.5 min` — family size and selection identical); paths; provenance-table caption `(forestsearch 0.2.0, read at render time) → (forestsearch 0.3.1, …)` and its default column `vi.grf.min: -0.2 → NULL` ×2 — the table echoes the *package default* read at render time (this doc passes no `vi.grf.min` and runs `use_grf = FALSE`, so the value is inert; NEWS.md L52-54 documents the default flip).
- **Attributed to the committed `run_loo` flip (`d1f0e88e`)** — baseline's LOO output removed, skip-note added: `LOO loaded from cache: ../../../quarto/GuoHe/cv_out_frozen_family.rds`, `LOO [stable]: out-of-sample harm HR 2.22 (1.17, 4.22) (n = 75 (10.9%)) …`, the Cross-Validation Summary and Cross-Validation Metrics gt tables (with their ~106 id-scoped CSS lines / 108 bare-brace lines), and the summary-row cell `2.22 (1.17, 4.22) → — (run_loo = FALSE)`; added `Leave-one-out skipped (run_loo = FALSE). …` and `Leave-one-out (N-fold, 11 → 102 workers)` label. Note the baseline's cache line names the **old focus-less cache key** (`cv_out_frozen_family.rds`) — pre-`48c419eb`, consistent with §2.
- **Attributed to post-baseline source commits** (`cf4d6432`/`48c419eb`) — new final section `1 Reproducibility payload` (+ its `Wrote payload to: /home/…` output), TOC/section renumbering (`1…8 → 0.1…0.8` with the new section as `1`), and 1,674 changed code-listing `<span id="cb…">` lines re-echoing the evolved source.
- **SUBSTANTIVE: none** — every estimate, CI, count and membership size printed by both versions is identical (ITT 686, harm subgroup n=75/HR 2.22 (1.17, 4.22), complement 611/0.61 (0.47, 0.79), GH family 116/sg_0001). Image blocks: 1 / 1, not compared.

### 5.4 analysis_gbsg_survival_multimethod — 196 hunks · SUBSTANTIVE 0 · CLAUSE-ORDER 0 · VOLATILE 196

The 90-line text residue (after removing gt-id CSS churn and structural HTML) is fully accounted for; families with examples:
- Render date `August 15, 2026 → August 30, 2026`; paths `/Users/… → /home/…`.
- Workers/topology: `Using 11 of 14 total cores → Using 102 of 128 total cores`; `Parallel config: workers = 11 , batch_size = 22` + `Batch 1..6 / 6` → `workers = 102 , batch_size = 123` + `Batch 1 / 1` (and the 42-candidate analog for the DINA-screen fit).
- Timing echoes: `Subgroup search completed in 0.05 → 0.12 minutes`; `Seconds and minutes forestsearch overall = 4.713 0.0785 → 14.375 0.2396` (and `2.855 0.0476 → 6.766 0.1128`); `ForestSearch completed in 7.4 → 30 seconds`; `Bootstrap completed in 10.1 → 5 minutes`; `DINA-mode bootstrap 1.2 → 2.6 minutes`; `GRF-mode bootstrap 3.8 → 3.2 minutes`.
- Alias announcements added ×4 (repro-report item 1, from `16e6bd96`, which postdates the 0.2.0 baseline): `sg_focus 'effMaxSG' resolves to canonical rule 'hrMaxSG' (aliases: effMaxSG).`
- Package-attach noise added (`Attaching package: 'data.table'` / `%notin%` masking) and the full sessionInfo block (R 4.5.2/aarch64-apple-darwin20/vecLib → R 4.6.1/x86_64-linux/reference BLAS 3.12.0; `forestsearch_0.2.0 → forestsearch_0.3.1`).
- Cross-BLAS 2-dp print boundary, **matching repro-report item 6 verbatim** — GRF screening-illustration leaf table, leaf.node 5 `control.mean`:
```
-#> 21         5         3.85       177.00       1.87     2
+#> 21         5         3.86       177.00       1.87     2
```
  (all other leaves and every downstream value identical; classified VOLATILE as a pre-attributed environment echo, same rationale as §5.2's items).
- **SUBSTANTIVE: none. CLAUSE-ORDER: none** — the headline `Subgroup identified: {er <= 0} & {size <= 35}` and every candidate-table clause string appear in identical order. Image blocks: 8 / 8, not compared.

### 5.x vi.grf.min NULL path, measured

`fs_dina` and `fs_grf` are the only in-scope calls that exercise the new `vi.grf.min = NULL` default (inventory §9.3). Old (0.2.0 baseline) vs new (0.3.1 render) selected-subgroup lines — **byte-identical**, no hunk touches them:

| fit | baseline (0.2.0) | 0.3.1 render |
|---|---|---|
| `fs_dina` | `SELECTED: {grade3 >=  1} & {pgr <= 10}  (n = 89, mean tau-hat = 0.0563)` (L12128); `DINA-selected subgroup: {grade3 >= 1} & {pgr <= 10}`; `Subgroup size: 89 (13.0% of ITT)` | same, verbatim (new HTML L12152, L12155, L12179) |
| `fs_grf` | `[forestsearch] GRF selection (subgroup_method = "grf") … grf_depth: 2 / dmin.grf: 0 / frac.tau: 0.8 / n.min: 60 / SELECTED (depth 1): er <= 0`; `GRF-selected subgroup: {er <= 0}` (L14004–14028) | same, verbatim (new HTML L14031–14055) |

This is consistent with the in-code comments in both fits: `subgroup_method = "dina"/"grf"` bypasses the screening block where `vi.grf.min` acts, so the NULL default is inert here. The four documents contain no other unpinned `vi.grf.min` call (frozen_family/template run `use_grf = FALSE`; every other fit pins `-0.2`).

## 6 Triage verdict

**PASSED.** 0 SUBSTANTIVE hunks in any document. CLAUSE-ORDER: 12 total (template 6, compare_all 1, frozen_family 5, multimethod 0), every one a reordering of an unchanged clause set with identical adjacent statistics. Four numeric print deltas exist, each reproducing an 08-27-attributed cross-BLAS artifact byte-for-byte against the same macOS 0.2.0 baselines (§5.2: 13th-digit estimates, 79.12→79.13 band prose; §5.4: leaf 5 `3.85→3.86`); classified VOLATILE with the rationale stated inline — flagged for Larry to overrule if he reads them as substantive.

## 7 Payload staging

- Pre-move `git status --porcelain` showed the fresh outputs as `??` (`gbsg/_payloads/` wholly untracked, `actg175/_payloads/analysis_actg175_continuous_compare_all/` untracked) and the template payload as ` M` (tracked, overwritten — expected).
- Moves: the three untracked dirs → `quarto/applications/{gbsg,actg175}/_payloads_2026-08-31/`; the now-empty `gbsg/_payloads/` removed by `rmdir`.
- Template: `cp -r` into the dated dir, then `git checkout -- quarto/applications/actg175/_payloads/template_actg175_continuous/` — path clean afterwards; restored file reads back `0.2.0 | 2026-08-17 09:17:37` (both vintages survive).
- Readback table (**G7 pass** — every staged payload `forestsearch_version == 0.3.1`):

| staged payload | version | built_at | est_scale |
|---|---|---|---|
| `gbsg/_payloads_2026-08-31/analysis_gbsg_survival_multimethod/…_payload.rds` | 0.3.1 | 2026-08-30 12:28:01 | hr |
| `gbsg/_payloads_2026-08-31/analysis_gbsg_survival_frozen_family/…_payload.rds` | 0.3.1 | 2026-08-30 12:15:57 | hr |
| `actg175/_payloads_2026-08-31/template_actg175_continuous/…_payload.rds` | 0.3.1 | 2026-08-30 12:06:05 | md |

- compare_all's two `.rds` (no version field, as expected):

| file | mtime | top-level names |
|---|---|---|
| `selected_subgroups_continuous.rds` | 2026-08-30 12:12:12 | analysis,combo,combo_label,sg_focus,selection_rule,outcome_type,effect_measure,subgroup,N_H,N_Hc,effect,Pcons,K,n_passed,on_frontier,error |
| `comparison_continuous.rds` | 2026-08-30 12:12:31 | combos,fs,ci_tab,plots,plot_grid,plot_combined,plot_combined_subsets,combined_skip_reason,console,diagnostics,errors |

- `git check-ignore -v` on both dated dirs: no output, exit 1 — **not ignored**; plain `git add` used (no `-f` needed).

## 8 Commits

| sha | message | contents |
|---|---|---|
| `728ccf5a` | `tasks: add applications render task (2026-08-31)` | task document (on arrival) |
| `d1f0e88e` | `gbsg frozen_family: run_loo -> FALSE (core analyses only; LOO/CV off per decision 2026-08-31)` | the single flag edit |
| `713cd93e` | `applications: 0.3.1 payloads under _payloads_2026-08-31 (dated; committed for cross-machine provenance)` | the two dated payload dirs |
| `e872763d` | `applications: 0.3.1 reference renders (actg175 template+compare_all; gbsg frozen_family+multimethod)` | the four `.html` (replacing the 0.2.0 baselines at HEAD; prior vintages remain at `51fa758d`/`43b051b6`/`9b8d92ae`) |
| (this report's commit — sha recorded in `git log` below) | `applications: 0.3.1 render-and-compare report` | this file |

Scratch dirs left uncommitted, as specified: `dev/tasks/_baseline_html_2026-08-31/` (HEAD-pinned baseline snapshots), `dev/tasks/_diffs_2026-08-31/` (the four text-layer diffs). **No push** — each commit is separate and droppable pre-push.

## 9 Deviations and notes

1. **Diff strip widened** from the task's `data:image` to all `data:` URIs (§5 preamble) — the narrower pattern left megabytes of base64 font/CSS payloads in the template diff; no rendered text is affected.
2. **frozen_family's baseline is not source-identical** with today's `.qmd` (baseline HTML at `9b8d92ae` 2026-08-13 predates the payload-export section, the keyed LOO cache, and the ordered `run_loo` flip). Its diff therefore carries attributed source-evolution and flag-flip hunks (§5.3), reported in their own groups rather than forced into VOLATILE; none is a result change.
3. **Four numeric print deltas classified VOLATILE by attribution** (see §6) — each matches an item the 08-27 repro report attributed to macOS/vecLib vs Linux/reference-BLAS arithmetic (items 4, 5, 6), against the same baselines. If Larry prefers the strict reading (any changed printed estimate ⇒ SUBSTANTIVE), the renders and payloads are committed and droppable pre-push.
4. **Template rendered in 11.7 s** (task expected "minutes") — verified genuine: payload `built_at` matches the render minute, version 0.3.1, no quarto freeze/cache configured.
5. **frozen_family well under estimate** (3m10 vs 10–30 min guessed) — the estimate priced the GH bootstrap conservatively and the old LOO; with `run_loo = FALSE`, GH (B=1000) took 0.5 min.
6. **Clause-order instability under 0.3.1 is real but presentation-level**: 12 CLAUSE-ORDER hunks across three documents (§5.1–5.3), always the same clause set with identical statistics. Downstream string-matching consumers (e.g. any future membership comparison keyed on definition strings) should canonicalize clause order first — the inventory (§4) already recommends membership-based comparison over string equality.
7. The multimethod render also refreshed the header date and sessionInfo; the tracked 0.2.0 stamps now exist only in the git history of the four `.html` (baseline commits listed in §2) and in the uncommitted `_baseline_html_2026-08-31/` snapshots.
