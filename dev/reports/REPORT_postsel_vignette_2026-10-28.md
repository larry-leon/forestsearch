# REPORT — post-selection vignette rewritten: the GBSG analysis of León & Anderson (2026a) with the full bootstrap illustrated (2026-10-28)

Task: `dev/tasks/cc_task_postsel_vignette_2026-10-28.md` (committed `c4ce38cb`). Record: `dev/reports/postsel_vignette_2026-10-28/`.
**Outcome:** every gate passed on its first evaluation. There were no deviations, and nothing under `R/`, `man/`, `DESCRIPTION`, `tests/` or `_pkgdown.yml` was touched. Not pushed.

## Provenance

- **Repo:** host pop-os, branch `feature/vignettes-0.4.1`.
  - HEAD at task start: `db5056f9`.
  - `git merge-base --is-ancestor db5056f9 HEAD`: pass.
- **Pre-existing untracked paths** (not staged): `vignettes/articles/methodology.html`, `vignettes/articles/paper_simulations.html`.
- **Input:** `~/Downloads/survival-post-selection_v2.qmd`, md5 `2f955e6926d231e4eae1bd191eeb2fa9`, 469 lines (GATE pass).
- **Install check (GATE 0.3 pass):**
  - Installed forestsearch is `0.4.0.9000`.
  - `git diff 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` is empty.
- **Environment:** R 4.6.1, survival 3.8.12, grf 2.6.1, future 1.80.0, doFuture 1.4.0, gt 1.4.0, knitr 1.52, quarto 1.10.18 (RStudio-bundled). `/proc/loadavg` read after the work: 0.91 / 0.95 / 0.64 on 128 threads.
- **Commits:**

| # | hash | message |
|---|---|---|
| 1 | `c4ce38cb` | dev: task — post-selection vignette (2026-10-28) |
| 2 | `0384adc0` | docs: post-selection vignette — the GBSG analysis of León & Anderson (2026a): three identifiers, multiplier resampling in full, the full bootstrap at a small B, mr_estimates_table() |
| 3 | `437eeefb` | docs: vignette — field bounds and the two-subgroup claim attributed to 2026b, the correction to 2026a |
| 4 | `525e3813` | docs: NEWS — post-selection vignette rewritten |
| 5 | (this report) | dev: report — post-selection vignette (2026-10-28) |

## §1 `vignettes/survival-post-selection.qmd`

- **Markdown gates (pass):** both greps are empty — the numbered-list grep and the `dev/|readRDS|read\.csv|load\(|\.\./` grep.
- **Render** (`render_postsel.log`; under `/usr/bin/time -v`, cap 20 min):
  - Exit 0.
  - Wall clock 2:50.61.
  - Peak RSS 977,668 kB.
  - HTML size 1,289,499 bytes.
- **HTML gates (pass):**

| pattern | count |
|---|---|
| `Error in` | 0 |
| `NaN` | 0 |
| `Warning` | 0 |
| `gt_table` | 8 (≥ 2 required) |

- **Selected regions, as rendered:**

| identifier | region | N | fit (s) | multiplier stage (s) | family |
|---|---|---|---|---|---|
| FS | `{er <= 0} & {size <= 35}` | 61 | 31.7 | 4.6 | 1744 |
| DINA | `{grade3 >= 1} & {pgr <= 10}` | 89 | 1.5 | 0.3 | 84 |
| GRF | `{er <= 0}` | 82 | 10.2 | 2.1 | 858 |

- **Other printed values:**
  - Whole trial: HR 0.695 (0.544, 0.888).
  - FS selection bias: 0.555 on the log scale.
  - FS p̂(Ĥ) = 0.237. The next four candidates have p̂ 0.140, 0.104, 0.058 and 0.052.

### Identity with stage 0 (GATE pass)

- **What ran:** `gate_identity.R` purled the vignette to `survival-post-selection_purled.R` and sourced it in a fresh `Rscript`. It then compared the results with `postsel_stage0_2026-10-28/fits/fit_M-{FS,DINA,GRF}.rds`.
- **Run cost:** wall 2:44.65 (cap 20 min); peak RSS 967,088 kB.
- **Result:** 0 gated mismatches out of 51 gated comparisons. Every gated quantity is `identical()` or `abs(diff) == 0`.
- **What that implies:** the arguments the vignette's calls drop are inert, and its explicit `ci_method` settings equal what the stage-0 calls resolved to: `field` for FS, `ij` for DINA and GRF.
- Full table: `identity.csv`.

| quantity | FS (vignette = stage 0) | DINA (vignette = stage 0) | GRF (vignette = stage 0) | equal |
|---|---|---|---|---|
| sg.harm | {er <= 0} & {size <= 35} | {grade3 >= 1} & {pgr <= 10} | {er <= 0} | TRUE ×3 |
| n_selected | 61 | 89 | 82 | TRUE ×3 |
| complement n | 625 | 597 | 604 | TRUE ×3 |
| naive est (lower, upper) | 2.53692 (1.24540, 5.16779) | 1.59936 (0.87961, 2.90806) | 1.95139 (1.03661, 3.67345) | TRUE ×9 |
| complement naive | 0.60794 (0.46860, 0.78871) | 0.62907 (0.48228, 0.82054) | 0.61500 (0.47232, 0.80078) | TRUE ×9 |
| debiased est (lower, upper) | 1.44395 (0.62061, 3.35958) | 1.11457 (0.45206, 2.74805) | 1.05385 (0.47459, 2.34013) | TRUE ×9 |
| complement debiased | 0.64410 (0.40160, 1.03304) | 0.65910 (0.39980, 1.08658) | 0.66089 (0.41025, 1.06464) | TRUE ×9 |
| n_family | 1744 | 84 | 858 | TRUE ×3 |
| selection_bias | 0.554570 | 0.357291 | 0.608769 | TRUE ×3 |
| ci_method (recorded) | field | ij | ij | TRUE ×3 |

- Values are shown at 5 dp here. The comparison itself was exact, at full precision.

**FS field quantities (recorded, not gated) — all `identical()` to M-FS** (`field_fs.csv`). Adding `field_recovery = TRUE` does not move the field bounds.

| quantity | vignette | stage 0 |
|---|---|---|
| `field$lower_1s` | 0.698976085274789 | 0.698976085274789 |
| `field$complement$upper_1s_s` | 0.805324673781695 | 0.805324673781695 |
| `field$joint_s$bonf_lower_H` | 0.572869542258108 | 0.572869542258108 |
| `field$joint_s$bonf_upper_Hc` | 0.835272174795449 | 0.835272174795449 |

**Membership recovery (FS, `recovery.csv`):**

| sens_H | ppv_H | q10 | q50 | q90 | share_equal_1 | n_used |
|---|---|---|---|---|---|---|
| 0.543 | 0.465 | 0 | 0.656 | 0.967 | 0.020 | 986 |

### Table S2 — naive and multiplier (GATE pass)

- **Check:** in all six rows, the `side_by_side` cells equal `table_s2`'s cells as strings, after the vignette's own 2-dp formatting. The `n (%)` cells agree as well. Full detail: `table_s2_check.csv`.

| method | SG | n | naive (vignette = S2) | multiplier (vignette = S2) | bootstrap, this run (recorded) | bootstrap, S2 at B = 1,000 |
|---|---|---|---|---|---|---|
| FS | H | 61 (9%) | 2.54 (1.25, 5.17) | 1.44 (0.62, 3.36) | 1.88 (0.01, 485.27) | 1.94 (0.89, 4.21) |
| FS | Hc | 625 (91%) | 0.61 (0.47, 0.79) | 0.64 (0.40, 1.03) | 0.63 (0.34, 1.20) | 0.63 (0.43, 0.93) |
| DINA | H | 89 (13%) | 1.60 (0.88, 2.91) | 1.11 (0.45, 2.75) | 0.83 (0.20, 3.45) | 0.89 (0.38, 2.07) |
| DINA | Hc | 597 (87%) | 0.63 (0.48, 0.82) | 0.66 (0.40, 1.09) | 0.66 (0.33, 1.35) | 0.66 (0.40, 1.11) |
| GRF | H | 82 (12%) | 1.95 (1.04, 3.67) | 1.05 (0.47, 2.34) | 1.29 (0.10, 15.93) | 1.04 (0.49, 2.24) |
| GRF | Hc | 604 (88%) | 0.61 (0.47, 0.80) | 0.66 (0.41, 1.06) | 0.73 (0.34, 1.57) | 0.66 (0.39, 1.13) |

- The rendered HTML carries the same bootstrap strings as the identity run: same seed, sequential plan.

### Bootstrap values and timings (`bootstrap.csv`, `cost.csv`, `fit_secs.csv`)

| identifier | B | H2 | H2 interval | Hc H2 | Hc interval | boot s (identity run) | boot s (render) |
|---|---|---|---|---|---|---|---|
| FS | 25 | 1.8786 | (0.0073, 485.27) | 0.6348 | (0.3368, 1.1967) | 53.04 | 54.3 |
| DINA | 25 | 0.8339 | (0.2014, 3.4535) | 0.6615 | (0.3251, 1.3462) | 30.50 | 30.6 |
| GRF | 10 | 1.2929 | (0.1049, 15.933) | 0.7347 | (0.3436, 1.5709) | 36.07 | 36.4 |

**Per-replicate cost.** Cost here is `bootstrap_seconds / B`, including the ~2.55 s of fixed per-call overhead measured in stage 0.

| identifier | this run (s/rep) | stage 0 total/5 (s/rep) | stage 0 mean iteration (s) | task projection (s) | this run (s) | MR stage (s) | projected at B = 1,000 (s) |
|---|---|---|---|---|---|---|---|
| FS | 2.12 | 2.19 | 1.673 | ~44 | 53.0 | 4.61 | 2122 |
| DINA | 1.22 | 0.85 | 0.332 | ~11 | 30.5 | 0.34 | 1220 |
| GRF | 3.61 | 3.99 | 3.482 | ~37 | 36.1 | 2.09 | 3606 |

- **Fit times in the identity run:** FS 31.66 s, DINA 1.34 s, GRF 10.36 s. Stage 0's were 31.83 / 1.73 / 11.28 s.

## §2 `vignettes/forestsearch.qmd`

- **Edit:** the attribution paragraph was replaced exactly as specified. The diff touches only those lines; it is a 6-line paragraph replaced by 7 lines.
- **Purl gate (pass):** `forestsearch_purl_before.R` and `forestsearch_purl_after.R` are byte-identical (`cmp`).
- **Render** (`render_forestsearch.log`): exit 0, wall 1:13.42, peak RSS 1,098,288 kB, HTML 1,285,738 bytes.
- **HTML gates (pass):**

| pattern | count |
|---|---|
| `1.276` | 4 |
| `1.076` | 2 |
| `1.111` | 2 |
| `<ol` | 0 |

## §3 `NEWS.md`

- The bullet was replaced exactly as specified. Committed alone.

## Deviations (verbatim)

None.

- One ordering note: the §2 and §3 file edits, and §2's purl check, were made while the §1 render was running. Each was committed only after the gates before it had passed, in the specified commit order. The §2 render ran after the §1 identity run, so no two timed runs overlapped.

## OPEN ITEMS

1. **The FS and GRF bootstrap intervals at small B are degenerate in the rendered table.**
   - FS Ĥ: 1.88 (0.01, 485.27). GRF Ĥ: 1.29 (0.10, 15.93).
   - The vignette says to read the interval columns as machinery output, not results, so this is within its framing. Still, a 485 upper bound printed in the side-by-side table may distract readers.
   - Option, for Larry's call: print only the bootstrap point estimate in `side_by_side` at this B. That would be a vignette change; not done here.
2. **The prose and the run disagree about how long the multiplier stage takes.**
   - The cost section says the manuscript reports "0.1 to 0.4 s for the multiplier stage". This run's `timing_seconds` are 4.61 / 0.34 / 2.09 s, and stage 0 recorded 4.68 / 0.49 / 2.41 s.
   - Likely cause: FS and GRF `timing_seconds` cover more than the manuscript's measure, e.g. the influence functions over 1,744 / 858 candidates, and for FS the field. This is unverified.
   - Proposed fix: make the sentence say what each number measures.
3. **DINA's bootstrap cost ran well over the task's projection:** 30.5 s against about 11 s, or 1.22 s per replicate against stage 0's 0.85. FS and GRF matched stage 0 within 10%. Possible explanations, neither checked:
   - Stage 0 had only 5 replicates, with iteration times from 0.10 to 0.58 s.
   - At B = 25 the regenerated DINA surfaces may include slower fits.
   - Recorded, not gated.
4. **The membership-recovery prose may need softening.** The prose reads "the typical re-selection retains much of the identified region".
   - What the run shows: q50 = 0.656 and sens_H = 0.543, but q10 = 0, so at least 10% of draws retain none of Ĥ. Only 2.0% of draws contain all of Ĥ.
   - The sentence holds at the median. Larry may want it to acknowledge the lower tail.
5. **The FS and DINA bootstrap point estimates at small B already sit near the published B = 1,000 values:** FS Ĥ 1.88 vs 1.94, DINA Ĥ 0.83 vs 0.89. GRF Ĥ (B = 10) does not: 1.29 vs 1.04. This is an observation only; exact reproduction was not sought.

## Copies to `~/Downloads/postsel_vignette_2026-10-28/`

(listing below)

```
total 2580
drwxrwxr-x  2 larryleon larryleon    4096 Oct 28 11:56 .
drwxr-xr-x 10 larryleon larryleon   20480 Oct 28 11:56 ..
-rw-rw-r--  1 larryleon larryleon 1285738 Oct 28 11:56 forestsearch.html
-rw-rw-r--  1 larryleon larryleon   10575 Oct 28 11:56 REPORT_postsel_vignette_2026-10-28.md
-rw-rw-r--  1 larryleon larryleon 1289499 Oct 28 11:56 survival-post-selection.html
-rw-r--r--  1 larryleon larryleon   22029 Oct 28 11:56 survival-post-selection.qmd
```

## Math, presentation and two-worker bootstrap pass (2026-10-28)

Task: `dev/tasks/cc_task_vignettes_math_2026-10-28_v2.md`. The task ran in three sessions. The first stopped at §0.3 because `forestsearch_v6.qmd` was not yet in `~/Downloads`. The second stopped twice at §1.4 and §2.4, on gate definitions; both were amended by Larry (deviations below).

### Commits

| # | hash | message |
|---|---|---|
| 1 | `aceb413b` | dev: task — vignette math v2 (2026-10-28) |
| 2 | `2f36d0ce` | docs: vignette — LaTeX math for the thresholds, the closed-form rate, the selection rule and the declaration statistic |
| 3 | `ce9cfb38` | docs: post-selection vignette — LaTeX math; bootstrap on two workers behind an n_workers knob; bootstrap column as point estimate at small B; cost and recovery sentences made exact |
| 4 | (this report) | dev: report — vignette math pass (2026-10-28) |

### §0

- **Branch and base:** branch `feature/vignettes-0.4.1`; HEAD at task start `e7f4fa6d`, which passes the ancestor check.
- **Untracked paths:** the pre-existing ones are `vignettes/articles/methodology.html` and `vignettes/articles/paper_simulations.html`. Neither was staged.
- **Input checksums:**
  - `forestsearch_v6.qmd`: `8ed0d81d373a4427a1a74a0b0dc9817f`, 662 lines.
  - `survival-post-selection_v5.qmd`: `8c56077e1fdbd0b1c13fbf4c5d082e15`, 494 lines.
- **Base checksums at HEAD:**
  - `vignettes/forestsearch.qmd`: `234322275f3df0ef4fee47c5509b81e5`.
  - `vignettes/survival-post-selection.qmd`: `2f955e6926d231e4eae1bd191eeb2fa9`.
- **Installed build:** `git diff 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` is empty.
- **Cores:** `future::availableCores()` returns 128, so `n_workers` resolves to 2 here.

### §1 `vignettes/forestsearch.qmd` ← v6

- **Prose only (pass):** `knitr::purl()` of the HEAD file and of v6 are byte-identical.
- **Markdown (pass):**
  - The numbered-list grep is empty.
  - β, Ĥ, Φ, κ̂ do not occur outside code chunks.
- **Render** (`render_forestsearch_math.log`; under `/usr/bin/time -v`):
  - Exit 0.
  - Wall clock 1:14.14.
  - Peak RSS 1,107,060 kB.
  - HTML size 1,373,170 bytes, against 1,285,738 last run.
  - The growth is the polyfill that Quarto embeds alongside the CDN MathJax reference (see deviation 1).
- **HTML gates as amended (pass; `gates_forestsearch_math.txt`):**

| check | count |
|---|---|
| pandoc plain `Error in` / `Warning` / `NaN` | 0 / 0 / 0 |
| pandoc plain `grep -cP '\bNA \('` | 0 |
| raw HTML outside `<script>`: `Error in` / `Warning` / `NaN` / `NA (` | 0 / 0 / 0 / 0 |
| `<ol` | 0 |
| `mathjax` (case-insensitive) | 4 |
| `class="math inline"` (`grep -o \| wc -l`) | 41 (≥ 20) |
| `class="math display"` (`grep -o \| wc -l`) | 2 (= 2) |
| `1.276` / `1.076` / `1.111` | 4 / 2 / 2 |

- **MathJax source:** the HTML references `https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml-full.js`. MathJax itself is not embedded.
- **Excerpts:** the pandoc plain text around the first display equation and around T(g) is in `excerpts_forestsearch_math.txt`. Pandoc shows the TeX source; the browser renders it. The first display equation reads:

```
p_{\mathrm{cons}}(g) \;=\; 2\,\Phi\!\left(\frac{\hat\beta(g) - \log
c_2}{\sigma_D(g)}\right) - 1, \] with \(\sigma_D(g)\) the scale of the
```

```
candidate's standardized statistic \[ T(g) \;=\; \frac{\hat\beta(g) -
\log c_2}{\sigma_D(g)}, \] with every candidate shifted to sit at
\(c_0\), under the law of the fit's own multiplier draws. The calibrated
rule admits a candidate when \(T(g) \ge \hat\kappa_\alpha(c_0)\); it
```

### §2 `vignettes/survival-post-selection.qmd` ← v5

- **Code change confined (pass; `gate_purl_postsel_math.txt`):** every line that differs between the purled scripts lies in the `bootstrap` chunk (5 lines) or the `three-columns` chunk (12 lines). The `fit-fs`, `fit-dina` and `fit-grf` chunks are identical, so the fit-identity result of `e7f4fa6d` stands.
  - The chunk of each differing line was found from `purl(documentation = 1)` chunk headers.
  - My first parse cut labels at the first hyphen; that was a bug in my own check script. It was fixed and the check rerun, with no change to the gate.
- **Markdown (pass):** both greps are empty.
- **Render** (`render_postsel_math.log`; under `/usr/bin/time -v`):
  - Exit 0.
  - Wall clock 2:07.47, against 2:50.61 sequential last run.
  - Peak RSS 992,776 kB.
  - HTML size 1,377,332 bytes.
- **HTML gates as amended (pass; `gates_postsel_math.txt`):**

| check | count |
|---|---|
| pandoc plain `Error in` / `Warning` / `NaN` | 0 / 0 / 0 |
| pandoc plain `grep -cP '\bNA \('` | 0 |
| raw HTML outside `<script>`: `Error in` / `Warning` / `NaN` | 0 / 0 / 0 |
| `gt_table` | 8 (≥ 2) |
| `mathjax` | 4 |
| `1.44 (0.62, 3.36)` / `0.64 (0.40, 1.03)` | 4 / 4 |
| `1.11 (0.45, 2.75)` / `0.66 (0.40, 1.09)` | 3 / 3 |
| `1.05 (0.47, 2.34)` / `0.66 (0.41, 1.06)` | 3 / 3 |

- **Bootstrap column:** "Full bootstrap estimate (B = 25 / 25 / 10)" holds bare numbers only (pass). In pandoc plain text it reads 1.88, 0.63, 0.83, 0.66, 1.29, 0.73.
- **Not a gate here:** the raw HTML outside `<script>` contains one `NA (`, inside "DINA (0.89 and 1.11)".
- **`boot_secs` on two workers:** FS 32.6, DINA 23.6, GRF 21.8 s. The previous run, sequential, was 53, 30.5 and 36 s.
- **Side-by-side H bootstrap estimates:** FS 1.88, DINA 0.83, GRF 1.29, identical to the previous run's.
- **Recovery sentence, as rendered:**

  > Read together: the median draw retains 66% of the region's members, a tenth of the draws retain 0% or less, and 2% of draws contain all of it. Exact re-selection is rare; what moves from draw to draw is the boundary. […]

- **Cost table, as rendered** (`excerpts_postsel_math.txt`):

| Method | MR (s) | B | bootstrap (s) | s per replicate | projected at B = 1,000 (s) |
|---|---|---|---|---|---|
| FS | 4.66 | 25 | 32.62 | 1.30 | 1305 |
| DINA | 0.34 | 25 | 23.65 | 0.95 | 946 |
| GRF | 2.12 | 10 | 21.80 | 2.18 | 2180 |

### §2.5 Plan invariance (recorded, not gated)

`plan_invariance.R` and its output `plan_invariance_output.txt`, run in a fresh `Rscript`.

- **Setup:**
  - The purled vignette is evaluated expression by expression until `fit_fs` exists. It asserts `fit_dina` does not.
  - Then `forestsearch_bootstrap_dofuture(fs.est = fit_fs, nb_boots = 25L, seed = 8316951L, mr_in_replicates = FALSE, ...)` runs once sequentially and once under `multisession` with 2 workers.
- **Wall clock:** sequential 53.2 s; two workers 32.3 s.
- **Estimates:** identical under both plans.
  - `H_estimates$H2` = 1.8786205217.
  - `Hc_estimates$H2` = 0.6348402476.
- **Per-replicate table (`results`, 25 × 34):** `identical()` is FALSE. The only differences are timing:
  - the columns `tmins_search` and `tmins_iteration`;
  - the `timing` attribute.

  With those two columns and the attributes dropped, the remaining 32 columns are `identical()` TRUE.
- **Conclusion:** the bootstrap is invariant to the plan apart from wall-clock fields. The vignette claims nothing about this.

### Deviations (verbatim approvals)

1. **§1.4 / §2.4 greps moved to rendered text; math spans counted by occurrence.**
   - Why: as written, `NaN` = 17 in `forestsearch.html`. All 17 matches are in a JavaScript polyfill (`isNaN(a)?NaN:…`) that Quarto embeds with MathJax under `embed-resources: true`. Also, `grep -c 'class="math inline"'` = 18 counts lines, while there are 41 occurrences.
   - Larry: "Approved: amend §1.4 and §2.4 so the Error in / Warning / NaN / NA ( greps run on `pandoc -t plain` of the HTML (record alongside that the raw HTML outside <script> also has 0), and count math spans with `grep -o ... | wc -l`; record the amendment as a deviation and resume from §1.4."
   - The §1 render was not repeated, because the qmd was unchanged since it was rendered. The amended gates were applied to that HTML.
2. **`NA (` anchored on a word boundary.**
   - Why: under deviation 1, §2.4 counted `NA (` = 1 in the plain text, from "DINA (0.89 and 1.11)".
   - Larry: "Approved with option 2: in §1.4 and §2.4 replace the `NA (` grep by `grep -cP '\bNA \('` on the pandoc plain text (0 required); record the change as a deviation and resume from §2.4."
3. **The §2.5 comparison ran twice.** The first run showed `results` differing only in timing fields. The second appended a comparison that drops the timing fields. Each run was under the 10-min cap, and the numbers above come from the second run.
4. **Render logs got new names** (`render_*_math.log`) so the previous run's logs were not overwritten.

### OPEN ITEMS

1. **The cost-table caption is now wrong.** It still says "Single-core timings from this run. The projection is linear in B." The bootstrap column now runs on `n_workers` = 2 here, so "s per replicate" and the B = 1,000 projection are two-worker figures. Proposed fix, for Larry's call:
   - say the bootstrap column is on `n_workers` workers; or
   - divide by the plan's worker count before projecting.

   The v5 text was placed unchanged.
2. **"0% or less" in the recovery sentence.** It renders as "a tenth of the draws retain 0% or less". It is accurate, since q10 = 0, but awkward. Possible wording: "at least a tenth of the draws retain none of it".
3. **The multiplier-stage seconds.** The rendered cost prose now says what each number measures (the whole multiplier stage, `timing_seconds`, against the whole bootstrap call). This closes item 2 of the previous pass.

## Copies to `~/Downloads/vignettes_math_2026-10-28/`

`forestsearch.html`, `survival-post-selection.html`, `REPORT_postsel_vignette_2026-10-28.md`.

## Caption and recovery-sentence fix (2026-10-28)

Task: `dev/tasks/cc_task_postsel_caption_2026-10-28.md` (commit `8de6cdaf`). Base HEAD `25318e8c`. Vignette commit `e3e18b07`.

### Gates (all passed)

- **§0.1:** branch `feature/vignettes-0.4.1`; `25318e8c` is an ancestor of HEAD. The two pre-existing untracked article `.html` files were left unstaged.
- **§0.3:** the input md5 is `660b2f87bd26edd0a5e1fdd9801fd09e` (499 lines; line 447 has 214 characters). The base `vignettes/survival-post-selection.qmd` md5 is `8c56077e1fdbd0b1c13fbf4c5d082e15`.
- **§0.4:** `git diff 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` is empty. `future::availableCores()` = 128, and `n_workers` resolves to 2.
- **§1.2 (purl):** the purled diff contains only the added `plan_label <-` line (purled line 190) and the four-line `sprintf(paste(…), plan_label)` caption (purled lines 194–197). All of it lies inside `## ----cost` (purled 182) and before `## ----phat` (purled 200). No other chunk differs.
- **§1.3 (markdown):** the numbered-list grep returns nothing; `Single-core` count in the qmd = 0.
- **§1.4 (render):**
  - exit 0.
  - In the plain text: `Error in` 0, `NaN` 0, `Warning` 0, `\bNA \(` 0.
  - In the raw HTML: `gt_table` 8, `mathjax` (case-insensitive) 4.
  - Multiplier strings in the plain text: `1.44 (0.62, 3.36)` 4, `0.64 (0.40, 1.03)` 4, `1.11 (0.45, 2.75)` 3, `0.66 (0.40, 1.09)` 3, `1.05 (0.47, 2.34)` 3, `0.66 (0.41, 1.06)` 3. All are ≥ 2.
  - Raw HTML: `the bootstrap on 2 workers` occurs 1 time; `Single-core` occurs 0 times.
  - Plain text: `0% or less` occurs 0 times; `at least a tenth of the draws retain none of them` occurs 1 time.
  - "Full bootstrap estimate" column, Ĥ rows: FS 1.88, DINA 0.83, GRF 1.29. These are unchanged from the two previous renders.

### Render

- **Wall clock:** 2:11.80 (`/usr/bin/time -v`). Max RSS 986,976 kB.
- **HTML size:** 1,378,503 bytes.
- **Log:** `postsel_vignette_2026-10-28/render_postsel_caption.log`.
- **`boot_secs`:** FS 34.4 / DINA 24.2 / GRF 22.6 s this run, against FS 32.6 / DINA 23.6 / GRF 21.8 s in the previous run. This is wall-clock jitter on two workers. The estimates are unchanged.

### Recorded excerpts (from `pandoc --wrap=none -t plain`)

Cost table and caption:

```
  -------------------------------------------------------------------------------------
  Method     MR (s)    B   bootstrap (s)   s per replicate   projected at B = 1,000 (s)
  -------- -------- ---- --------------- ----------------- ----------------------------
  FS           4.72   25           34.36              1.37                         1374

  DINA         0.35   25           24.22              0.97                          969

  GRF          2.13   10           22.61              2.26                         2261
  -------------------------------------------------------------------------------------

  : Wall-clock timings from this run: the multiplier stage under the fits’ sequential plan, the bootstrap on 2 workers. The projection is linear in B at that worker count. {.caption-top .table .table-sm .table-striped .small}
```

The paragraph after the cost table:

```
The multiplier correction is seconds; the bootstrap is minutes at the B its intervals need, and the ratio widens with B. The seconds are those of this machine and this package version: the MR column is the whole multiplier stage as the fit records it (timing_seconds); the bootstrap column is the whole forestsearch_bootstrap_dofuture() call on 2 workers, in wall-clock seconds, so the per-replicate figure and the projection fall as workers are added. The manuscript’s figures — 72 s to about 10 minutes at B = 1,000 on 11 workers, 0.1 to 0.4 s for the multiplier stage — were recorded under forestsearch 0.2.0 on another machine, and are not comparable second for second; the order of magnitude between the two routes is the point.
```

"Read together" paragraph:

```
Read together: the median draw retains 66% of the region’s members, at least a tenth of the draws retain none of them, and 2% of draws contain all of it. Exact re-selection is rare; what moves from draw to draw is the boundary. These diagnostics measure re-selection within the fixed candidate family under multiplier perturbation; the bootstrap’s and the cross-validation’s recovery rates (forestsearch_bootstrap_dofuture(), forestsearch_Kfold()) measure re-discovery from resampled data with the family rebuilt, which is stronger and costs a full re-run each time.
```

### Deviations

None. The task was executed as written.

### OPEN ITEMS

None. This commit closes both open items of the previous section:

1. The cost-table caption is now built from `n_workers`. It reads "…the bootstrap on 2 workers. The projection is linear in B at that worker count."
2. The recovery sentence now reads "at least a tenth of the draws retain none of them". It no longer says "0% or less".

## Copies to `~/Downloads/postsel_caption_2026-10-28/`

`survival-post-selection.html`, `REPORT_postsel_vignette_2026-10-28.md`.
