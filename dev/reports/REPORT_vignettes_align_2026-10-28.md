# REPORT — vignettes align (2026-10-28)

Task: `dev/tasks/cc_task_vignettes_align_2026-10-28.md`. Branch `feature/vignettes-0.4.1`; HEAD at start
`9567a568`. Record: `dev/reports/vignettes_align_2026-10-28/`. Not pushed.

Commits: `75bdac25` task · `946849ea` §1 · `055a7f7d` §2 · `0a73e7df` §3 · `a1393fe9` §4 · `1a88c295` §5 ·
this report §6.

## Provenance

- pop-os, 128 cores; R 4.6.1; forestsearch 0.4.0.9000 installed, Built 2026-10-27 21:01:28 UTC; quarto CLI
  1.10.18; pandoc 3.10 (`RSTUDIO_PANDOC`); pkgdown 2.2.1, knitr 1.52, rmarkdown 2.32 (`provenance.txt`).
- §0.4 NEWS-identity gate **failed** (installed `NEWS.md` lacks the `### Documentation` block added in
  `0de9a8e5`, 21:02 PDT, after the 14:01 PDT build) and was **waived** by Larry. Replacement gate, PASS:
  `git diff 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` is empty (`936b6734` = last commit before the
  build). Only `man/fit_causal_forest_glm.Rd` deleted since (`de298322`), for a function already absent
  from `R/` at `936b6734`.
- Pre-existing paths: `git status --short` clean. `vignettes/forestsearch.html` ignored;
  `vignettes/articles/forestsearch-2024.html` is **tracked** (Larry, `9567a568`), not untracked as the task
  expected — left as is.

## §1 Primary vignette introduction (v4)

- v4 md5 `fe9925e92b598a0f913af53446032ca7`, 626 lines — match.
- purl identity: HEAD and v4 purled scripts byte-identical (md5 `eb960012dfdc653c7e1d5796161fdac2`).
- Diff confined: hunk 1 at lines 26–114 (before `## The data`, line 116), hunk 2 at 590–610 (inside
  `## References`, line 588).
- Render: exit 0, wall-clock 1:13.49, max RSS 1.11 GB; HTML 1,283,200 bytes; `Error in` / `Warning` /
  `NaN` / `NA (` all 0.
- Table 1 block (`pandoc -t plain`, `Table 1` → `Declaration at a protected level`): identical to the
  committed render, 49 lines, md5 `5ba789c969677b725730a66b87b46841`.

## §2 `survival-post-selection.qmd`

- Edits A–I applied, each old string matched once. `Anderson (2026)` grep: none.
- Settings grep: as specified it matched lines 145/154 (`0.80 / 0.85` reading thresholds) — gate stopped;
  narrowed by Larry (see Deviations); narrowed grep: none.
- Render: exit 0, wall-clock 0:36.98 (fit 32.2 s), max RSS 0.82 GB; HTML 1,223,367 bytes;
  `Error in` / `Warning` / `NaN` 0.
- Rendered (3 dp): region `{er <= 0} & {pgr <= 26}`; unadjusted 2.222 / corrected 1.276 / field lower
  0.634 / field-s upper 0.804; `print(fs)` Bonferroni 0.509 / 0.840; IJ two-sided (0.570, 2.857).
- Full precision (`postsel_full_precision.R`: the vignette's own data + fit chunks, purled, re-run on the
  installed build; record only, not a gate) beside S3:

  | quantity | this fit | S3 |
  |---|---|---|
  | unadjusted H | 2.221839 | 2.221839 |
  | corrected H | 1.276181 | 1.276181 |
  | field lower H | 0.633561 | 0.633561 |
  | field-s upper Hc | 0.803570 | 0.803570 |
  | Bonferroni H / Hc | 0.508684 / 0.839594 | 0.508684 / 0.839594 |
  | corrected Hc | 0.644653 | 0.644653 |

  n_family 1744, field R_out/R_in 1000/500 — both as S3.
- **`field_recovery` reading:** every quantity equals S3 to 6 dp → the recovery diagnostics do not
  consume draws (they do not perturb the RNG stream that builds the bounds).
- No stage-0 counterpart, as rendered: p-hat(Hhat) = 0.005 over 1744 candidates; sens_H 0.572; ppv_H
  0.579; containment q10/q50/q90 0.000 / 0.733 / 1.000; share of draws containing all of Hhat 0.255
  (over 990 draws).

## §3 `methodology.Rmd`, `references.bib`

- Pointer note inserted before `## Introduction`.
- Header comment replaced in both listings (both carried the same "All forestsearch() parameters shown"
  wording).
- Deleted (4 lines): `max_n_confounders = 1000, # default 1000` and `vi.grf.min = -0.2, # default -0.2`,
  in each block.
- `# default` comments changed (both blocks; old → new):
  - `use_lasso`: `TRUE` → `FALSE` (`<-- FOCAL` kept)
  - `use_grf`: `TRUE` → `FALSE` (`<-- FOCAL` kept)
  - `max_subgroups_search`: `10` → `Inf` (value 10 unchanged)
  - `stop_threshold`: `0.95` → `pconsistency.threshold`
  - `outcome_type`: `"survival"` → `c("survival", "binary", "continuous", "count")`
  - `overdispersion`: `"none"` → `c("none", "quasi", "negbin")`
  - `grf_count_transform`: `"log"` → `c("log", "identity")`
  - `ps_adjust_method`: `"none"` → `c("none", "iptw", "dr_gcomp")`
- No listed argument was a non-formal (none deleted on that ground). `parallel_args` comment
  (`default plan="multisession"`) is outside the gated form; left.
- Gate (`methodology_listings_check.R`, installed formals): before 22 findings
  (`logs/methodology_check_before.log`), after ALL CHECKS PASS.
- Stale statements: minimum events → `10` / `d0.min = 10` / `d1.min = 10`; `fs.splits` default 1000 with
  the closed-form proviso. Defaults verified in `formals()` (`fs.splits` 1000, `d0.min`/`d1.min` 10,
  `consistency_method` first choice `"resample"`).
- Bib: neither arXiv ID present → appended `Leon_Anderson_2026a`, `Leon_Anderson_2026b` (after the
  trailing `@Comment{jabref-meta}`); both keys added to `nocite:` after `@Wager_2014`.
- Render: exit 0, wall-clock 0:03.03; HTML 3,340,667 bytes, untracked; References lists both 2026a
  (`10.48550/arXiv.2609.38361`) and 2026b (`10.48550/arXiv.2610.28379`).

## §4 `paper_simulations.Rmd`

- Pointer note inserted before `# Introduction`; `vi.grf.min = -0.2,` deleted from `define-params`.
- Grep `vi.grf.min|max_n_confounders|max.minutes`: one hit, line 55 — prose (the pointer note). None in an
  evaluated chunk.
- Render: exit 0, wall-clock 0:46.76 (within the 30-min cap), max RSS 0.50 GB (parent); 120 workers
  (`setup-parallel` unchanged); n_sims 44 (H1) + 44 (H0) as rendered; HTML 3,835,980 bytes, untracked; no
  deprecation message in the log.

## §5 NEWS

- Bullet appended under `### Documentation` in the development section.

## Deviations (verbatim)

- §0.4: "The §0.4 NEWS-identity gate is waived (Larry, 2026-10-28): the installed build is HEAD's R/ tree —
  record your 936b6734 check (no R/, DESCRIPTION, NAMESPACE, inst/ change since) in provenance.txt as the
  replacement gate and list the waiver under deviations; do not reinstall."
- §0.1: "vignettes/articles/forestsearch-2024.html is tracked by 9567a568 — leave it as it is."
- §2 gate: "Narrow the §2 grep gate to the pattern 'pconsistency.threshold = 0.80\|fs.splits\|draws = 500L'
  — the 0.80/0.85 reading thresholds on lines 145 and 154 stay as they are."
- §2 record: an extra ~33 s re-run of the vignette's fit (`postsel_full_precision.R`) to read the S3
  comparison at 6 dp; the render prints 3 dp. Not a gate.

## OPEN ITEMS

- `methodology.Rmd` lines 955 and 1064, the second comment line of each listing: "LASSO is enabled by
  default" and "Enable both LASSO and GRF (the package default for variable selection)". Both
  `use_lasso` and `use_grf` default to `FALSE`. Not in the task's edit list; left.
- `methodology.Rmd`: the four `match.arg` comments now show the full choice vector (that is the literal
  formal default). The effective default is the first element; if a reader-facing form such as
  `"survival"` is preferred, the check's definition must change with it.
- `paper_simulations.Rmd` `define-params` keeps `use_grf = TRUE`, `fs.splits = 400`,
  `d0.min = d1.min = 12` as 2024-workflow values (the note covers `fs.splits`).

## ~/Downloads/vignettes_align_2026-10-28/

```
total 9680
drwxrwxr-x 2 larryleon larryleon    4096 Oct 28 09:21 .
drwxr-xr-x 7 larryleon larryleon   20480 Oct 28 09:20 ..
-rw-rw-r-- 1 larryleon larryleon    1751 Oct 28 09:20 diffs.txt
-rw-rw-r-- 1 larryleon larryleon 1283200 Oct 28 09:20 forestsearch.html
-rw-rw-r-- 1 larryleon larryleon 3340667 Oct 28 09:20 methodology.html
-rw-r--r-- 1 larryleon larryleon   54158 Oct 28 09:20 methodology.Rmd
-rw-rw-r-- 1 larryleon larryleon 3835980 Oct 28 09:20 paper_simulations.html
-rw-rw-r-- 1 larryleon larryleon  107906 Oct 28 09:20 paper_simulations.Rmd
-rw-rw-r-- 1 larryleon larryleon    7446 Oct 28 09:21 REPORT_vignettes_align_2026-10-28.md
-rw-rw-r-- 1 larryleon larryleon 1223367 Oct 28 09:20 survival-post-selection.html
-rw-r--r-- 1 larryleon larryleon   10057 Oct 28 09:20 survival-post-selection.qmd
(listing taken before this section was filled in; the REPORT size shown is the pre-listing size)
```
