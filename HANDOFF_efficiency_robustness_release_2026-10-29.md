# HANDOFF — efficiency routes closed, robustness items in, 0.4.1 release branch ready (2026-10-29)

**For:** the next forestsearch chats — first the continued
vignette/article revisions, then the win-builder results and the 10
November submission. **Read with**
`HANDOFF_vignettes_release_0.4.1_2026-10-28.md` (references, rulings,
task-doc lessons, where things live — all still hold) and its `_v2`
(vignette state and the release plan, of which §1 and §3 are superseded
below). **Written by** chat at the close of the efficiency-routes chat
(2026-10-28/29), from the CC reports of that chat and the GitHub
`master` read read-only.

------------------------------------------------------------------------

## 1. Where things stand (2026-10-29, late morning)

- **`master` = `624bd457`** (pushed by Larry): the vignette work (one
  CRAN vignette `vignettes/forestsearch.qmd`; `survival-post-selection`
  as a pkgdown article), the lean fits A/F, the deprecations, the `xlog`
  fix, plus today’s three merges — robustness
  (`feature/robustness-0.4.1`, 10 commits), example timing
  (`feature/example-timing-041`, 5 commits), the README badge removal,
  and two check logs committed by CC.
- **`release/0.4.1` = `257204a0`** (+ a pending one-line
  `cran-comments.md` edit “about 50 s” → “about 40 s”): version 0.4.1;
  NEWS heading; the three roxygen fixes; lifecycle badge `stable`;
  `cran-comments.md` for 0.4.1 with the two win-builder lines left for
  Larry. **Final `R CMD check --as-cran --run-donttest`: Status OK,
  0/0/0**, 12 min 24 s (examples 41 s, tests 410 s with `NOT_CRAN=true`,
  vignette 73 s); the tarball ships only the `.qmd` vignette source.
  Record: `dev/reports/REPORT_release_041_2026-10-29.md`, log
  `dev/reports/release_041_2026-10-29/check_as_cran_release_2026-10-29.log`;
  the built tarball is in `~/Downloads`.
- **Remaining for 0.4.1 (Larry):** `devtools::check_win_release()` /
  `check_win_devel()` on `release/0.4.1` (queue took hours on 10-28; run
  early); fill the two win-builder lines of `cran-comments.md`; merge
  `release/0.4.1` → `master` in GitHub Desktop and push; **10
  November:** `devtools::submit_cran()` from `master`, confirm the email
  link; after acceptance tag `v0.4.1`, bump to `0.4.1.9000` with a new
  `# forestsearch (development version)` NEWS heading.
- `pandoc` 3.1.3 is now installed on the Pop!\_OS machine
  (`/usr/bin/pandoc`), so a CC-run check is like-for-like with an
  RStudio check (no README/NEWS NOTE).

## 2. Release-content freeze — what the vignette chat may and may not touch before 10 November

- **In the CRAN tarball** (any change re-opens the release: the final
  check must be re-run on the release branch and `cran-comments.md`
  re-stated): `vignettes/forestsearch.qmd`, `R/`, `man/`, `DESCRIPTION`,
  `NAMESPACE`, `NEWS.md`, `README.md`, `inst/`, `tests/`. **Treat these
  as frozen for 0.4.1** unless a defect is found; if the primary
  vignette must change before the 10th, the change goes on
  `release/0.4.1` (or on `master` before the release merge, followed by
  a re-merge) and ends with the one as-cran check.
- **Not in the tarball** (`.Rbuildignore`: `^vignettes/articles$`,
  `^_pkgdown\.yml$`, `^docs$`, `^quarto$`, `^dev$`, `^HANDOFF_.*\.md$`):
  the pkgdown articles, `_pkgdown.yml`, the quarto documents, dev
  records, handoffs. **Article revisions can proceed on a feature branch
  now and merge any time**; they reach the site through the pkgdown
  workflow on `master` and never touch the submission. Standing rules
  for them are in the pkgdown-articles notes: `.qmd`, based on the CRAN
  release, small simulation counts with the “real applications run
  5,000–10,000” note, rendered HTML kept locally and never staged by a
  task.
- Vignette/article items carried over, none blocking: the primary
  vignette’s “about a minute and a half on one core” against a measured
  1:14 (cosmetic; **in the tarball**, so it waits for the next release
  unless the release branch is reopened anyway); the `.gitignore`
  decision for `vignettes/articles/*.html` (`forestsearch-2024.html`
  tracked, three renders untracked and not ignored); the ACTG175 GLM
  article follows the release; the v1 handoff §6 items 2–11.

## 3. The efficiency question — closed for 0.4.1 (decision memo `decision_efficiency_routes_0.4.1_2026-10-28.md`, chat-side; the project KB was out of room)

- Re-based on the post-lean replicate (lean-fits Stage C: survival 2.02
  s, continuous 1.84 s, binary 2.75 s), the deferred consistency-stage
  items are **B+H 20.4% (survival), E+H 10.9% (continuous), G+H 8.0%
  (binary)** — all below the 30% bar; C (`multicore`) 2.6–8%; D
  (chunking) ≤ 14.6%; medians 5–6% (declined 09-02). High confidence:
  direct-timer method, which predicted A and F to within 0.07 s.
- Recorded, not proposed: B+H could clear only on a campaign shape with
  ≈ 210+ candidates entering the consistency stage (applications:
  120–136) — unmeasured; a zero-`R/` one-call measurement under a
  null/uniform-benefit DGM at c1 0.90 / c2 0.80 / p\* 0.90 would settle
  it, after the release. Search-stage bookkeeping (ceilings 19–30%,
  realizable share unmeasured, touches enumeration order — not for a
  patch release). MR inference is outside the survey’s scope and is, by
  subtraction, most of the vignette’s 74 s and of the application
  documents.
- Lesson kept: if B/E/G/H are ever done, do them as one “consistency
  pieces on vectors” refactor sharing the lean cores, not four items.

## 4. What changed in `R/` since the 10-28 v2 handoff, and why no analysis moves

| change | category | evidence of no change to conducted analyses |
|----|----|----|
| `survival (>= 3.5-8)` in `Imports` | DESCRIPTION only | no code runs differently |
| continuous IPTW weights fix (`.make_lm_estimator()`: `weights = sw`; `globalVariables("sw")`) | **changes behaviour** in one branch: `ps_adjust_method = "iptw"` + continuous, which failed every fit since 0.3.2 (`REPORT_fastpath_unadjusted_2026-09-02.md` F6) | branch unreachable in every committed analysis (`is.RCT = TRUE` throughout; grep of `quarto/applications`, `vignettes`, `tests`, `inst` finds only the new regression test); `else` branch byte-identical; the applications run the fast closure, a different function |
| lean Cox core `.fs_cox_two_group()` + once-per-session probe `.fs_cox_lean_ok()` + `.fs_cox_lean_conventions()` (ties/nocenter read from `coxph()`’s formals) | moves code; identical on the normal path; formula fallback with one message if `coxph.fit()` ever stops reproducing `coxph()` | Stage A (unmodified master) = committed lean-fits records = Stage C on every component, 3 fixtures, end to end and `nb = 20` bootstraps; `test-lean-fits.R` extended (probe passes, fallback identical) |
| [`run_null_calibration()`](https://larry-leon.github.io/forestsearch/reference/run_null_calibration.md) example at illustrative size (2 N values, 10 sims, no LASSO/GRF) | roxygen only | 51 s → 3.3 s; function untouched (it requires `n_sims >= 10`, two N values) |
| [`suggest_thresholds()`](https://larry-leon.github.io/forestsearch/reference/suggest_thresholds.md) example: binary call shown, not run | roxygen only | survival call live (3–4 s); binary call ~10 s was the last \>5 s example |
| [`mr_estimates_table()`](https://larry-leon.github.io/forestsearch/reference/mr_estimates_table.md) roxygen (`ci_method` default), [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md) example comment (article, not vignette), README badges (R-CMD-check removed — the workflow never existed, 404 at CRAN’s URL check; lifecycle `stable`) | documentation | — |

Records: `dev/reports/REPORT_robustness_041_2026-10-28.md` (+
`robustness_041_2026-10-28/`), `REPORT_example_timing_2026-10-29.md` (+
`example_timing_2026-10-29/`), `REPORT_release_041_2026-10-29.md` (+
`release_041_2026-10-29/`); task docs under `dev/tasks/` with the same
stems.

## 5. Open items after the release (recorded; a task only when Larry commissions)

1.  The consistency-stage pieces
    ([`.consistency_glm_pieces()`](https://larry-leon.github.io/forestsearch/reference/dot-consistency_glm_pieces.md),
    `R/consistency_resample.R`) take no weights for any
    `ps_adjust_method`: an IPTW analysis is weighted in the search stage
    and unweighted in the consistency screen (binary OR too). A method
    question.
2.  At `is.RCT = FALSE` with defaults (`ps_method = "grf"`,
    `ps_adjust_method = "none"`) a propensity forest is fitted per call
    and per bootstrap replicate and never consumed, and the rebuilt
    estimator closure is untagged (lean path lost). Also a stale design
    comment was corrected (`forestsearch_main.R` PS block: default is
    `"none"`).
3.  [`subgroup_glm()`](https://larry-leon.github.io/forestsearch/reference/subgroup_glm.md)
    (`R/run_subgroup_sims.R`) routes only the continuous case to a fast
    closure; one constructor shared with
    [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
    would close the drift (no committed consumer).
4.  With `use_lasso = TRUE` under a null DGM, replicates where LASSO
    selects a single covariate end in
    `get_FSdata: no candidate cut expressions … returning NULL result`
    (seen in the old `run_null_calibration` example; caught by the
    function; invisible in the fixed-family configuration). Looks like a
    defect of the LASSO-screened path when exactly one variable
    survives.
5.  The campaign-shape measurement of §3 (B+H).
6.  Hard-coded `1.96` in the GLM search-stage CI vs `qnorm(0.975)` on
    the Cox path; fixed 95% — a flexibility note; changing it moves
    displayed GLM bounds at the 4th digit.
7.  `.gitignore` for `vignettes/articles/*.html`; the primary vignette’s
    “minute and a half” sentence; v1 handoff §6 items 2–11.

## 6. Protocol clarifications from this chat (apply going forward)

- **Check policy, refined (Larry, 2026-10-29):** the concern was
  needless repetition, not CC running checks. CC runs no check mid-task
  or per edit; the full `R CMD check --as-cran` is the **last gate of
  any branch touching `R/`** (once, before Larry merges) and of the
  release task, and may be run by CC. Interim verification stays with
  the task’s own acceptance tests and machine-checkable post-conditions.
- Closing-block greps must exclude comment lines **and string literals**
  (CC’s open item (d) after three self-inflicted mismatches in the
  robustness task); never count a pattern the task document itself told
  CC to write into a comment or a message.
- `R CMD check` on a plain shell needs pandoc on `PATH` — now installed;
  `RSTUDIO_PANDOC` alone leaves the README/NEWS NOTE.
- CC stopping at a wrong-branch / unmerged-base gate and reporting is
  correct behaviour; the fix is Larry’s merge, then re-issue.
- A measurement stage (Stage A on the unmodified tree, installed first)
  is the baseline for any `R/`-touching task; identity is
  saved-against-saved, component by component, Stage C vs Stage A, with
  the committed earlier records compared informationally.
- Behaviour-invariance claims are stated in the task document’s category
  block and proven in the report from the evidence, never asserted from
  description.

## 7. Where everything lives

- Repo root: this handoff;
  `HANDOFF_vignettes_release_0.4.1_2026-10-28.md` and `_v2.md`;
  `HANDOFF_bootstrap_efficiency_2026-09-02.md`.
- `dev/tasks/`: `cc_task_robustness_041_2026-10-28.md`,
  `cc_task_example_timing_2026-10-29.md`,
  `cc_task_release_041_2026-10-29.md`. `dev/reports/`: the three reports
  and record directories above;
  `robustness_041_2026-10-28/check_as_cran_2026-10-29.log` (branch
  check, 3 NOTEs, pre-pandoc) and
  `example_timing_2026-10-29/check_as_cran_2026-10-29b.log` (master, 2
  NOTEs: version; `suggest_thresholds`; README 404 — all since
  resolved).
- Chat-side only (KB out of room):
  `decision_efficiency_routes_0.4.1_2026-10-28.md` (file card in the
  efficiency chat); `survey_bootstrap_efficiency_revisit_2026-10-25.md`,
  `REPORT_bootstrap_stage0*_2026-10-26.md` are in the KB.
- `~/Downloads`: the report/log/cran-comments copies and
  `forestsearch_0.4.1.tar.gz`.

## 8. Protocol carried forward

Three layers: Larry decides and approves · chat writes specs and task
documents · Claude Code does all repo operations. Task documents travel
`~/Downloads` → `dev/tasks/`, committed first and alone; kickoffs short,
inline, fenced, the absolute final content of a delivery; “clear CC
first, then paste” with every kickoff. Nothing points CC at project-KB
paths. Any task touching `R/` is called out as moves code / changes
behaviour / changes method. Gates stop the task; they do not stop to
ask. Compute is a separate go/no-go with wall-clock stated from a
measured record. Verify from source, never from description. Findings go
in the record; a task is proposed only when something blocks or Larry
commissions it. Larry pushes via GitHub Desktop; CC never pushes or
merges. Committed work is not re-run or re-verified. Text after `>` or
`❯` in a transcript is CC, not Larry.
