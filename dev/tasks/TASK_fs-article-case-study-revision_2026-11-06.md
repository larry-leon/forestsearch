# TASK_fs-article-case-study-revision_2026-11-06

**Repo:** forestsearch. **Branch:** `fs-article-case-study-revision` from `master` (after the
fs-articles-wrappers merge; G0 verifies). **Touches:** `vignettes/articles/mrct-consistency.qmd` and
its tracked render `vignettes/articles/mrct-consistency.html` (re-render) only. **No `R/`.**
**No simulation compute:** render-only; the article's in-call MR draws (step 3 and the §2.7 one-call
check) are its normal render cost — expect roughly 10 minutes. Before any edit, untracked local
backup copies of both files are made under `legacy_vignette/` (see Legacy copies below).

## Context (orientation, not evidence — verify from source)

Five prose edits adopted from an external (CoPilot) review of the application article, disposed by
Larry 2026-11-06; the rendered preview was approved in chat. The edit anchors in the companion file
were drafted against the article as it stood before the 2026-11-06 wrapper migration
(fs-articles-wrappers). The five target paragraphs are prose that the wrapper tasks did not list as
changed, but nothing is assumed: anchors apply as exact matches and any miss is a STOP, not a
judgment call.

**Companion file:** `EDITS_fs-article-case-study-revision_2026-11-06.md` — the five verbatim
old → new blocks. It travels with this doc.

## First action

Copy this task doc and the companion EDITS file from `~/Downloads` into `dev/tasks/` and commit them
together. Run `export GIT_TERMINAL_PROMPT=0` at session start (a credentials-failed
`pull`/`fetch`/`ls-remote` is recorded, not a stop).

## Legacy copies (before any edit — Larry's standing instruction, 2026-11-06)

Local copies of vignettes being re-written are kept in the repo, untracked, in date-named
subdirectories:

- At the repo root: `d="legacy_vignette/$(date +%d%b%Y)"; mkdir -p "$d"` (e.g. `07Nov2026`).
- Copy the files this task overwrites, as they stand at task start:
  `cp vignettes/articles/mrct-consistency.qmd vignettes/articles/mrct-consistency.html "$d"/`.
- Verify each copy against its source with `cmp` (byte-identical).
- Keep the directory **untracked and local**: append `legacy_vignette/` to `.git/info/exclude`
  (the local-only ignore). Do **not** modify `.gitignore` and do not `git add` anything under
  `legacy_vignette/`.
- Record the copy paths in NOTES.

## Edit-application convention

Locate each `old` block with whitespace-normalized matching (any run of spaces/newlines ≡ one
space); it must occur exactly once. Write the `new` text re-wrapped at the file's prevailing width
(~100 characters). Nothing outside the five blocks changes.

## Gates (stop-on-failure; a green run needs no chat round trip)

- **G0 — state and merge.** Working tree clean on `master`; record HEAD in NOTES.
  `git grep -q "mrct_trial_analysis" -- vignettes/articles/mrct-consistency.qmd` must succeed
  (the one-call section present ⇒ the wrappers merge landed on master). STOP if it fails or the
  tree is dirty.
- **G1 — anchors.** Every `old` block in the companion file matches exactly once under the
  convention above, and `fit_reg` is step 2's fit object (`grep -n "fit_reg" …qmd`). On any miss:
  STOP; NOTES names the edit and pastes the paragraph as it now stands.
- **G2 — apply.** Replace each old → new; `git diff` confirms no change outside the five edits.
- **G3 — call integrity.** The step-3 `forestsearch()` call chunk is byte-identical before and
  after (prose-only task).
- **G4 — semantics (STOP on mismatch).** From the package's own documentation (roxygen in `R/` or
  `man/`), confirm that `hr.consistency` with `pconsistency.threshold` means: the hazard-ratio
  criterion must hold at `hr.consistency` **in both halves of random data splits, in at least
  `pconsistency.threshold` of the splits**. If the documented meaning differs, STOP — edit E3 needs
  Larry. Subgroup-size / per-arm event minimums: if the identifier documents such defaults, keep
  E3's minimums sentence; if it does not, omit that one sentence and record the finding in NOTES
  (this either-way action is specified, not a gate failure).
- **G5 — render.** `quarto render vignettes/articles/mrct-consistency.qmd` green; every
  `stopifnot` (the §2.7 identity check included) passes; the step-2 CI appears in the HTML.
- **G6 — render delta.** Against the tracked baseline HTML: outside the five edited paragraphs and
  timing sentences (seconds/minutes), numeric tokens identical; inside them, the only value-bearing
  addition is the step-2 CI pair, which must equal `fit_reg`'s printed lower/upper from the step-2
  chunk output. Any other numeric change: STOP.
- **G7 — discipline.** The `.qmd` diff introduces no manuscript-value comparison and no
  computational-provenance note; the only new inline-R expressions are the two
  `fit_reg$conf.int` reads in E2.

## Close-out

- Commit the edits and the re-rendered tracked `mrct-consistency.html`; fast-forward merge to
  `master`; delete the branch (`git branch -r` checked first). Larry pushes.
- `git status --porcelain` empty at close: the `legacy_vignette/` copies stay local via
  `.git/info/exclude`; nothing under that directory is committed.
- `NOTES_fs-article-case-study-revision_2026-11-06.md` to `~/Downloads`: HEADs before/after, gate
  outcomes, the G4 finding (documented semantics wording; whether the minimums sentence was kept),
  the G6 delta listing, the legacy-copy paths, and — as a finding only, no action — the current
  state of `vignettes/articles/_pkgdown.yml` with respect to the two mrct articles.
