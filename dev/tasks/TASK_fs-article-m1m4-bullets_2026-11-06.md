# TASK_fs-article-m1m4-bullets_2026-11-06

**Repo:** forestsearch. **Branch:** `fs-article-m1m4-bullets` from `master`. **Touches:**
`vignettes/articles/mrct-consistency.qmd` and its tracked render
`vignettes/articles/mrct-consistency.html` (re-render) only. **No `R/`.** **No simulation
compute:** render-only; expect roughly 10 minutes.

## Context (orientation, not evidence)

One approved edit (rendered preview v2, approved by Larry 2026-11-06): the (M1)–(M4) assumptions
paragraph in §1.2 "The data-generating mechanism" becomes a bulleted list, with one added sentence
inside the (M2) bullet describing the region-dependent extension. The γ display paragraph that
follows ("Writing the effect with treatment-by-region interactions …") is untouched. The (M1)–(M4)
tags are kept verbatim — "(M2)" is cross-referenced by that following paragraph.

**Companion file:** `EDITS_fs-article-m1m4-bullets_2026-11-06.md` — the single old → new block.

## First action

Copy this task doc and the companion from `~/Downloads` into `dev/tasks/` and commit them together.
Run `export GIT_TERMINAL_PROMPT=0` at session start. `master` may be ahead of `origin/master`
(Larry pushes); that is not a stop.

## Legacy copies (standing instruction, before any edit)

- `d="legacy_vignette/$(date +%d%b%Y)"`; if `$d` already contains either file (today's earlier task
  wrote `legacy_vignette/06Nov2026/`), use the next free suffix instead: `${d}_2`, `${d}_3`, …
  Never overwrite an existing copy.
- `mkdir -p` the chosen directory; copy `vignettes/articles/mrct-consistency.qmd` and
  `vignettes/articles/mrct-consistency.html` into it as they stand at task start; verify each with
  `cmp`.
- The directory stays untracked and local: `legacy_vignette/` is already in `.git/info/exclude`
  (add it if missing; never touch `.gitignore`, never `git add` under it).
- Record the copy paths in NOTES.

## Edit-application convention

Locate the companion's `old` block with whitespace-normalized matching (any run of
spaces/newlines ≡ one space); it must occur exactly once. Write the `new` text with its list
structure as given, prose lines re-wrapped at the file's prevailing width (~100 characters).
Nothing else changes.

## Gates (stop-on-failure)

- **G0 — state.** Working tree clean on `master`; record HEAD in NOTES.
- **G1 — anchor.** The `old` block matches exactly once. On a miss: STOP; NOTES pastes the
  paragraph as it now stands.
- **G2 — apply.** `git diff` shows changes in this paragraph only; the following
  "Writing the effect with treatment-by-region interactions" paragraph is byte-identical.
- **G3 — code integrity.** Every code chunk byte-identical (prose-only task).
- **G4 — render.** `quarto render` green; every `stopifnot` passes.
- **G5 — render delta.** Against the tracked baseline HTML: outside the edited paragraph and
  timing positions, numeric tokens identical; inside it, the additions are the companion's
  verbatim text (the only numeric-bearing addition is the "0" of the γ null) — no computed value
  changes anywhere. The assumptions render as a list (`<ul>`/`<li>` present at that passage).
- **G6 — discipline.** No new inline R; no manuscript-value comparison; no provenance note.

## Close-out

- Commit the edit and the re-rendered tracked HTML; fast-forward merge to `master`; delete the
  branch (`git branch -r` checked first). Larry pushes.
- `git status --porcelain` empty at close (legacy copies local via `.git/info/exclude`).
- `NOTES_fs-article-m1m4-bullets_2026-11-06.md` to `~/Downloads`: HEADs, gate outcomes, the G5
  delta listing, the legacy-copy paths.
