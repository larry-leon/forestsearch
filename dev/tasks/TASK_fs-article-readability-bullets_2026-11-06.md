# TASK_fs-article-readability-bullets_2026-11-06

**Repo:** forestsearch. **Branch:** `fs-article-readability-bullets` from `master`. **Touches:**
`vignettes/articles/mrct-consistency.qmd` and its tracked render
`vignettes/articles/mrct-consistency.html` (re-render) only. **No `R/`.** **No simulation
compute:** render-only; expect roughly 10 minutes.

## Context (orientation, not evidence)

Seven readability restructures R1–R7, approved by Larry from a rendered preview (2026-11-06):
dense paragraphs become bulleted or numbered lists in the style of the committed (M1)–(M4) passage.
No numeric content changes; all values remain inline-R expressions; R5 only inserts paragraph
breaks and moves one sentence. Every anchor was verified verbatim against the current render.

**Companion file:** `EDITS_fs-article-readability-bullets_2026-11-06.md` — the seven old → new
blocks.

## First action

Copy this task doc and the companion from `~/Downloads` into `dev/tasks/` and commit them together.
Run `export GIT_TERMINAL_PROMPT=0` at session start. `master` may be ahead of `origin/master`
(Larry pushes); not a stop.

## Legacy copies (standing instruction, before any edit)

- `d="legacy_vignette/$(date +%d%b%Y)"`; if `$d` or a suffixed sibling already holds either file,
  use the next free directory (`${d}_2`, `${d}_3`, …) — never overwrite an existing copy
  (`06Nov2026/` and `06Nov2026_2/` exist from today's earlier tasks).
- Copy `vignettes/articles/mrct-consistency.qmd` and `vignettes/articles/mrct-consistency.html`
  into it as they stand at task start; verify each with `cmp`.
- Untracked and local: `legacy_vignette/` stays in `.git/info/exclude`; never touch `.gitignore`,
  never `git add` under it. Record the copy paths in NOTES.

## Edit-application convention

Locate each `old` block with whitespace-normalized matching (any run of spaces/newlines ≡ one
space); each must occur exactly once. Write the `new` text with its list structure as given, prose
lines re-wrapped at the file's prevailing width (~100 characters). Nothing outside the seven blocks
changes.

## Gates (stop-on-failure)

- **G0 — state.** Working tree clean on `master`; record HEAD in NOTES.
- **G1 — anchors.** Every `old` block matches exactly once. On any miss: STOP; NOTES names the
  block and pastes the passage as it now stands.
- **G2 — apply.** `git diff` confined to the seven passages. Byte-identical in particular: R1's
  following "The question is whether …" sentences; R2's companion-link sentence; R6's two lead
  sentences; R7's "Multiplier resampling (MR, …)" sentence and the *With your trial* aside.
- **G3 — code integrity.** Every code chunk byte-identical (prose-only task).
- **G4 — render.** `quarto render` green; every `stopifnot` passes.
- **G5 — render delta.** Against the tracked baseline HTML: outside the seven passages and timing
  positions, numeric tokens identical; inside them, additions are the companion's verbatim text
  only — no computed value changes (the bullets repeat threshold/label constants such as 0.90,
  0.80, 90%, 0.55, 2, 3, 20%, 1,000, 40, 11, 230, 0.20, 1.0). Each restructured passage renders as
  a list (`<ul>`/`<ol>` present at R1, R2, R3, R4 (two lists), R6, R7); R5 renders as three
  paragraphs.
- **G6 — discipline.** The multiset of inline-R expressions in the file is identical before and
  after (R5 moves its inline-R sentences intact; no block adds or removes any); no manuscript-value
  comparison; no provenance note.

## Close-out

- Commit the edits and the re-rendered tracked HTML; fast-forward merge to `master`; delete the
  branch (`git branch -r` checked first). Larry pushes.
- `git status --porcelain` empty at close.
- `NOTES_fs-article-readability-bullets_2026-11-06.md` to `~/Downloads`: HEADs, gate outcomes, the
  G5 delta listing, the legacy-copy paths.
