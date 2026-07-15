# CLAUDE.md — forestsearch project rules

This file encodes standing instructions for Claude Code when working
in the `forestsearch` repository. It is read automatically at the
start of every session and on every turn. You do not need to repeat
these rules in your prompts.

---

## 1. Source of truth

- The most recent version of a file uploaded or pasted in
  conversation, or the file at the exact path I cite, is the source
  of truth for that file.
- Do not pre-emptively re-read a file from project knowledge,
  uploaded assets, or earlier in the conversation unless I ask for
  it. If you are uncertain which version to work from, ask me.
- Before any file revision, state which version you are building on
  (path + short description) so I can confirm.
- After any successful file edit, treat earlier views of that file
  in your context as stale. Re-read the file before further edits
  to it.

## 2. Scope discipline

- Modify only what is explicitly requested. Do not opportunistically
  fix unrelated issues, "improve" surrounding code, or refactor for
  style while doing other work.
- If you notice a side issue while doing requested work, **flag it
  separately at the end of the turn**. Do not silently fix it. Wait
  for me to decide whether to address it.
- When in doubt about whether a change is in scope, ask.

## 3. Communication conventions

- When proposing multiple options for a decision, present them as a
  short numbered list with a clear recommendation. I will choose.
- For substantive technical claims (e.g.\ "this is a Cox
  non-collapsibility issue"), explain the reasoning briefly. Do not
  just assert.
- If you make a mistake, acknowledge it directly. Do not
  over-apologise. Move to the fix.
- Honest disagreement is welcome. If I propose a change you think
  is wrong, say so and why.

## 4. R package conventions

### File organisation
- Functions live in `R/<name>.R`, one major exported function per
  file, named after the function.
- Helpers used by a single main function live in the same file as
  that function, prefixed with `.`.
- Helpers shared across multiple functions live in a separate file
  with a descriptive name (e.g.\ `R/generate_aft_dgm_helpers.R`).

### Roxygen
- Match the roxygen style of existing siblings in the same file or
  in a closely related file. For new calibration utilities, use
  `R/sim_aft_gbsg.R::calibrate_k_inter()` as the reference.
- Every exported function needs `@title` (first line of the
  description), `@param` for every argument, `@return`, `@examples`
  (wrapped in `\dontrun{}` if they require fitted DGMs or
  simulations), `@seealso` for related functions, and `@export`.
- Use `@importFrom <pkg> <fn>` for any non-base function the
  implementation uses. Never use bare `pkg::fn()` without the
  matching `@importFrom`.

### API patterns
- Functions wrapping `generate_aft_dgm_flex()` should accept a
  `base_args` named list rather than `...`. This makes the call
  site explicit and avoids ambiguity about which arguments are
  being calibrated versus held fixed.
- Inside such wrappers, override arguments via
  `utils::modifyList(base_args, list(<override> = <value>))` — never
  via `c(base_args, list(<override> = <value>))`, which produces a
  duplicate-argument error when `base_args` already contains the
  key being overridden.
- Calibration functions that wrap a DGM constructor should accept a
  `verbose` flag (default `FALSE`) and emit diagnostic `message()`
  calls — not `cat()` — when `TRUE`.

### Validation
- Use `stopifnot()` at the top of exported functions to validate
  argument types, lengths, and basic constraints (positivity,
  ordering, etc.).
- When a function depends on a `seed` field in `base_args` for
  reproducibility, emit a `warning()` if it is `NULL`.

## 5. Vignette conventions

- The package vignettes use Quarto (`.qmd`), not R Markdown
  (`.Rmd`).
- The canonical comprehensive vignette is `extreme_subgroups.qmd`;
  the canonical glossary vignette is `treatment_effect_definitions.qmd`;
  the canonical biomarker reference is `biomarker_effects.qmd`.
- Use `\text{Extreme Value}` (not `\mathrm{Gumbel}(0,1)`) for the
  AFT residual distribution, to match the cross-vignette
  convention. Add the explicit density
  $f_\varepsilon(x) = \exp(x - e^x)$ on first introduction.
- Distinguish carefully between:
  - **Marginal Cox HR** — `dgm$hazard_ratios$overall`, the
    no-covariate Cox fit on stacked potential outcomes. Attenuated
    by non-collapsibility.
  - **Average hazard ratio (AHR)** — `dgm$hazard_ratios$AHR`,
    equals the patient-level individual HR under `model = "null"`.
  - **Individual / patient-level HR** — what each patient sees.
- Do not call the marginal Cox HR the "true HR." It isn't.
- Calibration targets the AHR by default
  (`calibrate_k_treat(..., use_ahr = TRUE)`) because the simulation
  Cox is stratified-by-grade, which estimates the patient-level
  quantity under proportional hazards.

## 6. Rendering and diagnostics

- Use `quarto render <file>.qmd` for vignettes. The
  `embed-resources: true` YAML option produces a single
  self-contained HTML; preserve this setting.
- After rendering, check the resulting HTML for `cell-output-error`
  and `cell-output-stderr` blocks and report any warnings or errors
  before declaring success.
- For data.table operations, prefer `data.table::set()` over
  `dt[[col]] <-` to avoid shallow-copy warnings.
- Long-running simulation chunks should be guarded with a small
  default `n_sims` (e.g.\ 100) and a comment indicating the
  publication-quality value (e.g.\ 10,000 or 20,000). Do not change
  the default upward without asking.

## 7. Git and branch operations

- The active development branch is `feature/glm-extension`. The
  trunk is `master`.
- Commit `R/`, `vignettes/`, `tests/`, `man/`, `NAMESPACE`, and
  `DESCRIPTION` changes. Do not commit `_files/` sidecar
  directories from Quarto renders, `.Rproj.user/`,
  `.RData`/`.Rhistory`, or `inst/doc/` (these should be
  `.gitignore`d).
- Propose commit messages but do not commit without explicit
  approval. Same for `git push`.
- For any operation that rewrites history (`rebase`, `commit
  --amend`, `reset --hard`), ask first and explain why.

## 8. Testing

- For new exported functions, propose a smoke test that confirms
  the function runs end-to-end on a small example. For calibration
  functions, the test should confirm convergence to the target.
- Tests go in `tests/testthat/test-<name>.R`, named after the
  function being tested.
- If a test fails, do not edit the test to make it pass without
  first confirming the failure is in the test (not the function).

## 9. Things to avoid

- Do not run `install.packages()` for non-trivial dependency chains
  without asking first.
- Do not run any operation that exceeds ~30 seconds of wall time
  without warning me first. For long simulation runs, suggest
  running them in a separate terminal and saving results to disk.
- Do not edit `DESCRIPTION` (version bump, dependency add) without
  asking.
- Do not delete files. If something looks unused, propose the
  deletion and wait.

## 10. Style preferences

- When providing files, place them where I can download them
  directly. Always offer to present the file at the end of a turn
  rather than only showing diffs.
- Provide source files for functions under discussion so I can copy
  to my local tree.
- For numbered options or tradeoff lists, give a clear
  recommendation. Default to action, not deferral.
- Prefer concise responses over comprehensive ones unless I ask for
  depth. The chat is iterative; we can drill in if needed.
