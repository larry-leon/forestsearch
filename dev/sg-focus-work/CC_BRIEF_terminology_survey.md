# CC BRIEF — terminology survey: "gate" and "tier" (READ-ONLY, SURVEY ONLY)

```
claude "Read dev/terminology-work/CC_BRIEF_terminology_survey.md and execute it."
```

**This brief changes nothing.** It produces an inventory and a *proposed* rename
map for maintainer approval. A second brief executes whatever is approved.

---

## 1. Objective

The package uses "gate" for at least four distinct things and "tier-1"/"tier-2"
for two named procedures. The maintainer's decisions:

* **"tier-1" becomes "full bootstrap (FB)"** and **"tier-2" becomes "multiplier
  resampling (MR)"**. These are the established manuscript abbreviations and say
  what the procedures are; "tier" encodes only an ordering.
* **"gate" is to be purged except where it is accurate.** See §3 --- one of the
  four senses genuinely is a gate and should keep the word.

Inventory every occurrence, classify it, and propose a mapping. Do not rename.

---

## 2. Hard constraints

* Modify nothing outside `dev/terminology-work/`.
* Hash the package sources before and after; report any change. Reuse
  `fs_hash_sources()` / `fs_guard_verify()` from
  `dev/efficiency-eval/R/00_guard.R`.
* `git status` at the end must show changes only under `dev/terminology-work/`.
* No `sed`, no bulk substitution, not even as a dry run that writes files.

---

## 3. The vocabulary this is measured against

Four things have been conflated. The distinction matters because a blanket purge
would delete a word that is doing real work.

| term | what it is | genuinely a "gate"? |
|------|-----------|--------------------|
| **Admissibility criteria** | which candidates are eligible: `hr.threshold`, `pconsistency.threshold`, `n.min`, `d0.min`, `d1.min`, `maxk`. Part of the identifier's definition. | **no** --- calling a threshold a gate makes `n.min` a gate too, which is not a useful sense |
| **Selection rule** | which admissible candidate is chosen: `sg_focus` | **no** |
| **Post-selection inference** | FB and MR de-biasing and interval estimation, applied to a completed analysis | **no** --- this is the bulk of what `debias_gate` currently does |
| **Harm confirmation** | the decision that a de-biased estimate still indicates harm | **yes** --- a pass/fail applied to a result |

So `debias_gate` names the fourth while chiefly performing the third, and
`gate_draws` counts MR draws rather than anything belonging to a gate. The `gate`
argument of `fs_debias_gate()` (values such as `"point"`) does appear to name the
harm-confirmation decision and may be correct as-is --- **verify rather than
assume.**

---

## 4. What to inventory

Every occurrence of `gate`, `Gate`, `dg` (as an identifier prefix or suffix),
`tier`, `Tier`, `tier1`, `tier2`, `t1`, `t2` across:

`R/`, `man/`, `NAMESPACE`, `tests/`, `vignettes/`, `quarto/` (including the
simulation `.qmd` files under `quarto/simulations/gbsg_redux/`), `NEWS.md`,
`README*`, `DESCRIPTION`, `_pkgdown.yml` if present, and `dev/`.

Beware false positives: `t1`/`t2` and `dg` will match unrelated identifiers.
Report matches you judged unrelated rather than silently dropping them.

For each occurrence record: file, line, the identifier or phrase, which of the
four senses in §3 it belongs to, and this classification:

| class | meaning | rename risk |
|-------|---------|-------------|
| **API** | exported name, formal argument, or list element a user sets or reads | **breaks user code** |
| **Internal** | non-exported function, local variable, internal list field | safe |
| **Docs** | roxygen text, vignette or `.qmd` prose, `man/` | safe |
| **History** | `NEWS.md` entries for already-released versions | **do not rewrite** |

---

## 5. Specific questions to answer

1. **Full API surface.** Which `gate`/`tier` names can a user set or read?
   Known starting points: `debias_gate`, `debias_gate_args`, `gate_draws`,
   `fs_debias_gate()`, `out$debias_gate`, and the `gate` / `t_gate` arguments.
   Find the rest, including anything in returned-object structure that user code
   or the simulation `.qmd` files index into.
2. **Which occurrences are FB and which are MR?** "tier-1" should map to FB and
   "tier-2" to MR, but confirm each site refers to the procedure the label
   claims. Report any that are ambiguous or mislabelled.
3. **Does `fs_debias_gate()` do FB, MR, or both?** The name suggests one thing;
   `debias_gate_args$draws` and `ci_method = "ij"` suggest MR. If the function
   spans both, a single rename cannot be right and the finding is structural.
4. **Is the `gate` argument accurately named?** Read what it does. If it selects
   among harm-confirmation decision rules, it keeps the word.
5. **Blast radius of each API rename.** For every API name, count call sites in
   `tests/`, `vignettes/`, `quarto/` --- including the six simulation cells and
   the GBSG analysis `.qmd`. This is what determines whether a rename needs a
   deprecation shim.
6. **`.fs_dg_*` prefix.** Enumerate these helpers and say what `dg` stands for at
   each. Propose a prefix consistent with the FB/MR vocabulary.

---

## 5A. Functional specification of the current "gate" (REQUIRED, do this before §6)

A rename cannot be judged without knowing what is being renamed. Produce a
complete, evidence-based description of the existing behaviour. Every answer must
cite file and line, and be **confirmed by running the code** where that is
possible --- not inferred from roxygen, which has already been found loose (the
`reselection_default` documentation at `fs_debias_gate_methods.R:156` says the
consistency engine hardcodes `"maxcons"`, while `forestsearch_main.R:2779`
derives it from `sg_focus`).

### 5A.1 What is computed

Walk `fs_debias_gate()` end to end and describe, in order: its inputs, what it
resamples, how many times, what it re-selects and by which rule, what estimates
it forms, and what it returns. Specifically:

* What quantity is de-biased, and on what scale (log-effect or ratio)?
* What interval is produced, and by which method (`ci_method = "ij"` and any
  alternatives)?
* What does the `gate` argument select among? Enumerate its permitted values and
  what each does.
* What is `t_gate`, what is its default when `NULL`, and what decision does it
  enter?
* What is `include_complement` for, and does it change the returned structure?
* How many draws by default, and what is `multiplier = "poisson"` doing?

### 5A.2 Is it always calculated?

Establish each of these by execution, not by reading defaults:

1. The default value of `debias_gate` in `forestsearch()`. Is it on or off out of
   the box?
2. Is it computed for all three engines --- consistency, DINA, GRF --- or only
   some? Note that `.fs_apply_debias_gate()` is called at
   `forestsearch_main.R:1808` and `:1940` while the consistency branch calls
   `fs_debias_gate()` directly at `:2770`. Explain why the paths differ.
3. What happens when **no subgroup is identified**? Does it run, skip, or error?
4. What happens under `sg_focus = "maxeff"`, where per-candidate consistency is
   never evaluated?
5. The consistency-branch call is wrapped in `tryCatch(...)` returning `NULL` with
   a `warning()`. **Under what conditions does that fire in practice?** Run at
   least one case where it does. A silent `NULL` that a user reads as "no harm
   confirmed" rather than "computation failed" is a reportable finding.
6. Is it computed inside the full-bootstrap replicates, inside cross-validation
   folds, or only once on the original analysis?

### 5A.3 How is it utilised

The maintainer's position is that this is **post-selection inference on a
completed analysis, not part of the identification algorithm**. Verify that the
code honours it:

1. Does anything downstream **branch** on the result --- does it alter `sg.harm`,
   `sg.harm.id`, the reported subgroup, or any estimate outside its own returned
   object? If it does anywhere, that is a defect and the headline finding of this
   survey.
2. Is the result purely reported, or does it feed a subsequent computation?
3. Does the re-selection replay the identifier faithfully? It must mirror **both**
   the selection rule (`sg_focus`) **and** the admissibility criteria in force.
   One known mismatch is already recorded: under `sg_focus = "maxeff"` the
   re-selection restricts to candidates satisfying $t_g$ --- which combines the
   effect **and** consistency floors --- while the identifier applies no
   consistency floor. Check every focus for the same class of mismatch and
   tabulate it.

### 5A.4 How is it output

1. Full field inventory of the returned object: name, type, meaning, and scale
   for every element, nested included.
2. Which of those fields does **user code** read? Check the six simulation cells
   in `quarto/simulations/gbsg_redux/`, the GBSG analysis `.qmd`, the vignettes,
   and `tests/`. The simulation recorder stores a `gate_ok` column --- trace what
   produces it. Any field user code indexes into is **API**, and its name is
   subject to the same deprecation question as a formal argument.
3. Does it appear in any `print()`, `summary()`, or plotting method? If so, the
   displayed labels are user-visible strings needing the same treatment.

### 5A.5 FB and MR: which code path is which

The maintainer is replacing "tier-1" with **full bootstrap (FB)** and "tier-2"
with **multiplier resampling (MR)**. Establish where each actually lives:

1. Which function or functions implement FB? (`nb_boots` in the simulation
   `.qmd` suggests a separate bootstrap path --- identify it.)
2. Which implement MR?
3. Does `fs_debias_gate()` perform MR only, or does it span both? **If it spans
   both, no single rename is correct and the function may need splitting** ---
   report that as a structural finding rather than proposing a name.
4. Are "tier-1"/"tier-2" used in the package at all, or only in the simulation
   `.qmd` files and the manuscripts? The scope of the rename depends on the
   answer.

Write this up as a standalone section a reader could use as a specification of
current behaviour, independent of the rename question. It is likely to be the
most durable part of this survey.

## 6. Proposed rename map

Produce a table: current name -> proposed name -> class -> call-site count ->
deprecation needed (Y/N) -> rationale.

**Propose; do not decide.** The maintainer fixes the target vocabulary. Where
several defensible options exist (e.g. `fs_mr_inference()` vs `fs_mr_debias()`,
`mr_inference` vs `run_mr` vs keeping a short form), list them with the tradeoff
rather than choosing.

For every API rename, state the deprecation mechanism you would use: accept the
old name with a `deprecated()` warning and forward it, versus a hard break with a
`NEWS.md` entry. Note that `forestsearch` 0.1.0 is on CRAN, so a hard break to
`debias_gate` invalidates existing user scripts.

Flag any rename that would collide with an existing name.

---

## 7. Deliverable

`dev/terminology-work/TERMINOLOGY_SURVEY.md`:

* guard verdict;
* full occurrence inventory, grouped by class, with counts per file;
* answers to the six questions in §5;
* the §5A functional specification, as a standalone section;
* the proposed rename map of §6;
* a **suggested execution order** --- which renames are independent, which must
  land together, and which require a deprecation shim;
* what you could not determine, and why.

## 8. Out of scope

* Any rename, in any file. Including "obviously safe" internal ones.
* Rewriting `NEWS.md` history for released versions.
* The `sg_focus` work (complete; Phases 1--5 committed) and Phase 6 (on hold).
* The `maxeff` admissibility mismatch recorded in
  `dev/sg-focus-work/article_alignment_note.qmd` --- a real defect, but a
  behaviour question rather than a naming one.
* Renaming anything in the manuscripts. FB/MR are already the manuscript
  abbreviations; this brief concerns the codebase.
