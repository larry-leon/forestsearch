# CC BRIEF — bug review of two revised `.qmd` files

```
claude "Read dev/review/CC_BRIEF_qmd_bug_review.md and execute it."
```

Targets, both in their current maintainer copies:

* `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd` (2507 lines, 69 chunks)
* `quarto/simulations/gbsg_redux/sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd` (1329 lines, 11 chunks)

**Neither has been rendered in its current state.** Both parse.

**Report only. Fix nothing, and do not edit either `.qmd`** — earlier parallel
edits forked these files and the reconciliation was expensive. Write only under
`dev/review/`.

---

## 1. Why this review exists

Both files went through several rounds of mechanical edits: a Tier→FB / gate→MR
vocabulary sweep, `sg_focus "eff"` → `"maxcons"`, symbol renames (`t1_*`→`fb_*`,
`t2_*`→`mr_*`), parameters pinned to their resolved values, a definition
reordering, an output-stem change, trimmed provenance tables, and a new
provenance block in the simulation's payload.

Four defects have already been found *after* the edits were declared safe:

| defect | how it escaped |
|---|---|
| `t2m` not renamed | `\bt2\b` cannot match inside `t2m` |
| `_t1_t2_` in the stem not renamed | `\b` does not match between `_` and `t` |
| `focus_tag` referenced `sg_focus` 60 lines before its definition | parse succeeds; only execution fails |
| `mr_in_replicates` pinned at some call sites, not all | a regex matched only the first occurrence per function name |

Every one is a class that **parsing cannot catch**. Look for more of the same
rather than reading for style.

---

## 2. Known-incomplete, confirm and quantify

`analysis_gbsg_survival_multimethod.qmd` pins `mr_in_replicates = FALSE` at
three call sites. The DINA and GRF bootstraps call the *same* function as the
main one, so they were missed. Report the exact count: how many calls to
`forestsearch_bootstrap_dofuture()`, `forestsearch_tenfold()` and
`forestsearch_Kfold()` exist, and which carry the argument. This is a no-op
either way — `FALSE` is the package default — so it is a completeness gap, not
a correctness one. Confirm that characterisation is right.

---

## 3. Checks to run — execution, not reading

### 3.1 Symbol resolution in definition order

For each file, extract the setup chunk and source it in a clean session
(stubbing `library()` and any `source()` of an external file, noting what you
stubbed). It must run to completion with no `object ... not found`. Then report
the resolved value of every knob the rest of the document consumes.

Then, for **every** chunk in order, list symbols read before they are assigned
anywhere earlier in the document. Forward references that happen to work because
a variable exists from an earlier session are the failure mode here.

### 3.2 Rename completeness

In both files, find:

* any surviving `t1`/`t2` symbol, including as a substring of a longer name
  (`t2m` is the precedent);
* any `fb_*` or `mr_*` name that is **assigned but never read**;
* any `fb_*` or `mr_*` name that is **read but never assigned**;
* the same for `gate`/`Tier` in any identifier, comment, heading, `cat()`
  string or table label. Note that `frontier`, `delegate*`, `gatekeeper`,
  `ungated` and `Negate` legitimately contain those substrings.

A read-never-assigned name yields a silent `NA` column, not an error.

### 3.3 `sprintf` and `cat` arity

Several format strings were edited for wording only. For every `sprintf()` and
`cat(sprintf())` in both files, check the number of format specifiers against
the number of arguments supplied. A mismatch is a runtime error or silent
truncation, and parsing will not catch it.

### 3.4 Field paths against real objects

Verify by introspection on a real fitted object and a real bootstrap object —
not by reading the source — that every extraction path used still resolves:

* multimethod: the payload assembly (`P$table` columns, `meta$timings$*`), the
  `fs_param_provenance()` read of `fit$args_call_all`, and the CV/LOO reads;
* simulation: `boot$H_estimates$H2` and siblings, and
  `fs.est$mr_inference$debiased$*`.

A wrong path returns `NULL`, which becomes `NA` downstream rather than erroring.

### 3.5 The new, never-executed code

Two blocks in the simulation have never run:

1. **The provenance block** in the batch save. Confirm every symbol it
   references resolves at that point (`sg_focus`, `focus_tag`,
   `consistency_method`, `stop_threshold`, `selection_rule`,
   `effect_neighborhood`, `n_workers`), that
   `paste(deparse(stop_threshold), collapse = " ")` yields `"NULL"` and not an
   error, and that `utils::packageVersion("forestsearch")` is reachable.
   Confirm the combine-path propagation using `%||%` returns `NA_character_`
   against a real pre-provenance batch `.rds`.
2. **`focus_tag`'s new branches.** It now folds `eff`/`hr`/`maxcons`,
   `effMaxSG`/`hrMaxSG`, `effMinSG`/`hrMinSG`. Exercise every legal `sg_focus`
   value through it and confirm `maxSG`/`minSG` stay distinct from
   `effMaxSG`/`effMinSG` — they are different rules, not spellings.

Also confirm the combine-mode agreement key vector is still
`c("n_sample", "nb_boots", "mr_draws", "subgroup_method", "seed_base")` with no
provenance field added. Adding one would resolve to `NA` for every existing
batch and `stop()` the pool.

### 3.6 Chunk options and eval guards

For every chunk with an `eval =` condition, confirm the guard variable is
defined before that chunk and that a `FALSE` branch leaves no downstream chunk
reading an undefined object. In the multimethod file, `run_cv` and `run_loo` are
`FALSE`, so confirm nothing after those blocks depends on objects they create.
The three bootstrap chunks carry `#| message: false`; confirm that suppresses
only console noise and no chunk relies on a message for output.

### 3.7 NULL semantics

`args_call_all[[n]] <- NULL` deletes rather than stores — this was a real
package defect. Check both files for any `[[<-` or `$<-` assignment whose
right-hand side can be `NULL`, where the intent is to store rather than remove.

---

## 4. Deliverable

`dev/review/QMD_BUG_REVIEW.md`:

1. **Verdict** — is either file unsafe to render? List anything that would
   error, and anything that would silently produce `NA` or a wrong label.
2. Per check in §3, per file: pass/fail with evidence.
3. §2 quantified.
4. Anything found this brief did not anticipate — that section is the point of
   the exercise.

Rank findings by whether they change a number, change displayed output, or
neither. Do not fix anything.
