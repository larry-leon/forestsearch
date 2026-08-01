---
title: "CC instructions — `mr_coverage_sweep_*` vocabulary + provenance pass (reorg Step 6, part 1)"
bibliography: []
---

# CC instructions — `mr_coverage_sweep_*` vocabulary + provenance pass

**Repo** `larry-leon/forestsearch` · **branch** `feature/mr-in-replicates` (base
`master`, not `main`) · **reference file**
`quarto/simulations/gbsg_redux/mr_coverage_sweep_h10_knoise0.qmd`.

## Where this file lives, and how to invoke it

Store at `dev/sim-structure/CC_SWEEP_PARITY_INSTRUCTIONS.md`, alongside
`SIM_REORG_PROPOSAL.md` and `rds_manifest.csv` — this pass implements Step 6
(part 1) of that proposal. Commit it on `feature/mr-in-replicates`.

**Single-line invocation:**

```
Read dev/sim-structure/CC_SWEEP_PARITY_INSTRUCTIONS.md and execute it end to end, exactly as written. Do not deviate from the scope constraints in it.
```

This document is designed to run **without manual intervention**. Section 1 is
an automated precondition check: if it passes, continue straight through to
Section 7; if it fails, halt and report — do not attempt to repair whatever it
flagged.

CC commits on `feature/mr-in-replicates` at the end and **does not push**.
Larry pushes.

## Hard constraint — the package codebase is read-only

**No file under `R/`, `man/`, `NAMESPACE`, `DESCRIPTION`, `tests/`, or `src/`
is to be created, edited, moved, or deleted.** The `forestsearch` R functions
are frozen for this pass. Files this pass may modify, in full:

```
quarto/simulations/gbsg_redux/mr_coverage_sweep_*.qmd
dev/sim-structure/CC_SWEEP_PARITY_INSTRUCTIONS.md   (this file, if it needs correcting)
```

Consequences that are easy to get wrong:

* `betaHhat_truth.R` is **read-only**. It is not read at all by this pass.
* **No `.rds` under `mr_sweep/<run_tag>/` is rewritten, migrated, or deleted.**
  This pass does not change the bundle schema (see Section 0, deferred item).
* Scratch files go in `/tmp`, never in the repo.
* Every command outside the sweep directory — `ls`, `git ls-files`, `git log`,
  `grep`, `diff` — is read-only by construction. A step that appears to require
  a write outside the two paths above is a wrong step: halt and report.
* Do **not** edit `.Rbuildignore`. Verify it excludes `dev/` (Section 7); if it
  does not, report it.

Scope gate, run before committing:

```bash
git status --porcelain \
  | grep -v '^.. quarto/simulations/gbsg_redux/mr_coverage_sweep_' \
  | grep -v '^.. dev/sim-structure/CC_SWEEP_PARITY_INSTRUCTIONS.md'
```

Any output means something out of scope moved. Halt and report what it was.

---

## 0. Scope, and what was verified

**In scope:** vocabulary parity with the maxcons template, and record-only
provenance in the bundle `meta`. Nothing else.

**Deliberately deferred — do not do these:**

* **The `t1_`/`t2_` → `fb_`/`mr_` column rename.** `t1_*` and `t2_*` are the
  column-name prefixes of the bundle `results` data frame (32 lines in the
  reference file), so renaming them is a **schema change to every `.rds`
  already committed under `mr_sweep/<run_tag>/`** — including the payloads the
  manuscript supplement reads live. Verified in R 4.3.3: after such a rename,
  `rH[["mr_H_lo"]]` on a legacy frame returns `NULL`, `.coverage_meta()` hits
  its `n_eff == 0L` guard and returns an all-`NA` row, which
  `drop_all_na = TRUE` discards — the cell vanishes from the grid with **no
  error and no warning**. Larry's decision: defer to the fold pass, which
  regenerates into fresh `run_tag` directories anyway. Section 3, B3 pins a
  comment on `.EST_KEYS` so nobody "tidies" it in the meantime.
* **The 8 → 1 structural fold.** Step 6 part 2. Do not add a `params:` block.

**Verified before drafting, not inferred from comments or filenames:**

| Claim | How verified |
|---|---|
| `sg_focus = "eff"` and `"maxcons"` are exact synonyms | `.normalize_sg_focus()` (`R/forestsearch_helpers.R`) maps both to `"hr"`; both were executed through `.normalize_sg_focus()` → `.fs_mr_reselection_from_focus()` for `engine = "consistency"` **and** `engine = "effect"` — identical output (`maxcons` / `maxeff` respectively). The change is inert on all three detector arms |
| Cell-to-cell delta is **2 lines** for seven cells, **3 for `h25`** | `diff` of `h10_knoise0` against `h075`, `h15`, `h20`: only `run_tag` and `target_hr_harm` differ. (`diff` prints 4 lines because it shows both sides.) `h25` additionally sets `n_sims <- 500L` against `1000L` everywhere else — a real configuration difference, encoded in its `run_tag` (`m1_h25_knoise0_s500`). `n_sims` is therefore an allowed delta in the Section 1 gate, cross-checked against the `_s<N>` token so an accidental change still fails. Section 1 re-checks all of this live, including the `knoise` cells |
| The sweep's `meta` records **no** DGM parameters | `run_cell()` writes `n_sample, n_sims, nb_boots, mr_draws, subgroup_method, sim_id_start, sim_id_end, seed_base, parallel_mode, elapsed_sec` — not `target_hr_harm`, not `k_random_noise`, not `run_tag`. Those two values (plus `n_sims` for `h25`) are the only things distinguishing the eight documents, so a bundle currently identifies its cell **only by the directory it sits in** — a filename, in the one directory known for lying filenames. Section 4 closes this, and it is a prerequisite for the fold |

---

## 1. Automated precondition gate

Write this to `/tmp/fs_preflight.R` (not the repo) and run it from the repo
root. It is read-only. **Exit status 0 ⇒ continue automatically to Section 2.
Non-zero ⇒ halt and report its output verbatim.**

```r
# /tmp/fs_preflight.R  --  read-only precondition check for the sweep parity pass
dir_sweep <- "quarto/simulations/gbsg_redux"
ref_name  <- "mr_coverage_sweep_h10_knoise0.qmd"
fail      <- character(0)

files <- sort(list.files(dir_sweep, pattern = "^mr_coverage_sweep_.*\\.qmd$",
                         full.names = TRUE))
if (!length(files)) stop("no mr_coverage_sweep_*.qmd found under ", dir_sweep)
ref <- file.path(dir_sweep, ref_name)
if (!file.exists(ref)) stop("reference file missing: ", ref)

# --- extract a top-level assignment, refusing ambiguity ----------------------
grab <- function(path, var) {
  ln <- grep(sprintf("^\\s*%s\\s*<-", var), readLines(path, warn = FALSE),
             value = TRUE)
  if (!length(ln))  return(NA_character_)
  if (length(ln) > 1L)
    stop(sprintf("%s: '%s' assigned %d times; cannot resolve",
                 basename(path), var, length(ln)))
  trimws(sub("#.*$", "", sub(sprintf("^\\s*%s\\s*<-", var), "", ln[1])))
}

tab <- do.call(rbind, lapply(files, function(f) data.frame(
  file           = basename(f),
  run_tag        = grab(f, "run_tag"),
  target_hr_harm = grab(f, "target_hr_harm"),
  k_random_noise = grab(f, "k_random_noise"),
  n_sims         = grab(f, "n_sims"),
  stringsAsFactors = FALSE)))

cat("Enumerated", nrow(tab), "sweep document(s):\n\n")
print(tab, row.names = FALSE)
cat("\n")

# --- filenames must agree with the parameters they advertise -----------------
# h075 -> 0.75, h10 -> 1.0, h15 -> 1.5 (2 digits => /10, 3 digits => /100).
for (i in seq_len(nrow(tab))) {
  cell <- sub("\\.qmd$", "", sub("^mr_coverage_sweep_", "", tab$file[i]))
  if (grepl("^h[0-9]+", cell)) {
    d <- sub("^h([0-9]+).*$", "\\1", cell)
    exp_hr <- if (nchar(d) >= 3L) as.numeric(d) / 100 else as.numeric(d) / 10
    got_hr <- suppressWarnings(as.numeric(tab$target_hr_harm[i]))
    if (!isTRUE(all.equal(exp_hr, got_hr)))
      fail <- c(fail, sprintf(
        "%s: filename implies target_hr_harm = %s but file sets %s",
        tab$file[i], format(exp_hr), tab$target_hr_harm[i]))
  }
  if (grepl("knoise[0-9]+", cell)) {
    exp_k <- as.numeric(sub(".*knoise([0-9]+).*", "\\1", cell))
    got_k <- suppressWarnings(as.numeric(tab$k_random_noise[i]))
    if (!isTRUE(all.equal(exp_k, got_k)))
      fail <- c(fail, sprintf(
        "%s: filename implies k_random_noise = %s but file sets %s",
        tab$file[i], format(exp_k), tab$k_random_noise[i]))
  }
}

# --- the only permitted cell-to-cell delta -----------------------------------
allowed <- c("run_tag", "target_hr_harm", "k_random_noise", "n_sims")

# n_sims must agree with the _s<N> token in run_tag.  n_sims is now an allowed
# cell-to-cell delta, so this is what keeps an accidental change from passing.
for (i in seq_len(nrow(tab))) {
  tg <- gsub('"', "", tab$run_tag[i])
  if (grepl("_s[0-9]+$", tg)) {
    exp_n <- as.numeric(sub(".*_s([0-9]+)$", "\\1", tg))
    got_n <- suppressWarnings(as.numeric(sub("L$", "", tab$n_sims[i])))
    if (!isTRUE(all.equal(exp_n, got_n)))
      fail <- c(fail, sprintf("%s: run_tag implies n_sims = %s but file sets %s",
                              tab$file[i], format(exp_n), tab$n_sims[i]))
  }
}

ref_lines <- readLines(ref, warn = FALSE)
for (f in setdiff(files, ref)) {
  fl <- readLines(f, warn = FALSE)
  if (length(fl) != length(ref_lines)) {
    fail <- c(fail, sprintf("%s: %d lines vs reference %d - structural drift",
                            basename(f), length(fl), length(ref_lines)))
    next
  }
  d <- which(fl != ref_lines)
  bad <- d[!vapply(d, function(i)
    any(vapply(allowed, function(v)
      grepl(sprintf("^\\s*%s\\s*<-", v), fl[i]), logical(1))), logical(1))]
  if (length(bad))
    fail <- c(fail, sprintf("%s: unexpected delta at line(s) %s: %s",
                            basename(f), paste(bad, collapse = ", "),
                            paste(trimws(fl[bad]), collapse = " | ")))
}

# --- the edits this pass will make must still be findable --------------------
must_find <- c("gate" = "gate",
               "run_tier1" = "run_tier1",
               "sg_focus_eff" = "sg_focus\\s*<-\\s*\"eff\"",
               "EST_KEYS" = "\\.EST_KEYS",
               "meta_block" = "elapsed_sec = elapsed")
for (f in files) {
  txt <- readLines(f, warn = FALSE)
  miss <- names(must_find)[!vapply(must_find, function(p)
    any(grepl(p, txt)), logical(1))]
  if (length(miss))
    fail <- c(fail, sprintf("%s: expected anchor(s) absent - file has moved on: %s",
                            basename(f), paste(miss, collapse = ", ")))
}

if (length(fail)) {
  cat("PREFLIGHT FAILED:\n"); cat(paste0("  - ", fail, collapse = "\n"), "\n")
  quit(status = 1L)
}
cat("PREFLIGHT OK - proceed to Section 2.\n")
```

```bash
Rscript /tmp/fs_preflight.R
```

Carry the printed table into the Section 7 report either way.

---

## 2. Group A — prose and comment vocabulary

Apply to every enumerated file. All five `gate` occurrences in the reference are
real; **search case-insensitively and for substrings** regardless, because a
prior residual check in this directory missed a lowercase `tier`:

```bash
cd quarto/simulations/gbsg_redux
grep -nic 'gate' mr_coverage_sweep_*.qmd      # expect 5 before, 0 after
grep -ni  'gate' mr_coverage_sweep_*.qmd \
  | grep -vi 'frontier\|delegate\|gatekeeper\|ungated\|Negate'
```

**A1** — intro prose (~L41)

```
each detected replicate records the **MR** de-biased gate
```
→
```
each detected replicate records the **MR** de-biased estimate and interval
```

**A2** — recorder prose (~L202)

```
single-cell recorder: identify Ĥ + naive + MR gate in one `forestsearch()` call,
```
→
```
single-cell recorder: identify Ĥ + naive + MR correction in one `forestsearch()` call,
```

**A3** — code comment (~L288)

```
  # ---- Identification + naive + MR gate (single forestsearch call) ----------
```
→
```
  # ---- Identification + naive + MR correction (single forestsearch call) ----
```

**A4** — plot-section prose (~L734)

```
C† (marginal subgroup HR) coverage of the **MR** gate on the **harm** block Ĥ,
```
→
```
C† (marginal subgroup HR) coverage of the **MR**-corrected interval on the **harm** block Ĥ,
```

**A5** — RESOLVED: no edit. The cross-reference at ~L72 is **correct as written**.

```
# (also sourced by the fs_t1_t2 FB engine); place it in this document's directory.
```

This step originally asserted that "the `fs_t1_t2` engine document was renamed by
the reorg" and asked for its successor. **That premise was false**, and the step
is retained only so the question is not reopened:

* The `fs_t1_t2` family was **never renamed away**. 84 `fs_t1_t2_*.qmd` are still
  tracked, roughly 40 of them in `quarto/simulations/gbsg_redux/` itself. Run
  `git ls-files '*fs_t1_t2*.qmd' | wc -l` **from the repo root** — a pathspec run
  from inside the sweep directory resolves relative to the working directory and
  returns 0, which is what made the family look deleted.
* There is **no rename path** from any `fs_t1_t2` document to any current
  `sim_fs_*.qmd`. `git log --follow` on the `sim_fs_*` files returns nothing;
  commit `98115f3` archived 12 `t1_t2` combine documents to `legacy/` and adopted
  the maxcons template as a **new** file. The relationship is adoption, not a
  rename, and it is many-to-many besides.

The line stays exactly as it is. **Do not re-litigate this in the fold pass** —
there is no stale filename here to fix.

**A6** — MR must not be described as re-selecting. MR de-biases against a fitted
candidate family; it does not identify or re-select a subgroup. After A1–A5:

```bash
grep -ni 're-select\|reselect\|re-selection' mr_coverage_sweep_*.qmd
```

**Expect 0 hits.** An earlier draft expected one survivor — the `sg_focus`
comment written in B1 — but the verbatim B1 text says "*derives MR's de-biasing
rule*" and contains no form of "re-select". Any hit at all is a prose bug:
report it with the line, do not fix it silently. Matches 5.1.

---

## 3. Group B — configuration vocabulary

**B1** — `sg_focus` spelling (~L127–128). Verified inert (Section 0).

```
consistency_method <- "resample"
sg_focus           <- "eff"     # consistency + "eff" ⇒ gate re-selection "maxcons"
```
→
```
consistency_method <- "resample"
# "maxcons" and "eff" are exact synonyms: .normalize_sg_focus() maps both to the
# canonical internal "hr".  "maxcons" is preferred because it names what the rule
# does.  forestsearch() derives MR's de-biasing rule from the normalized value and
# the engine -- "maxcons" for the consistency engine, "maxeff" for dina/grf -- so
# the spelling change is inert on all three detector arms.
sg_focus           <- "maxcons"
```

The replaced comment was also *incomplete*, not merely stale: it described only
the consistency arm, while this sweep runs `dina` and `grf` too.

**B2** — `run_tier1` → `run_fb` (~L265, ~L361). Local flag; no schema contact.

```
  run_tier1 <- is.numeric(nb_boots) && length(nb_boots) == 1L && nb_boots >= 1L
```
→
```
  run_fb <- is.numeric(nb_boots) && length(nb_boots) == 1L && nb_boots >= 1L
```

```
  if (run_tier1 && method %in% c("consistency", "dina", "grf")) {
```
→
```
  if (run_fb && method %in% c("consistency", "dina", "grf")) {
```

Residual: `grep -nic 'tier' mr_coverage_sweep_*.qmd` should leave only
`grf_selection <- "frontier"` (~L139) — a **protected word**. Do not touch it.

**B3** — pin the `.EST_KEYS` deferral (~L447). Values unchanged; comment added
so the deferred rename is not mistaken for an oversight and "tidied" into a
silent failure.

```
.EST_KEYS   <- c(oracle = "or", naive = "nv", FB = "t1", MR = "t2")
```
→
```
# DELIBERATE: the "t1"/"t2" values are the COLUMN PREFIXES of the bundle results
# frame, not vocabulary.  Every .rds already committed under mr_sweep/<run_tag>/
# carries them.  Renaming these to "fb"/"mr" without also migrating those bundles
# makes legacy cells vanish from the grid SILENTLY -- `[[` returns NULL on a
# missing column, .coverage_meta() absorbs it into an all-NA row, and
# drop_all_na discards it.  The rename is deferred to the fold pass, which
# regenerates into fresh run_tag directories.  Do not "tidy" this line.
.EST_KEYS   <- c(oracle = "or", naive = "nv", FB = "t1", MR = "t2")
```

---

## 4. Group C — record-only provenance in `run_cell()`

Replace the `meta` list in `run_cell()` (~L434–438) with:

```r
    meta = list(n_sample = n_sample, n_sims = n_sims, nb_boots = nb_boots,
                mr_draws = mr_draws, subgroup_method = method,
                sim_id_start = sim_id_start, sim_id_end = max(sim_ids),
                seed_base = seed_base, parallel_mode = parallel_mode,
                elapsed_sec = elapsed,
                # --- provenance, RECORD ONLY ------------------------------
                # Never added to the agreement check in accumulate_coverage():
                # that compares truth and seed_base only.  A missing field would
                # resolve to NA there and NA counts as a distinct value, so every
                # bundle written before this change would stop() the pool.  These
                # exist so a bundle can be identified after the fact -- which DGM
                # cell, which config, which build -- never as an agreement key.
                run_tag             = run_tag,
                dgm_model           = dgm_model,
                target_hr_harm      = target_hr_harm,
                k_random_noise      = k_random_noise,
                analysis_time       = analysis_time,
                cens_adjust         = cens_adjust,
                sg_focus            = sg_focus,
                consistency_method  = consistency_method,
                selection_rule      = selection_rule,
                forestsearch_version =
                  as.character(utils::packageVersion("forestsearch")),
                n_workers           = n_workers,
                r_version           = as.character(getRversion()),
                hostname            = unname(Sys.info()[["nodename"]]),
                built_at            = Sys.time()))
```

`target_hr_harm` and `k_random_noise` are what make a bundle self-describing
rather than directory-described, and are what the fold will key on. Every symbol
referenced is defined in the setup chunk, which precedes `run_cell()`.

**Do not** add any of these to the `require_same_truth` / `seed_base` comparison
in `accumulate_coverage()`. That check stays exactly as it is.

---

## 5. Verification

Per file, in order. Do not skip a rung.

### Hard rule for every execution-based step (5.3 and 5.4)

**Any step that *executes* document code — 5.3 as much as 5.4 — runs from a
scratch worktree, with `results_dir` redirected to a COPY of the bundles under
`/tmp`. Never point an executing step at the real `mr_sweep` tree.**

`run_sweep <- FALSE` is **not** a write guard. It gates per-cell *recomputation*
only. The `build-grid` chunk calls `save_coverage_grid()` unconditionally, so any
execution that reaches it overwrites
`mr_sweep/<run_tag>/mr_coverage_grid_<run_tag>.rds` — a **tracked** file — no
matter what `run_sweep` is set to.

This bit in 5.3, not 5.4: a chunk-slicing helper that split labels on `-` failed
to match `coverage-helpers`, the slice silently widened to the whole document,
`build-grid` ran, and the tracked grid `.rds` was overwritten. It was restored
byte-identical from `HEAD`, but a step advertised as read-only was not. Two
consequences, both mandatory:

* Redirect `results_dir` **immediately after the `setup` chunk**, before any
  later chunk can read it — `run-sweep` calls `dir.create(results_dir)` too:

  ```r
  results_dir <- file.path("/tmp/fs-sweep-scratch", run_tag)
  grid_path   <- file.path(results_dir, sprintf("mr_coverage_grid_%s.rds", run_tag))
  ```

* Assert the tree is clean after every executing step, and treat a hit as a
  failure of that step:

  ```bash
  git status --porcelain -- quarto/simulations/gbsg_redux/mr_sweep/   # expect empty
  ```

Also run executing steps from a directory where a stray `Rplots.pdf` does not
land in the repo, or delete it afterwards — base graphics opens a device.

**5.1 Residual vocabulary scan** — case-insensitive, substrings included.
`\bt2\b` does not match `t2m`, and `\b` does not match between `_` and `t`, so
`_t1_t2_` survives a word-boundary pattern. Use the unbounded forms:

```bash
grep -nic 'gate' mr_coverage_sweep_*.qmd          # expect 0
grep -ni  'tier' mr_coverage_sweep_*.qmd          # expect only grf "frontier"
grep -ni  're-select\|reselect' mr_coverage_sweep_*.qmd   # expect 0
grep -n   'sg_focus *<- *"eff"' mr_coverage_sweep_*.qmd   # expect 0
```

`t1`/`t2` will still appear — that is the deferred schema (B3) and is correct.

Both zero-expectations are exact, and each depends on the verbatim replacement
text elsewhere in this document; if either is edited, re-check the other:

* **`gate` → 0.** The Section 4 provenance comment ends "*never as an agreement
  key*" precisely so this scan stays clean. An earlier draft ended "*not to gate
  anything*", which left one hit per file and made a passing scan look like a
  failure.
* **`re-select|reselect` → 0**, not "only the B1 comment". The verbatim B1 text
  in Section 3 says "*derives MR's de-biasing rule*" and contains no such word.
  A6's reasoning is unchanged — MR de-biases, it does not re-select — the count
  is simply 0.

**5.2 Parse.** Extract and parse the R, per file:

```r
for (f in list.files(pattern = "^mr_coverage_sweep_.*\\.qmd$")) {
  p <- knitr::purl(f, output = tempfile(fileext = ".R"), quiet = TRUE)
  parse(p); cat("parse OK:", f, "\n")
}
```

**5.3 Execute the setup chunk in a clean session.** Parsing is not enough — a
symbol can be referenced before definition and still parse cleanly; it only
fails at execution. In a fresh R process, with `run_sweep <- FALSE`, source the
`setup` through `coverage-helpers` chunks of the reference file and confirm
`run_cell` and the three `*_from_bundle` functions exist and `.EST_KEYS` is
intact.

This executes document code, so the hard rule above applies in full: run it from
the worktree with `results_dir` on scratch, and confirm `mr_sweep/` is clean
afterwards. Slice the chunks by a label parse that splits on `,` and **not** on
`-` — chunk labels themselves contain `-` (`coverage-helpers`, `run-cell`,
`build-dgm`), and a `-` split silently widens the slice to the whole document:

```r
lab <- trimws(gsub("-+$", "", sub(",.*$", "", sub("^## -+", "", ln[hdr]))))
if (is.na(match("coverage-helpers", lab))) stop("chunk labels not resolved")
```

Fail loudly on an unresolved label. Do not let it fall through to a default of
"everything".

**5.4 Same-machine control render, reference file only.** `tolerance = 0` is not
a cross-machine test — worker count feeds the search batch size and changes
summation order. Same machine, same worker count, using a scratch worktree so
the live tree is never perturbed:

```bash
git worktree add /tmp/fs-control HEAD
# Render BOTH the worktree copy of mr_coverage_sweep_h10_knoise0.qmd and the
#   edited file, both with run_sweep <- FALSE, and both with results_dir
#   redirected to a scratch COPY of the bundles under /tmp -- one copy per side:
#     mkdir -p /tmp/fs-sweep-A/<run_tag> && cp <real>/*_mr_n*_res.rds /tmp/fs-sweep-A/<run_tag>/
#     (likewise /tmp/fs-sweep-B)
#   NOT at the real mr_sweep tree.  run_sweep = FALSE stops per-cell
#   recomputation but does NOT stop the grid write: build-grid calls
#   save_coverage_grid() unconditionally, so rendering against the real tree
#   overwrites the tracked mr_coverage_grid_<run_tag>.rds.  Same bundles in on
#   both sides, so the grid comparison is unaffected by the redirect.
git worktree remove /tmp/fs-control
git status --porcelain -- quarto/simulations/gbsg_redux/mr_sweep/   # expect empty
```

Diff `grid$coverage`, `grid$length`, `grid$manifest` at `tolerance = 0`,
excluding `built` and any `*_sec` column. **Expected: identical.** This pass
changes no read path, so anything other than identical means an edit landed
where it should not have — halt and report, do not rationalise a small
difference.

One control render on the reference is sufficient given the cell delta (2 lines,
3 for `h25`), provided 5.1–5.3 are clean on every file.

---

## 6. Carried forward to Step 6 part 2 (the fold)

Recorded so this pass is not redone:

* **The `t1_`/`t2_` → `fb_`/`mr_` rename**, with either a migrate-on-read shim
  or regeneration into fresh `run_tag` directories. See Section 0 for why it
  fails silently if done naively.
* **The `params:` block** is `target_hr_harm`, `k_random_noise`, `run_tag`,
  **`n_sims`** — confirmed live for all eight cells by the Section 1 gate. The
  delta is 2 lines for seven cells and 3 for `h25`, which runs `n_sims <- 500L`
  against `1000L` everywhere else. `n_sims` is **not** optional in this list:
  omit it and `h25` silently folds to 1000 replicates.
* **`run_tag` is not derivable** and must stay an explicit parameter, not a
  `sprintf()` of the other two: the reference cell's tag is
  `"m1_h10_knoise0new_s1000"`, with a `new` no formula produces.
* **Documents stay at top level.** `betaHhat_truth.R` is `source()`d from the
  working directory by ~70 documents; moving a `.qmd` into a subdirectory breaks
  that.
* **Do not revisit untracking renders** — cancelled by Larry; tracked stays
  tracked. **`git mv`, never delete.**

---

## 7. Commit and report

```bash
git status --porcelain \
  | grep -v '^.. quarto/simulations/gbsg_redux/mr_coverage_sweep_' \
  | grep -v '^.. dev/sim-structure/CC_SWEEP_PARITY_INSTRUCTIONS.md'
grep -n 'dev' .Rbuildignore        # verify only; do NOT edit
git add quarto/simulations/gbsg_redux/mr_coverage_sweep_*.qmd \
        dev/sim-structure/CC_SWEEP_PARITY_INSTRUCTIONS.md
git commit -m "sweep: FB/MR vocabulary parity and record-only provenance (Step 6 pt 1)"
```

**Do not push.**

Report these separately — do not fold them together, and do not append new
issues to an answer as though they were part of it:

1. The Section 1 preflight table, and its pass/fail.
2. The scope-gate output, plus a plain statement that nothing under `R/`,
   `man/`, `NAMESPACE`, `DESCRIPTION`, `tests/` or `src/` was modified and no
   `.rds` was rewritten or deleted.
3. The A5 outcome: resolved successor filename, or left unchanged and why.
4. Control-render result: identical, or the exact columns that moved.
5. Anything new that turned up, labelled **new finding**.
