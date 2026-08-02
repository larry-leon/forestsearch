---
title: "CC instructions — sweep template v2: FB/MR schema rename and engine-derived run_tag"
bibliography: []
---

# CC instructions — sweep template v2

**Repo** `larry-leon/forestsearch` · **branch** `feature/mr-in-replicates` ·
**single file**
`quarto/simulations/gbsg_redux/mr_coverage_sweep_h10_knoise0.qmd`.

## Invocation

```
Read dev/sim-structure/CC_TEMPLATE_V2_INSTRUCTIONS.md and execute it end to end, exactly as written. Do not deviate from the scope constraints in it.
```

Store at `dev/sim-structure/CC_TEMPLATE_V2_INSTRUCTIONS.md`. Section 1 is an
automated gate: pass ⇒ continue to Section 2; fail ⇒ halt and report, do not
repair. Commit at the end, do not push.

## Why this pass exists

All eight sweep cells are going to be re-run. That removes the constraint that
forced the schema rename to be deferred in Step 6 pt 1: with fresh `run_tag`
directories, no committed bundle is in the blast radius, so nothing can silently
drop. This pass hardens **one** document as the template; the other seven are
disposable copies of it and are **not** touched.

## Hard constraints

**No file under `R/`, `man/`, `NAMESPACE`, `DESCRIPTION`, `tests/`, `src/` is
created, edited, moved, or deleted.** The `forestsearch` R functions are frozen.

**Exactly two files may change:**

```
quarto/simulations/gbsg_redux/mr_coverage_sweep_h10_knoise0.qmd
dev/sim-structure/CC_TEMPLATE_V2_INSTRUCTIONS.md   (this file, if it needs correcting)
```

The other seven `mr_coverage_sweep_*.qmd` are **out of scope**. Do not apply
these edits to them.

**No `.rds` is written into the repo.** Section 3's smoke run uses a `/tmp`
scratch `results_dir`. Section 5 hard rule of the Step 6 instructions applies in
full: any step that *executes* document code runs from a scratch worktree with
`results_dir` redirected before any chunk after `setup`, because `build-grid`
calls `save_coverage_grid()` unconditionally — `run_sweep = FALSE` is not a
write guard.

Scope gate before committing:

```bash
git status --porcelain \
  | grep -v '^.. quarto/simulations/gbsg_redux/mr_coverage_sweep_h10_knoise0.qmd' \
  | grep -v '^.. dev/sim-structure/CC_TEMPLATE_V2_INSTRUCTIONS.md'
```

Any output beyond the six known pre-existing untracked files means something out
of scope moved. Halt and report.

## Line numbers are stale — use anchors

The line numbers below come from the **pre-Step-6** file and are approximate:
B1 added five lines near L128 and the provenance block added roughly twenty near
L434, so everything after those points has shifted. **Every edit is specified as
an exact old-string anchor with an expected occurrence count of 1.** Use the
`/tmp/fs_apply.R` validate-all-then-write machinery from the previous pass:
transform in memory, assert every anchor matched exactly once, write nothing
unless all edits validate.

---

## 0. What was verified

| Claim | How |
|---|---|
| The sweep is **MR-only as configured** | `nb_boots <- NULL` (L86); its only call site is `run_cell(m, n, nb_boots = nb_boots, ...)` in `run-sweep` — no cell overrides it. `run_fb` is therefore `FALSE` and the FB branch never executes. All four grid consumers default to `estimators = c("oracle","naive","MR")`. The Group C smoke test confirmed `nb_boots : NULL` in written `meta` |
| The FB machinery is **present and dormant**, not absent | Ten `t1_*` columns in `.na_record`, the bootstrap branch, and the `FB` entry in `.EST_KEYS`. One line (`nb_boots <- 500L`) activates it. **Retained deliberately** — this pass renames it, does not remove it |
| `rds_tag` drives **both** the filename and the grid glob | `cell_path()` builds `sprintf("%s_%s_n%d_res.rds", method_tag, rds_tag, n)`; `build-grid` globs `sprintf("*_%s_n*_res.rds", rds_tag)`. Deriving `rds_tag` from `nb_boots` therefore cannot desynchronise the two |
| Definition order permits the derivation | `n_sims` (L85), `nb_boots` (L86), `mr_draws` (L87) all precede `run_tag` (L97) and `rds_tag` (L99) |
| `meta` already writes `mr_draws` | `mr_draws = mr_draws` in `run_cell()`. The `gate_draws` key in the 189 committed bundles is **legacy payload only** — no document edit is needed, and fresh `run_tag` directories retire it |
| Exactly **one** `t1_`/`t2_` occurrence is not schema | L72, `fs_t1_t2` — a live cross-reference to a 60-document family that was never renamed (A5). A naive `t1_` → `fb_` prefix rename corrupts it to `fs_fb_t2`. E3 guards it explicitly |

---

## 1. Automated precondition gate

Write to `/tmp/fs_v2_preflight.R`, run from the repo root. Read-only.
**Exit 0 ⇒ continue. Non-zero ⇒ halt and report verbatim.**

```r
# /tmp/fs_v2_preflight.R -- read-only precondition check for the v2 template pass
f    <- "quarto/simulations/gbsg_redux/mr_coverage_sweep_h10_knoise0.qmd"
fail <- character(0)
txt  <- readLines(f, warn = FALSE)
one  <- function(p) sum(grepl(p, txt, fixed = TRUE))

# --- Step 6 pt 1 must already be in place ------------------------------------
if (one('sg_focus           <- "maxcons"') != 1L)
  fail <- c(fail, "Step 6 pt 1 not present: sg_focus is not \"maxcons\"")
if (!any(grepl("run_fb <- is.numeric", txt, fixed = TRUE)))
  fail <- c(fail, "Step 6 pt 1 not present: run_fb absent (still run_tier1?)")
if (!any(grepl("target_hr_harm      = target_hr_harm", txt, fixed = TRUE)))
  fail <- c(fail, "Step 6 pt 1 not present: provenance meta absent")

# --- anchors this pass edits must each appear exactly once -------------------
anchors <- c(
  run_tag  = 'run_tag <- "m1_h10_knoise0new_s1000"',
  rds_tag  = 'rds_tag     <- "mr"',
  est_keys = '.EST_KEYS   <- c(oracle = "or", naive = "nv", FB = "t1", MR = "t2")',
  truthvar = "    t1 <- bundles[[1]]$truth",
  # full E6 old-string: the short form "nb_boots = nb_boots," is a SUBSTRING that
  # also matches the two call sites (forestsearch_bootstrap_dofuture(...) and
  # run_cell(...)), so it can never be unique.  E6 itself targets the full line.
  nbmeta   = "meta = list(n_sample = n_sample, n_sims = n_sims, nb_boots = nb_boots,")
for (nm in names(anchors)) {
  n <- one(anchors[[nm]])
  if (n != 1L) fail <- c(fail, sprintf("anchor '%s' matched %d times (expected 1)", nm, n))
}

# --- the fs_t1_t2 cross-reference must be exactly where E3 expects it --------
n72 <- sum(grepl("fs_t1_t2", txt, fixed = TRUE))
if (n72 != 1L)
  fail <- c(fail, sprintf("fs_t1_t2 appears %d times (expected 1); E3's guard assumes 1", n72))

# --- schema token inventory: everything else must be renameable -------------
sch <- grep("t1_|t2_", txt, value = TRUE)
sch <- sch[!grepl("fs_t1_t2", sch, fixed = TRUE)]
cat("schema lines carrying t1_/t2_ (excluding the fs_t1_t2 cross-ref):",
    length(sch), "\n")
if (length(sch) < 20L)
  fail <- c(fail, sprintf("only %d schema lines found; file has moved on", length(sch)))

# --- do the figure fragments hardcode the old run_tag? -----------------------
# WARNING, NOT A FAILURE.  run_tag becomes derived in E1, so a fragment holding
# the literal old tag stops matching.  That was originally treated as a blocker.
# It is not, and the single-path lookup was wrong twice over:
#   * the path _sim_mr_coverage.qmd does not exist -- the fragments are
#     _sim_mr_coverage_h075.qmd, _sim_mr_coverage_h10_knoise.qmd, etc.  The old
#     check reported "not found" and passed for the wrong reason.  Hence a glob.
#   * _sim_mr_coverage_h10_knoise.qmd has NO includer anywhere in the repo, and
#     the sweep template has no {{< include >}} directive, so the shared-knitr
#     -environment path is not in play; it reconstructs the grid path from its own
#     hardcoded run_tags vector rather than reading the parent's grid_path.
#   * it is a THREE-peer overlay (knoise 0/3/6), not this document's fragment.
#   * as configured it reads payload_dir = "_payloads/", not mr_sweep/ --
#     read_from_sweep is FALSE -- and _payloads/ does not exist.
#   * a missing grid takes warning(...) + return(NULL): it degrades to a skipped
#     curve, not an error.
# So it cannot block E1.  It is still worth printing: after a v2 rename the tag
# list is internally mismatched (one v2 tag, two v1 peers), which renders an
# inconsistent overlay rather than an empty one.  Fragments are OUT OF SCOPE.
frags <- sort(Sys.glob("quarto/simulations/gbsg_redux/_sim_mr_coverage*.qmd"))
if (!length(frags)) {
  cat("NOTE: no _sim_mr_coverage*.qmd fragments found -- report this\n")
} else {
  hit <- frags[vapply(frags, function(p)
    any(grepl("m1_h10_knoise0new_s1000", readLines(p, warn = FALSE), fixed = TRUE)),
    logical(1))]
  cat("figure fragments scanned:", length(frags), "\n")
  if (length(hit)) {
    cat("WARNING: fragment(s) hardcode the old run_tag (not a blocker; see above):\n")
    cat(paste0("  - ", hit, collapse = "\n"), "\n")
  } else cat("no fragment hardcodes the old run_tag: OK\n")
}

if (length(fail)) {
  cat("PREFLIGHT FAILED:\n"); cat(paste0("  - ", fail, collapse = "\n"), "\n")
  quit(status = 1L)
}
cat("PREFLIGHT OK - proceed to Section 2.\n")
```

---

## 2. Edits

### E1 — engine-derived tags

`nb_boots` becomes the single source of truth for the engine token. MR-only and
FB+MR runs land in **separate `run_tag` directories**, so they cannot pool and
`file.exists()` resume cannot cross configurations.

Old:

```
run_tag <- "m1_h10_knoise0new_s1000" # ← bump this for a fresh, isolated run
results_dir <- file.path("mr_sweep", run_tag)   # per-run folder: bundles + grid live here
rds_tag     <- "mr"                        # per-cell engine token; the method tag (fs/dina/grf) LEADS the filename
```

New:

```
# FB dormancy is the single source of truth for the engine token.  It names the
# bundle files AND the build-grid glob, which both derive from rds_tag, so the
# two cannot fall out of sync.  MR-only and FB+MR runs therefore land in separate
# run_tag directories: they never pool, and file.exists() resume never crosses a
# configuration boundary.  Bump the v2 token for a fresh isolated run.
fb_active   <- is.numeric(nb_boots) && length(nb_boots) == 1L && nb_boots >= 1L
engine_tag  <- if (fb_active) "fbmr" else "mr"
run_tag     <- sprintf("m1_h10_knoise0_%s_v2_s%d", engine_tag, n_sims)
results_dir <- file.path("mr_sweep", run_tag)   # per-run folder: bundles + grid live here
rds_tag     <- engine_tag                  # per-cell engine token; the method tag (fs/dina/grf) LEADS the filename
```

Also update the convention line in the preceding comment block if it still reads
`Convention: <dgm>_s<n_sims>_d<mr_draws>.` — the convention is now
`<dgm>_<engine_tag>_v2_s<n_sims>`.

### E2 — `run_fb` reuses `fb_active`

Old: `  run_fb <- is.numeric(nb_boots) && length(nb_boots) == 1L && nb_boots >= 1L`

New: `  run_fb <- fb_active   # single source of truth, set in setup`

### E3 — schema rename, with the `fs_t1_t2` guard

Prefix rename `t1_` → `fb_`, `t2_` → `mr_`, then two fixups. **The L72
`fs_t1_t2` cross-reference must be protected first** — it names a live
60-document family that was never renamed, and a naive prefix rename corrupts it
to `fs_fb_t2`.

```r
txt <- paste(readLines(f, warn = FALSE), collapse = "\n")

# 1. protect the cross-reference (preflight asserted exactly one occurrence)
stopifnot(length(gregexpr("fs_t1_t2", txt, fixed = TRUE)[[1]]) == 1L)
txt <- gsub("fs_t1_t2", "\001XREF\001", txt, fixed = TRUE)

# 2. prefix rename
txt <- gsub("t1_", "fb_", txt, fixed = TRUE)
txt <- gsub("t2_", "mr_", txt, fixed = TRUE)

# 3. fixups: the maxcons template names the MR timer fit_mr_secs
txt <- gsub("mr_secs", "fit_mr_secs", txt, fixed = TRUE)

# 4. restore
txt <- gsub("\001XREF\001", "fs_t1_t2", txt, fixed = TRUE)
```

Resulting names, which must match the maxcons template exactly:
`fb_secs`, `fit_mr_secs`, `fb_err`, `fb_t0`, `mr_t0`,
`fb_H_{est,lo,hi,se}`, `fb_Hc_{est,lo,hi,se}`,
`mr_H_{est,lo,hi,se_ij}`, `mr_Hc_{est,lo,hi,se_ij}`, `mr_harm_flag`.

**Assert afterwards**: `fs_t1_t2` still present exactly once; zero occurrences of
`t1_` or `t2_` **after masking the cross-reference**; `fit_fit_mr_secs` absent
(double substitution).

The masking is not optional pedantry: `fs_t1_t2` *contains* the substring `t1_`,
so a literal "zero `t1_`" assert can never pass at the same time as "cross-ref
preserved" — the two conditions are contradictory as originally written, and the
apply step fails on a correct rename. Count on the masked text:

```r
masked <- gsub("fs_t1_t2", "", txt, fixed = TRUE)
stopifnot(count_fixed(masked, "t1_") == 0L, count_fixed(masked, "t2_") == 0L)
```

(`t2_` alone happens to reach 0 unmasked, because `fs_t1_t2` ends in `t2` with no
trailing underscore. Mask both anyway — relying on that asymmetry is brittle.)

### E4 — `.EST_KEYS`, and the deferral comment it carried

The B3 comment explaining why the rename was deferred is now obsolete — remove
it with the rename. Old (the comment block plus the assignment, as B3 left it):

```
.EST_KEYS   <- c(oracle = "or", naive = "nv", FB = "t1", MR = "t2")
```

New:

```
.EST_KEYS   <- c(oracle = "or", naive = "nv", FB = "fb", MR = "mr")
```

Delete the preceding `# DELIBERATE: ...` block in full — it documents a
constraint that no longer applies. If its exact text does not match what B3
wrote, report rather than guessing at a partial deletion.

### E5 — truth variable

Old:

```
    t1 <- bundles[[1]]$truth
```

New:

```
    truth_ref <- bundles[[1]]$truth
```

and the corresponding `all.equal(t1, ...)` → `all.equal(truth_ref, ...)`.

### E6 — `nb_boots` in `meta`

`NULL` survives `saveRDS` but drops out of any `Filter`/`unlist` keying path,
and `meta` is about to become a keying structure.

Old: `meta = list(n_sample = n_sample, n_sims = n_sims, nb_boots = nb_boots,`

New: `meta = list(n_sample = n_sample, n_sims = n_sims,`
     `            nb_boots = if (is.null(nb_boots)) NA_integer_ else nb_boots,`

### E7 — loud failure when FB is requested from an MR-only bundle

Without this, requesting FB from a dormant bundle yields all-`NA` rows that
`drop_all_na` discards — silently. Add to **each** of `coverage_from_bundle()`,
`length_from_bundle()`, `bias_from_bundle()`, immediately after `res` /
`bundle$results` is resolved:

```r
  if ("FB" %in% estimators &&
      is.na(bundle$meta$nb_boots %||% NA_integer_))
    stop("FB requested but this bundle was written with nb_boots dormant ",
         "(MR-only): the fb_* columns are all NA by construction.",
         call. = FALSE)
```

---

## 3. Verification

1. **Residual scan.**
   ```bash
   sed 's/fs_t1_t2//g' <file> | grep -c 't1_\|t2_'   # expect 0 (cross-ref masked)
   grep -c 'fs_t1_t2' <file>          # expect 1
   grep -n 'fit_fit_mr_secs' <file>   # expect none
   grep -n 'run_tier1\|gate' <file>   # expect none
   ```

   The `sed` mask is required: an unmasked `grep -c 't1_\|t2_'` returns 1, not 0,
   because the preserved `fs_t1_t2` cross-reference contains `t1_`. Same reason
   as the E3 assert above.
2. **Parse.** `knitr::purl()` then `parse()`.
3. **Clean-session setup execution.** Fresh R process, source `setup` through
   `coverage-helpers`. Confirm `fb_active` is `FALSE`, `engine_tag` is `"mr"`,
   `run_tag` is `"m1_h10_knoise0_mr_v2_s1000"`, and `rds_tag` is `"mr"`.
   Slice chunk labels by splitting on `,` **not** `-`; fail loudly on an
   unresolved label.
4. **Smoke run.** Scratch worktree, `results_dir` on `/tmp` **before any chunk
   after `setup`**, `n_sims <- 2L`, `methods_run <- "consistency"`,
   `n_grid_run <- 500L`, small `mr_draws`, `run_sweep <- TRUE`. Then `readRDS`
   the bundle and confirm:
   * filename is `fs_mr_n500_res.rds` (engine token `mr`, FB dormant)
   * `results` carries `fb_*` and `mr_*` columns, **no** `t1_*`/`t2_*`
   * `mr_H_est` / `mr_H_lo` / `mr_H_hi` are populated on detected replicates
   * `fb_*` columns are all `NA` (dormant, by design)
   * `meta` has 24 keys with `nb_boots` = `NA` (not `NULL`)
5. **FB-active tag check, no compute.** In the same clean session set
   `nb_boots <- 500L` and re-evaluate only the E1 block; confirm `engine_tag`
   becomes `"fbmr"` and `run_tag` becomes `"m1_h10_knoise0_fbmr_v2_s1000"`.
   Do **not** run a cell — FB is expensive and the point is the tag.
6. `git status --porcelain -- quarto/simulations/gbsg_redux/mr_sweep/` → empty.

No control render: there is no legacy grid to be invariant against, which is the
point of the fresh tag.

---

## 4. Report

Separately, not folded together:

1. Preflight output and pass/fail.
2. Scope gate, plus a plain statement that nothing under `R/`, `man/`,
   `NAMESPACE`, `DESCRIPTION`, `tests/`, `src/` was modified, no `.rds` was
   written into the repo, and the other seven sweep documents are untouched.
3. Smoke-run bundle: filename, column names, and printed `meta`.
4. The FB-active tag check result.
5. Anything new, labelled **new finding**.
