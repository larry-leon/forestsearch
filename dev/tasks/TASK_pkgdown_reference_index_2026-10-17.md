# TASK: Fix the pkgdown reference index so the site build passes

- Date: 2026-10-17
- Repo: `larry-leon/forestsearch`, branch `master`. The starting HEAD must be `9b4f1c6c` or a descendant of it.
- **Why.** Every GitHub Actions pkgdown run since the 0.4.0 release (2026-10-10 onward) has failed at `build_reference_index()`:
  ```
  In _pkgdown.yml, reference[15].contents[2] (default_fs_params) must be a known topic name or alias.
  ```
  The `reference:` index in `_pkgdown.yml` still reflects the March site. It names topics that no longer exist and omits most of the topics that 0.4.0 added. The last successful site build was on 2026-03-23.
- **Change.** `_pkgdown.yml` `reference:` section only. **No `R/`, `man/`, `NAMESPACE` or `DESCRIPTION` changes.**
- **Report.** `~/Downloads/REPORT_pkgdown_reference_index_2026-10-17.md`

## Standing rules

- Do not run `git fetch`, `git pull` or `git push`. Run each git command separately. Never stage pre-existing untracked files.
- Do NOT run `R CMD check`, testthat, `pkgdown::build_site()`, `build_articles()`, `build_reference()` with examples, or `clean_*()`. Do not install packages.
- Gates are stop-on-failure, not stop-to-ask.
- Hard cap: 20 minutes.

## Step 0 — Record the task

1. Copy this file from `~/Downloads/` to `dev/tasks/TASK_pkgdown_reference_index_2026-10-17.md`.
2. `git add` that path.
3. Commit with the message `dev/tasks: add pkgdown reference index task (2026-10-17)`.

## Step 1 — Diagnose with pkgdown itself (report section 1)

Using the installed pkgdown 2.2.1, run `pkgdown::check_pkgdown(".")` and record its full output.

Then compute two lists from `man/*.Rd` (`\name` and `\alias`, and `\keyword{internal}`) against the `contents:` entries of every `reference:` section. Treat selector entries (`starts_with(`, `matches(`, `has_keyword(`, `has_concept(`, `-...`) the way pkgdown does; if in doubt, call pkgdown's own internal selector code.

- **A. Listed but unknown.** Entries that match no topic name or alias. Expected: `default_fs_params`, `default_grf_params`.
- **B. Missing from the index.** Topics that are not `\keyword{internal}` and are not matched by any entry, whether by name, by alias, or by a selector.

Record both lists. For each name in A, also record:
- whether it exists in the installed CRAN 0.4.0 namespace (`exists(x, envir = asNamespace("forestsearch"))`)
- `git log -S'<name>' --oneline -- R | head -3`

## Step 2 — Edit `_pkgdown.yml` (report section 2)

1. **Remove** every entry in list A from its section. Leave no empty `contents:` behind.
2. **Add** every topic in list B to a new section appended as the **last** `reference:` section:
   ```yaml
     - title: "Additional functions (0.4.0)"
       desc: >
         Functions added in forestsearch 0.4.0, pending placement in the
         sections above.
       contents:
         - <topic 1>
         - <topic 2>
   ```
   List the topics alphabetically, one per line. Change no other line of `_pkgdown.yml`.
3. Gate: `git diff --stat` touches only `_pkgdown.yml`.
4. Gate: `yaml::read_yaml("_pkgdown.yml")` parses.

## Step 3 — Verify (report section 3)

1. Gate: `pkgdown::check_pkgdown(".")` reports no errors. Record its output.
2. Gate: recomputing lists A and B from Step 1 gives empty lists.
3. Gate: `pkgdown::build_reference_index(".")` completes without error. This writes only `docs/reference/index.html` and runs no examples; leave the output in place.

## Step 4 — Commit (report section 4)

1. `git add _pkgdown.yml`
2. Commit with the message `pkgdown: fix reference index for 0.4.0 (drop removed topics, add new ones)`.
3. Gate: `git show --stat HEAD` lists only `_pkgdown.yml`.

## Step 5 — Record only (report section 5)

These are recorded for Larry's decision. Do not act on either.

- **CLAUDE.md on the public site.** The CI log shows `Reading CLAUDE.md` → `Writing CLAUDE.html`, so `CLAUDE.md` is published on the public site. Record:
  - `git ls-files CLAUDE.md .claude`
  - whether `.claude/CLAUDE.md` exists
  - the first 10 lines of `CLAUDE.md`
- **Other top-level `.md` files.** List every other top-level `.md` file that pkgdown would render as a page, as shown by `pkgdown:::package_mds(".")` if that function exists.

## Report

Write `~/Downloads/REPORT_pkgdown_reference_index_2026-10-17.md` with sections 1–5. In the terminal, print one line only: the report path and overall PASS/FAIL.
