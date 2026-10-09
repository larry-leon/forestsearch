# REPORT — CRAN release, Stage 2: remote checks (forestsearch 0.4.0)

- Date: 2026-10-09
- Task: `dev/tasks/TASK_cran_release_remote_checks_2026-10-09.md` (committed alone as `95b1153d`).
- HEAD SHA at upload time: `95b1153d951f17c6d103ac285222a076927920d4` (branch `feature/glm-extension`).
- Preconditions: P1 clean tree; P2 `Version: 0.4.0`, Stage 1 report present; P3 `pandoc`, `quarto`, `pdflatex` on `PATH`, `devtools` 2.5.2 loadable. All passed.

## Build (Step 1)

- `R CMD build` (with vignettes) in a scratch directory outside the repo: exit 0.
- Tarball: `forestsearch_0.4.0.tar.gz`, 3,346,018 bytes.
- Deviation (permitted by Step 2's rule): devtools 2.5.2 `check_win_devel()`, `check_win_release()` and `check_mac_release()` accept only a package directory and call `pkgbuild::build()` into `tempdir()`. Each call was therefore made from the repo root with `pkg = "."` and **rebuilt** the package; the Step 1 tarball itself was not uploaded. Each rebuild succeeded and `git status --porcelain` stayed clean after each call.
- No interactive prompt occurred: in non-interactive `Rscript`, `confirm_maintainer_email()` returns without prompting.

## Upload outcomes (Step 2)

- **win-builder R-devel: uploaded** (exit 0).
- **win-builder R-release: uploaded** (exit 0).
- **mac-builder R-release: failed twice, HTTP 502 Bad Gateway; no results URL.** First attempt failed; Larry approved one retry after a 10-minute wait (retry started 2026-10-09 16:08:47 PDT), which failed identically. No further retry per that disposition.
- mac-builder results URL: none (both attempts failed before a URL was returned).
- win-builder results arrive by email to the maintainer address `larry.leon.05@post.harvard.edu` (~30–90 min).

### win-builder R-devel — console output (verbatim)

```
Building windows version of forestsearch (0.4.0)
ℹ Using R-devel with win-builder.r-project.org.
── R CMD build ─────────────────────────────────────────────────────────────────
* checking for file ‘/Users/larryleon/Documents/GitHub/forestsearch/DESCRIPTION’ ... OK
* preparing ‘forestsearch’:
* checking DESCRIPTION meta-information ... OK
* installing the package to build vignettes
* creating vignettes ... OK
* checking for LF line-endings in source and make files and shell scripts
* checking for empty or unneeded directories
Removed empty directory ‘forestsearch/tests/testthat/_snaps’
* building ‘forestsearch_0.4.0.tar.gz’

────────────────────────────────────────────────────────────────────────────────
ℹ Check <'larry.leon.05@post.harvard.edu'> for the results in 15-30 mins
  (~04:21 PM).
```

### win-builder R-release — console output (verbatim)

```
Building windows version of forestsearch (0.4.0)
ℹ Using R-release with win-builder.r-project.org.
── R CMD build ─────────────────────────────────────────────────────────────────
* checking for file ‘/Users/larryleon/Documents/GitHub/forestsearch/DESCRIPTION’ ... OK
* preparing ‘forestsearch’:
* checking DESCRIPTION meta-information ... OK
* installing the package to build vignettes
* creating vignettes ... OK
* checking for LF line-endings in source and make files and shell scripts
* checking for empty or unneeded directories
Removed empty directory ‘forestsearch/tests/testthat/_snaps’
* building ‘forestsearch_0.4.0.tar.gz’

────────────────────────────────────────────────────────────────────────────────
ℹ Check <'larry.leon.05@post.harvard.edu'> for the results in 15-30 mins
  (~04:22 PM).
```

### mac-builder R-release, attempt 1 — console output (verbatim)

```
Checking macOS version of forestsearch (0.4.0)
ℹ Using https://mac.r-project.org/macbuilder/submit.html.
── R CMD build ─────────────────────────────────────────────────────────────────
* checking for file ‘/Users/larryleon/Documents/GitHub/forestsearch/DESCRIPTION’ ... OK
* preparing ‘forestsearch’:
* checking DESCRIPTION meta-information ... OK
* installing the package to build vignettes
* creating vignettes ... OK
* checking for LF line-endings in source and make files and shell scripts
* checking for empty or unneeded directories
Removed empty directory ‘forestsearch/tests/testthat/_snaps’
* building ‘forestsearch_0.4.0.tar.gz’

Error in `httr2::req_perform()`:
! HTTP 502 Bad Gateway.
Backtrace:
     ▆
  1. └─devtools::check_mac_release(pkg = ".")
  2.   └─devtools:::check_mac(...)
  3.     └─httr2::req_perform(req)
  4.       └─httr2:::handle_resp(req, resp, error_call = error_call)
  5.         └─httr2:::resp_failure_cnd(req, resp, error_call = error_call)
  6.           ├─rlang::catch_cnd(...)
  7.           │ ├─rlang::eval_bare(...)
  8.           │ ├─base::tryCatch(...)
  9.           │ │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.           │ │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.           │ │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.           │ └─base::force(expr)
 13.           └─rlang::abort(...)
Execution halted
```

### mac-builder R-release, attempt 2 (after 10-minute wait) — console output (verbatim)

```
Checking macOS version of forestsearch (0.4.0)
ℹ Using https://mac.r-project.org/macbuilder/submit.html.
── R CMD build ─────────────────────────────────────────────────────────────────
* checking for file ‘/Users/larryleon/Documents/GitHub/forestsearch/DESCRIPTION’ ... OK
* preparing ‘forestsearch’:
* checking DESCRIPTION meta-information ... OK
* installing the package to build vignettes
* creating vignettes ... OK
* checking for LF line-endings in source and make files and shell scripts
* checking for empty or unneeded directories
Removed empty directory ‘forestsearch/tests/testthat/_snaps’
* building ‘forestsearch_0.4.0.tar.gz’

Error in `httr2::req_perform()`:
! HTTP 502 Bad Gateway.
Backtrace:
     ▆
  1. └─devtools::check_mac_release(pkg = ".")
  2.   └─devtools:::check_mac(...)
  3.     └─httr2::req_perform(req)
  4.       └─httr2:::handle_resp(req, resp, error_call = error_call)
  5.         └─httr2:::resp_failure_cnd(req, resp, error_call = error_call)
  6.           ├─rlang::catch_cnd(...)
  7.           │ ├─rlang::eval_bare(...)
  8.           │ ├─base::tryCatch(...)
  9.           │ │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.           │ │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.           │ │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.           │ └─base::force(expr)
 13.           └─rlang::abort(...)
Execution halted
```

## Follow-ups (out of scope here)

- macOS check outstanding: mac-builder returned 502 on both attempts; retry later or accept win-builder + local macOS results.
- Update `cran-comments.md` once win-builder emails arrive.
- CRAN submission itself (Larry).
