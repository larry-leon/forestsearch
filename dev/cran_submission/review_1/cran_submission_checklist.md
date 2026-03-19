# forestsearch CRAN Submission Checklist

## Current Status

- **Branch:** `cran-review-1`
- **Local check:** `devtools::check(args = "--as-cran")` — 0 errors, 0 warnings, 0 notes ✓
- **Exports:** 53 (down from 118)
- **All 8 CRAN reviewer concerns addressed**

---

## Step 1: Win-Builder Testing

Build the tarball from `cran-review-1`:

```r
devtools::build()
```

Upload `forestsearch_0.1.0.tar.gz` to **both**:

- R-release: <https://win-builder.r-project.org/upload.aspx>
- R-devel: <https://win-builder.r-project.org/upload.aspx>

Results arrive by email, typically within 30 minutes. Review for any
NOTEs, WARNINGs, or ERRORs not seen locally.

**If fixes are needed:** fix on `cran-review-1`, rebuild tarball,
resubmit to win-builder. Do not merge until win-builder passes.

---

## Step 2: Update cran-comments.md

After win-builder passes, update `cran-comments.md` to reflect the
test environments:

```markdown
## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Local: macOS Tahoe 26.3.1, R 4.5.2 (aarch64-apple-darwin20)
* win-builder: R-devel (yyyy-mm-dd rNNNNN ucrt)
* win-builder: R-release (R x.y.z)

## Downstream dependencies

This is a new submission. There are no reverse dependencies.

## Resubmission

This is a resubmission addressing all reviewer feedback from the
initial submission, including:

* Removed trailing spaces in DESCRIPTION
* Added @return documentation to all exported functions
* Removed examples from unexported/internal functions
* Replaced all \\dontrun{} with \\donttest{} or removed
* Ensured output_dir defaults to tempdir()
* Added on.exit(par(oldpar)) for all par() changes in functions
* Restored par() in vignette code chunks
* Replaced options(warn = -1) with suppressWarnings()
* Reduced exported functions from 118 to 53
```

Commit this update on `cran-review-1`.

---

## Step 3: Merge to Master

In GitHub Desktop:

1. Switch to `master` branch
2. Merge `cran-review-1` into `master`
3. Push `master`

This triggers GitHub Actions to rebuild the pkgdown site.

---

## Step 4: Verify pkgdown Site

Wait 2–3 minutes after push, then verify:

- <https://larry-leon.github.io/forestsearch/> — home page loads
- <https://larry-leon.github.io/forestsearch/reference/> — shows 53
  functions organized in 10 sections (plus Internal)
- <https://larry-leon.github.io/forestsearch/articles/forestsearch.html> —
  Getting Started vignette renders correctly

---

## Step 5: CRAN Submission

Submit at <https://cran.r-project.org/submit.html>:

1. Upload `forestsearch_0.1.0.tar.gz`
2. Paste `cran-comments.md` contents into the comments field
3. Confirm maintainer email
4. Submit

Expect a confirmation email. A CRAN volunteer typically reviews within
5–10 business days.

---

## Quick Reference: What Changed

| Area | Before | After |
|------|--------|-------|
| Exported functions | 118 | 53 |
| `\dontrun{}` blocks | 32 | 0 |
| `\donttest{}` blocks | — | 15 (all self-contained, fast) |
| Comment-only stub examples | Many | 0 |
| Legacy dead code functions | 4 | 0 |
| Missing `@return` tags | 13+ | 0 |
| `options(warn = -1)` | Present | Removed |
| Vignette `par()` restoration | `on.exit()` (broken in knitr) | Explicit `par(oldpar)` |
