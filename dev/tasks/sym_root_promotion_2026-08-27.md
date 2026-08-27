Save this message as dev/tasks/sym_root_promotion_2026-08-27.md and commit it, then execute. Base d884adbf. Gates are internal: stop only if one fails, otherwise finish and report once.

Rules: anchor by content, never line number. Write last. set -o pipefail, report ${PIPESTATUS[0]}. export PATH="/usr/lib/rstudio/resources/app/bin/quarto/bin:$PATH", sequential renders. Never push. Do not touch: dev/glm-continuous-sims/eval_bracket_direct.R, REPORT_blockers_1-3_ADDENDUM.md, a1_anchor_correction.py, the two frozen parent documents, the labelled provenance comment at driver site 1.

GOAL: the symmetric matrix root currently defined as .sym_root() inside the prediction document's {r anchor} chunk becomes an exported package function with a test asserting the property that motivated it. A document chunk cannot carry that test; that is the whole point of moving it.

1. NEW FILE R/fs_sym_root.R exporting fs_sym_root(S, scale = 2), body identical to the document's .sym_root: symmetrise as (S + t(S))/2, eigen(symmetric = TRUE), d <- sqrt(pmax(scale * values, 0)), return V %*% diag(d, nrow = length(d)) %*% t(V). Copy the body from the document rather than retyping it, and diff the two to prove they match before committing.

Roxygen, package conventions: markdown is ON (Roxygen: list(markdown = TRUE)) so write literal % < > & and never hand-escape Rd; @section titles plain text only. Document why the symmetric root and not V D^{1/2}: the eigenvector basis is not a continuous function of its input — signs flip and vectors rotate freely inside degenerate subspaces — so the asymmetric root is non-reproducible across arithmetically equivalent routes even though both satisfy R R' = scale * S. Include @return, @param, and a runnable @examples. Follow the file and roxygen style of R/fs_dgm_scale.R.

2. NEW FILE tests/testthat/test-fs_sym_root.R. Tests:
   - R %*% t(R) equals scale * S within tolerance, on a full-rank S and on a rank-deficient S (build one with rank 7 of 16 to match the case that motivated this).
   - the returned matrix is symmetric.
   - CONTINUITY: perturb S by 1e-12 and assert the root moves by less than 1e-6. Also assert the asymmetric root V D^{1/2} FAILS this on the same rank-deficient input — a test that only shows the fix works is weaker than one that shows the defect was real.
   - REPRODUCIBILITY: two arithmetically equivalent constructions of the same S give bit-identical roots.
   - scale is honoured (scale = 1 vs scale = 2).
   - a slightly negative eigenvalue from numerical noise is clamped to zero without error.

3. DOCUMENT: delete the .sym_root definition from {r anchor} and switch all five call sites to the package function, matching the document's existing style for calling package functions (check whether it uses library(forestsearch) with bare names or forestsearch:: prefixes, and match). Post-condition: grep -n "\.sym_root" on the .qmd returns ZERO hits — that is what rules out a missed site.

4. devtools::document(), then devtools::install() — R/ has changed, so load_all() is not sufficient. Then R CMD check on the package. Report the full result: errors, warnings, notes, each named. This is the check that was skipped on 08-26.

5. RENDER the prediction document and verify ZERO figure movement against a pre-edit baseline, tag-stripped whole-file diff, reported in full. The package function must be numerically identical to the chunk it replaces; any figure movement means it is not, and that is a stop-and-report.

Commits: one for R/ + tests, one for the document switch.

Report once: commit hashes, the body-diff proving R/ matches the document's original, the test results, the full R CMD check result, the grep proving zero .sym_root hits remain, render exit and wall-clock, the whole-file diff, git diff --stat d884adbf..HEAD, and git status --short.
