# Concluded QMD bug review — artefacts, not render targets

A completed review of simulation `.qmd` correctness, all of it from
2026-07-31:

| file | role |
|---|---|
| `CC_BRIEF_qmd_bug_review.md` | the brief that commissioned the review |
| `CC_BRIEF_qmd_fixes.md` | the brief for applying what it found |
| `CC_BRIEF_mr_ok_and_commit.md` | the sign-off brief |
| `QMD_BUG_REVIEW.md` / `.html` | the findings |
| `QMD_FIXES_APPLIED.md` / `.html` | what was changed in response |
| `analysis_gbsg_survival_multimethod.qmd` | a **subject** of the review |
| `sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd` | a **subject** of the review |

The two `.qmd` files are copies of documents *being reviewed*. They are
evidence of what the code looked like at review time, not documents this
directory is meant to render.

## The known-broken sourcer

`sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd:85` calls

```r
source("betaHhat_truth.R")
```

which does not resolve here. `dev/review/betaHhat_truth.R` has **never** existed
in any commit — it was not deleted, it was never there. The document was copied
in for inspection, not for execution, and it renders from its own home
directory where the dependency resolves.

**Do not restore a copy here.** See
`quarto/simulations/gbsg_redux/legacy/README.md` for why bare-path `source()`
calls make copied-for-inspection documents look broken when they are not.

Nothing in this directory should be updated, re-pointed, or re-rendered. The
review is concluded and its findings are the deliverable.
