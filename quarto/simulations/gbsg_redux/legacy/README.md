# Superseded — retained as history, these do not render

The `.qmd` files here are the **superseded `t1_t2` combine generation**. They
were archived out of the parent directory by `98115f3` ("adopt kappa-carrying
maxcons template; archive 12 t1_t2 combine docs"), which adopted the revised
`sim_fs_maxcons_fb_mr_*` template carrying the kappa / beta(H) treatment. That
commit moved files only — *"No file deleted; every move is a git mv."*

The `.html` beside each `.qmd` is its output **as rendered from the parent
directory, before the move**. That is the historical record; it is what these
documents produced when they were current.

## They no longer render, and that is expected

Each sources a sibling by bare relative path:

```r
source("betaHhat_truth.R")     # line 82
```

That resolved in `quarto/simulations/gbsg_redux/`, which has the file. It does
not resolve here, because the archiving move carried the documents but not the
dependency. `betaHhat_truth.R` has **never** existed at
`quarto/simulations/gbsg_redux/legacy/` in any commit.

**Do not restore a copy here.** These were authored against a different layout
and were correct in it. Adding a dependency they never had would resurrect a
configuration that never existed, and would put a current shim inside an
archive of a superseded generation. If one of these ever needs re-running, run
it from the parent directory where it was written, or port it to the current
template.

## The lesson, stated plainly

**A `git mv` of a `.qmd` that bare-path sources a sibling orphans the
dependency. Check sourcers before archiving moves.**

A bare `source("x.R")` resolves against the render working directory, so moving
the document changes what it resolves to — silently, with no diff to the
document and nothing for a path-based grep to find. Six documents here were
orphaned this way by a move that was otherwise clean and well documented.

The same pattern has now bitten three times in this repository: these six, the
`sim-check` top-level inputs, and `dev/identifier-alignment/sim_analyses/`,
where a document sourced a file that was simply absent. Before any move or
deletion involving a simulation document, grep for bare `source()` calls in the
directory rather than for the directory's path.
