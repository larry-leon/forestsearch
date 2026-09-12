#!/usr/bin/env python3
"""Transplant summary_dinamr.qmd -> summary_grfmr.qmd.

TASK_grfmr_campaign_2026-09-11, STAGE 3.  Globs and labels only, plus the
GRF-specific substitutions the task names:

  * admitted_n REPLACES n_family as the STRATIFIER in the strata section, with
    the reason stated in every caption it touches -- GRF's n_family is the
    outcome-independent enumerated pool.  n_family is KEPT in the descriptive
    tables, and an admitted_n descriptive table is added beside it.
  * the HR 1.00 null cells are removed from the cell list: they are not in this
    task and go to a follow-up session, so they are absent by design rather
    than deferred.
  * the framing is rewritten: this does not certify GRF, GRF is NOT
    FS-analogous, and every coverage column is the conditional-on-proposed-
    family estimand.

Everything else -- the p-hat strata and the joint count table, the error-SD
Gaussian reference beside the marginal one with both formulas, Wilson intervals
on every rate, marginal SD and error SD side by side, absolute coverage levels
for every product, and the FS comparator with its sg_focus / eps and the
confound sentence -- is carried through UNCHANGED.

Run from quarto/simulations/gbsg_020.
"""
import re, sys

SRC, DST = "summary_dinamr.qmd", "summary_grfmr.qmd"
s = open(SRC).read()
n_edits = 0
def sub1(old, new, why):
    global s, n_edits
    if s.count(old) != 1:
        sys.exit("EXPECTED EXACTLY ONE OCCURRENCE (%d found) for: %s" % (s.count(old), why))
    s = s.replace(old, new); n_edits += 1

# ---------------------------------------------------------------- 1. header
sub1('title: "DINA at the FS-analogous criterion: naive / field / field-s / IJ two-term across both prevalences"',
     'title: "GRF: naive / field / field-s / IJ two-term across both prevalences"',
     "title")
old_sub = s[s.index('subtitle: "'):s.index('\nauthor:')]
s = s.replace(old_sub,
  'subtitle: "Campaign grfmr (GRF frontier, effMaxSG eps 0.20, n 500 / 1000 / 1500 at HR 1.50 / 1.75, '
  '12.4% and 31% prevalence); TASK_grfmr_campaign_2026-09-11; transplanted from summary_dinamr.qmd"')
n_edits += 1

# framing block: everything from the first blockquote to the '# Bundles' header
i0 = s.index("> **This does not certify DINA.**")
i1 = s.index("# Bundles and the Gate 2 record")
s = s[:i0] + """> **This does not certify GRF.** **Every coverage number in this document is coverage of
> β(Ĥ) conditional on the proposed family, over selected replicates**, and every table says so
> in its caption. Nothing here is scored against an acceptance criterion; none was
> pre-registered and none is proposed. The FS numbers quoted throughout are read as they stand,
> as reference lines, not as a bar.

> **GRF is not "FS-analogous".** A GRF-to-FS or GRF-to-DINA comparison differs in **identifier**,
> in **family construction**, in **detection set**, and — **at the DR pre-filter only, not at
> admission** — in the **scale of the selection criterion**. `dmin.grf = 0.0` is a floor on the
> doubly-robust score (`grf_main.R:291`), which targets RMST for survival outcomes and is not
> alignable with the Cox log-HR floor that FS and DINA share. The floor that binds on the
> re-selection path is `hr.threshold = 0.90`, the same one DINA carries. Both statements belong
> beside every comparison in this document.

**The engine.** `subgroup_method = "grf"`, `grf_selection = "frontier"`,
`grf_select_statistic = "effect"`, `grf_depth = 2`, `dmin.grf = 0.0`, `sg_focus = "effMaxSG"` at
ε = 0.20. The frontier applies the shared inclusion band as a **filter with no empty-band
fallback**, where DINA uses it as a sort key. `FS_S7_ER_JCUTS` is **inert on GRF** — it feeds the
consistency-only `method_args` — so it was deliberately not set.

**`n_family` on GRF is the enumerated pool, not the qualified set.** `.grf_dr_candidates()`
(`R/grf_subgroup_labels.R:255–277`) enumerates from quantiles of X subject to `n_min`, so the pool
**does not depend on the outcome**: the cost probes measured it identical across the two
prevalences replicate by replicate. Stratifying on it would stratify on something
outcome-independent. **`admitted_n` — the count of enumerated candidates whose inferential effect
clears the resolved admission floor (`R/forestsearch_helpers.R:1654`) — is therefore the
stratifier used in @sec-strat**, and it is the analogue of DINA's family-size stratifier.
`n_family` is kept in the descriptive tables of @sec-family, beside it.

**Grid.** Twelve harm cells at ε = 0.20: **Block A** = 12.4% (`harm_z1_quantile` at the template
default 0.25) and **Block B** = 31% (`FS_S7_Z1Q = 0.60`), each at HR 1.50 and 1.75 for
n = 500 / 1000 / 1500. Seeds 8316951 + sim_id, sim_id 1–2000 as two seed-disjoint batches of
1,000, then combined. **The HR 1.00 null cells are not in this task** — they go to a follow-up
session, so they are absent by design and are not listed as deferred.

**Standing conventions.** Bounds read by location, never as significance at the null; Wilson
95% limits everywhere; marginal and error SDs side by side; NPV beside sens/spec/PPV;
winner-only and winner-floor excluded throughout.

""" + s[i1:]
n_edits += 1

# ---------------------------------------------------------------- 2. globs
sub1('''       file = sprintf("%sdina_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d%s_nb20_dinamr_combined_1_2000.rds",''',
     '''       file = sprintf("%sgrf_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d%s_nb20_grfmr_combined_1_2000.rds",''',
     "bundle glob")

# ---------------------------------------------------------------- 3. cells: harm only
sub1("""  lapply(list(c(1.00,500),c(1.00,1000),c(1.00,1500)),
         function(x) mk("A", 0.124, as.integer(x[2]), x[1], "null")),
  lapply(list(c(1.00,500),c(1.00,1000),c(1.00,1500)),
         function(x) mk("B", 0.31,  as.integer(x[2]), x[1], "null")))""",
"""  # The HR 1.00 null cells are NOT in this task (kickoff: "The HR 1.00 null
  # cells are NOT in this task").  They are removed from the cell list rather
  # than left to the absent-cell guard, so that the render reports them as out
  # of scope instead of as deferred.
  NULL)""",
     "drop the null cells")
sub1("cells <- c(\n", "cells <- Filter(Negate(is.null), c(\n", "Filter opener")
sub1("names(cells) <- vapply(cells, `[[`, \"\", \"label\")",
     "))\nnames(cells) <- vapply(cells, `[[`, \"\", \"label\")", "Filter closer")

sub1('cat(sprintf("Cells on disk: %d of %d.\\n", sum(present), length(cells)))',
     'cat(sprintf("Cells on disk: %d of %d harm cells (the HR 1.00 null cells are out of scope for this task).\\n", sum(present), length(cells)))',
     "cells-on-disk line")

# ---------------------------------------------------------------- 4. stratifier
sub1("""strata_of <- function(d, by) {
  v <- if (by == "K") d$n_family else d$p_hat_H""",
"""# THE GRF SUBSTITUTION.  by == "K" now reads admitted_n, NOT n_family.
# n_family on GRF is .grf_dr_candidates()'s enumeration from quantiles of X
# subject to n_min (R/grf_subgroup_labels.R:255-277), so it does not depend on
# the outcome and is near-constant across replicates AND across prevalences --
# stratifying on it would stratify on nothing the outcome touches.  admitted_n
# is the count that clears the resolved admission floor
# (R/forestsearch_helpers.R:1654) and is the analogue of DINA's family-size
# stratifier.  The p-hat stratification is UNCHANGED.
strata_of <- function(d, by) {
  v <- if (by == "K") d$admitted_n else d$p_hat_H""",
     "strata_of stratifier")
sub1("""                    format(min(v[!is.na(g) & g==t])), format(max(v[!is.na(g) & g==t])))) }
  if (by == "K") {
    out[[length(out)+1]] <- list(d = d[d$n_family == 1L, , drop = FALSE], ov = TRUE, lab = "K = 1 (overlapping)")
    out[[length(out)+1]] <- list(d = d[d$n_family <= 5L, , drop = FALSE], ov = TRUE, lab = "K <= 5 (overlapping)") }""",
"""                    format(min(v[!is.na(g) & g==t])), format(max(v[!is.na(g) & g==t])))) }
  if (by == "K") {
    out[[length(out)+1]] <- list(d = d[d$admitted_n == 1L, , drop = FALSE], ov = TRUE, lab = "admitted_n = 1 (overlapping)")
    out[[length(out)+1]] <- list(d = d[d$admitted_n <= 5L, , drop = FALSE], ov = TRUE, lab = "admitted_n <= 5 (overlapping)") }""",
     "overlapping rows")
sub1("""                    format(min(v[!is.na(g) & g==t])), format(max(v[!is.na(g) & g==t])))) }""",
"""                    format(min(v[!is.na(g) & g==t])), format(max(v[!is.na(g) & g==t])))) }""",
     "noop anchor") if False else None
s = s.replace('''lab = sprintf("%s T%d [%s, %s]", if (by == "K") "K" else "p-hat", t,''',
              '''lab = sprintf("%s T%d [%s, %s]", if (by == "K") "adm" else "p-hat", t,''')
n_edits += 1

# det_rows must require admitted_n finite as well
sub1("""  keep <- is.finite(d$betaHhat_H) & is.finite(d$betaHhat_Hc) & is.finite(d$n_family) &
          is.finite(d$p_hat_H) & is.finite(d$fld_H_lo1s) & is.finite(d$fld_Hc_up1s_s)""",
"""  keep <- is.finite(d$betaHhat_H) & is.finite(d$betaHhat_Hc) & is.finite(d$n_family) &
          is.finite(d$admitted_n) &
          is.finite(d$p_hat_H) & is.finite(d$fld_H_lo1s) & is.finite(d$fld_Hc_up1s_s)""",
     "det_rows requires admitted_n")

# the strata prose
sub1("""- **Proposed-family size.** Within-cell empirical tertiles of `n_family`, plus **two extra rows that
  deliberately overlap them**: `K = 1` (the family is a single candidate, so the "family" imposes no
  choice at all) and `K <= 5`. Overlapping rows are marked; they are not a partition and must not be
  summed.""",
"""- **Admitted-candidate count.** Within-cell empirical tertiles of **`admitted_n`**, plus **two extra
  rows that deliberately overlap them**: `admitted_n = 1` (one candidate clears the floor, so the
  admitted set imposes no choice at all) and `admitted_n <= 5`. Overlapping rows are marked; they are
  not a partition and must not be summed.
  **`admitted_n` and not `n_family`, deliberately.** `n_family` on GRF is the enumerated pool —
  `.grf_dr_candidates()` enumerates from quantiles of X subject to `n_min`
  (`R/grf_subgroup_labels.R:255–277`), so it does not depend on the outcome and is near-constant
  within a cell and identical across prevalences replicate by replicate. A stratification on it
  would not be a stratification on anything the outcome touches. `admitted_n` is the count whose
  inferential effect clears the resolved admission floor (`R/forestsearch_helpers.R:1654`) and is
  the analogue of the family-size stratifier the DINA summary uses. `n_family` is kept in the
  descriptive tables of @sec-family.""",
     "strata prose")

# ---------------------------------------------------------------- 5. captions
s = s.replace("EVERY PRODUCT by within-cell proposed-family-size (n_family) stratum",
              "EVERY PRODUCT by within-cell ADMITTED-CANDIDATE-COUNT (admitted_n) stratum -- admitted_n and NOT n_family, because n_family on GRF is the outcome-independent enumerated pool (.grf_dr_candidates() enumerates from quantiles of X subject to n_min) while admitted_n is the count clearing the resolved admission floor")
s = s.replace("K tertiles are a partition; the K = 1, K <= 5 and 'all detected' rows OVERLAP",
              "admitted_n tertiles are a partition; the admitted_n = 1, admitted_n <= 5 and 'all detected' rows OVERLAP")
n_edits += 2

# section titles that name the stratifier
s = s.replace("## Every product by proposed-family-size stratum {#sec-strat-k}",
              "## Every product by admitted-candidate-count stratum {#sec-strat-k}")
s = s.replace("## Family-size stratum by p&#770; bin: counts {#sec-strat-xtab}",
              "## Admitted-candidate-count stratum by p&#770; bin: counts {#sec-strat-xtab}")
s = s.replace("# Stratified tables: by proposed-family size and by p&#770; {#sec-strat}",
              "# Stratified tables: by admitted-candidate count and by p&#770; {#sec-strat}")
s = s.replace('skip_note("the family-size stratified tables")',
              'skip_note("the admitted_n stratified tables")')
n_edits += 4

# remaining n_family mentions inside the strata section's helpers/captions
s = s.replace('paste0("K T", tert_fs(d$n_family))', 'paste0("adm T", tert_fs(d$admitted_n))')
s = s.replace('K=paste0("K T",tert_fs(d$n_family))', 'K=paste0("adm T",tert_fs(d$admitted_n))')
n_edits += 1

# ---------------------------------------------------------------- 6. labels
s = s.replace("DINA", "GRF").replace("dinamr", "grfmr").replace("dina_", "grf_")
s = s.replace("summary_cert20.qmd", "summary_cert20.qmd")   # unchanged, kept explicit
n_edits += 1

open(DST, "w").write(s)
print("wrote %s (%d guarded edits + label sweep)" % (DST, n_edits))
