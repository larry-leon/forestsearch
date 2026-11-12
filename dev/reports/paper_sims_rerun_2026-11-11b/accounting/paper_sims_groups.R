# paper_sims_groups.R -- assigns each number of accounting_numbers.csv to the task's groups
# and writes the per-block markdown tables.
#   (i)   specificity and NPV of the complement: sens(Hc), ppv(Hc), Spec, NPV
#   (ii)  detection rate and PPV: Found / any(H) / Power / Type I error / empirical detection /
#         detection counts ("39 of 44", "39/44 (89%) estimable", "2 false detections"), PPV / ppv(H)
#   (iii) everything else. Run-dependent values (timing, render date, memory, session) are flagged.
# Usage: Rscript paper_sims_groups.R <dir with accounting_numbers.csv>
a <- commandArgs(TRUE); d <- a[1]
x <- read.csv(file.path(d, "accounting_numbers.csv"), stringsAsFactors = FALSE)
for (k in c("context", "after", "row_label", "column_label")) x[[k]][is.na(x[[k]])] <- ""
cl <- trimws(x$column_label); rl <- trimws(x$row_label); cx <- x$context; af <- x$after
is_cell <- x$kind == "table-row" & x$context != "row label"
g1 <- is_cell & (cl %in% c("Spec", "NPV") | rl %in% c("sens(Hc)", "ppv(Hc)")) |
      !is_cell & grepl("(Spec|NPV)\\s*[=:]\\s*$", cx)
g2 <- is_cell & (cl %in% c("Found", "PPV") | rl %in% c("any(H)", "ppv(H)")) |
      !is_cell & (grepl("(Power|PPV|Type I Error)\\s*[=:]\\s*$", cx) |
                  grepl("^Empirical (FS|FSlg|GRF):\\s*$", cx) |
                  grepl("^ of \\d+ simulations", af) | grepl("of \\d+ simulations \\($", cx) |
                  grepl("^ false detections", af) | grepl("\\(FS:\\s*$", cx) |
                  grepl("\\(FS: \\d+/\\d+ \\($", cx) | grepl("^ estimable", af))
x$group <- ifelse(g1, "(i)", ifelse(g2, "(ii)", "(iii)"))
run_blocks <- c("(title block)", "11 Computational Timing Summary", "11.0.1 Per-Simulation Timing",
                "17 Cleanup", "18 Session Info")
x$run_dependent <- x$group == "(iii)" & (x$block %in% run_blocks |
  grepl("minutes|seconds|sec/sim|Completed in|forestsearch overall|workers|elapsed|Time \\(", paste(x$old_line, x$new_line)) |
  grepl("^(R version|Running under|Platform|BLAS|LAPACK)", x$old_line) |
  grepl("^\\d{4}-\\d{2}-\\d{2}$", x$old_line))           # the render date (date: today)
# numbers typed into the article (Appendix B: the 2024 published tables, render_reference_table();
# section 15: worked examples in prose) -- reported, but not computed by the render
x$static <- grepl("^(15|16)(\\.|\\s)", x$block)
write.csv(x, file.path(d, "accounting_grouped.csv"), row.names = FALSE)

fmt <- function(v) ifelse(is.na(v), "—", v)
diffs <- function(s) ifelse(is.na(s$old_value) | is.na(s$new_value), "—",
                            ifelse(s$changed, formatC(s$difference, format = "fg", digits = 4, flag = "+"), "0"))
lab <- function(s) {
  ifelse(s$kind == "table-row",
         sprintf("T%s %s × %s", s$table, s$row_label, ifelse(s$context == "row label", "(row label)", s$column_label)),
         sprintf("%s: …%s [#]%s…", s$kind, substr(s$context, pmax(1, nchar(s$context) - 24), nchar(s$context)),
                 substr(s$after, 1, 18)))
}
md <- c()
for (b in unique(x$block)) {
  s <- x[x$block == b, ]
  md <- c(md, "", sprintf("#### %s", b), "",
          sprintf("%d numbers: (i) %d, (ii) %d, (iii) %d; changed %d.", nrow(s), sum(s$group == "(i)"),
                  sum(s$group == "(ii)"), sum(s$group == "(iii)"), sum(s$changed)))
  show <- s[s$group %in% c("(i)", "(ii)") | s$changed, ]
  if (nrow(show)) {
    show <- show[order(match(show$group, c("(i)", "(ii)", "(iii)"))), ]
    md <- c(md, "", "| Group | Metric and cell | 2026-11-11 | 2026-11-11b | Difference |", "|---|---|---|---|---|",
            sprintf("| %s%s | %s | %s | %s | %s |", show$group, ifelse(show$run_dependent, " run", ""),
                    gsub("\\|", "\\\\|", lab(show)), fmt(show$old), fmt(show$new), diffs(show)))
  }
  nun <- sum(s$group == "(iii)" & !s$changed)
  if (nun) md <- c(md, "", sprintf("(iii) unchanged in this block: %d numbers, not listed.", nun))
}
writeLines(md, file.path(d, "accounting_tables.md"))
for (gr in c("(i)", "(ii)", "(iii)")) {
  s <- x[x$group == gr, ]; ch <- s[s$changed & !is.na(s$difference), ]
  big <- if (nrow(ch)) ch[which.max(abs(ch$difference)), ] else NULL
  cat(sprintf("%s: %d numbers (%d computed by the render, %d typed in the source), %d changed%s%s\n", gr, nrow(s), sum(!s$static), sum(s$static), sum(s$changed),
              if (!is.null(big)) sprintf(", largest change %+.4g (%s, %s -> %s)", big$difference, lab(big), big$old, big$new) else "",
              if (gr == "(iii)") sprintf("; of the changed, %d run-dependent", sum(s$changed & s$run_dependent)) else ""))
}
