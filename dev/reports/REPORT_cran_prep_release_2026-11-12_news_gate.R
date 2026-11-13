# news_gate.R -- TASK_cran_prep_release_2026-11-12, Part 3 step 2 gate.
# Usage: Rscript --vanilla news_gate.R <NEWS.md before> <NEWS.md after>
#
# Independent of news_restructure.R (a different parser): the text above the
# "# forestsearch 0.4.0" heading is split into blank-line-separated paragraphs,
# each paragraph is split again before every line that starts with "* " (tight
# lists), and a bullet is a piece whose first character is "*".  Gate:
#   G1  the multiset of bullets before equals the multiset after, except
#       identical duplicates merged (listed);
#   G2  everything from "# forestsearch 0.4.0" to the end is byte-identical;
#   G3  above 0.4.0 the "after" file has exactly one version heading,
#       "# forestsearch 0.5.0", and the five group headings in the task's order;
#   G4  the opening summary (prose between the version heading and the first
#       group heading) is at most eight lines;
#   G5  the before file's prose that is not a bullet or a heading is listed
#       (it is replaced by the summary, not moved).
a <- commandArgs(TRUE); stopifnot(length(a) == 2)
rd <- function(f) rawToChar(readBin(f, "raw", file.size(f)))
B <- rd(a[1]); A <- rd(a[2])
Encoding(B) <- "UTF-8"; Encoding(A) <- "UTF-8"
split_at <- function(s) {
  k <- regexpr("\n# forestsearch 0.4.0\n", s, fixed = TRUE)
  stopifnot(k > 0)
  list(head = substr(s, 1, k), tail = substr(s, k + 1, nchar(s)))
}
b <- split_at(B); f <- split_at(A)
paras <- function(s) {
  p <- strsplit(s, "\n\n+")[[1]]
  p <- unlist(strsplit(p, "\n(?=\\* )", perl = TRUE))   # tight lists
  p <- sub("^\n+", "", sub("\n+$", "", p)); p[nzchar(p)]
}
pb <- paras(b$head); pa <- paras(f$head)
bul_b <- pb[startsWith(pb, "*")]; bul_a <- pa[startsWith(pa, "*")]

cat("before:", a[1], "\nafter: ", a[2], "\n\n")
cat("bullets above 0.4.0: before", length(bul_b), " after", length(bul_a), "\n")
tb <- table(bul_b); ta <- table(bul_a)
dup_b <- names(tb)[tb > 1]
only_b <- setdiff(names(tb), names(ta)); only_a <- setdiff(names(ta), names(tb))
common <- intersect(names(tb), names(ta))
count_diff <- common[tb[common] != ta[common]]
merged <- count_diff[ta[count_diff] == 1 & tb[count_diff] > 1]
other_diff <- setdiff(count_diff, merged)
cat("identical duplicates in before:", length(dup_b), "\n")
cat("merged (duplicates kept once):", length(merged), "\n")
for (m in merged) cat("  MERGED:", substr(m, 1, 90), "\n")
cat("bullets only before:", length(only_b), "\n")
for (m in only_b) cat("  ONLY BEFORE:", substr(m, 1, 90), "\n")
cat("bullets only after:", length(only_a), "\n")
for (m in only_a) cat("  ONLY AFTER:", substr(m, 1, 90), "\n")
cat("bullets with other count differences:", length(other_diff), "\n")
G1 <- !length(only_b) && !length(only_a) && !length(other_diff)
# line-level cross-check: every line of every bullet, as a multiset
lines_b <- sort(unlist(strsplit(bul_b, "\n"))); lines_a <- sort(unlist(strsplit(bul_a, "\n")))
cat("bullet lines: before", length(lines_b), " after", length(lines_a),
    " identical multisets:", identical(lines_b, lines_a), "\n")
G1 <- G1 && (identical(lines_b, lines_a) || length(merged) > 0)

G2 <- identical(charToRaw(b$tail), charToRaw(f$tail))
cat("\nG2 text from '# forestsearch 0.4.0' to the end byte-identical:", G2,
    " (", nchar(b$tail, "bytes"), "bytes )\n")

hl <- strsplit(f$head, "\n")[[1]]
vh <- grep("^# ", hl, value = TRUE); gh <- grep("^#{2,} ", hl, value = TRUE)
groups <- c("### Breaking changes", "### Bug fixes", "### New features",
            "### Superseded and deprecated", "### Documentation")
G3 <- identical(vh, "# forestsearch 0.5.0") && identical(gh, groups)
cat("\nG3 version headings above 0.4.0:", paste(vh, collapse = " | "), "\n")
cat("   group headings:", paste(gh, collapse = " | "), "\n   G3:", G3, "\n")

i1 <- which(hl == "# forestsearch 0.5.0"); i2 <- which(hl == groups[1])
summ <- hl[(i1 + 1):(i2 - 1)]; summ <- summ[nzchar(summ)]
G4 <- length(summ) >= 1 && length(summ) <= 8
cat("\nG4 opening summary lines:", length(summ), "(at most 8):", G4, "\n")
cat(paste0("   | ", summ), sep = "\n")

prose_b <- pb[!startsWith(pb, "*") & !startsWith(pb, "#")]
cat("\nG5 prose in before that is not a bullet or heading (replaced by the summary):",
    length(prose_b), "paragraph(s)\n")
for (p in prose_b) cat(paste0("   > ", strsplit(p, "\n")[[1]]), sep = "\n")

# group membership of each bullet in the after file (for the report)
cur <- NA_character_; grp <- character(0)
for (l in hl) { if (startsWith(l, "### ")) cur <- sub("^### ", "", l)
                if (startsWith(l, "* ")) grp <- c(grp, cur) }
cat("\nbullets per group (after):\n"); print(table(factor(grp, levels = sub("### ", "", groups))))

ok <- G1 && G2 && G3 && G4
cat("\nG1 multiset of bullets before == after (identical duplicates aside):", G1, "\n")
cat("NEWS GATE:", if (ok) "PASS" else "FAIL", "\n")
quit(status = if (ok) 0 else 1)
