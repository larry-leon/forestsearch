# §6 check for cc_task_oc_verification_repoint_2026-08-31: every number the rendered
# document shows for the alternative and null cells must agree with the corrected
# object at the document's precision, and no superseded (08-29) value may remain.
D <- "dev/glm-continuous-sims"
C  <- readRDS(file.path(D, "oc_wrapper_grid_corrected_2026-08-30.rds"))
G  <- readRDS(file.path(D, "oc_wrapper_grid_2026-08-29.rds")); N0 <- readRDS(file.path(D, "oc_wrapper_null_2026-08-29.rds"))
html <- readLines("quarto/simulations/actg175/continuous/oc_wrapper_verification.html", warn = FALSE)
txt <- paste(html, collapse = " ")
# code-tools: true embeds the .qmd source in the page; cut it off (it carries the chunk text, not rendered numbers)
.emb <- regexpr("quarto-embedded-source", txt, fixed = TRUE); if (.emb > 0) txt <- substr(txt, 1, .emb - 1)
txt <- gsub("<script.*?</script>", " ", txt, perl = TRUE); txt <- gsub("<style.*?</style>", " ", txt, perl = TRUE)
txt <- gsub("<[^>]+>", " ", txt); txt <- gsub("&nbsp;", " ", txt); txt <- gsub("\\s+", " ", txt)
# drop the folded code blocks' source (they carry literal numbers only as code, e.g. "1601" in the stopifnot path) -- keep prose + tables
.f  <- function(x, d = 3) formatC(x, format = "f", digits = d); .fm <- function(x) formatC(x, format = "d", big.mark = ",")
# a number is a hit only as a whole token: not inside a longer number (178 in 1784) and not as a longer decimal (70.9 in 70.90)
hits <- function(s) lengths(regmatches(txt, gregexpr(paste0("(?<![0-9.,])", gsub(".", "\\.", s, fixed = TRUE), "(?![0-9.])"), txt, perl = TRUE)))
rows <- c("det_rate", "EnH", "Esens", "Espec", "Eppv", "Enpv", "EbetaH", "Enaive_bias"); d <- c(3, 1, 3, 3, 3, 3, 2, 2)
expect <- list(); old <- list()
add <- function(lst, key, val) { lst[[key]] <- val; lst }
for (n in c(500, 700)) for (g in c("resample", "split")) {
  pn <- C$alt$runs[[sprintf("n%d_%s", n, g)]]$pred; po <- G$runs[[sprintf("n%d_%s", n, g)]]$pred
  for (i in seq_along(rows)) { expect <- add(expect, sprintf("alt n%d %s %s", n, g, rows[i]), .f(pn[[rows[i]]], d[i])); old <- add(old, sprintf("alt n%d %s %s", n, g, rows[i]), .f(po[[rows[i]]], d[i])) }
  expect <- add(expect, sprintf("alt n%d %s M", n, g), .fm(pn$M)); old <- add(old, sprintf("alt n%d %s M", n, g), .fm(po$M))
}
for (n in c(500, 700)) { expect <- add(expect, sprintf("family M %d", n), .fm(C$alt$families[[as.character(n)]]$M)); old <- add(old, sprintf("family M %d", n), .fm(G$families[[as.character(n)]]$M))
  k <- C$alt$families[[as.character(n)]]$counts; ko <- G$families[[as.character(n)]]$counts
  for (f in c("enumerated", "kept", "duplicate", "size")) { expect <- add(expect, sprintf("counts %d %s", n, f), as.character(k[[f]])); old <- add(old, sprintf("counts %d %s", n, f), as.character(ko[[f]])) } }
ivn <- C$alt$invert$table; ivo <- G$invert$table
for (i in seq_len(nrow(ivn))) { key <- sprintf("alt invert n%d %s %.2f", ivn$n[i], ivn$consistency_method[i], ivn$target[i])
  j <- which(ivo$n == ivn$n[i] & ivo$consistency_method == ivn$consistency_method[i] & ivo$target == ivn$target[i])
  expect <- add(expect, key, .f(ivn$value[i], 2)); old <- add(old, key, .f(ivo$c1[j], 2)) }
nr <- C$null$runs$resample$pred; ns <- C$null$runs$split$pred; or_ <- N0$runs$resample$pred; os <- N0$runs$split$pred
nd <- c(det_rate = 4, EnH = 1, Espec = 3, Enpv = 3, EbetaH = 3, Enaive_bias = 2)
for (q in names(nd)) { expect <- add(expect, paste("null resample", q), .f(nr[[q]], nd[[q]])); old <- add(old, paste("null resample", q), .f(or_[[q]], nd[[q]]))
                       expect <- add(expect, paste("null split", q), .f(ns[[q]], nd[[q]])); old <- add(old, paste("null split", q), .f(os[[q]], nd[[q]])) }
expect <- add(expect, "null M", .fm(C$null$family$M)); old <- add(old, "null M", .fm(N0$family$M))
expect <- add(expect, "null L_eff resample", .f(C$null$runs$resample$L_eff, 1)); old <- add(old, "null L_eff resample", .f(N0$runs$resample$L_eff, 1))
expect <- add(expect, "null L_eff split", .f(C$null$runs$split$L_eff, 1)); old <- add(old, "null L_eff split", .f(N0$runs$split$L_eff, 1))
nin <- C$null$invert$table; nio <- N0$invert$table
for (i in seq_len(nrow(nin))) { key <- sprintf("null invert %s %.2f", nin$consistency_method[i], nin$target[i]); j <- which(nio$consistency_method == nin$consistency_method[i] & nio$target == nin$target[i])
  expect <- add(expect, key, .f(nin$value[i], 2)); old <- add(old, key, .f(nio$value[j], 2)) }
for (n in c(500, 700)) { m <- C$alt$measured[[as.character(n)]]; for (i in seq_along(rows)) expect <- add(expect, sprintf("measured n%d %s", n, rows[i]), .f(m[[rows[i]]], d[i])) }
res <- data.frame(key = names(expect), corrected = unlist(expect), found = sapply(unlist(expect), hits) > 0, stringsAsFactors = FALSE, row.names = NULL)
res$old <- sapply(res$key, function(k) if (is.null(old[[k]])) NA_character_ else old[[k]])
res$old_differs <- !is.na(res$old) & res$old != res$corrected & !(res$old %in% res$corrected)   # an old string equal to some corrected value elsewhere cannot be distinguished by grep
res$old_hits <- NA_integer_; res$old_hits[res$old_differs] <- sapply(res$old[res$old_differs], hits)
cat(sprintf("corrected values checked: %d ; found in render: %d ; missing: %d\n", nrow(res), sum(res$found), sum(!res$found)))
if (any(!res$found)) print(res[!res$found, c("key", "corrected")], row.names = FALSE)
cat(sprintf("superseded values that differ from corrected at the document's precision: %d ; of those still present in the render: %d\n", sum(res$old_differs), sum(res$old_hits > 0, na.rm = TRUE)))
if (any(res$old_hits > 0, na.rm = TRUE)) print(res[which(res$old_hits > 0), c("key", "old", "corrected", "old_hits")], row.names = FALSE)
cat("raw greps (whole rendered text incl. folded code):\n")
for (s in c("1601", "1,601", "1784", "1,784", "1,890", "1,696", "0.9970", "91.90", "133.11", "grid_2026-08-29", "null_2026-08-29", "corrected_2026-08-30")) cat(sprintf("  %-22s %d\n", s, hits(s)))
i <- regexpr("The family throughout", txt); cat("historical note: ", substr(txt, i, i + 520), "\n")
print(res[, c("key", "corrected", "found", "old", "old_differs", "old_hits")], row.names = FALSE)
