# =============================================================================
# vi_grf_smoke_analysis_2026-08-30.R -- second pass over vi_grf_smoke_2026-08-30.rds
# -----------------------------------------------------------------------------
# The harness keyed the candidate family by its CUT-SET (which q-columns a row
# selects).  That key showed the families differing in most seeds even where
# every filter count and the winner were identical.  This pass rebuilds each
# candidate's SUBJECT MEMBERSHIP from its cut labels and compares the families
# as sets of subject sets -- the comparison primitive the task asks for -- and
# explains the cut-set differences (nested cuts of one variable surviving the
# order-dependent redundancy walk in one order but not the other).
# No forestsearch() calls; only the generators, to rebuild the data per seed.
# Writes dev/glm-continuous-sims/vi_grf_smoke_analysis_2026-08-30.rds
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-synthetic-dgm.R")
N <- 500L
src <- readLines("dev/glm-continuous-sims/vi_grf_smoke_2026-08-30.R")
i0 <- grep("^# ---- fixture generators", src); i1 <- grep("^# ---- the forestsearch\\(\\) call", src) - 1L
source(textConnection(paste(src[i0:i1], collapse = "\n")))          # gen_F1..gen_F8, GEN
S <- readRDS("dev/glm-continuous-sims/vi_grf_smoke_2026-08-30.rds")
res <- S$res; pairs <- S$pairs

# a cut label -> logical membership on df
lab_memb <- function(df, lb) {
  if (grepl("<=", lb, fixed = TRUE)) {
    v <- trimws(strsplit(lb, "<=", fixed = TRUE)[[1]]); df[[v[1]]] <= as.numeric(v[2])
  } else as.numeric(as.character(df[[lb]])) == 1
}
# a cut-set key "q3.1&q21.0" -> subject index set, given the label vector (q index = position)
key_memb <- function(key, df, labels) {
  parts <- strsplit(key, "&", fixed = TRUE)[[1]]
  m <- rep(TRUE, nrow(df))
  for (p in parts) {
    q <- as.integer(sub("^q([0-9]+)\\.[01]$", "\\1", p)); dir <- sub("^q[0-9]+\\.([01])$", "\\1", p)
    x <- lab_memb(df, labels[q]); m <- m & (if (dir == "1") x else !x)
  }
  paste(which(m), collapse = ",")
}
nested_same_var <- function(key, labels) {
  parts <- strsplit(key, "&", fixed = TRUE)[[1]]
  if (length(parts) != 2L) return(FALSE)
  v <- vapply(parts, function(p) { q <- as.integer(sub("^q([0-9]+)\\.[01]$", "\\1", p)); sub("\\s*<=.*$", "", labels[q]) }, "")
  v[1] == v[2]
}

out <- list()
for (i in seq_len(nrow(res))) {
  r <- res[i, ]; pk <- paste(r$fixture, r$focus, r$pcons, r$seed, sep = "_"); p <- pairs[[pk]]
  g <- GEN[[r$fixture]](r$seed); df <- g$df
  ma <- vapply(p$cand_a, key_memb, "", df = df, labels = p$conf_a)
  mb <- vapply(p$cand_b, key_memb, "", df = df, labels = p$conf_b)
  only_a <- setdiff(p$cand_a, p$cand_b); only_b <- setdiff(p$cand_b, p$cand_a)
  out[[i]] <- data.frame(
    fixture = r$fixture, focus = r$focus, pcons = r$pcons, seed = r$seed, class = r$class,
    n_rows_a = length(p$cand_a), n_rows_b = length(p$cand_b),
    cutset_symdiff = length(only_a) + length(only_b),
    cutset_symdiff_nested = sum(vapply(only_a, nested_same_var, TRUE, labels = p$conf_a)) +
                            sum(vapply(only_b, nested_same_var, TRUE, labels = p$conf_b)),
    n_memb_a = length(unique(ma)), n_memb_b = length(unique(mb)),
    memb_family_same = setequal(unique(ma), unique(mb)),
    memb_symdiff = length(setdiff(unique(ma), unique(mb))) + length(setdiff(unique(mb), unique(ma))),
    n_q_a = length(unique(unlist(lapply(strsplit(p$cand_a, "&"), function(z) sub("\\.[01]$", "", z))))),
    n_q_b = length(unique(unlist(lapply(strsplit(p$cand_b, "&"), function(z) sub("\\.[01]$", "", z))))),
    n_total_same = r$n_total_a == r$n_total_b, n_passed_same = r$n_passed_a == r$n_passed_b,
    n_eval_same = r$n_eval_a == r$n_eval_b,
    stringsAsFactors = FALSE)
}
A <- do.call(rbind, out)
cat("== candidate family: cut-set rows vs subject-membership sets, by fixture ==\n")
print(aggregate(cbind(n_rows_a, n_rows_b, cutset_symdiff, cutset_symdiff_nested, n_memb_a, n_memb_b,
                      memb_family_same, memb_symdiff, n_q_a, n_q_b, n_total_same, n_passed_same, n_eval_same) ~ fixture + focus + pcons,
                data = A, FUN = mean), digits = 3)
cat("\nfraction of cut-set differences that are nested same-variable pairs, overall (non-F4):",
    with(A[A$fixture != "F4", ], sum(cutset_symdiff_nested) / sum(cutset_symdiff)), "\n")
cat("membership families identical in", sum(A$memb_family_same[A$fixture != "F4"]), "of", sum(A$fixture != "F4"), "non-F4 pairs\n")
cat("\n== F4: what differs ==\n")
f4 <- res[res$fixture == "F4", ]
print(f4[, c("pcons", "seed", "class", "n_harm_a", "n_harm_b", "effect_a", "effect_b", "n_cand_a", "n_cand_b", "n_total_a", "n_total_b", "n_fits_a", "n_fits_b", "secs_a", "secs_b")], digits = 4, row.names = FALSE)
cat("\nF4 q-columns present in the family: -0.2", unique(A$n_q_a[A$fixture == "F4"]), " NULL", unique(A$n_q_b[A$fixture == "F4"]), "\n")
f4p <- pairs[grep("^F4_maxeffCons_0.5_", names(pairs))]
cat("F4 seed 1: sg -0.2 =", f4p[[1]]$sg_a, "| NULL =", f4p[[1]]$sg_b, "\n")
qa <- sort(unique(unlist(lapply(strsplit(f4p[[1]]$cand_a, "&"), function(z) as.integer(sub("^q([0-9]+)\\..*$", "\\1", z))))))
cat("F4 seed 1: cut labels surviving the cap at -0.2:", paste(f4p[[1]]$conf_a[qa], collapse = " | "), "\n")
cat("\n== timing: -0.2 minus NULL, total vs direct forest, by fixture (maxeffCons, pcons 0.5) ==\n")
tm <- aggregate(cbind(secs_diff, forest_secs, resid = secs_diff - forest_secs) ~ fixture + focus + pcons, data = res, FUN = mean)
tm_se <- aggregate(cbind(secs_diff) ~ fixture + focus + pcons, data = res, FUN = function(v) stats::sd(v) / sqrt(length(v)))
tm$secs_diff_se <- tm_se$secs_diff
print(tm, digits = 3)
saveRDS(list(A = A, timing = tm), "dev/glm-continuous-sims/vi_grf_smoke_analysis_2026-08-30.rds")
cat("written: dev/glm-continuous-sims/vi_grf_smoke_analysis_2026-08-30.rds\n")
