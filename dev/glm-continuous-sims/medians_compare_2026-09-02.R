# =============================================================================
# medians_compare_2026-09-02.R -- section 5 equality gates for the medians-on-
# survivors task (dev/tasks/cc_task_medians_on_survivors_2026-09-02.md)
# -----------------------------------------------------------------------------
# Loads medians_baseline_2026-09-02.rds (Stage A, 0.3.3) and
# medians_postchange_2026-09-02.rds (Stage C, 0.3.4) and compares with
# identical() under the reused volatile-field exclusion (already applied by
# the battery's pruner).  Gates:
#   E-1  F2, F5           E-2  F1, F3
#   E-3  Eties, Enamed, Efinite, Ezero, Edegen
#   E-4  the 20-replicate gbsg bootstrap payload
# On any failure: case, component, first differing values verbatim, both
# sides, full precision.  Exit status 1 on failure.
# =============================================================================
D <- "dev/glm-continuous-sims"
a <- readRDS(file.path(D, "medians_baseline_2026-09-02.rds"))
b <- readRDS(file.path(D, "medians_postchange_2026-09-02.rds"))
cat(sprintf("baseline:   pkg %s @ %s\npostchange: pkg %s @ %s\n",
            a$pkg_version, a$git_head, b$pkg_version, b$git_head))

first_diff <- function(a, b, path = "") {
  if (identical(a, b)) return(NULL)
  if (is.list(a) && is.list(b) && !is.data.frame(a) && !is.data.frame(b) &&
      identical(names(a), names(b)) && length(a) == length(b)) {
    for (i in seq_along(a)) {
      d <- first_diff(a[[i]], b[[i]], paste0(path, "$", if (!is.null(names(a))) names(a)[i] else i))
      if (!is.null(d)) return(d)
    }
    return(list(path = paste0(path, " [attributes differ]"), a = attributes(a), b = attributes(b)))
  }
  if (is.data.frame(a) && is.data.frame(b) && identical(names(a), names(b))) {
    for (cn in names(a)) if (!identical(a[[cn]], b[[cn]])) {
      j <- which(a[[cn]] != b[[cn]] | xor(is.na(a[[cn]]), is.na(b[[cn]])))[1]
      return(list(path = paste0(path, "$", cn, if (!is.na(j)) paste0("[", j, "]")),
                  a = if (!is.na(j)) a[[cn]][j] else a[[cn]], b = if (!is.na(j)) b[[cn]][j] else b[[cn]]))
    }
    return(list(path = paste0(path, " [frame attributes differ]"), a = attributes(a), b = attributes(b)))
  }
  list(path = path, a = a, b = b)
}

GATES <- list(`E-1` = c("F2", "F5"),
              `E-2` = c("F1", "F3"),
              `E-3` = c("Eties", "Enamed", "Efinite", "Ezero", "Edegen"))
fails <- 0L
matrix_rows <- list()
for (g in names(GATES)) for (cs in GATES[[g]]) {
  ra <- a$results[[cs]]; rb <- b$results[[cs]]
  for (comp in c("pruned", "warnings", "counters", "sg.harm")) {
    ok <- identical(ra[[comp]], rb[[comp]])
    matrix_rows[[length(matrix_rows) + 1L]] <-
      data.frame(gate = g, case = cs, component = comp, identical = ok)
    if (!ok) {
      fails <- fails + 1L
      d <- first_diff(ra[[comp]], rb[[comp]], paste0(cs, "$", comp))
      cat(sprintf("\n** GATE %s FAIL: %s / %s -- first diff at %s\n", g, cs, comp, d$path))
      cat("baseline value:\n");   print(d$a, digits = 22)
      cat("postchange value:\n"); print(d$b, digits = 22)
    }
  }
}
ok_boot <- identical(a$boot$pruned, b$boot$pruned)
matrix_rows[[length(matrix_rows) + 1L]] <-
  data.frame(gate = "E-4", case = "gbsg_bootstrap_20rep", component = "pruned", identical = ok_boot)
if (!ok_boot) {
  fails <- fails + 1L
  d <- first_diff(a$boot$pruned, b$boot$pruned, "boot$pruned")
  cat(sprintf("\n** GATE E-4 FAIL: bootstrap payload -- first diff at %s\n", d$path))
  cat("baseline value:\n");   print(d$a, digits = 22)
  cat("postchange value:\n"); print(d$b, digits = 22)
}

cat("\nequality matrix:\n")
print(do.call(rbind, matrix_rows), row.names = FALSE)
cat(sprintf("\nself-consistency  A: F2 %s, Eties %s, boot %s | C: F2 %s, Eties %s, boot %s\n",
            a$selfcheck$F2, a$selfcheck$Eties, a$boot_selfcheck_ok,
            b$selfcheck$F2, b$selfcheck$Eties, b$boot_selfcheck_ok))
cat(sprintf("bootstrap per-replicate: A %.2f s -> C %.2f s (B = 1000 projection: %.1f -> %.1f min)\n",
            a$boot$per_rep, b$boot$per_rep, a$boot$per_rep * 1000 / 60, b$boot$per_rep * 1000 / 60))
cat("\nsurvfit dispatch counts (elision proof, section 6):\n")
for (cs in names(a$report)) {
  cat(sprintf("  %-8s A: n=%4d [%s] | C: n=%4d [%s]\n", cs,
              a$report[[cs]]$survfit$n, paste(names(a$report[[cs]]$survfit$by),
                                              a$report[[cs]]$survfit$by, collapse = ", "),
              b$report[[cs]]$survfit$n, paste(names(b$report[[cs]]$survfit$by),
                                              b$report[[cs]]$survfit$by, collapse = ", ")))
}
cat(sprintf("\nwalls (s): %s\n", paste(vapply(names(a$report), function(cs)
  sprintf("%s %.2f->%.2f", cs, a$report[[cs]]$wall, b$report[[cs]]$wall), character(1)), collapse = " | ")))
if (fails > 0L) { cat(sprintf("\n%d GATE FAILURE(S)\n", fails)); quit(status = 1L) }
cat("\nALL EQUALITY GATES PASS\n")
