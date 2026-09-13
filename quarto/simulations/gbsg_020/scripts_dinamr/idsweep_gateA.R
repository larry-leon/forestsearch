# Gate A -- alignment, PER RUN, stop-on-failure (TASK_idsweep_2026-09-12).
#
# Resolves, for one render, every value the gate names and prints it:
#   MR off, with `_nomr` in the stem and the bundle filename;
#   sg_focus as intended (and the stem tag fs_focus_tag() gives it);
#   eps 0.20 on effMaxSG / effMinSG, and FS_S7_NBHD UNSET otherwise (the
#     template then forwards its 0.10 default, which those foci never read);
#   on GRF: grf_selection = "frontier", frontier_rule as intended, eps 0.20.
#
# Each value is read from as many independent places as exist:
#   0. the environment the driver passed to render.sh (set and unset lists);
#   1. the render log (RC = 0);
#   2. the rendered document -- the template's own audit lines "Output stem:",
#      "Run config:" and "Template knobs:" (the echoed format strings, which
#      contain %s, are skipped);
#   3. the bundle's meta;
#   4. GRF: frontier_rule is not recorded anywhere in the bundle.  It is the
#      switch at R/forestsearch_main.R (sg_focus -> frontier_rule), so it is
#      resolved by evaluating THAT switch -- from the installed package's
#      forestsearch() body AND from the repository source -- on the normalized
#      focus, and the two must agree.  admitted_n finite on the bundle is the
#      frontier/effect path's own output (see gate3.R).
#
# usage: Rscript idsweep_gateA.R <out-basename> <engine> <sg_focus> <z1q|-> <n> <hr> "<env set>" "<env unset>"
#   IDSWEEP_TAG (default idsweep) and IDSWEEP_REPS (default 500) exist only so
#   the gate can be exercised on the committed pBoc smoke before the sweep.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 6L)
OUT <- args[1]; E <- args[2]; FO <- args[3]; Z <- args[4]
N <- as.integer(args[5]); H <- as.numeric(args[6])
KN <- if (length(args) >= 7L) strsplit(trimws(args[7]), " +")[[1]] else character(0)
UN <- if (length(args) >= 8L) strsplit(trimws(args[8]), " +")[[1]] else character(0)
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
REPO    <- normalizePath(file.path(QMD_DIR, "..", "..", ".."))
TAG     <- Sys.getenv("IDSWEEP_TAG", unset = "idsweep")
REPS    <- as.integer(Sys.getenv("IDSWEEP_REPS", unset = "500"))
WORKERS <- 12L
band    <- FO %in% c("effMaxSG", "effMinSG")
z1q     <- if (Z == "-") 0.25 else as.numeric(Z)
EPS     <- if (band) 0.20 else 0.10

fail <- 0L; pass <- 0L
say <- function(ok, msg) { ok <- isTRUE(ok)
  if (ok) pass <<- pass + 1L else fail <<- fail + 1L
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", msg)) }
kv <- function(k, v) cat(sprintf("  %-26s %s\n", k, paste(format(v), collapse = " ")))

mt   <- if (E == "consistency") "fs" else E
ftag <- forestsearch::fs_focus_tag(E, FO)
stem <- sprintf("%s_%s_fb_mr_field_m1_h%03d_knoise0_n%d%s%s_nomr_%s", mt, ftag, round(100 * H), N,
                if (abs(z1q - 0.25) < 1e-12) "" else sprintf("_z1q%02d", round(100 * z1q)),
                if (band) "_nb20" else "", TAG)
bun  <- file.path(QMD_DIR, "results", sprintf("%s_res_1_%d.rds", stem, REPS))
html <- file.path(QMD_DIR, paste0(OUT, ".html"))
lg   <- file.path(SCRATCH, "logs", paste0(OUT, ".log"))
cat(sprintf("=== GATE A: %s  (engine %s, sg_focus %s, z1q %s, n %d, HR %.2f) ===\n", OUT, E, FO, Z, N, H))
kv("expected stem", stem)

## 0. the environment the driver passed
cat("\n--- 0. driver environment ---\n")
kv("set", KN); kv("unset", UN)
has  <- function(x) x %in% KN
unst <- function(v) v %in% UN
say(has("FS_S7_MR=FALSE"), "FS_S7_MR=FALSE passed")
say(has(paste0("FS_S7_FOCUS=", FO)), sprintf("FS_S7_FOCUS=%s passed", FO))
if (band) {
  say(has("FS_S7_NBHD=0.20") && !unst("FS_S7_NBHD"), "band focus: FS_S7_NBHD=0.20 passed")
} else {
  say(unst("FS_S7_NBHD") && !any(grepl("^FS_S7_NBHD=", KN)), "non-band focus: FS_S7_NBHD unset")
}
if (E == "consistency") say(unst("FS_S7_METHOD") && !any(grepl("^FS_S7_METHOD=", KN)), "consistency: FS_S7_METHOD unset")
if (E != "consistency") say(has(paste0("FS_S7_METHOD=", E)), sprintf("FS_S7_METHOD=%s passed", E))
if (Z == "-") say(unst("FS_S7_Z1Q") && !any(grepl("^FS_S7_Z1Q=", KN)), "12.4%: FS_S7_Z1Q unset")
if (Z != "-") say(has(paste0("FS_S7_Z1Q=", Z)), sprintf("31%%: FS_S7_Z1Q=%s passed", Z))
say(unst("FS_S7_ER_JCUTS") && !any(grepl("^FS_S7_ER_JCUTS=", KN)), "FS_S7_ER_JCUTS unset")

## 1. render log
cat("\n--- 1. render log ---\n")
wl <- if (file.exists(lg)) grep("^WALL_SECONDS=", readLines(lg, warn = FALSE), value = TRUE) else character(0)
kv("log", if (length(wl)) tail(wl, 1) else "<no WALL_SECONDS line>")
say(length(wl) && grepl(" RC=0 ", tail(wl, 1)), "render RC = 0")

## 2. rendered document
cat("\n--- 2. rendered document (the template's audit lines) ---\n")
ht <- if (file.exists(html)) paste(readLines(html, warn = FALSE), collapse = "\n") else ""
say(nzchar(ht), sprintf("document present: %s", basename(html)))
grab <- function(pat) { m <- regmatches(ht, gregexpr(pat, ht))[[1]]; m <- m[!grepl("%s", m, fixed = TRUE)]
  if (length(m)) m[1] else NA_character_ }
st <- grab("Output stem: [^<\n]*"); rc <- grab("Run config: [^<\n]*"); tk <- grab("Template knobs: [^<\n]*")
kv("Output stem", st); kv("Run config", rc); kv("Template knobs", tk)
say(identical(sub("^Output stem: ", "", st), stem), "document's stem == expected stem")
say(grepl("_nomr_", st, fixed = TRUE), "`_nomr` in the document's stem")
tok <- if (is.na(tk)) character(0) else strsplit(sub("^Template knobs: ", "", tk), " ")[[1]]
tok <- tok[grepl("=", tok)]
K <- setNames(sub("^[^=]*=", "", tok), sub("=.*$", "", tok))
want <- c(method = E, hr = sprintf("%.2f", H), n = as.character(N), z1q = sprintf("%.2f", z1q),
          focus = FO, nbhd = sprintf("%.2f", EPS), er_jcuts = "10", fb_mode = "none",
          campaign = TAG, run_mode = "batch", quickrun = "FALSE", mr_inference = "FALSE")
for (k in names(want)) say(identical(unname(K[k]), unname(want[k])),
                           sprintf("audit %s = %s (intended %s)", k, K[k], want[k]))
say(grepl(sprintf("sim_id 1-%d", REPS), tk, fixed = TRUE), sprintf("audit sim_id 1-%d", REPS))
cat("  inert with MR off, reported only: ")
cat(paste(sprintf("%s=%s", names(K), K)[names(K) %in% c("ci_method", "uniform", "field_complement",
      "field_decompose", "field_scale_complement", "field_recovery", "return_reselection", "ij_residual")],
      collapse = " "), "\n")
if (E == "grf") say(grepl("method=grf/frontier", rc, fixed = TRUE), "Run config reads method=grf/frontier")
say(grepl(sprintf("workers=%d", WORKERS), rc, fixed = TRUE), sprintf("Run config workers=%d", WORKERS))

## 3. bundle
cat("\n--- 3. bundle meta ---\n")
say(file.exists(bun), sprintf("bundle present: results/%s", basename(bun)))
say(grepl("_nomr_", basename(bun), fixed = TRUE), "`_nomr` in the bundle filename")
if (file.exists(bun)) {
  b <- readRDS(bun); m <- b$meta; r <- b$results
  for (k in c("mr_inference", "subgroup_method", "sg_focus", "focus_tag", "effect_neighborhood",
              "selection_rule", "campaign_tag", "n_sample", "target_hr_harm", "harm_z1_quantile",
              "harm_prevalence_super", "n_sims", "sim_id_start", "sim_id_end", "seed_base",
              "n_workers", "fb_mode", "er_jcuts", "forestsearch_version", "hostname"))
    kv(k, if (is.null(m[[k]])) "<absent>" else m[[k]])
  say(identical(m$mr_inference, FALSE), "meta$mr_inference is FALSE")
  say(identical(m$subgroup_method, E), sprintf("meta$subgroup_method = %s", m$subgroup_method))
  say(identical(m$sg_focus, FO), sprintf("meta$sg_focus = %s (intended %s)", m$sg_focus, FO))
  say(identical(m$focus_tag, ftag), sprintf("meta$focus_tag = %s (fs_focus_tag: %s)", m$focus_tag, ftag))
  say(isTRUE(all.equal(m$effect_neighborhood, EPS)),
      sprintf("meta$effect_neighborhood = %.2f (%s)", m$effect_neighborhood,
              if (band) "eps 0.20 on a band focus" else "template default 0.10 forwarded; FS_S7_NBHD unset; not read by this focus"))
  say(identical(m$campaign_tag, TAG), sprintf("meta$campaign_tag = %s", m$campaign_tag))
  say(identical(as.integer(m$n_sample), N) && isTRUE(all.equal(m$target_hr_harm, H)) &&
      isTRUE(all.equal(m$harm_z1_quantile, z1q)), "meta n / HR / z1q as intended")
  say(identical(as.integer(m$n_sims), REPS) && identical(as.integer(m$sim_id_start), 1L) &&
      identical(as.integer(m$sim_id_end), REPS), sprintf("meta sim_id 1-%d, one batch", REPS))
  say(identical(as.integer(m$seed_base), 8316951L), "meta$seed_base = 8316951")
  say(identical(as.integer(m$n_workers), WORKERS), sprintf("meta$n_workers = %d", m$n_workers))
  say(identical(nrow(r), REPS), sprintf("bundle rows = %d", nrow(r)))
}

## 4. GRF
if (E == "grf") {
  cat("\n--- 4. GRF: grf_selection, frontier_rule, eps ---\n")
  src <- readLines(file.path(QMD_DIR, "sim_fs_maxeffCons_fb_mr_field_m1_template.qmd"), warn = FALSE)
  lit <- function(nm) { i <- grep(sprintf("^%s\\s*<-", nm), src)[1]
    list(line = i, val = eval(parse(text = sub("#.*$", "", src[i])))) }
  gs <- lit("grf_selection"); gss <- lit("grf_select_statistic")
  kv("template grf_selection", sprintf("\"%s\" (template :%d)", gs$val, gs$line))
  kv("template grf_select_stat.", sprintf("\"%s\" (template :%d)", gss$val, gss$line))
  say(identical(gs$val, "frontier"), "grf_selection = \"frontier\" (source literal)")
  say(identical(gss$val, "effect"), "grf_select_statistic = \"effect\" (source literal)")
  find_switch <- function(expr, var) { out <- NULL
    walk <- function(e) {
      if (!is.null(out) || !is.call(e)) return(invisible(NULL))
      if (length(e) == 3L && is.symbol(e[[1L]]) && as.character(e[[1L]]) %in% c("<-", "=") &&
          is.symbol(e[[2L]]) && identical(as.character(e[[2L]]), var) &&
          is.call(e[[3L]]) && is.symbol(e[[3L]][[1L]]) && identical(as.character(e[[3L]][[1L]]), "switch")) {
        out <<- e[[3L]]; return(invisible(NULL)) }
      for (i in seq_along(e)) tryCatch(walk(e[[i]]), error = function(err) NULL)
      invisible(NULL) }
    walk(expr); out }
  resolve <- function(sw, f) { if (is.null(sw)) return(NA_character_); sw[[2L]] <- f
    tryCatch(eval(sw, baseenv()), error = function(e) paste("ERROR:", conditionMessage(e))) }
  nf <- forestsearch:::.normalize_sg_focus(FO)
  fr_inst <- resolve(find_switch(body(forestsearch::forestsearch), "frontier_rule"), nf)
  ex <- parse(file.path(REPO, "R", "forestsearch_main.R"), keep.source = FALSE)
  sw_src <- NULL; for (e in ex) { sw_src <- find_switch(e, "frontier_rule"); if (!is.null(sw_src)) break }
  fr_src <- resolve(sw_src, nf)
  intended <- c(effMaxSG = "effMaxSG", effMinSG = "effMinSG", maxSG = "maxSG", minSG = "minSG",
                maxeffCons = "eff", maxeff = "eff")[[FO]]
  kv("normalized sg_focus", nf)
  kv("frontier_rule (installed)", fr_inst); kv("frontier_rule (source)", fr_src)
  say(identical(fr_inst, intended) && identical(fr_src, intended),
      sprintf("frontier_rule = %s, installed and source agree (intended %s)", fr_inst, intended))
  if (file.exists(bun)) {
    say(sum(is.finite(r$admitted_n)) > 0L,
        sprintf("admitted_n finite on %d of %d rows -- the effect/frontier re-selection ran",
                sum(is.finite(r$admitted_n)), nrow(r)))
    say(isTRUE(all.equal(m$effect_neighborhood, EPS)),
        sprintf("GRF eps = %.2f%s", m$effect_neighborhood,
                if (band) " (set explicitly; GRF's own default is 0.10)" else " (default; not read by this rule)"))
  }
}
if (E == "dina") {
  src <- readLines(file.path(QMD_DIR, "sim_fs_maxeffCons_fb_mr_field_m1_template.qmd"), warn = FALSE)
  i <- grep("^dina_select_statistic\\s*<-", src)[1]
  cat(sprintf("\n--- 4. DINA ---\n  template :%d  %s\n", i, trimws(src[i])))
}

## 5. thread variables (render.sh exports them for every render)
rs <- readLines(file.path(SCRATCH, "render.sh"), warn = FALSE)
say(all(c("export VECLIB_MAXIMUM_THREADS=1", "export OMP_NUM_THREADS=1", "export OPENBLAS_NUM_THREADS=1") %in% trimws(rs)),
    "render.sh exports VECLIB_MAXIMUM_THREADS / OMP_NUM_THREADS / OPENBLAS_NUM_THREADS = 1")

cat(sprintf("\nGATE A %s: %d passes, %d failures -- %s\n", OUT, pass, fail, if (fail) "FAIL: STOP" else "PASS"))
if (fail) quit(status = 1L)
