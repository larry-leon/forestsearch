# Per-cell numbers for the grfmr report, EXTRACTED BY RE-EXECUTING
# summary_grfmr.qmd's OWN CHUNK CODE.
#
# The chunks named below are pulled verbatim out of the .qmd and eval'd, so
# every number this prints is the number the rendered document shows -- not a
# reimplementation that could drift from it.  The only thing added is markdown
# formatting of the resulting data frames.
#
# EVERY coverage figure is coverage of beta(Hhat), computed OVER THE DETECTED
# REPLICATES of the cell it names.  Nothing certifies GRF; no acceptance
# criterion is applied; no recommendation is made.
QMD <- normalizePath(file.path(Sys.getenv("DINAMR_QMD_DIR", unset = ".."),
                               "summary_grfmr.qmd"), mustWork = TRUE)
src <- readLines(QMD, warn = FALSE)
setwd(dirname(QMD))          # the chunks read results/ relative to the .qmd
chunk <- function(nm) {
  st <- grep(sprintf("^```\\{r %s[,}]", nm), src)
  stopifnot(length(st) == 1L)
  en <- st + which(src[(st+1):length(src)] == "```")[1]
  paste(src[(st+1):(en-1)], collapse = "\n") }
run <- function(nm) eval(parse(text = chunk(nm)), envir = globalenv())

suppressPackageStartupMessages({library(forestsearch); library(ggplot2)})
# Only the DEFINITION chunks are re-executed.  The display chunks (inventory,
# cov, strat-k, ...) are deliberately NOT run: they call kbl()/kable_styling()
# to render, which would print their own tables here.  Their computations are
# repeated below from the same definitions, which is what makes these numbers
# the rendered document's numbers.
sink(tempfile())                     # setup's own cat() of the cell inventory
for (nm in c("setup","cov-fns","wilson-fn","strat-fns")) run(nm)
sink()

fmt <- function(x, d = 4) formatC(x, format = "f", digits = d)
md <- function(df, cols, digits = NULL, caption = NULL) {
  df <- df[, cols, drop = FALSE]
  if (!is.null(caption)) cat("\n", caption, "\n\n", sep = "")
  for (j in names(df)) {
    v <- df[[j]]
    if (is.factor(v)) { df[[j]] <- as.character(v); next }
    if (!is.numeric(v)) next
    # integers (counts, n, n_sample) print as integers, not as 1997.0000
    d <- if (!is.null(digits)) digits else if (isTRUE(all(v == round(v), na.rm = TRUE))) 0 else 4
    df[[j]] <- formatC(v, format = "f", digits = d) }
  cat("|", paste(names(df), collapse = " | "), "|\n")
  cat("|", paste(rep("---", ncol(df)), collapse = " | "), "|\n")
  for (i in seq_len(nrow(df)))
    cat("|", paste(as.character(unlist(df[i, ])), collapse = " | "), "|\n") }

## ---------------------------------------------------------------- 1. standard
COV <- do.call(rbind, lapply(labs, function(k) { r <- B[[k]]$results
  rbind(cbind(cell = k, block_grp = BLK[[k]], n_sample = NN[[k]], hr_cell = HR[[k]],
              block = "Hhat (lower)",   cov_block(r, "H",  "lower")),
        cbind(cell = k, block_grp = BLK[[k]], n_sample = NN[[k]], hr_cell = HR[[k]],
              block = "Hhat^c (upper)", cov_block(r, "Hc", "upper"))) }))
COV$construction <- factor(COV$construction, levels = c("naive","field","field-s","IJ two-term"))
COV <- COV[order(match(COV$cell, labs), COV$block, COV$construction), ]
cat("\n## 1. STANDARD TABLE -- every construction, both blocks, every cell\n")
for (blk in c("A","B")) {
  sub <- COV[COV$block_grp == blk, ]
  md(sub, c("cell","block","construction","n","bias_log","sd_emp","sd_err","b","b_err",
            "se_mean","r","se_over_sd_err","cov1","cov1_wilson_lo","cov1_wilson_hi",
            "cov2","cov1_ref","cov1_ref_err"),
     caption = sprintf("### Block %s (%s prevalence)", blk, if (blk=="A") "12.4%" else "31%")) }

## ---------------------------------------------------------------- 2. products + joint
cat("\n\n## 2. EVERY PRODUCT incl. the Bonferroni joint, with the IJ miss split by side\n")
md(VSN, c("cell","n_eval","harm_field","harm_field_lo","harm_field_hi",
          "comp_field","comp_field_lo","comp_field_hi",
          "comp_field_s","comp_field_s_lo","comp_field_s_hi"),
   caption = "### 2a. field and field-s, one-sided on the exposed side [Wilson]")
md(VSN, c("cell","n_eval","ij2_H","ij2_H_lo","ij2_H_hi","miss_below","miss_above",
          "ij2_Hc","ij2_Hc_lo","ij2_Hc_hi"),
   caption = "### 2b. IJ two-term two-sided, with the Hhat miss split by side")
md(VSN, c("cell","n_eval","joint_bonf","joint_bonf_lo","joint_bonf_hi",
          "joint_s_bonf","joint_s_bonf_lo","joint_s_bonf_hi"),
   caption = "### 2c. Bonferroni joint (harm lower, complement upper) [Wilson]")

## ---------------------------------------------------------------- 3. strata
SK <- do.call(rbind, lapply(labs, function(k) { d <- det_rows(k)
  do.call(rbind, lapply(strata_of(d, "K"), function(S) prod_rows(S$d, k, S$lab, S$ov))) }))
SP <- do.call(rbind, lapply(labs, function(k) { d <- det_rows(k)
  do.call(rbind, lapply(strata_of(d, "p"), function(S) prod_rows(S$d, k, S$lab, S$ov))) }))
fld <- function(X) X[X$block == "Hhat" & X$product == "field", ]
cat("\n\n## 3. THE FIELD LOWER BOUND by stratum\n")
md(fld(SK), c("cell","stratum","overlapping","n","coverage","wilson_lo","wilson_hi",
              "retained_bias_log"),
   caption = "### 3a. by admitted_n tertile (the stratifier)")
md(fld(SP), c("cell","stratum","overlapping","n","coverage","wilson_lo","wilson_hi",
              "retained_bias_log"),
   caption = "### 3b. by p-hat bin (the FS summaries' within-cell tertiles)")

## ---------------------------------------------------------------- 4. joint counts
XT <- do.call(rbind, lapply(labs, function(k) { d <- det_rows(k)
  gk <- tert_fs(d$admitted_n); gp <- tert_fs(d$p_hat_H)
  base <- do.call(rbind, lapply(sort(unique(gk[!is.na(gk)])), function(a)
    do.call(rbind, lapply(sort(unique(gp[!is.na(gp)])), function(b)
      data.frame(cell = k, adm_stratum = sprintf("adm T%d", a),
                 p_hat_bin = sprintf("p-hat T%d", b),
                 count = sum(!is.na(gk) & !is.na(gp) & gk == a & gp == b),
                 stringsAsFactors = FALSE)))))
  extra <- do.call(rbind, lapply(list(c("admitted_n = 1","1"), c("admitted_n <= 5","5")), function(x)
    do.call(rbind, lapply(sort(unique(gp[!is.na(gp)])), function(b)
      data.frame(cell = k, adm_stratum = paste0(x[1], " (overlapping)"),
                 p_hat_bin = sprintf("p-hat T%d", b),
                 count = sum(d$admitted_n <= as.integer(x[2]) &
                             (if (x[2]=="1") d$admitted_n == 1L else TRUE) &
                             !is.na(gp) & gp == b), stringsAsFactors = FALSE)))))
  rbind(base, extra) }))
XW <- reshape(XT, idvar = c("cell","adm_stratum"), timevar = "p_hat_bin", direction = "wide")
names(XW) <- sub("^count\\.", "", names(XW))
cat("\n\n## 4. JOINT COUNTS -- admitted_n stratum by p-hat bin\n")
md(XW, names(XW), digits = 0, caption = "")

## ---------------------------------------------------------------- 5. location
cat("\n\n## 5. BOUND LOCATION ON THE HR SCALE\n")
LOC <- do.call(rbind, lapply(labs, function(k) { d <- det_rows(k)
  data.frame(cell = k, n_eval = nrow(d),
             med_bound = stats::median(d$fld_H_lo1s),
             med_theta = stats::median(d$betaHhat_H),
             med_est2  = stats::median(d$fld_H_est2),
             med_naive = stats::median(d$nv_H_est),
             gap_diff  = stats::median(d$fld_H_lo1s) - stats::median(d$betaHhat_H),
             gap_ratio = stats::median(d$fld_H_lo1s) / stats::median(d$betaHhat_H),
             paired_ratio = stats::median(d$fld_H_lo1s / d$betaHhat_H),
             planted_marg_H = unname(B[[k]]$truth$marg_H),
             share_ge_100 = mean(d$fld_H_lo1s >= 1.00),
             s100_lo = wil(mean(d$fld_H_lo1s >= 1.00), nrow(d))[1],
             s100_hi = wil(mean(d$fld_H_lo1s >= 1.00), nrow(d))[2],
             share_ge_125 = mean(d$fld_H_lo1s >= 1.25),
             s125_lo = wil(mean(d$fld_H_lo1s >= 1.25), nrow(d))[1],
             s125_hi = wil(mean(d$fld_H_lo1s >= 1.25), nrow(d))[2],
             stringsAsFactors = FALSE) }))
md(LOC, names(LOC), caption = "")
saveRDS(list(COV=COV, VSN=VSN, SK=SK, SP=SP, XW=XW, LOC=LOC),
        file.path("scripts_dinamr","grfmr_tables.rds"))
