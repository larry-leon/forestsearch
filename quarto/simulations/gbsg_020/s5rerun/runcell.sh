#!/bin/bash
# Section 5 full re-run (TASK_section5_full_rerun_2026-09-24), Mac Studio runner.
# Mirrors mrs5sweep/runcell.sh (itself the p12x20 / cert20 knob set) with these
# deviations only: Mac paths; the installed library (R_LIBS unset) instead of a
# scratch R_LIBS; worker count from $S5_WORKERS (template caps at physical-1 = 13);
# MODE/START/NSIMS are arguments so the same runner serves batch and combine.
# usage: runcell.sh <tag> <cellid> <HR> <N> <MODE> <START> <NSIMS> [Z1Q]
set -u
TAG=$1; CELL=$2; HR=$3; N=$4; MODE=$5; START=$6; NS=$7; Z=${8:-}
: "${S5_WORKERS:?S5_WORKERS must be set}"
ZK=(); [ -n "$Z" ] && ZK=(FS_S7_Z1Q=$Z)
cd /Users/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020
export SP=$PWD/s5rerun
unset R_LIBS
mkdir -p $SP/logs
OUT=${TAG}_${CELL}_${MODE}_${START}
Rscript -e '
library(future); W <- as.integer(Sys.getenv("S5_WORKERS")); plan(multisession, workers = W)
fs <- lapply(seq_len(W), function(i) future({ f <- get("fs_mr_inference", asNamespace("forestsearch"))
  c(Sys.getpid(), find.package("forestsearch"), "pconsistency.digits" %in% names(formals(f)),
    any(grepl(".fs_pcons_eff", deparse(f), fixed = TRUE)), packageDescription("forestsearch")$Built) }, seed = TRUE))
v <- do.call(rbind, lapply(fs, value)); plan(sequential)
lib <- unique(v[,2]); ok <- length(unique(v[,1])) >= min(W, 2) && length(lib) == 1 && all(v[,3] == "TRUE") && all(v[,4] == "TRUE")
cat(sprintf("WORKER_BUILD_CHECK pids=%d lib=%s built=%s fix_formal=%s pcons_eff=%s OK=%s\n",
  length(unique(v[,1])), paste(lib, collapse="|"), paste(unique(v[,5]), collapse="|"), all(v[,3]=="TRUE"), all(v[,4]=="TRUE"), ok))
if (!ok) quit(status = 3)' 2>&1 | tee $SP/logs/${OUT}.buildcheck | grep WORKER_BUILD_CHECK
[ ${PIPESTATUS[0]} -eq 0 ] || { echo "BUILD CHECK FAILED"; exit 3; }
UNSET=(); for v in $(env | grep -o '^FS_S7_[A-Z0-9_]*'); do UNSET+=(-u "$v"); done
env ${UNSET[@]+"${UNSET[@]}"} FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected \
  FS_S7_FIELD_DECOMP=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=$TAG FS_S7_RETURN_RESEL=TRUE \
  ${ZK[@]+"${ZK[@]}"} FS_S7_WORKERS=$S5_WORKERS FS_S7_HR=$HR FS_S7_N=$N FS_S7_MODE=$MODE FS_S7_START=$START FS_S7_NSIMS=$NS \
  zsh scripts_dinamr/render.sh $OUT
