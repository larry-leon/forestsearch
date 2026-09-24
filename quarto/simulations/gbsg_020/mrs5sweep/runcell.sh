#!/bin/bash
# usage: runcell.sh <tag> <cellid> <HR> <N> <NSIMS> [Z1Q] [LIBNAME]
set -u
TAG=$1; CELL=$2; HR=$3; N=$4; NS=$5; Z=${6:-}; LIBN=${7:-Rlib_head}
ZK=(); [ -n "$Z" ] && ZK=(FS_S7_Z1Q=$Z)
cd /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020
export SP=$PWD/$TAG PATH=/usr/lib/rstudio/resources/app/bin/quarto/bin:$PATH
export R_LIBS=/tmp/claude-1000/-home-larryleon-Documents-GitHub-forestsearch/ab0a4121-64a0-4c61-972b-998577a64b16/scratchpad/$LIBN
mkdir -p $SP/logs
OUT=${TAG}_${CELL}_batch_1
# Gate: assert from inside workers that the HEAD build is what loads
export LIBN; WANT=TRUE; [ "$LIBN" = Rlib_pre ] && WANT=FALSE; export WANT
Rscript -e '
library(future); plan(multisession, workers = 64)
fs <- lapply(1:64, function(i) future({ f <- get("fs_mr_inference", asNamespace("forestsearch"))
  c(Sys.getpid(), find.package("forestsearch"), "pconsistency.digits" %in% names(formals(f)),
    any(grepl(".fs_pcons_eff", deparse(f), fixed = TRUE))) }, seed = TRUE))
v <- do.call(rbind, lapply(fs, value)); plan(sequential)
lib <- unique(v[,2]); ok <- length(unique(v[,1])) >= 32 && length(lib) == 1 &&
  grepl(Sys.getenv("LIBN"), lib) && all(v[,3] == Sys.getenv("WANT")) && all(v[,4] == Sys.getenv("WANT"))
cat(sprintf("WORKER_BUILD_CHECK pids=%d lib=%s fix_formal=%s pcons_eff=%s OK=%s\n",
  length(unique(v[,1])), paste(lib, collapse="|"), all(v[,3]=="TRUE"), all(v[,4]=="TRUE"), ok))
if (!ok) quit(status = 3)' 2>&1 | tee $SP/logs/${OUT}.buildcheck | grep WORKER_BUILD_CHECK
[ ${PIPESTATUS[0]} -eq 0 ] || { echo "BUILD CHECK FAILED"; exit 3; }
UNSET=(); for v in $(env | grep -o '^FS_S7_[A-Z0-9_]*'); do UNSET+=(-u "$v"); done
env "${UNSET[@]}" FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected \
  FS_S7_FIELD_DECOMP=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=$TAG FS_S7_RETURN_RESEL=TRUE \
  "${ZK[@]}" FS_S7_WORKERS=64 FS_S7_HR=$HR FS_S7_N=$N FS_S7_MODE=batch FS_S7_START=1 FS_S7_NSIMS=$NS \
  bash scripts_dinamr/render.sh $OUT
