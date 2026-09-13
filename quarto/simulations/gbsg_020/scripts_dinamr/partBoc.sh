#!/bin/zsh
# Part B OC smoke (TASK_partB_enabling_2026-09-12, Part 2).  Runs ONLY after
# Gate T3 passes.
#
# One cell: 12.4% (FS_S7_Z1Q unset), HR 1.50, n 500.  30 replicates, sim_id
# 1-30, MR OFF (FS_S7_MR=FALSE), FB none, 12 workers; render.sh sets the three
# thread variables to 1.  Campaign tag `pBoc` (no campaign glob contains it),
# and FS_S7_MR=FALSE adds `_nomr` to every stem.
#
# Sixteen runs -- the full engine x criterion structure:
#   consistency : effMaxSG effMinSG maxeffCons maxeff maxSG minSG
#   dina, grf   : effMaxSG effMinSG maxSG minSG maxeffCons
# On DINA and GRF `maxeffCons` stands for eff / maxeff / maxeffCons: all three
# resolve to the plain effect argmax there and stem-tag as `eff`
# (fs_focus_tag.R:76-86).  eps 0.20 (FS_S7_NBHD) is set for effMaxSG and
# effMinSG ONLY -- the template now rejects it for any other focus -- which
# also sets it explicitly on GRF, whose own default is 0.10.  GRF's rule is
# carried by frontier_rule under the template literal grf_selection =
# "frontier" (template :503; forestsearch_main.R:2412-2415).
#
# Cap: 90 min wall.  A run is not started once 85 min have elapsed; anything
# skipped is printed as SKIPPED.
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"   # render.sh reads $SP
cd $SP
mkdir -p logs
CAP=$((85 * 60))
BASE=(FS_S7_MR=FALSE FS_S7_FB=none FS_S7_CAMPAIGN=pBoc FS_S7_WORKERS=12
      FS_S7_HR=1.50 FS_S7_N=500 FS_S7_MODE=batch FS_S7_START=1 FS_S7_NSIMS=30)
T0=$(date +%s)
run () {  # $1 = engine (consistency|dina|grf)  $2 = sg_focus
  local E=$1 F=$2 now=$(date +%s)
  if (( now - T0 > CAP )); then echo "SKIPPED (cap): $E $F"; return 0; fi
  local KN=($BASE FS_S7_FOCUS=$F)
  local UNSET=(-u FS_S7_Z1Q -u FS_S7_ER_JCUTS)
  [[ "$F" == effMaxSG || "$F" == effMinSG ]] && KN+=(FS_S7_NBHD=0.20) || UNSET+=(-u FS_S7_NBHD)
  [[ "$E" == consistency ]] && UNSET+=(-u FS_S7_METHOD) || KN+=(FS_S7_METHOD=$E)
  local OUT="partBoc_${E}_${F}"
  echo "=== RUN START: $E $F $(date) ==="
  env $UNSET $KN ./render.sh $OUT || { echo "RUN FAILED: $E $F"; return 1; }
}
for F in effMaxSG effMinSG maxeffCons maxeff maxSG minSG; do run consistency $F; done
for E in dina grf; do
  for F in effMaxSG effMinSG maxSG minSG maxeffCons; do run $E $F; done
done
echo "OC SMOKE COMPLETE $(date) total_wall=$(( $(date +%s) - T0 ))s"
