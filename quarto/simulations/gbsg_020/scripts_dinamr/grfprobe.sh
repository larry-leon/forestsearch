#!/bin/zsh
# Part B GRF cost probes (TASK_dinamr_blockC_grfprobe_2026-09-11).
# Five 36-replicate corners, 12 workers, tag `grfprobe`.
#
# COST MEASUREMENT ONLY.  No coverage claim, no acceptance criterion, no
# recommendation, and NOT a GRF campaign.  GRF is not "FS-analogous": a
# GRF-to-FS or GRF-to-DINA comparison differs in identifier, family
# construction, detection set, selection criterion AND the scale of the
# selection criterion.
#
# dmin.grf = 0.0 is the template's own value (line 506) and is Larry's decision
# of 2026-09-11; it is not overridable from the environment, so this driver
# does not set it.  The rationale, AND the qualification tracing the path adds
# to it, are recorded in full at the head of grfprobe.R: dmin.grf is a DR-score
# PRE-FILTER (grf_main.R:291, grf_subg_harm_glm.R:523 ->
# grf_subgroup_labels.R:358), while the BINDING effect-scale floor on the
# re-selection path is hr.threshold = 0.90 -- the same floor DINA carries
# (forestsearch_helpers.R:1632-1635, forestsearch_main.R:2026-2030, template
# line 531).  So "GRF's floor is not alignable with DINA's" is true of
# dmin.grf and only of dmin.grf.
#
# Peak memory: each render runs under /usr/bin/time -l, whose "maximum resident
# set size" line is appended to the probe's own log, and a 5-second sampler
# records the summed RSS of the whole quarto process tree beside it.
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"   # override with DINAMR_SCRATCH
   # EXPORTED: render.sh reads $SP for its log directory and its header
   # states the caller exports what it needs.  Unexported, a clean-environment
   # run resolved $SP to empty and died on `mkdir -p /logs`.
cd $SP
mkdir -p logs
sample () {  # $1 = label; writes total-RSS samples (KB) until killed
  while :; do
    ps -Ao rss,command | grep -E "quarto|/R$|exec/R" | grep -v grep \
      | awk -v t="$(date +%s)" '{s+=$1} END {print t, s+0}' >> logs/rss_$1.txt
    sleep 5
  done
}
run () {  # $1=z1q ("" = unset/12.4%)  $2=n  $3=hr  $4=label
  local Z=$1 N=$2 H=$3 L=$4
  : > logs/rss_$L.txt
  sample $L & local SPID=$!
  local KN=(FS_S7_METHOD=grf FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE
            FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE
            FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=grfprobe
            FS_S7_HR=$H FS_S7_N=$N FS_S7_NSIMS=36 FS_S7_START=1 FS_S7_WORKERS=12)
  if [[ -z "$Z" ]]; then
    env -u FS_S7_Z1Q -u FS_S7_ER_JCUTS $KN \
      /usr/bin/time -l ./render.sh grfprobe_$L 2>> logs/time_$L.txt
  else
    env -u FS_S7_ER_JCUTS $KN FS_S7_Z1Q=$Z \
      /usr/bin/time -l ./render.sh grfprobe_$L 2>> logs/time_$L.txt
  fi
  kill $SPID 2>/dev/null
  echo "GRF PROBE DONE: $L $(date)"
}
run ""    500 1.50 g_p124_h150_n500
run ""   1500 1.50 g_p124_h150_n1500
run 0.60  500 1.50 g_p31_h150_n500
run 0.60 1500 1.50 g_p31_h150_n1500
run ""    500 1.00 g_p124_h100_n500
echo "ALL GRF PROBES DONE $(date)"
