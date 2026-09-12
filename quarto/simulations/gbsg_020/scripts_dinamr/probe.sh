#!/bin/zsh
# Gate 1 cost probes: 36 replicates per (prevalence, n, HR) corner, 12 workers.
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"   # override with DINAMR_SCRATCH
   # EXPORTED: render.sh reads $SP for its log directory and its header
   # states the caller exports what it needs.  Unexported, a clean-environment
   # run resolved $SP to empty and died on `mkdir -p /logs`.
cd $SP
run () {  # $1=z1q ("" = unset/12.4%)  $2=n  $3=hr  $4=label
  local Z=$1 N=$2 H=$3 L=$4
  if [[ -z "$Z" ]]; then
    env -u FS_S7_Z1Q -u FS_S7_ER_JCUTS \
      FS_S7_METHOD=dina FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE \
      FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE \
      FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=dinamrprobe \
      FS_S7_HR=$H FS_S7_N=$N FS_S7_NSIMS=36 FS_S7_START=1 FS_S7_WORKERS=12 \
      ./render.sh probe_$L
  else
    env -u FS_S7_ER_JCUTS \
      FS_S7_METHOD=dina FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE \
      FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE \
      FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=dinamrprobe \
      FS_S7_Z1Q=$Z FS_S7_HR=$H FS_S7_N=$N FS_S7_NSIMS=36 FS_S7_START=1 FS_S7_WORKERS=12 \
      ./render.sh probe_$L
  fi
  echo "PROBE DONE: $L"
}
run ""    500  1.50 p124_h150_n500
run ""   1000  1.50 p124_h150_n1000
run ""   1500  1.50 p124_h150_n1500
run 0.60  500  1.50 p31_h150_n500
run 0.60 1000  1.50 p31_h150_n1000
run 0.60 1500  1.50 p31_h150_n1500
run ""    500  1.00 p124_h100_n500
run ""   1500  1.00 p124_h100_n1500
run 0.60  500  1.00 p31_h100_n500
run 0.60 1500  1.00 p31_h100_n1500
echo "ALL PROBES DONE"
