#!/bin/bash
# Run-length sensitivity of the e-detector (mc_arl.R). Waits for the overdispersion run (NBDONE in
# results/mc_negbin.log) so the two do not compete for the 4 cores; skips finished scenarios.
# Log: results/mc_arl.log (ARLDONE at the end); summary in results/mc_arl_summary.log.
cd "$(dirname "$0")"
until grep -q NBDONE results/mc_negbin.log 2>/dev/null; do sleep 60; done
for sc in nc_log nc_gom 1c 2c 3c; do
  [ -f results/mc_arl_${sc}.rds ] && continue
  echo "$(date '+%F %T') start $sc"
  Rscript mc_arl.R $sc 4
done
Rscript summarise_arl.R > results/mc_arl_summary.log 2>&1
echo ARLDONE
