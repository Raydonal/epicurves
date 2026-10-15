#!/bin/bash
# Overdispersion Monte Carlo (mc_negbin.R): no-change Logistic/Gompertz, phi = 3 and 10, J = 60.
# Skips configurations whose .rds already exists. Log: results/mc_negbin.log (NBDONE at the end).
cd "$(dirname "$0")"
for phi in 3 10; do for cv in log gom; do
  [ -f results/mc_negbin_${cv}_phi${phi}.rds ] && continue
  echo "$(date '+%F %T') start $cv phi=$phi"
  Rscript mc_negbin.R $cv $phi 60 4
done; done
echo NBDONE
