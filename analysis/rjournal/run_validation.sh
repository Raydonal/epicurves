#!/bin/bash
cd "$(dirname "$0")"
for m in seird sirdv; do
  [ -f results/validate_$m.rds ] && continue
  echo "$(date '+%F %T') start $m"
  Rscript validate_generic.R $m 100 4
done
echo ALLDONE
