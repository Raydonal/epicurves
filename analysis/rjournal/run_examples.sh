#!/bin/bash
cd "$(dirname "$0")"
while ! grep -q ALLDONE results/validation_run.log; do sleep 60; done
echo "$(date '+%F %T') start example_seird"
Rscript example_seird.R > results/example_seird.log 2>&1
echo "$(date "+%F %T") start example_paraiba"
Rscript example_paraiba.R > results/example_paraiba.log 2>&1
echo "$(date "+%F %T") start timing"
Rscript timing.R > results/timing.log 2>&1
echo EXAMPLESDONE
