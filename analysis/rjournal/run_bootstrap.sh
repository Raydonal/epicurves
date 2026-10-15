#!/bin/bash
# Runs example_bootstrap.R (B = 200) and ends with BOOTDONE.
cd "$(dirname "$0")"
echo "$(date '+%F %T') start example_bootstrap"
Rscript example_bootstrap.R > results/example_bootstrap.log 2>&1
echo BOOTDONE
