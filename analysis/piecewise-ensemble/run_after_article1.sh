#!/bin/bash
# Espera a grade do Artigo 1 (resume_all.sh) terminar e roda o Monte Carlo do Artigo 2 (J=60),
# pulando cenários já concluídos com J=60.
cd "$(dirname "$0")"
while pgrep -f "resume_all.sh|sim_study.R" > /dev/null; do sleep 120; done
echo "$(date '+%F %T') Article 1 grid finished; starting Article 2 MC"
for sc in 1c 2c 3c nc_log nc_gom 1 2 3; do
  f=results/mc_${sc}.rds
  if [ -f $f ] && Rscript -e "quit(status = as.integer(length(readRDS('$f')\$res) < 60))"; then echo "skip $sc"; continue; fi
  echo "$(date '+%F %T') start mc $sc"; Rscript mc_study.R $sc 60 4
done
echo ALLDONE
