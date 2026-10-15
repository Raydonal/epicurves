# Chance baseline for detection rates: cuts that each rule places in the no-change series,
# scored against the true change points of the continuous scenarios (tolerance 7 days).
suppressMessages(library(epicurves))
source("scenarios.R")
rules <- c("thesis_flat", "cum_bic", "inc_bic", "edetect_z0", "edetect_z1")
nc <- c(readRDS("results/mc_nc_log.rds")$res, readRDS("results/mc_nc_gom.rds")$res)
for (s in 1:3) {
  tc <- scenario_curve(s, continuous = TRUE)$cuts
  r <- sapply(rules, function(rl) mean(sapply(nc, function(x) sapply(tc, function(c) any(abs(x[[rl]]$cuts - c) <= 7)))))
  cat(sprintf("scenario %dc cuts (%s): %s\n", s, paste(tc, collapse = " "),
              paste(sprintf("%s=%.2f", rules, r), collapse = "  ")))
}
