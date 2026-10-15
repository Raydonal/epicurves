# Summary of the Paraíba parametric bootstrap (paraiba_bootstrap.R) against the Wald-type intervals
# of paraiba_report.R (Table 3). Percentile intervals, and basic intervals 2*est - quantiles.
fs <- list.files("results", "^paraiba_boot_[0-9]+\\.rds$", full.names = TRUE)
reps <- unlist(lapply(fs, readRDS), recursive = FALSE)
ok <- vapply(reps, function(x) x[["ok"]] == 1, logical(1))
M <- do.call(rbind, reps[ok])
rep_tab <- readRDS("results/paraiba_report.rds")$table
q <- c("inf_period_days", "ifr", "p", "attack", "i0N", "tau2", "tau2_d", paste0("R_block", 1:10))
est <- rep_tab[q, "estimate"]
lo <- apply(M[, q], 2, quantile, 0.025); hi <- apply(M[, q], 2, quantile, 0.975)
out <- data.frame(estimate = est, wald_lower = rep_tab[q, "lower"], wald_upper = rep_tab[q, "upper"],
                  boot_mean = colMeans(M[, q]), boot_median = apply(M[, q], 2, median),
                  pct_lower = lo, pct_upper = hi, basic_lower = 2 * est - hi, basic_upper = 2 * est - lo,
                  row.names = q)
out$width_ratio <- (out$pct_upper - out$pct_lower) / (out$wald_upper - out$wald_lower)
cat("replicates:", length(reps), " successful fits:", sum(ok), " nonzero convergence codes:",
    sum(M[, "conv"] != 0), "\n")
cat("median refit time (s):", median(M[, "secs"]), "\n")
cat("true attack rate in the simulated epidemics: median", signif(median(M[, "attack_true"]), 3),
    " range", signif(range(M[, "attack_true"]), 3), "\n")
print(signif(out, 3))
saveRDS(list(table = out, reps = M), "results/paraiba_bootstrap_summary.rds")
