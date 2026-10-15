# Summary of the overdispersion Monte Carlo (results/mc_negbin_*.rds): per rule and configuration,
# share of series with at least one change point and mean number of change points (all are false,
# since the series have no change), and the timing of the first false alarm of the sub-Gaussian
# detector with z = 1. Output: results/mc_negbin_summary.{log,rds}.
dir <- "results"
cfg <- expand.grid(curve = c("log", "gom"), phi = c(3, 10), stringsAsFactors = FALSE)
rules <- c("pl_poisson", "pl_negbin", "ed_poisson_z1", "ed_subgaussian_z0", "ed_subgaussian_z1")
out <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(i) {
  r <- readRDS(file.path(dir, sprintf("mc_negbin_%s_phi%g.rds", cfg$curve[i], cfg$phi[i])))$res
  do.call(rbind, lapply(rules, function(rl) {
    k <- vapply(r, function(z) length(z[[rl]]$cuts), numeric(1))
    data.frame(curve = cfg$curve[i], phi = cfg$phi[i], rule = rl, J = length(r),
               any_cut = mean(k > 0), mean_cuts = mean(k))
  }))
}))
print(out, digits = 3, row.names = FALSE)
first <- unlist(lapply(Sys.glob(file.path(dir, "mc_negbin_*.rds")), function(f)
  unlist(lapply(readRDS(f)$res, function(z) if (length(z$ed_subgaussian_z1$alarms)) min(z$ed_subgaussian_z1$alarms)))))
cat(sprintf("\nsub-Gaussian z = 1: %d series with a false alarm; first alarm by day 30 (first 10 monitored days) in %d (%.2f)\n",
            length(first), sum(first <= 30), mean(first <= 30)))
saveRDS(list(table = out, first_alarm_sg_z1 = first), file.path(dir, "mc_negbin_summary.rds"))
