# Ratio of the median Wald standard error to the Monte Carlo SD (unconstrained scale), and bias in SD units.
tf <- list(lambda = log, gamma = log, alpha = log, p = function(v) qlogis(v / 0.999))
for (f in sort(list.files("results", "^sim_s.*rds$", full.names = TRUE))) {
  x <- readRDS(f); tr <- attr(x, "truth")
  r <- sapply(names(tf), function(pp) { u <- x[[paste0("u_", pp)]]; se <- x[[paste0("se_", pp)]]; ok <- is.finite(se)
    c(sdratio = median(se[ok]) / sd(u[ok]), zbias = (mean(u[ok]) - tf[[pp]](tr[[pp]])) / sd(u[ok])) })
  cat(sprintf("%-28s SE/SD: %s | bias/SD: %s | mean tries %.2f\n", basename(f),
      paste(sprintf("%.2f", r[1, ]), collapse = " "), paste(sprintf("%.2f", r[2, ]), collapse = " "), mean(x$tries)))
}
