# Paraíba: sensitivity of the time-varying SIRD-GA fit to the block length of lambda(t).
# Blocks of L = 2, 3, 4, 5, 6, 8 weeks over the 40 complete epidemiological weeks (13-52 of 2020),
# block of week k = ceiling(k / L) as in paraiba_tv.R, so when L does not divide 40 the last block
# is shorter (L = 3: 14 blocks, the last of 1 week; L = 6: 7 blocks, the last of 4 weeks).
# Each L is fitted from the three generic starts of paraiba_tv.R (one core per start); the best
# optimum is reported. Number of parameters k = (#blocks) + 6 (gamma, alpha, i0, p, tau2, tau2_d).
# AIC = -2 ll + 2k; BIC = -2 ll + k log(n) with n = 80 scalar observations (40 weeks x 2 series).
# Output: results/paraiba_blocks_L<L>.rds (one per L, skipped if present) and
# results/paraiba_blocks.rds (summary table).
# Usage: Rscript paraiba_blocks.R [cores]
source("regime_helpers.R")
args <- commandArgs(TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 3L
dat <- paraiba_weekly(); n <- length(dat$cases); N <- dat$N
stopifnot(n == 40)
Ls <- c(2, 3, 4, 5, 6, 8)
for (L in Ls) {
  out_file <- sprintf("results/paraiba_blocks_L%d.rds", L)
  if (file.exists(out_file)) next
  blk <- ceiling(seq_len(n) / L)
  t0 <- proc.time()[3]
  fits <- parallel::mclapply(generic_starts_tv(max(blk)), function(st)
    tryCatch(sird_ga_fit_tv(dat$cases, dat$deaths, N, blk, "both", st, Delta = 1, nsub = 20L),
             error = function(e) NULL), mc.cores = cores)
  saveRDS(list(L = L, blk = blk, fits = fits), out_file)
  cat(format(Sys.time()), sprintf("L = %d: %d blocks, logliks %s (%.0fs)\n", L, max(blk),
      paste(sprintf("%.2f", sapply(fits, function(f) if (is.null(f)) NA else f$loglik)), collapse = " "),
      proc.time()[3] - t0))
}

tab <- do.call(rbind, lapply(Ls, function(L) {
  x <- readRDS(sprintf("results/paraiba_blocks_L%d.rds", L))
  ll <- sapply(x$fits, function(f) if (is.null(f)) NA else f$loglik)
  f <- x$fits[[which.max(ll)]]
  B <- max(x$blk); k <- B + 6; nobs <- 80
  est <- derived_from_par(f$par, f$sc_names, dat$cases, dat$deaths, N, x$blk)
  ev <- min(eigen((f$hessian + t(f$hessian)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  data.frame(L = L, blocks = B, last_block_weeks = sum(x$blk == B), k = k,
             loglik = f$loglik, AIC = -2 * f$loglik + 2 * k, BIC = -2 * f$loglik + k * log(nobs),
             ll_spread_starts = diff(range(ll, na.rm = TRUE)),
             p = est[["p"]], ifr = est[["ifr"]], inf_period_days = est[["inf_period_days"]],
             attack = est[["attack"]], i0N = est[["i0N"]], tau2 = est[["tau2"]],
             tau2_d = est[["tau2_d"]], R_block1 = est[["R_block1"]],
             conv = f$convergence, hess_min_eig = ev)
}))
tab$dAIC <- tab$AIC - min(tab$AIC); tab$dBIC <- tab$BIC - min(tab$BIC)
options(width = 200)
cat("\nBlock-length sensitivity (n = 80 for BIC)\n")
print(format(tab, digits = 4), row.names = FALSE)
saveRDS(tab, "results/paraiba_blocks.rds")
cat("BLOCKSDONE\n")
