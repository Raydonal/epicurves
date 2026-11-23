# Paraíba block-length sensitivity, second pass. paraiba_blocks.R fitted each block length from
# three generic starts; for L = 2 and L = 5 those starts disagreed, so the reported maxima may be
# local. Here each L is refitted from the optima of the other block lengths: the weekly lambda path
# of the other fit is mapped to the new blocks by averaging log lambda over the weeks of each block
# (exact when the other blocks are unions of the new ones, e.g. 4 -> 2), the other parameters are
# copied. The best of generic and warm starts is reported, with the same columns as paraiba_blocks.R
# and a check that nested models are ordered (L = 2 >= 4 >= 8, 2 >= 6, 3 >= 6).
# Output: results/paraiba_blocks_warm.rds.   Usage: Rscript paraiba_blocks_warm.R [cores]
source("regime_helpers.R")
args <- commandArgs(TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 3L
dat <- paraiba_weekly(); n <- length(dat$cases); N <- dat$N
Ls <- c(2, 3, 4, 5, 6, 8)
gen <- lapply(Ls, function(L) readRDS(sprintf("results/paraiba_blocks_L%d.rds", L)))
names(gen) <- Ls
best_of <- function(fits) { ll <- sapply(fits, function(f) if (is.null(f)) -Inf else f$loglik); fits[[which.max(ll)]] }
gbest <- lapply(gen, function(x) best_of(x$fits))

out_file <- "results/paraiba_blocks_warm_fits.rds"
if (!file.exists(out_file)) {
  jobs <- expand.grid(to = Ls, from = Ls); jobs <- jobs[jobs$to != jobs$from, ]
  fits <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
    L <- jobs$to[j]; src <- gbest[[as.character(jobs$from[j])]]
    blk <- ceiling(seq_len(n) / L)
    wl <- log(src$lambda[src$blocks])                       # weekly log lambda of the source fit
    st <- c(list(lambda = exp(tapply(wl, blk, mean))), as.list(src$theta[c("gamma", "alpha", "i0", "p", "tau2", "tau2_d")]))
    tryCatch(sird_ga_fit_tv(dat$cases, dat$deaths, N, blk, "both", st, Delta = 1, nsub = 20L),
             error = function(e) NULL)
  }, mc.cores = cores)
  saveRDS(list(jobs = jobs, fits = fits), out_file)
}
w <- readRDS(out_file)
cat("Warm-start log-likelihoods (rows: fitted L, columns: source L; generic best in column 'generic')\n")
M <- matrix(NA, length(Ls), length(Ls), dimnames = list(paste0("L", Ls), paste0("from", Ls)))
for (j in seq_len(nrow(w$jobs))) M[paste0("L", w$jobs$to[j]), paste0("from", w$jobs$from[j])] <-
  if (is.null(w$fits[[j]])) NA else w$fits[[j]]$loglik
print(round(cbind(generic = sapply(gbest, `[[`, "loglik"), M), 2))

tab <- do.call(rbind, lapply(Ls, function(L) {
  cand <- c(gen[[as.character(L)]]$fits, w$fits[w$jobs$to == L])
  src <- c(paste0("generic", 1:3), paste0("warm_from_L", w$jobs$from[w$jobs$to == L]))
  ll <- sapply(cand, function(f) if (is.null(f)) -Inf else f$loglik)
  f <- cand[[which.max(ll)]]; blk <- ceiling(seq_len(n) / L)
  B <- max(blk); k <- B + 6; nobs <- 80
  est <- derived_from_par(f$par, f$sc_names, dat$cases, dat$deaths, N, blk)
  ev <- min(eigen((f$hessian + t(f$hessian)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  data.frame(L = L, blocks = B, last_block_weeks = sum(blk == B), k = k, best_from = src[which.max(ll)],
             loglik = f$loglik, loglik_generic = max(ll[1:3]),
             AIC = -2 * f$loglik + 2 * k, BIC = -2 * f$loglik + k * log(nobs),
             p = est[["p"]], ifr = est[["ifr"]], inf_period_days = est[["inf_period_days"]],
             attack = est[["attack"]], i0N = est[["i0N"]], tau2 = est[["tau2"]], tau2_d = est[["tau2_d"]],
             R_block1 = est[["R_block1"]], conv = f$convergence, hess_min_eig = ev)
}))
tab$dAIC <- tab$AIC - min(tab$AIC); tab$dBIC <- tab$BIC - min(tab$BIC)
options(width = 220)
cat("\nBlock-length sensitivity, best of generic and warm starts (n = 80 for BIC)\n")
print(format(tab, digits = 4), row.names = FALSE)
ll <- setNames(tab$loglik, tab$L)
cat("\nnesting checks (finer >= coarser):",
    sprintf("2>=4 %s, 4>=8 %s, 2>=8 %s, 2>=6 %s, 3>=6 %s", ll["2"] >= ll["4"] - 1e-6, ll["4"] >= ll["8"] - 1e-6,
            ll["2"] >= ll["8"] - 1e-6, ll["2"] >= ll["6"] - 1e-6, ll["3"] >= ll["6"] - 1e-6), "\n")
saveRDS(list(table = tab, warm_ll = M), "results/paraiba_blocks_warm.rds")
cat("BLOCKSWARMDONE\n")
