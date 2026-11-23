# Paraíba: refined profile log-likelihood of the case reporting probability p (4-week blocks).
# The earlier profile (paraiba_tv_profile.R, results/paraiba_tv_profile_p.rds) used 10 grid values,
# each fitted from the overall optimum only; neighbouring fitted infectious periods and R_t jump,
# which suggests different local maxima. Here, on a finer grid (18 log-spaced values between 0.02
# and 0.8, plus 0.066 and 0.10), every grid value is fitted with p fixed from
#   (a) the three generic starts of paraiba_tv.R (phase "generic", grid values in parallel),
#   (b) a forward sweep: warm start from the best fit so far at the previous (smaller) grid value,
#   (c) a backward sweep: warm start from the best fit so far at the next (larger) grid value;
# the two sweeps run in parallel. The profile is the best of the five fits.
# Profile interval: points with ll >= max - 1.92 (max over the profile and the unconstrained
# optimum of paraiba_tv.rds), limits by linear interpolation of ll between grid values. The same
# rule is applied to the old profile for comparison.
# Output: results/paraiba_profile_p_refined_{generic,fwd,bwd}.rds (phases, skipped if present),
# results/paraiba_profile_p_refined.rds.   Usage: Rscript paraiba_profile_p_refined.R [cores]
source("regime_helpers.R")
args <- commandArgs(TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 3L
r <- readRDS("results/paraiba_tv.rds")
fbest <- r$res[[which.max(sapply(r$res, `[[`, "loglik"))]]
cases <- r$cases; deaths <- r$deaths; N <- r$N; blk <- r$blk; B <- max(blk)
grid <- sort(unique(signif(c(exp(seq(log(0.02), log(0.8), length.out = 18)), 0.066, 0.10), 4)))
G <- length(grid)
cat("grid (", G, "):", grid, "\n")

fit_p <- function(pp, st) tryCatch(
  sird_ga_fit_tv(cases, deaths, N, blk, "both", st, fixed = list(p = pp), Delta = 1, nsub = 20L),
  error = function(e) NULL)
warm <- function(f) c(list(lambda = f$lambda), as.list(f$theta[c("gamma", "alpha", "i0", "tau2", "tau2_d")]))
row_of <- function(f, pp, start) {
  if (is.null(f)) return(NULL)        # no fit (first point of the forward sweep, last of the backward)
  g <- f$theta[["gamma"]] + f$theta[["alpha"]]
  data.frame(p = pp, start = start, loglik = f$loglik, conv = f$convergence,
             inf_period_days = 7 / g, ifr = f$theta[["alpha"]] / g,
             attack = 1 - f$kf$x[nrow(f$kf$x), 1], i0N = f$theta[["i0"]] * N,
             tau2 = f$theta[["tau2"]], tau2_d = f$theta[["tau2_d"]],
             t(setNames(tapply(f$R_t, blk, mean), paste0("R_block", seq_len(B)))))
}
better <- function(a, b) if (is.null(a)) b else if (is.null(b)) a else if (b$loglik > a$loglik) b else a

# phase 1: generic starts
gfile <- "results/paraiba_profile_p_refined_generic.rds"
if (!file.exists(gfile)) {
  jobs <- expand.grid(i = seq_len(G), s = 1:3)
  sts <- generic_starts_tv(B)
  res <- parallel::mclapply(seq_len(nrow(jobs)), function(j)
    fit_p(grid[jobs$i[j]], sts[[jobs$s[j]]]), mc.cores = cores)
  saveRDS(list(jobs = jobs, fits = res), gfile)
  cat(format(Sys.time()), "generic phase done\n")
}
gen <- readRDS(gfile)
gen_best <- lapply(seq_len(G), function(i) {
  b <- NULL
  for (j in which(gen$jobs$i == i)) b <- better(b, gen$fits[[j]])
  b
})

# phase 2: forward and backward sweeps (in parallel)
sweep <- function(dir) {
  ord <- if (dir == "fwd") seq_len(G) else rev(seq_len(G))
  out <- vector("list", G); prev <- NULL
  for (i in ord) {
    f <- if (is.null(prev)) NULL else fit_p(grid[i], warm(prev))
    out[i] <- list(f)                 # keeps NULL entries (out[[i]] <- NULL would drop them)
    prev <- better(gen_best[[i]], f)
  }
  out
}
sfile <- c(fwd = "results/paraiba_profile_p_refined_fwd.rds", bwd = "results/paraiba_profile_p_refined_bwd.rds")
todo <- names(sfile)[!file.exists(sfile)]
if (length(todo)) {
  sw <- parallel::mclapply(todo, sweep, mc.cores = min(cores, length(todo)))
  for (k in seq_along(todo)) saveRDS(sw[[k]], sfile[[todo[k]]])
  cat(format(Sys.time()), "sweeps done\n")
}
fwd <- readRDS(sfile[["fwd"]]); bwd <- readRDS(sfile[["bwd"]])
# files written before the fix above lost their NULL entry (first of fwd, last of bwd): restore
if (length(fwd) < G) fwd <- c(list(NULL), fwd)
if (length(bwd) < G) length(bwd) <- G

# all starts, and the profile
lab <- c("generic_D5", "generic_D7", "generic_D10")
all <- do.call(rbind, c(
  lapply(seq_len(nrow(gen$jobs)), function(j) row_of(gen$fits[[j]], grid[gen$jobs$i[j]], lab[gen$jobs$s[j]])),
  lapply(seq_len(G), function(i) row_of(fwd[[i]], grid[i], "forward")),
  lapply(seq_len(G), function(i) row_of(bwd[[i]], grid[i], "backward"))))
all <- all[order(all$p, all$start), ]
prof <- do.call(rbind, lapply(split(all, all$p), function(d) {
  d <- d[is.finite(d$loglik), ]
  b <- d[which.max(d$loglik), ]
  b$ll_spread_starts <- diff(range(d$loglik)); b$n_within_0.1 <- sum(d$loglik > max(d$loglik) - 0.1)
  b
}))
rownames(prof) <- NULL

prof_interval <- function(p, ll, llmax) {
  cut <- llmax - 1.92
  inside <- ll >= cut
  if (!any(inside)) return(c(lower = NA, upper = NA))
  i1 <- min(which(inside)); i2 <- max(which(inside))
  interp <- function(a, b) { # linear in (p, ll) between grid points a (outside) and b (inside)
    p[a] + (cut - ll[a]) * (p[b] - p[a]) / (ll[b] - ll[a]) }
  lo <- if (i1 > 1) interp(i1 - 1, i1) else NA      # NA: interval reaches the grid edge
  hi <- if (i2 < length(p)) interp(i2 + 1, i2) else NA
  c(lower = lo, upper = hi, n_inside = sum(inside), contiguous = all(inside[i1:i2]))
}
llmax <- max(c(prof$loglik, fbest$loglik))
old <- readRDS("results/paraiba_tv_profile_p.rds")
ci_new <- prof_interval(prof$p, prof$loglik, llmax)
ci_old <- prof_interval(old$p, old$loglik, max(c(old$loglik, fbest$loglik)))

options(width = 220)
cat(sprintf("\nunconstrained optimum (paraiba_tv.rds): ll = %.3f, p = %.4f\n", fbest$loglik, fbest$theta[["p"]]))
cat("\nAll starts (loglik, infectious period, IFR, attack, R_t block 1 and 10, convergence)\n")
print(format(all[, c("p", "start", "loglik", "conv", "inf_period_days", "ifr", "attack", "R_block1", "R_block10")],
             digits = 5), row.names = FALSE)
cat("\nRefined profile (best of 5 starts per p)\n")
print(format(prof[, c("p", "start", "loglik", "ll_spread_starts", "n_within_0.1", "conv", "inf_period_days", "ifr",
                      "attack", "i0N", "tau2", "R_block1", "R_block10")], digits = 5), row.names = FALSE)
cat("\nOld profile (paraiba_tv_profile_p.rds) vs refined at the old grid values\n")
cmp <- merge(old[, c("p", "loglik", "inf_period", "R1")],
             prof[, c("p", "loglik", "inf_period_days", "R_block1")], by = "p", all.x = TRUE,
             suffixes = c("_old", "_new"))
print(format(cmp, digits = 5), row.names = FALSE)
cat(sprintf("\n95%% profile interval for p (cut max - 1.92 = %.3f, linear interpolation)\n", llmax - 1.92))
cat(sprintf("  refined: [%.4f, %.4f]  (%d grid values inside, contiguous = %s)\n",
            ci_new[["lower"]], ci_new[["upper"]], ci_new[["n_inside"]], as.logical(ci_new[["contiguous"]])))
cat(sprintf("  old    : [%.4f, %.4f]  (%d grid values inside)\n",
            ci_old[["lower"]], ci_old[["upper"]], ci_old[["n_inside"]]))
saveRDS(list(grid = grid, all_starts = all, profile = prof, ci_refined = ci_new, ci_old = ci_old,
             llmax = llmax, comparison_old = cmp), "results/paraiba_profile_p_refined.rds")
cat("PROFILEDONE\n")
