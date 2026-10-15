# Accuracy of the transition quantities and of the log-likelihood as functions of the number of
# integration steps per observation interval.
# Reference: joint RK4 integration of (x, Phi, V) with 4000 steps per interval.
# Error of Q_k (and A_k): max over k of ||Q_k - Q_ref||_F / ||Q_ref||_F (Frobenius norm).
# (a) Scenario 1 rates, daily intervals, 100 days, i0 = 10/20000: RK4 and Simpson quadrature
#     (sird_ga_Qk, lsoda with rtol = atol = 1e-6, on the substep grid).
# (b) Paraíba fitted rates on the weekly scale (first-block lambda), Delta = 1 week, 40 weeks.
# (c) Change of the maximised log-likelihood of the Paraíba fit and of a Scenario 1 data set.
suppressMessages(library(epicurves))
slice <- function(M, k) if (length(dim(M)) == 3) M[, , k] else matrix(M[, k], 3)
relerr <- function(get, ref, n) max(sapply(1:n, function(k) {
  R <- slice(ref, k); norm(get(k) - R, "F") / norm(R, "F") }))

cat("(a) Scenario 1, daily intervals\n")
la <- 0.6; ga <- 0.3; al <- 0.1; i0 <- 10 / 20000; x0 <- c(1 - i0, i0, 0); n <- 100
ref <- sird_ga_transition_cpp(n, 1, x0, la, ga, al, 4000)
for (ns in c(2, 5, 10, 20, 40)) {
  m <- sird_ga_transition_cpp(n, 1, x0, la, ga, al, ns)
  cat(sprintf("  RK4     steps=%3d  Q rel.err=%.2e  A rel.err=%.2e\n", ns,
              relerr(function(k) slice(m$M, k), ref$M, n), relerr(function(k) slice(m$A, k), ref$A, n)))
}
for (ns in c(10, 20, 40, 80)) {
  tr <- sird_ga_transition(0:n, x0, list(lambda = la, gamma = ga, alpha = al), n_sub = ns)
  cat(sprintf("  Simpson steps=%3d  Q rel.err=%.2e\n", ns,
              relerr(function(k) sird_ga_Qk(tr, k, 1), ref$M, n)))
}

cat("(b) Paraiba fitted rates, weekly intervals\n")
r <- readRDS("results/paraiba_tv.rds"); f <- r$res[[which.max(sapply(r$res, `[[`, "loglik"))]]
th <- f$theta; i0 <- th[["i0"]]; x0 <- c(1 - i0, i0, 0); n <- length(r$cases)
ref <- sird_ga_transition_cpp(n, 1, x0, f$lambda[1], th[["gamma"]], th[["alpha"]], 4000)
for (ns in c(2, 5, 10, 20, 40)) {
  m <- sird_ga_transition_cpp(n, 1, x0, f$lambda[1], th[["gamma"]], th[["alpha"]], ns)
  cat(sprintf("  RK4     steps=%3d  Q rel.err=%.2e  A rel.err=%.2e\n", ns,
              relerr(function(k) slice(m$M, k), ref$M, n), relerr(function(k) slice(m$A, k), ref$A, n)))
}

cat("(c) log-likelihood at the fitted parameters versus 4000 steps\n")
thl <- as.list(th); thl$s0 <- 1 - thl$i0; thl$r0 <- 0; thl$lambda <- f$lambda[1]
llp <- function(ns) sird_ga_loglik_inc(thl, r$cases, r$deaths, r$N, "both", 1, ns,
                                       lambda_t = f$lambda[r$blk])
l_ref <- llp(4000)
for (ns in c(2, 5, 10, 20, 40)) cat(sprintf("  Paraiba     steps=%3d  ll - ll_ref = %.2e\n", ns, llp(ns) - l_ref))
set.seed(1)
si <- sird_ga_simulate_inc(20000, 20000 - 10, 10, 0, 0.6, 0.3, 0.1, p = 0.5, tau2 = 0.16, q = 1,
                           tau2_d = 0.1, Tmax = 100)
th1 <- list(lambda = 0.6, gamma = 0.3, alpha = 0.1, s0 = 1 - 10 / 20000, i0 = 10 / 20000, r0 = 0,
            p = 0.5, tau2 = 0.16, q = 1, tau2_d = 0.1)
lls <- function(ns) sird_ga_loglik_inc(th1, si$cases, si$deaths, 20000, "both", 1, ns)
l_ref <- lls(4000)
for (ns in c(2, 5, 10, 20, 40)) cat(sprintf("  Scenario 1  steps=%3d  ll - ll_ref = %.2e\n", ns, lls(ns) - l_ref))
