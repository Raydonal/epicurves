# Shared helpers for the application-regime studies (paraiba_blocks.R, paraiba_floor.R,
# sim_paraiba_regime.R). Sourced, not run on its own.
#
# - generic_starts_tv(): the three generic starts of paraiba_tv.R (mean infectious period
#   5, 7, 10 days; R0 1.2; IFR 0.02; p 0.2; tau2 20; tau2_d 2; i0 5e-6).
# - sird_ga_fit_tv_floor(): copy of epicurves::sird_ga_fit_tv() whose log-likelihood calls the
#   compiled filter with an explicit `var_floor` (lower bound, in persons, on the expected
#   increment that enters the observation variance R_k). With var_floor = 1e-3 it is the
#   package function.
# - wald_derived(): Wald-type intervals for the derived quantities of Table 3, computed as in
#   paraiba_report.R (4000 draws from N(par, H^{-1}) on the unconstrained scale, mapped through
#   the model, 2.5% and 97.5% quantiles), plus Wald intervals for the unconstrained parameters.
suppressMessages(library(epicurves))

generic_starts_tv <- function(B) {
  lapply(c(5, 7, 10), function(D) {
    rem <- 7 / D
    list(lambda = rep(1.2 * rem, B), gamma = 0.98 * rem, alpha = 0.02 * rem,
         i0 = 5e-6, p = 0.2, tau2 = 20, tau2_d = 2)
  })
}

loglik_floor <- function(th, cases, deaths, N, Delta, nsub, lambda_t, var_floor, full = FALSE) {
  bad <- any(c(th$lambda, th$gamma, th$alpha, th$i0, th$s0) <= 0) ||
    th$p <= 0 || th$p > 1 || th$q <= 0 || th$q > 1 ||
    th$tau2 < 0 || th$tau2_d < 0 || (th$s0 + th$i0 + th$r0) > 1 + 1e-9
  if (bad) return(if (full) list(loglik = -Inf) else -Inf)
  out <- epicurves:::sird_ga_loglik_inc_cpp(as.numeric(cases), as.numeric(deaths), N, Delta,
                                            c(th$s0, th$i0, th$r0), th$lambda, th$gamma, th$alpha,
                                            th$p, th$tau2, th$q, th$tau2_d, use = 3L,
                                            nsub = as.integer(nsub), var_floor = var_floor,
                                            lambda_t = as.numeric(lambda_t))
  if (full) out else out$loglik
}

# Same optimisation sequence as epicurves::sird_ga_fit_tv (BFGS, Nelder-Mead, BFGS; optimHess),
# use = "both", no fixed parameters.
sird_ga_fit_tv_floor <- function(cases, deaths, N, blocks, start, var_floor = 1e-3,
                                 q = 1, p_max = 0.999, Delta = 1, nsub = 20L, maxit = 3000) {
  B <- max(blocks)
  sc_names <- c("gamma", "alpha", "i0", "p", "tau2", "tau2_d")
  fwd <- list(gamma = log, alpha = log, tau2 = log, tau2_d = log, i0 = stats::qlogis,
              p = function(x) stats::qlogis(x / p_max))
  bwd <- list(gamma = exp, alpha = exp, tau2 = exp, tau2_d = exp, i0 = stats::plogis,
              p = function(v) p_max * stats::plogis(v))
  unpack <- function(v) {
    lam <- exp(v[seq_len(B)])
    th <- list(q = q, r0 = 0, lambda = lam[1])
    for (j in seq_along(sc_names)) th[[sc_names[j]]] <- bwd[[sc_names[j]]](unname(v[B + j]))
    th$s0 <- 1 - th$i0 - th$r0
    list(th = th, lam = lam)
  }
  nll <- function(v) {
    u <- unpack(v)
    ll <- loglik_floor(u$th, cases, deaths, N, Delta, nsub, u$lam[blocks], var_floor)
    if (!is.finite(ll)) 1e10 else -ll
  }
  v0 <- c(log(start$lambda), vapply(sc_names, function(nm) fwd[[nm]](start[[nm]]), 1))
  f <- stats::optim(v0, nll, method = "BFGS", control = list(maxit = maxit, reltol = 1e-12))
  f2 <- stats::optim(f$par, nll, method = "Nelder-Mead", control = list(maxit = maxit, reltol = 1e-12))
  if (f2$value < f$value) f <- f2
  f3 <- stats::optim(f$par, nll, method = "BFGS", control = list(maxit = maxit, reltol = 1e-12))
  if (f3$value < f$value) f <- f3
  u <- unpack(f$par)
  H <- tryCatch(stats::optimHess(f$par, nll), error = function(e) NULL)
  kf <- loglik_floor(u$th, cases, deaths, N, Delta, nsub, u$lam[blocks], var_floor, full = TRUE)
  s_t <- kf$x[-1, 1]
  R_t <- u$lam[blocks] * s_t / (u$th$gamma + u$th$alpha)
  list(lambda = u$lam, theta = unlist(u$th[c("gamma", "alpha", "i0", "p", "tau2", "tau2_d", "q")]),
       loglik = -f$value, par = f$par, hessian = H, R_t = R_t, kf = kf, blocks = blocks,
       sc_names = sc_names, fixed = list(), convergence = f$convergence, var_floor = var_floor)
}

# Derived quantities of Table 3 from an unconstrained parameter vector (as in paraiba_report.R).
derived_from_par <- function(v, sc_names, cases, deaths, N, blk, q = 1, p_max = 0.999) {
  B <- max(blk)
  th <- list(q = q, r0 = 0, lambda = exp(v[1]))
  bw <- list(gamma = exp, alpha = exp, tau2 = exp, tau2_d = exp, i0 = plogis,
             p = function(z) p_max * plogis(z))
  for (j in seq_along(sc_names)) th[[sc_names[j]]] <- bw[[sc_names[j]]](unname(v[B + j]))
  th$s0 <- 1 - th$i0 - th$r0
  lam <- exp(v[seq_len(B)])
  kf <- sird_ga_loglik_inc(th, cases, deaths, N, "both", 1, 20L, full = TRUE,
                           lambda_t = lam[blk])
  s <- kf$x[-1, 1]
  Rt <- tapply(lam[blk] * s / (th$gamma + th$alpha), blk, mean)
  c(inf_period_days = 7 / (th$gamma + th$alpha), ifr = th$alpha / (th$gamma + th$alpha),
    p = th$p, i0N = th$i0 * N, tau2 = th$tau2, tau2_d = th$tau2_d,
    attack = 1 - s[length(s)], setNames(Rt, paste0("R_block", seq_len(B))))
}

# Point estimates and nominal 95% intervals. Derived quantities: quantiles of `ndraw` draws
# (paraiba_report.R), and also estimate +- 1.96 * SD(draws). Unconstrained: par +- 1.96 * SE.
wald_derived <- function(f, cases, deaths, N, blk, ndraw = 4000, seed = 1) {
  est <- derived_from_par(f$par, f$sc_names, cases, deaths, N, blk)
  H <- f$hessian
  ev <- if (is.null(H) || any(!is.finite(H))) NA_real_ else
    min(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  na <- rep(NA_real_, length(est)); names(na) <- names(est)
  nap <- rep(NA_real_, length(f$par))
  out <- list(est = est, lower = na, upper = na, sd = na, min_eig = ev, pd = isTRUE(ev > 0),
              par = f$par, se_par = nap)
  if (!out$pd) return(out)
  V <- tryCatch(solve(H), error = function(e) NULL)
  L <- if (is.null(V)) NULL else tryCatch(chol((V + t(V)) / 2), error = function(e) NULL)
  if (is.null(L)) { out$pd <- FALSE; return(out) }
  out$se_par <- sqrt(diag(V))
  set.seed(seed)
  draws <- t(replicate(ndraw, derived_from_par(f$par + drop(rnorm(length(f$par)) %*% L),
                                               f$sc_names, cases, deaths, N, blk)))
  out$lower <- apply(draws, 2, quantile, 0.025, na.rm = TRUE)
  out$upper <- apply(draws, 2, quantile, 0.975, na.rm = TRUE)
  out$sd <- apply(draws, 2, sd, na.rm = TRUE)
  out
}

# Weekly Paraíba data, as in paraiba_tv.R.
paraiba_weekly <- function() {
  data(paraiba_covid19, package = "epicurves", envir = environment())
  data(paraiba_population, package = "epicurves", envir = environment())
  d <- paraiba_covid19
  wk <- d$epi_week; keep <- wk %in% names(which(table(wk) == 7))
  list(cases = as.numeric(tapply(d$new_cases[keep], wk[keep], sum)),
       deaths = as.numeric(tapply(d$new_deaths[keep], wk[keep], sum)),
       weeks = as.numeric(names(table(wk[keep]))), N = paraiba_population)
}
