#' Simulate daily reported cases and deaths from the SIRD jump process
#'
#' Runs the exact Gillespie simulator [sird_gillespie()] and turns the
#' latent daily increments into surveillance counts: reported cases are
#' `Binomial(new infections, p)` plus `N(0, tau2 * new infections)`
#' measurement error, reported deaths are `Binomial(new deaths, q)` plus
#' `N(0, tau2_d * new deaths)`. Counts are rounded and truncated at zero.
#'
#' @inheritParams sird_gillespie
#' @param p,tau2 Case reporting probability and case measurement-error scale.
#' @param q,tau2_d Death reporting probability and death measurement-error
#'   scale.
#' @return A list with `latent` (the Gillespie output), `cases` and
#'   `deaths` (numeric vectors of length `Tmax / Delta`, one count per
#'   interval `(t_{k-1}, t_k]`), `new_inf` and `new_deaths` (latent
#'   increments).
#' @export
sird_ga_simulate_inc <- function(N, S0, I0, R0, lambda, gamma, alpha,
                                 p, tau2, q = 1, tau2_d = 0, Tmax, Delta = 1) {
  latent <- sird_gillespie(N, S0, I0, R0, lambda, gamma, alpha, Tmax, Delta)
  new_inf <- -diff(latent$S)
  new_deaths <- diff(latent$D)
  obs <- function(x, pr, t2) {
    o <- stats::rbinom(length(x), size = pmax(round(x), 0), prob = pr) +
      stats::rnorm(length(x), 0, sqrt(pmax(t2 * x, 0)))
    pmax(round(o), 0)
  }
  list(latent = latent, cases = obs(new_inf, p, tau2),
       deaths = obs(new_deaths, q, tau2_d),
       new_inf = new_inf, new_deaths = new_deaths)
}

#' SIRD-GA log-likelihood for daily reported cases and/or deaths
#'
#' Kalman-filter approximate log-likelihood of the SIRD-GA state-space model
#' when the observations are daily *increments* (new reported cases, new
#' reported deaths), as recorded by routine surveillance systems. Daily
#' increments are linear in the augmented state \eqn{Z_k=(X_k, X_{k-1})}, so
#' the linear-Gaussian structure is preserved. Transition quantities are
#' computed in compiled code (RK4 on the joint system for the ODE solution,
#' the resolvent and the Lyapunov equation of the state-noise covariance).
#'
#' @param theta Named vector with `lambda`, `gamma`, `alpha`, `s0`, `i0`,
#'   `r0`, `p`, `tau2`, `q`, `tau2_d`.
#' @param cases,deaths Daily counts for intervals `1, ..., n` (`NA`
#'   allowed). Either may be `NULL` when `use` excludes it.
#' @param N Population size.
#' @param use Which streams enter the likelihood: `"both"`, `"cases"` or
#'   `"deaths"`.
#' @param Delta Observation interval.
#' @param nsub RK4 steps per interval.
#' @param full If `TRUE`, return the full list (fitted one-step-ahead
#'   means and standard deviations, ODE path) instead of the scalar.
#' @param lambda_t Optional vector of transmission rates, one per
#'   observation interval; when given it overrides `theta[["lambda"]]`
#'   (piecewise-constant transmission).
#' @export
sird_ga_loglik_inc <- function(theta, cases = NULL, deaths = NULL, N,
                               use = c("both", "cases", "deaths"),
                               Delta = 1, nsub = 10L, full = FALSE, lambda_t = NULL) {
  use <- match.arg(use)
  code <- switch(use, cases = 1L, deaths = 2L, both = 3L)
  th <- as.list(theta)
  n <- max(length(cases), length(deaths))
  if (is.null(cases)) cases <- rep(NA_real_, n)
  if (is.null(deaths)) deaths <- rep(NA_real_, n)
  bad <- any(c(th$lambda, th$gamma, th$alpha, th$i0, th$s0) <= 0) ||
    th$p <= 0 || th$p > 1 || th$q <= 0 || th$q > 1 ||
    th$tau2 < 0 || th$tau2_d < 0 || (th$s0 + th$i0 + th$r0) > 1 + 1e-9
  if (bad) return(if (full) list(loglik = -Inf) else -Inf)
  out <- sird_ga_loglik_inc_cpp(as.numeric(cases), as.numeric(deaths), N, Delta,
                                c(th$s0, th$i0, th$r0), th$lambda, th$gamma, th$alpha,
                                th$p, th$tau2, th$q, th$tau2_d, use = code,
                                nsub = as.integer(nsub),
                                lambda_t = if (is.null(lambda_t)) numeric(0) else as.numeric(lambda_t))
  if (full) out else out$loglik
}

#' Fast Kalman log-likelihood for prevalence observations
#'
#' Compiled counterpart of [sird_ga_loglik()] (same model and observation
#' equation, RK4 transition quantities).
#'
#' @inheritParams sird_ga_loglik
#' @param nsub RK4 steps per interval.
#' @export
sird_ga_loglik_prev <- function(theta, y, N, Delta = 1, nsub = 10L) {
  th <- as.list(theta)
  if (any(c(th$lambda, th$gamma, th$alpha, th$tau2, th$i0, th$s0) <= 0) ||
      th$p <= 0 || th$p >= 1) return(-Inf)
  sird_ga_loglik_prev_cpp(as.numeric(y), N, Delta, c(th$s0, th$i0, th$r0),
                          th$lambda, th$gamma, th$alpha, th$p, th$tau2,
                          nsub = as.integer(nsub))$loglik
}

#' Maximum likelihood for SIRD-GA from daily cases and/or deaths
#'
#' Nelder-Mead maximisation of [sird_ga_loglik_inc()] on an unconstrained
#' scale (log rates, logit probabilities and fractions), from one or several
#' starting points; the best optimum is polished by a final BFGS run.
#' Parameters not listed in `free` are held at their values in `fixed`.
#'
#' @param cases,deaths,N,use,Delta,nsub See [sird_ga_loglik_inc()].
#' @param starts A list of named lists of starting values (one per start).
#' @param free Character vector of parameters to estimate, a subset of
#'   `c("lambda", "gamma", "alpha", "i0", "s0", "p", "tau2", "q", "tau2_d")`.
#' @param fixed Named list with values for the remaining parameters (`r0`
#'   is always derived as `1 - s0 - i0` when `s0` is free, otherwise taken
#'   from `fixed`).
#' @param p_max Upper bound for `p`.
#' @param maxit Iteration limit for each Nelder-Mead run.
#' @param hessian Logical; return the numerical Hessian (unconstrained
#'   scale) at the optimum.
#' @return A list with `theta`, `loglik`, `starts` (data frame of the
#'   optimum reached from each start), `hessian`, `free` and the
#'   transformation functions.
#' @export
sird_ga_fit_inc <- function(cases = NULL, deaths = NULL, N,
                            use = c("both", "cases", "deaths"),
                            starts,
                            free = c("lambda", "gamma", "alpha", "i0", "p", "tau2", "tau2_d"),
                            fixed = list(s0 = NULL, r0 = 0, q = 1, p = 1, tau2_d = 0),
                            p_max = 0.999, Delta = 1, nsub = 10L, maxit = 4000,
                            hessian = TRUE) {
  use <- match.arg(use)
  if (!is.list(starts[[1]])) starts <- list(starts)
  tf <- list(
    lambda = c(log, exp), gamma = c(log, exp), alpha = c(log, exp),
    tau2 = c(log, exp), tau2_d = c(log, exp),
    i0 = c(stats::qlogis, stats::plogis), s0 = c(stats::qlogis, stats::plogis),
    q = c(stats::qlogis, stats::plogis),
    p = c(function(x) stats::qlogis(x / p_max), function(v) p_max * stats::plogis(v))
  )
  build <- function(v) {
    th <- fixed
    for (j in seq_along(free)) th[[free[j]]] <- tf[[free[j]]][[2]](unname(v[j]))
    if (is.null(th$s0)) th$s0 <- 1 - th$i0 - th$r0
    else if ("s0" %in% free) th$r0 <- max(1 - th$s0 - th$i0, 0)
    if (is.null(th$r0)) th$r0 <- 0
    th
  }
  nll <- function(v) {
    th <- build(v)
    if (th$s0 + th$i0 > 1) return(1e10)
    ll <- sird_ga_loglik_inc(th, cases, deaths, N, use, Delta, nsub)
    if (!is.finite(ll)) 1e10 else -ll
  }
  runs <- lapply(starts, function(st) {
    v0 <- vapply(free, function(nm) tf[[nm]][[1]](st[[nm]]), numeric(1))
    f <- tryCatch(stats::optim(v0, nll, method = "Nelder-Mead",
                               control = list(maxit = maxit, reltol = 1e-10)),
                  error = function(e) NULL)
    if (is.null(f)) return(NULL)
    f2 <- tryCatch(stats::optim(f$par, nll, method = "Nelder-Mead",
                                control = list(maxit = maxit, reltol = 1e-12)),
                   error = function(e) f)
    if (f2$value <= f$value) f2 else f
  })
  ok <- !vapply(runs, is.null, logical(1))
  vals <- vapply(runs, function(r) if (is.null(r)) Inf else r$value, numeric(1))
  best <- runs[[which.min(vals)]]
  pol <- tryCatch(stats::optim(best$par, nll, method = "BFGS",
                               control = list(maxit = 500, reltol = 1e-12)),
                  error = function(e) best)
  if (pol$value < best$value) best <- pol
  H <- NULL
  if (hessian) H <- tryCatch(stats::optimHess(best$par, nll), error = function(e) NULL)
  th <- unlist(build(best$par))
  th <- th[c("lambda", "gamma", "alpha", "s0", "i0", "r0", "p", "tau2", "q", "tau2_d")]
  start_tab <- do.call(rbind, lapply(which(ok), function(j) {
    data.frame(start = j, loglik = -runs[[j]]$value,
               t(unlist(build(runs[[j]]$par))[c(free)]))
  }))
  list(theta = th, loglik = -best$value, par = best$par, hessian = H,
       starts = start_tab, free = free, build = build, nll = nll,
       settings = list(N = N, use = use, Delta = Delta, nsub = nsub, p_max = p_max))
}

#' Profile log-likelihood of one SIRD-GA parameter
#'
#' Maximises the log-likelihood over the remaining free parameters on a
#' grid of fixed values of `param`, warm-starting each grid point from the
#' previous optimum.
#'
#' @param fit Output of [sird_ga_fit_inc()].
#' @param cases,deaths Data used for `fit`.
#' @param param Name of the profiled parameter (one of `fit$free`).
#' @param grid Values at which to profile.
#' @param maxit Nelder-Mead iteration limit per grid point.
#' @return Data frame with `value`, `loglik` and the re-optimised parameters.
#' @export
sird_ga_profile <- function(fit, cases = NULL, deaths = NULL, param, grid,
                            maxit = 3000) {
  st <- as.list(fit$theta)
  fixed <- st
  free <- setdiff(fit$free, param)
  out <- vector("list", length(grid))
  start <- st
  for (g in seq_along(grid)) {
    fixed[[param]] <- grid[g]
    f <- sird_ga_fit_inc(cases, deaths, fit$settings$N, use = fit$settings$use,
                         starts = list(start), free = free, fixed = fixed,
                         p_max = fit$settings$p_max, Delta = fit$settings$Delta,
                         nsub = fit$settings$nsub, maxit = maxit, hessian = FALSE)
    out[[g]] <- data.frame(value = grid[g], loglik = f$loglik, t(f$theta))
    start <- as.list(f$theta)
  }
  do.call(rbind, out)
}

#' SIRD-GA with piecewise-constant transmission, fitted to cases and deaths
#'
#' Maximum likelihood for the SIRD-GA model in which the transmission rate is
#' constant within blocks of observation intervals,
#' \eqn{\lambda(t)=\lambda_b} for \eqn{t} in block \eqn{b}, while the removal
#' rates \eqn{\gamma}, \eqn{\alpha}, the initial infectious fraction, the case
#' reporting probability and the two noise scales are shared. Parameters
#' listed in `fixed` are held at the given values (e.g. a profile over `p`).
#'
#' @param cases,deaths,N,use,Delta,nsub See [sird_ga_loglik_inc()].
#' @param blocks Integer vector, one entry per observation interval, giving
#'   the block of each interval (`1, ..., B`).
#' @param start Named list of starting values for `gamma`, `alpha`, `i0`,
#'   `p`, `tau2`, `tau2_d` and a vector `lambda` of length `B`.
#' @param fixed Named list of parameters held fixed (any of `gamma`,
#'   `alpha`, `i0`, `p`, `tau2`, `tau2_d`); `q` (default 1) is always fixed.
#' @param q Death reporting probability (fixed).
#' @param p_max Upper bound for `p`.
#' @param maxit Iteration limit.
#' @return A list with `lambda` (per block), `theta`, `loglik`, `par`,
#'   `hessian`, `R_t` (per interval, \eqn{\lambda(t)s(t)/(\gamma+\alpha)}) and
#'   the full Kalman output `kf`.
#' @export
sird_ga_fit_tv <- function(cases = NULL, deaths = NULL, N, blocks,
                           use = c("both", "cases", "deaths"), start, fixed = list(),
                           q = 1, p_max = 0.999, Delta = 1, nsub = 10L, maxit = 3000) {
  use <- match.arg(use)
  B <- max(blocks)
  sc_names <- setdiff(c("gamma", "alpha", "i0", "p", "tau2", "tau2_d"), names(fixed))
  if (use == "cases") sc_names <- setdiff(sc_names, "tau2_d")
  if (use == "deaths") sc_names <- setdiff(sc_names, c("p", "tau2"))
  fwd <- list(gamma = log, alpha = log, tau2 = log, tau2_d = log, i0 = stats::qlogis,
              p = function(x) stats::qlogis(x / p_max))
  bwd <- list(gamma = exp, alpha = exp, tau2 = exp, tau2_d = exp, i0 = stats::plogis,
              p = function(v) p_max * stats::plogis(v))
  unpack <- function(v) {
    lam <- exp(v[seq_len(B)])
    th <- list(q = q, r0 = 0, p = 1, tau2 = 1, tau2_d = 1, lambda = lam[1])
    for (nm in names(fixed)) th[[nm]] <- fixed[[nm]]
    for (j in seq_along(sc_names)) th[[sc_names[j]]] <- bwd[[sc_names[j]]](unname(v[B + j]))
    th$s0 <- 1 - th$i0 - th$r0
    list(th = th, lam = lam)
  }
  nll <- function(v) {
    u <- unpack(v)
    ll <- sird_ga_loglik_inc(u$th, cases, deaths, N, use, Delta, nsub,
                             lambda_t = u$lam[blocks])
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
  kf <- sird_ga_loglik_inc(u$th, cases, deaths, N, use, Delta, nsub, full = TRUE,
                           lambda_t = u$lam[blocks])
  s_t <- kf$x[-1, 1]
  R_t <- u$lam[blocks] * s_t / (u$th$gamma + u$th$alpha)
  list(lambda = u$lam, theta = unlist(u$th[c("gamma", "alpha", "i0", "p", "tau2", "tau2_d", "q")]),
       loglik = -f$value, par = f$par, hessian = H, R_t = R_t, kf = kf, blocks = blocks,
       sc_names = sc_names, fixed = fixed, convergence = f$convergence)
}
