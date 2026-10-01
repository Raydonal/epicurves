#' Maximum (approximate) likelihood estimation for the SIRD-GA model
#'
#' Estimates `theta = (lambda, gamma, alpha, i0, p, tau2)` (and optionally
#' `r0`) by Nelder-Mead maximisation of the Kalman-filter approximate
#' log-likelihood ([sird_ga_loglik()]), following the estimation strategy
#' described in the SIRD-GA manuscript. Positivity and unit-interval
#' constraints are enforced through a log/logit reparametrisation, so the
#' optimiser itself runs unconstrained.
#'
#' @param y Numeric vector of observed infectious counts, `y_0, ..., y_n`.
#' @param N Population size.
#' @param Delta Time step between observations.
#' @param start Named numeric vector/list with starting values for
#'   `lambda`, `gamma`, `alpha`, `i0`, `p`, `tau2` (and `r0` if
#'   `estimate_r0 = TRUE`).
#' @param r0_fixed Value at which `r0` is held fixed when
#'   `estimate_r0 = FALSE` (default `0`, appropriate for data starting at
#'   or near the onset of the epidemic).
#' @param estimate_r0 Logical; estimate `r0` jointly instead of fixing it.
#' @param p_max Upper bound enforced on the reporting probability `p`
#'   (default `0.999`). Real daily case series are rarely reported with
#'   near-100% completeness; without an upper cap strictly below 1, the
#'   optimiser can drive `p -> 1` and `tau2` to a degenerate extreme
#'   (perfect reporting, all residual variance absorbed by the
#'   measurement-noise term), a boundary artefact rather than a genuine
#'   optimum. Setting `p_max` around 0.95-0.98 is recommended for real
#'   data; simulation studies where the reporting rate is itself the
#'   target of recovery can use `p_max = 0.999`.
#' @param n_sub,Q0_diag Passed to [sird_ga_loglik()].
#' @param control List of control parameters passed to [stats::optim()].
#' @param hessian Logical; compute the numerical Hessian in the
#'   unconstrained parametrisation (used only for a quick asymptotic
#'   diagnostic — see [sird_ga_bootstrap()] for the inferential confidence
#'   intervals reported in the manuscript).
#'
#' @return A list with elements `theta` (named vector of MLEs in natural
#'   scale), `loglik`, `convergence`, `optim_fit` (raw [stats::optim()]
#'   output in the unconstrained parametrisation), and `free_se` (delta-free
#'   standard errors on the unconstrained scale, or `NULL`).
#' @export
sird_ga_fit <- function(y, N, Delta = 1, start, r0_fixed = 0,
                         estimate_r0 = FALSE, p_max = 0.999,
                         n_sub = 10L, Q0_diag = 1e-8,
                         control = list(maxit = 3000, reltol = 1e-10),
                         hessian = TRUE) {
  start <- as.list(start)

  to_free <- function(th) {
    v <- c(log(th$lambda), log(th$gamma), log(th$alpha),
           stats::qlogis(th$i0), stats::qlogis(th$p / p_max), log(th$tau2))
    if (estimate_r0) {
      r0_frac <- th$r0 / (1 - th$i0)
      v <- c(v, stats::qlogis(max(min(r0_frac, 1 - 1e-8), 1e-8)))
    }
    v
  }
  from_free <- function(v) {
    lambda <- exp(v[1]); gamma <- exp(v[2]); alpha <- exp(v[3])
    i0 <- stats::plogis(v[4]); p <- p_max * stats::plogis(v[5]); tau2 <- exp(v[6])
    r0 <- if (estimate_r0) (1 - i0) * stats::plogis(v[7]) else r0_fixed
    s0 <- 1 - i0 - r0
    c(lambda = lambda, gamma = gamma, alpha = alpha,
      s0 = s0, i0 = i0, r0 = r0, p = p, tau2 = tau2)
  }

  negloglik <- function(v) {
    theta <- from_free(v)
    ll <- sird_ga_loglik(theta, y = y, N = N, Delta = Delta,
                          n_sub = n_sub, Q0_diag = Q0_diag)
    if (!is.finite(ll)) return(1e10)
    -ll
  }

  v0 <- to_free(start)
  fit <- stats::optim(v0, negloglik, method = "Nelder-Mead",
                       control = control, hessian = hessian)

  theta_hat <- from_free(fit$par)
  free_se <- NULL
  if (hessian && !is.null(fit$hessian)) {
    free_se <- tryCatch(sqrt(diag(solve(fit$hessian))), error = function(e) NULL)
  }

  list(theta = theta_hat, loglik = -fit$value, convergence = fit$convergence,
       optim_fit = fit, free_se = free_se,
       to_free = to_free, from_free = from_free,
       settings = list(N = N, Delta = Delta, n_sub = n_sub, Q0_diag = Q0_diag))
}

#' Parametric bootstrap confidence intervals for the SIRD-GA MLE
#'
#' Two bootstrap schemes are supported. `method = "gillespie"` simulates `B`
#' independent epidemic trajectories and surveillance records from the
#' fitted parameters via the exact jump-process simulator
#' ([sird_ga_simulate_obs()]) and re-estimates the model on each replicate —
#' appropriate when extinction risk is low. `method = "residual"` resamples
#' the standardised one-step-ahead Kalman residuals of the fitted model with
#' replacement, reconstructs bootstrap observation series around the fitted
#' predictive mean, and re-estimates the model on each — appropriate for
#' real-data settings near the epidemic threshold, where re-simulating the
#' jump process can produce early extinction.
#'
#' @param fit Object returned by [sird_ga_fit()].
#' @param y Original observed series used to produce `fit`.
#' @param N,Delta Passed through to [sird_ga_loglik()] / the simulator.
#' @param method `"gillespie"` or `"residual"`.
#' @param B Number of bootstrap replicates.
#' @param level Confidence level.
#' @param S0,I0,R0 Initial head counts, required for `method = "gillespie"`
#'   (derived from `fit$theta` and `N` if omitted).
#' @return A list with `boot_theta` (`B x p` matrix of bootstrap estimates)
#'   and `ci` (matrix of percentile confidence limits, one row per
#'   parameter).
#' @export
sird_ga_bootstrap <- function(fit, y, N, Delta = 1,
                               method = c("gillespie", "residual"),
                               B = 200L, level = 0.95,
                               S0 = NULL, I0 = NULL, R0 = NULL) {
  method <- match.arg(method)
  theta <- fit$theta
  settings <- fit$settings

  refit_one <- function(y_rep, start) {
    tryCatch({
      f <- sird_ga_fit(y_rep, N = N, Delta = Delta, start = start,
                        n_sub = settings$n_sub, Q0_diag = settings$Q0_diag,
                        hessian = FALSE)
      f$theta
    }, error = function(e) rep(NA_real_, length(theta)))
  }

  boot_theta <- matrix(NA_real_, B, length(theta))
  colnames(boot_theta) <- names(theta)

  if (method == "gillespie") {
    if (is.null(S0)) S0 <- round(N * theta[["s0"]])
    if (is.null(I0)) I0 <- round(N * theta[["i0"]])
    if (is.null(R0)) R0 <- round(N * theta[["r0"]])
    Tmax <- (length(y) - 1L) * Delta
    for (b in seq_len(B)) {
      sim <- sird_ga_simulate_obs(N, S0, I0, R0,
                                   theta[["lambda"]], theta[["gamma"]], theta[["alpha"]],
                                   theta[["p"]], theta[["tau2"]], Tmax, Delta)
      boot_theta[b, ] <- refit_one(sim$y, as.list(theta))
    }
  } else {
    kf <- sird_ga_loglik(theta, y, N, Delta, settings$n_sub, settings$Q0_diag,
                          keep_states = TRUE)
    M <- theta[["p"]] * kf$predicted[, "i"]
    resid <- y - M
    resid_c <- resid - mean(resid)
    for (b in seq_len(B)) {
      y_rep <- pmax(M + sample(resid_c, length(resid_c), replace = TRUE), 0)
      boot_theta[b, ] <- refit_one(y_rep, as.list(theta))
    }
  }

  alpha2 <- (1 - level) / 2
  ci <- t(apply(boot_theta, 2, stats::quantile,
                probs = c(alpha2, 1 - alpha2), na.rm = TRUE))
  colnames(ci) <- c("lower", "upper")
  list(boot_theta = boot_theta, ci = ci, method = method, B = B, level = level)
}

#' Multi-start SIRD-GA fitting, with a practical-identifiability summary
#'
#' Fits the SIRD-GA model from several starting points spanning a range of
#' assumed mean infectious periods, and reports all converged optima sorted
#' by log-likelihood. This is useful because, on real surveillance data, the
#' SIRD-GA log-likelihood can have a near-flat ridge along which
#' `lambda` trades off against `gamma + alpha` while `R0` and the shape of
#' the fitted curve stay essentially unchanged — a known weak-identifiability
#' feature of SIR/SIRD models fit to incidence-only series (the removal
#' rate is identified mainly through its sum with the transmission rate and
#' the curve's inflection timing, not independently). Reporting the single
#' best-fitting optimum without checking for this ridge can be misleading;
#' this function makes the ridge (if present) visible.
#'
#' @inheritParams sird_ga_fit
#' @param r0_guess Starting guess for \eqn{R_0} (used to scale `lambda`
#'   across starts, together with each candidate infectious period).
#' @param infectious_period_days Numeric vector of candidate mean
#'   infectious periods \eqn{1/(\gamma+\alpha)} (in days) to start from
#'   (default spans 2-10 days, covering typical respiratory-infection
#'   priors).
#' @param alpha_frac Fraction of the total removal rate `gamma + alpha`
#'   assigned to `alpha` at each starting point (default 0.1).
#' @param p0,tau2_0,i0_0 Starting values shared across all starts for `p`,
#'   `tau2` and `i0`.
#' @param ... Further arguments passed to [sird_ga_fit()].
#' @return A list with `results` (a list of [sird_ga_fit()] outputs, one per
#'   start, in the same order as `infectious_period_days`), `summary` (a
#'   data frame with the starting infectious period, final log-likelihood,
#'   `R0` and key parameter estimates for each start, sorted best-first) and
#'   `best` (the single best-fitting result).
#' @export
sird_ga_fit_multistart <- function(y, N, Delta = 1, r0_guess = 1.1,
                                    infectious_period_days = c(2, 3, 4, 5, 7, 10),
                                    alpha_frac = 0.1, p0 = 0.9, tau2_0 = 1,
                                    i0_0 = 1e-4, ...) {
  results <- vector("list", length(infectious_period_days))
  for (j in seq_along(infectious_period_days)) {
    removal <- 1 / infectious_period_days[j]
    start <- list(lambda = r0_guess * removal, gamma = (1 - alpha_frac) * removal,
                   alpha = alpha_frac * removal, i0 = i0_0, p = p0, tau2 = tau2_0)
    results[[j]] <- tryCatch(sird_ga_fit(y, N, Delta = Delta, start = start, hessian = FALSE, ...),
                              error = function(e) NULL)
  }

  summary_df <- do.call(rbind, lapply(seq_along(results), function(j) {
    r <- results[[j]]
    if (is.null(r)) {
      return(data.frame(infectious_period_start = infectious_period_days[j],
                          loglik = NA_real_, R0 = NA_real_,
                          lambda = NA_real_, gamma = NA_real_, alpha = NA_real_,
                          p = NA_real_, tau2 = NA_real_))
    }
    data.frame(infectious_period_start = infectious_period_days[j], loglik = r$loglik,
               R0 = sird_R0(r$theta[["lambda"]], r$theta[["gamma"]], r$theta[["alpha"]]),
               lambda = r$theta[["lambda"]], gamma = r$theta[["gamma"]],
               alpha = r$theta[["alpha"]], p = r$theta[["p"]], tau2 = r$theta[["tau2"]])
  }))
  summary_df <- summary_df[order(-summary_df$loglik), ]

  best_idx <- which.max(vapply(results, function(r) if (is.null(r)) -Inf else r$loglik, numeric(1)))
  list(results = results, summary = summary_df, best = results[[best_idx]])
}
