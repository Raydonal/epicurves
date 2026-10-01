#' Logistic growth curve
#' @param t Time.
#' @param y0 Cumulative count at the start of the fitting window.
#' @param K Final epidemic size (carrying capacity).
#' @param gamma Growth rate.
#' @return Numeric vector, same length as `t`.
#' @export
logistic_curve <- function(t, y0, K, gamma) {
  K / (1 + (K / y0 - 1) * exp(-gamma * t))
}

#' Gompertz growth curve
#' @inheritParams logistic_curve
#' @export
gompertz_curve <- function(t, y0, K, gamma) {
  K * exp(-exp(-gamma * t) * log(K / y0))
}

#' Richards growth curve
#'
#' @inheritParams logistic_curve
#' @param beta Shape parameter.
#'
#' @details
#' Parametrised so that `richards_curve(0, y0, K, gamma, beta) == y0` and
#' `richards_curve(t, ...) -> K` as `t -> Inf`, consistent with the
#' Logistic and Gompertz curves in this package and with the standard
#' generalised-logistic form (Richards, 1959). Note this uses
#' `exp(-beta * gamma * t)` (decaying exponential); the thesis chapter this
#' package accompanies transcribes the equation with `exp(+beta * gamma *
#' t)`, which drives the curve to 0 rather than to `K` as `t` grows — a
#' sign slip relative to the standard form, corrected here (and in the
#' revised manuscript).
#' @export
richards_curve <- function(t, y0, K, gamma, beta) {
  K * (1 + ((K / y0)^beta - 1) * exp(-beta * gamma * t))^(-1 / beta)
}

#' Fit a single growth-curve model by nonlinear least squares
#'
#' Fits [logistic_curve()], [gompertz_curve()] or [richards_curve()] to a
#' cumulative-incidence series by nonlinear least squares
#' ([minpack.lm::nlsLM()]), which is markedly more robust to the
#' identification issues of these S-shaped curves than base [stats::nls()].
#'
#' @param t Numeric vector of time indices (`1, 2, ..., n`, relative to the
#'   start of the fitting window).
#' @param y Numeric vector of cumulative counts, same length as `t`.
#' @param model One of `"logistic"`, `"gompertz"`, `"richards"`.
#' @param beta_fixed If not `NULL` (Richards only), the shape parameter is
#'   held fixed at this value instead of estimated (improves convergence on
#'   short/fast-saturating phases, as used throughout the simulation study).
#' @param start Optional named list overriding the default starting values
#'   `y0`, `K`, `gamma` (and `beta` for Richards).
#' @param lower,upper Optional named lists overriding the default box
#'   constraints.
#' @return An object of class `growth_fit`: a list with `model`, the fitted
#'   `nls` object, `coef` (named coefficient vector), `fitted` (fitted
#'   values at `t`), `resid`, `mse` (mean squared error, Equation for EQM in
#'   the manuscript), `aic`, and `n_par` (number of estimated parameters,
#'   used by the AIC-based ensemble weights).
#' @export
fit_growth_model <- function(t, y, model = c("logistic", "gompertz", "richards"),
                              beta_fixed = 1.5, start = NULL, lower = NULL, upper = NULL) {
  model <- match.arg(model)
  n <- length(y)

  defaults <- list(
    start = list(y0 = max(y[1L], 1), K = max(y) * 1.5, gamma = 0.05, beta = beta_fixed %||% 1),
    lower = list(y0 = 1e-3, K = max(y) * 1.001, gamma = 1e-4, beta = 0.05),
    upper = list(y0 = max(y), K = max(y) * 50, gamma = 2, beta = 10)
  )
  st <- utils::modifyList(defaults$start, start %||% list())
  lo <- utils::modifyList(defaults$lower, lower %||% list())
  up <- utils::modifyList(defaults$upper, upper %||% list())

  dat <- data.frame(t = t, y = y)

  if (model == "logistic") {
    frm <- y ~ logistic_curve(t, y0, K, gamma)
    pnames <- c("y0", "K", "gamma")
  } else if (model == "gompertz") {
    frm <- y ~ gompertz_curve(t, y0, K, gamma)
    pnames <- c("y0", "K", "gamma")
  } else {
    if (!is.null(beta_fixed)) {
      frm <- stats::as.formula(sprintf("y ~ richards_curve(t, y0, K, gamma, %.10f)", beta_fixed))
      pnames <- c("y0", "K", "gamma")
    } else {
      frm <- y ~ richards_curve(t, y0, K, gamma, beta)
      pnames <- c("y0", "K", "gamma", "beta")
    }
  }

  # Growth-curve objectives are prone to flat, boundary-hugging local optima
  # when the optimiser is started far from the true curvature (common when
  # this function is called hundreds of times inside bootstrap/grid-search
  # loops on short, noisy windows). We therefore try a small grid of
  # (K, gamma) starting multipliers and keep the converged fit with the
  # lowest residual sum of squares, rather than trusting a single start.
  K_mult <- c(1.1, 1.5, 3)
  gamma0 <- c(0.02, 0.08)
  starts_grid <- expand.grid(K_mult = K_mult, gamma0 = gamma0)

  best_fit <- NULL
  best_rss <- Inf
  for (r in seq_len(nrow(starts_grid))) {
    st_r <- st
    st_r$K <- max(y) * starts_grid$K_mult[r]
    st_r$gamma <- starts_grid$gamma0[r]
    cand <- tryCatch(
      minpack.lm::nlsLM(
        frm, data = dat,
        start = st_r[pnames], lower = unlist(lo[pnames]), upper = unlist(up[pnames]),
        control = minpack.lm::nls.lm.control(maxiter = 500)
      ),
      error = function(e) NULL
    )
    if (!is.null(cand)) {
      rss <- sum(stats::residuals(cand)^2)
      if (is.finite(rss) && rss < best_rss) {
        best_rss <- rss
        best_fit <- cand
      }
    }
  }
  if (is.null(best_fit)) stop("fit_growth_model: no starting value converged for model '", model, "'")
  fit <- best_fit

  fitted_vals <- stats::fitted(fit)
  resid <- y - fitted_vals
  mse <- mean(resid^2)
  n_par <- length(pnames) + (!is.null(beta_fixed) && model == "richards") * 0L
  sigma2 <- sum(resid^2) / n
  loglik <- -0.5 * n * (log(2 * pi) + log(sigma2) + 1)
  aic <- 2 * (n_par + 1) - 2 * loglik # +1 for the estimated residual variance

  structure(list(model = model, fit = fit, coef = stats::coef(fit),
                  beta_fixed = if (model == "richards") beta_fixed else NA,
                  fitted = fitted_vals, resid = resid, mse = mse, aic = aic,
                  n_par = n_par, t = t, y = y),
             class = "growth_fit")
}

#' @export
predict.growth_fit <- function(object, newdata = NULL, ...) {
  t_new <- if (is.null(newdata)) object$t else newdata
  cf <- as.list(object$coef)
  switch(object$model,
    logistic  = logistic_curve(t_new, cf$y0, cf$K, cf$gamma),
    gompertz  = gompertz_curve(t_new, cf$y0, cf$K, cf$gamma),
    richards  = richards_curve(t_new, cf$y0, cf$K, cf$gamma,
                                 if (is.na(object$beta_fixed)) cf$beta else object$beta_fixed)
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a
