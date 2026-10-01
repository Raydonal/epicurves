#' Inverse-criterion ensemble weights
#'
#' Computes weights proportional to the inverse of each model's goodness-of-
#' fit criterion (mean squared error or AIC), so that better-fitting models
#' receive more weight, following Chowell et al. (2021) as adapted in the
#' manuscript.
#'
#' @param crit Named numeric vector of per-model criterion values (MSE or
#'   AIC; lower is better for both).
#' @return Named numeric vector of weights summing to 1.
#' @export
ensemble_weights <- function(crit) {
  inv <- 1 / crit
  inv / sum(inv)
}

#' Fit the full ensemble of growth-curve models to one data window
#'
#' Fits Logistic, Gompertz and Richards curves ([fit_growth_model()]) to a
#' single window of cumulative counts, computes MSE- and AIC-based weights,
#' and returns the ensemble (weighted-average) curve together with the
#' individual fits.
#'
#' @param t,y As in [fit_growth_model()].
#' @param criterion `"mse"` (default, matching the main text) or `"aic"`
#'   (matching the manuscript's comparative analysis).
#' @param beta_fixed Shape parameter passed to the Richards fit (see
#'   [fit_growth_model()]).
#' @return An object of class `ensemble_fit`: a list with `models` (named
#'   list of three `growth_fit` objects), `weights`, `criterion`,
#'   `fitted` (ensemble fitted values at `t`), `mse` (ensemble MSE) and
#'   `t`, `y`.
#' @export
fit_ensemble <- function(t, y, criterion = c("mse", "aic"), beta_fixed = 1.5) {
  criterion <- match.arg(criterion)
  models <- list(
    gompertz = fit_growth_model(t, y, "gompertz"),
    richards = fit_growth_model(t, y, "richards", beta_fixed = beta_fixed),
    logistic = fit_growth_model(t, y, "logistic")
  )
  crit_vals <- vapply(models, function(m) if (criterion == "mse") m$mse else m$aic, numeric(1))
  w <- ensemble_weights(crit_vals)

  fitted_mat <- vapply(models, function(m) m$fitted, numeric(length(t)))
  fitted_ens <- as.vector(fitted_mat %*% w)
  mse_ens <- mean((y - fitted_ens)^2)

  structure(list(models = models, weights = w, criterion = criterion,
                  crit_vals = crit_vals, fitted = fitted_ens, mse = mse_ens,
                  t = t, y = y),
             class = "ensemble_fit")
}

#' Predict the ensemble curve at arbitrary times
#' @param object An `ensemble_fit` object.
#' @param newdata Numeric vector of times (defaults to `object$t`).
#' @export
predict.ensemble_fit <- function(object, newdata = NULL, ...) {
  t_new <- newdata %||% object$t
  preds <- vapply(object$models, function(m) predict(m, newdata = t_new), numeric(length(t_new)))
  as.vector(preds %*% object$weights)
}

#' Parametric (Poisson) bootstrap confidence band for an ensemble fit
#'
#' Implements the bootstrap procedure of the manuscript: bootstrap
#' pseudo-observations are drawn as \eqn{y_k^{(b)} \sim
#' \mathrm{Poisson}(f_{ens}(t_k))}, the full ensemble (all three growth
#' models, re-weighted) is refit to each replicate, and the fitted/
#' forecast curve is recomputed. Percentile confidence limits are taken
#' across the `B` replicates at every requested time point.
#'
#' @param fit An `ensemble_fit` object.
#' @param newdata Numeric vector of times at which to build the band
#'   (defaults to `fit$t`; include out-of-window times for a forecast
#'   band).
#' @param B Number of bootstrap replicates.
#' @param level Confidence level.
#' @param beta_fixed Passed to [fit_ensemble()] for each replicate.
#' @return A list with `boot_curves` (`B x length(newdata)` matrix) and
#'   `ci` (data frame with `time`, `lower`, `upper`, `fit`).
#' @export
ensemble_bootstrap_ci <- function(fit, newdata = NULL, B = 500L, level = 0.95,
                                   beta_fixed = 1.5) {
  newdata <- newdata %||% fit$t
  boot_curves <- matrix(NA_real_, B, length(newdata))

  for (b in seq_len(B)) {
    y_rep <- stats::rpois(length(fit$t), lambda = pmax(fit$fitted, 1e-6))
    rep_fit <- tryCatch(fit_ensemble(fit$t, y_rep, criterion = fit$criterion,
                                       beta_fixed = beta_fixed),
                          error = function(e) NULL)
    if (!is.null(rep_fit)) {
      boot_curves[b, ] <- tryCatch(predict(rep_fit, newdata = newdata),
                                     error = function(e) NA_real_)
    }
  }

  alpha2 <- (1 - level) / 2
  lims <- apply(boot_curves, 2, stats::quantile, probs = c(alpha2, 1 - alpha2), na.rm = TRUE)
  ci <- data.frame(time = newdata, lower = lims[1L, ], upper = lims[2L, ],
                    fit = predict(fit, newdata = newdata))
  list(boot_curves = boot_curves, ci = ci, B = B, level = level)
}
