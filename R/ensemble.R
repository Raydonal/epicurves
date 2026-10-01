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
#' @param ... Unused; present for S3 consistency with [stats::predict()].
#' @export
predict.ensemble_fit <- function(object, newdata = NULL, ...) {
  t_new <- newdata %||% object$t
  preds <- vapply(object$models, function(m) predict(m, newdata = t_new), numeric(length(t_new)))
  as.vector(preds %*% object$weights)
}

#' Bootstrap confidence band for an ensemble fit
#'
#' Two resampling schemes are supported. `method = "poisson"` implements
#' the parametric bootstrap of the manuscript's simulation study: bootstrap
#' pseudo-observations are drawn as \eqn{y_k^{(b)} \sim
#' \mathrm{Poisson}(f_{ens}(t_k))}, appropriate when the data-generating
#' process is itself (approximately) Poisson, as in the simulated growth
#' curves used there. `method = "residual"` instead resamples the
#' standardised fitted residuals with replacement and adds them back onto
#' the fitted curve, which is appropriate for real surveillance counts: the
#' Poisson scheme assumes \eqn{\mathrm{Var}(y_k)=f_{ens}(t_k)}, which is
#' routinely violated by over-dispersed daily death/case counts (reporting
#' batches, weekday effects), producing artificially narrow bands (see the
#' real-data application in the manuscript, where this is checked
#' directly). In both cases the full ensemble (all three growth models,
#' re-weighted) is refit to each replicate, and percentile confidence
#' limits are taken across the `B` replicates at every requested time
#' point.
#'
#' @param fit An `ensemble_fit` object.
#' @param newdata Numeric vector of times at which to build the band
#'   (defaults to `fit$t`; include out-of-window times for a forecast
#'   band).
#' @param method `"poisson"` (default) or `"residual"`.
#' @param B Number of bootstrap replicates.
#' @param level Confidence level.
#' @param beta_fixed Passed to [fit_ensemble()] for each replicate.
#' @return A list with `boot_curves` (`B x length(newdata)` matrix) and
#'   `ci` (data frame with `time`, `lower`, `upper`, `fit`).
#' @export
ensemble_bootstrap_ci <- function(fit, newdata = NULL, method = c("poisson", "residual"),
                                   B = 500L, level = 0.95, beta_fixed = 1.5) {
  method <- match.arg(method)
  newdata <- newdata %||% fit$t
  boot_curves <- matrix(NA_real_, B, length(newdata))
  resid <- fit$y - fit$fitted

  for (b in seq_len(B)) {
    y_rep <- switch(method,
      poisson  = stats::rpois(length(fit$t), lambda = pmax(fit$fitted, 1e-6)),
      residual = pmax(fit$fitted + sample(resid, length(resid), replace = TRUE), 0)
    )
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
  list(boot_curves = boot_curves, ci = ci, B = B, level = level, method = method)
}

#' Bootstrap confidence intervals for each growth model's parameters
#'
#' Resamples the ensemble fit (see [ensemble_bootstrap_ci()] for the choice
#' between `"poisson"` and `"residual"` resampling) and, instead of (or in
#' addition to) the ensemble curve, collects each of the three growth
#' models' own parameter estimates (\eqn{Y_0}, \eqn{\gamma}, \eqn{K}) and
#' weight on every replicate, returning percentile confidence limits for
#' each. This is what the manuscript's per-segment parameter tables report.
#'
#' @inheritParams ensemble_bootstrap_ci
#' @return A data frame with one row per model x parameter, giving the
#'   point estimate, lower/upper confidence limits, and number of
#'   successful replicates (`B` minus failed refits).
#' @export
ensemble_bootstrap_params <- function(fit, method = c("poisson", "residual"),
                                       B = 200L, level = 0.95, beta_fixed = 1.5) {
  method <- match.arg(method)
  resid <- fit$y - fit$fitted
  model_names <- names(fit$models)
  par_names <- c("y0", "gamma", "K")

  boot_par <- array(NA_real_, dim = c(B, length(model_names), length(par_names)),
                     dimnames = list(NULL, model_names, par_names))
  boot_w <- matrix(NA_real_, B, length(model_names), dimnames = list(NULL, model_names))

  for (b in seq_len(B)) {
    y_rep <- switch(method,
      poisson  = stats::rpois(length(fit$t), lambda = pmax(fit$fitted, 1e-6)),
      residual = pmax(fit$fitted + sample(resid, length(resid), replace = TRUE), 0)
    )
    rep_fit <- tryCatch(fit_ensemble(fit$t, y_rep, criterion = fit$criterion,
                                       beta_fixed = beta_fixed),
                          error = function(e) NULL)
    if (!is.null(rep_fit)) {
      for (m in model_names) {
        cf <- rep_fit$models[[m]]$coef
        boot_par[b, m, ] <- c(cf[["y0"]], cf[["gamma"]], cf[["K"]])
      }
      boot_w[b, ] <- rep_fit$weights[model_names]
    }
  }

  alpha2 <- (1 - level) / 2
  out <- do.call(rbind, lapply(model_names, function(m) {
    do.call(rbind, lapply(par_names, function(p) {
      vals <- boot_par[, m, p]
      lims <- stats::quantile(vals, probs = c(alpha2, 1 - alpha2), na.rm = TRUE)
      data.frame(model = m, parameter = p,
                 estimate = fit$models[[m]]$coef[[p]],
                 lower = lims[1L], upper = lims[2L],
                 n_ok = sum(!is.na(vals)))
    }))
  }))
  w_rows <- data.frame(model = model_names, parameter = "weight",
                        estimate = fit$weights[model_names],
                        lower = apply(boot_w, 2, stats::quantile, probs = alpha2, na.rm = TRUE),
                        upper = apply(boot_w, 2, stats::quantile, probs = 1 - alpha2, na.rm = TRUE),
                        n_ok = colSums(!is.na(boot_w)))
  rownames(out) <- NULL
  rownames(w_rows) <- NULL
  out <- rbind(out, w_rows)
  out$method <- method
  out
}
