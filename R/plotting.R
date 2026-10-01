#' Okabe-Ito colorblind-safe palette used throughout this package's figures
#' @keywords internal
epicurves_palette <- c(
  gompertz = "#E69F00", richards = "#0072B2", logistic = "#D55E00",
  ensemble = "#000000", observed = "#009E73", fitted = "#CC79A7"
)

#' Plot an ensemble fit against the observed data
#'
#' @param fit An `ensemble_fit` object.
#' @param show_models Logical; overlay the three individual growth-curve
#'   fits in addition to the ensemble curve.
#' @return A `ggplot` object.
#' @export
plot_ensemble_fit <- function(fit, show_models = TRUE) {
  df <- data.frame(t = fit$t, y = fit$y, ensemble = fit$fitted)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$t)) +
    ggplot2::geom_point(ggplot2::aes(y = .data$y), color = epicurves_palette[["observed"]],
                         size = 1.3, alpha = 0.7)

  if (show_models) {
    for (nm in names(fit$models)) {
      df[[nm]] <- fit$models[[nm]]$fitted
      p <- p + ggplot2::geom_line(ggplot2::aes_string(y = nm), data = df,
                                   color = epicurves_palette[[nm]], linewidth = 0.6, alpha = 0.8)
    }
  }
  p + ggplot2::geom_line(ggplot2::aes(y = .data$ensemble), color = epicurves_palette[["ensemble"]],
                          linewidth = 1.1) +
    ggplot2::labs(x = "Time (days)", y = "Cumulative count") +
    ggplot2::theme_classic(base_size = 13)
}

#' Plot an MSE (or AIC) profile from a partition-point grid search
#'
#' Reproduces the manuscript's EQM-profile figures.
#'
#' @param search Object returned by [find_partition_point()].
#' @param true_d Optional true partition point (vertical reference line).
#' @export
plot_partition_profile <- function(search, true_d = NULL) {
  df <- data.frame(d = search$grid, mse = search$mse_profile)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$d, y = .data$mse)) +
    ggplot2::geom_line(color = epicurves_palette[["ensemble"]], linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = search$d_hat, color = epicurves_palette[["richards"]],
                         linetype = "dashed", linewidth = 0.7) +
    ggplot2::labs(x = "Candidate partition point d", y = "Pooled MSE") +
    ggplot2::theme_classic(base_size = 13)
  if (!is.null(true_d)) {
    p <- p + ggplot2::geom_vline(xintercept = true_d, color = epicurves_palette[["gompertz"]],
                                  linetype = "dotted", linewidth = 0.7)
  }
  p
}

#' Plot a fitted piecewise ensemble curve, with confirmed partition points
#'
#' @param fit A `piecewise_fit` object.
#' @param ci Optional data frame with columns `time`, `lower`, `upper`
#'   (e.g. built by calling [ensemble_bootstrap_ci()] per segment and
#'   row-binding the results) to draw as a shaded confidence band.
#' @param forecast_from Optional time beyond which a dashed vertical line
#'   marks the start of an out-of-sample forecast window.
#' @export
plot_piecewise_fit <- function(fit, ci = NULL, forecast_from = NULL) {
  df <- data.frame(t = fit$t, y = fit$y, fitted = fit$fitted)
  p <- ggplot2::ggplot()
  if (!is.null(ci)) {
    p <- p + ggplot2::geom_ribbon(data = ci, ggplot2::aes(x = .data$time, ymin = .data$lower,
                                                           ymax = .data$upper),
                                   fill = epicurves_palette[["richards"]], alpha = 0.2)
  }
  p <- p +
    ggplot2::geom_point(data = df, ggplot2::aes(x = .data$t, y = .data$y),
                         color = epicurves_palette[["observed"]], size = 1, alpha = 0.6) +
    ggplot2::geom_line(data = df, ggplot2::aes(x = .data$t, y = .data$fitted),
                        color = epicurves_palette[["ensemble"]], linewidth = 1) +
    ggplot2::geom_vline(xintercept = fit$cuts, color = epicurves_palette[["gompertz"]],
                         linetype = "dashed", linewidth = 0.5)
  if (!is.null(forecast_from)) {
    p <- p + ggplot2::geom_vline(xintercept = forecast_from, color = "grey50", linetype = "dotted")
  }
  p + ggplot2::labs(x = "Time (days)", y = "Cumulative count") +
    ggplot2::theme_classic(base_size = 13)
}

#' Plot observed vs. fitted SIRD-GA incidence, with an optional bootstrap band
#'
#' @param y Observed counts.
#' @param kf Output of [sird_ga_loglik()] called with `keep_states = TRUE`.
#' @param ci Optional data frame (`time`, `lower`, `upper`) on the
#'   head-count scale, e.g. built from [sird_ga_bootstrap()] replicate
#'   trajectories.
#' @param cumulative Logical; plot cumulative (default `FALSE`, i.e. daily)
#'   incidence.
#' @export
plot_sird_ga_fit <- function(y, kf, ci = NULL, cumulative = FALSE) {
  yy <- if (cumulative) cumsum(y) else y
  ff <- if (cumulative) cumsum(kf$fitted) else kf$fitted
  df <- data.frame(t = kf$times, observed = yy, fitted = ff)

  p <- ggplot2::ggplot()
  if (!is.null(ci)) {
    p <- p + ggplot2::geom_ribbon(data = ci, ggplot2::aes(x = .data$time, ymin = .data$lower,
                                                           ymax = .data$upper),
                                   fill = epicurves_palette[["richards"]], alpha = 0.2)
  }
  p + ggplot2::geom_line(data = df, ggplot2::aes(x = .data$t, y = .data$observed,
                                                   color = "Observed"), linewidth = 0.7) +
    ggplot2::geom_line(data = df, ggplot2::aes(x = .data$t, y = .data$fitted,
                                                 color = "SIRD-GA"), linewidth = 0.9) +
    ggplot2::scale_color_manual(name = "", values = c(Observed = epicurves_palette[["observed"]],
                                                        `SIRD-GA` = epicurves_palette[["ensemble"]])) +
    ggplot2::labs(x = "Time (days)", y = if (cumulative) "Cumulative cases" else "Daily cases") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.85, 0.85))
}
