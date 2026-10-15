#' Maximum likelihood for a compartmental model under the Gaussian approximation
#'
#' Maximises [ga_loglik()] over the parameters in `start`, on an
#' unconstrained scale (logarithms of positive parameters, logits of those
#' in `unit`), by a Nelder-Mead search (relative tolerance `1e-8`) from each
#' starting point followed by BFGS refinement of the best optimum. Parameters in `fixed` are held at their
#' values.
#'
#' @inheritParams ga_loglik
#' @param init A function of the parameter list returning the named vector
#'   of initial fractions, e.g.
#'   `function(p) c(S = 1 - p$i0, I = p$i0, R = 0, D = 0)`.
#' @param start Named list of starting values of the free parameters, or a
#'   list of such lists (several starts). A parameter listed in `blocks`
#'   takes a vector with one value per block (a scalar is recycled).
#' @param fixed Named list of parameters held fixed.
#' @param unit Names of free parameters constrained to `(0, 1)`. By
#'   default, the reporting probabilities named in `obs` and any free
#'   parameter whose name ends in `0` (initial fractions such as `i0`).
#' @param blocks Optional named list giving, for a time-varying rate
#'   parameter, an integer vector with the block (`1, ..., B`) of each
#'   observation interval; the parameter then takes one value per block.
#' @param hessian Logical; compute the numerical Hessian at the optimum.
#' @param maxit Iteration limit of each Nelder-Mead run.
#' @return An object of class `ga_fit` with elements `coef` (estimates on the
#'   natural scale, block-specific values named `name[b]`), `loglik`,
#'   `convergence` (the `optim` code of the final run), `par` and `hessian`
#'   (unconstrained scale), `vcov` (on that scale),
#'   `starts` (log-likelihood reached from each start) and the inputs needed
#'   by [ga_profile()], [confint.ga_fit()] and [predict.ga_fit()].
#' @export
ga_fit <- function(model, data, obs, N, init, start, fixed = list(), unit = NULL,
                   blocks = NULL, Delta = 1, nsub = 10L, hessian = TRUE, maxit = 1500) {
  if (!is.list(start[[1]])) start <- list(start)
  free <- names(start[[1]])
  n <- nrow(as.matrix(data))
  if (is.null(unit)) {
    pr <- unlist(lapply(obs, function(o) if (is.character(o$prob)) o$prob))
    unit <- union(intersect(pr, free), grep("0$", free, value = TRUE))
  }
  B <- vapply(free, function(p) if (p %in% names(blocks)) max(blocks[[p]]) else 1L, 1L)
  for (p in names(blocks)) stopifnot(length(blocks[[p]]) == n)
  lab <- unlist(lapply(free, function(p) if (B[[p]] > 1) sprintf("%s[%d]", p, seq_len(B[[p]])) else p))
  fwd <- function(p, v) if (p %in% unit) stats::qlogis(v) else log(v)
  bwd <- function(p, z) if (p %in% unit) stats::plogis(z) else exp(z)
  pack <- function(st) unlist(lapply(free, function(p) fwd(p, rep_len(st[[p]], B[[p]]))))
  unpack <- function(z) {
    th <- fixed; pos <- 0L
    for (p in free) { th[[p]] <- bwd(p, unname(z[pos + seq_len(B[[p]])])); pos <- pos + B[[p]] }
    th
  }
  ll_of <- function(th) {
    varying <- if (length(blocks)) lapply(stats::setNames(names(blocks), names(blocks)),
                                          function(p) th[[p]][blocks[[p]]])
    th_c <- th; for (p in names(blocks)) th_c[[p]] <- th[[p]][1]
    ga_loglik(model, th_c, init(th), data, obs, N, Delta, nsub, varying = varying)
  }
  nll <- function(z) { v <- tryCatch(ll_of(unpack(z)), error = function(e) -Inf); if (is.finite(v)) -v else 1e10 }
  runs <- lapply(start, function(st)
    stats::optim(pack(st), nll, method = "Nelder-Mead", control = list(maxit = maxit, reltol = 1e-8)))
  vals <- vapply(runs, `[[`, 1, "value")
  best <- runs[[which.min(vals)]]
  pol <- tryCatch(stats::optim(best$par, nll, method = "BFGS", control = list(maxit = 500, reltol = 1e-12)),
                  error = function(e) best)
  if (pol$value < best$value) best <- pol
  names(best$par) <- lab
  H <- if (hessian) tryCatch(stats::optimHess(best$par, nll), error = function(e) NULL)
  V <- if (!is.null(H)) tryCatch({ v <- solve(H); if (all(diag(v) > 0)) v else NULL }, error = function(e) NULL)
  if (!is.null(V)) dimnames(V) <- list(lab, lab)
  th <- unpack(best$par)
  coef <- unlist(lapply(free, function(p) th[[p]])); names(coef) <- lab
  structure(list(coef = coef, loglik = -best$value, convergence = best$convergence,
                 par = best$par, hessian = H, vcov = V,
                 starts = data.frame(start = seq_along(runs), loglik = -vals),
                 theta = th, free = free, unit = unit, blocks = blocks, fixed = fixed,
                 call_args = list(model = model, data = data, obs = obs, N = N, init = init,
                                  Delta = Delta, nsub = nsub),
                 ll_of = ll_of, unpack = unpack, nll = nll),
            class = "ga_fit")
}

#' @export
print.ga_fit <- function(x, ...) {
  cat("Gaussian-approximation fit of a", length(x$call_args$model$compartments), "compartment model\n")
  cat("log-likelihood: ", format(x$loglik, digits = 8), " (", nrow(x$starts), " start(s))\n", sep = "")
  ci <- confint(x)
  print(data.frame(estimate = signif(x$coef, 4), lower = signif(ci[, 1], 4), upper = signif(ci[, 2], 4)))
  invisible(x)
}

#' @export
coef.ga_fit <- function(object, ...) object$coef

#' @export
logLik.ga_fit <- function(object, ...) {
  structure(object$loglik, df = length(object$par), class = "logLik")
}

#' Wald confidence intervals for a `ga_fit`
#'
#' Intervals are computed on the unconstrained scale of the optimisation and
#' transformed back, so they respect the parameter bounds.
#' @param object A [ga_fit()] result.
#' @param parm Ignored (all parameters are returned).
#' @param level Confidence level.
#' @param ... Unused.
#' @export
confint.ga_fit <- function(object, parm, level = 0.95, ...) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  se <- if (!is.null(object$vcov)) sqrt(diag(object$vcov)) else rep(NA_real_, length(object$par))
  base <- sub("\\[.*$", "", names(object$par))
  tr <- function(p, v) if (p %in% object$unit) stats::plogis(v) else exp(v)
  out <- cbind(mapply(tr, base, object$par - z * se), mapply(tr, base, object$par + z * se))
  dimnames(out) <- list(names(object$par), sprintf("%.1f %%", 100 * c(1 - level, 1 + level) / 2))
  out
}

#' One-step-ahead predictions from a `ga_fit`
#'
#' @param object A [ga_fit()] result.
#' @param ... Unused.
#' @return A list with `fitted` and `sd` (matrices of predicted means and
#'   standard deviations of the counts, one column per stream) and `x` (path
#'   of the mean, as fractions).
#' @export
predict.ga_fit <- function(object, ...) {
  a <- object$call_args; th <- object$theta
  varying <- if (length(object$blocks)) lapply(stats::setNames(names(object$blocks), names(object$blocks)),
                                               function(p) th[[p]][object$blocks[[p]]])
  th_c <- th; for (p in names(object$blocks)) th_c[[p]] <- th[[p]][1]
  ga_loglik(a$model, th_c, a$init(th), a$data, a$obs, a$N, a$Delta, a$nsub, varying = varying, full = TRUE)
}

#' Profile log-likelihood for one parameter of a `ga_fit`
#'
#' @param fit A [ga_fit()] result.
#' @param param Name of a free, scalar parameter.
#' @param grid Values at which the profile is computed. Each point is
#'   warm-started from the previous one, so an ordered grid that starts near
#'   the estimate works best.
#' @param maxit Nelder-Mead iteration limit per point.
#' @return A data frame with the grid value, the profile log-likelihood and
#'   the other estimates.
#' @export
ga_profile <- function(fit, param, grid, maxit = 3000) {
  a <- fit$call_args
  st <- fit$theta[setdiff(fit$free, param)]
  out <- vector("list", length(grid))
  for (g in seq_along(grid)) {
    fx <- fit$fixed; fx[[param]] <- grid[g]
    f <- ga_fit(a$model, a$data, a$obs, a$N, a$init, start = st, fixed = fx, unit = fit$unit,
                blocks = fit$blocks, Delta = a$Delta, nsub = a$nsub, hessian = FALSE, maxit = maxit)
    out[[g]] <- data.frame(value = grid[g], loglik = f$loglik, t(f$coef), check.names = FALSE)
    st <- f$theta[setdiff(fit$free, param)]
  }
  do.call(rbind, out)
}

#' Plot the one-step-ahead predictions of a `ga_fit`
#'
#' One panel per observation stream: observed counts (points), predicted
#' means (line) and predicted means plus and minus two standard deviations
#' (band), all from the Kalman filter at the estimates.
#'
#' @param x A [ga_fit()] result.
#' @param labels Optional panel labels, one per stream.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
plot.ga_fit <- function(x, labels = NULL, ...) {
  pr <- predict(x)
  Y <- as.matrix(x$call_args$data)
  if (is.null(labels)) labels <- vapply(x$call_args$obs, function(o) paste(o$type, "of", o$target), "")
  df <- do.call(rbind, lapply(seq_len(ncol(Y)), function(j)
    data.frame(time = seq_len(nrow(Y)) * x$call_args$Delta, observed = Y[, j], fitted = pr$fitted[, j],
               lower = pmax(pr$fitted[, j] - 2 * pr$sd[, j], 0), upper = pr$fitted[, j] + 2 * pr$sd[, j],
               stream = factor(labels[j], levels = labels))))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), fill = "grey85") +
    ggplot2::geom_point(ggplot2::aes(y = .data$observed), color = epicurves_palette[["observed"]], size = 1) +
    ggplot2::geom_line(ggplot2::aes(y = .data$fitted), color = epicurves_palette[["ensemble"]], linewidth = 0.7) +
    ggplot2::facet_wrap(~ stream, scales = "free_y") +
    ggplot2::labs(x = "Time", y = "Count") +
    ggplot2::theme_classic(base_size = 12)
}
