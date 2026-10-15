#' Parametric bootstrap for a `ga_fit`
#'
#' Simulates `B` data sets from the fitted model with [ga_simulate()], that
#' is, from the Markov jump process and the observation model assumed by
#' [ga_loglik()], refits each by [ga_fit()] started at the estimates, and
#' records a statistic of every refit. Missing values in the original data
#' are kept missing in the simulated data sets, and rates that vary between
#' blocks are simulated with their block values.
#'
#' Each replicate draws its own seed from the current random-number stream,
#' so the result is reproducible under [set.seed()] (or `seed`) and does not
#' depend on `cores`.
#'
#' @param fit A [ga_fit()] result.
#' @param B Number of bootstrap replicates.
#' @param statistic Function of a `ga_fit` returning a named numeric vector,
#'   applied to `fit` and to every refit. The default returns the estimates.
#' @param init_counts Named vector of initial counts per compartment, summing
#'   to the population size of the fit. By default the initial fractions at
#'   the estimates times the population size, rounded, with the rounding
#'   remainder added to the largest compartment.
#' @param cores Number of processes for [parallel::mclapply()], which forks;
#'   use 1 on Windows.
#' @param seed Optional seed set before the replicate seeds are drawn.
#' @param maxit Nelder-Mead iteration limit of each refit.
#' @return An object of class `ga_boot` with elements `t0` (statistic at
#'   `fit`), `t` (matrix with one row per replicate, `NA` where the
#'   simulation or the refit failed), `loglik` and `convergence` (the `optim`
#'   code of each refit) and `init_counts`.
#' @seealso [confint.ga_boot()], [ga_profile()].
#' @examples
#' \donttest{
#' sir <- ga_model(c("S", "I", "R"), list(
#'   infection = ga_event("S", "I", "lambda", by = "I"),
#'   removal   = ga_event("I", "R", "gamma")))
#' obs <- list(cases = ga_incidence("infection", prob = "p", noise = "tau2"))
#' set.seed(1)
#' sim <- ga_simulate(sir, list(lambda = 0.5, gamma = 0.2, p = 0.5, tau2 = 0.1),
#'                    c(S = 4990, I = 10, R = 0), obs, n = 60)
#' init <- function(p) c(S = 1 - p$i0, I = p$i0, R = 0)
#' fit <- ga_fit(sir, sim$data, obs, 5000, init,
#'               start = list(lambda = 0.4, gamma = 0.25, i0 = 0.005, p = 0.4, tau2 = 0.2))
#' bt <- ga_bootstrap(fit, B = 20, seed = 2)
#' confint(bt)
#' }
#' @export
ga_bootstrap <- function(fit, B = 200, statistic = stats::coef, init_counts = NULL, cores = 1L,
                         seed = NULL, maxit = 1500) {
  a <- fit$call_args; th <- fit$theta
  Y0 <- as.matrix(a$data); n <- nrow(Y0)
  if (is.null(init_counts)) {
    fr <- a$init(th)[a$model$compartments]; fr <- fr / sum(fr)
    init_counts <- round(a$N * fr)
    i <- which.max(init_counts); init_counts[i] <- init_counts[i] + a$N - sum(init_counts)
    if (any(init_counts == 0 & fr > 0))
      warning("a compartment with a positive initial fraction has no individuals after rounding: ",
              paste(names(fr)[init_counts == 0 & fr > 0], collapse = ", "))
  }
  if (abs(sum(init_counts) - a$N) > 0.5) stop("init_counts must sum to the population size of the fit")
  varying <- if (length(fit$blocks)) lapply(stats::setNames(names(fit$blocks), names(fit$blocks)),
                                            function(p) th[[p]][fit$blocks[[p]]])
  th_c <- th; for (p in names(fit$blocks)) th_c[[p]] <- th[[p]][1]
  t0 <- statistic(fit)
  if (!is.null(seed)) set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  one <- function(b) {
    set.seed(seeds[b])
    sim <- ga_simulate(a$model, th_c, init_counts, a$obs, n, a$Delta, varying = varying)
    Y <- sim$data; Y[is.na(Y0)] <- NA
    f <- ga_fit(a$model, Y, a$obs, a$N, a$init, start = th[fit$free], fixed = fit$fixed,
                unit = fit$unit, blocks = fit$blocks, Delta = a$Delta, nsub = a$nsub,
                hessian = FALSE, maxit = maxit)
    s <- statistic(f)
    if (length(s) != length(t0)) stop("statistic returned a vector of a different length")
    list(t = unname(s), loglik = f$loglik, convergence = f$convergence)
  }
  res <- if (cores > 1L) parallel::mclapply(seq_len(B), function(b) tryCatch(one(b), error = function(e) NULL),
                                            mc.cores = cores, mc.preschedule = FALSE)
         else lapply(seq_len(B), function(b) tryCatch(one(b), error = function(e) NULL))
  ok <- vapply(res, is.list, logical(1))
  tt <- matrix(NA_real_, B, length(t0), dimnames = list(NULL, names(t0)))
  ll <- cv <- rep(NA_real_, B)
  for (b in which(ok)) { tt[b, ] <- res[[b]]$t; ll[b] <- res[[b]]$loglik; cv[b] <- res[[b]]$convergence }
  structure(list(t0 = t0, t = tt, loglik = ll, convergence = cv, init_counts = init_counts, B = B),
            class = "ga_boot")
}

#' Bootstrap confidence intervals for a `ga_boot`
#'
#' Percentile intervals use the quantiles of the replicates; basic intervals
#' reflect them about the original estimate. Replicates with a missing value
#' of a quantity are dropped for that quantity.
#'
#' @param object A [ga_bootstrap()] result.
#' @param parm Names of the quantities; all by default.
#' @param level Confidence level.
#' @param type `"percentile"` or `"basic"`.
#' @param ... Unused.
#' @export
confint.ga_boot <- function(object, parm, level = 0.95, type = c("percentile", "basic"), ...) {
  type <- match.arg(type)
  if (missing(parm)) parm <- names(object$t0)
  pr <- c(1 - level, 1 + level) / 2
  q <- t(vapply(parm, function(p) stats::quantile(object$t[, p], pr, na.rm = TRUE, names = FALSE), numeric(2)))
  if (type == "basic") q <- 2 * object$t0[parm] - q[, 2:1, drop = FALSE]
  dimnames(q) <- list(parm, sprintf("%.1f %%", 100 * pr))
  q
}

#' @export
print.ga_boot <- function(x, ...) {
  ok <- stats::complete.cases(x$t)
  cat("Parametric bootstrap of a Gaussian-approximation fit:", x$B, "replicates,",
      sum(!ok), "failed,", sum(x$convergence[ok] != 0), "with a nonzero convergence code\n")
  ci <- confint(x)
  print(data.frame(estimate = signif(x$t0, 4), boot_mean = signif(colMeans(x$t, na.rm = TRUE), 4),
                   boot_sd = signif(apply(x$t, 2, stats::sd, na.rm = TRUE), 4),
                   lower = signif(ci[, 1], 4), upper = signif(ci[, 2], 4)))
  invisible(x)
}
