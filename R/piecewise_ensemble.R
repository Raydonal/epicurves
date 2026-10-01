#' Mean squared error of a split at a single candidate partition point
#'
#' Fits the growth-curve ensemble separately to the two sides of a
#' candidate split `d` and returns the pooled MSE, as used by the grid
#' search in [find_partition_point()].
#'
#' @param t,y Time index and cumulative counts of the (open) window being
#'   searched, both starting at 1.
#' @param d Candidate partition point, `1 <= d < length(y)`.
#' @param criterion,beta_fixed Passed to [fit_ensemble()].
#' @return A list with `mse` (pooled MSE, `NA` if either side fails to
#'   converge) and `fit1`, `fit2` (the two `ensemble_fit` objects, or
#'   `NULL` on failure).
#' @export
split_mse <- function(t, y, d, criterion = "mse", beta_fixed = 1.5) {
  n <- length(y)
  idx1 <- seq_len(d); idx2 <- (d + 1L):n
  fit1 <- tryCatch(fit_ensemble(t[idx1], y[idx1], criterion, beta_fixed), error = function(e) NULL)
  fit2 <- tryCatch(fit_ensemble(t[idx2] - d, y[idx2], criterion, beta_fixed), error = function(e) NULL)
  if (is.null(fit1) || is.null(fit2)) return(list(mse = NA_real_, fit1 = NULL, fit2 = NULL))
  sse <- sum((y[idx1] - fit1$fitted)^2) + sum((y[idx2] - fit2$fitted)^2)
  list(mse = sse / n, fit1 = fit1, fit2 = fit2)
}

#' Grid search for a single structural change (partition) point
#'
#' Scans a grid of candidate partition points with increments of size 1,
#' as described in the manuscript, and selects the point minimising the
#' pooled mean squared error, subject to a minimum number of observations
#' on each side.
#'
#' @param t,y Time index and cumulative counts of the open window.
#' @param min_obs Minimum number of observations required on each side of
#'   the candidate split (default 12, as in Simulation Scenarios 1 and 3 of
#'   the manuscript; use 15 to match Scenario 2).
#' @param criterion,beta_fixed Passed to [fit_ensemble()].
#' @return A list with `grid` (candidate points searched), `mse_profile`
#'   (MSE at each grid point), `d_hat` (selected partition point), and
#'   `mse_split` (MSE at `d_hat`).
#' @export
find_partition_point <- function(t, y, min_obs = 12L, criterion = "mse", beta_fixed = 1.5) {
  n <- length(y)
  grid <- seq.int(min_obs, n - min_obs)
  if (length(grid) < 1L) {
    return(list(grid = integer(0), mse_profile = numeric(0), d_hat = NA_integer_, mse_split = NA_real_))
  }
  mse_profile <- vapply(grid, function(d) split_mse(t, y, d, criterion, beta_fixed)$mse, numeric(1))
  best <- which.min(mse_profile)
  list(grid = grid, mse_profile = mse_profile,
       d_hat = grid[best], mse_split = mse_profile[best])
}

#' Sequential piecewise ensemble model
#'
#' Reproduces the sequential partition-point identification algorithm of
#' the manuscript: starting from an initial data window, candidate
#' partition points are searched by grid search
#' ([find_partition_point()]); the window is then advanced by
#' `update_every` observations and, at each update, the most recent
#' (still-unsplit) segment is re-examined for a new structural change,
#' which is confirmed only if the pooled MSE of the best split is at least
#' `threshold` lower than the MSE of leaving that segment unsplit.
#'
#' @param y Numeric vector of cumulative counts (the full observed series).
#' @param t Optional numeric time index (defaults to `seq_along(y)`).
#' @param window Size of the initial data window used to search for the
#'   first partition point (default 90, as in the manuscript).
#' @param update_every Number of observations by which the window advances
#'   at each sequential update (default 20).
#' @param min_obs Minimum admissible segment length (default 12; use 15 to
#'   match Simulation Scenario 2 of the manuscript).
#' @param criterion,beta_fixed Passed to [fit_ensemble()].
#' @param threshold Relative MSE reduction required to confirm a new
#'   partition point under `confirm_criterion = "mse_threshold"` (default
#'   0.5, i.e. the split must at least halve the MSE relative to the
#'   unsplit fit — this is the rule described in the manuscript).
#' @param confirm_criterion How a candidate split is confirmed.
#'   `"mse_threshold"` (default) reproduces the manuscript's rule exactly.
#'   `"aic"` instead compares a changepoint-penalised AIC for the split
#'   against the unsplit fit — splitting doubles the number of growth-curve
#'   parameters (three models re-fit on each side) and adds one degree of
#'   freedom for the partition-point location itself, and `"aic"` is the
#'   more conservative choice when a flat MSE-reduction threshold
#'   over-segments noisy real series (see the package vignette for a
#'   worked comparison on both simulated and real surveillance data).
#'   `"permutation"` instead performs a formal hypothesis test
#'   ([permutation_test_partition()]): the candidate split is confirmed
#'   only if its MSE reduction is unlikely under a null, unsplit growth
#'   process, at level `alpha`, properly accounting for the multiple
#'   comparisons implicit in scanning the whole candidate grid. This is
#'   the criterion we recommend when statistical significance, not just a
#'   descriptive rule, is required (e.g. for direct comparison against
#'   nonparametric permutation-based detectors such as
#'   \code{ecp::e.divisive}); it is markedly more expensive
#'   computationally than the other two.
#' @param alpha Significance level used only when
#'   `confirm_criterion = "permutation"`.
#' @param B Number of permutation-test null replicates, used only when
#'   `confirm_criterion = "permutation"`; passed to
#'   [permutation_test_partition()].
#' @return An object of class `piecewise_fit`: a list with `cuts` (sorted
#'   confirmed partition points), `segments` (list of `ensemble_fit`
#'   objects, one per final segment, each carrying the original time
#'   range in `$t_range`), `fitted` (full fitted curve aligned with `y`),
#'   and `history` (list of grid-search diagnostics at every update step,
#'   suitable for reproducing the manuscript's MSE-profile figures).
#' @export
fit_piecewise_ensemble <- function(y, t = NULL, window = 90L, update_every = 20L,
                                    min_obs = 12L, criterion = "mse",
                                    beta_fixed = 1.5, threshold = 0.5,
                                    confirm_criterion = c("mse_threshold", "aic", "permutation"),
                                    alpha = 0.05, B = 99L) {
  confirm_criterion <- match.arg(confirm_criterion)
  n_total <- length(y)
  t <- t %||% seq_len(n_total)
  cuts <- integer(0)
  history <- list()
  end <- min(window, n_total)

  repeat {
    open_start <- if (length(cuts) == 0L) 1L else cuts[length(cuts)] + 1L
    open_idx <- open_start:end
    t_open <- t[open_idx] - t[open_start] + 1L
    y_open <- y[open_idx]
    n_open <- length(open_idx)

    step_log <- list(open_range = c(open_start, end))

    if (n_open >= 2L * min_obs) {
      search <- find_partition_point(t_open, y_open, min_obs, criterion, beta_fixed)
      nosplit <- tryCatch(fit_ensemble(t_open, y_open, criterion, beta_fixed), error = function(e) NULL)
      mse_nosplit <- if (is.null(nosplit)) Inf else nosplit$mse
      step_log$search <- search
      step_log$mse_nosplit <- mse_nosplit

      confirmed <- FALSE
      if (!is.na(search$mse_split)) {
        if (confirm_criterion == "mse_threshold") {
          confirmed <- search$mse_split <= (1 - threshold) * mse_nosplit
        } else if (confirm_criterion == "aic") {
          # BIC-style comparison: splitting re-fits all three growth models
          # on each side (doubling their parameters) and adds one degree of
          # freedom for the partition-point location itself. The BIC
          # penalty (k * log(n), rather than AIC's 2k) is deliberately
          # stricter, since segments here are short (12-90 observations)
          # relative to the 9 growth-curve parameters being fit, a regime
          # where AIC under-penalises and a flat MSE-reduction threshold
          # over-segments (see the package vignette).
          k_single <- sum(vapply(nosplit$models, `[[`, numeric(1), "n_par")) + 1L
          bic_nosplit <- n_open * log(mse_nosplit) + k_single * log(n_open)
          bic_split <- n_open * log(search$mse_split) + 2 * k_single * log(n_open) + log(n_open)
          confirmed <- is.finite(bic_split) && bic_split < bic_nosplit
        } else {
          pt <- permutation_test_partition(t_open, y_open, min_obs, B = B, beta_fixed = beta_fixed)
          step_log$permutation_test <- pt
          confirmed <- pt$p_value < alpha
        }
      }
      if (confirmed) {
        new_cut <- open_start - 1L + search$d_hat
        cuts <- c(cuts, new_cut)
        step_log$new_cut <- new_cut
      }
    }
    history[[length(history) + 1L]] <- step_log

    if (end >= n_total) break
    end <- min(end + update_every, n_total)
  }

  bounds <- c(0L, cuts, n_total)
  segments <- vector("list", length(bounds) - 1L)
  fitted <- rep(NA_real_, n_total)
  for (s in seq_along(segments)) {
    idx <- (bounds[s] + 1L):bounds[s + 1L]
    t_rel <- t[idx] - t[idx[1L]] + 1L
    ef <- fit_ensemble(t_rel, y[idx], criterion, beta_fixed)
    ef$t_range <- range(t[idx])
    segments[[s]] <- ef
    fitted[idx] <- ef$fitted
  }

  structure(list(cuts = cuts, segments = segments, fitted = fitted,
                  t = t, y = y, history = history,
                  window = window, update_every = update_every,
                  min_obs = min_obs, criterion = criterion, threshold = threshold),
             class = "piecewise_fit")
}

#' Predict from a piecewise ensemble fit
#'
#' @param object A `piecewise_fit` object.
#' @param newdata Numeric vector of times (original time scale); times
#'   beyond the last confirmed segment's final observation are forecast
#'   using that segment's ensemble fit.
#' @param ... Unused; present for S3 consistency with [stats::predict()].
#' @export
predict.piecewise_fit <- function(object, newdata, ...) {
  vapply(newdata, function(tt) {
    s <- min(sum(tt > object$cuts) + 1L, length(object$segments))
    seg <- object$segments[[s]]
    t_rel <- tt - seg$t_range[1L] + 1L
    predict(seg, newdata = t_rel)
  }, numeric(1))
}

#' Best single-model split MSE at a candidate partition point (fast proxy)
#'
#' A single-growth-model analogue of [split_mse()], roughly three times
#' cheaper (one model fit per side instead of three), used internally by
#' [permutation_test_partition()] to generate the null distribution of the
#' grid-maximum MSE-reduction statistic at a cost that keeps a full
#' permutation test tractable. The model used is `null_model` (default
#' `"gompertz"`), fixed across the whole null simulation, not re-selected
#' per replicate.
#' @inheritParams split_mse
#' @param null_model One of `"logistic"`, `"gompertz"`, `"richards"`.
#' @keywords internal
split_mse_single <- function(t, y, d, null_model = "gompertz", beta_fixed = 1.5) {
  n <- length(y)
  idx1 <- seq_len(d); idx2 <- (d + 1L):n
  f1 <- tryCatch(fit_growth_model(t[idx1], y[idx1], null_model, beta_fixed = beta_fixed),
                  error = function(e) NULL)
  f2 <- tryCatch(fit_growth_model(t[idx2] - d, y[idx2], null_model, beta_fixed = beta_fixed),
                  error = function(e) NULL)
  if (is.null(f1) || is.null(f2)) return(NA_real_)
  (sum((y[idx1] - f1$fitted)^2) + sum((y[idx2] - f2$fitted)^2)) / n
}

#' Permutation test for a single candidate structural change point
#'
#' Tests whether the observed reduction in pooled MSE from the best
#' candidate split in [find_partition_point()] is larger than what a
#' single homogeneous (unsplit) growth process would produce by chance,
#' properly accounting for the multiple comparisons implicit in scanning
#' the whole candidate grid: the null distribution is built from the
#' *grid-maximum* MSE reduction in each of `B` parametric-bootstrap
#' replicates simulated from the unsplit null fit, following the same
#' logic as the permutation test of \code{ecp::e.divisive}
#' (Matteson & James, 2013) but specialised to parametric growth curves
#' rather than a general multivariate distributional distance. For
#' computational tractability, the null replicates' grid search uses a
#' single growth model ([split_mse_single()]) rather than the full
#' three-model ensemble used for the reported fit; this is an
#' approximation, disclosed here rather than left implicit.
#'
#' @param t,y Time index and cumulative counts of the open window.
#' @param min_obs Minimum segment length (as in [find_partition_point()]).
#' @param B Number of parametric-bootstrap null replicates.
#' @param null_model Growth model used for the (cheaper) null-replicate
#'   grid searches; defaults to whichever of the three models best fits
#'   the unsplit window.
#' @param beta_fixed Passed to the growth-curve fits.
#' @param use_ensemble_null Logical. If `TRUE`, the null replicates are
#'   generated from, and re-fit with, the full three-model ensemble
#'   (matching the observed statistic exactly, at roughly three times the
#'   cost); if `FALSE` (default), the cheaper single-model proxy of
#'   [split_mse_single()] is used. We found the single-model proxy
#'   under-powered in simulation (it can fail to confirm genuine change
#'   points detected by the full ensemble); `use_ensemble_null = TRUE` is
#'   recommended whenever the extra cost is affordable.
#' @return A list with `p_value`, `observed_reduction` (MSE_nosplit minus
#'   MSE at the best split), `null_reductions` (the `B` simulated null
#'   maxima) and `d_hat` (the best split, as in [find_partition_point()]).
#' @export
permutation_test_partition <- function(t, y, min_obs = 12L, B = 99L,
                                        null_model = NULL, beta_fixed = 1.5,
                                        use_ensemble_null = TRUE) {
  n <- length(y)
  nosplit <- fit_ensemble(t, y, "mse", beta_fixed)
  search <- find_partition_point(t, y, min_obs, "mse", beta_fixed)
  observed_reduction <- nosplit$mse - search$mse_split
  grid <- seq.int(min_obs, n - min_obs)

  if (use_ensemble_null) {
    # Null replicates simulated from, and re-fit with, the same full
    # ensemble used for the observed statistic: an exact (if ~3x costlier)
    # analogue of the single-model proxy below.
    null_reductions <- vapply(seq_len(B), function(b) {
      y_star <- stats::rpois(n, lambda = pmax(nosplit$fitted, 1e-6))
      nosplit_star <- tryCatch(fit_ensemble(t, y_star, "mse", beta_fixed), error = function(e) NULL)
      if (is.null(nosplit_star)) return(NA_real_)
      profile <- vapply(grid, function(d) split_mse(t, y_star, d, "mse", beta_fixed)$mse, numeric(1))
      if (all(is.na(profile))) return(NA_real_)
      nosplit_star$mse - min(profile, na.rm = TRUE)
    }, numeric(1))
    null_model <- "ensemble"
  } else {
    if (is.null(null_model)) {
      null_model <- names(nosplit$models)[which.min(vapply(nosplit$models, `[[`, numeric(1), "mse"))]
    }
    null_fit <- fit_growth_model(t, y, null_model, beta_fixed = beta_fixed)
    null_reductions <- vapply(seq_len(B), function(b) {
      y_star <- stats::rpois(n, lambda = pmax(null_fit$fitted, 1e-6))
      # Refit the null model to y_star itself (not to the original data) so
      # that the "no-split" reference absorbs the same refitting advantage
      # the split profile gets; otherwise min(profile) is compared against
      # the true generating variance rather than an equally-overfit
      # baseline, which makes any split look spuriously better.
      nosplit_star <- tryCatch(fit_growth_model(t, y_star, null_model, beta_fixed = beta_fixed),
                                 error = function(e) NULL)
      if (is.null(nosplit_star)) return(NA_real_)
      profile <- vapply(grid, function(d) split_mse_single(t, y_star, d, null_model, beta_fixed), numeric(1))
      if (all(is.na(profile))) return(NA_real_)
      nosplit_star$mse - min(profile, na.rm = TRUE)
    }, numeric(1))
  }

  p_value <- (sum(null_reductions >= observed_reduction, na.rm = TRUE) + 1) / (sum(!is.na(null_reductions)) + 1)
  list(p_value = p_value, observed_reduction = observed_reduction,
       null_reductions = null_reductions, d_hat = search$d_hat, null_model = null_model)
}
