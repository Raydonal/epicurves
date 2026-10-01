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
#'   partition point (default 0.5, i.e. the split must at least halve the
#'   MSE relative to the unsplit fit).
#' @return An object of class `piecewise_fit`: a list with `cuts` (sorted
#'   confirmed partition points), `segments` (list of `ensemble_fit`
#'   objects, one per final segment, each carrying the original time
#'   range in `$t_range`), `fitted` (full fitted curve aligned with `y`),
#'   and `history` (list of grid-search diagnostics at every update step,
#'   suitable for reproducing the manuscript's MSE-profile figures).
#' @export
fit_piecewise_ensemble <- function(y, t = NULL, window = 90L, update_every = 20L,
                                    min_obs = 12L, criterion = "mse",
                                    beta_fixed = 1.5, threshold = 0.5) {
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

      if (!is.na(search$mse_split) && search$mse_split <= (1 - threshold) * mse_nosplit) {
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
#' @export
predict.piecewise_fit <- function(object, newdata, ...) {
  vapply(newdata, function(tt) {
    s <- min(sum(tt > object$cuts) + 1L, length(object$segments))
    seg <- object$segments[[s]]
    t_rel <- tt - seg$t_range[1L] + 1L
    predict(seg, newdata = t_rel)
  }, numeric(1))
}
