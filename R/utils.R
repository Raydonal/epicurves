#' Centered moving-average smoother
#'
#' A simple centered rolling mean, used to attenuate day-of-week reporting
#' artefacts (lower weekend counts, Monday catch-up spikes) in daily
#' COVID-19 surveillance series before SIRD-GA estimation. Edge windows are
#' truncated (not padded), so the smoothed series has the same length as
#' the input.
#'
#' @param x Numeric vector.
#' @param window Window width in days (default 7). Must be odd for a
#'   symmetric window; even values are incremented by 1 with a warning.
#' @return A numeric vector of the same length as `x`.
#' @export
smooth_series <- function(x, window = 7L) {
  if (window %% 2L == 0L) {
    window <- window + 1L
    warning("`window` must be odd for a centered average; using ", window)
  }
  half <- (window - 1L) %/% 2L
  n <- length(x)
  vapply(seq_len(n), function(i) {
    lo <- max(1L, i - half); hi <- min(n, i + half)
    mean(x[lo:hi])
  }, numeric(1))
}
