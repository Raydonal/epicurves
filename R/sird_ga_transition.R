#' Solve the SIRD-GA state-space transition quantities
#'
#' Jointly solves the deterministic SIRD trajectory \eqn{x(t)} and the
#' fundamental (sensitivity) matrix of the linearised dynamics in a
#' **single** [deSolve::ode()] call spanning the whole observation window,
#' using a *reset event* that sets the fundamental matrix back to the
#' identity at the start of every observation interval. This gives the
#' same locally-anchored resolvent \eqn{\Phi(s, t_{k-1})} as restarting the
#' solver once per interval (and the same numerical stability — the matrix
#' never has to span more than one interval, which matters for long
#' series), while paying the fixed per-call overhead of [deSolve::ode()]
#' only once instead of `n` times. This is the single biggest cost in
#' repeated likelihood evaluation (MLE optimisation, bootstrap), so the
#' single-call form is the default implementation.
#'
#' @param times Numeric vector of observation times, `t_0 < t_1 < ... < t_n`
#'   (typically `0:n`).
#' @param x0 Initial state `c(s0, i0, r0)` at `times[1]`.
#' @param parms Named list/vector with `lambda`, `gamma`, `alpha`.
#' @param n_sub Number of internal subdivisions per observation interval,
#'   used for the [sird_ga_Qk()] quadrature. Defaults to 10 (even, required
#'   by Simpson's rule).
#' @return A list of class `sird_ga_transition` with components `x_obs`
#'   (matrix of states at `t_0, ..., t_n`), `segments` (one list per
#'   interval, with fine-grid `t`, `x` and `P`, `P[,,1]` always the
#'   identity), and `parms`.
#' @export
sird_ga_transition <- function(times, x0, parms, n_sub = 10L) {
  stopifnot(length(times) >= 2L, all(diff(times) > 0))
  if (n_sub %% 2L == 1L) n_sub <- n_sub + 1L # Simpson needs an even count

  n <- length(times) - 1L

  fine <- times[1L]
  for (k in seq_len(n)) {
    seg <- seq(times[k], times[k + 1L], length.out = n_sub + 1L)
    fine <- if (k == 1L) seg else c(fine, seg[-1L])
  }
  # position of each t_k (k = 0, ..., n) within `fine`
  obs_index <- c(1L, cumsum(rep(n_sub, n)) + 1L)

  rhs <- function(t, y, parms) {
    x <- y[1:3]
    P <- matrix(y[4:12], nrow = 3L)
    dx <- sird_drift(x, parms)
    J <- sird_jacobian(x, parms)
    dP <- J %*% P
    list(c(dx, as.vector(dP)))
  }
  reset_P <- function(t, y, parms) {
    y[4:12] <- as.vector(diag(3))
    y
  }

  y0 <- c(x0, as.vector(diag(3)))
  event_times <- times[-1L][-n] # reset at the start of every interval after the first (t_1, ..., t_{n-1})
  ev <- if (length(event_times) > 0L) list(func = reset_P, time = event_times) else NULL

  out <- deSolve::ode(y = y0, times = fine, func = rhs, parms = parms,
                       method = "lsoda", events = ev)
  out <- as.matrix(out)

  x_obs <- out[obs_index, 2:4, drop = FALSE]
  colnames(x_obs) <- c("s", "i", "r")

  segments <- vector("list", n)
  for (k in seq_len(n)) {
    rows <- obs_index[k]:obs_index[k + 1L]
    x_mat <- out[rows, 2:4, drop = FALSE]
    P_arr <- array(0, dim = c(3L, 3L, length(rows)))
    for (j in seq_along(rows)) P_arr[, , j] <- matrix(out[rows[j], 5:13], nrow = 3L)
    # deSolve records the pre-event value at an event time, so the row at t_{k-1}
    # holds Phi(t_{k-1}, t_{k-2}) rather than the reset value
    P_arr[, , 1L] <- diag(3)
    segments[[k]] <- list(t = out[rows, 1], x = x_mat, P = P_arr)
  }

  structure(list(x_obs = x_obs, segments = segments, parms = parms, times = times),
             class = "sird_ga_transition")
}

#' Resolvent matrix \eqn{A_{k-1} = \Phi(t_k, t_{k-1})} of interval `k`
#' @param tr Object returned by [sird_ga_transition()].
#' @param k Integer interval index, `1, ..., n`.
#' @export
sird_ga_Ak <- function(tr, k) {
  seg <- tr$segments[[k]]
  seg$P[, , dim(seg$P)[3L]]
}

#' Drift increment \eqn{F_k} of interval `k`
#' @inheritParams sird_ga_Ak
#' @export
sird_ga_Fk <- function(tr, k) {
  seg <- tr$segments[[k]]
  x_s <- seg$x[1L, ]; x_t <- seg$x[nrow(seg$x), ]
  A <- sird_ga_Ak(tr, k)
  x_t - as.vector(A %*% x_s)
}

#' State-noise covariance \eqn{Q_k} of interval `k`
#'
#' Computes \eqn{Q_k = N^{-1}\int_{t_{k-1}}^{t_k} \Phi(t_k,s)\,\Sigma(x(s))\,
#' \Phi(t_k,s)^\top ds} by Simpson's rule, with
#' \eqn{\Phi(t_k,s) = \Phi(t_k,t_{k-1})\,\Phi(s,t_{k-1})^{-1}}, using the
#' locally-anchored fundamental matrix stored in `tr$segments[[k]]`.
#' @inheritParams sird_ga_Ak
#' @param N Population size.
#' @export
sird_ga_Qk <- function(tr, k, N) {
  seg <- tr$segments[[k]]
  m <- nrow(seg$x) - 1L
  h <- (seg$t[length(seg$t)] - seg$t[1L]) / m
  A <- sird_ga_Ak(tr, k)

  w <- rep(1, m + 1L)
  if (m >= 2L) {
    w[seq(2L, m, by = 2L)] <- 4
    w[seq(3L, m - 1L, by = 2L)] <- 2
  }
  acc <- matrix(0, 3L, 3L)
  for (j in seq_len(m + 1L)) {
    Phi_ts <- A %*% solve(seg$P[, , j])
    Sigma_s <- sird_diffusion(seg$x[j, ], tr$parms)
    acc <- acc + w[j] * (Phi_ts %*% Sigma_s %*% t(Phi_ts))
  }
  (h / 3) * acc / N
}
