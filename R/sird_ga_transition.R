#' Solve the SIRD-GA state-space transition quantities
#'
#' For each observation interval \eqn{[t_{k-1}, t_k]}, jointly solves the
#' deterministic SIRD trajectory \eqn{x(t)} and the local fundamental
#' (sensitivity) matrix \eqn{\Phi(s, t_{k-1})} of the linearised dynamics,
#' resetting the fundamental matrix to the identity at the start of each
#' interval. Resetting locally (rather than accumulating a single
#' fundamental matrix from \eqn{t_0}) keeps the matrix well-conditioned even
#' over long series, since \eqn{\Phi(t,0)} can otherwise span many orders of
#' magnitude once the epidemic has grown and decayed.
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
#'   interval, with fine-grid `t`, `x` and `P`), and `parms`.
#' @export
sird_ga_transition <- function(times, x0, parms, n_sub = 10L) {
  stopifnot(length(times) >= 2L, all(diff(times) > 0))
  if (n_sub %% 2L == 1L) n_sub <- n_sub + 1L # Simpson needs an even count

  n <- length(times) - 1L
  rhs <- function(t, y, parms) {
    x <- y[1:3]
    P <- matrix(y[4:12], nrow = 3L)
    dx <- sird_drift(x, parms)
    J <- sird_jacobian(x, parms)
    dP <- J %*% P
    list(c(dx, as.vector(dP)))
  }

  x_obs <- matrix(NA_real_, n + 1L, 3L)
  colnames(x_obs) <- c("s", "i", "r")
  x_obs[1L, ] <- x0

  segments <- vector("list", n)
  x_cur <- x0
  for (k in seq_len(n)) {
    tg <- seq(times[k], times[k + 1L], length.out = n_sub + 1L)
    y0 <- c(x_cur, as.vector(diag(3)))
    out <- deSolve::ode(y = y0, times = tg, func = rhs, parms = parms, method = "lsoda")
    out <- as.matrix(out)
    x_mat <- out[, 2:4, drop = FALSE]
    P_arr <- array(0, dim = c(3L, 3L, nrow(out)))
    for (j in seq_len(nrow(out))) P_arr[, , j] <- matrix(out[j, 5:13], nrow = 3L)

    segments[[k]] <- list(t = tg, x = x_mat, P = P_arr)
    x_cur <- x_mat[nrow(x_mat), ]
    x_obs[k + 1L, ] <- x_cur
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
