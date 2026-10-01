#' Solve the Gaussian-approximation transition quantities for a general jump model
#'
#' Dimension-generic analogue of [sird_ga_transition()]: identical
#' event-reset construction (a single [deSolve::ode()] call per model,
#' with the fundamental matrix reset to the identity at the start of each
#' observation interval), but driven by a [new_jump_model()] object of any
#' dimension `d = length(model$compartments)` instead of the hard-coded
#' `d = 3` SIRD case. [sird_ga_transition()] remains the validated,
#' hand-derived implementation used for the main SIRD-GA results; this
#' generic engine is what [general_ga_fit()] uses for extensions (SEIRD,
#' SIRDV, ...).
#'
#' @param times,x0,n_sub As in [sird_ga_transition()].
#' @param model A `jump_model` (see [new_jump_model()]).
#' @param parms Named list/vector of the model's rate parameters.
#' @export
general_ga_transition <- function(times, x0, model, parms, n_sub = 10L) {
  stopifnot(length(times) >= 2L, all(diff(times) > 0))
  if (n_sub %% 2L == 1L) n_sub <- n_sub + 1L
  d <- length(x0)
  n <- length(times) - 1L

  rhs <- function(t, y, parms) {
    x <- y[1:d]
    P <- matrix(y[(d + 1):(d + d * d)], nrow = d)
    dx <- jm_drift(model, x, parms)
    J <- jm_jacobian(model, x, parms)
    dP <- J %*% P
    list(c(dx, as.vector(dP)))
  }

  x_obs <- matrix(NA_real_, n + 1L, d)
  colnames(x_obs) <- model$compartments
  x_obs[1L, ] <- x0

  segments <- vector("list", n)
  x_cur <- x0
  for (k in seq_len(n)) {
    tg <- seq(times[k], times[k + 1L], length.out = n_sub + 1L)
    y0 <- c(x_cur, as.vector(diag(d)))
    out <- deSolve::ode(y = y0, times = tg, func = rhs, parms = parms, method = "lsoda")
    out <- as.matrix(out)
    x_mat <- out[, 2:(d + 1L), drop = FALSE]
    P_arr <- array(0, dim = c(d, d, nrow(out)))
    for (j in seq_len(nrow(out))) P_arr[, , j] <- matrix(out[j, (d + 2):(d + 1L + d * d)], nrow = d)

    segments[[k]] <- list(t = tg, x = x_mat, P = P_arr)
    x_cur <- x_mat[nrow(x_mat), ]
    x_obs[k + 1L, ] <- x_cur
  }

  structure(list(x_obs = x_obs, segments = segments, model = model, parms = parms,
                  times = times, d = d),
             class = "general_ga_transition")
}

#' @describeIn general_ga_transition Resolvent matrix of interval `k`.
#' @param tr Object returned by [general_ga_transition()].
#' @param k Integer interval index.
#' @export
general_ga_Ak <- function(tr, k) {
  seg <- tr$segments[[k]]
  seg$P[, , dim(seg$P)[3L]]
}

#' @describeIn general_ga_transition Drift increment of interval `k`.
#' @export
general_ga_Fk <- function(tr, k) {
  seg <- tr$segments[[k]]
  x_s <- seg$x[1L, ]; x_t <- seg$x[nrow(seg$x), ]
  A <- general_ga_Ak(tr, k)
  x_t - as.vector(A %*% x_s)
}

#' @describeIn general_ga_transition State-noise covariance of interval `k` (Simpson's rule).
#' @param N Population size.
#' @export
general_ga_Qk <- function(tr, k, N) {
  seg <- tr$segments[[k]]
  m <- nrow(seg$x) - 1L
  h <- (seg$t[length(seg$t)] - seg$t[1L]) / m
  A <- general_ga_Ak(tr, k)
  d <- tr$d

  w <- rep(1, m + 1L)
  if (m >= 2L) {
    w[seq(2L, m, by = 2L)] <- 4
    w[seq(3L, m - 1L, by = 2L)] <- 2
  }
  acc <- matrix(0, d, d)
  for (j in seq_len(m + 1L)) {
    Phi_ts <- A %*% solve(seg$P[, , j])
    Sigma_s <- jm_diffusion(tr$model, seg$x[j, ], tr$parms)
    acc <- acc + w[j] * (Phi_ts %*% Sigma_s %*% t(Phi_ts))
  }
  (h / 3) * acc / N
}

#' Kalman-filter log-likelihood for a general jump model
#'
#' Dimension-generic analogue of [sird_ga_loglik()]. A single compartment
#' is observed (binomial under-reporting plus Gaussian measurement error,
#' as in the SIRD-GA observation model), identified by name in
#' `obs_compartment`.
#'
#' @param theta Named numeric vector: the model's rate parameters, its
#'   initial condition (one named entry per `model$compartments`), `p`
#'   and `tau2`.
#' @param y,N,Delta,n_sub,Q0_diag,keep_states As in [sird_ga_loglik()].
#' @param model A `jump_model`.
#' @param obs_compartment Name of the observed compartment (must be one
#'   of `model$compartments`).
#' @export
general_ga_loglik <- function(theta, y, N, model, obs_compartment,
                               Delta = 1, n_sub = 10L, Q0_diag = 1e-8,
                               keep_states = FALSE) {
  d <- length(model$compartments)
  x0 <- theta[model$compartments]
  p <- theta[["p"]]; tau2 <- theta[["tau2"]]
  rate_names <- setdiff(names(theta), c(model$compartments, "p", "tau2"))
  parms <- as.list(theta[rate_names])

  bad <- any(unlist(parms) <= 0) || p <= 0 || p >= 1 || tau2 <= 0 ||
    x0[[obs_compartment]] <= 0 || any(x0 < 0) || sum(x0) > 1 + 1e-6
  if (bad) return(if (keep_states) list(loglik = -Inf) else -Inf)

  n <- length(y) - 1L
  times <- (0:n) * Delta
  y_frac <- y / N

  tr <- tryCatch(general_ga_transition(times, x0, model, parms, n_sub = n_sub),
                  error = function(e) NULL)
  obs_col <- which(model$compartments == obs_compartment)
  if (is.null(tr) || anyNA(tr$x_obs) || any(tr$x_obs[, obs_col] <= 0)) {
    return(if (keep_states) list(loglik = -Inf) else -Inf)
  }

  B <- matrix(0, 1L, d); B[1L, obs_col] <- p
  Xbar <- x0
  Sbar <- diag(Q0_diag, d)
  loglik <- 0
  filtered <- matrix(NA_real_, n + 1L, d)
  predicted <- matrix(NA_real_, n + 1L, d)

  step_update <- function(Xhat, Shat, obs_val, y_obs) {
    Rk <- (1 / N) * (p * (1 - p) + tau2) * obs_val
    nu <- y_obs - as.vector(B %*% Xhat)
    Gamma <- as.vector(B %*% Shat %*% t(B)) + Rk
    if (!is.finite(Gamma) || Gamma <= 0) return(NULL)
    ll <- -0.5 * (log(2 * pi) + log(Gamma) + nu^2 / Gamma)
    K <- (Shat %*% t(B)) / Gamma
    list(ll = ll, Xbar = Xhat + as.vector(K * nu), Sbar = Shat - K %*% B %*% Shat)
  }

  up <- step_update(Xbar, Sbar, tr$x_obs[1L, obs_col], y_frac[1L])
  if (is.null(up)) return(if (keep_states) list(loglik = -Inf) else -Inf)
  loglik <- loglik + up$ll; Xbar <- up$Xbar; Sbar <- up$Sbar
  filtered[1L, ] <- Xbar; predicted[1L, ] <- x0

  for (k in seq_len(n)) {
    A <- general_ga_Ak(tr, k); Fk <- general_ga_Fk(tr, k); Qk <- general_ga_Qk(tr, k, N)
    Xhat <- Fk + as.vector(A %*% Xbar)
    Shat <- A %*% Sbar %*% t(A) + Qk
    up <- step_update(Xhat, Shat, tr$x_obs[k + 1L, obs_col], y_frac[k + 1L])
    if (is.null(up)) return(if (keep_states) list(loglik = -Inf) else -Inf)
    loglik <- loglik + up$ll; Xbar <- up$Xbar; Sbar <- up$Sbar
    filtered[k + 1L, ] <- Xbar; predicted[k + 1L, ] <- Xhat
  }

  if (!keep_states) return(as.numeric(loglik))
  colnames(filtered) <- colnames(predicted) <- model$compartments
  list(loglik = as.numeric(loglik), filtered = filtered, predicted = predicted,
       fitted = N * p * predicted[, obs_col], times = times, N = N)
}
