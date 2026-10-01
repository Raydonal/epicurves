#' Kalman-filter log-likelihood of the SIRD-GA state-space model
#'
#' Runs the predict/update Kalman recursion of the linear-Gaussian
#' state-space representation of the SIRD-GA model and returns the
#' approximate log-likelihood of the observed (partially reported, noisy)
#' infectious counts.
#'
#' @param theta Named numeric vector with elements `lambda`, `gamma`,
#'   `alpha`, `s0`, `i0`, `r0`, `p`, `tau2` (all strictly positive except
#'   `r0`, which may be zero; `p` in `(0,1)`).
#' @param y Numeric vector of observed infectious **head counts**
#'   \eqn{y_0, y_1, \ldots, y_n} (not normalised; normalisation by `N` is
#'   performed internally, matching the \eqn{O_N(t_k)=O(t_k)/N} observation
#'   model of the manuscript).
#' @param N Population size.
#' @param Delta Time step between consecutive observations (default `1`).
#' @param n_sub Internal quadrature resolution, passed to
#'   [sird_ga_transition()].
#' @param Q0_diag Variance assigned to each component of the initial-state
#'   distribution (small, reflects near-certainty about `x0`).
#' @param keep_states Logical; if `TRUE`, also return the filtered state
#'   trajectory (in normalised units) and the predictive mean on the
#'   original (head-count) scale.
#'
#' @details
#' The state-space system is
#' \deqn{X_k = F_k(\eta,\Delta) + A_{k-1}(\eta,\Delta) X_{k-1} + U_k, \quad
#'       U_k \sim N_3(0, Q_k(\eta,\Delta)),}
#' \deqn{Y_k/N = p \, X_{k,2} + V_k, \quad V_k \sim N(0, R_k(\theta)),}
#' with \eqn{X_k=(S_k,I_k,R_k)/N}, observation operator `B = (0, p, 0)` and
#' \eqn{R_k(\theta) = N^{-1}(p(1-p)+\tau^2) \, i(\eta,t_k)}, following the
#' observation model of Narci et al. (2021) extended to the three-event
#' SIRD jump process.
#'
#' @return If `keep_states` is `FALSE` (default), the scalar approximate
#'   log-likelihood. Otherwise a list with elements `loglik`, `filtered`
#'   (matrix of filtered means \eqn{\bar X_k}, normalised units), `predicted`
#'   (matrix of one-step-ahead predicted means \eqn{\hat X_k}, normalised
#'   units) and `fitted` (predictive mean of `y`, on the original head-count
#'   scale, i.e. `N * p * predicted[, "i"]`).
#' @export
sird_ga_loglik <- function(theta, y, N, Delta = 1, n_sub = 10L,
                            Q0_diag = 1e-8, keep_states = FALSE) {
  lambda <- theta[["lambda"]]; gamma <- theta[["gamma"]]; alpha <- theta[["alpha"]]
  s0 <- theta[["s0"]]; i0 <- theta[["i0"]]; r0 <- theta[["r0"]]
  p <- theta[["p"]]; tau2 <- theta[["tau2"]]

  bad <- any(c(lambda, gamma, alpha, tau2) <= 0) || p <= 0 || p >= 1 ||
    i0 <= 0 || s0 <= 0 || (s0 + i0 + r0) > 1 + 1e-6
  if (bad) return(if (keep_states) list(loglik = -Inf) else -Inf)

  n <- length(y) - 1L
  times <- (0:n) * Delta
  parms <- list(lambda = lambda, gamma = gamma, alpha = alpha)
  x0 <- c(s0, i0, r0)
  y_frac <- y / N

  tr <- tryCatch(sird_ga_transition(times, x0, parms, n_sub = n_sub),
                  error = function(e) NULL)
  if (is.null(tr) || anyNA(tr$x_obs) || any(tr$x_obs[, "i"] <= 0)) {
    return(if (keep_states) list(loglik = -Inf) else -Inf)
  }

  B <- matrix(c(0, p, 0), nrow = 1L)
  Xbar <- x0
  Sbar <- diag(Q0_diag, 3L)
  loglik <- 0
  filtered <- matrix(NA_real_, n + 1L, 3L)
  predicted <- matrix(NA_real_, n + 1L, 3L)

  step_update <- function(Xhat, Shat, i_t, y_obs) {
    Rk <- (1 / N) * (p * (1 - p) + tau2) * i_t
    nu <- y_obs - as.vector(B %*% Xhat)
    Gamma <- as.vector(B %*% Shat %*% t(B)) + Rk
    if (!is.finite(Gamma) || Gamma <= 0) return(NULL)
    ll <- -0.5 * (log(2 * pi) + log(Gamma) + nu^2 / Gamma)
    K <- (Shat %*% t(B)) / Gamma
    Xbar_new <- Xhat + as.vector(K * nu)
    Sbar_new <- Shat - K %*% B %*% Shat
    list(ll = ll, Xbar = Xbar_new, Sbar = Sbar_new)
  }

  up <- step_update(Xbar, Sbar, tr$x_obs[1L, "i"], y_frac[1L])
  if (is.null(up)) return(if (keep_states) list(loglik = -Inf) else -Inf)
  loglik <- loglik + up$ll
  Xbar <- up$Xbar; Sbar <- up$Sbar
  filtered[1L, ] <- Xbar; predicted[1L, ] <- x0

  for (k in seq_len(n)) {
    A <- sird_ga_Ak(tr, k)
    Fk <- sird_ga_Fk(tr, k)
    Qk <- sird_ga_Qk(tr, k, N)

    Xhat <- Fk + as.vector(A %*% Xbar)
    Shat <- A %*% Sbar %*% t(A) + Qk
    i_t <- tr$x_obs[k + 1L, "i"]

    up <- step_update(Xhat, Shat, i_t, y_frac[k + 1L])
    if (is.null(up)) return(if (keep_states) list(loglik = -Inf) else -Inf)
    loglik <- loglik + up$ll
    Xbar <- up$Xbar; Sbar <- up$Sbar
    filtered[k + 1L, ] <- Xbar; predicted[k + 1L, ] <- Xhat
  }

  if (!keep_states) return(as.numeric(loglik))
  colnames(filtered) <- colnames(predicted) <- c("s", "i", "r")
  list(loglik = as.numeric(loglik), filtered = filtered, predicted = predicted,
       fitted = N * p * predicted[, "i"], times = times, N = N)
}
