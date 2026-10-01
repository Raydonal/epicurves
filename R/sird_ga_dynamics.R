#' Deterministic SIRD drift vector
#'
#' Right-hand side of the normalised SIRD ODE system
#' \eqn{\dot s = -\lambda s i}, \eqn{\dot i = \lambda s i - (\gamma+\alpha) i},
#' \eqn{\dot r = \gamma i}.
#'
#' @param state Numeric vector `c(s, i, r)`.
#' @param parms Named numeric vector/list with elements `lambda`, `gamma`, `alpha`.
#' @return Numeric vector of length 3, `b(x)`.
#' @export
sird_drift <- function(state, parms) {
  s <- state[1L]; i <- state[2L]; r <- state[3L]
  lambda <- parms[["lambda"]]; gamma <- parms[["gamma"]]; alpha <- parms[["alpha"]]
  c(
    -lambda * s * i,
    lambda * s * i - (gamma + alpha) * i,
    gamma * i
  )
}

#' Jacobian of the SIRD drift vector
#'
#' Analytic Jacobian \eqn{\nabla_x b(x)} of [sird_drift()], used to build the
#' state-transition (resolvent) matrix \eqn{\Phi(t,s)} of the SIRD-GA model.
#'
#' @inheritParams sird_drift
#' @return A 3x3 numeric matrix.
#' @export
sird_jacobian <- function(state, parms) {
  s <- state[1L]; i <- state[2L]
  lambda <- parms[["lambda"]]; gamma <- parms[["gamma"]]; alpha <- parms[["alpha"]]
  matrix(
    c(-lambda * i, -lambda * s, 0,
      lambda * i,  lambda * s - gamma - alpha, 0,
      0,           gamma, 0),
    nrow = 3L, byrow = TRUE
  )
}

#' Diffusion covariance matrix of the SIRD jump process
#'
#' Symmetric matrix \eqn{\Sigma(x)} from the second-order moments of the jump
#' rates, obtained from the three elementary events (infection, recovery,
#' disease-induced death).
#'
#' @inheritParams sird_drift
#' @return A 3x3 symmetric numeric matrix.
#' @export
sird_diffusion <- function(state, parms) {
  s <- state[1L]; i <- state[2L]
  lambda <- parms[["lambda"]]; gamma <- parms[["gamma"]]; alpha <- parms[["alpha"]]
  li <- lambda * s * i
  matrix(
    c(li,        -li,                 0,
      -li,  li + (gamma + alpha) * i, -gamma * i,
      0,         -gamma * i,          gamma * i),
    nrow = 3L, byrow = TRUE
  )
}

#' Diffusion "square-root" factor of the SIRD jump process
#'
#' Lower-triangular-like matrix \eqn{\sigma(x)} such that
#' \eqn{\sigma(x)\sigma(x)^\top = \Sigma(x)}, driven by a 3-dimensional
#' Brownian motion with one independent noise channel per elementary event.
#'
#' @inheritParams sird_drift
#' @return A 3x3 numeric matrix.
#' @export
sird_sigma <- function(state, parms) {
  s <- state[1L]; i <- state[2L]
  lambda <- parms[["lambda"]]; gamma <- parms[["gamma"]]; alpha <- parms[["alpha"]]
  ga <- gamma + alpha
  matrix(
    c(sqrt(lambda * s * i),              0,                      0,
      -sqrt(lambda * s * i),             sqrt(ga * i),           0,
      0,                        -gamma * sqrt(i / ga),  sqrt(gamma * alpha * i / ga)),
    nrow = 3L, byrow = TRUE
  )
}

#' Basic reproduction number of the SIRD model
#'
#' @param lambda,gamma,alpha Transmission, recovery and disease-induced
#'   mortality rates.
#' @return Numeric scalar \eqn{R_0 = \lambda / (\gamma+\alpha)}.
#' @export
sird_R0 <- function(lambda, gamma, alpha) {
  lambda / (gamma + alpha)
}

#' Solve the deterministic normalised SIRD ODE system
#'
#' @param times Numeric vector of output times (must include `t0`).
#' @param x0 Initial state `c(s0, i0, r0)`.
#' @param parms Named list/vector with `lambda`, `gamma`, `alpha`.
#' @return A data frame with columns `time`, `s`, `i`, `r`, `d` (data frame as
#'   returned by [deSolve::ode()], with the implicit death fraction added).
#' @export
sird_ode_solve <- function(times, x0, parms) {
  rhs <- function(t, state, parms) list(sird_drift(state, parms))
  out <- deSolve::ode(y = stats::setNames(x0, c("s", "i", "r")),
                       times = times, func = rhs, parms = parms,
                       method = "lsoda")
  out <- as.data.frame(out)
  out$d <- 1 - out$s - out$i - out$r
  out
}
