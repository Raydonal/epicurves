#' Exact stochastic simulation of the SIRD jump process (Gillespie algorithm)
#'
#' Simulates one exact trajectory of the continuous-time, density-dependent
#' Markov jump process underlying the SIRD-GA model, using the direct
#' Gillespie method, and returns the compartment counts recorded at an
#' equally-spaced observation grid.
#'
#' @param N Population size (integer head count).
#' @param S0,I0,R0 Initial compartment counts (integers, `S0+I0+R0 <= N`).
#' @param lambda,gamma,alpha Transmission, recovery and disease-induced
#'   mortality rates.
#' @param Tmax Final observation time.
#' @param Delta Spacing of the recorded observation grid (default `1`).
#' @param max_events Safety cap on the number of simulated events.
#' @return A data frame with columns `time`, `S`, `I`, `R`, `D`, recorded on
#'   the grid `seq(0, Tmax, by = Delta)` (state held constant between
#'   events, as is standard for jump processes).
#' @export
sird_gillespie <- function(N, S0, I0, R0, lambda, gamma, alpha, Tmax,
                            Delta = 1, max_events = 50L * N) {
  t <- 0; S <- S0; I <- I0; R <- R0
  grid <- seq(0, Tmax, by = Delta)
  rec <- matrix(NA_real_, length(grid), 4L)
  colnames(rec) <- c("S", "I", "R", "D")
  g_idx <- 1L
  while (g_idx <= length(grid) && grid[g_idx] <= t) {
    rec[g_idx, ] <- c(S, I, R, N - S - I - R)
    g_idx <- g_idx + 1L
  }

  n_events <- 0L
  while (t < Tmax && I > 0 && n_events < max_events) {
    r_inf <- lambda * S * I / N
    r_rec <- gamma * I
    r_death <- alpha * I
    r_tot <- r_inf + r_rec + r_death
    if (r_tot <= 0) break
    t <- t + stats::rexp(1, rate = r_tot)
    n_events <- n_events + 1L

    while (g_idx <= length(grid) && grid[g_idx] < t) {
      rec[g_idx, ] <- c(S, I, R, N - S - I - R)
      g_idx <- g_idx + 1L
    }
    if (g_idx > length(grid)) break

    u <- stats::runif(1, 0, r_tot)
    if (u < r_inf) {
      S <- S - 1L; I <- I + 1L
    } else if (u < r_inf + r_rec) {
      I <- I - 1L; R <- R + 1L
    } else {
      I <- I - 1L
    }
  }
  while (g_idx <= length(grid)) {
    rec[g_idx, ] <- c(S, I, R, N - S - I - R)
    g_idx <- g_idx + 1L
  }
  data.frame(time = grid, rec)
}

#' Simulate partially-observed SIRD-GA surveillance data
#'
#' Wraps [sird_gillespie()] with the two-part observation-error model
#' (binomial under-reporting plus Gaussian measurement error) used
#' throughout the SIRD-GA manuscript.
#'
#' @inheritParams sird_gillespie
#' @param p Reporting probability in `(0,1)`.
#' @param tau2 Measurement-error variance parameter.
#' @return A list with `latent` (output of [sird_gillespie()]) and `y`
#'   (numeric vector of simulated observations, one per grid time).
#' @export
sird_ga_simulate_obs <- function(N, S0, I0, R0, lambda, gamma, alpha,
                                  p, tau2, Tmax, Delta = 1) {
  latent <- sird_gillespie(N, S0, I0, R0, lambda, gamma, alpha, Tmax, Delta)
  I_t <- latent$I
  O1 <- stats::rbinom(length(I_t), size = pmax(round(I_t), 0), prob = p)
  O2 <- stats::rnorm(length(I_t), mean = 0, sd = sqrt(pmax(tau2 * I_t, 0)))
  y <- pmax(O1 + O2, 0)
  list(latent = latent, y = y)
}
