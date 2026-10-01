#' Exact stochastic simulation of a general jump model (Gillespie algorithm)
#'
#' Dimension-generic analogue of [sird_gillespie()], driven by a
#' [new_jump_model()] of any dimension.
#'
#' @param model A `jump_model`.
#' @param N Population size.
#' @param x0_counts Named integer vector of initial head counts, one per
#'   `model$compartments`.
#' @param parms Named list of rate parameters.
#' @param Tmax,Delta As in [sird_gillespie()].
#' @param max_events Safety cap on the number of simulated events.
#' @return A data frame with column `time` and one column per compartment
#'   in `model$compartments` (head counts), recorded on
#'   `seq(0, Tmax, by = Delta)`.
#' @export
general_gillespie <- function(model, N, x0_counts, parms, Tmax, Delta = 1,
                               max_events = 50L * N) {
  d <- length(model$compartments)
  state <- x0_counts[model$compartments]
  t <- 0
  grid <- seq(0, Tmax, by = Delta)
  rec <- matrix(NA_real_, length(grid), d)
  colnames(rec) <- model$compartments
  g_idx <- 1L
  while (g_idx <= length(grid) && grid[g_idx] <= t) { rec[g_idx, ] <- state; g_idx <- g_idx + 1L }

  n_events <- 0L
  while (t < Tmax && n_events < max_events) {
    x_frac <- state / N
    rates <- vapply(model$events, function(ev) ev$rate(x_frac, parms), numeric(1)) * N
    r_tot <- sum(rates)
    if (r_tot <= 0 || !is.finite(r_tot)) break
    t <- t + stats::rexp(1, rate = r_tot)
    n_events <- n_events + 1L
    while (g_idx <= length(grid) && grid[g_idx] < t) { rec[g_idx, ] <- state; g_idx <- g_idx + 1L }
    if (g_idx > length(grid)) break
    ev_idx <- sample.int(length(model$events), 1L, prob = rates / r_tot)
    state <- pmax(state + model$events[[ev_idx]]$jump, 0)
  }
  while (g_idx <= length(grid)) { rec[g_idx, ] <- state; g_idx <- g_idx + 1L }
  data.frame(time = grid, rec)
}

#' Simulate partially-observed surveillance data from a general jump model
#' @inheritParams general_gillespie
#' @param obs_compartment Compartment observed (with under-reporting/noise).
#' @param p,tau2 As in [sird_ga_simulate_obs()].
#' @export
general_ga_simulate_obs <- function(model, N, x0_counts, parms, obs_compartment,
                                     p, tau2, Tmax, Delta = 1) {
  latent <- general_gillespie(model, N, x0_counts, parms, Tmax, Delta)
  obs_counts <- latent[[obs_compartment]]
  O1 <- stats::rbinom(length(obs_counts), size = pmax(round(obs_counts), 0), prob = p)
  O2 <- stats::rnorm(length(obs_counts), mean = 0, sd = sqrt(pmax(tau2 * obs_counts, 0)))
  list(latent = latent, y = pmax(O1 + O2, 0))
}

#' Maximum-likelihood estimation for a general jump model
#'
#' Dimension-generic analogue of [sird_ga_fit()]. All explicit
#' compartments except one (`free_init`) are held fixed at `start`'s
#' values for the initial condition; by default only the observed
#' compartment's initial fraction is estimated, with all others fixed
#' (mirroring [sird_ga_fit()]'s treatment of `r0`).
#'
#' @param y,N,Delta As in [sird_ga_fit()].
#' @param model A `jump_model`.
#' @param obs_compartment Name of the observed compartment.
#' @param start Named list with one entry per rate parameter, one entry
#'   per compartment in `model$compartments` (the initial condition), `p`
#'   and `tau2`.
#' @param free_init Character vector of compartment names whose initial
#'   fraction is estimated (default: just `obs_compartment`); all other
#'   compartments are held fixed at their `start` value, and the
#'   "susceptible-like" first compartment is always treated as `1 - sum(
#'   other explicit compartments)` so the initial condition stays on the
#'   simplex.
#' @param p_max,n_sub,Q0_diag,control,hessian As in [sird_ga_fit()].
#' @export
general_ga_fit <- function(y, N, model, obs_compartment, Delta = 1, start,
                            free_init = obs_compartment, p_max = 0.999,
                            n_sub = 10L, Q0_diag = 1e-8,
                            control = list(maxit = 3000, reltol = 1e-10),
                            hessian = FALSE) {
  start <- as.list(start)
  comp <- model$compartments
  rate_names <- setdiff(names(start), c(comp, "p", "tau2"))
  fixed_init <- setdiff(comp[-1], free_init) # first compartment is always derived

  to_free <- function(th) {
    v <- c(vapply(rate_names, function(nm) log(th[[nm]]), numeric(1)),
           vapply(free_init, function(nm) stats::qlogis(th[[nm]]), numeric(1)),
           stats::qlogis(th$p / p_max), log(th$tau2))
    names(v) <- c(rate_names, free_init, "p", "tau2")
    v
  }
  from_free <- function(v) {
    rates <- exp(v[rate_names]); names(rates) <- rate_names
    init_free <- stats::plogis(v[free_init]); names(init_free) <- free_init
    p <- p_max * stats::plogis(v[["p"]]); tau2 <- exp(v[["tau2"]])
    init <- numeric(length(comp)); names(init) <- comp
    for (nm in fixed_init) init[nm] <- start[[nm]]
    for (nm in free_init) init[nm] <- init_free[nm]
    init[comp[1]] <- 1 - sum(init[comp[-1]])
    c(rates, init, p = p, tau2 = tau2)
  }

  negloglik <- function(v) {
    theta <- from_free(v)
    ll <- general_ga_loglik(theta, y, N, model, obs_compartment, Delta, n_sub, Q0_diag)
    if (!is.finite(ll)) return(1e10)
    -ll
  }

  v0 <- to_free(start)
  fit <- stats::optim(v0, negloglik, method = "Nelder-Mead", control = control, hessian = hessian)
  list(theta = from_free(fit$par), loglik = -fit$value, convergence = fit$convergence,
       optim_fit = fit, settings = list(N = N, Delta = Delta, n_sub = n_sub, Q0_diag = Q0_diag,
                                          model = model, obs_compartment = obs_compartment))
}
