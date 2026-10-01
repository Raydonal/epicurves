test_that("jm_drift/jm_diffusion for sird_model() match the hand-derived SIRD functions", {
  parms <- list(lambda = 0.4, gamma = 0.15, alpha = 0.08)
  state <- c(s = 0.7, i = 0.2, r = 0.1)
  model <- sird_model()
  # unname(): sird_drift()'s internal arithmetic on named scalars incidentally
  # propagates element names that differ from jm_drift()'s; only the numeric
  # values are meant to be compared here.
  expect_equal(unname(jm_drift(model, state, parms)), unname(sird_drift(state, parms)), tolerance = 1e-10)
  expect_equal(unname(jm_diffusion(model, state, parms)), unname(sird_diffusion(state, parms)), tolerance = 1e-10)
  expect_equal(unname(jm_jacobian(model, state, parms)), unname(sird_jacobian(state, parms)), tolerance = 1e-5)
})

test_that("general_ga_loglik(sird_model()) matches sird_ga_loglik() on identical data", {
  skip_on_cran()
  set.seed(42)
  N <- 20000
  theta_full <- c(lambda = 0.4, gamma = 0.15, alpha = 0.08, s0 = 0.999, i0 = 0.001,
                   r0 = 0, p = 0.6, tau2 = 0.2)
  sim <- sird_ga_simulate_obs(N, round(N * 0.999), round(N * 0.001), 0,
                               theta_full[["lambda"]], theta_full[["gamma"]], theta_full[["alpha"]],
                               theta_full[["p"]], theta_full[["tau2"]], Tmax = 60, Delta = 1)
  ll_sird <- sird_ga_loglik(theta_full, sim$y, N, Delta = 1, n_sub = 6)

  theta_general <- c(lambda = 0.4, gamma = 0.15, alpha = 0.08, s = 0.999, i = 0.001,
                      r = 0, p = 0.6, tau2 = 0.2)
  ll_general <- general_ga_loglik(theta_general, sim$y, N, sird_model(), "i", Delta = 1, n_sub = 6)

  # numerical (central-difference) vs analytic Jacobian -> small, bounded discrepancy
  expect_equal(ll_general, ll_sird, tolerance = 1e-3)
})

test_that("SEIRD-GA preserves mass balance (s+e+i+r+d = 1, d implicit, non-increasing with deaths)", {
  N <- 20000
  model <- seird_model()
  parms <- list(lambda = 0.5, sigma = 0.3, gamma = 0.15, alpha = 0.05)
  tr <- general_ga_transition(0:60, c(s = 0.999, e = 0, i = 0.001, r = 0), model, parms, n_sub = 6)
  total <- rowSums(tr$x_obs)
  expect_true(all(total <= 1 + 1e-8))
  expect_true(all(diff(total) <= 1e-8)) # implicit death fraction is non-decreasing
})

test_that("SEIRD-GA recovers simulated parameters to reasonable accuracy", {
  skip_on_cran()
  set.seed(7)
  N <- 20000
  model <- seird_model()
  parms <- list(lambda = 0.5, sigma = 0.3, gamma = 0.15, alpha = 0.05)
  sim <- general_ga_simulate_obs(model, N, c(s = round(N * 0.999), e = 0, i = round(N * 0.001), r = 0),
                                   parms, obs_compartment = "i", p = 0.6, tau2 = 0.2, Tmax = 80, Delta = 1)

  fit <- general_ga_fit(sim$y, N, model, "i", Delta = 1,
                         start = list(lambda = 0.4, sigma = 0.2, gamma = 0.1, alpha = 0.1,
                                       s = 0.999, e = 0, i = 0.002, r = 0, p = 0.5, tau2 = 0.3),
                         free_init = "i", n_sub = 6, control = list(maxit = 1000, reltol = 1e-9))

  expect_equal(unname(fit$theta[["lambda"]]), 0.5, tolerance = 0.3)
  R0_hat <- fit$theta[["lambda"]] / (fit$theta[["gamma"]] + fit$theta[["alpha"]])
  R0_true <- 0.5 / (0.15 + 0.05)
  expect_equal(unname(R0_hat), R0_true, tolerance = 0.3)
})
