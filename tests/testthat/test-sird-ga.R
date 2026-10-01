test_that("sird_ga_loglik is finite at a sensible parameter value and -Inf outside the parameter space", {
  set.seed(10)
  N <- 20000
  theta <- c(lambda = 0.4, gamma = 0.15, alpha = 0.08, s0 = 0.999, i0 = 0.001,
             r0 = 0, p = 0.6, tau2 = 0.2)
  sim <- sird_ga_simulate_obs(N, S0 = round(N * 0.999), I0 = round(N * 0.001), R0 = 0,
                               lambda = theta[["lambda"]], gamma = theta[["gamma"]],
                               alpha = theta[["alpha"]], p = theta[["p"]], tau2 = theta[["tau2"]],
                               Tmax = 60, Delta = 1)
  ll <- sird_ga_loglik(theta, sim$y, N, Delta = 1, n_sub = 6)
  expect_true(is.finite(ll))

  bad_theta <- theta; bad_theta[["p"]] <- 1.5
  expect_equal(sird_ga_loglik(bad_theta, sim$y, N), -Inf)
})

test_that("sird_ga_fit_multistart returns one result per starting infectious period, sorted by loglik", {
  skip_on_cran()
  set.seed(12)
  N <- 20000
  theta_true <- c(lambda = 0.4, gamma = 0.15, alpha = 0.08, s0 = 0.999, i0 = 0.001,
                   r0 = 0, p = 0.6, tau2 = 0.2)
  sim <- sird_ga_simulate_obs(N, S0 = round(N * 0.999), I0 = round(N * 0.001), R0 = 0,
                               lambda = theta_true[["lambda"]], gamma = theta_true[["gamma"]],
                               alpha = theta_true[["alpha"]], p = theta_true[["p"]],
                               tau2 = theta_true[["tau2"]], Tmax = 80, Delta = 1)
  ms <- sird_ga_fit_multistart(sim$y, N, Delta = 1, r0_guess = 1.5,
                                infectious_period_days = c(3, 5, 8), n_sub = 6,
                                control = list(maxit = 500, reltol = 1e-8))
  expect_equal(nrow(ms$summary), 3L)
  expect_true(all(diff(ms$summary$loglik) <= 0)) # sorted best-first (non-increasing)
  expect_identical(ms$best$loglik, max(ms$summary$loglik, na.rm = TRUE))
})

test_that("sird_ga_fit recovers R0 and lambda to reasonable accuracy on simulated data", {
  skip_on_cran()
  set.seed(11)
  N <- 20000
  theta_true <- c(lambda = 0.4, gamma = 0.15, alpha = 0.08, s0 = 0.999, i0 = 0.001,
                   r0 = 0, p = 0.6, tau2 = 0.2)
  sim <- sird_ga_simulate_obs(N, S0 = round(N * 0.999), I0 = round(N * 0.001), R0 = 0,
                               lambda = theta_true[["lambda"]], gamma = theta_true[["gamma"]],
                               alpha = theta_true[["alpha"]], p = theta_true[["p"]],
                               tau2 = theta_true[["tau2"]], Tmax = 100, Delta = 1)

  fit <- sird_ga_fit(sim$y, N, Delta = 1,
                      start = list(lambda = 0.3, gamma = 0.1, alpha = 0.1,
                                   i0 = 0.002, p = 0.5, tau2 = 0.3),
                      n_sub = 6, hessian = FALSE,
                      control = list(maxit = 800, reltol = 1e-8))

  R0_true <- sird_R0(theta_true[["lambda"]], theta_true[["gamma"]], theta_true[["alpha"]])
  R0_hat <- sird_R0(fit$theta[["lambda"]], fit$theta[["gamma"]], fit$theta[["alpha"]])
  expect_equal(R0_hat, R0_true, tolerance = 0.3)
  expect_equal(unname(fit$theta[["lambda"]]), unname(theta_true[["lambda"]]), tolerance = 0.3)
})
