test_that("R0 formula matches lambda / (gamma + alpha)", {
  expect_equal(sird_R0(0.4, 0.15, 0.08), 0.4 / 0.23)
})

test_that("the SIRD ODE conserves mass (s + i + r + d = 1)", {
  parms <- list(lambda = 0.4, gamma = 0.15, alpha = 0.08)
  out <- sird_ode_solve(times = seq(0, 100, by = 1), x0 = c(s = 0.999, i = 0.001, r = 0),
                         parms = parms)
  expect_true(all(abs(out$s + out$i + out$r + out$d - 1) < 1e-8))
  expect_true(all(out$s >= -1e-8 & out$i >= -1e-8 & out$r >= -1e-8))
})

test_that("sigma is a valid square-root factor of Sigma", {
  parms <- list(lambda = 0.4, gamma = 0.15, alpha = 0.08)
  state <- c(s = 0.7, i = 0.2, r = 0.1)
  Sigma <- sird_diffusion(state, parms)
  sigma <- sird_sigma(state, parms)
  expect_equal(sigma %*% t(sigma), Sigma, tolerance = 1e-10)
  expect_equal(Sigma, t(Sigma), tolerance = 1e-10) # symmetric
})

test_that("sird_jacobian matches a numerical Jacobian of sird_drift", {
  parms <- list(lambda = 0.4, gamma = 0.15, alpha = 0.08)
  state <- c(s = 0.7, i = 0.2, r = 0.1)
  analytic <- sird_jacobian(state, parms)
  h <- 1e-6
  numeric_J <- matrix(0, 3, 3)
  for (j in 1:3) {
    up <- state; up[j] <- up[j] + h
    dn <- state; dn[j] <- dn[j] - h
    numeric_J[, j] <- (sird_drift(up, parms) - sird_drift(dn, parms)) / (2 * h)
  }
  expect_equal(analytic, numeric_J, tolerance = 1e-5)
})
