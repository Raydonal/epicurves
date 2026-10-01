test_that("growth curves satisfy C(0) = y0 and C(t) -> K", {
  y0 <- 10; K <- 1000; gamma <- 0.1; beta <- 1.5
  expect_equal(logistic_curve(0, y0, K, gamma), y0, tolerance = 1e-8)
  expect_equal(gompertz_curve(0, y0, K, gamma), y0, tolerance = 1e-8)
  expect_equal(richards_curve(0, y0, K, gamma, beta), y0, tolerance = 1e-6)

  expect_equal(logistic_curve(500, y0, K, gamma), K, tolerance = 1e-3)
  expect_equal(gompertz_curve(500, y0, K, gamma), K, tolerance = 1e-3)
  expect_equal(richards_curve(500, y0, K, gamma, beta), K, tolerance = 1e-3)
})

test_that("growth curves are monotonically increasing for y0 < K", {
  t <- seq(0, 100, length.out = 50)
  expect_true(all(diff(logistic_curve(t, 10, 1000, 0.1)) >= 0))
  expect_true(all(diff(gompertz_curve(t, 10, 1000, 0.1)) >= 0))
  expect_true(all(diff(richards_curve(t, 10, 1000, 0.1, 1.5)) >= 0))
})

test_that("fit_growth_model recovers known parameters from low-noise data", {
  set.seed(1)
  t <- 1:60
  y_true <- gompertz_curve(t, y0 = 20, K = 2000, gamma = 0.08)
  y <- y_true + rnorm(60, 0, 3)
  fit <- fit_growth_model(t, y, "gompertz")
  expect_equal(unname(fit$coef["K"]), 2000, tolerance = 0.1)
  expect_equal(unname(fit$coef["gamma"]), 0.08, tolerance = 0.3)
  expect_lt(fit$mse, 50)
})

test_that("ensemble weights sum to one and are all non-negative", {
  set.seed(2)
  t <- 1:60
  y_true <- richards_curve(t, y0 = 1374, K = 2200, gamma = 0.03, beta = 1.5)
  y <- y_true + rnorm(60, 0, 5)
  ef <- fit_ensemble(t, y, criterion = "mse", beta_fixed = 1.5)
  expect_equal(sum(ef$weights), 1, tolerance = 1e-8)
  expect_true(all(ef$weights >= 0))
  # With y0/K only ~1.6, Logistic/Gompertz/Richards(beta=1.5) are nearly
  # indistinguishable over this window, so no single model is guaranteed to
  # dominate; we only require that Richards (the true generator) is not
  # driven to a negligible weight by a spurious fit failure.
  expect_gt(ef$weights[["richards"]], 0.15)
})

test_that("ensemble clearly favours the generating model when curve shapes are more distinct", {
  set.seed(3)
  t <- 1:80
  y_true <- richards_curve(t, y0 = 50, K = 5000, gamma = 0.07, beta = 4)
  y <- y_true + rnorm(80, 0, 10)
  ef <- fit_ensemble(t, y, criterion = "mse", beta_fixed = 4)
  expect_equal(unname(which.max(ef$weights)),
               unname(which(names(ef$models) == "richards")))
})
