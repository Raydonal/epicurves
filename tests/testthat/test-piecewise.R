simulate_scenario1 <- function(seed = 123) {
  set.seed(seed)
  C <- numeric(200)
  C[1:70]    <- logistic_curve(1:70, y0 = 10,   K = 1200, gamma = 0.08)
  C[71:95]   <- gompertz_curve(1:25, y0 = 812,  K = 1500, gamma = 0.07)
  C[96:135]  <- richards_curve(1:40, y0 = 1374, K = 2200, gamma = 0.03, beta = 1.5)
  C[136:200] <- gompertz_curve(1:65, y0 = 1973, K = 2600, gamma = 0.02)
  lambda_t <- pmax(c(C[1], diff(C)), 0.01)
  cumsum(stats::rpois(200, lambda_t))
}

test_that("fit_piecewise_ensemble recovers the manuscript's Scenario 1 change points", {
  skip_on_cran()
  y <- simulate_scenario1()
  pw <- fit_piecewise_ensemble(y, window = 90, update_every = 20, min_obs = 12,
                                criterion = "mse", beta_fixed = 1.5, threshold = 0.5)
  # Thesis reports d1=72, d2=95, d3=132 against true values 70, 95, 135
  # (within a few days of truth); we require the same order of accuracy.
  true_cuts <- c(70, 95, 135)
  first_three <- pw$cuts[seq_len(min(3, length(pw$cuts)))]
  expect_equal(length(first_three), 3L)
  expect_true(all(abs(first_three - true_cuts) <= 5))
})

test_that("the AIC/BIC split criterion is at least as conservative as the raw MSE threshold", {
  skip_on_cran()
  y <- simulate_scenario1(seed = 123)
  pw_mse <- fit_piecewise_ensemble(y, min_obs = 12, confirm_criterion = "mse_threshold")
  pw_aic <- fit_piecewise_ensemble(y, min_obs = 12, confirm_criterion = "aic")
  # The BIC-style penalty on the doubled parameter count should never
  # confirm more splits than the flat 50%-MSE-reduction rule.
  expect_lte(length(pw_aic$cuts), length(pw_mse$cuts))
  # It must still find the three genuine change points.
  expect_true(all(abs(pw_aic$cuts[1:3] - c(70, 95, 135)) <= 5))
})

test_that("predict.piecewise_fit returns a value for every requested time", {
  y <- simulate_scenario1()
  pw <- fit_piecewise_ensemble(y, min_obs = 12)
  pred <- predict(pw, newdata = c(1, 50, 100, 150, 200))
  expect_length(pred, 5)
  expect_true(all(is.finite(pred)))
})
