test_that("compiled prevalence likelihood agrees with the R implementation", {
  set.seed(1)
  N <- 20000
  sim <- sird_ga_simulate_obs(N, 19980, 20, 0, 0.6, 0.3, 0.1, p = 0.5, tau2 = 0.16, Tmax = 60)
  th <- c(lambda = 0.6, gamma = 0.3, alpha = 0.1, s0 = 0.999, i0 = 0.001, r0 = 0, p = 0.5, tau2 = 0.16)
  expect_equal(sird_ga_loglik_prev(th, sim$y, N, nsub = 20),
               sird_ga_loglik(th, sim$y, N, n_sub = 40), tolerance = 1e-3)
})

test_that("case-only likelihood depends on gamma and alpha only through their sum", {
  set.seed(2)
  N <- 20000
  si <- sird_ga_simulate_inc(N, 19990, 10, 0, 0.6, 0.3, 0.1, p = 0.5, tau2 = 0.16,
                             q = 1, tau2_d = 0.1, Tmax = 80)
  base <- list(lambda = 0.6, s0 = 1 - 10 / N, i0 = 10 / N, r0 = 0, p = 0.5, tau2 = 0.16,
               q = 1, tau2_d = 0.1)
  ll <- function(a, use) {
    th <- base; th$alpha <- a; th$gamma <- 0.4 - a
    sird_ga_loglik_inc(th, si$cases, si$deaths, N, use = use)
  }
  lc <- vapply(c(0.01, 0.1, 0.3), ll, numeric(1), use = "cases")
  expect_lt(diff(range(lc)), 1e-6)
  lb <- vapply(c(0.01, 0.1, 0.3), ll, numeric(1), use = "both")
  expect_gt(diff(range(lb)), 10)
})

test_that("growth_incidence integrates to the curve increment", {
  mu <- growth_incidence(1:400, "logistic", y0 = 10, K = 1000, gamma = 0.1)
  expect_equal(sum(mu), 1000 - 10, tolerance = 1e-3)
  mu <- growth_incidence(1:2000, "richards", y0 = 10, K = 1000, gamma = 0.05, beta = 1.5)
  expect_equal(sum(mu), 1000 - 10, tolerance = 1e-2)
})

test_that("Akaike weights are a probability vector and favour the generating curve", {
  set.seed(3)
  mu <- growth_incidence(1:80, "gompertz", y0 = 20, K = 3000, gamma = 0.05)
  x <- rpois(80, mu)
  e <- fit_ensemble_incidence(x)
  expect_equal(sum(e$weights), 1)
  expect_true(all(e$weights >= 0))
  expect_equal(names(which.max(e$weights)), "gompertz")
})

test_that("edetect_segment runs and returns ordered cuts", {
  skip_if_not_installed("edetect")
  skip_on_cran()
  set.seed(4)
  mu <- c(growth_incidence(1:60, "logistic", 10, 1500, 0.1),
          growth_incidence(1:60, "logistic", 10, 4000, 0.12))
  x <- rpois(120, mu)
  seg <- edetect_segment(x, d_min = 20, arl = 500, z = 1)
  expect_true(all(diff(seg$cuts) > 0))
  expect_equal(nrow(seg$segments), length(seg$cuts) + 1L)
})
