sir_boot_setup <- function() {
  sir <- ga_model(c("S", "I", "R"), list(
    infection = ga_event("S", "I", "lambda", by = "I"),
    removal   = ga_event("I", "R", "gamma")))
  obs <- list(cases = ga_incidence("infection", prob = "p", noise = "tau2"),
              prev = ga_prevalence("R", prob = 1, noise = 0.1))
  set.seed(11)
  sim <- ga_simulate(sir, list(lambda = 0.5, gamma = 0.2, p = 0.5, tau2 = 0.1),
                     c(S = 1990, I = 10, R = 0), obs, n = 40)
  dat <- sim$data; dat[, 2] <- NA
  init <- function(p) c(S = 1 - p$i0, I = p$i0, R = 0)
  fit <- ga_fit(sir, dat, obs, 2000, init, fixed = list(tau2 = 0.1),
                start = list(lambda = 0.4, gamma = 0.25, i0 = 0.005, p = 0.4), maxit = 300)
  list(model = sir, obs = obs, fit = fit)
}

test_that("ga_simulate with constant varying rates reproduces the constant simulation", {
  s <- sir_boot_setup()
  th <- list(lambda = 0.5, gamma = 0.2, p = 0.5, tau2 = 0.1)
  set.seed(3); a <- ga_simulate(s$model, th, c(S = 990, I = 10, R = 0), s$obs, n = 20)
  set.seed(3); b <- ga_simulate(s$model, th["gamma" != names(th)], c(S = 990, I = 10, R = 0), s$obs, n = 20,
                                varying = list(gamma = rep(0.2, 20)))
  expect_identical(a$data, b$data)
})

test_that("ga_bootstrap is reproducible, keeps missing values and gives intervals", {
  s <- sir_boot_setup()
  stat <- function(f) c(coef(f), R0 = f$theta$lambda / f$theta$gamma)
  b1 <- ga_bootstrap(s$fit, B = 4, statistic = stat, seed = 5, maxit = 300)
  b2 <- ga_bootstrap(s$fit, B = 4, statistic = stat, seed = 5, maxit = 300)
  expect_s3_class(b1, "ga_boot")
  expect_identical(b1$t, b2$t)
  skip_on_os("windows")
  b3 <- ga_bootstrap(s$fit, B = 4, statistic = stat, seed = 5, maxit = 300, cores = 2)
  expect_identical(b1$t, b3$t)
  expect_equal(dim(b1$t), c(4L, 5L))
  expect_equal(colnames(b1$t), names(stat(s$fit)))
  expect_equal(sum(b1$init_counts), 2000)
  expect_true(all(is.finite(b1$t)))
  ci <- confint(b1)
  expect_true(all(ci[, 1] <= ci[, 2]))
  cb <- confint(b1, "R0", type = "basic")
  expect_equal(unname(cb[1, ]), unname(2 * b1$t0[["R0"]] - rev(confint(b1, "R0")[1, ])))
  expect_error(ga_bootstrap(s$fit, B = 1, init_counts = c(S = 990, I = 10, R = 0)), "sum")
})

test_that("ga_bootstrap handles block-wise rates", {
  s <- sir_boot_setup()
  blk <- rep(1:2, each = 20)
  set.seed(21)
  sim <- ga_simulate(s$model, list(gamma = 0.2, p = 0.5, tau2 = 0.1), c(S = 1990, I = 10, R = 0),
                     s$obs[1], n = 40, varying = list(lambda = c(0.6, 0.3)[blk]))
  init <- function(p) c(S = 1 - p$i0, I = p$i0, R = 0)
  fit <- ga_fit(s$model, sim$data, s$obs[1], 2000, init, fixed = list(tau2 = 0.1), blocks = list(lambda = blk),
                start = list(lambda = c(0.5, 0.5), gamma = 0.25, i0 = 0.005, p = 0.4), maxit = 300)
  b <- ga_bootstrap(fit, B = 3, seed = 1, maxit = 300)
  expect_equal(colnames(b$t), c("lambda[1]", "lambda[2]", "gamma", "i0", "p"))
  expect_true(all(is.finite(b$t)))
})
