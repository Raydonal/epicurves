sird_gm <- function() ga_model(c("S", "I", "R", "D"), list(
  infection = ga_event("S", "I", "lambda", by = "I"),
  recovery  = ga_event("I", "R", "gamma"),
  death     = ga_event("I", "D", "alpha")))

test_that("ga_model builds the stoichiometry matrix and parameter list", {
  m <- sird_gm()
  expect_s3_class(m, "ga_model")
  expect_equal(unname(m$V[, "infection"]), c(-1, 1, 0, 0))
  expect_equal(unname(m$V[, "death"]), c(0, -1, 0, 1))
  expect_equal(m$params, c("lambda", "gamma", "alpha"))
  expect_error(ga_model(c("S", "I"), list(x = ga_event("S", "Q", "k"))), "unknown compartment")
})

test_that("generic engine reproduces the hand-derived SIRD incidence likelihood", {
  set.seed(1)
  N <- 20000
  si <- sird_ga_simulate_inc(N, N - 10, 10, 0, 0.6, 0.3, 0.1, p = 0.5, tau2 = 0.16,
                             q = 1, tau2_d = 0.1, Tmax = 80)
  th <- list(lambda = 0.6, gamma = 0.3, alpha = 0.1, s0 = 1 - 10 / N, i0 = 10 / N, r0 = 0,
             p = 0.5, tau2 = 0.16, q = 1, tau2_d = 0.1)
  obs <- list(ga_incidence("infection", "p", "tau2"), ga_incidence("death", "q", "tau2_d"))
  x0 <- c(S = 1 - 10 / N, I = 10 / N, R = 0, D = 0)
  dat <- cbind(si$cases, si$deaths)
  expect_equal(ga_loglik(sird_gm(), th, x0, dat, obs, N),
               sird_ga_loglik_inc(th, si$cases, si$deaths, N, "both"), tolerance = 1e-9)
  dc <- dat; dc[, 2] <- NA
  expect_equal(ga_loglik(sird_gm(), th, x0, dc, obs, N),
               sird_ga_loglik_inc(th, si$cases, si$deaths, N, "cases"), tolerance = 1e-9)
  lam <- rep(c(0.6, 0.45), each = 40)
  expect_equal(ga_loglik(sird_gm(), th, x0, dat, obs, N, varying = list(lambda = lam)),
               sird_ga_loglik_inc(th, si$cases, si$deaths, N, "both", lambda_t = lam), tolerance = 1e-9)
})

test_that("prevalence observations give a finite likelihood and a mass-conserving mean path", {
  set.seed(2)
  N <- 10000
  th <- list(lambda = 0.5, gamma = 0.2, alpha = 0.05, s0 = 0.999, i0 = 0.001, r0 = 0, p = 0.6, tau2 = 0.3)
  y <- round(seq(5, 400, length.out = 40))
  x0 <- c(S = 0.999, I = 0.001, R = 0, D = 0)
  out <- ga_loglik(sird_gm(), th, x0, matrix(y), list(ga_prevalence("I", "p", "tau2")), N, full = TRUE)
  expect_true(is.finite(out$loglik))
  expect_equal(ncol(out$x), 4L)
  expect_equal(unname(rowSums(out$x)), rep(1, 41), tolerance = 1e-10)
})

test_that("ga_simulate conserves the population and records event counts", {
  set.seed(3)
  m <- sird_gm()
  obs <- list(ga_incidence("infection", 1, 0), ga_prevalence("I", 1, 0))
  s <- ga_simulate(m, list(lambda = 0.6, gamma = 0.3, alpha = 0.1), c(S = 990, I = 10, R = 0, D = 0), obs, n = 30)
  expect_equal(unname(rowSums(s$states)), rep(1000, 31))
  expect_equal(sum(s$latent[, 1]), unname(s$states[1, "S"] - s$states[31, "S"]))
  expect_equal(unname(s$latent[, 2]), unname(s$states[-1, "I"]))
  expect_equal(s$data, s$latent, ignore_attr = TRUE)  # p = 1, no noise
})

test_that("ga_fit recovers SEIRD parameters from simulated cases and deaths", {
  skip_on_cran()
  seird <- ga_model(c("S", "E", "I", "R", "D"), list(
    infection   = ga_event("S", "E", "lambda", by = "I"),
    progression = ga_event("E", "I", "sigma"),
    recovery    = ga_event("I", "R", "gamma"),
    death       = ga_event("I", "D", "alpha")))
  N <- 50000
  obs <- list(ga_incidence("progression", "p", "tau2"), ga_incidence("death", 1, "tau2_d"))
  truth <- list(lambda = 0.5, sigma = 0.25, gamma = 0.2, alpha = 0.05, p = 0.4, tau2 = 0.2, tau2_d = 0.1)
  set.seed(3)
  sim <- ga_simulate(seird, truth, c(S = N - 20, E = 0, I = 20, R = 0, D = 0), obs, n = 120)
  f <- ga_fit(seird, sim$data, obs, N, init = function(p) c(S = 1 - p$i0, E = 0, I = p$i0, R = 0, D = 0),
              start = c(lapply(truth, function(v) v * 1.2), i0 = 20 / N))
  ci <- confint(f)
  for (p in c("lambda", "sigma", "gamma", "alpha", "p"))
    expect_true(ci[p, 1] < truth[[p]] && truth[[p]] < ci[p, 2], info = p)
})
