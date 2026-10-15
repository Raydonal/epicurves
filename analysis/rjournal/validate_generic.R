# Simulation check of the generic Gaussian-approximation engine (ga_model / ga_fit) on models
# that the hand-derived SIRD code does not cover.
# Usage: Rscript validate_generic.R <seird|sirdv> [J] [cores]
suppressMessages(library(epicurves))
args <- commandArgs(trailingOnly = TRUE)
mod <- args[1]; J <- if (length(args) >= 2) as.integer(args[2]) else 100L
cores <- if (length(args) >= 3) as.integer(args[3]) else 4L
N <- 50000; n <- 120; I0 <- 20
if (mod == "seird") {
  model <- ga_model(c("S", "E", "I", "R", "D"), list(
    infection   = ga_event("S", "E", "lambda", by = "I"),
    progression = ga_event("E", "I", "sigma"),
    recovery    = ga_event("I", "R", "gamma"),
    death       = ga_event("I", "D", "alpha")))
  obs <- list(ga_incidence("progression", "p", "tau2"), ga_incidence("death", 1, "tau2_d"))
  truth <- list(lambda = 0.5, sigma = 0.25, gamma = 0.2, alpha = 0.05, p = 0.4, tau2 = 0.2, tau2_d = 0.1)
  init <- function(p) c(S = 1 - p$i0, E = 0, I = p$i0, R = 0, D = 0)
  x0 <- c(S = N - I0, E = 0, I = I0, R = 0, D = 0)
} else {
  model <- ga_model(c("S", "I", "R", "D", "V"), list(
    infection   = ga_event("S", "I", "lambda", by = "I"),
    recovery    = ga_event("I", "R", "gamma"),
    death       = ga_event("I", "D", "alpha"),
    vaccination = ga_event("S", "V", "nu")))
  obs <- list(ga_incidence("infection", "p", "tau2"), ga_incidence("death", 1, "tau2_d"),
              ga_incidence("vaccination", 1, "tau2_v"))
  truth <- list(lambda = 0.45, gamma = 0.2, alpha = 0.05, nu = 0.005, p = 0.5, tau2 = 0.2,
                tau2_d = 0.1, tau2_v = 0.5)
  init <- function(p) c(S = 1 - p$i0, I = p$i0, R = 0, D = 0, V = 0)
  x0 <- c(S = N - I0, I = I0, R = 0, D = 0, V = 0)
}
truth_all <- c(truth, i0 = I0 / N)

one <- function(j) {
  set.seed(1000L * match(mod, c("seird", "sirdv")) + j)
  tries <- 0L
  repeat {
    tries <- tries + 1L
    sim <- ga_simulate(model, truth, x0, obs, n = n)
    if (sum(sim$latent[, 1]) >= 0.05 * N) break   # discard minor outbreaks
  }
  st <- lapply(truth_all, function(v) v * exp(stats::rnorm(1, 0, 0.25)))
  st$p <- min(st$p, 0.95)
  t0 <- proc.time()[3]
  f <- tryCatch(ga_fit(model, sim$data, obs, N, init, start = st), error = function(e) NULL)
  secs <- proc.time()[3] - t0
  if (is.null(f)) return(NULL)
  ci <- confint(f)
  data.frame(rep = j, tries = tries, secs = secs, loglik = f$loglik, param = names(f$coef),
             truth = unlist(truth_all[names(f$coef)]), est = f$coef, lower = ci[, 1], upper = ci[, 2],
             se_ok = !is.null(f$vcov), row.names = NULL)
}
res <- parallel::mclapply(seq_len(J), one, mc.cores = cores, mc.preschedule = FALSE)
res <- do.call(rbind, res[!vapply(res, function(r) is.null(r) || inherits(r, "try-error"), logical(1))])
saveRDS(res, sprintf("results/validate_%s.rds", mod))
cat("done", mod, length(unique(res$rep)), "reps\n")
