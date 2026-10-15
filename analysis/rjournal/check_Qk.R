# Accuracy of the transition quantities of the general engine as a function of the number of RK4
# steps per observation interval (nsub), for the SEIRD model of example_seird.R at its true values.
# Reference: nsub = 4000. Error of Q_k (and A_k): max over k of ||Q_k - Q_ref||_F / ||Q_ref||_F.
# Also the change of the log-likelihood of the simulated data set at the true values.
suppressMessages(library(epicurves))
ex <- readRDS("results/example_seird.rds")
seird <- ga_model(c("S", "E", "I", "R", "D"), list(
  infection   = ga_event("S", "E", "lambda", by = "I"),
  progression = ga_event("E", "I", "sigma"),
  recovery    = ga_event("I", "R", "gamma"),
  death       = ga_event("I", "D", "alpha")))
obs <- list(cases = ga_incidence("progression", prob = "p", noise = "tau2"),
            deaths = ga_incidence("death", prob = 1, noise = "tau2_d"))
tr <- ex$truth; N <- ex$N; n <- nrow(ex$sim$data)
init <- c(S = N - 20, E = 0, I = 20, R = 0, D = 0) / N
sp <- epicurves:::ga_spec(seird, obs, seird$params)
th <- matrix(unlist(tr[seird$params]), nrow = 1, dimnames = list(NULL, seird$params))
x0 <- c(init, numeric(length(sp$counters)))
run <- function(ns) epicurves:::ga_transition_cpp(n, 1, x0, th, sp$V, sp$pidx, sp$r1, sp$r2, nsub = ns)
relerr <- function(a, b) max(sapply(seq_len(dim(b)[3]), function(k)
  norm(a[, , k] - b[, , k], "F") / norm(b[, , k], "F")))
ref <- run(4000L)
ll <- function(ns) ga_loglik(seird, tr, init, ex$sim$data, obs, N, nsub = ns)
ll_ref <- ll(4000L)
steps <- c(2L, 5L, 10L, 20L, 40L)
tab <- do.call(rbind, lapply(steps, function(ns) { m <- run(ns)
  data.frame(nsub = ns, Q = relerr(m$M, ref$M), A = relerr(m$A, ref$A), dll = ll(ns) - ll_ref) }))
tab$ratio_Q <- c(NA, head(tab$Q, -1) / tab$Q[-1])
print(tab, digits = 3)
saveRDS(tab, "results/check_Qk.rds")
