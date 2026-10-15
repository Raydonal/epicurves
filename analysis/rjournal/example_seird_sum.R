# Complement to example_seird.R: cases-only fits of the SEIRD example with the removal rates
# fixed, (i) two splits of the same sum gamma + alpha = 0.25 and (ii) a grid of sums.
suppressMessages(library(epicurves))
r <- readRDS("results/example_seird.rds")
fc <- r$fit_cases
seird <- fc$call_args$model; obs <- fc$call_args$obs; init <- fc$call_args$init
cases_only <- r$sim$data; cases_only[, 2] <- NA
st <- list(lambda = 0.5, sigma = 0.25, i0 = 4e-4, p = 0.4, tau2 = 0.2)
fitfix <- function(g, a) ga_fit(seird, cases_only, obs, r$N, init, start = st,
                                fixed = list(gamma = g, alpha = a, tau2_d = 0.1), hessian = FALSE)
split1 <- fitfix(0.20, 0.05); split2 <- fitfix(0.10, 0.15)
sums <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
prof_sum <- data.frame(sum = sums, loglik = sapply(sums, function(s) fitfix(0.8 * s, 0.2 * s)$loglik))
out <- list(split1 = split1$loglik, split2 = split2$loglik,
            split1_coef = coef(split1), split2_coef = coef(split2), prof_sum = prof_sum,
            max_cases = fc$loglik)
saveRDS(out, "results/example_seird_sum.rds")
print(out)
