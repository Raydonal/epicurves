# Growth-curve part of the R Journal article: daily COVID-19 deaths in Paraiba (2020) divided into
# phases by offline penalised likelihood and by the sequential e-detector, with negative binomial
# ensembles per phase.
suppressMessages(library(epicurves))
data(paraiba_covid19)
d <- paraiba_covid19
x <- d$new_deaths
t0 <- proc.time()[3]
seg <- bic_segment(x, d_min = 12, family = "negbin", step = 5)
secs_bic <- proc.time()[3] - t0
t0 <- proc.time()[3]
det <- edetect_segment(x, d_min = 20, arl = 1000, refit_every = 7, chart = "subgaussian", z = 1)
secs_det <- proc.time()[3] - t0
bounds <- seg$segments
ens <- lapply(seq_len(nrow(bounds)), function(j) fit_ensemble_incidence(x[bounds$start[j]:bounds$end[j]], family = "negbin"))
saveRDS(list(x = x, dates = d$date, seg = seg, det = det, ens = ens, secs_bic = secs_bic, secs_det = secs_det),
        "results/example_paraiba.rds")
cat("penalised likelihood cuts:", format(d$date[seg$cuts]), " (", round(secs_bic), "s)\n")
cat("e-detector cuts:", format(d$date[det$cuts]), " alarms:", format(d$date[det$alarms]), " (", round(secs_det), "s)\n")
for (j in seq_along(ens)) { cat("phase", j, ":", format(d$date[bounds$start[j]]), "-", format(d$date[bounds$end[j]]), " weights "); print(round(ens[[j]]$weights, 3)) }
