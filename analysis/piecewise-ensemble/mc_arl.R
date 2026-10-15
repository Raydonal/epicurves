# Sensitivity of the e-detector (Poisson, z = 1) to the average run length to false alarm.
# Reuses the series of results/mc_<scenario>.rds (same data as Table 1, where arl = 1000) and
# reruns edetect_segment() with arl = 100, 200 and 365; the arl = 1000 results are copied over.
# Usage: Rscript mc_arl.R <scenario: 1c|2c|3c|nc_log|nc_gom> <cores>
suppressMessages(library(epicurves))
args <- commandArgs(trailingOnly = TRUE)
sc <- args[1]; cores <- as.integer(args[2])
dir <- "results"
m <- readRDS(file.path(dir, sprintf("mc_%s.rds", sc)))
arls <- c(100, 200, 365)
timed <- function(expr) { t0 <- proc.time()[3]; v <- tryCatch(expr, error = function(e) NULL); list(v = v, s = proc.time()[3] - t0) }
one <- function(r) {
  out <- list(x = r$x, edetect_arl1000 = r$edetect_z1)
  for (a in arls) {
    e <- timed(edetect_segment(r$x, d_min = 20, arl = a, refit_every = 7, z = 1))
    out[[paste0("edetect_arl", a)]] <- list(cuts = e$v$cuts, alarms = e$v$alarms, secs = e$s)
  }
  out
}
res <- parallel::mclapply(m$res, one, mc.cores = cores, mc.preschedule = FALSE)
saveRDS(list(res = res, truth = m$truth, scenario = sc, arls = c(arls, 1000)),
        file.path(dir, sprintf("mc_arl_%s.rds", sc)))
cat("done", sc, "\n")
