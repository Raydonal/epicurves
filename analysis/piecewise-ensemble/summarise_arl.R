# Summary of the run-length sensitivity (results/mc_arl_*.rds), scored as in summarise_mc.R (tol = 7).
# Writes results/mc_arl_summary.rds; printed output goes to results/mc_arl_summary.log via run_arl.sh.
tol <- 7
dir <- "results"
summ <- function(res, tc, rule) {
  per <- lapply(res, function(r) {
    o <- r[[rule]]; ec <- if (is.null(o$cuts)) numeric() else o$cuts
    det <- vapply(tc, function(c) any(abs(ec - c) <= tol), logical(1))
    false <- sum(vapply(ec, function(e) !any(abs(tc - e) <= tol), logical(1)))
    delay <- if (length(tc)) vapply(tc, function(c) {
      k <- which(abs(ec - c) <= tol); if (length(k)) o$alarms[k[1]] - c else NA_real_ }, numeric(1)) else NA_real_
    list(failed = is.null(o$cuts), n = length(ec), det = det, false = false, delay = delay)
  })
  n <- vapply(per, `[[`, numeric(1), "n")
  data.frame(rule = rule, failed = sum(vapply(per, `[[`, logical(1), "failed")), mean_cuts = mean(n),
             any_false = mean(vapply(per, `[[`, numeric(1), "false") > 0),
             false_per_rep = mean(vapply(per, `[[`, numeric(1), "false")),
             detect = if (length(tc)) mean(unlist(lapply(per, `[[`, "det"))) else NA,
             median_delay = suppressWarnings(median(unlist(lapply(per, `[[`, "delay")), na.rm = TRUE)),
             mean_delay = suppressWarnings(mean(unlist(lapply(per, `[[`, "delay")), na.rm = TRUE)))
}
all <- list()
for (sc in c("nc_log", "nc_gom", "1c", "2c", "3c")) {
  f <- file.path(dir, sprintf("mc_arl_%s.rds", sc)); if (!file.exists(f)) next
  m <- readRDS(f)
  s <- do.call(rbind, lapply(paste0("edetect_arl", m$arls), summ, res = m$res, tc = m$truth$cuts))
  s <- cbind(scenario = sc, arl = m$arls, s); all[[sc]] <- s
  cat(sprintf("\n== scenario %s, J = %d, true cuts: %s\n", sc, length(m$res), paste(m$truth$cuts, collapse = " ")))
  print(s[, -1], digits = 3, row.names = FALSE)
}
saveRDS(do.call(rbind, all), file.path(dir, "mc_arl_summary.rds"))
