# LaTeX rows for the Monte Carlo table of the piecewise ensemble article, from results/mc_summary_tol7.rds
# (run summarise_mc.R first). Columns: mean number of cuts, share of series with the true number,
# share of true change points detected within 7 days, false cuts per series, median detection delay
# (detectors only), median seconds per series.
s <- readRDS("results/mc_summary_tol7.rds")
lab_rule <- c(thesis_flat = "Proportional, cumulative", cum_bic = "BIC, cumulative",
              inc_bic = "Penalised likelihood", edetect_z0 = "E-detector, $z=0$",
              edetect_z1 = "E-detector, $z=1$")
lab_sc <- c(nc_log = "Logistic, no change", nc_gom = "Gompertz, no change",
            `1c` = "1, continuous (3)", `2c` = "2, continuous (4)", `3c` = "3, continuous (5)",
            `1` = "1, with spikes (3)", `2` = "2, with spikes (4)", `3` = "3, with spikes (5)")
f <- function(x, d) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))
out <- character()
for (sc in names(lab_sc)) {
  d <- s[s$scenario == sc, ]
  if (!nrow(d)) next
  for (i in seq_len(nrow(d))) {
    first <- if (i == 1) sprintf("\\multirow{5}{*}{\\parbox{2.6cm}{%s}}", lab_sc[[sc]]) else ""
    out <- c(out, paste0(first, " & ", lab_rule[[d$rule[i]]], " & ", f(d$mean_cuts[i], 2), " & ",
                         f(d$exact[i], 2), " & ", f(d$detect[i], 2), " & ", f(d$false_per_rep[i], 2), " & ",
                         f(d$median_delay[i], 0), " & ", f(d$median_secs[i], 0), "\\\\"))
  }
  out <- c(out, "\\midrule")
}
out <- out[-length(out)]
cat(out, sep = "\n")
cat(sprintf("\n%% failures: %d\n", sum(s$failed)))
