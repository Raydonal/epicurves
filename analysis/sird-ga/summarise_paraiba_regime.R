# Summary of sim_paraiba_regime.R, per scenario (A, B, C):
#  - multimodality: share of replicates whose three generic starts differ by > 0.5 in loglik,
#    median and 90th percentile of the spread of p-hat among the starts, share where the best start
#    is not unique within 0.5; convergence codes; share of positive-definite Hessians;
#  - unconstrained parameters: bias, Monte Carlo SD, median Wald SE, coverage of par +- 1.96 SE;
#  - Table 3 quantities: bias, relative bias, RMSE, Monte Carlo SD, median SD of the 4000 draws,
#    coverage of the draw-quantile interval of paraiba_report.R (truth inside [2.5%, 97.5%]) and of
#    estimate +- 1.96 SD(draws). Truth: deterministic path at the true parameters (for C, the
#    q = 1 pseudo-truth); attack rate and block R_t also against the simulated path.
# Usage: Rscript summarise_paraiba_regime.R [dir]   -> results/paraiba_regime_summary.{rds,log via stdout}
args <- commandArgs(TRUE)
dir <- if (length(args) >= 1) args[1] else "results"
options(width = 200)
out <- list()
for (sc in c("A", "B", "C")) {
  files <- list.files(dir, sprintf("^regime_%s_[0-9]+\\.rds$", sc), full.names = TRUE)
  if (!length(files)) next
  tr <- readRDS(file.path(dir, sprintf("regime_%s_truth.rds", sc)))
  reps <- unlist(lapply(files, readRDS), recursive = FALSE)
  reps <- Filter(function(z) is.list(z) && isTRUE(z$ok), reps)
  R <- length(reps)
  cat(sprintf("\n================ Scenario %s: p = %.3f, q = %.2f, replicates with a fit = %d\n",
              sc, tr$p, tr$q, R))
  # multimodality
  spread_ll <- sapply(reps, function(z) diff(range(z$starts$loglik, na.rm = TRUE)))
  spread_p <- sapply(reps, function(z) diff(range(z$starts$p, na.rm = TRUE)))
  nfail <- sapply(reps, function(z) sum(is.na(z$starts$loglik)))
  conv <- sapply(reps, `[[`, "conv"); pd <- sapply(reps, `[[`, "pd")
  multi <- data.frame(share_ll_spread_gt_0.5 = mean(spread_ll > 0.5), median_ll_spread = median(spread_ll),
                      share_p_spread_gt_0.01 = mean(spread_p > 0.01), median_p_spread = median(spread_p),
                      q90_p_spread = unname(quantile(spread_p, 0.9)), failed_starts = sum(nfail),
                      share_conv0 = mean(conv == 0), share_hessian_pd = mean(pd),
                      median_secs = median(sapply(reps, `[[`, "secs")))
  cat("\nMultimodality / numerics\n"); print(format(multi, digits = 3), row.names = FALSE)
  cat("convergence codes of best fit:", paste(names(table(conv)), table(conv), sep = ":", collapse = " "), "\n")
  # among replicates where starts disagree, p-hat of each start
  dis <- which(spread_ll > 0.5)
  if (length(dis)) {
    cat("replicates with ll spread > 0.5 (first 10): start p-hat / loglik\n")
    for (i in head(dis, 10)) cat(sprintf("  rep %d: p = %s | ll = %s\n", reps[[i]]$b,
      paste(sprintf("%.3f", reps[[i]]$starts$p), collapse = " "),
      paste(sprintf("%.2f", reps[[i]]$starts$loglik), collapse = " ")))
  }
  # unconstrained parameters
  P <- t(sapply(reps, `[[`, "par")); SE <- t(sapply(reps, `[[`, "se_par"))
  ok <- pd
  un <- do.call(rbind, lapply(seq_along(tr$par_true), function(j) {
    t0 <- tr$par_true[[j]]; e <- P[, j]; s <- SE[ok, j]
    data.frame(param = names(tr$par_true)[j], true = t0, mean = mean(e), bias = mean(e) - t0,
               mc_sd = sd(e), median_se = median(s, na.rm = TRUE),
               cover95 = mean(abs(P[ok, j] - t0) <= 1.96 * s, na.rm = TRUE), n = sum(is.finite(s)))
  }))
  cat("\nUnconstrained parameters (Wald, PD Hessians only for SE/coverage)\n")
  print(format(un, digits = 3), row.names = FALSE)
  # derived quantities
  E <- t(sapply(reps, `[[`, "est")); Lo <- t(sapply(reps, `[[`, "lower"))
  Up <- t(sapply(reps, `[[`, "upper")); SD <- t(sapply(reps, `[[`, "sd"))
  qrow <- function(name, t0, j) {
    e <- E[, j]; okj <- ok & is.finite(Lo[, j])
    data.frame(quantity = name, true = t0, mean = mean(e), median = median(e), bias = mean(e) - t0,
               rel_bias = (mean(e) - t0) / t0, rmse = sqrt(mean((e - t0)^2)), mc_sd = sd(e),
               median_draw_sd = median(SD[okj, j]),
               cover_quantile = mean(Lo[okj, j] <= t0 & t0 <= Up[okj, j]),
               cover_sym = mean(abs(E[okj, j] - t0) <= 1.96 * SD[okj, j]),
               below = mean(Up[okj, j] < t0), above = mean(Lo[okj, j] > t0), n = sum(okj))
  }
  qn <- colnames(E)
  der <- do.call(rbind, lapply(seq_along(qn), function(j) qrow(qn[j], tr$derived_true[[qn[j]]], j)))
  att_sim <- sapply(reps, `[[`, "attack_sim"); Rsim <- t(sapply(reps, `[[`, "R_block_sim"))
  ja <- which(qn == "attack")
  extra <- {
    e <- E[, ja]; okj <- ok & is.finite(Lo[, ja])
    data.frame(quantity = "attack_vs_simulated", true = mean(att_sim), mean = mean(e), median = median(e),
               bias = mean(e - att_sim), rel_bias = mean(e - att_sim) / mean(att_sim),
               rmse = sqrt(mean((e - att_sim)^2)), mc_sd = sd(e - att_sim),
               median_draw_sd = median(SD[okj, ja]),
               cover_quantile = mean(Lo[okj, ja] <= att_sim[okj] & att_sim[okj] <= Up[okj, ja]),
               cover_sym = mean(abs(e[okj] - att_sim[okj]) <= 1.96 * SD[okj, ja]),
               below = mean(Up[okj, ja] < att_sim[okj]), above = mean(Lo[okj, ja] > att_sim[okj]), n = sum(okj))
  }
  Rext <- do.call(rbind, lapply(grep("^R_block", qn), function(j) {
    t0 <- Rsim[, qn[j]]; e <- E[, j]; okj <- ok & is.finite(Lo[, j])
    data.frame(quantity = paste0(qn[j], "_vs_simulated"), true = mean(t0), mean = mean(e), median = median(e),
               bias = mean(e - t0), rel_bias = mean(e - t0) / mean(t0), rmse = sqrt(mean((e - t0)^2)),
               mc_sd = sd(e - t0), median_draw_sd = median(SD[okj, j]),
               cover_quantile = mean(Lo[okj, j] <= t0[okj] & t0[okj] <= Up[okj, j]),
               cover_sym = mean(abs(e[okj] - t0[okj]) <= 1.96 * SD[okj, j]),
               below = mean(Up[okj, j] < t0[okj]), above = mean(Lo[okj, j] > t0[okj]), n = sum(okj))
  }))
  der <- rbind(der, extra, Rext)
  if (sc == "C") {
    jn <- which(qn == "ifr")
    der <- rbind(der, qrow("ifr_vs_nominal_(q=0.8_truth)", tr$derived_nominal[["ifr"]], jn))
  }
  cat("\nTable 3 quantities (truth = deterministic path", if (sc == "C") "of the q = 1 pseudo-truth" else "", ")\n")
  print(format(der, digits = 3), row.names = FALSE)
  out[[sc]] <- list(R = R, multimodality = multi, unconstrained = un, derived = der,
                    spread_ll = spread_ll, spread_p = spread_p)
}
saveRDS(out, file.path(dir, "paraiba_regime_summary.rds"))
