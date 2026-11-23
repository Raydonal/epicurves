# Summary of sim_generic_starts.R: the summaries of summarise_sim.R for the generic-start fits,
# side by side with the original perturbed-truth fits of the same replicates, and the agreement
# between the two (same data checked by ll_check_diff).
# Usage: Rscript summarise_generic_starts.R [dir]
source("summarise_sim.R")
args <- commandArgs(TRUE)
dir <- if (length(args) >= 1) args[1] else "results"
options(width = 200)
out <- list()
for (f in sort(list.files(dir, "^sim_generic_s[0-9]+_N[0-9]+\\.rds$", full.names = TRUE))) {
  g <- readRDS(f)
  key <- sub("^sim_generic_(s[0-9]+_N[0-9]+)\\.rds$", "\\1", basename(f))
  o <- readRDS(file.path("results", sprintf("sim_%s.rds", key)))
  gok <- g[!g$failed & is.finite(g$loglik), ]
  tmp <- tempfile(fileext = ".rds"); attributes_keep <- attributes(g)[c("truth", "scenario")]
  gg <- gok; attr(gg, "truth") <- attributes_keep$truth; attr(gg, "scenario") <- attributes_keep$scenario
  saveRDS(gg, tmp)
  sg <- summ(tmp)
  so <- summ(file.path("results", sprintf("sim_%s.rds", key)))
  cat(sprintf("\n================ %s: generic starts, %d of %d replicates fitted (failed %d); original: %d\n",
              key, nrow(gok), nrow(g), sum(g$failed), nrow(o)))
  cat("same-data check, |loglik at original estimate on regenerated data - stored loglik|: max =",
      signif(max(abs(g$ll_check_diff), na.rm = TRUE), 3), " (n compared =", sum(is.finite(g$ll_check_diff)), ")\n")
  cat("\nGeneric starts (best of three)\n"); print(sg, digits = 3)
  cat("\nOriginal single start from perturbed truth\n"); print(so, digits = 3)
  m <- merge(gok, o, by = "rep", suffixes = c("_g", "_o"))
  dll <- m$loglik_g - m$loglik_o
  sp <- apply(gok[, paste0("ll_start", 1:3)], 1, function(v) diff(range(v, na.rm = TRUE)))
  ps <- apply(gok[, paste0("p_start", 1:3)], 1, function(v) diff(range(v, na.rm = TRUE)))
  agree <- data.frame(n_common = nrow(m), median_dll = median(dll),
                      share_abs_dll_lt_0.01 = mean(abs(dll) < 0.01),
                      share_generic_better_gt_0.5 = mean(dll > 0.5),
                      share_original_better_gt_0.5 = mean(dll < -0.5),
                      max_abs_diff_p = max(abs(m$p_g - m$p_o)),
                      share_starts_ll_spread_gt_0.5 = mean(sp > 0.5),
                      median_p_spread_starts = median(ps), share_hessian_pd = mean(gok$hessian_pd),
                      median_secs = median(gok$secs))
  cat("\nAgreement generic vs original (dll = loglik generic - loglik original)\n")
  print(format(agree, digits = 3), row.names = FALSE)
  out[[key]] <- list(generic = sg, original = so, agreement = agree, dll = dll)
}
saveRDS(out, file.path(dir, "sim_generic_summary.rds"))
