# LaTeX tables of the SIRD-GA simulation study from results/sim_s*_N*.rds (uses summarise_sim.R).
source("summarise_sim.R")
sc_order <- c(1, 2, 3, 4, 5); Ns <- c(5000, 20000, 100000)
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
rows_a <- rows_b <- character()
disc <- c()
for (sc in sc_order) for (N in Ns) {
  f <- sprintf("results/sim_s%d_N%d.rds", sc, N)
  if (!file.exists(f)) next
  x <- readRDS(f); s <- summ(f); g <- function(p, col) s[s$param == p, col]
  lab <- if (N == Ns[1]) sprintf("\\multirow{3}{*}{%d}", sc) else ""
  Nlab <- gsub(",", "\\\\,", formatC(N, format = "d", big.mark = ","))
  pars <- c("lambda", "gamma", "alpha", "p", "R0", "ifr")
  rows_a <- c(rows_a, paste0(lab, " & ", Nlab, " & ",
    paste(sapply(pars, function(p) sprintf("%s (%s)", fmt(g(p, "mean")), fmt(g(p, "rmse")))), collapse = " & "), " \\\\"))
  cov <- sapply(c("lambda", "gamma", "alpha", "p"), function(p) fmt(g(p, "cover95"), 2))
  pdisc <- 100 * sum(x$tries - 1) / sum(x$tries); disc <- c(disc, setNames(pdisc, paste(sc, N)))
  rows_b <- c(rows_b, paste0(lab, " & ", Nlab, " & ", paste(cov, collapse = " & "), " & ",
    g("lambda", "n_se"), " & ", fmt(pdisc, 1), " & ", round(median(x$secs)), " \\\\"))
  if (N == Ns[3] && sc != tail(sc_order, 1)) { rows_a <- c(rows_a, "\\midrule"); rows_b <- c(rows_b, "\\midrule") }
}
truth_line <- function(sc) {
  x <- readRDS(sprintf("results/sim_s%d_N5000.rds", sc)); t <- attr(x, "truth")
  sprintf("%d: $\\lambda=%g$, $\\gamma=%g$, $\\alpha=%g$, $p=%g$, $R_0=%s$, IFR $=%s$", sc, t[["lambda"]], t[["gamma"]],
          t[["alpha"]], t[["p"]], fmt(t[["lambda"]] / (t[["gamma"]] + t[["alpha"]]), 2),
          fmt(t[["alpha"]] / (t[["gamma"]] + t[["alpha"]]), 3))
}
cat("% Table A rows: mean (RMSE)\n", paste(rows_a, collapse = "\n"), "\n\n% Table B rows: coverage, n with SE, % discarded, median secs\n",
    paste(rows_b, collapse = "\n"), "\n\n% true values\n", paste(sapply(sc_order[file.exists(sprintf("results/sim_s%d_N5000.rds", sc_order))], truth_line), collapse = "\n"), "\n", sep = "")
cat(sprintf("\n%% overall discarded draws: %.2f%% (range %.1f-%.1f%% across cells)\n",
            100 * sum(sapply(list.files("results", "^sim_s.*rds$", full.names = TRUE), function(f) sum(readRDS(f)$tries - 1))) /
              sum(sapply(list.files("results", "^sim_s.*rds$", full.names = TRUE), function(f) sum(readRDS(f)$tries))),
            min(disc), max(disc)))
