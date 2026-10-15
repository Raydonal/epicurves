# Cumulative curves of the three simulation scenarios (thesis, Ch. 3).
# continuous = FALSE reproduces the thesis construction: each phase restarts its curve from a
# prescribed y0, so the cumulative curve jumps at every change point and the daily mean has a
# one-day spike there. continuous = TRUE keeps the same curve families, K and gamma but starts
# each phase from the cumulative value at the end of the previous one (no jump, no spike).
phase_spec <- function(sc) switch(as.character(sc),
  `1` = list(y0 = 10, n = c(70, 25, 40, 65), model = c("L", "G", "R", "G"),
             y0s = c(10, 812, 1374, 1973), K = c(1200, 1500, 2200, 2600),
             gamma = c(0.08, 0.07, 0.03, 0.02)),
  `2` = list(y0 = 20, n = c(38, 56, 27, 25, 54), model = c("G", "R", "L", "G", "G"),
             y0s = c(20, 590, 8670, 17430, 22261), K = c(2440, 14040, 24500, 24560, 25520),
             gamma = c(0.03, 0.054, 0.055, 0.05, 0.047)),
  `3` = list(y0 = 26, n = c(40, 49, 39, 17, 21, 34), model = c("G", "L", "L", "R", "R", "G"),
             y0s = c(26, 1210, 10259, 20907, 23423, 24956),
             K = c(1685, 19337, 25790, 25600, 25600, 25620),
             gamma = c(0.05, 0.056, 0.047, 0.032, 0.033, 0.046)))
scenario_curve <- function(sc, continuous = FALSE) {
  p <- phase_spec(sc)
  crv <- function(m, u, y0, K, g) switch(m, L = logistic_curve(u, y0 = y0, K = K, gamma = g),
                                         G = gompertz_curve(u, y0 = y0, K = K, gamma = g),
                                         R = richards_curve(u, y0 = y0, K = K, gamma = g, beta = 1.5))
  C <- numeric(); prev <- p$y0
  for (j in seq_along(p$n)) {
    start <- if (continuous) prev else p$y0s[j]
    cc <- crv(p$model[j], 1:p$n[j], start, p$K[j], p$gamma[j])
    C <- c(C, cc); prev <- cc[length(cc)]
  }
  cuts <- cumsum(p$n)[-length(p$n)]
  # daily mean: first difference of the cumulative curve, with C(0) = y0
  list(C = C, mu = pmax(diff(c(p$y0, C)), 0.01), cuts = cuts, y0 = p$y0)
}
