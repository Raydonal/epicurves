# Figure: profile log-likelihood of alpha, cases only vs cases and deaths (results/profile_alpha_sim.rds).
r <- readRDS("results/profile_alpha_sim.rds")$profile
figdir <- if (dir.exists("../CSDA_submission/figures")) "../CSDA_submission/figures" else "figures"
dir.create(figdir, showWarnings = FALSE)
pdf(file.path(figdir, "profile_alpha_sim.pdf"), width = 5.5, height = 3.6)
par(mar = c(4, 4.5, 0.5, 0.5), las = 1)
both <- r[r$data == "cases and deaths", ]; cas <- r[r$data == "cases only", ]
plot(both$value, both$loglik - both$max, type = "b", pch = 19, ylim = c(-30, 1),
     xlab = expression(alpha~"(fixed)"), ylab = "Profile log-likelihood minus maximum")
lines(cas$value, cas$loglik - cas$max, type = "b", pch = 1, lty = 2)
abline(v = 0.1, col = "grey50", lty = 3); abline(h = -1.92, col = "grey50", lty = 3)
legend("bottomright", c("Cases and deaths", "Cases only"), pch = c(19, 1), lty = c(1, 2), bty = "n")
dev.off()
print(cbind(both[, c("value")], round(both$loglik - both$max, 2), round(cas$loglik - cas$max, 3)))
