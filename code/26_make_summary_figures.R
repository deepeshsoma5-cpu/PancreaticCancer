cv <- read.csv("results/07_cv_results.csv")
nn <- read.csv("results/09_nn_cv_results.csv")
scores <- read.csv("results/25_fresh_scores.csv")

cv_labels <- c("Cox-LASSO", "Random survival forest", "Neural network")
cv_mean <- c(
  mean(cv$cindex[cv$model == "Cox-LASSO"]),
  mean(cv$cindex[cv$model == "RandomSurvivalForest"]),
  mean(nn$cindex)
)
cv_sd <- c(
  sd(cv$cindex[cv$model == "Cox-LASSO"]),
  sd(cv$cindex[cv$model == "RandomSurvivalForest"]),
  sd(nn$cindex)
)

png("results/26_cv_comparison.png", width = 1800, height = 1200, res = 200)
par(mar = c(5, 10, 3, 2), cex.lab = 1.1, cex.axis = 1.05)
y <- rev(seq_along(cv_labels))
plot(cv_mean, y, xlim = c(0.35, 0.85), ylim = c(0.5, 3.5), yaxt = "n",
     xlab = "Mean cross-validated concordance index (C-index)",
     ylab = "", pch = 19, cex = 1.5, col = "#1f5a85",
     main = "Model comparison in TCGA-PAAD (10 repeated 5-fold CV)")
axis(2, at = y, labels = cv_labels, las = 1)
segments(cv_mean - cv_sd, y, cv_mean + cv_sd, y, lwd = 2, col = "#6b7280")
segments(cv_mean - cv_sd, y - 0.08, cv_mean - cv_sd, y + 0.08, col = "#6b7280")
segments(cv_mean + cv_sd, y - 0.08, cv_mean + cv_sd, y + 0.08, col = "#6b7280")
abline(v = 0.5, lty = 2, col = "#9ca3af")
text(cv_mean, y + 0.18, sprintf("%.3f +/- %.3f SD", cv_mean, cv_sd),
     cex = 0.85, pos = 4, offset = 0.4)
dev.off()

keep <- scores$cohort != "GSE57495 (all 63, sensitivity)"
val <- scores[keep, ]
val$label <- c("GSE62452 (previously viewed)", "GSE21501 (fresh)",
               "GSE57495 (fresh; screened)")
colors <- c("#667085", "#1f5a85", "#1f5a85")
png("results/26_validation_forest.png", width = 2000, height = 1500, res = 200)
layout(matrix(c(1, 2), ncol = 1), heights = c(1, 1))
par(mar = c(4, 12, 3, 2), cex.lab = 1.05, cex.axis = 0.95)
y <- rev(seq_len(nrow(val)))
plot(val$cindex, y, xlim = c(0.35, 0.75), ylim = c(0.5, nrow(val) + 0.5),
     yaxt = "n", xlab = "C-index (95% bootstrap percentile interval)",
     ylab = "", pch = 19, cex = 1.4, col = colors,
     main = "Locked 8-gene signature in external cohorts")
axis(2, at = y, labels = paste0(val$label, "  (n=", val$n, ", events=", val$events, ")"),
     las = 1)
segments(val$ci_lo, y, val$ci_hi, y, lwd = 2, col = colors)
segments(val$ci_lo, y - 0.08, val$ci_lo, y + 0.08, col = colors)
segments(val$ci_hi, y - 0.08, val$ci_hi, y + 0.08, col = colors)
abline(v = 0.5, lty = 2, col = "#9ca3af")
par(mar = c(5, 12, 2, 2), cex.lab = 1.05, cex.axis = 0.95)
plot(val$hr_sd, y, log = "x", xlim = c(0.65, 2.5),
     ylim = c(0.5, nrow(val) + 0.5), yaxt = "n",
     xlab = "Hazard ratio per 1 SD higher risk score (95% Cox CI)",
     ylab = "", pch = 19, cex = 1.4, col = colors)
axis(2, at = y, labels = rep("", nrow(val)))
segments(val$hr_lo, y, val$hr_hi, y, lwd = 2, col = colors)
segments(val$hr_lo, y - 0.08, val$hr_lo, y + 0.08, col = colors)
segments(val$hr_hi, y - 0.08, val$hr_hi, y + 0.08, col = colors)
abline(v = 1, lty = 2, col = "#9ca3af")
dev.off()
