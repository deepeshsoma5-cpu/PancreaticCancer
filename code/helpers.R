# small functions used by the lock and validate scripts
library(survival)
library(glmnet)

cidx <- function(risk, t, s) {
  unname(concordance(Surv(t, s) ~ risk, reverse = TRUE)$concordance)
}

uni_p <- function(x, t, s) {
  apply(x, 2, function(g) {
    f <- tryCatch(suppressWarnings(coxph(Surv(t, s) ~ g)), error = function(e) NULL)
    if (is.null(f)) 1 else summary(f)$coefficients[1, 5]
  })
}

# lasso path, pick the first lambda that gives lo to hi genes
fit_sig <- function(x, y, lo = 5, hi = 10) {
  fit <- glmnet(x, y, family = "cox")
  ix <- which(fit$df >= lo & fit$df <= hi)
  if (!length(ix)) ix <- which(fit$df >= lo)
  li <- if (length(ix)) ix[1] else length(fit$df)
  cf <- as.matrix(coef(fit, s = fit$lambda[li]))
  keep <- cf[, 1] != 0
  list(genes = rownames(cf)[keep], beta = cf[keep, 1])
}

stage_to_num <- function(x) {
  st <- sub("^Stage ", "", x)
  ifelse(grepl("^IV", st), 4,
  ifelse(grepl("^III", st), 3,
  ifelse(grepl("^II", st), 2,
  ifelse(grepl("^I", st), 1, NA))))
}

# kaplan meier plot, df needs time_months, status, grp (TRUE = high risk)
save_km <- function(df, file) {
  suppressMessages(library(survminer))
  km <- survfit(Surv(time_months, status) ~ grp, data = df)
  p <- ggsurvplot(km, data = df, pval = TRUE, risk.table = TRUE,
                  legend.labs = c("low risk", "high risk"), xlab = "months")
  png(file, width = 1600, height = 1300, res = 200, type = "cairo")
  print(p)
  dev.off()
}

logrank_p <- function(df) {
  s <- survdiff(Surv(time_months, status) ~ grp, data = df)
  1 - pchisq(s$chisq, 1)
}

hr_ci <- function(df) {
  s <- summary(coxph(Surv(time_months, status) ~ grp, data = df))
  round(s$conf.int[1, c(1, 3, 4)], 2)
}
