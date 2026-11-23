options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"), Ncpus = 2)
for (p in c("randomForestSRC", "survminer", "timeROC")) {
  tryCatch(install.packages(p), error = function(e) cat("FAILED:", p, conditionMessage(e), "\n"))
}
for (p in c("GEOquery", "TCGAbiolinks")) {
  tryCatch(BiocManager::install(p, ask = FALSE, update = FALSE), error = function(e) cat("FAILED:", p, conditionMessage(e), "\n"))
}
