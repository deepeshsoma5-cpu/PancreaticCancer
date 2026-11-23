options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"), Ncpus = 2)
for (p in c("mvtnorm", "plyr", "pec", "timereg", "timeROC")) {
  tryCatch(install.packages(p, dependencies = TRUE), error = function(e) cat("FAILED:", p, conditionMessage(e), "\n"))
}
cat("timeROC installed:", requireNamespace("timeROC", quietly = TRUE), "\n")
