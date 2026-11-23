# 05_inspect_gse62452.R: download the external cohort and look at its clinical columns
suppressMessages(library(GEOquery))
gse <- getGEO("GSE62452", GSEMatrix = TRUE, destdir = "data/")
eset <- gse[[1]]
pd <- Biobase::pData(eset)
cat("\n===== GSE62452 =====\n")
cat("Samples:", ncol(eset), " Features (probes):", nrow(eset), "\n")
cat("Platform:", Biobase::annotation(eset), "\n")
cat("\nClinical columns:\n"); print(colnames(pd))
cols <- grep("characteristics|tissue|survival|death|status|stage", colnames(pd), value = TRUE, ignore.case = TRUE)
cat("\nFirst 5 rows of the relevant columns:\n")
print(head(pd[, cols, drop = FALSE], 5))
for (cc in cols) { cat("\n--", cc, "--\n"); print(head(sort(table(pd[[cc]]), decreasing = TRUE), 8)) }
