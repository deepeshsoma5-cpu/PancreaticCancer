# finish the inspection that got cut off (gse57495), then check gse79668 events and files
library(GEOquery)
library(Biobase)

g <- getGEO("GSE57495", GSEMatrix = TRUE, destdir = "data/")
e <- g[[1]]
pd <- pData(e)
cat("\n##### GSE57495 | platform:", annotation(e), "| samples:", ncol(e), "| probes:", nrow(e), "\n")
cat("annotation columns:", paste(colnames(fData(e)), collapse = ", "), "\n")
cols <- grep("characteristics|survival|status|stage|tissue|source", colnames(pd), value = TRUE, ignore.case = TRUE)
for (cc in cols) {
  cat("\n--", cc, "--\n")
  print(head(sort(table(pd[[cc]]), decreasing = TRUE), 8))
}

g2 <- getGEO("GSE79668", GSEMatrix = TRUE, destdir = "data/")
pd2 <- pData(g2[[1]])
sc <- grep("survival status", colnames(pd2), value = TRUE, ignore.case = TRUE)[1]
tcol <- grep("survival time", colnames(pd2), value = TRUE, ignore.case = TRUE)[1]
cat("\n##### GSE79668 | columns used:", sc, "|", tcol, "\n")
st <- tolower(sub(".*: ", "", as.character(pd2[[sc]])))
days <- suppressWarnings(as.numeric(sub(".*: ", "", as.character(pd2[[tcol]]))))
cat("status counts:\n")
print(table(st))
cat("events (dead):", sum(st == "dead"), " censored:", sum(st != "dead"), "\n")
cat("median survival time (months):", round(median(days / 30.44, na.rm = TRUE), 1), "\n")
cat("patients with a time:", sum(!is.na(days)), "\n")

cat("\nsupplementary files listed for GSE79668:\n")
sf <- tryCatch(getGEOSuppFiles("GSE79668", baseDir = "data/", fetch_files = FALSE), error = function(e) NULL)
print(sf)
