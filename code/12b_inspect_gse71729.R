# just gse71729 this time, short output
library(GEOquery)
library(Biobase)
g <- getGEO("GSE71729", GSEMatrix = TRUE, destdir = "data/")
for (i in seq_along(g)) {
  e <- g[[i]]
  pd <- pData(e)
  cat("\nplatform:", annotation(e), " samples:", ncol(e), "\n")
  cols <- grep("characteristics|tissue|survival|death|status|stage|time|event|grade|subtype",
               colnames(pd), value = TRUE, ignore.case = TRUE)
  cols <- setdiff(cols, c("status"))
  cat("\nclinical columns:\n")
  print(cols)
  for (cc in cols) {
    cat("\n--", cc, "--\n")
    print(head(sort(table(pd[[cc]]), decreasing = TRUE), 6))
  }
}
