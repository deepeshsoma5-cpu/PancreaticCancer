# just looking at gse71729 and gse21501 before deciding anything
library(GEOquery)
library(Biobase)

for (id in c("GSE71729", "GSE21501")) {
  cat("\n\n##########", id, "##########\n")
  g <- getGEO(id, GSEMatrix = TRUE, destdir = "data/")
  for (i in seq_along(g)) {
    e <- g[[i]]
    pd <- pData(e)
    cat("\nplatform:", annotation(e), " samples:", ncol(e), " probes:", nrow(e), "\n")
    cat("title:", experimentData(e)@title, "\n")
    cols <- grep("characteristics|tissue|survival|death|status|stage|time|event|source",
                 colnames(pd), value = TRUE, ignore.case = TRUE)
    cat("\nclinical columns:\n")
    print(cols)
    cat("\nfirst 3 rows:\n")
    print(head(pd[, cols, drop = FALSE], 3))
    for (cc in cols) {
      cat("\n--", cc, "--\n")
      print(head(sort(table(pd[[cc]]), decreasing = TRUE), 6))
    }
  }
}
