# metadata only: sample counts, platform, and which fields might hold survival. nothing is scored
suppressMessages(library(GEOquery))
for (id in c("GSE57495", "GSE78229", "GSE79668", "GSE28735")) {
  cat("\n#####", id, "\n")
  g <- tryCatch(getGEO(id, GSEMatrix = TRUE, getGPL = FALSE), error = function(e) { cat("failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(g)) next
  for (k in seq_along(g)) {
    p <- pData(g[[k]])
    cat("platform:", annotation(g[[k]]), "| samples:", ncol(exprs(g[[k]])), "| genes/probes:", nrow(exprs(g[[k]])), "\n")
    ch <- grep("characteristics|title|source", names(p), value = TRUE)
    print(head(p[, ch, drop = FALSE], 3))
  }
}
