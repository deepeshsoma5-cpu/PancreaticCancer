pkgs <- c("survival","glmnet","randomForestSRC","survminer","timeROC","TCGAbiolinks","GEOquery","limma","edgeR")
for (p in pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("%-18s %s\n", p, ifelse(ok, "OK", "MISSING")))
}
