# overlap final: drop every moffitt sample whose best gse21501 match has identical time + status
# also check gse62452 vs the other two (secondary test must not share patients with training)
mo <- readRDS("data/gse71729_tumors.rds")
st <- readRDS("data/gse21501_tumors.rds")

best_match <- function(a, b) {
  g <- intersect(rownames(a$expr), rownames(b$expr))
  r <- cor(a$expr[g, ], b$expr[g, ], use = "pairwise.complete.obs", method = "spearman")
  j <- apply(r, 1, which.max)
  data.frame(a = rownames(r), b = colnames(r)[j], r = round(r[cbind(seq_len(nrow(r)), j)], 3),
             same_fu = a$clin[rownames(r), "time_months"] == b$clin[colnames(r)[j], "time_months"] &
                       a$clin[rownames(r), "status"] == b$clin[colnames(r)[j], "status"])
}

d <- best_match(mo, st)
drop <- d$a[d$same_fu]
writeLines(drop, "results/overlap_drop_gse71729.txt")
cat("dropping", length(drop), "of", nrow(mo$clin), "gse71729 samples; left for training:", nrow(mo$clin) - length(drop), "\n")
cat("events left:", sum(mo$clin[setdiff(rownames(mo$clin), drop), "status"]), "\n")

f <- "data/gse62452_tumors.rds"
if (file.exists(f)) {
  ex <- readRDS(f)
  for (nm in c("mo", "st")) {
    o <- best_match(ex, get(nm))
    cat("\ngse62452 vs", nm, ": identical-followup best matches =", sum(o$same_fu), "of", nrow(o), "\n")
    print(head(o[o$same_fu, ][order(-o$r[o$same_fu]), ], 10))
  }
} else cat("\nfile", f, "not found, check the ls output above and edit the name\n")
