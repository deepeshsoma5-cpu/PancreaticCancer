# are any gse71729 patients also in gse21501? same lab, same reference rna, so worth checking
# only looks at expression, no outcome is used to decide anything
mo <- readRDS("data/gse71729_tumors.rds")
st <- readRDS("data/gse21501_tumors.rds")
g <- intersect(rownames(mo$expr), rownames(st$expr))
cat("genes in common:", length(g), "\n")

r <- cor(mo$expr[g, ], st$expr[g, ], use = "pairwise.complete.obs")
best <- apply(r, 1, max)
j <- apply(r, 1, which.max)

cat("\nbest match correlation for each gse71729 sample:\n")
print(round(quantile(best, c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)), 3))
for (th in c(0.8, 0.9, 0.95)) cat("samples with a match >=", th, ":", sum(best >= th), "\n")

pairs <- data.frame(moffitt = rownames(r), match = colnames(r)[j], r = round(best, 3),
                    t_mo = mo$clin[rownames(r), "time_months"], s_mo = mo$clin[rownames(r), "status"],
                    t_st = st$clin[colnames(r)[j], "time_months"], s_st = st$clin[colnames(r)[j], "status"])
pairs <- pairs[order(-pairs$r), ]
cat("\ntop 10 pairs (do the survival times and status match?):\n")
print(head(pairs, 10), row.names = FALSE)

dup <- names(best)[best >= 0.9]
writeLines(dup, "results/14_dropped_moffitt_samples.txt")
cat("\ndropping", length(dup), "gse71729 samples from training\n")
