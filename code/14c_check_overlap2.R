# overlap check v2: corr threshold was too strict (platform diffs), so also use mutual best match + identical follow-up
# no signature or model performance used here
mo <- readRDS("data/gse71729_tumors.rds")
st <- readRDS("data/gse21501_tumors.rds")
g <- intersect(rownames(mo$expr), rownames(st$expr))
r <- cor(mo$expr[g, ], st$expr[g, ], use = "pairwise.complete.obs", method = "spearman")

j <- apply(r, 1, which.max)          # best gse21501 for each moffitt
k <- apply(r, 2, which.max)          # best moffitt for each gse21501
best <- r[cbind(seq_len(nrow(r)), j)]
sec <- apply(r, 1, function(x) sort(x, decreasing = TRUE)[2])
mutual <- k[j] == seq_len(nrow(r))

t_mo <- mo$clin[rownames(r), "time_months"]; s_mo <- mo$clin[rownames(r), "status"]
t_st <- st$clin[colnames(r)[j], "time_months"]; s_st <- st$clin[colnames(r)[j], "status"]
same <- t_mo == t_st & s_mo == s_st

d <- data.frame(moffitt = rownames(r), match = colnames(r)[j], r = round(best, 3),
                gap = round(best - sec, 3), mutual = mutual, same_fu = same)
cat("mutual best matches:", sum(mutual), "\n")
cat("identical time+status among best matches:", sum(same), "of", nrow(d), "\n")
cat("identical time+status AND mutual:", sum(same & mutual), "\n")

# chance level: how often does a random gse21501 patient share time+status with a moffitt patient?
set.seed(1)
nul <- replicate(2000, {
  x <- sample(nrow(st$clin), nrow(d), replace = TRUE)
  mean(t_mo == st$clin$time_months[x] & s_mo == st$clin$status[x])
})
cat("expected fraction by chance:", round(mean(nul), 4), " observed:", round(mean(same), 3), "\n")

drop <- d$moffitt[same & mutual]
cat("\nflagged as same patient:", length(drop), "\n")
print(head(d[same & mutual, ][order(-d$r[same & mutual]), ], 15))
writeLines(drop, "results/overlap_drop_gse71729.txt")
cat("\nunflagged samples that still have identical fu (look at these by eye):\n")
print(d[same & !mutual, ])
