# all pairs duplicate audit across the training cohort and the scored cohorts
# thresholds are the ones committed in docs/plan_fresh_cohort.md
# z score inside each cohort, 3000 most variable shared genes
# flagged if best r >= 0.8, or r >= 0.5 with survival time within 0.5 months and the same status
tc <- readRDS("data/tcga_paad_pdac.rds")
rd <- function(f) {
  x <- readRDS(f)
  list(expr = x$expr, clin = x$clin[colnames(x$expr), ])
}
co <- list(
  TCGA = list(expr = tc$logcpm, clin = tc$clin[colnames(tc$logcpm), ]),
  GSE62452 = rd("data/gse62452_tumors.rds"),
  GSE21501 = rd("data/gse21501_tumors.rds"),
  GSE57495 = rd("data/gse57495_tumors.rds"),
  GSE71729 = rd("data/gse71729_tumors.rds"))

screen <- function(N, R) {
  g <- intersect(rownames(N$expr), rownames(R$expr))
  rk <- (rank(-apply(N$expr[g, ], 1, mad)) + rank(-apply(R$expr[g, ], 1, mad))) / 2
  g <- names(sort(rk))[1:min(3000, length(g))]
  zN <- scale(t(N$expr[g, ]))
  zR <- scale(t(R$expr[g, ]))
  ok <- apply(zN, 2, function(v) all(is.finite(v))) & apply(zR, 2, function(v) all(is.finite(v)))
  r <- cor(t(zN[, ok]), t(zR[, ok]))
  same <- outer(N$clin$time_months, R$clin$time_months, function(a, b) abs(a - b) <= 0.5) &
          outer(N$clin$status, R$clin$status, "==")
  rs <- r
  rs[!same] <- -1
  best <- apply(r, 1, max)
  bests <- apply(rs, 1, max)
  data.frame(sample = rownames(r), match = colnames(r)[apply(r, 1, which.max)],
             best_r = best, best_r_survmatch = bests,
             flagged = best >= 0.8 | bests >= 0.5, stringsAsFactors = FALSE)
}

nm <- names(co)
flags <- list()
cnt <- matrix(NA_integer_, length(nm), length(nm), dimnames = list(screened = nm, against = nm))
for (a in nm) {
  for (b in nm) {
    if (a == b) next
    s <- screen(co[[a]], co[[b]])
    cnt[a, b] <- sum(s$flagged)
    if (any(s$flagged)) {
      f <- s[s$flagged, c("sample", "match", "best_r", "best_r_survmatch")]
      f$cohort <- a
      f$against <- b
      flags[[paste(a, b)]] <- f
    }
  }
}
cat("patients flagged (rows) when screened against each other cohort (columns):\n")
print(cnt)
fl <- if (length(flags)) do.call(rbind, flags) else data.frame()
write.csv(fl, "results/41_audit_flags_all.csv", row.names = FALSE)

# who would be dropped from the three scored cohorts in the sensitivity rescore
scored <- c("GSE62452", "GSE21501", "GSE57495")
drops <- data.frame(cohort = character(), sample = character(), against = character(),
                    best_r = numeric(), best_r_survmatch = numeric(), stringsAsFactors = FALSE)
if (nrow(fl)) {
  for (i in seq_along(scored)) {
    x <- scored[i]
    earlier <- c("TCGA", scored[seq_len(i - 1)])
    d <- fl[fl$cohort == x & fl$against %in% earlier, , drop = FALSE]
    if (nrow(d)) {
      drops <- rbind(drops, data.frame(cohort = x, sample = d$sample, against = d$against,
                                       best_r = d$best_r, best_r_survmatch = d$best_r_survmatch,
                                       stringsAsFactors = FALSE))
    }
  }
}
if (nrow(drops)) drops <- drops[!duplicated(paste(drops$cohort, drops$sample)), ]
write.csv(drops, "results/41_audit_drops.csv", row.names = FALSE)
cat("\nscored-cohort patients flagged against TCGA or an earlier scored cohort:", nrow(drops), "\n")
if (nrow(drops)) print(drops, row.names = FALSE)
