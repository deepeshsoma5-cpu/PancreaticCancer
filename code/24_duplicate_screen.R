# duplicate screen with the thresholds committed in docs/plan_fresh_cohort.md
# r >= 0.8 alone, or r >= 0.5 with survival time within 0.5 months and the same status
tc <- readRDS("data/tcga_paad_pdac.rds")
sig <- read.csv("results/10_locked_signature.csv")
get <- function(f) {
  x <- readRDS(f)
  list(expr = x$expr, clin = x$clin[colnames(x$expr), ])
}
refs <- list(
  tcga = list(expr = tc$logcpm, clin = tc$clin[colnames(tc$logcpm), ]),
  gse62452 = get("data/gse62452_tumors.rds"),
  gse21501 = get("data/gse21501_tumors.rds"),
  gse71729 = get("data/gse71729_tumors.rds"))
news <- list(gse57495 = get("data/gse57495_tumors.rds"), gse79668 = get("data/gse79668_tumors.rds"))

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
  data.frame(sample = rownames(r), match = colnames(r)[apply(r, 1, which.max)], best_r = best,
             match_surv = colnames(r)[apply(rs, 1, which.max)], best_r_survmatch = bests,
             flagged = best >= 0.8 | bests >= 0.5, stringsAsFactors = FALSE)
}

for (nm in names(news)) {
  N <- news[[nm]]
  others <- c(refs, news[setdiff(names(news), nm)])
  flagged <- character(0)
  allres <- list()
  cat("\n##########", nm, "##########\n")
  for (rn in names(others)) {
    s <- screen(N, others[[rn]])
    cat(sprintf("vs %-9s | best r: median %.2f, max %.2f | flagged %d of %d\n",
                rn, median(s$best_r), max(s$best_r), sum(s$flagged), nrow(s)))
    if (any(s$flagged)) {
      print(head(s[s$flagged, c("sample", "match", "best_r", "best_r_survmatch")], 5), row.names = FALSE)
    }
    flagged <- union(flagged, s$sample[s$flagged])
    s$against <- rn
    allres[[rn]] <- s
  }
  write.csv(do.call(rbind, allres), paste0("results/24_screen_", nm, ".csv"), row.names = FALSE)
  writeLines(flagged, paste0("results/24_dropped_", nm, ".txt"))

  keep <- setdiff(colnames(N$expr), flagged)
  clin <- N$clin[keep, ]
  miss <- setdiff(sig$gene, rownames(N$expr))
  ok <- nrow(clin) >= 50 && sum(clin$status) >= 30 && length(miss) == 0
  cat(sprintf("\n%s after screen: %d patients, %d events, missing signature genes: %s -> %s\n",
              nm, nrow(clin), sum(clin$status), if (length(miss)) paste(miss, collapse = ", ") else "none",
              if (ok) "PASSES the inclusion rule" else "FAILS the inclusion rule, skip and document"))
  saveRDS(list(expr = N$expr[, keep, drop = FALSE], clin = clin, passes = ok),
          paste0("data/fresh_", nm, ".rds"))
}
