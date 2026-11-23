# quick check on two things
# 1. are the genes that keep getting picked actually expressed much?
# 2. how well do stage + age alone predict survival? (clinical baseline)
library(survival)
set.seed(42)

tc <- readRDS("data/tcga_paad_pdac.rds")
S <- read.csv("results/07_cv_signature_genes.csv")

# --- 1. expression level of the most picked genes ---
tab <- sort(table(S$gene), decreasing = TRUE)
top <- names(tab)[1:15]
top <- top[top %in% rownames(tc$logcpm)]
avg <- rowMeans(tc$logcpm)
pct <- ecdf(avg)(avg[top]) * 100
out <- data.frame(gene = top, times_picked = as.integer(tab[top]),
                  mean_logcpm = round(avg[top], 2), expr_percentile = round(pct))
cat("\n===== most picked genes and how expressed they are =====\n")
print(out, row.names = FALSE)
cat("\nmedian expression percentile of these genes:", median(out$expr_percentile), "\n")

# --- 2. clinical only baseline (stage + age), same kind of cv ---
cl <- tc$clin
st <- sub("^Stage ", "", cl$ajcc_pathologic_stage)
cl$stage_num <- ifelse(grepl("^IV", st), 4,
                ifelse(grepl("^III", st), 3,
                ifelse(grepl("^II", st), 2,
                ifelse(grepl("^I", st), 1, NA))))
cl$age <- cl$age_at_diagnosis / 365.25
cl <- cl[complete.cases(cl[, c("time_months", "status", "stage_num", "age")]), ]
cat("\npatients with stage and age:", nrow(cl), "\n")

cidx <- function(risk, t, s) {
  unname(concordance(Surv(t, s) ~ risk, reverse = TRUE)$concordance)
}

scores <- c()
for (r in 1:10) {
  f <- integer(nrow(cl))
  for (v in unique(cl$status)) {
    i <- which(cl$status == v)
    f[i] <- sample(rep_len(1:5, length(i)))
  }
  for (k in 1:5) {
    tr <- cl[f != k, ]
    te <- cl[f == k, ]
    m <- coxph(Surv(time_months, status) ~ stage_num + age, data = tr)
    risk <- predict(m, newdata = te, type = "lp")
    scores <- c(scores, cidx(risk, te$time_months, te$status))
  }
}
cat(sprintf("\n===== clinical baseline (stage + age) =====\nmean c-index %.3f | sd %.3f\n",
            mean(scores), sd(scores)))
