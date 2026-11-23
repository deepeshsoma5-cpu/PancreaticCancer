# 07_model_selection.R
# Repeated 5-fold CV on TCGA-PAAD ONLY. Compares Cox-LASSO vs random survival forest.
# All supervised steps (gene filter, LASSO, forest importance) happen inside each training fold.
suppressMessages({ library(survival); library(glmnet); library(randomForestSRC) })
set.seed(42)
N_REPEATS <- as.integer(Sys.getenv("N_REPEATS", "3"))
K <- 5; N_VAR <- 2000; P_CUT <- 0.01

# ---- 1. Load, align, restrict to shared genes ----
tc <- readRDS("data/tcga_paad_pdac.rds")
ex <- readRDS("data/gse62452_tumors.rds")
shared <- intersect(rownames(tc$logcpm), rownames(ex$expr))

# Unsupervised variance filter (no survival info used): top N_VAR genes by MAD in TCGA
mads <- apply(tc$logcpm[shared, ], 1, mad)
top  <- names(sort(mads, decreasing = TRUE))[1:N_VAR]

# z-score each gene within its own cohort
Xt <- scale(t(tc$logcpm[top, ]))
Xe <- scale(t(ex$expr[top, ]))
ok <- apply(Xt, 2, function(v) all(is.finite(v))) & apply(Xe, 2, function(v) all(is.finite(v)))
Xt <- Xt[, ok]; Xe <- Xe[, ok]

# formula-safe gene names (keep a lookup to the real symbols)
orig <- colnames(Xt)
safe <- make.names(orig, unique = TRUE)
colnames(Xt) <- safe; colnames(Xe) <- safe
gene_map <- setNames(orig, safe)

clin <- tc$clin[match(rownames(Xt), rownames(tc$clin)), ]
time <- clin$time_months; status <- clin$status
saveRDS(list(Xt = Xt, Xe = Xe, clin_t = clin, clin_e = ex$clin, gene_map = gene_map),
        "data/model_input.rds")
cat("Training patients:", nrow(Xt), "| candidate genes:", ncol(Xt), "\n")

# ---- 2. Helpers ----
cidx <- function(risk, t, s) {
  unname(concordance(Surv(t, s) ~ risk, reverse = TRUE)$concordance)
}
uni_p <- function(X, t, s) {
  apply(X, 2, function(g) {
    f <- tryCatch(suppressWarnings(coxph(Surv(t, s) ~ g)), error = function(e) NULL)
    if (is.null(f)) 1 else summary(f)$coefficients[1, 5]
  })
}
make_folds <- function(s, K) {
  f <- integer(length(s))
  for (v in unique(s)) { i <- which(s == v); f[i] <- sample(rep_len(1:K, length(i))) }
  f
}

# ---- 3. Repeated cross-validation ----
res <- list(); sigs <- list(); n <- 0
for (r in 1:N_REPEATS) {
  fold <- make_folds(status, K)
  for (k in 1:K) {
    tr <- which(fold != k); te <- which(fold == k)
    Xtr <- Xt[tr, ]; Xte <- Xt[te, ]
    ytr <- Surv(time[tr], status[tr])

    # supervised gene filter, training fold only
    p <- uni_p(Xtr, time[tr], status[tr])
    sel <- names(p)[p < P_CUT]
    if (length(sel) < 20) sel <- names(sort(p))[1:20]

    # --- Cox-LASSO, aiming for a 5-10 gene signature ---
    fit <- glmnet(Xtr[, sel], ytr, family = "cox")
    ix <- which(fit$df >= 5 & fit$df <= 10)
    if (!length(ix)) ix <- which(fit$df >= 5)
    li <- if (length(ix)) ix[1] else length(fit$df)
    lam <- fit$lambda[li]
    risk_l <- as.numeric(predict(fit, newx = Xte[, sel], s = lam, type = "link"))
    cf <- coef(fit, s = lam); g_l <- rownames(cf)[as.numeric(cf) != 0]
    n <- n + 1
    res[[n]] <- data.frame(rep = r, fold = k, model = "Cox-LASSO",
                           cindex = cidx(risk_l, time[te], status[te]), n_genes = length(g_l))
    sigs[[n]] <- data.frame(model = "Cox-LASSO", gene = g_l)

    # --- Random survival forest: importance on top 100 filtered genes, refit on top 10 ---
    sel_rf <- head(names(sort(p[sel])), 100)
    d0 <- data.frame(time = time[tr], status = status[tr], Xtr[, sel_rf, drop = FALSE])
    rf0 <- rfsrc(Surv(time, status) ~ ., data = d0, ntree = 300, importance = "permute", seed = -42)
    g_r <- names(sort(rf0$importance, decreasing = TRUE))[1:10]
    d1 <- data.frame(time = time[tr], status = status[tr], Xtr[, g_r, drop = FALSE])
    rf1 <- rfsrc(Surv(time, status) ~ ., data = d1, ntree = 300, seed = -42)
    risk_r <- predict(rf1, newdata = data.frame(Xte[, g_r, drop = FALSE]))$predicted
    n <- n + 1
    res[[n]] <- data.frame(rep = r, fold = k, model = "RandomSurvivalForest",
                           cindex = cidx(risk_r, time[te], status[te]), n_genes = length(g_r))
    sigs[[n]] <- data.frame(model = "RandomSurvivalForest", gene = g_r)
  }
  cat("Finished repeat", r, "of", N_REPEATS, "\n")
}

# ---- 4. Results ----
R <- do.call(rbind, res); S <- do.call(rbind, sigs)
write.csv(R, "results/07_cv_results.csv", row.names = FALSE)
S$gene <- gene_map[S$gene]
write.csv(S, "results/07_cv_signature_genes.csv", row.names = FALSE)

cat("\n===== CROSS-VALIDATED C-INDEX (TCGA only) =====\n")
for (m in unique(R$model)) {
  v <- R$cindex[R$model == m]
  cat(sprintf("%-22s mean %.3f | sd %.3f | median genes %d\n", m, mean(v), sd(v),
              median(R$n_genes[R$model == m])))
}
w <- reshape(R[, c("rep", "fold", "model", "cindex")], idvar = c("rep", "fold"),
             timevar = "model", direction = "wide")
d <- w[["cindex.Cox-LASSO"]] - w[["cindex.RandomSurvivalForest"]]
cat(sprintf("\nMean difference (LASSO - forest): %.3f | LASSO won %d of %d folds\n",
            mean(d), sum(d > 0), length(d)))
cat("\nMost frequently selected genes (stability):\n")
tab <- sort(table(S$gene), decreasing = TRUE)
print(head(tab, 15))
