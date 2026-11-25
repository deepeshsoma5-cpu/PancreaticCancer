# build gse57495 and gse79668 as gene level cohorts. no scoring here
library(GEOquery)
library(Biobase)
library(edgeR)

tc <- readRDS("data/tcga_paad_pdac.rds")
ref <- rownames(tc$logcpm)
sig <- read.csv("results/10_locked_signature.csv")

clean_sym <- function(x) trimws(sub(" ?(///|//) .*", "", as.character(x)))

# pick the annotation column that matches the most tcga gene names, one row per symbol
to_genes <- function(ex, fd) {
  sc <- sapply(colnames(fd), function(cn) length(intersect(unique(clean_sym(fd[[cn]])), ref)))
  col <- names(which.max(sc))
  cat("annotation column used:", col, "(matches", max(sc), "tcga genes)\n")
  sym <- clean_sym(fd[[col]])
  ok <- !is.na(sym) & sym != "" & sym != "---" & sym != "NA"
  ex <- ex[ok, , drop = FALSE]
  sym <- sym[ok]
  v <- apply(ex, 1, var, na.rm = TRUE)
  v[is.na(v)] <- -Inf
  o <- order(v, decreasing = TRUE)
  ex <- ex[o, , drop = FALSE]
  sym <- sym[o]
  d <- duplicated(sym)
  ex <- ex[!d, , drop = FALSE]
  rownames(ex) <- sym[!d]
  ex
}

report <- function(name, ex, clin) {
  cat("\n=====", name, "=====\n")
  cat("patients:", nrow(clin), " events:", sum(clin$status == 1), " censored:", sum(clin$status == 0), "\n")
  cat("median time (months), events:", round(median(clin$time_months[clin$status == 1]), 1),
      "| censored:", round(median(clin$time_months[clin$status == 0]), 1), "\n")
  cat("genes:", nrow(ex), " shared with tcga:", length(intersect(rownames(ex), ref)), "\n")
  miss <- setdiff(sig$gene, rownames(ex))
  cat("signature genes missing:", if (length(miss)) paste(miss, collapse = ", ") else "none", "\n")
}

# ---- gse57495 ----
g1 <- getGEO("GSE57495", GSEMatrix = TRUE, destdir = "data/")
e1 <- g1[[1]]
pd1 <- pData(e1)
tm1 <- suppressWarnings(as.numeric(as.character(pd1[["overall survival (month):ch1"]])))
st1 <- ifelse(toupper(as.character(pd1[["vital.status:ch1"]])) == "DEAD", 1, 0)
keep1 <- !is.na(tm1) & tm1 > 0
ex1 <- exprs(e1)[, keep1]
if (max(ex1, na.rm = TRUE) > 50) { ex1 <- log2(ex1 + 1); cat("applied log2 to gse57495\n") }
ex1 <- to_genes(ex1, fData(e1))
clin1 <- data.frame(sample = colnames(ex1), time_months = tm1[keep1], status = st1[keep1],
                    stage = as.character(pd1[["Stage:ch1"]][keep1]))
rownames(clin1) <- clin1$sample
saveRDS(list(expr = ex1, clin = clin1), "data/gse57495_tumors.rds")
report("GSE57495", ex1, clin1)

# ---- gse79668 (rna-seq counts from the supplementary file) ----
g2 <- getGEO("GSE79668", GSEMatrix = TRUE, destdir = "data/")
pd2 <- pData(g2[[1]])
sc2 <- grep("survival status", colnames(pd2), value = TRUE, ignore.case = TRUE)[1]
tc2 <- grep("survival time", colnames(pd2), value = TRUE, ignore.case = TRUE)[1]
days <- suppressWarnings(as.numeric(sub(".*: ", "", as.character(pd2[[tc2]]))))
st2 <- ifelse(tolower(sub(".*: ", "", as.character(pd2[[sc2]]))) == "dead", 1, 0)

f <- file.path("data", "GSE79668", "GSE79668_51_tumors_sharedgenecounts.txt.gz")
if (!file.exists(f)) getGEOSuppFiles("GSE79668", baseDir = "data/", makeDirectory = TRUE)
cn <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
cat("\ncounts file dims:", dim(cn), "\n")
print(cn[1:3, 1:min(5, ncol(cn))])
num <- sapply(cn, is.numeric)
cands <- names(cn)[!num]
cat("non numeric columns:", cands, "\n")
score <- sapply(cands, function(cc) length(intersect(unique(clean_sym(cn[[cc]])), ref)))
if (!length(cands) || max(score) < 8000) {
  stop("could not find gene symbols in the counts file, paste me the head printed above")
}
gid <- clean_sym(cn[[names(which.max(score))]])
m <- as.matrix(cn[, num, drop = FALSE])
good <- !is.na(gid) & gid != ""
m <- rowsum(m[good, , drop = FALSE], gid[good])
ord <- match(colnames(m), pd2$title)
if (any(is.na(ord))) {
  cat("sample names did not match the series titles, assuming the same order\n")
  stopifnot(ncol(m) == nrow(pd2))
  ord <- seq_len(ncol(m))
}
days <- days[ord]
st2 <- st2[ord]
keep2 <- !is.na(days) & days > 0 & colSums(m) > 0
m <- m[rowSums(m > 0) >= 5, keep2, drop = FALSE]
y <- calcNormFactors(DGEList(m))
ex2 <- cpm(y, log = TRUE, prior.count = 1)
clin2 <- data.frame(sample = colnames(ex2), time_months = days[keep2] / 30.44, status = st2[keep2])
rownames(clin2) <- clin2$sample
saveRDS(list(expr = ex2, clin = clin2), "data/gse79668_tumors.rds")
report("GSE79668", ex2, clin2)
