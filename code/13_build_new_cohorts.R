# build gse71729 (moffitt) and gse21501 as gene level cohorts
# tumors with survival only, one row per gene symbol
library(GEOquery)
library(Biobase)

symbol_col <- function(fd) {
  hit <- grep("^gene_?symbol$|^symbol$|gene symbol", colnames(fd), ignore.case = TRUE, value = TRUE)
  if (!length(hit)) stop("no symbol column, columns are: ", paste(colnames(fd), collapse = ", "))
  hit[1]
}

# probes -> genes, keeping the most variable probe for each symbol
to_genes <- function(ex, fd) {
  sym <- trimws(sub(" /// .*", "", as.character(fd[[symbol_col(fd)]])))
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

tc <- readRDS("data/tcga_paad_pdac.rds")

summarize <- function(name, ex, clin) {
  cat("\n=====", name, "=====\n")
  cat("patients:", nrow(clin), " events:", sum(clin$status == 1), " censored:", sum(clin$status == 0), "\n")
  cat("median time, events:", round(median(clin$time_months[clin$status == 1]), 1),
      "| censored:", round(median(clin$time_months[clin$status == 0]), 1), "\n")
  cat("genes:", nrow(ex), " shared with tcga:", length(intersect(rownames(ex), rownames(tc$logcpm))), "\n")
}

# ---- gse71729 ----
g1 <- getGEO("GSE71729", GSEMatrix = TRUE, destdir = "data/")
e1 <- g1[[1]]
pd <- pData(e1)
need <- c("tissue type:ch2", "survival_months:ch2", "death_event_1death_0censor:ch2")
stopifnot(all(need %in% colnames(pd)))
time <- suppressWarnings(as.numeric(pd[["survival_months:ch2"]]))
status <- suppressWarnings(as.numeric(pd[["death_event_1death_0censor:ch2"]]))
keep <- (pd[["tissue type:ch2"]] %in% "Primary") & !is.na(time) & !is.na(status) & time > 0
cat("platform:", annotation(e1), "\n")

ex1 <- exprs(e1)[, keep]
if (max(ex1, na.rm = TRUE) > 50) { ex1 <- log2(ex1 + 1); cat("applied log2\n") }
ex1 <- to_genes(ex1, fData(e1))
clin1 <- data.frame(sample = colnames(ex1), time_months = time[keep], status = status[keep])
sub_col <- "tumor_subtype_0na_1classical_2basal:ch2"
if (sub_col %in% colnames(pd)) clin1$subtype <- pd[[sub_col]][keep]
rownames(clin1) <- clin1$sample
saveRDS(list(expr = ex1, clin = clin1), "data/gse71729_tumors.rds")
summarize("GSE71729 (primary tumors with survival)", ex1, clin1)

# ---- gse21501 ----
g2 <- getGEO("GSE21501", GSEMatrix = TRUE, destdir = "data/")
e2 <- g2[[1]]
pd2 <- pData(e2)
stopifnot(all(c("os time:ch2", "os event:ch2") %in% colnames(pd2)))
time2 <- suppressWarnings(as.numeric(pd2[["os time:ch2"]]))
status2 <- suppressWarnings(as.numeric(pd2[["os event:ch2"]]))
keep2 <- !is.na(time2) & !is.na(status2) & time2 > 0
cat("\nplatform:", annotation(e2), "\n")

ex2 <- exprs(e2)[, keep2]
if (max(ex2, na.rm = TRUE) > 50) { ex2 <- log2(ex2 + 1); cat("applied log2\n") }
ex2 <- to_genes(ex2, fData(e2))
clin2 <- data.frame(sample = colnames(ex2), time_months = time2[keep2], status = status2[keep2])
rownames(clin2) <- clin2$sample
saveRDS(list(expr = ex2, clin = clin2), "data/gse21501_tumors.rds")
summarize("GSE21501 (tumors with survival)", ex2, clin2)
