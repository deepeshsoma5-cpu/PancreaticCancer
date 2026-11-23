# build gse71729 (moffitt) and gse21501 as gene level cohorts
# tumors with survival only, one row per gene symbol
# the symbol column is picked by testing every annotation column against tcga gene names
library(GEOquery)
library(Biobase)

tc <- readRDS("data/tcga_paad_pdac.rds")
ref <- rownames(tc$logcpm)

get_geo <- function(id) {
  f <- file.path("data", paste0(id, "_series_matrix.txt.gz"))
  for (i in 1:5) {
    r <- tryCatch({
      if (file.exists(f)) getGEO(filename = f, destdir = "data/", GSEMatrix = TRUE)
      else getGEO(id, GSEMatrix = TRUE, destdir = "data/")
    }, error = function(e) { cat("try", i, "failed for", id, ":", conditionMessage(e), "\n"); NULL })
    if (!is.null(r)) {
      if (inherits(r, "ExpressionSet")) return(r)
      return(r[[1]])
    }
    Sys.sleep(10)
  }
  stop("could not get ", id, " after 5 tries")
}

clean_sym <- function(x) trimws(sub(" ?(///|//) .*", "", as.character(x)))

# which annotation column has the most gene names that tcga also has?
symbol_col <- function(fd) {
  sc <- sapply(colnames(fd), function(cn) length(intersect(unique(clean_sym(fd[[cn]])), ref)))
  cat("annotation columns and how many tcga genes they match:\n")
  print(sc)
  if (max(sc) < 5000) stop("no annotation column matches tcga gene names well")
  names(which.max(sc))
}

# probes -> genes, keeping the most variable probe for each symbol
to_genes <- function(ex, fd) {
  col <- symbol_col(fd)
  sym <- clean_sym(fd[[col]])
  cat("using annotation column:", col, "| examples:", paste(head(sym[sym != ""], 6), collapse = ", "), "\n")
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

summarize <- function(name, ex, clin) {
  shared <- length(intersect(rownames(ex), ref))
  cat("\n=====", name, "=====\n")
  cat("patients:", nrow(clin), " events:", sum(clin$status == 1), " censored:", sum(clin$status == 0), "\n")
  cat("median time, events:", round(median(clin$time_months[clin$status == 1]), 1),
      "| censored:", round(median(clin$time_months[clin$status == 0]), 1), "\n")
  cat("genes:", nrow(ex), " shared with tcga:", shared, "\n")
  if (shared < 8000) cat("WARNING: very few shared genes, the gene names probably didnt map right\n")
}

# ---- gse71729 ----
e1 <- get_geo("GSE71729")
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
e2 <- get_geo("GSE21501")
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

# ---- all four cohorts together ----
e5 <- readRDS("data/gse62452_tumors.rds")
all4 <- Reduce(intersect, list(rownames(tc$logcpm), rownames(ex1), rownames(ex2), rownames(e5$expr)))
cat("\ngenes present in all four cohorts:", length(all4), "\n")
