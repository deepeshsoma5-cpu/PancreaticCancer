# 06_build_gse62452.R: external validation cohort (tumors only, gene-level matrix)
suppressMessages({ library(GEOquery); library(Biobase) })

gse  <- getGEO("GSE62452", GSEMatrix = TRUE, destdir = "data/")
eset <- gse[[1]]
pd   <- pData(eset)
cat("Platform:", annotation(eset), "\n")
cat("\nSeries title:", experimentData(eset)@title, "\n")
cat("\nSeries description (check how survival status is coded):\n")
cat(substr(experimentData(eset)@abstract, 1, 1500), "\n")

# ---- keep tumors with valid survival ----
tumor  <- pd[["tissue:ch1"]] == "Pancreatic tumor"
time   <- suppressWarnings(as.numeric(pd[["survival months:ch1"]]))
status <- suppressWarnings(as.numeric(pd[["survival status:ch1"]]))
keep   <- tumor & !is.na(time) & !is.na(status) & time > 0
cat("\nTumors:", sum(tumor), " | with valid survival:", sum(keep), "\n")

clin <- data.frame(sample = rownames(pd)[keep], time_months = time[keep], status = status[keep],
                   stage = pd[["Stage:ch1"]][keep], stringsAsFactors = FALSE)
if ("grading:ch1" %in% names(pd)) clin$grade <- pd[["grading:ch1"]][keep]
rownames(clin) <- clin$sample

# ---- expression: log2 check, probe -> gene symbol ----
ex <- exprs(eset)[, clin$sample]
if (max(ex, na.rm = TRUE) > 50) { ex <- log2(ex + 1); cat("Applied log2 transform\n") }

fd  <- fData(eset)
sym <- NULL
if ("gene_assignment" %in% names(fd)) {
  sym <- sapply(strsplit(as.character(fd$gene_assignment), " // "), function(x) if (length(x) >= 2) x[2] else NA)
} else {
  cand <- grep("symbol", names(fd), ignore.case = TRUE, value = TRUE)
  if (length(cand)) sym <- as.character(fd[[cand[1]]])
}
if (is.null(sym)) stop("No gene symbol column found. Columns: ", paste(names(fd), collapse = ", "))
sym <- trimws(sub(" /// .*", "", sym))

ok  <- !is.na(sym) & sym != "" & sym != "---"
ex  <- ex[ok, ]; sym <- sym[ok]
ord <- order(rowMeans(ex), decreasing = TRUE)
ex  <- ex[ord, ]; sym <- sym[ord]
dup <- duplicated(sym)
ex  <- ex[!dup, ]; rownames(ex) <- sym[!dup]

saveRDS(list(expr = ex, clin = clin), "data/gse62452_tumors.rds")

# ---- summary ----
tc     <- readRDS("data/tcga_paad_pdac.rds")
shared <- intersect(rownames(ex), rownames(tc$logcpm))
cat("\n===== EXTERNAL COHORT (GSE62452, tumors) =====\n")
cat("Patients (n):", nrow(clin), "\n")
cat("Status 1:", sum(clin$status == 1), " Status 0:", sum(clin$status == 0), "\n")
cat("Median months, status 1:", round(median(clin$time_months[clin$status == 1]), 1),
    "| status 0:", round(median(clin$time_months[clin$status == 0]), 1), "\n")
cat("Genes (unique symbols):", nrow(ex), "\n")
cat("Genes shared with TCGA:", length(shared), "\n")
cat("\nStage:\n"); print(table(clin$stage))
