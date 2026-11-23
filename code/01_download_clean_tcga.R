# 01_download_clean_tcga.R
# Downloads TCGA-PAAD, keeps one primary tumor per patient, removes non-ductal
# histologies, builds the survival outcome, filters low-expressed genes.

suppressMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(edgeR)
})
set.seed(42)

# ---- 1. Download RNA-seq counts (about 185 samples) ----
query <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query, method = "api", files.per.chunk = 20)
se <- GDCprepare(query)

cd <- as.data.frame(colData(se))
cat("\nSamples downloaded:", ncol(se), "\n")
cat("\nSample types:\n");  print(table(cd$sample_type))
cat("\nDiagnoses (before filtering):\n"); print(table(cd$primary_diagnosis))

# ---- 2. Keep primary tumors only ----
keep <- cd$sample_type == "Primary Tumor"

# ---- 3. Remove non-ductal histologies (neuroendocrine, acinar, cystic, etc.) ----
non_pdac <- "neuroendocrine|islet|acinar|cystadenocarcinoma|mucinous|solid pseudopapillary|carcinoid"
keep <- keep & !grepl(non_pdac, cd$primary_diagnosis, ignore.case = TRUE)
se <- se[, keep]
cd <- as.data.frame(colData(se))

# ---- 4. One sample per patient ----
patient <- substr(colnames(se), 1, 12)
first <- !duplicated(patient)
se <- se[, first]
cd <- as.data.frame(colData(se))
cd$patient <- patient[first]

# ---- 5. Survival outcome (time in months, status 1 = died) ----
cd$status <- ifelse(cd$vital_status == "Dead", 1, 0)
days <- ifelse(cd$status == 1, cd$days_to_death, cd$days_to_last_follow_up)
cd$time_months <- days / 30.44

ok <- !is.na(cd$time_months) & cd$time_months > 0
cat("\nRemoved for missing/zero survival time:", sum(!ok), "\n")
se <- se[, ok]
cd <- cd[ok, ]

# ---- 6. Expression: protein-coding genes, gene symbols, filter low expression ----
counts <- assay(se, "unstranded")
rd <- as.data.frame(rowData(se))
if ("gene_type" %in% names(rd)) {
  pc <- rd$gene_type == "protein_coding"
  counts <- counts[pc, ]
  rd <- rd[pc, ]
}
rownames(counts) <- make.unique(rd$gene_name)
counts <- counts[order(rowMeans(counts), decreasing = TRUE), ]
counts <- counts[!duplicated(sub("\\.[0-9]+$", "", rownames(counts))), ]

y <- DGEList(counts = counts)
y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
logcpm <- cpm(y, log = TRUE, prior.count = 1)

# ---- 7. Clinical table (only columns that exist) ----
want <- c("patient", "time_months", "status", "age_at_diagnosis", "gender",
          "ajcc_pathologic_stage", "primary_diagnosis", "tumor_grade")
clin <- cd[, intersect(want, names(cd))]
rownames(clin) <- colnames(logcpm)

# ---- 8. Save ----
saveRDS(list(counts = y$counts, logcpm = logcpm, clin = clin), "data/tcga_paad_clean.rds")

# ---- 9. Summary for the Methods and logbook ----
cat("\n===== FINAL COHORT =====\n")
cat("Patients (n):", nrow(clin), "\n")
cat("Deaths:", sum(clin$status), " Censored:", sum(clin$status == 0), "\n")
cat("Median follow-up (months):", round(median(clin$time_months), 1), "\n")
cat("Genes kept after filtering:", nrow(logcpm), "\n")
cat("\nDiagnoses (after filtering):\n"); print(table(clin$primary_diagnosis))
if ("ajcc_pathologic_stage" %in% names(clin)) { cat("\nStage:\n"); print(table(clin$ajcc_pathologic_stage)) }
