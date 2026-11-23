# 04_filter_tcga.R: keep ductal adenocarcinoma only
d <- readRDS("data/tcga_paad_clean.rds")
bad <- "basal cell|melanoma|undifferentiated|mixed subtypes"
keep <- !grepl(bad, d$clin$primary_diagnosis, ignore.case = TRUE)
cat("Removed", sum(!keep), "non-PDAC / ambiguous samples\n")
d$clin <- d$clin[keep, ]
d$logcpm <- d$logcpm[, keep]
d$counts <- d$counts[, keep]
saveRDS(d, "data/tcga_paad_pdac.rds")
cat("\n===== TCGA-PAAD PDAC COHORT =====\n")
cat("Patients (n):", nrow(d$clin), "\n")
cat("Deaths:", sum(d$clin$status), " Censored:", sum(d$clin$status == 0), "\n")
cat("Median observed time (months):", round(median(d$clin$time_months), 1), "\n")
print(table(d$clin$primary_diagnosis))
