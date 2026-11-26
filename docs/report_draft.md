# Predicting Pancreatic Cancer Survival and Its Biological Correlates from Gene Expression Using Machine Learning and a Prognostic Gene Signature

**Research report draft · November 25, 2026**

## Abstract

Pancreatic ductal adenocarcinoma (PDAC) has poor survival, motivating molecular
prognostic tools that can be tested across patient cohorts. We developed an
eight-gene expression score using TCGA-PAAD and assessed it without refitting
in public pancreatic cancer cohorts. In repeated five-fold cross-validation
(10 repeats), Cox-LASSO, random survival forest, and a small neural network had
mean C-indices of 0.615, 0.612, and 0.608, respectively; the neural network
did not improve prediction. A pooled TCGA/Moffitt training attempt also failed
to improve external performance, so the TCGA-only signature remained primary.
The locked score had a C-index of 0.565 (95% bootstrap CI 0.472–0.654) in
GSE21501 and 0.558 (0.449–0.664) in duplicate-screened GSE57495. Neither
fresh cohort met the pre-specified replication rule: the C-index interval had
to exclude 0.5 and the hazard ratio per standard deviation (SD) had to exceed
1. The previously viewed GSE62452 reference cohort met that rule (C-index
0.592, 0.506–0.678; HR/SD 1.40, 1.07–1.82), but is not fresh validation.
These results support further independent testing, not clinical use.

## Introduction

Pancreatic cancer has a five-year relative survival rate of only 13% in the
United States. Only about 17% of patients are diagnosed while the tumor is
still local, and five-year survival falls from 44% for local disease to 3%
for distant disease [1]. Diagnosis is often late because symptoms are vague,
and survival varies substantially among patients. Better prognostic tools may
help clinicians identify changing levels of risk and assign treatment
accordingly.

Gene-expression signatures are compact sets of genes whose activity correlates
with cancer outcomes. Prior pancreatic cancer studies have developed
Cox-based, immune-related, exosomal, and metastasis-related signatures, with
varying external validation [2–4]. Many approaches rely on a single modeling
method and do not connect model performance to biological correlates. This
project compares three survival-learning methods and evaluates whether a
locked signature retains prognostic value outside its training cohort.

The study had two questions: whether a compact expression signature can
predict overall survival, and whether a signature built in TCGA-PAAD remains
associated with survival and biological programs in independent pancreatic
cancer cohorts. The signature was fitted only in TCGA, locked before external
scoring, and not retuned on test cohorts.

The primary hypothesis was that higher locked-signature risk would be
associated with shorter overall survival in independent cohorts. For a cohort
to count as a replication, its 95% bootstrap confidence interval (CI) for the
C-index had to exclude 0.5, and the score's hazard ratio (HR) per SD had to be
greater than 1. The signature, duplicate screen, eligibility criteria, and
replication rule were committed before the new-cohort score was calculated.

## Methods

### Public datasets and cohort eligibility

TCGA-PAAD RNA-sequencing counts were obtained from the Genomic Data Commons.
One primary tumor per patient was retained, non-ductal or ambiguous histologies
were excluded, and patients without a positive survival time were removed.
The resulting PDAC training set contained 145 patients. GSE62452 was used as
an external reference cohort; GSE21501 and GSE57495 were the fresh test
cohorts. GSE71729 (Moffitt) was considered for a pooled-training experiment.
The TCGA, GSE62452, GSE21501, GSE71729, GSE57495, and GSE79668 records are
listed in References [1–6].

The fresh-cohort plan required primary tumors, overall survival, at least 50
patients and 30 events after duplicate screening, and measurement of all eight
signature genes. Candidate datasets were screened before score calculation.
GSE78229 was excluded because its sample IDs and clinical survival
observations matched GSE62452. GSE28735 did not meet the 50-patient minimum.
GSE79668 had 51 tumors and 45 events before screening; two samples were
flagged, leaving 49 patients, and SLC16A2 was missing, so it failed the rule
and was not scored. E-MTAB-6134 and ICGC PACA-CA/PACA-AU were not added:
the ArrayExpress download route and ICGC download route could not be verified
before scoring and are not silently treated as eligible cohorts.

### Expression processing and model development

For TCGA, protein-coding counts were filtered with edgeR's `filterByExpr`,
normalized using TMM, and represented as log-CPM. For model comparison, genes
were restricted to those shared with GSE62452; genes with average TCGA
expression below the median were excluded, and the 2,000 genes with greatest
median absolute deviation (MAD) were retained. Gene values were standardized
within each dataset. The same median-expression filter and 2,000-gene MAD
feature limit were used in the pooled-training attempt.

Cox-LASSO and random survival forest were compared using 10 repeats of
status-stratified five-fold cross-validation in TCGA. Within each training
fold, Cox-LASSO used univariate Cox screening (p < 0.01; the top 20 genes were
used if fewer passed) and a penalty chosen to retain approximately 5–10
features. The forest ranked genes by permutation importance and used the top
10. A separate neural network used a 16-unit hidden layer, dropout, Cox loss,
and early stopping on a validation subset of each training fold. CV C-indices
are reported as means and standard deviations across folds, not as independent
patient-level confidence intervals.

The final eight-gene Cox-LASSO model was fit using TCGA only and locked before
external scoring. `glmnet`'s Cox implementation used its default Breslow tie
handling. The locked risk score was computed by standardizing each of the
eight genes within a cohort and taking the sum of expression z-scores
multiplied by the fixed TCGA coefficients. Neither gene selection nor
coefficients were refit in validation cohorts.

| Gene | Locked coefficient |
|---|---:|
| ARNTL2 | 0.1788 |
| CYP27A1 | -0.1055 |
| MET | 0.0499 |
| SIPA1L2 | -0.0263 |
| DPYD | 0.0192 |
| SLC16A2 | -0.0161 |
| MICAL1 | -0.0148 |
| CAV2 | 0.0117 |

### Duplicate screening and validation

Before any fresh-cohort scoring, each new sample was compared with TCGA,
GSE62452, GSE21501, GSE71729, and the other candidate cohort. For each
pairwise comparison, the 3,000 shared genes with the strongest mean MAD rank
were standardized within cohort and correlated. A sample was flagged if its
best expression correlation was at least 0.8, or if correlation was at least
0.5 and the matched survival time differed by no more than 0.5 months with
identical event status. Flagged patients were removed; a cohort falling below
the size/event criteria or missing a signature gene was excluded without
rescue.

Discrimination was summarized by Harrell's C-index with 2,000 patient-level
bootstrap resamples for percentile 95% CIs. The HR per SD was estimated using
a univariable Cox model with the cohort-standardized risk score. The
pre-specified replication definition was a C-index CI wholly above or below
0.5 **and** an HR point estimate greater than 1; it did not include a separate
p-value cutoff. Fixed-effect HR pooling was an exploratory, post hoc analysis.

### Secondary and exploratory analyses

For expression-level comparison, TCGA counts were analyzed with voom/limma
and microarray cohorts with limma. The correlation of gene-level t-statistics
between TCGA and GSE62452 was calculated. Hallmark pathway activity was tested
with two complementary approaches. The primary gene-set method used `camera()`
with estimated inter-gene correlation rather than the lower default assumption
used by an initial `cameraPR` run. The complementary method calculated a mean
z-score for each Hallmark set in each patient and correlated that score with
the continuous locked-risk score using Spearman correlation, with
Benjamini–Hochberg correction within cohort. A 32-set adjusted score analysis
that additionally modeled MCP-counter stromal and immune proxies was retained
only as a sensitivity analysis because of likely over-adjustment.

Immune estimates used MCP-counter only. Checkpoint genes, a cytolytic score,
Moffitt subtype, TCGA mutations, and immune/stromal adjustments were
exploratory. These analyses did not alter the locked signature. No single-cell
RNA-seq or Human Protein Atlas claims are included; these analyses were not
performed and are reserved for future work.

## Results

### Model comparison and the pooled-training attempt

Across 50 held-out folds, mean C-indices were 0.615 (SD 0.083) for Cox-LASSO,
0.612 (0.073) for random survival forest, and 0.608 (0.079) for the neural
network. Cox-LASSO exceeded the forest in 26 of 50 folds; the mean LASSO–
forest difference was 0.003. The neural network's mean was 0.007 below
Cox-LASSO, so it was not adopted.

An overlap screen found that 86 of 123 GSE71729 patients shared identical
survival time and event status with GSE21501, despite weak expression
correlations. Only 37 patients and 26 events remained eligible for the pooled
training attempt. That model used TCGA plus the non-overlapping Moffitt
patients (182 total, 113 events) and a seven-gene signature. Its
cross-validated C-index was 0.588 in TCGA and 0.502 in Moffitt. On GSE21501,
the pooled score's C-index was 0.551 (95% CI 0.474–0.631), compared with
0.565 (0.480–0.655) for the TCGA-only score; on GSE62452 the respective
estimates were 0.548 and 0.592. The pooled approach did not improve
generalization and was not selected as primary.

### Locked-signature results

The TCGA-only score had an apparent training C-index of 0.706; this is
optimistic because it is measured in the data used to fit the signature.
External-cohort results were:

| Cohort | Role | Patients / events | C-index (95% bootstrap CI) | HR per SD (95% Cox CI) | p | Meets rule? |
|---|---|---:|---:|---:|---:|---|
| GSE62452 | Previously viewed reference | 65 / 49 | 0.592 (0.506–0.678) | 1.40 (1.07–1.82) | 0.013 | Yes; not fresh |
| GSE21501 | Fresh | 102 / 66 | 0.565 (0.472–0.654) | 1.22 (0.94–1.58) | 0.137 | No |
| GSE57495, screened | Fresh | 59 / 38 | 0.558 (0.449–0.664) | 1.31 (0.90–1.90) | 0.156 | No |
| GSE57495, all samples | Sensitivity only | 63 / 42 | 0.571 (0.465–0.671) | 1.33 (0.93–1.89) | 0.116 | No |

Neither fresh cohort met the pre-specified replication rule: both C-index
intervals included 0.5. Both HR point estimates were above 1, but their
confidence intervals included 1. The GSE62452 reference passed the operational
rule, but it had been viewed before the fresh-cohort plan and is not
independent confirmation.

As a post hoc summary, fixed-effect pooling of HRs across GSE62452, GSE21501,
and screened GSE57495 gave HR/SD 1.30 (95% CI 1.10–1.54; p = 0.0017;
heterogeneity p = 0.77). Pooling only the two fresh cohorts gave 1.25
(1.01–1.54; p = 0.042; heterogeneity p = 0.75). These pooled estimates are
not a substitute for the pre-specified cohort-by-cohort replication test.

### Expression, immune, subtype, and mutation findings

Gene-level t-statistics were correlated between TCGA and GSE62452
(Spearman rho = 0.768). Of 5,925 genes with BH-adjusted p < 0.05 in TCGA,
3,041 also had nominal p < 0.05 and the same direction in GSE62452.

The initial `cameraPR` run reported concordant Hallmark sets, but its assumed
inter-gene correlation (0.01) was lower than the estimated median correlation
(0.064 in TCGA and 0.081 in GSE62452). With `camera()` and estimated
correlation, no Hallmark set met the BH < 0.05, same-direction rule in TCGA
and GSE62452; none replicated in GSE21501 either. The complementary, unadjusted
Hallmark score analysis found 17 same-direction associations after cohort-wise
BH correction, led by mitotic spindle,
glycolysis, G2M checkpoint, E2F targets, mTORC1 signaling, and MYC targets V1.
These score-test results are exploratory and do not replace the more
conservative gene-set test. TGF-beta signaling and hypoxia did not meet the
unadjusted score-test replication criteria. The model therefore does **not**
claim replicated pathway enrichment from the primary `camera()` analysis.
The adjusted score-test sensitivity result is not reported as replication.

In TCGA, nine of ten MCP-counter populations had lower estimated abundance in
the high-risk group after BH correction; only endothelial cells also met the
correction threshold in GSE62452. Checkpoint-gene and cytolytic-score findings
did not replicate. Adjustment for immune and stromal proxies did not remove
the risk association (TCGA HR/SD 1.98 unadjusted vs 2.14 adjusted; GSE62452
1.40 vs 1.40). These are bulk-expression estimates, not direct cell counts
or evidence of immune function.

In the Moffitt cohort, code 2 was labeled basal-like using marker expression.
Across all 123 patients, basal-like tumors had higher risk scores (AUC 0.70)
and worse survival (HR 1.85, p = 0.010). Because 86 Moffitt samples overlapped
GSE21501, this is not independent validation; after removing those samples,
the subtype-risk AUC was 0.61 (p = 0.284) in 37 patients.

Among 134 TCGA patients with mutation data and risk scores, KRAS mutations
occurred in 82.4% of high-risk versus 48.5% of low-risk tumors; TP53
mutations occurred in 75.0% versus 50.0%. Associations persisted after
adjustment for the immune proxy (KRAS OR 2.40, p = 0.0010; TP53 OR 2.04,
p = 0.0030). Mutation counts were also higher in the high-risk group, but the
count is not a calibrated mutations-per-megabase tumor mutational burden
measure.

## Discussion

The main result is deliberately modest: a TCGA-derived eight-gene score had
some prognostic signal in the previously viewed GSE62452 cohort, but it did
not meet the pre-specified replication rule in either fresh cohort. The
positive HR point estimates in GSE21501 and GSE57495 are directionally
consistent with the hypothesis, but their C-index intervals overlap chance
and their HR intervals include 1. The evidence is therefore insufficient to
claim reliable transportability.

The validation workflow prioritized protection against duplicate patients
and post hoc model tuning. It screened sample-level expression together with
survival observations, excluded ineligible cohorts before scoring, kept the
signature fixed, and reported non-replication rather than retraining on test
data. The pooled-training attempt and neural-network comparison also show why
more model complexity or more training records should not be assumed to
improve external performance.

Exploratory analyses suggest that risk is associated with both tumor
expression programs and tissue composition, and that high-risk tumors more
often carry KRAS and TP53 mutations. These observations are hypothesis
generating. The pathway results depended on the gene-set method, immune
deconvolution is indirect, and the Moffitt subtype analysis is confounded by
patient overlap. None establishes a causal mechanism.

## Limitations

The fresh validation cohorts were small, with 38 and 66 events, so estimates
are imprecise. Cohorts differ in platform, processing, clinical annotation,
and case selection. Within-cohort standardization may reduce scale
differences but does not remove all batch or population effects. Survival
definitions and missing covariates may differ across studies. Only one
previously viewed reference and two fresh cohorts met or nearly met the
evaluation requirements. A C-index does not establish calibration, clinical
utility, or benefit from treatment decisions. The TCGA apparent C-index is
optimistic. Several secondary analyses were exploratory and involved multiple
comparisons; their findings require independent confirmation. The score is a
research result, not a diagnostic or clinical decision tool.

## Conclusion

The locked TCGA-only eight-gene score did not meet the pre-specified
replication definition in two fresh PDAC cohorts. The data support further
testing with larger, genuinely independent cohorts and pre-specified
calibration and clinical-utility analyses; they do not support clinical use
or a claim of validated prognostic performance.

## Figures

1. **Cross-validation model comparison:** [results/26_cv_comparison.png](../results/26_cv_comparison.png). Points show mean C-index and whiskers show fold-to-fold SD, not confidence intervals.
2. **Locked-signature external evaluation:** [results/26_validation_forest.png](../results/26_validation_forest.png). C-index intervals are bootstrap percentile CIs; HR intervals are Cox-model CIs. GSE62452 is labeled as previously viewed.
3. **Kaplan–Meier display for the earlier external reference:** [results/10_km_external.png](../results/10_km_external.png). This is secondary to the continuous-score metrics and is not fresh validation.

## References

1. Gene Expression Omnibus. **GSE62452:** Microarray gene-expression profiles of 69 pancreatic tumors and 61 adjacent non-tumor tissue from patients with pancreatic ductal adenocarcinoma. [NCBI GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62452).
2. Gene Expression Omnibus. **GSE21501:** A six-gene signature predicts survival of patients with localized pancreatic ductal adenocarcinoma. [NCBI GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21501).
3. Gene Expression Omnibus. **GSE71729:** Virtual microdissection of pancreatic ductal adenocarcinoma reveals tumor and stroma subtypes. [NCBI GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71729).
4. Gene Expression Omnibus. **GSE57495:** Development of a prognostic gene signature in pancreatic cancer. [NCBI GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57495).
5. Cancer Genome Atlas Research Network. Integrated Genomic Characterization of Pancreatic Ductal Adenocarcinoma. *Cancer Cell*. 2017;32(2):185–203.e13. [doi:10.1016/j.ccell.2017.07.007](https://doi.org/10.1016/j.ccell.2017.07.007).
6. Gene Expression Omnibus. **GSE79668:** RNA-sequencing of human pancreatic adenocarcinoma cancer tissues (screened and excluded before scoring). [NCBI GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79668).
7. Cox DR. Regression models and life-tables. *J R Stat Soc Series B*. 1972;34(2):187–220. [doi:10.1111/j.2517-6161.1972.tb00899.x](https://doi.org/10.1111/j.2517-6161.1972.tb00899.x).
8. Simon N, Friedman J, Hastie T, Tibshirani R. Regularization paths for Cox's proportional hazards model via coordinate descent. *J Stat Softw*. 2011;39(5):1–13. [doi:10.18637/jss.v039.i05](https://doi.org/10.18637/jss.v039.i05).
9. Becht E, et al. Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. *Genome Biol*. 2016;17:218. [doi:10.1186/s13059-016-1070-5](https://doi.org/10.1186/s13059-016-1070-5).
10. Liberzon A, et al. The Molecular Signatures Database Hallmark Gene Set Collection. *Cell Syst*. 2015;1(6):417–425. [doi:10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004).

## Acknowledgments and data-use note

Generative AI (Copilot SDK in VS Code) assisted with code edits, organizing
analysis outputs, and drafting this report. The student researcher should
verify every interpretation, citation, and final submission requirement and
remains responsible for the submitted work. No new participants were
recruited or contacted; the analyses used previously collected public human
research data. The fair's SRC/IRB must determine the applicable paperwork
before submission.
