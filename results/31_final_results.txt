# Final results (generated from results/*.csv, do not edit by hand)

## Locked signature (TCGA only, 8 genes, positive coef = higher expression, higher risk)

| gene | coef |
|---|---|
| ARNTL2 | 0.1788 |
| CYP27A1 | -0.1055 |
| MET | 0.0499 |
| SIPA1L2 | -0.0263 |
| DPYD | 0.0192 |
| SLC16A2 | -0.0161 |
| MICAL1 | -0.0148 |
| CAV2 | 0.0117 |

## TCGA cross-validated C-index (screening inside each fold)

| run | model | mean_cindex | sd |
|---|---|---|---|
| first run (variance filter only) | Cox-LASSO | 0.537 | 0.072 |
| first run (variance filter only) | RandomSurvivalForest | 0.558 | 0.068 |
| final run (median expression filter added) | Cox-LASSO | 0.615 | 0.083 |
| final run (median expression filter added) | RandomSurvivalForest | 0.612 | 0.073 |
| final run (median expression filter added) | NeuralNetwork | 0.608 | 0.079 |

## External validation of the locked signature (each cohort scored once)

| cohort | n | events | cindex_95ci | hr_per_sd_95ci | p | meets_rule |
|---|---|---|---|---|---|---|
| GSE62452 (viewed before) | 65 | 49 | 0.592 (0.506-0.678) | 1.40 (1.07-1.82) | 0.013 | TRUE |
| GSE21501 (fresh) | 102 | 66 | 0.565 (0.472-0.654) | 1.22 (0.94-1.58) | 0.137 | FALSE |
| GSE57495 (fresh, after screen) | 59 | 38 | 0.558 (0.449-0.664) | 1.31 (0.90-1.90) | 0.156 | FALSE |
| GSE57495 (all 63, sensitivity) | 63 | 42 | 0.571 (0.465-0.671) | 1.33 (0.93-1.89) | 0.116 | FALSE |

## Pooled hazard ratio (post hoc, fixed effect)

- All three external cohorts: HR per SD 1.30 (1.10-1.54), p = 0.0017, heterogeneity p = 0.77, patients 226, events 153
- Fresh cohorts only: HR per SD 1.25 (1.01-1.54), p = 0.0424, heterogeneity p = 0.75, patients 161, events 104

## Hallmark pathways associated with the risk score (patient level)

Replicated pathways (TCGA FDR < 0.05, external meta FDR < 0.05, same sign in all four cohorts): 23 of 50.

| pathway | direction | rho_tcga | rho_external_meta | fdr_tcga | fdr_meta |
|---|---|---|---|---|---|
| G2M_CHECKPOINT | higher in high risk | 0.62 | 0.61 | 8.8e-16 | 7.0e-24 |
| E2F_TARGETS | higher in high risk | 0.59 | 0.60 | 5.6e-14 | 2.4e-23 |
| MTORC1_SIGNALING | higher in high risk | 0.56 | 0.55 | 3.8e-12 | 1.8e-18 |
| GLYCOLYSIS | higher in high risk | 0.63 | 0.55 | 8.8e-16 | 2.1e-18 |
| MYC_TARGETS_V1 | higher in high risk | 0.50 | 0.47 | 1.1e-09 | 3.2e-13 |
| MITOTIC_SPINDLE | higher in high risk | 0.64 | 0.47 | 3.4e-16 | 7.8e-13 |
| MYOGENESIS | lower in high risk | -0.30 | -0.43 | 6.8e-04 | 1.5e-10 |
| PANCREAS_BETA_CELLS | lower in high risk | -0.35 | -0.40 | 6.3e-05 | 1.9e-09 |
| DNA_REPAIR | higher in high risk | 0.36 | 0.38 | 4.1e-05 | 3.0e-08 |
| PROTEIN_SECRETION | higher in high risk | 0.40 | 0.34 | 3.0e-06 | 9.3e-07 |
| UNFOLDED_PROTEIN_RESPONSE | higher in high risk | 0.41 | 0.34 | 2.1e-06 | 9.3e-07 |
| ANDROGEN_RESPONSE | higher in high risk | 0.27 | 0.34 | 2.6e-03 | 9.3e-07 |
| ESTROGEN_RESPONSE_LATE | higher in high risk | 0.27 | 0.32 | 2.6e-03 | 3.4e-06 |
| MYC_TARGETS_V2 | higher in high risk | 0.39 | 0.31 | 7.3e-06 | 7.6e-06 |
| KRAS_SIGNALING_DN | lower in high risk | -0.38 | -0.30 | 1.3e-05 | 1.1e-05 |
| HYPOXIA | higher in high risk | 0.32 | 0.28 | 2.8e-04 | 8.8e-05 |
| TGF_BETA_SIGNALING | higher in high risk | 0.44 | 0.27 | 2.0e-07 | 1.2e-04 |
| PI3K_AKT_MTOR_SIGNALING | higher in high risk | 0.37 | 0.27 | 2.5e-05 | 1.4e-04 |
| APOPTOSIS | higher in high risk | 0.23 | 0.26 | 1.3e-02 | 2.1e-04 |
| ESTROGEN_RESPONSE_EARLY | higher in high risk | 0.28 | 0.24 | 1.6e-03 | 5.8e-04 |
| BILE_ACID_METABOLISM | lower in high risk | -0.28 | -0.23 | 1.7e-03 | 8.9e-04 |
| NOTCH_SIGNALING | higher in high risk | 0.22 | 0.22 | 1.4e-02 | 1.7e-03 |
| P53_PATHWAY | higher in high risk | 0.22 | 0.22 | 1.6e-02 | 1.9e-03 |

## What these results support

- The locked signature points the same way in every external cohort. Only GSE62452 met the pre-registered rule individually.
- The effect is modest (C-index about 0.56 to 0.59). It is not a clinically usable predictor.
- Pathway results are associations at the patient level, replicated across cohorts. They are not causal and tumor purity is an unresolved confounder.

## What they do not support

- Any C-index near 0.7, any claim that pooled training or the neural network helped, any claim that the signature beats stage, any use of the earlier cameraPR p-values as evidence, any causal wording.
