# How to reproduce

Environment: GitHub Codespace (Ubuntu 24.04), R with the packages listed in results/session_info.txt, and a Python virtual environment (.venv) with the packages in results/python_packages.txt (PyTorch CPU, used only for the neural network in code/09b_nn_cv.py).

If R package installs fail because of conda libraries, run first: export PATH=$(echo "$PATH" | tr ':' '\n' | grep -v '^/opt/conda' | paste -sd:); unset PKG_CONFIG_PATH LD_LIBRARY_PATH

Data are public and de-identified: TCGA-PAAD from the GDC through TCGAbiolinks, and the GEO series GSE62452, GSE21501, GSE57495, GSE71729 and GSE79668 (the last through its supplementary counts file). Files in data/ are not committed and are rebuilt by the scripts. Checksums of the processed inputs used for the final results are in results/data_checksums.txt.

Scripts in code/ are numbered in the order they were written and run. Main steps: 01 and 04 (TCGA download and filtering), 06 and 13 and 23 (external cohorts), 07 and 09a and 09b (model comparison), 10a and 10b (lock and first external test), 24 (duplicate screen), 25 (final scoring), 30 and 31 (pathways and results summary), 40 to 43 (integrity audit, extended analyses, model card). Other numbers are exploratory analyses.

Decisions and dates are in docs/logbook.md and docs/plan_*.md. docs/INTEGRITY_REPORT.md shows the commit order and which parts were specified before or after the results were known. Random seeds are fixed in each script.
