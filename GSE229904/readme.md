# GSE229904 Analysis Pipeline

This repository contains the R scripts used to clean, process, and analyze the [GSE229904 dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229904), which focuses on prostate cancer patient data. [Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10253008/)

## Overview

The analysis is divided into two main scripts:

### 1. Data Cleaning and Preparation
- Standardizes and tidies both phenotypic (`feno`) and genotypic (`genes`) datasets.
- Extracts patient IDs and reorganizes key columns.
- Filters for RNA-Seq data from primary tumor samples.
- Saves intermediate tables as `.rds` files for later use.

### 2. Exploratory Analysis and Modeling
- **Descriptive Statistics**: Calculates means, medians, and summaries of key clinical variables.
- **Associations**: Explores relationships between Gleason scores, PSA levels, age, and other variables using jitter plots and boxplots.
- **GT Summary Tables**: Generates visual summaries (`gtsummary`) to understand patient characteristics and stratify by recurrence and risk group.
- **Patient Consistency Check**: Validates that the number of unique patients matches the original publication.
- **Logistic Regression**: Assesses how clinical and gene expression variables relate to biochemical recurrence.
- **Gene Expression Analysis**: Performs differential gene expression (DGE) analysis between recurrence groups using the Wilcoxon test.
- **Visualization**: Includes plots like histograms, volcano plots, regression coefficient bar plots, and boxplots to enhance interpretability.

## Key Outcomes
- Identified significant phenotypic and genotypic factors associated with biochemical recurrence.
- Visualized the differential expression of genes and their statistical significance.
- Ensured data consistency and reproducibility through saved RDS outputs and versioned plots.

## Files
- `Script_1_Cleaning.R`: Data cleaning and preparation
- `Script_2_Analysis.R`: Full statistical analysis and modeling
- `Plots/`: Folder containing all generated figures
- `Tablas.rds/`: Processed tables in `.rds` format
- `Tablas.gt/`: GT summary tables exported as images

---

Feel free to adapt this template if you split the analysis into more scripts later. Would you like me to auto-generate this as a file in Markdown format?
