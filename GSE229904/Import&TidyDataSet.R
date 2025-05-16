# Set working directory to the folder where your dataset is located
# Replace "your/dataset/folder/path" with your actual path or use setwd() interactively
# setwd("your/dataset/folder/path")


# packages     !!!!!!
library(tidyverse)


# Import gene expression data (TPM normalized counts)
genes <- read.delim("GSE229904/Data/GSE229904_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")

# Import phenotypic/clinical metadata
# The first column is set as row names
feno <- read.csv("GSE229904/Data/GSE229904_Annotations.csv", row.names = 1)


# Remove columns that start with the word "character"
feno <- feno %>%
  dplyr::select(!starts_with("character"))

# Clean column names: remove the ".ch1" suffix at the end of column names
colnames(feno) <- gsub("\\.ch1$", "", colnames(feno))

# Convert relevant columns to factors
# You can customize this list depending on your dataset
cols_to_factor <- c("gleason.primary", "gleason.secondary", "gleason.sum", "risk.group",  
                    "source_name", "description", "library_selection", "library_strategy",  
                    "biochemical.recurrence", "cm", "isup", "lapc", "pn", "pt",  
                    "residual_tumor", "tmprss2.erg")

# Apply factor conversion
valid_cols <- intersect(cols_to_factor, colnames(feno))
feno[valid_cols] <- lapply(feno[valid_cols], factor)


# -----------------------------
# Check data types and redundancy
# -----------------------------

# Check if certain columns are numeric (should return a numeric vector preview)
feno$age
feno$gleason.primary
feno$gleason.secondary
feno$gleason.sum
feno$isup
feno$ln.percent
feno$postop.psa
feno$preop.psa

# Check for possible redundancy or invariant columns
# These tables allow you to see if the columns contain repeated or identical values
table(feno$description)
table(feno$library_selection)
table(feno$library_strategy)
table(feno$cm)
table(feno$lapc)
table(feno$residual_tumor)
table(feno$risk.group)

