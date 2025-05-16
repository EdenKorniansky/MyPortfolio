#You have your datasets feno and genes imported and tidied
#Now we focus on the analisis of this data

############################################
#Analysis of feno's stats
# Calculate the mean age and round the result
mean_age <- mean(feno$age, na.rm = TRUE)
round(mean_age)

# Display a summary of the age data (including min, median, mean, max, etc.)
summary(feno$age)

# Calculate the median age
median_age <- median(feno$age, na.rm = TRUE)
median_age

# Calculate the mean postoperative PSA and round the result
mean_psa_postop <- mean(feno$postop.psa, na.rm = TRUE)
round(mean_psa_postop)

# Calculate the mean preoperative PSA and round the result
mean_psa_preop <- mean(feno$preop.psa, na.rm = TRUE)
round(mean_psa_preop)



###########################################################################
# Specific questions to understand some characteristics of our patients


#Quantity of primary tumor tissue and normal prostatic tissue?
tipo_tejido <- table(feno$source_name_ch1)
tipo_tejido
# There is more primary tumor tissue than normal prostatic tissue

# How many patients of each age have each type of primary/secondary Gleason score?
table(feno$gleason.primary, feno$age)
table(feno$gleason.secondary, feno$age)
table(feno$gleason.sum, feno$age)
table(feno$biochemical.recurrence)
table(feno$tmprss2.erg)

# Age histogram
library(ggplot2)
feno$age

ggplot(feno, aes(x = age)) +
  geom_histogram(fill = "blue", color = "black") +
  labs(title = "Age Histogram", x = "Age", y = "Frequency") +
  theme_minimal()

ggsave("Plots/HistEdades.png", bg = "white", height = 4, width = 6, units = "in")

# Is there an association between the primary Gleason pattern and preoperative PSA? (jitter plot)
ggplot(feno, aes(x = gleason.primary, y = preop.psa)) +
  geom_jitter(width = 0.1)

# Is there an association between the secondary Gleason pattern and postoperative PSA? (jitter plot)
ggplot(feno, aes(x = gleason.secondary, y = postop.psa)) +
  geom_jitter(width = 0.1)

# Is there an association between the secondary Gleason pattern and preoperative PSA? (jitter plot)
ggplot(feno, aes(x = gleason.secondary, y = preop.psa)) +
  geom_jitter(width = 0.1)

# Is there an association between the primary Gleason pattern and postoperative PSA? (jitter plot)
ggplot(feno, aes(x = gleason.primary, y = postop.psa)) +
  geom_jitter(width = 0.1)




##########################################################################
library(skimr)
library(gtsummary)
library(tidyverse)
library(gt)
skim(feno)

# For getting a better summary of patient characteristics:
# General summary provided by R using table feno
tab1 <- feno %>% select(-title, -geo_accession, -gleason.primary, -gleason.secondary) %>% 
  tbl_summary()

gtsave(data = as_gt(tab1), filename = "Tablas.gt/prueba.gt.feno.png")
tab1




# Since I found irrelevant columns in my gt table, I decided to keep the relevant ones in a new list called lista_variables:

# See all column names in feno
colnames(feno)

# Create a list of all variables that contain relevant information from feno
variables <- colnames(feno)[c(3:9, 12:22)]
print(variables)

# Store these variables in a list
lista_variables <- as.list(variables)
print(lista_variables)




# Save phenotypic and genotypic tables as RDS files for later use

# Replace the path below with the folder where you want to save your data
saveRDS(feno, "Tablas.rds/PhenotypicTable.rds")
saveRDS(genes, "Tablas.rds/GenotypicTable.rds")






# Create a new table called feno1 from feno, keeping only the rows corresponding to RNA-Seq and primary tumor samples. So, pay attention in columns "library_strategy" and "source_name_ch1"). From now on, use this filtered data using these new tables.

feno1 <- feno %>%
  dplyr::filter(library_strategy == "RNA-Seq" & source_name_ch1 == "primary tumor")

# Adapt the genes table to match feno1; the result is called genes1
genes1 <- genes %>%
  dplyr::select('GeneID', any_of(feno1$geo_accession))



#Since multiple samples were collected from each patient, it was useful to create a new column called $ID.paciente, which extracts the patient ID from the first three digits of the feno$title column. This allowed us to confirm that the number of unique patients in the dataset matches the number reported in the original publication.


# Create a new column with patient ID, extracting the first 1 to 3 digits (optionally preceded by a letter)
feno$ID.paciente <- str_extract(feno$title, "^[A-Za-z]?[0-9]{1,3}")

# Create a second column with the remaining part of the title, representing the sample or analysis type
feno$ID.paciente.resto <- str_remove(feno$title, "^[A-Za-z]?[0-9]{1,3}")

# Reorder columns for easier visualization and sort by patient ID and sample type
feno <- feno %>%
  select(title, ID.paciente, ID.paciente.resto, source_name_ch1, library_strategy, everything()) %>%
  arrange(ID.paciente, source_name_ch1)



# Create a summary table of feno1 using gtsummary, excluding irrelevant columns
tab2 <- feno1 %>%
  select(-title, -geo_accession, -gleason.primary, -gleason.secondary, -ID.paciente, -ID.paciente.resto,
         -source_name_ch1, -library_strategy, -description, -library_selection) %>%
  tbl_summary()

# Save the summary table as a .png file
gtsave(data = as_gt(tab2), filename = "Tablas.gt/prueba.gt.feno1.png")


# Same gt summary table as before, but now stratified by a specific variable: biochemical.recurrence
tab3 <- feno1 %>%
  select(-title, -geo_accession, -gleason.primary, -gleason.secondary, -ID.paciente, -ID.paciente.resto,
         -source_name_ch1, -library_strategy, -description, -library_selection) %>%
  tbl_summary(by = biochemical.recurrence)

# Save the stratified summary table as a .png file
gtsave(data = as_gt(tab3), filename = "Tablas.gt/prueba.gt.biochrecu.png")
tab3


# Same gt summary table as before, but now stratified by a specific variable: risk.group
tab4=feno1 %>% select(-title,-geo_accession,-gleason.primary,-gleason.secondary, -ID.paciente, -ID.paciente.resto,-source_name_ch1, -library_strategy, -description, -library_selection) %>% 
  tbl_summary(by=risk.group)

# Save the stratified summary table as a .png file
gtsave(data = as_gt(tab4), filename = "Tablas.gt/prueba.gt.riskgroup.png")
tab4



saveRDS(feno1, file = "Tablas.rds/PhenotypicTable1.rds")
saveRDS(genes1, file = "Tablas.rds/GenotypicTable1.rds")





# Boxplots to determine the age distribution of patients who experience recurrence

# 1. Evaluated using biochemical recurrence status
library(ggpubr)
ggplot(data = feno1 %>% filter(biochemical.recurrence != "NA"), 
       aes(x = biochemical.recurrence, y = age)) +
  geom_boxplot(fill = "violet") +
  labs(title = "Age distribution by biochemical recurrence status",
       x = "Biochemical recurrence", 
       y = "Age (years)") +
  theme_minimal() +
  scale_y_continuous(breaks = seq(floor(min(feno1$age, na.rm = TRUE) / 5) * 5, 
                                  ceiling(max(feno1$age, na.rm = TRUE) / 5) * 5, 
                                  by = 5)) +
  stat_compare_means(method = "t.test", label = "p.format")
ggsave("Plots/boxplot_age_biochrecu.png", width = 8, height = 6, dpi = 300)

# 2. Evaluated based on postoperative PSA levels
# Need to determine a PSA threshold that is clinically relevant for recurrence.
# For now, we use an arbitrary cutoff of 0.30
ggplot(feno1, aes(x = postop.psa >= 0.3, y = age)) +
  geom_boxplot(fill = "violet") +
  scale_x_discrete(labels = c(">0.3", "<=0.3")) +
  labs(title = "Age distribution by postoperative PSA levels",
       x = "Postoperative PSA levels", 
       y = "Age (years)") +
  theme_minimal()
ggsave("Plots/boxplot_edades_postopPSA.png", width = 8, height = 6, dpi = 300)

# 3. Evaluated based on residual tumor status (R0 vs R1)
ggplot(feno1, aes(x = r.residual_tumor, y = age)) +
  geom_boxplot(fill = "violet") +
  labs(title = "Age distribution by residual tumor status", 
       x = "Residual tumor", 
       y = "Age (years)") +
  theme_minimal()
ggsave("Plots/boxplot_edades_residualtumor.png", width = 8, height = 6, dpi = 300)

# Is there an association between age and preoperative PSA? (scatter plot)
library(ggplot2)

ggplot(feno1, aes(x = age, y = preop.psa, color = gleason.sum)) +
  geom_jitter(alpha = 0.6, width = 0.2, height = 0.2) +
  geom_smooth(method = "lm", color = "black", se = TRUE, fill = "gray70") + 
  labs(title = "Relationship between age and preoperative PSA levels",
       x = "Age (years)", 
       y = "Preoperative PSA levels",
       color = "Gleason Score") +
  theme_minimal() +
  scale_color_brewer(palette = "Paired") 

ggsave("Plots/jitterplot_edades_PSApreoperatorio.png", width = 8, height = 6, dpi = 300)





#Logistic regresion

# This step is FUNDAMENTAL because it defines "brf" as improvement and "bcr" as deterioration
table(feno1$biochemical.recurrence)
class(feno1$biochemical.recurrence)
feno1$biochemical.recurrence <- factor(feno1$biochemical.recurrence, levels = c("brf", "bcr"))


# Logistic regression between biochemical recurrence and preoperative PSA
modelo <- glm(biochemical.recurrence ~ preop.psa, data = feno1, family = binomial)
summary(modelo)

# Logistic regression between biochemical recurrence and Gleason score
modelo2 <- glm(biochemical.recurrence ~ gleason.sum, data = feno1, family = binomial)
summary(modelo2)



#To perform a logistic regression between gene expression levels and the variable biochemical.recurrence, we first need to combine the phenotypic and genotypic datasets.
#We start by selecting only the genes of interest: AR, YWHAZ, NDRG1, APOE, CRIP2, and


# Selected genes by their GeneID
muestra <- c(367, 7534, 10397, 348, 1397, 10631)

# Find indices of the genes of interest
indices <- which(genes1$GeneID %in% muestra)

# Subset the gene expression table using those indices
nueva_tabla <- genes1[indices,]
View(nueva_tabla)

# Rename gene IDs to gene names
nueva_tabla[1,1] <- "YWHAZ"
nueva_tabla[2,1] <- "NDRG1"
nueva_tabla[3,1] <- "POSTN"
nueva_tabla[4,1] <- "CRIP2"
nueva_tabla[5,1] <- "APOE"
nueva_tabla[6,1] <- "AR"

# Set gene names as row names and remove GeneID column
row.names(nueva_tabla) <- nueva_tabla$GeneID
nueva_tabla = nueva_tabla[-1]

# Transpose to have patients in rows and genes in columns
genes_translocados <- t(nueva_tabla)
View(genes_translocados)

# Save the transposed table
write.csv(genes_translocados, "your/dataset/folder/path/genes_translocados.csv", row.names = TRUE)
saveRDS(genes_translocados, "Tablas.rds/nueva_tabla_translocada.rds")



#merge the clinical and molecular data using `merge`. I use `geo_accession` as the reference column that contains the patient IDs.
tabla_unificada <- merge.data.frame(feno1, genes_translocados, by.x = "geo_accession", by.y = "row.names", all = TRUE)
View(tabla_unificada)

#save the unified table as an RDS file
saveRDS(tabla_unificada, "Tablas.rds/tabla_unificada.rds")




###############################################################################
# Logistic regression between biochemical recurrence and selected genes. Then creat a table with the estimate and p-value for each gene with the aim of making a plot.

library(stats)

# Select the names of the genes
genes <- colnames(tabla_unificada)[23:28]
genes

# Initialize a list to store the models
modelos <- list()
print(colnames(tabla_unificada))

# Loop to fit a logistic regression model for each gene individually
for (gen in genes) {
  formula <- as.formula(paste("biochemical.recurrence ~", gen))  # Convert string to formula
  modelos[[gen]] <- glm(formula, data = tabla_unificada, family = binomial)
  
  # Show model summary
  print(summary(modelos[[gen]]))
}


# Create a table with three columns: Gene, Estimate, and P-Value

# Initialize an empty list to store the results
resultados <- list()

# Iterate over the genes to extract coefficient and p-value
for (gen in genes) {
  modelo <- modelos[[gen]]
  resumen <- summary(modelo)
  
  # Extract the coefficient (Estimate) and p-value (Pr(>|z|))
  coeficiente <- coef(resumen)[gen, "Estimate"]
  p_valor <- coef(resumen)[gen, "Pr(>|z|)"]
  
  # Save to the list
  resultados[[gen]] <- data.frame(Gene = gen, Estimate = coeficiente, P_Value = p_valor)
}

# Convert the list to a final data frame
tabla_gen_est_pv <- do.call(rbind, resultados)  # Combine list of data frames into one
rownames(tabla_gen_est_pv) <- NULL  # Remove row names if necessary


# Plotting
library(ggplot2)

# Create the plot
ggplot(tabla_gen_est_pv, aes(x = Estimate, y = reorder(Gene, Estimate), fill = P_Value)) +
  geom_col(width = 0.4, color = "black") +  # Narrower bars
  scale_fill_gradient(low = "black", high = "white", name = "P-Value") +  # Fill gradient
  scale_x_continuous(limits = c(-0.05, 0.05)) +
  theme_minimal() +  # Clean theme
  labs(x = "Estimate", y = NULL, title = "Regression Coefficients and Significance") +
  theme(axis.text.y = element_blank()) +  # Remove Y axis labels
  geom_text(aes(label = Gene), hjust = ifelse(tabla_gen_est_pv$Estimate > 0, -0.2, 1.2), size = 4) +  # Gene names next to bars
  theme(legend.position = "right")  
# Place legend to the right

# Save the plot
ggsave("Plots/tablaGenesEstimatePvalue.png", width = 7, height = 5, bg = "white")               


###################################################################################
# Logistic regression between biochemical.recurrence and clinical data + plots

# Initialize the variable 'datos' with the clinical data OF INTEREST
datos <- c("age", "cm", "gleason.sum", "isup", "lapc","ln.percent","pn", "postop.psa", "preop.psa", "pt", "r.residual_tumor", "risk.group","tmprss2.erg")
library(stats)

# Initialize the list of models
modelos_clinicos <- list()

# Loop to fit logistic regression models for each variable individually
for (dato in datos) {
  formula <- as.formula(paste("biochemical.recurrence ~", dato))  # Convert string to formula
  modelos_clinicos[[dato]] <- glm(formula, data = tabla_unificada, family = binomial)
  
  # Show model summary
  print(summary(modelos_clinicos[[dato]]))
}

# Initialize an empty list to store results
# Create a list to store the results
resultados <- list()

# Iterate over the variables and extract coefficient and p-value
for (dato in datos) {
  modelo <- modelos_clinicos[[dato]]
  resumen <- summary(modelo)
  coeficientes <- coef(resumen)  # Extract the coefficients table
  
  # Search for the correct coefficient name (handling inconsistent names)
  coef_nombres <- rownames(coeficientes)
  coef_indice <- grep(paste0("^", dato), coef_nombres)  # Find the correct row
  
  if (length(coef_indice) > 0) {  # Check if coefficient was found
    coeficiente <- coeficientes[coef_indice[1], "Estimate"]
    p_valor <- coeficientes[coef_indice[1], "Pr(>|z|)"]
    
    # Save into the list
    resultados[[dato]] <- data.frame(Dato = dato, Estimate = coeficiente, P_Value = p_valor)
  } else {
    cat("⚠️ Warning: Coefficient not found for", dato, "\n")
  }
}

# Convert the list into a final data frame
tabla_final <- do.call(rbind, resultados)
rownames(tabla_final) <- NULL  # Remove row names if they're causing issues
view(tabla_final)

# Plot
library(ggplot2)

ggplot(tabla_final, aes(x = Estimate, y = reorder(Dato, Estimate), fill = P_Value)) +
  geom_col(width = 0.6, color = "black") +  # Add black border for better visibility
  scale_x_continuous(expand = c(0, 0)) + # make sure the X scale fits the values tightly
  scale_fill_gradient(low = "black", high = "white", name = "P-Value") +  # Color gradient
  theme_minimal() +  # Clean theme
  labs(x = "Estimate", y = NULL, title = "Regression Coefficients of Clinical Variables") +
  theme(axis.text.y = element_text(size = 10)) +  # Restore Y axis labels
  theme(legend.position = "right")  # Place legend on the right

# Save the plot
ggsave(filename = "Plots/RegresionLogistica_BRC_DatosClinicos.png", width = 7, height = 5, bg = "white")

# Same plot as above but highlighting p-value differences more clearly with colors

# Plot with p-value >= 0.05 in white and <0.05 in darker/lighter color
# p >= 0.05 is assigned as NA, then use na.value = "white" in the color scale
# p < 0.05 uses the actual P_Value for coloring
# because scale_fill_gradient() requires numeric values

tabla_final$color <- ifelse(tabla_final$P_Value >= 0.05, NA, tabla_final$P_Value)

ggplot(tabla_final, aes(x = Estimate, y = reorder(Dato, Estimate), fill = color)) +
  geom_col(width = 0.6, color = "black") +  # Black border for all bars
  scale_x_continuous(expand = c(0, 0)) +  # Adjust X scale
  scale_fill_gradient(low = "pink", high = "red", name = "P-Value", na.value = "white") +  
  theme_minimal() +  
  labs(x = "Estimate", y = NULL, title = "Regression Coefficients of Clinical Variables") +
  theme(axis.text.y = element_text(size = 10), legend.position = "right")

# Save the plot
ggsave(filename = "Plots/RegresionLogistica_BRC_DatosClinicosConColores.png", width = 14, height = 5, bg = "white")


view(feno1)


###################################################################################
# Differential gene expression analysis (DGE)
# I analyze biochemical recurrence against the whole gene expression data

# I must use the table 'genes1'(for continue analizing cases with primary tumour and RNAseq) which keeps static, and 'feno1' transposed with only the 'geo_accession' variable

# Rows contain gene names, columns contain sample names
genesDGE <- genes1
View(genesDGE)
rownames(genesDGE) <- as.character(genesDGE$GeneID)  # gene names
genesDGE <- genesDGE[,-1]

# Keep only sample names and biochemical recurrence
fenoDGE <- data.frame(muestra = rownames(feno1), biochemical.recurrence = feno1$biochemical.recurrence)
View(fenoDGE)

saveRDS(fenoDGE, "Tablas.rds/fenoDGE.rds")
saveRDS(genesDGE,"Tablas.rds/genesDGE.rds")

# We need an expression matrix called 'count_norm' (already normalized) and a vector indicating the condition/group of each patient, called 'conditions', which should be a factor

count_norm <- genesDGE
conditions <- fenoDGE$biochemical.recurrence

# Ensure that count_norm is a data frame
count_norm <- as.data.frame(count_norm)

# Run the Wilcoxon rank-sum test for each gene
pvalues <- sapply(1:nrow(count_norm), function(i){
  data <- cbind.data.frame(gene = as.numeric(t(count_norm[i,])), conditions)
  p <- wilcox.test(gene~conditions, data)$p.value
  return(p)
})
# cbind.data.frame(...): creates a data frame with two columns:
# gene: gene expression
# conditions: condition of each sample

# Common adjustment when performing multiple statistical tests to reduce false positives
fdr <- p.adjust(pvalues, method = "fdr")

# Calculate the fold-change for each gene

conditionsLevel <- levels(conditions)
dataConBRF <- count_norm[,c(which(conditions==conditionsLevel[1]))]
dataConBCR <- count_norm[,c(which(conditions==conditionsLevel[2]))]

# Fold-change is calculated as the ratio of the means
foldChanges <- log2(rowMeans(dataConBCR)/rowMeans(dataConBRF))
# Calculates log2 of the fold-change gene by gene
# If the result is >0, the gene is overexpressed in BCR
# If <0, the gene is underexpressed in BCR

# Output results based on the FDR threshold 0.05
outRst <- data.frame(log2foldChange = foldChanges, pValues = pvalues, FDR = fdr)
rownames(outRst) <- rownames(count_norm)
outRst <- na.omit(outRst)
view(outRst)

saveRDS(outRst, "Tablas.rds/outRst.rds")

# How many genes have FDR less than 0.05?

# Genes with FDR < 0.05
significativos_005 <- which(fdr < 0.05, arr.ind = TRUE)
cantidad_significativos_005 <- length(significativos_005) 
cantidad_significativos_005 #1253

# Which are the 20 genes with the lowest FDR? --> the most significant ones!!!
orden_fdr <- order(fdr, na.last = NA)
veinte_menores_FDR <- orden_fdr[1:20]
veinte_menores_FDR

# Which are the 20 genes with the highest absolute log2 fold-change?
orden_logFC <- order(abs(foldChanges), decreasing = TRUE, na.last = NA)
veinte_mayores_LFC <- orden_logFC[1:20]
veinte_mayores_LFC

# Summarize these results: volcano plot

library(ggplot2)
library(ggrepel)  # To add labels without overlapping

# Create a column to classify genes based on FDR and logFC
outRst$Significativo <- ifelse(outRst$FDR < 0.05 & abs(outRst$log2foldChange) > 1, 
                               "Significant", "Not significant")

# Get the 20 most significant genes, i.e., with the lowest FDR
top_genes <- rownames(outRst)[orden_fdr[1:20]]

# Volcano plot
outRst$gene <- rownames(outRst)

volcanoplot <- ggplot(outRst, aes(x = log2foldChange,
                                  y = -log10(FDR),
                                  color = Significativo,
                                  text = paste0("Gene: ", gene,
                                                "\nlog2FC: ", round(log2foldChange, 2),
                                                "\nFDR: ", signif(FDR, 3)))) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray")) +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-log10(FDR)") +
  theme_minimal() +
  geom_text_repel(data = outRst[top_genes, ],
                  aes(label = gene),
                  size = 3, max.overlaps = 10)

ggsave("Plots/volcano_plot.png", width = 8, height = 6, dpi = 300)

library(plotly)
ggplotly(volcanoplot, tooltip = "text")

