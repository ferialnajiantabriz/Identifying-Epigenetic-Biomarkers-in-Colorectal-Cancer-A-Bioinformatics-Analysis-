if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")# Install and load the annotation package
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("karyoploteR")
BiocManager::install("DMRcate")
BiocManager::install("pheatmap")
BiocManager::install("genomation")
BiocManager::install("methylKit")
BiocManager::install("GenomicRanges", force = TRUE)
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
BiocManager::install("WGCNA" )
BiocManager::install("flashClust")
install.packages("ggforce")
install.packages("VennDiagram")
install.packages("ggvenn")
install.packages("WGCNA")
install.packages("withr")
install.packages("withr")
install.packages("WGCNA", force = TRUE)
install.packages("writexl")
install.packages("RIdeogram")
install.packages("ggrepel")
install.packages("readr")
library(dplyr)
library(karyoploteR)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg19)
library(VennDiagram)
library(ggvenn)
library(WGCNA)
library(flashClust)
library(curl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(methylKit)
library(genomation)
library(readr)
library(ggrepel)
library(pheatmap)
library(RIdeogram)
library(writexl)
library(DMRcate)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(GEOquery)
library(limma)
library(minfi)
library(ggplot2)
library(Gviz)
library(GenomicRanges)
library(RColorBrewer)
library(GEOquery)
library(edgeR)
library(limma)
library(ggplot2)
library(withr)

#############################
#Downloading the GEO Dataset
gse <- getGEO("GSE75546", GSEMatrix = TRUE)
str(gse)

# Extracting methylation data and metadata
meth_data <- exprs(gse[[1]]) 
meth_metadata <- pData(gse[[1]])
head(meth_data)
head(meth_metadata)

# Saving to files
write.csv(meth_data, file = "GSE75546_methylation_data.csv")
write.csv(meth_metadata, file = "GSE75546_sample_metadata.csv")
saveRDS(meth_data, file = "GSE75546_methylation_data.RDS")
saveRDS(meth_metadata, file = "GSE75546_sample_metadata.RDS")
# Create a zip folder containing the CSV file
zip("Project supplementary files.zip", files = "GSE75546_methylation_data.csv")
zip("Project supplementary files.zip", files = "GSE75546_sample_metadata.csv")

#Cleaning the data and filtering probes 
#Removing probes with missing values
probe_na <- rowSums(is.na(meth_data))
methylation_data <- meth_data[probe_na == 0, ]

# filitering probes on chromosome X and Y
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
keep_probes <- !(rownames(methylation_data) %in% ann450k$Name[ann450k$chr %in% c("chrX", "chrY")])
methylation_data <- methylation_data[keep_probes, ]

#Making sure there are no probes remaining on chromosome X and Y
x_y_probes <- ann450k$Name[ann450k$chr %in% c("chrX", "chrY")]
remaining_probes <- rownames(methylation_data) %in% x_y_probes
cat("Number of chrX/Y probes remaining:", sum(remaining_probes), "\n")

#filtering probes overlapping with SNPs
no_snp_probes <- ann450k$Name[is.na(ann450k$Probe_rs)]
methylation_data <- methylation_data[rownames(methylation_data) %in% no_snp_probes, ]

# Converting beta values to M-values for differential analysis
m_values <- t(apply(methylation_data, 1, function(x) log2(x / (1 - x))))

#savng the data
saveRDS(m_values, file = "GSE75546_m_values.RDS")
saveRDS(methylation_data, file = "GSE75546_b_values.RDS")

#Metadata Processing and Filtering
#Loading the metadata
metadata <- read.csv("GSE75546_sample_metadata.csv")

#Grouping variables to normal and cancer 
metadata$sample_type <- ifelse(metadata$source_name_ch1 == "adjacent normal tissue", "normal", "rectal cancer")

# Filtering the metadata for the samples
metadata <- metadata[metadata$sample_type %in% c("normal", "rectal cancer"), ]

# Assigning Row Names to the Metadata
row.names(metadata) <- metadata$geo_accession

#Loading the RDS files 
bval <- readRDS("GSE75546_b_values.RDS")
mval <- readRDS("GSE75546_m_values.RDS")

# Ensuring the data samples in the metadata and methylation data match 
matching_indices <- match(row.names(metadata), colnames(bval))
bval <- bval[, matching_indices, drop = FALSE]
mval <- mval[, matching_indices, drop = FALSE]
metadata <- metadata[!is.na(matching_indices), ]

# Differentially methlated CPG sites analysis
# Creating a design matrix required for linear modelling
metadata$sample_type <- factor(metadata$sample_type, levels = c("normal", "rectal cancer"))
design <- model.matrix(~ sample_type, data = metadata)

# Performing Differential Methylation Analysis using Linear Modeling and Empirical Bayes Moderation
fit <- lmFit(mval, design)
fit2 <- eBayes(fit)

# Extracting the top differentially methylated probes
deff.meth <- topTable(fit2, coef = 2, sort.by = "p", number = nrow(mval), adjust.method = "BY")

# Boxplot for top CpGs to visualize the differential methylation between normal and rectal cancer samples
png("boxplot_deff_meth.png", width = 3000, height = 3000)
par(mfrow = c(2, 5), mar = c(5, 4, 2, 1))

plot_custom_cpg <- function(bval, cpg, pheno) {
  beta_values <- bval[cpg, ]
  boxplot(beta_values ~ pheno, 
          ylab = "Beta values", xlab = "", cex.main = 3, cex.lab = 3, cex.axis = 3,
          main = cpg, col = c("skyblue", "pink"),
          xaxt = 'n')
  mtext(side = 1, at = 1:2, text = c("Normal", "Rectal Cancer"), line = 2.5, cex = 2)
}
sapply(rownames(deff.meth)[1:10], function(cpg) {
  plot_custom_cpg(bval = bval, cpg = cpg, pheno = metadata$sample_type)
})

dev.off()
# Volcano Plot of DMCs
dat <- data.frame(
  foldchange = deff.meth$logFC,
  logPvalue = -log10(deff.meth$P.Value)
)
dat$category <- "Non-significant"
dat$category[dat$foldchange > 0.4 & dat$logPvalue > -log10(0.05)] <- "Hypermethylation"
dat$category[dat$foldchange < -0.4 & dat$logPvalue > -log10(0.05)] <- "Hypomethylation"
###########################################
deff.meth$ID <- rownames(deff.meth)
deff.meth <- deff.meth[, !colnames(deff.meth) %in% "ID"]
# Define thresholds for categorization
logFC_threshold <- 0.4  # Match the threshold in 'dat'
p_value_threshold <- 0.05

# Add a 'Category' column based on conditions
deff.meth$Category <- ifelse(
  deff.meth$logFC > logFC_threshold & deff.meth$P.Value < p_value_threshold, "Hypermethylated",
  ifelse(
    deff.meth$logFC < -logFC_threshold & deff.meth$P.Value < p_value_threshold, "Hypomethylated",
    "Non-significant"
  )
)

# Write the 'deff.meth' data frame with the 'Category' column to a CSV file
write.csv(deff.meth, file = "differential_methylation_results.csv", row.names = TRUE)
zip("Project supplementary files.zip", files = "differential_methylation_results.csv")

# View the updated data frame
head(deff.meth)

##################################################
cols <- c("Non-significant" = "grey", "Hypermethylation" = "red", "Hypomethylation" = "blue")

ggplot(data = dat, aes(x = foldchange, y = logPvalue, color = category)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_colour_manual(values = cols, name = "Methylation Status") +
  geom_vline(xintercept = c(-0.4, 0.4), colour = "#990000", linetype = "dashed") +
  xlab("Log Fold Change") +
  ylab("-log10 p-value") +
  theme_bw() +
  theme(legend.position = "right")

ggsave("volcano_plot3.png", width = 10, height = 8)

#Setting the thresholds
p_value_threshold <- 0.05
logFC_threshold <- 0.4

#Filtering Differentially Methylated Sites:
total_screened <- nrow(deff.meth)
differentially_methylated <- deff.meth[deff.meth$P.Value < p_value_threshold & abs(deff.meth$logFC) > logFC_threshold, ]
total_differentially_methylated <- nrow(differentially_methylated)
hypermethylated <- nrow(differentially_methylated[differentially_methylated$logFC > logFC_threshold, ])
hypomethylated <- nrow(differentially_methylated[differentially_methylated$logFC < -logFC_threshold, ])

#Printing Results
cat("Total CpG sites screened:", total_screened, "\n")
cat("Total differentially methylated CpG sites:", total_differentially_methylated, "\n")
cat("Hypermethylated CpG sites:", hypermethylated, "\n")
cat("Hypomethylated CpG sites:", hypomethylated, "\n")



# Identifying Differentially Methylated Region

# 'metadata$sample_type' contains group information (normal vs rectal cancer)
pheno <- metadata$sample_type
pheno
# Creating the design matrix based on normal and rectal cancer groups
design <- model.matrix(~ pheno)
head(design)
# Setting the annotation for differential analysis 
myAnnotation <- cpg.annotate(object = mval, datatype = "array", 
                             what = "M", 
                             analysis.type = "differential", 
                             design = design, 
                             contrasts = FALSE, 
                             coef = 2,  # The coefficient for rectal cancer in the design matrix
                             arraytype = "450K", 
                             fdr = 0.001)  # False discovery rate threshold

# Checking the structure of the annotation object
str(myAnnotation)

# Performing DMR analysis (Differentially Methylated Regions)
DMRs <- dmrcate(myAnnotation, lambda = 1000, C = 2)

# Setting the ExperimentHub cache directory to desired location
Sys.setenv(EXPERIMENT_HUB_CACHE = "C:/Users/olada/OneDrive/dmcs")

# Extracting the results as genomic ranges
results.ranges <- extractRanges(DMRs)

# View the DMR results
print(results.ranges)

results_df <- as.data.frame(results.ranges)
print(results_df)
# Write to Excel file
write_xlsx(results_df, "results_ranges.xlsx")
zip("Project supplementary files.zip", files = "results_ranges.xlsx")
################################################################
# Add methylation type based on meandiff
mcols(results.ranges)$type <- ifelse(
  mcols(results.ranges)$meandiff > 0, 
  "hyper", 
  "hypo"
)

# Check the counts of hyper and hypo regions
table(mcols(results.ranges)$type)

# Trim ranges to fit within chromosome lengths
library(BSgenome.Hsapiens.UCSC.hg19)
seqinfo(results.ranges) <- seqinfo(Hsapiens)
results.ranges <- trim(results.ranges)
results.ranges <- keepStandardChromosomes(results.ranges, pruning.mode = "coarse")
# Remove unused levels from seqlevels
seqlevels(results.ranges) <- seqlevelsInUse(results.ranges)

# Verify unique seqnames again
unique(seqnames(results.ranges))

# Check for any out-of-bound ranges
out_of_bounds <- results.ranges[
  start(results.ranges) > seqlengths(results.ranges)[as.character(seqnames(results.ranges))] |
    end(results.ranges) > seqlengths(results.ranges)[as.character(seqnames(results.ranges))]
]
print(out_of_bounds)  # Should return an empty GRanges object

# Save high-resolution karyogram plot
png("DMR_karyogram.png", width = 4000, height = 3000, res = 400)

# Create the karyogram plot
kp <- plotKaryotype(genome = "hg19")
kpPlotRegions(
  kp,
  data = results.ranges,
  col = ifelse(mcols(results.ranges)$type == "hyper", "red", "blue"),
  border = NA  # Optional: Remove borders for cleaner visualization
)

dev.off()

#####################################################################
#Differentially expressed Genes
## Extracting Expression Data and Metadata
# Specify the GEO dataset ID
rectal_cancer <- "GSE75548"  # GEO dataset for rectal cancer samples

#Download the dataset from GEO
gse <- getGEO(rectal_cancer)  # Retrieve the dataset as an object
length(gse)  # Check the number of series (platforms) in the dataset

#Extract the first series from the dataset
gse <- gse[[1]]  # Select the first series for analysis

#Overview of the dataset
gse  # Print a summary of the dataset (e.g., number of samples, metadata)

# Access metadata (phenotypic information)
pheno_metadata <- pData(gse)  # Extract sample annotations
head(pheno_metadata)  # Display the first few rows of the sample metadata

#Access feature metadata (probe or gene information)
feature_metadata <- fData(gse)  # Extract annotations for features (e.g., probes or genes)
head(feature_metadata)  # Display the first few rows of the feature metadata

#Access the expression/methylation data matrix
data_matrix <- exprs(gse)  # Extract the numerical expression or beta-value matrix
dim(data_matrix)  # Check the dimensions of the matrix (features x samples)
head(data_matrix)  # Display the first few rows of the data matrix

#Summarize the data values
summary(data_matrix)  # Get summary statistics (min, max, median, etc.) for the data matrix
write.csv(pheno_metadata, file = "rectal_cancer_metadata.csv")  # Save metadata to CSV
write.csv(feature_metadata, file = "rectal_cancer_feature_metadata.csv")  # Save feature metadata to CSV
write.csv(data_matrix, file = "rectal_cancer_data_matrix.csv")  # Save data ma

zip("Project supplementary files.zip", files = "rectal_cancer_metadata.csv")
zip("Project supplementary files.zip", files = "rectal_cancer_feature_metadata.csv")
zip("Project supplementary files.zip", files = "rectal_cancer_data_matrix.csv")
#Boxplot to visualize data distributions
boxplot(data_matrix, outline = FALSE, 
        main = "Distribution of Expression/Beta Values Across Samples", 
        ylab = "Expression/Beta Values", 
        xlab = "Samples", 
        col = "lightblue")  # Create a boxplot to examine sample-wise data distributions
# Select relevant columns and rename them
sampleInfo <- dplyr::select(pheno_metadata, source_name_ch1, characteristics_ch1) %>%
  rename(group = source_name_ch1, patient = characteristics_ch1)

# View the updated dataset
sampleInfo
# Convert the group column to a factor with the correct levels
sampleInfo$group <- factor(sampleInfo$group, levels = c("adjacent normal tissue", "rectal cancer"))

# Now check the levels
levels(sampleInfo$group)


# Calculate the correlation matrix for the expression data
corMatrix <- cor(data_matrix, use = "c")

# Plot the heatmap of the correlation matrix
pheatmap(corMatrix)

# Align row names of metadata with the correlation matrix sample names
rownames(sampleInfo) <- colnames(corMatrix)

# Plot the heatmap with sample annotations
pheatmap(corMatrix,
         annotation_col = sampleInfo)

# Perform PCA on the transposed expression data
pca <- prcomp(t(data_matrix))

# Visualize the PCA results using ggplot2
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y = PC2, col = group, label = paste("Patient", patient))) + 
  geom_point() + 
  geom_text_repel()

# Combine feature data with expression data
combined_data <- cbind(feature_metadata, data_matrix)
head(combined_data)
write_csv(combined_data, file = "Expression_combined_data.csv")
zip("Project supplementary files.zip", files = "Expression_combined_data.csv")
# Select specific columns and rename `ID` to `illumna_ID`
features_selected <- feature_metadata %>%
  dplyr::select(ID, ILMN_Gene, Source_Reference_ID, Symbol, Entrez_Gene_ID, Chromosome, Cytoband) %>%
  rename(illumna_ID = ID)

# Combine again with expression data
combined_data <- cbind(features_selected, data_matrix)
head(combined_data)
# Save the combined data to CSV
write.csv(combined_data, file = "Expression_combined_data.csv", row.names = TRUE)
zip("Project supplementary files.zip", files = "Expression_combined_data.csv")

#Differential Expression
#Create the design matrix
design <- model.matrix(~0 + sampleInfo$group)
colnames(design) <- c("Normal", "Rectal_Cancer")
design

#Filter the data based on expression levels
summary(data_matrix)
cutoff <- median(data_matrix)
is_expressed <- data_matrix > cutoff
keep <- rowSums(is_expressed) > 2
table(keep)
gse <- gse[keep,]

#Perform differential expression analysis
fit <- lmFit(data_matrix, design)
contrast <- makeContrasts(Rectal_Cancer - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

#Extract top genes
top_genes <- topTable(fit2, number = 10, coef = 1)
head(top_genes)

#Get all results and add Gene_ID column
deg_result <- topTable(fit2, number = Inf)
deg_result <- tibble::rownames_to_column(deg_result, "Gene_ID")
nrow(deg_result)
str(deg_result)

#Filter the DEGs based on significance and fold change
p_value_threshold <- 0.05
logFC_threshold <- 0.5
DEGs <- deg_result[deg_result$P.Value < p_value_threshold & abs(deg_result$logFC) > logFC_threshold, ]
nrow(DEGs)
# Add a 'Category' column based on conditions
deg_result$Category <- ifelse(
  deg_result$logFC > logFC_threshold & deg_result$P.Value < p_value_threshold, "Upregulated",
  ifelse(
    deg_result$logFC < -logFC_threshold & deg_result$P.Value < p_value_threshold, "Downregulated",
    "Non-significant"
  )
)
# Ensure the `illumna_ID` column is aligned with the `Gene_ID` column
colnames(combined_data)[colnames(combined_data) == "illumna_ID"] <- "Gene_ID"

# Merge `deg_result` with `combined_data` to add the `Symbol` column
deg_result_with_symbols <- merge(deg_result, combined_data[, c("Gene_ID", "Symbol")], 
                                 by = "Gene_ID", all.x = TRUE)

# Save the updated results to a CSV file
write.csv(deg_result_with_symbols, file = "DEG_results_with_symbols.csv", row.names = TRUE)

# View the first few rows of the updated data
head(deg_result_with_symbols)

# Save the results with categories to a CSV file
write.csv(deg_result, file = "deg_result_with_symbols.CSV", row.names = TRUE)
zip("Project supplementary files.zip", files = "deg_result_with_symbols.CSV")
head(deg_result)
# Optional: Display a summary of the categories
table(deg_result$Category)
#Count upregulated and downregulated DEGs
upregulated_DEGs <- nrow(DEGs[DEGs$logFC > logFC_threshold, ])
downregulated_DEGs <- nrow(DEGs[DEGs$logFC < -logFC_threshold, ])

#Print summary statistics
cat("Total DEGs:", nrow(DEGs), "\n")
cat("Upregulated DEGs:", upregulated_DEGs, "\n")
cat("Downregulated DEGs:", downregulated_DEGs, "\n")


#Visualize DEGs with a Volcano Plot
deg_result$Category <- ifelse(
  deg_result$logFC > logFC_threshold & deg_result$P.Value < p_value_threshold, "Upregulated",
  ifelse(deg_result$logFC < -logFC_threshold & deg_result$P.Value < p_value_threshold, "Downregulated", "Non-significant")
)

ggplot(deg_result, aes(x = logFC, y = -log10(P.Value), color = Category)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Non-significant" = "gray")) +
  geom_hline(yintercept = -log10(p_value_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "Log2 Fold Change",
       y = "-log10 P-value") +
  theme_minimal(base_size = 15) +
  theme(plot.background = element_rect(fill = "white"),
        legend.position = "right")
ggsave("DEG_volcano_plot_new1.png", width = 10, height = 8)

#Filter the differential expression results
filtered_DEGs <- filter(deg_result_with_symbols, adj.P.Val < 0.05, abs(logFC) > 0.5)

#Check the number of filtered DEGs
nrow(filtered_DEGs)

#Save the filtered DEGs to a CSV file
write_csv(filtered_DEGs, file = "filtered_DEGs_results.csv")
zip("Project supplementary files.zip", files = "filtered_DEGs_results.csv")

#Ensure `filtered_DEGs` and `full_output_gene` have common identifiers
if (!"Gene_ID" %in% colnames(filtered_DEGs) || !"illumna_ID" %in% colnames(deg_result)) {
  stop("Common identifier (Gene_ID or illumna_ID) is missing in one of the dataframes.")
}

#Merge `filtered_DEGs` with `full_output_gene` to add gene symbols
filtered_DEGs <- merge(
  filtered_DEGs,
  combined_data[, c("illumna_ID", "Symbol")],  # Select relevant columns
  by.x = "Gene_ID",                              # Match `Gene_ID` in filtered_DEGs
  by.y = "illumna_ID",                           # Match `illumna_ID` in full_output_gene
  all.x = TRUE                                   # Keep all rows in filtered_DEGs
)
#Ensure filtered_DEGs contains Gene Symbols
if (!"Symbol" %in% colnames(filtered_DEGs)) {
  stop("The Symbol column is missing in the filtered_DEGs dataframe.")
}
colnames(combined_data)
colnames(filtered_DEGs)

#Replace row names of filtered_DEGs with gene symbols
filtered_DEGs$Symbol <- make.unique(ifelse(filtered_DEGs$Symbol == "", filtered_DEGs$Gene_ID, filtered_DEGs$Symbol))

#Subset the expression matrix using Gene_IDs and update row names with Symbols
deg_ids <- filtered_DEGs$Gene_ID  # Illumina IDs for subsetting rows
gene_matrix <- data_matrix[deg_ids, ]  # Subset rows of the expression matrix
rownames(gene_matrix) <- filtered_DEGs$Symbol  # Replace row names with gene symbols

#Generate Heatmap
png("DEG_Heatmap_GeneSymbols_new1.png", width = 4000, height = 12000, res = 400)

pheatmap(
  gene_matrix,
  scale = "row",               # Scale rows for better visualization
  cluster_cols = TRUE,         # Cluster columns (samples)
  cluster_rows = TRUE,         # Cluster rows (genes)
  show_rownames = TRUE,        # Display gene symbols as row names
  show_colnames = TRUE,        # Display sample names as column names
  fontsize_row = 3,            # Adjust row font size for better visibility
  fontsize_col = 8,            # Adjust column font size for better visibility
  main = "Heatmap of DEGs"  # Add a title to the heatmap
)

dev.off()

dev.list()
########################################


##########. MRGS. #############


#
#
###genes impacted by methylation(methylationDS).
results_df <- readRDS("/Users/ferialnajiantabriz/Desktop/Bioproject/results_df.RDS")

full_results <- readRDS("/Users/ferialnajiantabriz/Desktop/Bioproject/full_results.RDS") ####genes with altered expression(RNADS)
head(full_results)
View(methyl_regulated_genes_df)

# Check if 'overlapping.genes' column exists in `results_df`
if (!"overlapping.genes" %in% colnames(results_df)) {
  stop("The `results_df` does not contain the 'overlapping.genes' column.")
}

# Check if 'Symbol' column exists in `full_results`
if (!"Symbol" %in% colnames(full_results)) {
  stop("The `full_results` does not contain the 'Symbol' column.")
}

# Split overlapping.genes into separate gene names and remove repetetive genes
split_genes <- unique(unlist(strsplit(results_df$overlapping.genes, ", ")))

# Ensure the split_genes are of character type for comparison
split_genes <- as.character(split_genes)

# Find common genes between the split list and full_results
#it identifies the genes that are both differentially methylated and differentially expressed.
common_genes <- intersect(split_genes, full_results$Symbol)
# Filter the full_results for the common genes
#filter the RNA (gene expression) dataset (full_results) to keep only the rows that contain the common genes we identified. 
methyl_regulated_genes <- full_results[full_results$Symbol %in% common_genes, ]

# Display the result
#display the MRGs (Methylation-Regulated Genes) table first few rows
if (nrow(methyl_regulated_genes) > 0) {
  head(methyl_regulated_genes)
} else {
  cat("No methyl-regulated genes found with the current criteria.")
}

methyl_regulated_genes_df2 <- data.frame(methyl_regulated_genes)

head(methyl_regulated_genes_df)

#  save the DataFrame to a CSV file
write.csv(methyl_regulated_genes_df, "Methyl_Regulated_Genes.csv", row.names = FALSE)


head(full_results)
View(methyl_regulated_genes_df)
###################################################KEGG


###############################################################
#
#
#
#
#
#
#Function enrichment analysis
#
#
#
#
#
#
# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"), ask = FALSE)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Load methyl_regulated_genes_df from your CSV file
methyl_regulated_genes_df <- read.csv("Methyl_Regulated_Genes.csv")

# Ensure that the 'Symbol' column exists in methyl_regulated_genes_df
if (!"Symbol" %in% colnames(methyl_regulated_genes_df)) {
  stop("The `methyl_regulated_genes_df` does not contain the 'Symbol' column.")
}

# Extract the gene symbols for enrichment analysis
gene_symbols <- methyl_regulated_genes_df$Symbol

# Map gene symbols to Entrez IDs
gene_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Identify and display any unmapped genes
unmapped_genes <- setdiff(gene_symbols, gene_ids$SYMBOL)
if (length(unmapped_genes) > 0) {
  cat("Warning: The following gene symbols could not be mapped to Entrez IDs and will be excluded from the analysis:\n")
  print(unmapped_genes)
} else {
  cat("All gene symbols were successfully mapped.\n")
}

# Use only the successfully mapped Entrez IDs
mapped_entrez_ids <- gene_ids$ENTREZID

#  Perform KEGG Pathway Enrichment Analysis
kegg_results <- enrichKEGG(
  gene = mapped_entrez_ids,
  organism = "hsa",  
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

#Perform GO Enrichment Analysis
go_results <- enrichGO(
  gene = mapped_entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",  # Includes BP, MF, and CC categories
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

#  Display Results
cat("\nKEGG Enrichment Results:\n")
print(head(kegg_results))

cat("\nGO Enrichment Results:\n")
print(head(go_results))

# Visualize KEGG and GO Results (Optional)
# Dotplot visualization for KEGG results
if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
  dotplot(kegg_results, showCategory = 10) + ggtitle("KEGG Pathway Enrichment")
} else {
  cat("No significant KEGG enrichment results to display.\n")
}

# Dotplot visualization for GO results
if (!is.null(go_results) && nrow(go_results) > 0) {
  dotplot(go_results, showCategory = 10) + ggtitle("GO Enrichment Analysis")
} else {
  cat("No significant GO enrichment results to display.\n")
}

######################################
####################. WCGNA ########################################

######################################### WGCNA #####################




if(!require("BiocManager", quietly=TRUE)) install.packages("BiocManager")
packages <- c("GEOquery", "WGCNA", "pROC", "limma")
for (pkg in packages) {
  if(!require(pkg, character.only=TRUE, quietly=TRUE)) {
    BiocManager::install(pkg, ask=FALSE, update=TRUE)
    library(pkg, character.only=TRUE)
  }
}

allowWGCNAThreads()

if(!dir.exists("WGCNA_Results")) dir.create("WGCNA_Results")

## Step 1: Load and Preprocess Gene Expression Data
file_path <- "/Users/ferialnajiantabriz/Desktop/gse_full_output_gene.csv"
expression_data <- read.csv(file_path, check.names = FALSE)

# Ensure gene symbols as rownames
if ("Symbol" %in% colnames(expression_data)) {
  rownames(expression_data) <- make.unique(expression_data$Symbol)
  expression_data <- expression_data[, !colnames(expression_data) %in% "Symbol"]
} else {
  stop("No 'Symbol' column found in expression data. Please ensure gene symbols are available.")
}

# Remove gene annotation columns 
metadata_cols <- c("Entrez_Gene_ID", "Chromosome", "Cytoband")
metadata_cols <- metadata_cols[metadata_cols %in% colnames(expression_data)]
datExpr <- expression_data[, !colnames(expression_data) %in% metadata_cols]

# Transpose: samples as rows, genes as columns
datExpr <- t(datExpr)

# Save initial processed expression matrix
write.csv(datExpr, file = "WGCNA_Results/initial_datExpr.csv")


## Step 2: Download and Align Metadata from GEO (GSE75548)

gse_id <- "GSE75548"
gse <- getGEO(gse_id, GSEMatrix = TRUE)
if (length(gse) > 1) {
  # If multiple platforms, pick the first
  gse <- gse[[1]]
} else {
  gse <- gse[[1]]
}

meth_metadata <- pData(gse)

# Align samples in datExpr with metadata
common_samples <- intersect(rownames(meth_metadata), rownames(datExpr))
if (length(common_samples) == 0) {
  stop("No overlapping samples found. Ensure sample names in datExpr match GEO sample names.")
}

datExpr <- datExpr[common_samples, , drop=FALSE]
meth_metadata <- meth_metadata[common_samples, , drop=FALSE]

write.csv(meth_metadata, "WGCNA_Results/meth_metadata_subset.csv")


## Step 3: Define the Trait of Interest (Rectal Cancer vs Normal)

trait_col <- "tissue subtype:ch1"
if(!trait_col %in% colnames(meth_metadata)) {
  stop(paste("Trait column", trait_col, "not found in metadata. Check metadata columns."))
}

trait_factor <- factor(meth_metadata[[trait_col]])
levels(trait_factor) 

# Rectal cancer = 1, adjacent normal = 0
trait_numeric <- ifelse(trait_factor == "rectal cancer", 1, 0)
traits <- data.frame(Condition = trait_numeric)
rownames(traits) <- rownames(meth_metadata)

write.csv(traits, file="WGCNA_Results/traits.csv")


## Step 4: Data Quality Check and Filtering

# Filter out low variance genes
varThreshold <- 1
gene_variances <- apply(datExpr, 2, var)
highVarGenes <- gene_variances > varThreshold
datExpr <- datExpr[, highVarGenes]

write.csv(datExpr, file="WGCNA_Results/datExpr_highVar.csv")

# Check for good samples and genes
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  traits <- traits[rownames(datExpr), , drop=FALSE]
}

# Sample clustering to detect outliers
sampleTree <- hclust(dist(datExpr), method = "average")

pdf("WGCNA_Results/Step4_SampleClustering.pdf")
plot(sampleTree, main="Sample Clustering to detect outliers (Rectal Cancer)", sub="", xlab="")
dev.off()


## Step 5: Pick a Soft-thresholding Power

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")

pdf("WGCNA_Results/Step5_ScaleFreeTopology.pdf")
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     main = "Scale independence",
     type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")
abline(h=0.8, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     main = "Mean connectivity",
     type="n")
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels=powers, cex=0.9, col="red")
dev.off()

chosen_power <- sft$powerEstimate
if (is.na(chosen_power)) {
  chosen_power <- 6  # Fallback if needed
}

write.table(chosen_power, file="WGCNA_Results/chosen_power.txt", row.names=FALSE, col.names=FALSE)


## Step 6: Construct the Co-expression Network and Identify Modules

maxSize <- ncol(datExpr)
net <- blockwiseModules(datExpr, 
                        power = chosen_power,
                        TOMType = "signed", 
                        minModuleSize = 100, 
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        maxBlockSize = maxSize,
                        saveTOMs=TRUE,
                        saveTOMFileBase="WGCNA_Results/TOM",
                        verbose = 3)

moduleColors <- net$colors
MEs <- net$MEs

pdf("WGCNA_Results/Step6_ModuleDendrogram.pdf")
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors,
                    "Module colors",
                    main = "Gene dendrogram and module colors (Rectal Cancer)")
dev.off()

write.csv(moduleColors, file="WGCNA_Results/moduleColors.csv")
write.csv(MEs, file="WGCNA_Results/moduleEigengenes.csv")


## Step 7: Relate Modules to Rectal Cancer Trait

moduleTraitCor <- cor(MEs, traits$Condition, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

pdf("WGCNA_Results/Step7_ModuleTraitHeatmap.pdf")
textMatrix <- paste(signif(moduleTraitCor,2), "\n(",
                    signif(moduleTraitPvalue,1),")", sep="")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "Rectal Cancer Condition",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               main = "Module-trait relationships (Rectal Cancer)")
dev.off()

write.csv(cbind(moduleTraitCor, moduleTraitPvalue), file="WGCNA_Results/moduleTraitCorrelations.csv")

# Identify key modules: p<0.05 and |cor|>0.5
significantModules <- names(which(moduleTraitPvalue < 0.1 & abs(moduleTraitCor) > 0.3))
write.csv(significantModules, file="WGCNA_Results/significantModules.csv", row.names=FALSE)


## Step 8: Calculate Gene Significance (GS) and Module Membership (MM)

GS <- as.data.frame(cor(datExpr, traits$Condition, use = "p"))
colnames(GS) <- "GS"
GS.pvalue <- as.data.frame(corPvalueStudent(as.matrix(GS), nrow(datExpr)))

MM <- as.data.frame(sapply(MEs, function(me) cor(datExpr, me, use="p")))
MM.pvalue <- as.data.frame(sapply(MEs, function(me) corPvalueStudent(cor(datExpr, me, use="p"), nrow(datExpr))))

write.csv(GS, file="WGCNA_Results/geneSignificance.csv")
write.csv(GS.pvalue, file="WGCNA_Results/geneSignificance_pvalues.csv")
write.csv(MM, file="WGCNA_Results/moduleMembership.csv")
write.csv(MM.pvalue, file="WGCNA_Results/moduleMembership_pvalues.csv")


## Step 9: Identify Trait-Associated Genes from Key Modules
## Relaxed Criteria: |GS|>0.4, |MM|>0.6, p ≤ 0.1

key_module_genes <- list()
for (mod in significantModules) {
  inModule <- moduleColors == as.numeric(mod)
  modGenes <- names(which(inModule))
  
  column_ME <- paste0("ME", mod)
  modMM <- MM[modGenes, column_ME, drop=FALSE]
  modMM_p <- MM.pvalue[modGenes, column_ME, drop=FALSE]
  
  modGS <- GS[modGenes, "GS", drop=FALSE]
  modGS_p <- GS.pvalue[modGenes, "GS", drop=FALSE]
  
  # Relaxed thresholds
  selected <- which(abs(modGS$GS) > 0.4 & abs(modMM[,1]) > 0.6 & 
                      modGS_p[,1] < 0.1 & modMM_p[,1] < 0.1)
  
  selected_genes <- modGenes[selected]
  key_module_genes[[mod]] <- selected_genes
  
  cat("Module:", mod, "- Selected genes:", length(selected_genes), "\n")
}

# Rectal Cancer Genes (RCGs)
RCGs <- unique(unlist(key_module_genes))
write.csv(RCGs, file="RCGs2.csv", row.names = FALSE)

# Check summaries again after selection
cat("Summary of GS:\n")
print(summary(GS$GS))

cat("Summary of MM:\n")
print(summary(MM))
print(significantModules)

## Step 10: Load MRGs from CSV and Intersect with RCGs

# Read MRGs from the provided CSV table (make sure this file exists as a proper CSV)
MRGs_data <- read.csv("MRGs_table.csv", stringsAsFactors = FALSE)
MRGs <- unique(MRGs_data$Symbol)  # Extract Symbol column as MRG list

# Rectal Cancer-associated MRGs
RC_MRGs <- intersect(RCGs, MRGs)
write.csv(RC_MRGs, file="WGCNA_Results/RC_MRGs.csv", row.names=FALSE)


## Step 11: ROC Analysis for Diagnostic Genes

if (length(RC_MRGs) > 0) {
  auc_values <- c()
  pdf("WGCNA_Results/Step11_ROC_Curves.pdf")
  par(mfrow=c(2,2))
  
  for (gene in RC_MRGs) {
    if (gene %in% colnames(datExpr)) {
      gene_exp <- datExpr[, gene]
      roc_obj <- roc(traits$Condition, gene_exp, plot=FALSE, quiet=TRUE)
      auc_val <- auc(roc_obj)
      auc_values[gene] <- auc_val
      # Plot ROC curve
      plot(roc_obj, main=paste("ROC for", gene, "\nAUC =", round(auc_val,3)))
    }
  }
  
  dev.off()
  
  # Filter by AUC > 0.8
  diagnostic_genes <- names(auc_values[auc_values > 0.8])
  write.csv(diagnostic_genes, file="WGCNA_Results/diagnostic_genes.csv", row.names=FALSE)
  
  if (length(diagnostic_genes) > 0) {
    cat("Diagnostic genes (AUC > 0.8):", diagnostic_genes, "\n")
  } else {
    cat("No RC-MRGs with AUC > 0.8 found.\n")
  }
} else {
  cat("No RC-MRGs identified. Check MRG list or criteria.\n")
}

## Step 12: Final Notes

cat("WGCNA analysis for rectal cancer completed. Results saved in 'WGCNA_Results' directory.\n")

