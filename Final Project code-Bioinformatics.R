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
#MRGs
getwd()
results_ranges <- read.csv("C:/Users/olada/OneDrive/Documents/results_ranges.csv", header = TRUE)
# read the gene BED file
transcriptBED=system.file("extdata", "refseq.hg18.bed.txt", 
                          package = "methylKit")
gene.obj=readTranscriptFeatures(transcriptBED)

# annotate differentially methylated CpGs with 
# promoter/exon/intron using annotation data
#
annotation <- annotateWithGeneParts(as(results_ranges,"GRanges"),gene.obj)
# Perform the annotation
annotation <- annotateWithGeneParts(as(results_ranges, "GRanges"), gene.obj)


# Save the annotation summary to a text file
sink("annotation_summary.txt")

cat("Summary of AnnotationByGeneParts Object\n")
cat("========================================\n\n")

cat("Distance to TSS:\n")
print(head(annotation@dist.to.TSS))

cat("\nPercentage of Overlaps:\n")
print(annotation@annotation)

cat("\nPercentage of Overlaps with Precedence:\n")
print(annotation@precedence)

cat("\nNumber of Overlapping Features:\n")
print(annotation@num.annotation)

cat("\nNumber of Overlapping Features with Precedence:\n")
print(annotation@num.precedence)

cat("\nPercentage of Annotation Boundaries with Overlap:\n")
print(annotation@perc.of.OlapFeat)

# Close the sink
sink()

cat("Annotation summary saved to 'annotation_summary.txt'.\n")
zip("Project supplementary files.zip", files = 'annotation_summary.txt')

# Split overlapping.genes into separate gene names
split_genes <- unique(unlist(strsplit(results_ranges$overlapping.genes, ", ")))
# Find the number of unique genes in split_genes
length(split_genes)
nrow(filtered_DEGs)

# Ensure the split_genes are of character type for comparison
split_genes <- as.character(split_genes)
split_genes
# Find common genes between the split list and filtered_DEGs$Symbol
common_genes <- intersect(split_genes, filtered_DEGs$Symbol)
length(common_genes)
# Filter the full_results for the common genes
methyl_regulated_genes <- filtered_DEGs[filtered_DEGs$Symbol %in% common_genes, ]
length(methyl_regulated_genes)
head(methyl_regulated_genes)
dim(methyl_regulated_genes)
# Create a DataFrame with the methyl-regulated genes
methyl_regulated_genes_df <- data.frame(methyl_regulated_genes)
nrow(methyl_regulated_genes_df)
# Display the first few rows of the DataFrame
head(methyl_regulated_genes_df)

# Optionally, save the DataFrame to a CSV file
write.csv(methyl_regulated_genes_df, "Methyl_Regulated_Genes.csv", row.names = TRUE)
zip("Project supplementary files.zip", files = "Methyl_Regulated_Genes.csv")
# Define sets for the Venn diagram
set_split_genes <- split_genes
set_filtered_DEGs <- filtered_DEGs$Symbol
# Make sure both sets are vectors
set_split_genes <- as.vector(set_split_genes)
set_filtered_DEGs <- as.vector(set_filtered_DEGs)
# Combine into a list
venn_data <- list(
  "Methylation" = set_split_genes,
  "Expression" = set_filtered_DEGs
)

# Create and save the Venn diagram
venn_plot <- venn.diagram(
  x = venn_data,
  category.names = c("Methylation", "Expression"),
  filename = NULL,  # Capture the plot as a grob
  col = "black",    # Border color of circles
  fill = c("skyblue", "pink"), # Fill colors for the sets
  alpha = c(0.5, 0.5),         # Transparency for the sets
  cex = 1.5,                   # Font size for counts
  cat.cex = 1.5,               # Font size for category labels
  cat.pos = c(-15, 15),        # Adjust label positions
  cat.dist = c(0.05, 0.05),    # Adjust label distance
  main = "Venn Diagram of Methyl Regulated Genes",
  main.cex = 2,                # Title font size
  lty = "solid",               # Line type for circles
  ext.text = FALSE,            # Hide intersection count
  fill.intersect = "purple"    # Add a distinct color to the intersection
)

# Save the plot as a PNG file
png("VennDiagram_IntersectionColor.png", width = 2000, height = 2000, res = 300)
grid.draw(venn_plot)
dev.off()

#pathways

# Example of ID conversion if needed
library(org.Hs.eg.db)
gene_symbols <- methyl_regulated_genes_df$Symbol

# Replace with your gene symbols
entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
gene_list <- na.omit(as.integer(entrez_ids))
print(gene_list)


missing_genes <- gene_symbols[is.na(entrez_ids)]
if (length(missing_genes) > 0) {
  cat("Unmapped gene symbols:", missing_genes, "\n")
}


# Keep only genes with valid Entrez IDs
valid_entrez_ids <- na.omit(entrez_ids)

# Proceed with analysis using `valid_entrez_ids`
print(valid_entrez_ids)
# Check KEGG enrichment result for GNG7
print(kegg_results_df$geneID)  # Look for "2788" or "GNG7" in the gene IDs
# Verify if GNG7 is in your gene list for KEGG analysis
print(valid_entrez_ids %in% "2788")
kegg_results_df[kegg_results_df$geneID == "2788", ]

#KEGG pathway

kegg_results <- enrichKEGG(
  gene = valid_entrez_ids,
  organism = 'hsa',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.5,
  minGSSize = 3,
  maxGSSize = 2000
)
png("KEGG_Pathway_Enrichment_Dotplot.png", width = 3000, height = 5000, res = 300)
dotplot(kegg_results, showCategory = nrow(as.data.frame(kegg_results))) +
  ggtitle("KEGG Pathway Enrichment Analysis")
dev.off()
print(methyl_regulated_genes)
print(kegg_results)
View(kegg_results_df)
# Map Entrez IDs to Gene Symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = valid_entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
# Convert kegg_results to a data frame
kegg_results_df <- as.data.frame(kegg_results)

# Add Gene Symbols to KEGG results (use the data frame for this)
kegg_results_df$GeneSymbols <- sapply(kegg_results_df$geneID, function(gene_list) {
  genes <- unlist(strsplit(gene_list, "/"))
  symbols <- gene_symbols[genes]
  paste(symbols, collapse = "/")
})
print(kegg_results_df)
# Save the KEGG results with gene symbols to CSV
write.csv(as.data.frame(kegg_results_df), "KEGG_Pathway_Enrichment_Results.csv", row.names = TRUE)
zip("Project supplementary files.zip", files = "KEGG_Pathway_Enrichment_Results.csv")

#Gene Ontology
go_results <- enrichGO(
  gene = valid_entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.5
)
go_results


png("Gene_Ontology_Dotplot1.png", width = 3000, height = 5000, res = 300)
dotplot(go_results, showCategory = 4) + ggtitle("GO Biological Processes Enrichment")
dev.off()
go_results_df <- as.data.frame(go_results)
# Add Gene Symbols to KEGG results (use the data frame for this)
go_results_df$GeneSymbols <- sapply(go_results$geneID, function(gene_list) {
  genes <- unlist(strsplit(gene_list, "/"))
  symbols <- gene_symbols[genes]
  paste(symbols, collapse = "/")
})
View(go_results_df)
write.csv(as.data.frame(go_results_df), "GO_enrichment_results.csv", row.names = TRUE)
zip("Project supplementary files.zip", files = "GO_enrichment_results.csv")
