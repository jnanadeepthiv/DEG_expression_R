#Assignment 5 


# Downloading libraries required
# -------------------------
#required_pkgs <- c("readxl","tidyverse","pheatmap","RColorBrewer","ggplot2")
#installed <- installed.packages()[, "Package"]
#for(p in required_pkgs){
  #if(! p %in% installed) install.packages(p, repos = "https://cloud.r-project.org")
#}

# Loading the libraries
library(readxl)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)


# Part 1a - Load files 

gene_expr_path <- "Gene_Expression_Data.xlsx"
gene_info_path <- "Gene_Information.csv"
sample_info_path <- "Sample_Information.tsv"

# Read gene expression 
gene_expression <- read_excel(gene_expr_path)
gene_information <- read.csv(gene_info_path, stringsAsFactors = FALSE)
sample_information <- read.delim(sample_info_path, sep = "\t", stringsAsFactors = FALSE)

# Quick checks
cat("Dimensions: gene_expression =", dim(gene_expression), "\n")
cat("Dimensions: gene_information =", dim(gene_information), "\n")
cat("Dimensions: sample_information =", dim(sample_information), "\n")

cat("Column names in sample_information:\n")
colnames(sample_information)

# Optional: view first few rows to see the data
head(sample_information)

# 2) Part 1b - Detecting column names

# Ensure Probe_ID column exists in gene_expression
probe_name <- c("Probe_ID")


# Find sample ID and phenotype columns in sample_information
pheno_cols <- "group"



# Make a small mapping dataframe
sample_map <- sample_information %>%
  dplyr::select(phenotype = !!sym(pheno_cols)) %>%
  mutate(phenotype = tolower(as.character(phenotype))) # lower-case for consistent matching



# View first few phenotypes
cat("\nFirst few phenotypes:\n")
head(sample_map$phenotype)


# Rename sample columns in gene_expression using phenotype
# Get sample columns (everything except Probe_ID)
expr_sample_cols <- setdiff(names(gene_expression), "Probe_ID")

# Get the phenotype for each sample from sample_map
# Assume sample_map$phenotype is in the same order as gene_expression columns
phenotypes <- sample_map$phenotype[1:length(expr_sample_cols)]

# Create new column names: phenotype_1, phenotype_2, ...
new_colnames <- paste0(phenotypes, "_", ave(phenotypes, phenotypes, FUN = seq_along))

# Apply new column names
colnames(gene_expression) <- c("Probe_ID", new_colnames)

# Check first few columns
cat("Renamed sample columns:\n")
print(colnames(gene_expression)[1:min(12, ncol(gene_expression))])


# Part_1c -  Split into Tumor and Normal based on suffix "_tumor" or "_normal" in column names


sample_cols <- setdiff(names(gene_expression), "Probe_ID")

tumor_cols <- sample_cols[grepl("^tumor", sample_cols)]   # starts with "tumor"
normal_cols <- sample_cols[grepl("^normal", sample_cols)] # starts with "normal"

if(length(tumor_cols) < 1 | length(normal_cols) < 1){
  warning("Could not find both tumor and normal sample columns using suffix matching. Inspect column names.")
}

tumor_data <- gene_expression %>% dplyr::select(Probe_ID, all_of(tumor_cols))
normal_data <- gene_expression %>% dplyr::select(Probe_ID, all_of(normal_cols))

head(tumor_data)
head(normal_data)


# part_1d - Compute average expression per probe in each group

# Convert sample columns to numeric (coerce safely)
tumor_matrix <- tumor_data %>% column_to_rownames("Probe_ID") %>% mutate_all(as.numeric)
normal_matrix <- normal_data %>% column_to_rownames("Probe_ID") %>% mutate_all(as.numeric)

tumor_avg <- rowMeans(tumor_matrix, na.rm = TRUE)
normal_avg <- rowMeans(normal_matrix, na.rm = TRUE)

# Combine into a results table
results <- tibble(
  Probe_ID = names(tumor_avg),
  Tumor_Avg = as.numeric(tumor_avg),
  Normal_Avg = as.numeric(normal_avg)
)
head(results)

# Part_1e
# Compute fold changes- Standard safe log2 ratio: log2((Tumor + eps)/(Normal + eps))


results <- results %>%
  mutate(
    Log2_Fold_Change = (Tumor_Avg - Normal_Avg) / (Normal_Avg )
  )
# Show first rows
head(results)
dim(results)

# Part_1f
# Merge with gene_information and identify magnitude > 5

# Ensure gene_information has a Probe_ID column
ginfo_possible <- intersect(c("Probe_ID"), names(gene_information))
if(length(ginfo_possible)==0) {
  stop("gene_information must contain a Probe ID column (e.g. 'Probe_ID').")
}
gene_information <- gene_information %>% rename(Probe_ID = !!sym(ginfo_possible[1]))

merged_results <- results %>% left_join(gene_information, by = "Probe_ID")

# Add absolute fold change
merged_results <- merged_results %>% mutate(Abs_Log2_FC = abs(Log2_Fold_Change))

# Filter significant genes where |log2_fc| > 5
sig_threshold <- 5
significant_genes <- merged_results %>% filter(!is.na(Abs_Log2_FC) & Abs_Log2_FC > sig_threshold)

cat("Number of significant genes with |log2 FC| >", sig_threshold, ":", nrow(significant_genes), "\n")

# Part_1g
# Add column indicating higher expressed in Tumor or Normal

significant_genes <- significant_genes %>%
  mutate(Higher_Expression = ifelse(Log2_Fold_Change > 0, "Tumor", "Normal"))

head(significant_genes)
dim(significant_genes)


# Task 2

# Part2a - Exploratory Data Analysis (EDA)
# Part2b - DEG chromosome wise_histogram





# Ensure Chromosome column exists
if(!"Chromosome" %in% names(significant_genes)){
  significant_genes$Chromosome <- "unknown"
}

# Part 2b - Histogram of DEGs per chromosome
# Count DEGs per chromosome
degs_by_chrom <- table(significant_genes$Chromosome)
degs_df <- as.data.frame(degs_by_chrom)
colnames(degs_df) <- c("Chromosome", "Count")

# Reorder chromosomes by descending DEG count
degs_df$Chromosome <- factor(degs_df$Chromosome,
                             levels = degs_df$Chromosome[order(degs_df$Count, decreasing = TRUE)])

# Plot histogram (bar chart)
p_chrom <- ggplot(degs_df, aes(x = Chromosome, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Number of Significant DEGs per Chromosome",
       x = "Chromosome",
       y = "Number of DEGs") +
  theme_minimal()

print(p_chrom)


# Part 2c - Histogram by sample type (Tumor vs Normal)
# Reorder Chromosome factor in original data for consistency
significant_genes$Chromosome <- factor(significant_genes$Chromosome,
                                       levels = names(sort(table(significant_genes$Chromosome),
                                                           decreasing = TRUE)))

# Plot dodged bar chart
p_chrom_type <- ggplot(significant_genes, aes(x = Chromosome, fill = Higher_Expression)) +
  geom_bar(position = "dodge") +
  labs(title = "Significant DEGs by Chromosome and Expression (Tumor vs Normal)",
       x = "Chromosome",
       y = "Count",
       fill = "Higher Expression") +
  theme_minimal()

print(p_chrom_type)


# Part 2d - Bar chart: percentage of DEGs upregulated (Tumor) vs downregulated (Normal)
up_count <- sum(significant_genes$Higher_Expression == "Tumor", na.rm = TRUE)
down_count <- sum(significant_genes$Higher_Expression == "Normal", na.rm = TRUE)
total_sig <- up_count + down_count

# Create dataframe for percentage plot
regulation_df <- tibble(
  Regulation = c("Upregulated_in_Tumor", "Downregulated_in_Tumor"),
  Count = c(up_count, down_count),
  Percent = c(up_count / total_sig * 100, down_count / total_sig * 100)
)

# Plot percentages
p_percent <- ggplot(regulation_df, aes(x = Regulation, y = Percent, fill = Regulation)) +
  geom_col(color = "black", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Percent)), vjust = -0.3) +
  labs(title = "Percentage of Significant DEGs Up/Down in Tumor",
       y = "Percentage (%)", x = "") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10), expand = expansion(mult = c(0, 0.02))) +
  theme_minimal() +
  theme(legend.position = "none")

print(p_percent)


# Part 2e - Heatmap (raw data from part 1b)
#     Use the original gene_expression numeric matrix but transpose for pheatmap
# -------------------------
# Prepare expression matrix for pheatmap: rows = samples, cols = probes (we will sample for speed if huge)
expr_mat <- gene_expression %>% column_to_rownames("Probe_ID") %>% mutate_all(as.numeric)
# pheatmap prefers samples as rows, probes as cols -> transpose
expr_for_heatmap <- t(expr_mat)

# Extract sample phenotype from column names
sample_phenotype <- tibble(
  sample = rownames(expr_for_heatmap),
  phenotype = ifelse(grepl("^tumor", rownames(expr_for_heatmap)), "Tumor", "Normal")
) %>% column_to_rownames("sample")

# Draw heatmap
pheatmap(expr_for_heatmap,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_row = sample_phenotype,
         main = "Heatmap: Samples (rows) x Probes (cols) - row-scaled")


#  Part 2f - Clustermap 
# Clustermap for all genes
# Use all genes from the expression matrix
expr_all <- t(expr_mat)  # transpose so samples are rows, genes are columns

# Extract sample phenotype
sample_phenotype_all <- tibble(
  sample = rownames(expr_all),
  phenotype = ifelse(grepl("^tumor", rownames(expr_all)), "Tumor", "Normal")
) %>% column_to_rownames("sample")

# Draw clustermap for all genes
pheatmap(expr_all,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_row = sample_phenotype_all,
         main = "Clustermap: All Genes")

#  Part 2g - 
#    Chromosome 11 contains the highest number of differentially expressed genes (DEGs), 
#    while most other chromosomes have comparatively fewer. DEGs are notably concentrated on chromosomes 
#    1, 2, 3, 4, 5, and X.

#    The dataset contains 45 probes exhibiting substantial fold changes (|FC| > 5), 
#    all of which are upregulated in Tumor samples. This explains why 100% of the DEGs are classified as 
#    upregulated in tumors.

#    All identified DEGs show tumor-specific overexpression, confirming a 
#    consistent upregulation pattern across Tumor samples in the 1g dataset.

#    The heatmap and clustermap visualizations reveal groups of genes with similar 
#    expression trends, highlighting co-expression among DEGs.

#     The clustermap shows distinct clusters of samples with similar gene expression profiles, 
#    emphasizing the differences between Tumor and Normal samples.
