# Differential Gene Expression Analysis and Visualization

## Author Information
**Programmer:** Jnana Deepthi Vishnumolakala  
**Language:** R (4.4.1) 
**Date Submitted:** November 23rd , 2025  
**Purpose:**  
This script performs a comprehensive analysis of gene expression data, including exploratory data analysis (EDA) and visualization. The primary objectives are to identify differentially expressed genes (DEGs) between Tumor and Normal samples, analyze their distribution across chromosomes, and generate heatmaps and clustermaps to explore co-expression patterns.

---

## Required Files
The following input files are required to run the analysis:

| File Name                     | Description                                                                 |
|--------------------------------|-----------------------------------------------------------------------------|
| `Gene_Expression_Data.xlsx`    | Contains raw expression values for all probes across all samples.           |
| `Gene_Information.csv`         | Contains metadata about each gene, such as chromosome location and annotation. |
| `Sample_Information.tsv`       | Contains sample phenotypes (Tumor or Normal) corresponding to columns in the expression data. |

---

## Required Libraries / Software
To run this analysis, the following R packages are required:

- `readxl` – for reading Excel files  
- `tidyverse` – for data manipulation and plotting (`dplyr`, `ggplot2`, `tibble`)  
- `pheatmap` – for heatmaps and clustermaps  
- `RColorBrewer` – for color palettes  

**Optional:** Installing packages via `install.packages()` if not already available.

```r
# Install required packages (if not already installed)
required_pkgs <- c("readxl","tidyverse","pheatmap","RColorBrewer")
installed <- installed.packages()[, "Package"]
for(p in required_pkgs){
  if(! p %in% installed) install.packages(p)
}

# Load libraries
library(readxl)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

---

## **Required Software / Packages**
- **R** (Base installation; no external packages needed)  
- **RStudio** (optional but recommended)

---


## Overview of the Analysis

The script performs a comprehensive analysis of gene expression data, including differential expression and visualization. Key steps include:

1. **Data Loading:**  
   Loads gene expression data, gene annotation, and sample metadata.

2. **Sample Renaming:**  
   Renames sample columns based on phenotype (Tumor/Normal) and appends unique suffixes to ensure uniqueness.

3. **Data Splitting & Aggregation:**  
   Splits the dataset into Tumor and Normal groups and computes the average expression for each probe in both groups.

4. **Differential Expression Analysis:**  
   Calculates fold changes between Tumor and Normal samples and identifies differentially expressed genes (DEGs) with an absolute fold change > 5.  
   Each DEG is annotated as **higher expressed in Tumor** or **higher expressed in Normal**.

5. **Visualization:**  
   Generates informative plots for exploratory data analysis:
   - **Histogram** of DEGs by chromosome  
   - **Histogram** of DEGs by chromosome segregated by sample type (Tumor vs Normal)  
   - **Bar chart** showing the percentage of DEGs upregulated in Tumor vs downregulated  
   - **Heatmap** of the top 50 probes across all samples  
   - **Clustermap** highlighting co-expression patterns among top DEGs
