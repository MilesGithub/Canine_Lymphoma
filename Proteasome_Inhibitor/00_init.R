
project_path<-""
setwd(project_path)

# Load packages
packages <- c("BiocManager", "dplyr", "ggplot2", "gridExtra", "ggpubr", "survival", "survminer", "DESeq2")

load_packages <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% rownames(installed.packages())) {
      print(paste0("Package already installed: ", pkg))
    } else if (pkg %in% BiocManager::available()) {
      BiocManager::install(pkg)
      print(paste0("Installed from Bioconductor: ", pkg))
    } else {
      install.packages(pkg)
      print(paste0("Installed from CRAN: ", pkg))
    }
  }
  library(pkg, character.only = TRUE)
  print(paste0("Loaded: ", pkg))
}

lapply(packages, load_packages)

# Load data
gene_expression_count<-read.csv('gene_expression_count.csv')

CanFam3_df_sub<-CanFam3_df_sub[CanFam3_df_sub$gene_biotype  == "protein_coding",]

sample_ids<-c("001_S1","006_S6","012_S9","016_S10","009_S8","022_S12","003A_S3","005_S5","026A_S17","027A_S19","024A_S13","025A_S15","033A_S29","035A_S33","037A_S37")
