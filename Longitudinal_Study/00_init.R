
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
CanFam3_df_sub<-CanFam3_df_sub[CanFam3_df_sub$gene_biotype  == "protein_coding",]

#
sample_ids<-c("003A_S3","024A_S13","025A_S15","026A_S17","027A_S19","028A_S21","029A_S23","030A_S25","031A_S27",
              "033A_S29","034A_S31","035A_S33","036A_S35","037A_S37","038A_S39","003B_S4","024B_S14","025B_S16",
              "026B_S18","027B_S20","028B_S22","029B_S24","030B_S26","032B_S28","033B_S30","034B_S32","035B_S34",
              "036B_S36","037B_S38","038B_S40")

sample_ids_time_point_1 <- c("003A_S3", "024A_S13", "025A_S15", "026A_S17", "027A_S19", "028A_S21", "029A_S23", "030A_S25", "031A_S27", "033A_S29", "034A_S31", "035A_S33", "036A_S35", "037A_S37", "038A_S39")
sample_ids_time_point_2 <- c("003B_S4", "024B_S14", "025B_S16", "026B_S18", "027B_S20", "028B_S22", "029B_S24", "030B_S26", "032B_S28", "033B_S30", "034B_S32", "035B_S34", "036B_S36", "037B_S38", "038B_S40")

time_point <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
donor <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)


