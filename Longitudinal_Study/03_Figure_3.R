
#Figure 3B Cluster Group 1 vs Cluster Group 2 Gene Expression Box Plots
genes <- list(
  list(gene_id = 'ENSCAFG00000009669', gene_name = 'MDM4', fdr = 'FDR=2.471411e-37'),
  list(gene_id = 'ENSCAFG00000000418', gene_name = 'MDM2', fdr = 'FDR=3.966593e-84'),
  
  list(gene_id = 'ENSCAFG00000029853', gene_name = 'CDKN1A', fdr = 'FDR=3.010075e-83'),
  list(gene_id = 'ENSCAFG00000004436', gene_name = 'RB1', fdr = 'FDR=1.011493e-70'),
  
  list(gene_id = 'ENSCAFG00000009532', gene_name = 'NRAS', fdr = 'FDR=3.893408e-14'),
  list(gene_id = 'ENSCAFG00000003907', gene_name = 'BRAF', fdr = 'FDR=5.510416e-43'),
  list(gene_id = 'ENSCAFG00000019138', gene_name = 'MAP2K2', fdr = 'FDR=3.527271e-209'),
  list(gene_id = 'ENSCAFG00000015421', gene_name = 'MAPK1', fdr = 'FDR=1.107599e-15'),
  list(gene_id = 'ENSCAFG00000001086', gene_name = 'MYC', fdr = 'FDR=1.577550e-41'),
  
  list(gene_id = 'ENSCAFG00000011212', gene_name = 'PIK3CA', fdr = 'FDR=4.867437e-103'),
  list(gene_id = 'ENSCAFG00000015670', gene_name = 'PTEN', fdr = 'FDR=1.781399e-116'),
  list(gene_id = 'ENSCAFG00000005388', gene_name = 'AKT2', fdr = 'FDR=3.887699e-42'),
  list(gene_id = 'ENSCAFG00000016648', gene_name = 'MTOR', fdr = 'FDR=2.368915e-38'),
  list(gene_id = 'ENSCAFG00000005965', gene_name = 'FOXO1', fdr = 'FDR=1.039396e-59'),
  
  list(gene_id = 'ENSCAFG00000005449', gene_name = 'TGFBR2', fdr = 'FDR=2.658831e-21'),
  list(gene_id = 'ENSCAFG00000017567', gene_name = 'SMAD2', fdr = 'FDR=1.595678e-13'),
  list(gene_id = 'ENSCAFG00000017388', gene_name = 'SMAD3', fdr = 'FDR=7.895353e-55'),
  list(gene_id = 'ENSCAFG00000000164', gene_name = 'SMAD4', fdr = 'FDR=2.441427e-61'),
  
  list(gene_id = 'ENSCAFG00000010730', gene_name = 'NFKB1', fdr = 'FDR=3.421835e-06'),
  list(gene_id = 'ENSCAFG00000010155', gene_name = 'NFKB2', fdr = 'FDR=1.781399e-116'),
  list(gene_id = 'ENSCAFG00000013334', gene_name = 'RELA', fdr = 'FDR=4.933817e-67'),
  list(gene_id = 'ENSCAFG00000004597', gene_name = 'RELB', fdr = 'FDR=6.627144e-65'),
  list(gene_id = 'ENSCAFG00000005526', gene_name = 'IKBKB', fdr = 'FDR=1.046494e-120'),
  list(gene_id = 'ENSCAFG00000013418', gene_name = 'NFKBIA', fdr = 'FDR=1.162952e-71'),
  
  list(gene_id = 'ENSCAFG00000020355', gene_name = 'CTCF', fdr = 'FDR=1.058449e-116'),
  list(gene_id = 'ENSCAFG00000001835', gene_name = 'ABCB1', fdr = 'FDR=5.489665e-05')
)

# Initialize plots list
plots <- list()

# Function to create plots
create_plot <- function(gene, rna_seq_counts, sample_ids_time_point_2, cluster_group_df) {
  gene_id <- gene$gene_id
  gene_name <- gene$gene_name
  title <- gene$fdr
  
  # Subset the RNA sequence counts for the specific gene
  rna_seq_counts_sub <- rna_seq_counts[rownames(rna_seq_counts) == gene_id, , drop = FALSE]
  rna_seq_counts_sub_t <- data.frame(t(rna_seq_counts_sub))
  rna_seq_counts_sub_t$sample_id <- rownames(rna_seq_counts_sub_t)
  
  # Subset for time point 2
  rna_seq_counts_sub_t <- rna_seq_counts_sub_t[rna_seq_counts_sub_t$sample_id %in% sample_ids_time_point_2, ]
  
  # Assign group labels based on sample ID and group data
  rna_seq_counts_sub_t$group <- "Group 1"
  group_2_samples <- cluster_group_df[cluster_group_df$group == 1 & cluster_group_df$time_point == 2, ]$sample_id
  rna_seq_counts_sub_t[rna_seq_counts_sub_t$sample_id %in% group_2_samples, ]$group <- "Group 2"
  
  # Adjust column names
  colnames(rna_seq_counts_sub_t) <- c('counts', 'sample_id', 'group')
  
  # Log-transform counts
  rna_seq_counts_sub_t$counts <- log1p(rna_seq_counts_sub_t$counts)
  
  # Create boxplot with jitter
  p <- ggplot(rna_seq_counts_sub_t) +
    geom_boxplot(aes(x = as.factor(group), y = counts, group = as.factor(group)), size = 1, outlier.shape = NA) +
    geom_jitter(aes(x = as.factor(group), y = counts), position = position_jitter(0.2), size = 3) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"),
          axis.text.x = element_blank(), axis.text.y = element_text(size = 19), 
          axis.title.y = element_text(size = 30), axis.title.x = element_text(size = 29),
          plot.title = element_text(size = 22)) +
    labs(y = "", x = gene_name, title = title)
  
  return(p)
}

# Loop through genes and create plots
for (i in seq_along(genes)) {
  plots[[i]] <- create_plot(genes[[i]], rna_seq_counts, sample_ids_time_point_2, cluster_group_df)
}

# Save all plots to a PDF
output_path <- "D:/projects/002_Coomber_lab/canine_lymphoma_time_series/images/Figure_3B.pdf"
pdf(output_path, width = 26.5, height = 7)
grid.arrange(grobs = plots, ncol = 9, nrow = 3)
dev.off()
