
generate_colData <- function(group_df, group_condition, time_condition = NULL) {
  if (is.null(time_condition)) {
    data <- group_df[group_df$group == group_condition, ]
  } else {
    data <- group_df[group_df$time_point == time_condition, ]
  }
  colData <- data.frame(data$sample_id, as.factor(data$time_point), as.factor(data$donor))
  colnames(colData) <- c('sample_id', 'time_point', 'donor')
  rownames(colData) <- colData$sample_id
  return(colData)
}

run_deseq_analysis <- function(colData, rna_seq_counts, CanFam3_df, design_formula) {
  countData <- rna_seq_counts[, colnames(rna_seq_counts) %in% rownames(colData)]
  countData <- countData[, rownames(colData)]
  countData <- countData[rownames(countData) %in% CanFam3_df$gene_id, ]
  all(rownames(colData) %in% colnames(countData)) # Check rownames
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = design_formula)
  dds <- estimateSizeFactors(dds)
  idx <- rowSums(counts(dds, normalized = TRUE) >= 5) >= 3
  dds <- dds[idx, ]
  dds <- DESeq(dds)
  results <- data.frame(results(dds, alpha = 0.05))
  results$ensembl_id <- rownames(results)
  return(results)
}

merge_and_filter_results <- function(deseq_results, CanFam3_df, cgc_genes) {
  CanFam3_df_sub <- CanFam3_df[CanFam3_df$gene_id %in% deseq_results$ensembl_id, ]
  CanFam3_df_sub <- CanFam3_df_sub[!duplicated(CanFam3_df_sub$gene_id) & !is.na(CanFam3_df_sub$gene_name), ]
  deseq_results <- deseq_results[deseq_results$ensembl_id %in% CanFam3_df_sub$gene_id, ]

  deseq_results <- deseq_results[order(deseq_results$ensembl_id), ]
  CanFam3_df_sub <- CanFam3_df_sub[order(CanFam3_df_sub$gene_id), ]
  deseq_results$gene_name <- CanFam3_df_sub$gene_name
  deseq_results$log_padj <- -log(deseq_results$padj)

  deseq_results <- deseq_results[deseq_results$baseMean > 50, ]
  
  deseq_results <- deseq_results$gene_name %in% cgc_genes
  
  deseq_results <- deseq_results[deseq_results$padj < 0.05 & !is.na(deseq_results$gene_name), ]
  
  return(deseq_results)
}

# Group 1 - timepoint 1 vs timepoint 2 (g1_t1_vs_g1_t2)
colData_g1 <- generate_colData(group_df, group_condition = 1)
deseq_results_g1_t1_vs_g1_t2 <- run_deseq_analysis(colData_g1, rna_seq_counts, CanFam3_df, design_formula = ~ donor + time_point)
write.table(deseq_results_g1_t1_vs_t2, 'deseq_results_g1_t1_vs_g1_t2.csv', sep=",")
g1_t1_vs_g1_t2_filtered <- merge_and_filter_results(deseq_results_g1_t1_vs_g1_t2, CanFam3_df, cgc_genes)

# Group 2 - timepoint 1 vs timepoint 2 (g2_t1_vs_g2_t2)
colData_g2 <- generate_colData(group_df, group_condition = 0)
deseq_results_g2_t1_vs_g2_t2 <- run_deseq_analysis(colData_g2, rna_seq_counts, CanFam3_df, design_formula = ~ donor + time_point)
write.table(g2_t1_vs_g2_t2_filtered, 'g2_t1_vs_g2_t2_filtered.csv', sep=",")
g2_t1_vs_g2_t2_filtered <- merge_and_filter_results(deseq_results_g2_t1_vs_g2_t2, CanFam3_df, cgc_genes)

# Group 1 vs Group 2 - timepoint 1 (g1_t1_vs_g2_t1)
colData_g1_g2_t1 <- generate_colData(group_df, group_condition = NULL, time_condition = 1)
deseq_results_g1_t1_vs_g2_t1 <- run_deseq_analysis(colData_g1_g2_t1, rna_seq_counts, CanFam3_df, design_formula = ~ group)
write.table(g1_t2_vs_g2_t2_filtered, 'g1_t2_vs_g2_t2_filtered.csv', sep=",")
g1_t1_vs_g2_t1_filtered <- merge_and_filter_results(deseq_results_g1_t1_vs_g2_t1, CanFam3_df, cgc_genes)

# Group 1 vs Group 2 - timepoint 2 (g1_t2_vs_g2_t2)
colData_g1_g2_t2 <- generate_colData(group_df, group_condition = NULL, time_condition = 2)
deseq_results_g1_t2_vs_g2_t2 <- run_deseq_analysis(colData_g1_g2_t2, rna_seq_counts, CanFam3_df, design_formula = ~ group)
write.table(g1_t2_vs_g2_t2_filtered, 'g1_t2_vs_g2_t2_filtered.csv', sep=",")
g1_t2_vs_g2_t2_filtered <- merge_and_filter_results(deseq_results_g1_t2_vs_g2_t2, CanFam3_df, cgc_genes)

# Write DESeq results tables
write.table(g1_t2_vs_g2_t2_filtered, 'g1_t2_vs_g2_t2_filtered.csv', sep=",")
write.table(g1_t1_vs_g2_t1_filtered, 'g1_t1_vs_g2_t1_filtered.csv', sep=",")
write.table(g2_t1_vs_g2_t2_filtered, 'g2_t1_vs_g2_t2_filtered.csv',sep=",")
write.table(g1_t1_vs_g1_t2_filtered,'g1_t1_vs_g1_t2_filtered.csv',sep=",")


#Figure 2A DEG Volcano Plots
create_volcano_plot <- function(data, title, y_limits, y_breaks, x_limits, x_breaks) {
  ggplot(data, aes(x = as.numeric(as.character(log2FoldChange)), 
                   y = as.numeric(as.character(log_padj)))) +
    geom_point(aes(color = ifelse(log_padj > 3, "#4d4d4d", "#b3b3b3")), alpha = 1.2) +
    scale_color_identity() +
    labs(title = title) +
    geom_hline(yintercept = 3, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    scale_y_continuous(limits = y_limits, breaks = seq(y_limits[1], y_limits[2], y_breaks)) +
    scale_x_continuous(limits = x_limits, breaks = seq(x_limits[1], x_limits[2], x_breaks)) +
    xlab("Log2 Fold Change") +
    ylab("-log(FDR adjusted p-value)") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16))
}

p1 <- create_volcano_plot(deseq_results_g1_t1_vs_g1_t2, "g1_t1_vs_g1_t2", c(0, 15), 5, c(-5, 5), 5)
p2 <- create_volcano_plot(deseq_results_g2_t1_vs_g2_t2, "g2_t1_vs_g2_t2", c(0, 15), 5, c(-5, 5), 5)
p3 <- create_volcano_plot(deseq_results_g1_t1_vs_g2_t1, "g1_t1_vs_g2_t1", c(0, 800), 100, c(-15, 15), 5)
p4 <- create_volcano_plot(deseq_results_g1_t2_vs_g2_t2, "g1_t2_vs_g2_t2", c(0, 800), 100, c(-15, 15), 5)

# Save Figure PDF
combined_plot <- ggarrange(plotlist = list(p1,p2,p3,p4), ncol = 4)
output_file_name <- paste0(project_path,"/images/Figure_2A.pdf")
ggsave(output_file_name, combined_plot, width = 5, height = 4)


#Figure 2B DEG Bar Plots
create_DEG_df <- function(comparison_labels, deg_counts) {
  data.frame(comparison = factor(comparison_labels, levels = comparison_labels), 
             num_DEGs = deg_counts)
}

create_DEG_bar_plot <- function(deg_df, y_label = "Differentially Expressed Genes") {
  ggplot(deg_df, aes(x = comparison, y = num_DEGs)) +
    geom_bar(stat = "identity") +
    labs(x = "Category", y = y_label) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(color = "#2d3f44"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Increase x-axis text size
          axis.text.y = element_text(size = 14),  # Increase y-axis text size
          axis.title.y = element_text(size = 12),  # Increase y-axis title size
          axis.title.x = element_blank())
}

g1_t1_vs_g2_t1_filtered_up <- g1_t1_vs_g2_t1_filtered[g1_t1_vs_g2_t1_filtered$log2FoldChange > 0, ]
g1_t1_vs_g2_t1_filtered_down <- g1_t1_vs_g2_t1_filtered[g1_t1_vs_g2_t1_filtered$log2FoldChange < 0, ]

g1_t2_vs_g2_t2_filtered_up <- g1_t2_vs_g2_t2_filtered[g1_t2_vs_g2_t2_filtered$log2FoldChange > 0, ]
g1_t2_vs_g2_t2_filtered_down <- g1_t2_vs_g2_t2_filtered[g1_t2_vs_g2_t2_filtered$log2FoldChange < 0, ]

g1_t1_vs_g1_t2_filtered_up <- g1_t1_vs_g1_t2_filtered[g1_t1_vs_g1_t2_filtered$log2FoldChange > 0, ]
g1_t1_vs_g1_t2_filtered_down <- g1_t1_vs_g1_t2_filtered[g1_t1_vs_g1_t2_filtered$log2FoldChange < 0, ]

g2_t1_vs_g2_t2_filtered_up <- g2_t1_vs_g2_t2_filtered[g2_t1_vs_g2_t2_filtered$log2FoldChange > 0, ]
g2_t1_vs_g2_t2_filtered_down <- g2_t1_vs_g2_t2_filtered[g2_t1_vs_g2_t2_filtered$log2FoldChange < 0, ]


comparison_1 <- c('g1_t1_vs_g2_t1_up', 'g1_t1_vs_g2_t1_down', 'g1_t2_vs_g2_t2_up', 'g1_t2_vs_g2_t2_down')
DEG_counts_1 <- c(nrow(g1_t1_vs_g2_t1_filtered_up), nrow(g1_t1_vs_g2_t1_filtered_down),
                  nrow(g1_t2_vs_g2_t2_filtered_up), nrow(g1_t2_vs_g2_t2_filtered_down))

DEG_df_1 <- create_DEG_df(comparison_1, DEG_counts_1)
p1 <- create_DEG_bar_plot(DEG_df_1)

combined_1 <- g1_t1_vs_g2_t1_filtered_up[g1_t1_vs_g2_t1_filtered_up$gene_name %in% g1_t2_vs_g2_t2_filtered_down$gene_name, ]
combined_2 <- g1_t1_vs_g2_t1_filtered_up[g1_t1_vs_g2_t1_filtered_up$gene_name %in% g1_t2_vs_g2_t2_filtered_down$gene_name, ]

week_6_only_up <- g1_t2_vs_g2_t2_filtered_up[!(g1_t2_vs_g2_t2_filtered_up$gene_name %in% g1_t1_vs_g2_t1_results_filtered_up$gene_name), ]
week_6_only_down <- g1_t2_vs_g2_t2_filtered_down[!(g1_t2_vs_g2_t2_filtered_down$gene_name %in% g1_t1_vs_g2_t1_filtered_down$gene_name), ]

# Write results tables
write.table(week_6_only_up, 'week_6_only_up.csv', sep = ",")
write.table(week_6_only_down, 'week_6_only_down.csv', sep = ",")
write.table(combined_1, 'combined_up_same_FC.csv', sep = ",")
write.table(combined_2, 'combined_down_same_FC.csv', sep = ",")

combined_up_same_FC <- 3874
combined_down_same_FC <- 3236
combined_diff <- 0
week_6_only_up <- 203
week_6_only_down <- 102

comparison_2 <- c('combined_up_same_FC', 'combined_down_same_FC', 'combined_different_FC', 'week_6_only_up', 'week_6_only_down')
DEG_counts_2 <- c(combined_up_same_FC, combined_down_same_FC, combined_diff, week_6_only_up, week_6_only_down)

DEG_df_2 <- create_DEG_df(comparison_2, DEG_counts_2)
p2 <- create_bar_plot(DEG_df_2)


comparison_3 <- c('g1_t1_vs_g1_t2_up', 'g1_t1_vs_g1_t2_down', 'g2_t1_vs_g2_t2_up', 'g2_t1_vs_g2_t2_down')
DEG_counts_3 <- c(nrow(g1_t1_vs_g1_t2_filtered_up), nrow(g1_t1_vs_g1_t2_filtered_down),
                  nrow(g2_t1_vs_g2_t2_filtered_up), nrow(g2_t1_vs_g2_t2_filtered_down))

DEG_df_3 <- create_DEG_df(comparison_3, DEG_counts_3)
p3 <- create_bar_plot(DEG_df_3)

# Save Figure PDF
combined_plot <- ggarrange(plotlist = list(p1,p2,p3), ncol = 3)
output_file_name <- paste0(project_path,"/images/Figure_2B.pdf")
ggsave(output_file_name, combined_plot, width = 5, height = 4)


#Figure 2C Enriched Term Bar Plots
create_enrichment_bar_plot <- function(comparison, num_enriched_terms, y_limit, y_breaks) {
  
  enrichment_BP <- data.frame(comparison = factor(comparison, levels = comparison), num_enriched_terms = num_enriched_terms)
  
  ggplot(enrichment_BP, aes(x = comparison, y = num_enriched_terms)) +
    geom_bar(stat = "identity") +
    labs(x = "Category", y = "Enriched GO:BP terms") +
    scale_y_continuous(limits = c(0, y_limit), breaks = seq(0, y_limit, y_breaks)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(color = "#2d3f44"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 12),
          axis.title.x = element_blank())
}


comparison_data <- list(
  list(comparison = c('GO:MF','GO:BP','GO:CC','TRANSFAC'), num_enriched_terms = c(23, 265, 59, 110), y_limit = 1200, y_breaks = 200),
  list(comparison = c('GO:MF','GO:BP','GO:CC','TRANSFAC'), num_enriched_terms = c(156, 740, 223, 1091), y_limit = 1200, y_breaks = 200),
  list(comparison = c('GO:MF','GO:BP','GO:CC','TRANSFAC'), num_enriched_terms = c(22, 279, 54, 6), y_limit = 300, y_breaks = 50),
  list(comparison = c('GO:MF','GO:BP','GO:CC','TRANSFAC'), num_enriched_terms = c(56, 225, 91, 282), y_limit = 300, y_breaks = 50)
)

plots <- lapply(comparison_data, function(data) {
  create_enrichment_bar_plot(data$comparison, data$num_enriched_terms, data$y_limit, data$y_breaks)
})

# Save Figure PDF
combined_plot <- ggarrange(plotlist = plots, ncol = 4)
output_file_tsne <- paste0(project_path,"/images/Figure_2C.pdf")
ggsave(output_file_tsne, combined_plot, width = 5, height = 4)



