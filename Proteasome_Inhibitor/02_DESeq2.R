
prepare_colData <- function(countData, sample_ids, clinical_sub) {
  colData <- data.frame(sample_id = colnames(countData), group = 0)
  colData <- colData[colData$sample_id %in% sample_ids, ]
  clinical_sub_sub <- clinical_sub[clinical_sub$group == 1, ]
  colData$group[colData$sample_id %in% clinical_sub_sub$sample_id] <- 1
  colData$group <- as.factor(colData$group)
  rownames(colData) <- colnames(countData)
  return(colData)
}

# DESeq2 analysis
run_deseq2_analysis <- function(countData, colData) {
  countData <- countData[, colnames(countData) %in% rownames(colData)]
  countData <- countData[, rownames(colData)]
  
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ group)
  dds <- estimateSizeFactors(dds)
  idx <- rowSums(counts(dds, normalized = TRUE) >= 5) >= 3
  dds <- dds[idx, ]
  dds <- DESeq(dds)
  return(results(dds, alpha = 0.05))
}

annotate_deseq_results <- function(res, CanFam3_df_sub) {
  res_df <- as.data.frame(res[order(res$padj), ])
  res_df$gene_id <- rownames(res_df)
  CanFam3_df_sub <- CanFam3_df_sub[!duplicated(CanFam3_df_sub$gene_id), ]
  annotated_res <- merge(res_df, CanFam3_df_sub, by = "gene_id", all.x = TRUE)
  return(annotated_res[!is.na(annotated_res$gene_name) & annotated_res$padj < 0.05, ])
}

plot_volcano <- function(dseq_sig_genes) {
  dseq_sig_genes$neg_log <- -log(dseq_sig_genes$padj)
  
  p <- ggplot(dseq_sig_genes, aes(x = as.numeric(log2FoldChange), y = as.numeric(neg_log))) +
    geom_point(aes(color = neg_log > 3), alpha = 1.2, show.legend = FALSE) +
    scale_color_manual(values = c("#aaaaaa", "#333333")) +
    geom_hline(yintercept = 3, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    xlab("Log2 Fold Change") +
    ylab("-log(FDR adjusted p-value)") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14)) +
    theme(plot.margin = margin(1, 1, 1, 1))
  
  p <- p + geom_text(data = subset(dseq_sig_genes, neg_log > 3), 
                     aes(label = gene_name), vjust = -0.5, size = 3)
  return(p)
}

colData <- prepare_colData(gene_expression_count, patients, clinical_sub)
res <- run_deseq2_analysis(gene_expression_count, colData)

dseq_sig_genes <- annotate_deseq_results(res, CanFam3_df_sub)

volcano_plot <- plot_volcano(dseq_sig_genes)

ggsave(paste0(project_path,"/images/DESeq2_volcano.pdf"), plot = volcano_plot, width = 6, height = 5)

