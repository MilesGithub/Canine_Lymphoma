
prepare_clinical_sub <- function(clinical_sub, gene_expression_FPKM) {
  gene_expression_FPKM_t <- t(gene_expression_FPKM)
  gene_expression_FPKM_t_sub <- gene_expression_FPKM_t[rownames(gene_expression_FPKM_t) %in% clinical_sub$sample_id, ]
  gene_expression_FPKM_t_sub <- data.frame(gene_expression_FPKM_t_sub)
  gene_expression_FPKM_t_sub$patient <- rownames(gene_expression_FPKM_t_sub)
  gene_expression_FPKM_t_sub <- gene_expression_FPKM_t_sub[order(gene_expression_FPKM_t_sub$patient), ]
  clinical_sub <- cbind(clinical_sub, gene_expression_FPKM_t_sub)
  return(clinical_sub)
}

compute_differential_expression <- function(df_sub, df_sub_1, df_sub_2) {
  proteasome_subunits$log2fc <- 0
  proteasome_subunits$wilcoxon_pval <- 1
  proteasome_subunits$neg_log <- 0
  for (i in 1:nrow(proteasome_subunits)) {
    print(i)
    row <- proteasome_subunits[i, ]
    ensembl_dog_id <- row$ensembl_dog_id
    
    if (ensembl_dog_id %in% colnames(df_sub)) {
      low <- df_sub_1[, colnames(df_sub_1) == ensembl_dog_id, drop = FALSE]
      high <- df_sub_2[, colnames(df_sub_2) == ensembl_dog_id, drop = FALSE]
      wilcoxon_pval <- wilcox.test(high, low)$p.value
      proteasome_subunits[i, ]$wilcoxon_pval <- wilcoxon_pval
      fc <- foldchange(mean(high), mean(low))
      log2fc <- foldchange2logratio(fc)
      proteasome_subunits[i, ]$log2fc <- log2fc
    }
  }
  
  proteasome_subunits <- proteasome_subunits[order(proteasome_subunits$wilcoxon_pval), ]
  proteasome_subunits <- proteasome_subunits[!duplicated(proteasome_subunits$hgnc_symbol), ]
  proteasome_subunits$wilcoxon_qval <- p.adjust(proteasome_subunits$wilcoxon_pval)
  proteasome_subunits$neg_log <- -log(proteasome_subunits$wilcoxon_qval)
  
  return(proteasome_subunits)
}

plot_proteasome_subunits <- function(proteasome_subunits) {
  p <- ggplot(proteasome_subunits, aes(x = as.numeric(as.character(log2fc)), y = as.numeric(as.character(neg_log)))) +
    geom_point(alpha = 1.2) +
    labs(title = "") +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 1)) +
    geom_hline(yintercept = 3, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    scale_shape_manual(values = c(17, 19)) +
    xlab("Log2 Fold Change") +
    ylab("-log(FDR adjusted p-value)") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14)) +
    geom_text_repel(aes(label = hgnc_symbol), box.padding = 0.5) +
    theme(plot.margin = margin(1, 1, 1, 1))  # Adjusting plot margins
  
  return(p)
}


clinical_sub_gene_expression <- prepare_clinical_sub(clinical_sub, gene_expression_FPKM)

df_sub <- clinical_sub_gene_expression
df_sub_1 <- df_sub[df_sub$group == 0, ]
df_sub_2 <- df_sub[df_sub$group == 1, ]

proteasome_subunits <- compute_differential_expression(df_sub, df_sub_1, df_sub_2)

p1 <- plot_proteasome_subunits(proteasome_subunits)

ggsave(paste0(project_path,"/images/proteasome_subunits.pdf"), plot = p1, width = 5, height = 4)
