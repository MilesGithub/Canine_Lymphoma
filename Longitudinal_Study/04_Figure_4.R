
#Figure 4B Longitudinal Box Plots
genes <- list(
  c('ENSCAFG00000001835', 'ABCB1', '2.5', '7', '3', '7', 'FDR=0.0006083375'),
  c('ENSCAFG00000009638', 'ABCG2', '3.5', '6', '4', '6', 'FDR=0.0196082114'),
  c('ENSCAFG00000018208', 'ABCC1', '3.5', '8', '4', '8', 'FDR=0.6303075'),
  c('ENSCAFG00000002065', 'KIT', '2.5', '6', '3', '6', 'FDR=3.901365e-03'),
  c('ENSCAFG00000030851', 'HLF', '2.5', '6', '2', '6', 'FDR=1.286682e-04'),
  c('ENSCAFG00000018752', 'IL7R', '2.5', '8', '3', '8', 'FDR=2.121596e-03'),
  c('ENSCAFG00000000482', 'FLT4', '0', '6', '0', '6', 'FDR=1.116118e-02'),
  c('ENSCAFG00000006716', 'FLT3', '3', '6.5', '3', '6', 'FDR=0.008814308'),
  
  c('ENSCAFG00000000280', 'CDK4', '4.5', '9', '4', '9', 'FDR=0.0207515352'),
  c('ENSCAFG00000016714', 'TP53', '4', '9', '4', '9', 'FDR=0.02092605'),
  c('ENSCAFG00000006383', 'BRCA2', '5.8', '10', '6', '10', 'FDR=0.0226009487'),
  c('ENSCAFG00000002659', 'MSH2', '5.8', '10', '6', '10', 'FDR=0.0079303416'),
  c('ENSCAFG00000002664', 'MSH6', '6.5', '10', '7', '10', 'FDR=0.0121722313'),
  c('ENSCAFG00000004052', 'MCM2', '5.5', '9', '6', '9', 'FDR=0.0007213392'),
  c('ENSCAFG00000001709', 'MCM5', '5', '9', '5', '9', 'FDR=0.0046822300'),
  c('ENSCAFG00000014748', 'MCM7', '5.5', '9', '6', '9', 'FDR=0.0048874509')
)

plots <- list()
t <- 0

for (gene in genes) {
  gene_id <- gene[1]
  gene_name <- gene[2]
  min_1 <- as.numeric(gene[3])
  max_1 <- as.numeric(gene[4])
  min_2 <- as.numeric(gene[5])
  max_2 <- as.numeric(gene[6])
  
  rna_seq_counts_sub <- rna_seq_counts[rownames(rna_seq_counts) == gene_id, , drop = FALSE]
  rna_seq_counts_sub_t <- data.frame(t(rna_seq_counts_sub))
  rna_seq_counts_sub_t$sample_id <- rownames(rna_seq_counts_sub_t)
  rna_seq_counts_sub_t <- rna_seq_counts_sub_t[rna_seq_counts_sub_t$sample_id %in% cluster_group_df[cluster_group_df$group == 0, ]$sample_id, ]
  
  rna_seq_counts_sub_t$group <- "Day 0"
  rna_seq_counts_sub_t[rna_seq_counts_sub_t$sample_id %in% cluster_group_df[cluster_group_df$time_point == 2, ]$sample_id, ]$group <- "Week 6"
  
  colnames(rna_seq_counts_sub_t) <- c('counts', 'sample_id', 'group')
  rna_seq_counts_sub_t$counts <- log1p(rna_seq_counts_sub_t$counts)
  rna_seq_counts_sub_t$sample_id <- as.character(rna_seq_counts_sub_t$sample_id)
  rna_seq_counts_sub_t$donor <- rep(1:9, each = 2)
  
  p1 <- ggplot(rna_seq_counts_sub_t) +
    geom_boxplot(aes(x = as.factor(group), y = counts, group = as.factor(group)), size = 0.5) +
    geom_point(aes(x = as.factor(group), y = counts), size = 1.5, color = "#4d4d4d") +
    geom_line(aes(x = as.factor(group), y = counts, group = as.factor(donor)), size = 0.5, color = "#4d4d4d") +
    scale_y_continuous(limits = c(min_1, max_1), breaks = seq(min_2, max_2, 1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"),
          axis.text.y = element_text(size = 12), axis.text.x = element_blank(), axis.title.x = element_text(size = 16),
          plot.title = element_blank()) +
    labs(y = "", x = gene_name)
  
  t <- t + 1
  plots[[t]] <- p1
}

pdf("D:/projects/002_Coomber_lab/canine_lymphoma_time_series/images/Figure_4B.pdf", width = 12, height = 3.5)
grid.arrange(grobs = plots, ncol = 8, nrow = 2)
dev.off()
