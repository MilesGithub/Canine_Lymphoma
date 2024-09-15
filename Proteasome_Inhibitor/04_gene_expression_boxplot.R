
extract_gene_expression <- function(gene_id, gene_name, gene_expression_FPKM, clinical_sub) {
  sub <- gene_expression_FPKM[rownames(gene_expression_FPKM) == gene_id, ]
  sub <- t(sub)
  sub <- data.frame(sub[rownames(sub) %in% clinical_sub$sample_id, ])
  sub$patient <- rownames(sub)
  sub <- sub[order(sub$patient), ]
  sub <- data.frame(sub[, 1])
  colnames(sub) <- gene_name
  clinical_sub <- cbind(clinical_sub, sub)
  return(clinical_sub)
}

create_boxplot <- function(clinical_sub, gene_name, group_var = "group", color_var = "colour", my_comparisons) {
  p<-ggplot(clinical_sub, aes_string(x = group_var, y = gene_name)) +
    xlab("") +
    ylab("") +
    ggtitle(gene_name) +
    stat_compare_means(comparisons = my_comparisons) +
    geom_boxplot(outlier.colour = NA, outlier.shape = NA, outlier.size = 0) +
    geom_jitter(shape = 16, position = position_jitter(0.2), size = 3, color = clinical_sub[[color_var]]) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"),
          plot.title = element_text(size = 12), axis.text.y = element_text(size = 10))
  return(p)
}

gene_list <- list(
  list(id = 'ENSCAFG00000004116', name = 'PSMB1'),
  list(id = 'ENSCAFG00000003543', name = 'PSMB2'),
  list(id = 'ENSCAFG00000016587', name = 'PSMB3'),
  list(id = 'ENSCAFG00000011295', name = 'PSMB5'),
  list(id = 'ENSCAFG00000009353', name = 'PSMC1'),
  list(id = 'ENSCAFG00000004204', name = 'PSMC2'),
  list(id = 'ENSCAFG00000014484', name = 'ETV4'),
  list(id = 'ENSCAFG00000013420', name = 'ETV5'),
  list(id = 'ENSCAFG00000004658', name = 'CCR4'),
  list(id = 'ENSCAFG00000017146', name = 'CXCR3')
)


for (gene in gene_list) {
  clinical_sub <- extract_gene_expression(gene$id, gene$name, gene_expression_FPKM, clinical_sub)
}

my_comparisons <- list(c("1", "0"))
clinical_sub$colour <- ifelse(clinical_sub$group == 1, 'black', '#bbbbbb')
clinical_sub$group <- as.factor(clinical_sub$group)

gene_names_to_plot1 <- c("PSMB1", "PSMB3", "PSMC1", "PSMC2")
gene_names_to_plot2 <- c("ETV4", "ETV5", "CCR4", "CXCR3")

plot_list1 <- lapply(gene_names_to_plot1, function(gene) create_boxplot(clinical_sub, gene, my_comparisons = my_comparisons))
plot_list2 <- lapply(gene_names_to_plot2, function(gene) create_boxplot(clinical_sub, gene, my_comparisons = my_comparisons))

p1 <- ggarrange(plotlist = plot_list1, ncol = 4, nrow = 1)
ggsave(paste0(project_path,"/images/proteasome_subunit_boxplots.pdf"), plot = p1, width = 10, height = 2)

p2 <- ggarrange(plotlist = plot_list2, ncol = 2, nrow = 2)
ggsave(paste0(project_path,"/images/gene_boxplots.pdf"), plot = p2, width = 5, height = 4)

