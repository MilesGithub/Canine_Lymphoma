
process_rna_seq_data <- function(rna_seq_data, sample_ids_subset) {
  rna_seq_data %>% 
    select(all_of(sample_ids_subset)) %>% 
    filter(rownames(.) %in% CanFam3_df_sub$gene_id) %>% 
    t() %>% 
    as_tibble(rownames = "sample_id") %>% 
    mutate(across(everything(), as.numeric)) %>% 
    drop_na() %>% 
    mutate(across(everything(), log1p))
}

create_tsne_plot <- function(tsne_data, sample_ids_subset, title, group_condition) {
  tsne_result <- Rtsne(as.matrix(tsne_data), perplexity = 4)
  tsne_df <- tibble(sample_id = factor(sample_ids_subset), 
                    X = tsne_result$Y[, 1], 
                    Y = tsne_result$Y[, 2]) %>% 
    mutate(group = factor(if_else(group_condition(.), 1, 0)))
  
  ggplot(tsne_df, aes(x = X, y = Y, color = group)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("1" = "black", "0" = "#999999")) +
    labs(title = title, x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"),
          axis.text = element_text(size = 12), axis.title = element_text(size = 16),
          plot.title = element_text(size = 18)) +
    guides(color = FALSE)
}

#Figure 1A t-SNE
rna_seq_counts_t1 <- process_rna_seq_data(rna_seq_FPKM, sample_ids_time_point_1)
p1 <- create_tsne_plot(rna_seq_counts_t1, sample_ids_time_point_1, "Day 0", ~ .x$Y < 0)

rna_seq_counts_t2 <- process_rna_seq_data(rna_seq_counts, sample_ids_time_point_2)
p2 <- create_tsne_plot(rna_seq_counts_t2, sample_ids_time_point_2, "Week 6", ~ .x$X > 0)

combined_plot <- ggarrange(p1, p2, ncol = 2)

output_file_tsne <- paste0(project_path,"/images/Figure_1A.pdf")
ggsave(output_file_tsne, combined_plot, width = 5, height = 4)


#Figure 1B KM-plot
clinical_sub <- clinical %>% 
  filter(sample_id %in% sample_ids_time_point_1) %>% 
  mutate(group = factor(tsne_df_t1$group))

fit <- survfit(Surv(progression_free_survival.days., is_progressed) ~ group, data = clinical_sub)

p_survival <- ggsurvplot(fit, palette = c("#B3B3B3", "#333333"), 
                         pval = TRUE, pval.method = TRUE, log.rank.weights = "1", 
                         surv.median.line = 'hv', risk.table = FALSE, 
                         xlab = "Progression Free Survival (days)", 
                         legend = "none", 
                         pval.size = 5, 
                         pval.coord = c(300, 0.7),
                         pval.method.coord = c(300, 0.8),
                         size = 2)

output_file_km <- paste0(project_path,"/images/Figure_1B.pdf")
ggsave(output_file_km, p_survival, width = 5, height = 4)

#Create data.frame with cluster group annotation
cluster_group_df<-data.frame(sample_id=sample_ids,donor=donor,time_point=time_point, group=c(tsne_df_t1$group,tsne_df_t1$group))

