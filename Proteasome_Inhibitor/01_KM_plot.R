
prepare_clinical_data <- function(clinical, sample_ids, gene_expression, gene_column_name) {
  clinical_sub <- clinical[clinical$sample_id %in% sample_ids, ]
  gene_expression_sub <- data.frame(gene_expression[rownames(gene_expression) %in% sample_ids, ])
  colnames(gene_expression_sub) <- gene_column_name
  gene_expression_sub$sample_ids <- rownames(gene_expression_sub)
  gene_expression_sub <- gene_expression_sub[order(gene_expression_sub$patient), ]
  
  clinical_sub <- cbind(clinical_sub, gene_expression_sub)
  
  # Classify patients into high and low expression groups based on median
  clinical_sub$group <- ifelse(clinical_sub[[gene_column_name]] >= median(clinical_sub[[gene_column_name]]), 1, 0)
  
  return(clinical_sub)
}

# Kaplan-Meier plot
plot_km_survival <- function(clinical_sub, survival_time_col, survival_status_col, group_col) {
  formula <- as.formula(paste0("Surv(", survival_time_col, ", ", survival_status_col, ") ~ ", group_col))
  fit <- survfit(formula, data = clinical_sub)
  fit$call$formula <- formula

  p1 <- ggsurvplot(
    fit,
    data = clinical_sub,
    pval = TRUE,
    xlab = "Progression Free Survival (days)",
    ylab = "Survival Probability",
    legend = "none",  # Remove legend
    palette = c("grey", "black"),  # Set colors to grey and black
    pval.size = 6,  # Adjust p-value text size
    pval.coord = c(300, 0.82),  # Set p-value location
    pval.method = TRUE,  # Show p-value method
    pval.method.coord = c(300, 0.9),  # Set method p-value location
    size = 1.5,
    font.x = c(18),
    font.y = c(18),
    font.tickslab = c(14)
  )
  
  return(p1)
}

clinical_sub <- prepare_clinical_data(clinical, sample_ids, gene_expression, "PSMB1")
km_plot <- plot_km_survival(clinical_sub, "progression_free_survival.days.", "is_progressed", "group")

# Save KM plot to PDF
filename<-paste0(project_path,"/images/KM_plot.pdf")
pdf(filename, width = 5, height = 4)
print(km_plot)
dev.off()
