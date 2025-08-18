# Helper Functions for TermineR Analysis Pipeline
# These functions reduce code redundancy across analysis scripts

#' Create standardized volcano plot
#' @param data Data frame with differential analysis results
#' @param contrast_name Name of the contrast to plot
#' @param fc_threshold Fold change threshold
#' @param pval_threshold P-value threshold
#' @return ggplot object
create_volcano_plot <- function(data, contrast_name, fc_threshold = 2.5, pval_threshold = 0.05) {
  data %>%
    filter(contrast == contrast_name) %>%
    mutate(
      regulation = case_when(
        adj.P.Val < pval_threshold & logFC > log2(fc_threshold) ~ "upregulated",
        adj.P.Val < pval_threshold & logFC < -log2(fc_threshold) ~ "downregulated",
        TRUE ~ "not_significant"
      ),
      neg_log10_pval = -log10(adj.P.Val)
    ) %>%
    ggplot(aes(x = logFC, y = neg_log10_pval, color = regulation)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c(
      "upregulated" = "red",
      "downregulated" = "blue", 
      "not_significant" = "gray"
    )) +
    geom_vline(xintercept = c(-log2(fc_threshold), log2(fc_threshold)), 
               linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(pval_threshold), 
               linetype = "dashed", alpha = 0.5) +
    labs(
      title = paste("Volcano Plot -", contrast_name),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_minimal()
}

#' Perform PCA analysis with standardized output
#' @param se_object SummarizedExperiment object
#' @param experimental_design Data frame with sample annotation
#' @param title_suffix Suffix for plot title
#' @return List with plot, data, and variance explained
perform_pca_analysis <- function(se_object, experimental_design, title_suffix = "") {
  pca_matrix <- assay(se_object) %>% t() %>% na.omit()
  
  pca_result <- prcomp(pca_matrix, scale. = TRUE)
  
  pca_df <- data.frame(
    sample = rownames(pca_result$x),
    PC1 = pca_result$x[,1],
    PC2 = pca_result$x[,2]
  ) %>%
    left_join(experimental_design, by = "sample")
  
  var_exp <- summary(pca_result)$importance[2,]
  
  pca_plot <- pca_df %>%
    ggplot(aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.68) +
    labs(
      title = paste("PCA Analysis", title_suffix),
      x = paste("PC1 -", round(var_exp[1] * 100, 1), "% variance"),
      y = paste("PC2 -", round(var_exp[2] * 100, 1), "% variance"),
      color = "Condition"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(list(plot = pca_plot, data = pca_df, variance = var_exp))
}

#' Create identification statistics plot
#' @param data Annotated peptide data
#' @param experimental_design Sample annotation
#' @param instrument Instrument prefix
#' @return List with plot and summary data
create_identification_plots <- function(data, experimental_design, instrument) {
  ids_per_sample <- data %>%
    select(starts_with(instrument)) %>%
    pivot_longer(cols = starts_with(instrument), names_to = "sample", values_to = "value") %>%
    group_by(sample) %>%
    summarise(n_identifications = sum(!is.na(value)), .groups = 'drop') %>%
    left_join(experimental_design, by = "sample")
  
  p1 <- ids_per_sample %>%
    ggplot(aes(x = condition, y = n_identifications, fill = condition)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.7) +
    labs(title = "Number of Identifications by Condition",
         x = "Condition", y = "Number of Identifications") +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(list(plot = p1, data = ids_per_sample))
}

#' Generate limma topTables for multiple contrasts
#' @param fit_object limma fit object
#' @param contrasts_vector Named vector of contrasts
#' @param row_data Data frame with feature annotations
#' @return List of topTables
generate_toptables <- function(fit_object, contrasts_vector, row_data) {
  get_top_table <- function(x) {
    topTable(fit_object, coef = x, number = Inf, adjust.method = "BH") %>%
      rownames_to_column(var = "nterm_modif_peptide") %>%
      left_join(row_data, by = "nterm_modif_peptide")
  }
  
  topTables <- map(contrasts_vector, get_top_table)
  names(topTables) <- names(contrasts_vector)
  
  return(topTables)
}

#' Safe RDS file reading with existence check
#' @param file_path Path to RDS file
#' @param computation_function Function to run if file doesn't exist
#' @param message_text Custom message to display
#' @return Object from RDS file or computation
safe_rds_read <- function(file_path, computation_function = NULL, message_text = NULL) {
  if (file.exists(file_path)) {
    if (!is.null(message_text)) {
      message(paste("Loading cached data:", message_text))
    }
    return(readRDS(file_path))
  } else {
    if (!is.null(computation_function)) {
      if (!is.null(message_text)) {
        message(paste("Computing:", message_text))
      }
      result <- computation_function()
      saveRDS(result, file_path)
      return(result)
    } else {
      stop(paste("File does not exist and no computation function provided:", file_path))
    }
  }
}

#' Create regulation heatmap
#' @param results_data Data frame with differential analysis results
#' @param top_n Number of top features to include
#' @param fc_threshold Fold change threshold
#' @param pval_threshold P-value threshold
#' @return pheatmap object or NULL
create_regulation_heatmap <- function(results_data, top_n = 50, fc_threshold = 2.5, pval_threshold = 0.05) {
  top_features <- results_data %>%
    filter(adj.P.Val < pval_threshold, abs(logFC) > log2(fc_threshold)) %>%
    arrange(adj.P.Val) %>%
    slice_head(n = top_n)
  
  if (nrow(top_features) == 0) {
    message("No significantly regulated features found for heatmap.")
    return(NULL)
  }
  
  heatmap_matrix <- top_features %>%
    select(nterm_modif_peptide, contrast, logFC) %>%
    pivot_wider(names_from = contrast, values_from = logFC, values_fill = 0) %>%
    column_to_rownames("nterm_modif_peptide") %>%
    as.matrix()
  
  pheatmap(
    heatmap_matrix,
    scale = "none",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation", 
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = paste("Top", min(top_n, nrow(top_features)), "Regulated Features"),
    fontsize_row = 8,
    fontsize_col = 10
  )
}
