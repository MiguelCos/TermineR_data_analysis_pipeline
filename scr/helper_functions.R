# Helper Functions for TermineR Analysis Pipeline
# These functions reduce code redundancy across analysis scripts

is_pipeline_absolute_path <- function(path_value) {
  is.character(path_value) &&
    length(path_value) == 1 &&
    !is.na(path_value) &&
    grepl("^(?:[A-Za-z]:[/\\\\]|[/\\\\]{2}|~[/\\\\])", path_value)
}

resolve_pipeline_path <- function(path_value) {
  if (is.null(path_value) || !is.character(path_value) || length(path_value) != 1 || is.na(path_value) || path_value == "") {
    return(path_value)
  }

  if (is_pipeline_absolute_path(path_value)) {
    return(path_value)
  }

  here::here(path_value)
}

resolve_pipeline_config_paths <- function(config) {
  path_keys <- c(
    "diann_report_location",
    "fragpipe_parent_dir",
    "fragpipe_lf_parent_dir",
    "fragpipe_lf_annotation",
    "fragpipe_hl_file",
    "fragpipe_hl_annotation",
    "spectronaut_report",
    "location_annotation",
    "fasta_location",
    "uniprot_annotation",
    "protein_annotation_path",
    "rds_dir",
    "results_dir"
  )

  for (key in intersect(names(config), path_keys)) {
    config[[key]] <- resolve_pipeline_path(config[[key]])
  }

  config
}

simplify_pipeline_config <- function(x) {
  if (!is.list(x)) {
    return(x)
  }

  if (length(x) == 0) {
    return(x)
  }

  simplified <- lapply(x, simplify_pipeline_config)

  if (all(vapply(simplified, function(item) !is.list(item) && length(item) == 1, logical(1)))) {
    return(unlist(simplified, recursive = FALSE, use.names = TRUE))
  }

  simplified
}

load_pipeline_config <- function(stage, config_path = here::here("config.yaml")) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read config.yaml. Install it with install.packages('yaml').")
  }

  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path)
  }

  raw_config <- yaml::read_yaml(config_path)
  common_config <- raw_config$common %||% list()
  stage_configs <- raw_config$stages %||% list()

  if (!stage %in% names(stage_configs)) {
    stop(
      "Unknown pipeline stage: ", stage,
      ". Expected one of: ", paste(names(stage_configs), collapse = ", ")
    )
  }

  stage_config <- stage_configs[[stage]] %||% list()
  merged_config <- utils::modifyList(common_config, stage_config, keep.null = TRUE)
  merged_config <- simplify_pipeline_config(merged_config)
  merged_config <- resolve_pipeline_config_paths(merged_config)
  merged_config$stage <- stage

  merged_config
}

load_stage_parameters <- function(stage, config_path = here::here("config.yaml"), envir = parent.frame()) {
  config <- load_pipeline_config(stage = stage, config_path = config_path)
  list2env(config, envir = envir)
  invisible(config)
}

default_pipeline_config_path <- function() {
  candidate_paths <- c(
    here::here("r_project_files", "config.yaml"),
    here::here("config.yaml")
  )

  existing_path <- candidate_paths[file.exists(candidate_paths)][1]

  if (is.na(existing_path)) {
    stop(
      "Config file not found. Expected one of: ",
      paste(candidate_paths, collapse = ", ")
    )
  }

  existing_path
}

normalize_pipeline_experimental_design <- function(experimental_design) {
  if (!is.data.frame(experimental_design)) {
    stop("experimental_design must be a data.frame")
  }

  if ("replicate" %in% names(experimental_design) && !"bio_replicate" %in% names(experimental_design)) {
    experimental_design$bio_replicate <- experimental_design$replicate
  }

  if ("bio_replicate" %in% names(experimental_design) && !"replicate" %in% names(experimental_design)) {
    experimental_design$replicate <- experimental_design$bio_replicate
  }

  if ("replicate" %in% names(experimental_design)) {
    experimental_design$replicate <- as.character(experimental_design$replicate)
  }

  if ("bio_replicate" %in% names(experimental_design)) {
    experimental_design$bio_replicate <- as.character(experimental_design$bio_replicate)
  }

  if ("measurement_batch" %in% names(experimental_design)) {
    experimental_design$measurement_batch <- as.factor(experimental_design$measurement_batch)
  }

  if ("plex" %in% names(experimental_design)) {
    experimental_design$plex <- as.factor(experimental_design$plex)
  }

  if ("measurement_batch" %in% names(experimental_design) && !"plex" %in% names(experimental_design)) {
    experimental_design$plex <- experimental_design$measurement_batch
  }

  if ("plex" %in% names(experimental_design) && !"measurement_batch" %in% names(experimental_design)) {
    experimental_design$measurement_batch <- experimental_design$plex
  }

  experimental_design
}

load_pipeline_experimental_design <- function(location_annotation, search_data_colnames, instrument) {
  experimental_design <- readr::read_delim(location_annotation, show_col_types = FALSE)
  experimental_design <- dplyr::filter(
    experimental_design,
    sample %in% stringr::str_subset(search_data_colnames, paste0("^", instrument))
  )

  normalize_pipeline_experimental_design(experimental_design)
}

filter_pipeline_experimental_design <- function(experimental_design,
                                                observed_samples,
                                                replicate_column = "replicate") {
  experimental_design <- normalize_pipeline_experimental_design(experimental_design)
  experimental_design <- dplyr::filter(experimental_design, sample %in% observed_samples)

  if (!replicate_column %in% names(experimental_design)) {
    return(experimental_design)
  }

  replicate_values <- experimental_design[[replicate_column]]
  replicate_values <- replicate_values[!is.na(replicate_values) & replicate_values != ""]

  if (!any(duplicated(replicate_values))) {
    message(
      "No repeated groups detected in '", replicate_column,
      "'. Skipping replicate-group filtering."
    )
    return(experimental_design)
  }

  experimental_design <- dplyr::group_by(
    experimental_design,
    dplyr::across(dplyr::all_of(replicate_column))
  )
  experimental_design <- dplyr::filter(experimental_design, dplyr::n() >= 2)
  dplyr::ungroup(experimental_design)
}

resolve_limma_blocking_values <- function(experimental_design,
                                          blocking_column = "replicate",
                                          dataset_name = NULL) {
  experimental_design <- normalize_pipeline_experimental_design(experimental_design)

  if (!blocking_column %in% names(experimental_design)) {
    stop(
      "Blocking is enabled but column '", blocking_column,
      "' was not found in the experimental design",
      if (!is.null(dataset_name)) paste0(" for dataset '", dataset_name, "'"),
      ". Available columns: ",
      paste(names(experimental_design), collapse = ", ")
    )
  }

  blocking_values <- experimental_design[[blocking_column]]

  if (any(is.na(blocking_values) | blocking_values == "")) {
    stop(
      "Blocking column '", blocking_column,
      "' contains missing or blank values",
      if (!is.null(dataset_name)) paste0(" for dataset '", dataset_name, "'"),
      "."
    )
  }

  block_sizes <- table(blocking_values)

  if (!any(block_sizes >= 2)) {
    stop(
      "Blocking is enabled but column '", blocking_column,
      "' does not contain repeated groups",
      if (!is.null(dataset_name)) paste0(" for dataset '", dataset_name, "'"),
      ". Disable blocking or provide a repeated-measures column."
    )
  }

  factor(blocking_values)
}

#' Smart data loader that handles different search engine outputs
#' @param search_engine Type of search engine ("diann", "fragpipe_tmt", "fragpipe_lf", "fragpipe_heavy_light", "spectronaut")
#' @param file_paths Named list of file paths for the specific search engine
#' @param parameters Named list of search engine specific parameters
#' @return Data frame with processed search engine data
load_search_engine_data <- function(search_engine, file_paths, parameters = list()) {
  
  result <- switch(search_engine,
    "diann" = {
      # Handle both parquet and TSV files
      file_ext <- tools::file_ext(file_paths$report)
      
      if (file_ext == "parquet") {
        diann_adapter_parquet(
          path_to_file = file_paths$report,
          proteotypic = parameters$proteotypic %||% TRUE,
          summarization = parameters$summarization %||% "SUM"
        )
      } else if (file_ext == "tsv") {
        result <- diann_adapter(
          path_to_file = file_paths$report,
          proteotypic = parameters$proteotypic %||% TRUE,
          summarization = parameters$summarization %||% "SUM"
        )
        # Clean column names for TSV files
        if (!is.null(parameters$instrument)) {
          result <- result %>% 
            rename_with(~ str_remove(., "_.*"), starts_with(parameters$instrument))
        }
        result
      } else {
        stop("Unsupported DIA-NN file format. Use .parquet or .tsv files.")
      }
    },
    
    "fragpipe_tmt" = {
      fragpipe_adapter(
        parent_dir = file_paths$parent_dir,
        ref_sample = parameters$ref_sample %||% "norm",
        grouping_var = parameters$grouping_var %||% "nterm_modif_peptide",
        min_purity = parameters$min_purity %||% 0.5,
        tmt_delta = parameters$tmt_delta %||% "229"
      )
    },
    
    "fragpipe_lf" = {
      fragpipe_lf_adapter(
        parent_dir = file_paths$parent_dir,
        annotation_file_path = file_paths$annotation,
        grouping_var = parameters$grouping_var %||% "nterm_modif_peptide"
      )
    },
    
    "fragpipe_heavy_light" = {
      fragpipe_dda_heavy_light_adapter(
        combined_modified_peptide_file = file_paths$combined_file,
        annotation_file_path = file_paths$annotation,
        grouping_var = parameters$grouping_var %||% "nterm_modif_peptide"
      )
    },
    
    "spectronaut" = {
      spectronaut_adapter(
        path_to_file = file_paths$report,
        proteotypic = parameters$proteotypic %||% TRUE
      )
    },
    
    stop("Unsupported search engine: ", search_engine)
  )
  
  return(result)
}

#' Create identification statistics plot adapted for different search engines
#' @param data Processed search engine data
#' @param experimental_design Sample annotation
#' @param search_engine Type of search engine used
#' @return List with plot and summary data
create_identification_plots_adaptive <- function(data, experimental_design, search_engine) {
  
  # Get quantitative columns based on search engine
  quant_cols <- switch(search_engine,
    "diann" = names(data)[str_detect(names(data), "^[A-Z]{2}\\d+|^\\d+")],  # Sample columns
    "fragpipe_tmt" = names(data)[str_detect(names(data), "TMT|^\\d+")],      # TMT channels
    "fragpipe_lf" = names(data)[!names(data) %in% c("nterm_modif_peptide", "nterm_modif", "peptide", "protein")],
    "fragpipe_heavy_light" = names(data)[str_detect(names(data), "heavy|light")],
    "spectronaut" = names(data)[!names(data) %in% c("nterm_modif_peptide", "nterm_modif", "peptide", "protein")],
    names(data)[sapply(data, is.numeric)]  # Fallback to numeric columns
  )
  
  if (length(quant_cols) == 0) {
    warning("No quantitative columns found for ", search_engine)
    return(NULL)
  }
  
  ids_per_sample <- data %>%
    select(all_of(quant_cols)) %>%
    pivot_longer(cols = everything(), names_to = "sample", values_to = "value") %>%
    group_by(sample) %>%
    summarise(n_identifications = sum(!is.na(value)), .groups = 'drop')
  
  # Join with experimental design if available
  if (nrow(experimental_design) > 0 && "sample" %in% names(experimental_design)) {
    ids_per_sample <- ids_per_sample %>%
      left_join(experimental_design, by = "sample")
    
    x_var <- if("condition" %in% names(experimental_design)) "condition" else "sample"
    fill_var <- if("condition" %in% names(experimental_design)) "condition" else NULL
  } else {
    x_var <- "sample"
    fill_var <- NULL
  }
  
  p1 <- ids_per_sample %>%
    ggplot(aes(x = .data[[x_var]], y = n_identifications)) +
    {if(!is.null(fill_var)) geom_boxplot(aes(fill = .data[[fill_var]])) else geom_boxplot()} +
    {if(!is.null(fill_var)) geom_jitter(width = 0.2, alpha = 0.7) else geom_point()} +
    labs(
      title = paste("Number of Identifications by", str_to_title(x_var), "-", str_to_upper(search_engine)),
      x = str_to_title(x_var), 
      y = "Number of Identifications"
    ) +
    theme_minimal() +
    {if(x_var == "sample") theme(axis.text.x = element_text(angle = 45, hjust = 1)) else theme()} +
    {if(is.null(fill_var)) theme(legend.position = "none") else theme()}
  
  return(list(plot = p1, data = ids_per_sample))
}

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

#' Complete partial missing values with robust sequential imputation
#' @param counts Numeric feature-by-sample matrix with finite observations and NAs
#' @return Complete matrix in the input order, preserving all observed values
impute_partial_missing_values <- function(counts) {
  if (!anyNA(counts)) {
    return(counts)
  }

  if (!requireNamespace("rrcovNA", quietly = TRUE)) {
    stop("Package 'rrcovNA' is required to impute remaining partial missing values. Install it before rerunning imputation.")
  }

  imputation_result <- tryCatch(
    rrcovNA::impSeqRob(counts),
    error = function(error) {
      stop("Partial-missing imputation failed: ", conditionMessage(error), call. = FALSE)
    }
  )

  # rrcovNA returns xseq on its EM fallback and x on its sequential path.
  imputed_counts <- imputation_result[["x"]]
  if (is.null(imputed_counts)) {
    imputed_counts <- imputation_result[["xseq"]]
  }

  if (!is.matrix(imputed_counts) || !is.numeric(imputed_counts) ||
      !identical(dim(imputed_counts), dim(counts)) ||
      !setequal(rownames(imputed_counts), rownames(counts)) ||
      !setequal(colnames(imputed_counts), colnames(counts))) {
    stop("Partial-missing imputation returned incompatible matrix dimensions or feature/sample names.")
  }

  imputed_counts <- imputed_counts[rownames(counts), colnames(counts), drop = FALSE]
  missing_values <- is.na(counts)
  if (any(!is.finite(imputed_counts[missing_values]))) {
    stop("Partial-missing imputation left missing or non-finite values. Check the retained data before rerunning.")
  }

  counts[missing_values] <- imputed_counts[missing_values]
  return(counts)
}

#' Impute missing values and return imputation summary
#' @param se SummarizedExperiment with a numeric, log2-scale matrix assay named 'counts'
#' @param min_fraction_condition Maximum missing fraction allowed within a condition
#'   before using minimum-probability imputation (numeric 0-1). Features exceeding
#'   this fraction in every condition are excluded. The historical argument name
#'   is retained for compatibility; fractions count samples, not biological IDs.
#' @param tune_sigma Non-negative multiplier for the median feature standard deviation
#' @param tune_quantile Sample-specific low quantile for minimum-probability sampling (0-1)
#' @param seed Non-negative integer random seed, default 1. NULL uses the caller's
#'   random stream. With a seed, the caller's R random state is restored on exit.
#' @return List with the imputed SummarizedExperiment, imputation summary table,
#'   and sparse-feature exclusion metadata. Total_Samples is the missingness
#'   denominator; Total_Replicates records distinct biological replicate IDs.
terminer_imputation <- function(se,
                                min_fraction_condition = 0.5,
                                tune_sigma = 1,
                                tune_quantile = 1e-7,
                                seed = 1L) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("SummarizedExperiment required")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("tibble required")

  parameters <- list(
    min_fraction_condition = min_fraction_condition,
    tune_sigma = tune_sigma,
    tune_quantile = tune_quantile
  )
  for (parameter_name in names(parameters)) {
    value <- parameters[[parameter_name]]
    upper_limit <- if (parameter_name == "tune_sigma") Inf else 1
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || value > upper_limit) {
      stop(parameter_name, " must be a finite numeric scalar between 0 and ", upper_limit, ".")
    }
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed < 0 || seed > .Machine$integer.max || seed != floor(seed))) {
    stop("seed must be NULL or a non-negative integer no greater than ", .Machine$integer.max, ".")
  }

  counts <- SummarizedExperiment::assay(se, "counts")
  if (!is.matrix(counts) || !is.numeric(counts)) {
    stop("Assay 'counts' must be a numeric matrix in the provided SummarizedExperiment.")
  }
  if (nrow(counts) == 0L || ncol(counts) == 0L) {
    stop("Assay 'counts' must contain at least one feature and one sample.")
  }
  if (any(!is.finite(counts) & !is.na(counts))) {
    stop("Assay 'counts' contains infinite values. Supply finite log2 abundances or NA.")
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
      anyNA(rownames(counts)) || anyNA(colnames(counts)) ||
      any(rownames(counts) == "") || any(colnames(counts) == "") ||
      anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
    stop("Assay 'counts' must have unique, non-empty feature and sample names.")
  }
  n_features_input <- nrow(counts)

  # Prepare sample metadata and require exact alignment with the assay.
  col_ann <- as.data.frame(SummarizedExperiment::colData(se))
  if (!"sample" %in% names(col_ann)) {
    col_ann$sample <- colnames(counts)
  } else {
    col_ann$sample <- as.character(col_ann$sample)
    missing_sample_names <- is.na(col_ann$sample) | col_ann$sample == ""
    col_ann$sample[missing_sample_names] <- colnames(counts)[missing_sample_names]
  }
  samples <- col_ann$sample
  if (!identical(samples, colnames(counts))) {
    stop("colData$sample must match the assay column names in the same order.")
  }
  if (!"condition" %in% names(col_ann) || anyNA(col_ann$condition) ||
      any(as.character(col_ann$condition) == "")) {
    stop("colData must include a 'condition' column without missing or blank values.")
  }
  if (!"replicate" %in% names(col_ann)) {
    if ("bio_replicate" %in% names(col_ann)) {
      col_ann$replicate <- col_ann$bio_replicate
    } else if ("Replicate" %in% names(col_ann)) {
      col_ann$replicate <- col_ann$Replicate
    } else {
      col_ann$replicate <- samples
    }
  }

  # Characterize missingness at the sample level within each condition.
  condition_sizes <- col_ann %>%
    dplyr::group_by(condition) %>%
    dplyr::summarise(
      Total_Samples = dplyr::n(),
      Total_Replicates = dplyr::n_distinct(replicate),
      .groups = "drop"
    )

  quant_peptide_data_long <- as.data.frame(counts) %>%
    tibble::rownames_to_column(var = "nterm_modif_peptide") %>%
    tidyr::pivot_longer(cols = -nterm_modif_peptide, names_to = "sample", values_to = "Abundance") %>%
    dplyr::left_join(col_ann, by = "sample")

  peptide_missingness <- quant_peptide_data_long %>%
    dplyr::group_by(nterm_modif_peptide, condition) %>%
    dplyr::summarise(
      Num_Quantified_per_cond = sum(!is.na(Abundance)),
      Num_Missing_per_cond = sum(is.na(Abundance)),
      .groups = "drop"
    ) %>%
    dplyr::left_join(condition_sizes, by = "condition") %>%
    dplyr::mutate(
      Proportion_Missing = Num_Missing_per_cond / Total_Samples,
      Missingness_Category = dplyr::case_when(
        Proportion_Missing > min_fraction_condition ~ "Total_Missing",
        Proportion_Missing > 0 ~ "Partial_Missing",
        TRUE ~ "Complete"
      )
    ) %>%
    dplyr::group_by(nterm_modif_peptide) %>%
    dplyr::mutate(all_conditions_missing = all(Proportion_Missing > min_fraction_condition)) %>%
    dplyr::ungroup()

  sparse_feature_missingness <- peptide_missingness %>%
    dplyr::filter(all_conditions_missing)
  sparse_features <- unique(sparse_feature_missingness$nterm_modif_peptide)
  n_sparse_features_excluded <- length(sparse_features)
  message(
    "Sparse-feature filter excluded ", n_sparse_features_excluded,
    " of ", n_features_input,
    " features that exceeded the missingness threshold in every condition."
  )

  sparse_feature_summary_table <- sparse_feature_missingness %>%
    dplyr::group_by(nterm_modif_peptide) %>%
    dplyr::summarise(
      n_conditions = dplyr::n_distinct(condition),
      min_quantified_per_condition = if (dplyr::n() > 0L) min(Num_Quantified_per_cond) else NA_integer_,
      max_quantified_per_condition = if (dplyr::n() > 0L) max(Num_Quantified_per_cond) else NA_integer_,
      mean_missing_fraction = mean(Proportion_Missing),
      .groups = "drop"
    )
  sparse_feature_annotations <- as.data.frame(
    SummarizedExperiment::rowData(se)[sparse_features, , drop = FALSE]
  )
  if (nrow(sparse_feature_annotations) > 0 && !"nterm_modif_peptide" %in% names(sparse_feature_annotations)) {
    sparse_feature_annotations$nterm_modif_peptide <- rownames(sparse_feature_annotations)
  }

  # Check the retained data before estimating any imputation distributions.
  retained_features <- setdiff(rownames(counts), sparse_features)
  if (length(retained_features) == 0L) {
    stop("No features remain after the sparse-feature filter. Check the missingness threshold and input data.")
  }
  retained_counts <- counts[retained_features, , drop = FALSE]
  empty_samples <- colnames(retained_counts)[colSums(!is.na(retained_counts)) == 0L]
  if (length(empty_samples) > 0L) {
    stop(
      "Cannot impute samples with no observed values after sparse-feature filtering: ",
      paste(empty_samples, collapse = ", "), ". Check or exclude these samples before rerunning."
    )
  }

  quant_peptide_data_long <- quant_peptide_data_long %>%
    dplyr::filter(nterm_modif_peptide %in% retained_features) %>%
    dplyr::left_join(peptide_missingness, by = c("nterm_modif_peptide", "condition")) %>%
    dplyr::mutate(abundance_imputed = Abundance, imputation_method = "not_imputed")

  if (!is.null(seed) && anyNA(retained_counts)) {
    had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_random_seed) {
      previous_random_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (had_random_seed) {
        assign(".Random.seed", previous_random_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  # Estimate the spread using features observed in at least half of all samples.
  total_missing <- is.na(quant_peptide_data_long$Abundance) &
    quant_peptide_data_long$Missingness_Category == "Total_Missing"
  if (any(total_missing)) {
    peptide_missingness_overall <- quant_peptide_data_long %>%
      dplyr::group_by(nterm_modif_peptide) %>%
      dplyr::summarise(
        Num_Quantified_overall = sum(!is.na(Abundance)),
        Num_Missing_overall = sum(is.na(Abundance)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        Proportion_Missing_overall = Num_Missing_overall /
          (Num_Quantified_overall + Num_Missing_overall)
      )
    sd_features <- peptide_missingness_overall %>%
      dplyr::filter(Proportion_Missing_overall <= 0.5, Num_Quantified_overall >= 2L) %>%
      dplyr::pull(nterm_modif_peptide)
    feature_wise_sd <- quant_peptide_data_long %>%
      dplyr::filter(nterm_modif_peptide %in% sd_features) %>%
      dplyr::group_by(nterm_modif_peptide) %>%
      dplyr::summarise(sd = stats::sd(Abundance, na.rm = TRUE), .groups = "drop")
    sd_median <- stats::median(feature_wise_sd$sd, na.rm = TRUE) * tune_sigma
    if (!is.finite(sd_median) || sd_median < 0) {
      stop("Cannot estimate a finite imputation standard deviation from features observed in at least half of the samples and at least two samples.")
    }

    min_quantile_sample <- quant_peptide_data_long %>%
      dplyr::group_by(sample) %>%
      dplyr::summarise(
        min_per_sample = stats::quantile(Abundance, probs = tune_quantile, na.rm = TRUE),
        .groups = "drop"
      )
    quant_peptide_data_long <- quant_peptide_data_long %>%
      dplyr::left_join(min_quantile_sample, by = "sample")
    imputation_means <- quant_peptide_data_long$min_per_sample[total_missing]
    if (any(!is.finite(imputation_means))) {
      stop("Cannot estimate a finite sample-specific imputation floor.")
    }

    # Draw exactly one Gaussian value for each missing observation.
    imputed_values <- stats::rnorm(n = sum(total_missing), mean = imputation_means, sd = sd_median)
    if (any(!is.finite(imputed_values))) {
      stop("Minimum-probability imputation produced non-finite values. Check tune_sigma and the input abundances.")
    }
    quant_peptide_data_long$abundance_imputed[total_missing] <- imputed_values
    quant_peptide_data_long$imputation_method[total_missing] <- "minProb_dist"
  }

  # Only remaining missing observations require the robust imputation step.
  partial_missing_counts <- quant_peptide_data_long %>%
    dplyr::select(nterm_modif_peptide, sample, abundance_imputed) %>%
    tidyr::pivot_wider(names_from = sample, values_from = abundance_imputed) %>%
    dplyr::select(nterm_modif_peptide, dplyr::all_of(samples)) %>%
    tibble::column_to_rownames(var = "nterm_modif_peptide") %>%
    as.matrix()
  imputed_counts <- impute_partial_missing_values(partial_missing_counts)

  # Label methods only after successful replacement and retain the audit columns.
  imputed_data_long <- as.data.frame(imputed_counts) %>%
    tibble::rownames_to_column(var = "nterm_modif_peptide") %>%
    tidyr::pivot_longer(cols = -nterm_modif_peptide, names_to = "sample", values_to = "Abundance") %>%
    dplyr::left_join(
      dplyr::select(quant_peptide_data_long, -Abundance),
      by = c("nterm_modif_peptide", "sample")
    ) %>%
    dplyr::mutate(
      imputation_method = dplyr::if_else(is.na(abundance_imputed), "impSeqRob", imputation_method)
    )
  imputation_summary_table <- imputed_data_long %>%
    dplyr::select(
      nterm_modif_peptide, sample,
      dplyr::any_of(c(
        "sample_name", "Num_Quantified_per_cond", "Num_Missing_per_cond",
        "Proportion_Missing", "Total_Samples", "Total_Replicates",
        "Missingness_Category", "imputation_method"
      ))
    )

  imputed_counts <- imputed_counts[retained_features, samples, drop = FALSE]
  se_imputed <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = imputed_counts),
    colData = SummarizedExperiment::colData(se),
    rowData = SummarizedExperiment::rowData(se)[retained_features, , drop = FALSE]
  )

  return(list(
    se = se_imputed,
    imputation_summary_table = imputation_summary_table,
    n_sparse_features_excluded = n_sparse_features_excluded,
    sparse_features = sparse_features,
    sparse_feature_summary_table = sparse_feature_summary_table,
    sparse_feature_missingness = sparse_feature_missingness,
    sparse_feature_annotations = sparse_feature_annotations
  ))
}
