suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

load_imputation_helpers <- function() {
  helpers <- new.env(parent = globalenv())
  sys.source(file.path("..", "..", "scr", "helper_functions.R"), envir = helpers)
  helpers
}

make_imputation_counts <- function(n_features = 60L, n_samples = 6L) {
  withr::with_seed(123, matrix(
    rnorm(n_features * n_samples, mean = 20, sd = 2),
    nrow = n_features,
    dimnames = list(paste0("peptide_", seq_len(n_features)), paste0("sample_", seq_len(n_samples)))
  ))
}

make_imputation_se <- function(counts, replicates = seq_len(ncol(counts))) {
  SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
      sample = colnames(counts),
      sample_name = paste("Sample", seq_len(ncol(counts))),
      condition = rep(c("A", "B"), each = ncol(counts) / 2),
      bio_replicate = replicates
    ),
    rowData = S4Vectors::DataFrame(
      annotation = paste("Annotation", seq_len(nrow(counts))),
      row.names = rownames(counts)
    )
  )
}
