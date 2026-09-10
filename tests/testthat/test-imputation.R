test_that("the SD estimate excludes peptides missing in more than half of samples", {
  helpers <- load_imputation_helpers()
  counts <- make_imputation_counts()
  counts[1:40, 1:4] <- NA_real_
  input <- make_imputation_se(counts)
  result <- suppressMessages(helpers$terminer_imputation(input, seed = 12L))

  # Only the 20 complete peptides qualify: the other 40 are missing 4/6 samples.
  expected_sd <- median(apply(counts[41:60, ], 1, sd))
  floors <- apply(counts, 2, quantile, probs = 1e-7, na.rm = TRUE)
  expected_values <- withr::with_seed(12L, rnorm(120, mean = rep(floors[1:3], 40), sd = expected_sd))
  expect_equal(as.vector(t(assay(result$se)[1:40, 1:3])), expected_values)
  expect_true(all(is.finite(assay(result$se))))
  expect_identical(assay(result$se)[!is.na(counts)], counts[!is.na(counts)])
  expect_equal(sum(result$imputation_summary_table$imputation_method == "minProb_dist"), 120)
  expect_equal(sum(result$imputation_summary_table$imputation_method == "impSeqRob"), 40)
})

test_that("technical repeats use sample counts without changing biological replicate metadata", {
  helpers <- load_imputation_helpers()
  counts <- make_imputation_counts(n_samples = 8L)
  counts[1, 1:4] <- NA_real_
  counts[2, c(1, 2, 5, 6)] <- NA_real_
  input <- make_imputation_se(counts, replicates = rep(c(1, 1, 2, 2), 2))
  result <- suppressMessages(helpers$terminer_imputation(input))
  summary <- result$imputation_summary_table
  missing_group <- filter(summary, nterm_modif_peptide == "peptide_1", sample == "sample_1")
  expect_equal(missing_group$Proportion_Missing, 1)
  expect_equal(missing_group$Total_Samples, 4L)
  expect_equal(missing_group$Total_Replicates, 2L)
  expect_equal(missing_group$Missingness_Category, "Total_Missing")
  expect_equal(missing_group$imputation_method, "minProb_dist")
  expect_true("peptide_2" %in% rownames(result$se))
  expect_true(all(summary$Proportion_Missing >= 0 & summary$Proportion_Missing <= 1))
  expect_identical(colData(result$se), colData(input))
})

test_that("empty samples are rejected before robust imputation, including after filtering", {
  helpers <- load_imputation_helpers()
  helpers$impute_partial_missing_values <- function(counts) stop("Unexpected robust imputation")
  counts <- make_imputation_counts()
  counts[, 1] <- NA_real_
  expect_error(suppressMessages(helpers$terminer_imputation(make_imputation_se(counts))),
               "no observed values.*sample_1")
  counts[1, 1] <- 20
  counts[1, 2:6] <- NA_real_
  expect_error(suppressMessages(helpers$terminer_imputation(make_imputation_se(counts))),
               "no observed values after sparse-feature filtering: sample_1")
})

test_that("missing rrcovNA cannot produce an apparently imputed result", {
  helpers <- load_imputation_helpers()
  helpers$requireNamespace <- function(package, ...) {
    if (package == "rrcovNA") return(FALSE)
    base::requireNamespace(package, ...)
  }
  counts <- make_imputation_counts()
  counts[1, 1] <- NA_real_
  expect_error(suppressMessages(helpers$terminer_imputation(make_imputation_se(counts))),
               "rrcovNA.*required.*remaining partial missing")
})

test_that("complete small inputs bypass robust imputation and keep their metadata", {
  helpers <- load_imputation_helpers()
  helpers$requireNamespace <- function(package, ...) {
    if (package == "rrcovNA") stop("Unexpected robust dependency lookup")
    base::requireNamespace(package, ...)
  }
  for (n_features in c(1L, 2L)) {
    counts <- make_imputation_counts(n_features)
    input <- make_imputation_se(counts)
    expect_no_warning(result <- suppressMessages(helpers$terminer_imputation(input)))
    expect_identical(assay(result$se), counts)
    expect_identical(rowData(result$se), rowData(input))
    expect_identical(colData(result$se), colData(input))
    expect_true(all(result$imputation_summary_table$imputation_method == "not_imputed"))
    expect_equal(nrow(result$sparse_feature_summary_table), 0L)
  }
})

test_that("minimum-probability sampling consumes one Gaussian draw per replacement", {
  helpers <- load_imputation_helpers()
  counts <- make_imputation_counts()
  counts[1:10, 1:3] <- NA_real_
  input <- make_imputation_se(counts)
  withr::local_seed(92L)
  expected_seed <- withr::with_seed(92L, { rnorm(30); .Random.seed })
  result <- suppressMessages(helpers$terminer_imputation(input, seed = NULL))
  expect_identical(.Random.seed, expected_seed)
  expect_equal(sum(result$imputation_summary_table$imputation_method == "minProb_dist"), 30)
  expect_true(all(is.finite(assay(result$se))))
})

test_that("seeded mixed imputation is reproducible and restores the caller random stream", {
  helpers <- load_imputation_helpers()
  counts <- make_imputation_counts()
  counts[1, 1:4] <- NA_real_
  counts[2, 1] <- NA_real_
  input <- make_imputation_se(counts)
  withr::local_seed(1234L)
  previous_seed <- .Random.seed
  first <- suppressMessages(helpers$terminer_imputation(input, seed = 7L))
  expect_identical(.Random.seed, previous_seed)
  runif(5)
  second <- suppressMessages(helpers$terminer_imputation(input, seed = 7L))
  third <- suppressMessages(helpers$terminer_imputation(input, seed = 8L))
  expect_identical(assay(first$se), assay(second$se))
  expect_false(identical(assay(first$se)[1, 1:3], assay(third$se)[1, 1:3]))
})

test_that("seeded imputation restores RNG state even when imputation fails", {
  helpers <- load_imputation_helpers()
  helpers$requireNamespace <- function(package, ...) {
    if (package == "rrcovNA") return(FALSE)
    base::requireNamespace(package, ...)
  }
  counts <- make_imputation_counts()
  counts[1, 1:4] <- NA_real_
  withr::local_seed(98L)
  previous_seed <- .Random.seed
  expect_error(suppressMessages(helpers$terminer_imputation(make_imputation_se(counts))), "rrcovNA")
  expect_identical(.Random.seed, previous_seed)
})

test_that("sparse exclusions preserve surviving row order and annotation", {
  helpers <- load_imputation_helpers()
  counts <- make_imputation_counts()
  counts[1, ] <- NA_real_
  counts[2, c(1, 2, 4, 5)] <- NA_real_
  input <- make_imputation_se(counts)
  result <- suppressMessages(helpers$terminer_imputation(input))
  expect_identical(rownames(result$se), rownames(input)[-(1:2)])
  expect_identical(rowData(result$se), rowData(input)[-(1:2), , drop = FALSE])
  expect_setequal(result$sparse_features, rownames(input)[1:2])
  expect_equal(result$n_sparse_features_excluded, 2L)
  expect_equal(nrow(result$sparse_feature_summary_table), 2L)
  expect_setequal(rownames(result$sparse_feature_annotations), rownames(input)[1:2])
  expect_true(all(result$imputation_summary_table$imputation_method == "not_imputed"))
})

test_that("an empty retained dataset produces an actionable error", {
  helpers <- load_imputation_helpers()
  counts <- make_imputation_counts()
  counts[,] <- NA_real_
  expect_error(suppressMessages(helpers$terminer_imputation(make_imputation_se(counts))),
               "No features remain after the sparse-feature filter")
})

test_that("undefined spread is rejected only when minimum-probability imputation needs it", {
  helpers <- load_imputation_helpers()
  counts <- make_imputation_counts(6L)
  counts[,] <- NA_real_
  diag(counts) <- 20
  input <- make_imputation_se(counts)
  expect_error(suppressMessages(helpers$terminer_imputation(input, min_fraction_condition = 0.8)),
               "Cannot estimate a finite imputation standard deviation")
})

test_that("invalid parameters and infinite input are rejected early", {
  helpers <- load_imputation_helpers()
  counts <- make_imputation_counts()
  input <- make_imputation_se(counts)
  expect_error(helpers$terminer_imputation(input, min_fraction_condition = 1.1), "min_fraction_condition")
  expect_error(helpers$terminer_imputation(input, tune_quantile = NA_real_), "tune_quantile")
  expect_error(helpers$terminer_imputation(input, tune_sigma = -1), "tune_sigma")
  expect_error(helpers$terminer_imputation(input, seed = 1.2), "seed must be NULL")
  expect_error(helpers$terminer_imputation(input, seed = c(1, 2)), "seed must be NULL")
  counts[1, 1] <- Inf
  expect_error(helpers$terminer_imputation(make_imputation_se(counts)), "contains infinite values")
})

test_that("sample and condition metadata are validated before joining", {
  helpers <- load_imputation_helpers()
  input <- make_imputation_se(make_imputation_counts())
  colData(input)$sample <- rev(colData(input)$sample)
  expect_error(helpers$terminer_imputation(input), "sample must match.*same order")
  colData(input)$sample <- colnames(input)
  colData(input)$condition[1] <- NA_character_
  expect_error(helpers$terminer_imputation(input), "condition.*without missing or blank")
})

test_that("non-syntactic feature and sample names survive imputation", {
  helpers <- load_imputation_helpers()
  counts <- make_imputation_counts()
  rownames(counts)[1] <- "Acetyl peptide+1"
  colnames(counts)[1:2] <- c("sample 1", "run-2")
  counts[1, 1] <- NA_real_
  input <- make_imputation_se(counts)
  result <- suppressMessages(helpers$terminer_imputation(input))
  expect_identical(dimnames(assay(result$se)), dimnames(counts))
  expect_identical(assay(result$se)[!is.na(counts)], counts[!is.na(counts)])
})

test_that("the EM fallback preserves feature identities and observed values", {
  helpers <- load_imputation_helpers()
  counts <- make_imputation_counts()
  counts[cbind(seq_len(nrow(counts)), rep(seq_len(ncol(counts)), length.out = nrow(counts)))] <- NA_real_
  counts[1, 4] <- NA_real_
  input <- make_imputation_se(counts)
  result <- suppressMessages(helpers$terminer_imputation(input))
  expect_identical(dimnames(assay(result$se)), dimnames(counts))
  expect_true(all(is.finite(assay(result$se))))
  expect_identical(assay(result$se)[!is.na(counts)], counts[!is.na(counts)])
  expect_equal(sum(result$imputation_summary_table$imputation_method == "impSeqRob"), sum(is.na(counts)))
})

test_that("robust backend output is aligned by names and only missing cells are replaced", {
  helpers <- load_imputation_helpers()
  counts <- make_imputation_counts()
  counts[1, 1] <- NA_real_
  backend_result <- counts + 100
  backend_result[1, 1] <- 15
  backend_result <- backend_result[rev(rownames(counts)), rev(colnames(counts))]
  local_mocked_bindings(impSeqRob = function(x) list(xseq = backend_result), .package = "rrcovNA")
  result <- helpers$impute_partial_missing_values(counts)
  expect_identical(dimnames(result), dimnames(counts))
  expect_identical(result[!is.na(counts)], counts[!is.na(counts)])
  expect_equal(result[1, 1], 15)
})

test_that("invalid or incomplete robust backend output stops with a clear error", {
  helpers <- load_imputation_helpers()
  counts <- make_imputation_counts()
  counts[1, 1] <- NA_real_
  backend_result <- list(x = counts)
  local_mocked_bindings(impSeqRob = function(x) backend_result, .package = "rrcovNA")
  expect_error(helpers$impute_partial_missing_values(counts), "left missing or non-finite values")
  backend_result$x[1, 1] <- Inf
  expect_error(helpers$impute_partial_missing_values(counts), "left missing or non-finite values")
  backend_result <- list(x = counts[-1, , drop = FALSE])
  expect_error(helpers$impute_partial_missing_values(counts), "incompatible matrix dimensions")
})

test_that("robust backend failures retain their cause in the pipeline error", {
  helpers <- load_imputation_helpers()
  counts <- make_imputation_counts()
  counts[1, 1] <- NA_real_
  local_mocked_bindings(impSeqRob = function(x) stop("singular covariance"), .package = "rrcovNA")
  expect_error(helpers$impute_partial_missing_values(counts),
               "Partial-missing imputation failed: singular covariance")
})

test_that("a seeded call restores an initially absent R random state", {
  helpers <- load_imputation_helpers()
  counts <- make_imputation_counts()
  counts[1, 1:3] <- NA_real_
  withr::local_seed(42L)
  rm(".Random.seed", envir = .GlobalEnv)
  result <- suppressMessages(helpers$terminer_imputation(make_imputation_se(counts)))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_true(all(is.finite(assay(result$se))))
})
