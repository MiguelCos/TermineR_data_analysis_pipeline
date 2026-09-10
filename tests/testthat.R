# Run from the repository root with: Rscript --vanilla tests/testthat.R
testthat::test_dir(file.path("tests", "testthat"), reporter = "summary", stop_on_failure = TRUE)
