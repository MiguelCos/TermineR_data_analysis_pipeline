# Copilot Instructions for TermineR Data Analysis Pipeline

## Big picture architecture
- This repo is an R/Quarto pipeline with three staged reports: `terminer_exploratory_analysis.qmd` → `terminer_inferential_analysis.qmd` → `terminer_results_visualization.qmd`.
- Exploratory analysis loads search-engine outputs (DIA-NN, FragPipe TMT/LF/HL, Spectronaut), annotates with TermineR, and writes cached `rds/` artifacts used by downstream scripts.
- Inferential analysis consumes cached RDS objects for limma-based differential analysis and GO enrichment; visualization consumes those results to build publication plots.

## Key files and where logic lives
- `terminer_exploratory_analysis.qmd`: parameter block, adapter selection via `search_engine`, and initial caching (`*_data.rds`).
- `scr/helper_functions.R`: shared helpers (`load_search_engine_data()`, PCA/volcano helpers, `safe_rds_read()` cache wrapper).
- `scr/modif_diann_adapter.R`: DIA-NN parquet/TSV adapter logic used when `search_engine == "diann"`.
- `terminer_inferential_analysis.qmd` and `terminer_results_visualization.qmd`: downstream stages (assume inputs from `rds/` and outputs in `results/`).

## Data flow & conventions
- File paths are set with `here()` and assume repo-root execution; keep new paths aligned with `data/`, `rds/`, and `results/`.
- `search_engine` drives a `switch()` in the exploratory report; keep new adapters consistent with this pattern.
- DIA-NN: `.parquet` vs `.tsv` is auto-detected; TSV requires column renaming based on `instrument` prefix.
- FragPipe TMT: expects `mix_*/psm.tsv` under `data/` plus an `experiment_annotation.tsv` (see `data/example_fragpipe_tmt_search/`).
- Experimental annotation must contain `sample`, `sample_name`, `condition`, `bio_replicate`; sample IDs must match search output columns.
- RDS caching is standard: read if present, compute + `write_rds()` otherwise; deleting `rds/` forces recomputation.

## Developer workflow (R/Quarto)
- Run the pipeline in order: exploratory → inferential → results visualization.
- Quarto chunks use `cache: true` and suppress messages; keep chunk names stable to preserve caching.
- External deps include TermineR, diann, Bioconductor packages (limma, clusterProfiler, SummarizedExperiment, etc.) as listed in `README.md`.

## Project-specific patterns
- Prefer helper functions in `scr/helper_functions.R` over duplicating adapter logic.
- Use `safe_rds_read()` for any expensive computation to align with existing caching style.
- Results are TSVs in `results/`; intermediate objects in `rds/`.