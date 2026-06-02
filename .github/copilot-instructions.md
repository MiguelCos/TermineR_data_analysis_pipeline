# Copilot Instructions for TermineR Data Analysis Pipeline

## Start Here
- Use [README.md](../README.md) for dependency installation, supported search engines, and input file requirements.
- Treat [config.yaml](../config.yaml) as the source of truth for pipeline parameters. Prefer changing config values over hardcoding notebook variables.
- Configuration lookup is not root-only: `default_pipeline_config_path()` in [scr/helper_functions.R](../scr/helper_functions.R) checks `r_project_files/config.yaml` before the repo-root `config.yaml`.

## Pipeline Shape
- This repo is a four-stage R/Quarto workflow and should be handled in order: `terminer_data_preparation.qmd` -> `terminer_exploratory_analysis.qmd` -> `terminer_inferential_analysis.qmd` -> `terminer_results_visualization.qmd`.
- The notebooks are the execution surface. There is no CI, test harness, or scripted build pipeline in the repo.
- `save_key_objects.R` is a side utility for exporting inferential objects, not a required stage.

## Where Behavior Lives
- [scr/helper_functions.R](../scr/helper_functions.R): shared config loading, path resolution, experimental-design normalization, caching helpers, and plotting utilities.
- [scr/modif_diann_adapter.R](../scr/modif_diann_adapter.R): DIA-NN-specific loading and column handling.
- [terminer_data_preparation.qmd](../terminer_data_preparation.qmd): primary search-engine switch and initial cache creation.
- Downstream notebooks mostly consume cached `.rds` objects and `results/` tables rather than reimplementing loaders.

## Repo Conventions
- Keep paths repo-relative through `here()` unless the config intentionally uses absolute paths.
- `search_engine` is the main routing knob. Supported values are `diann`, `fragpipe_tmt`, `fragpipe_lf`, `fragpipe_heavy_light`, and `spectronaut`; keep adapter changes aligned with that switch-based pattern.
- `instrument` is used as a prefix filter for quantitative sample columns. If it does not match the actual sample columns, sample selection can fail quietly.
- Experimental annotations must align exactly with search-output sample names. The helper layer normalizes legacy `replicate` and `bio_replicate` columns, so prefer reusing that behavior instead of open-coding column fixes.
- Cached filenames are coupled to both `pre_fix` and `search_engine`. Changing either breaks downstream cache reads until artifacts are regenerated.
- Quarto chunk caching is enabled broadly with stable chunk names. Avoid renaming chunks unless you want to invalidate cached execution.

## Editing Guidance
- For new parameters, extend [config.yaml](../config.yaml) under `common` or the matching `stages.<stage>` block, then load them through `load_stage_parameters()`.
- Prefer existing helpers over duplicate notebook logic, especially for config resolution, sample/annotation normalization, and cache access.
- If a change affects adapter outputs or cache naming, verify both [terminer_data_preparation.qmd](../terminer_data_preparation.qmd) and [terminer_exploratory_analysis.qmd](../terminer_exploratory_analysis.qmd).
- If notebook output looks stale, inspect `rds/` before assuming the analysis code is wrong.