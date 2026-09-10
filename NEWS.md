# TermineR Data Analysis Pipeline News

## Unreleased

### Fixed
- Corrected overall imputation missingness to use sample counts and prevented repeated biological replicate IDs from distorting condition-level missingness and sparse-feature filtering.
- Added validation for empty retained data, unobserved samples, undefined imputation distributions, and unsuccessful robust imputation; method labels now describe successful replacements.
- Reduced minimum-probability sampling to one Gaussian draw per missing observation and skipped robust imputation for complete matrices.

### Added
- Added `imputation_seed` (default `1`) to make fresh imputation runs reproducible while restoring the caller's R random state.
- Added `Total_Samples` to imputation summaries and synthetic regression tests for imputation behavior and failure cases.

### Changed
- Centralized pipeline parameters in `config.yaml` and updated notebooks to load stage-specific settings through shared helper functions.
- Updated neo-termini annotation calls to rely on the current `TermineR::annotate_neo_termini()` TargetP-aware processing annotations.
- Removed the external TargetP2 post-processing step and the obsolete TargetP2 file path configuration.
- Standardized cached artifact naming around `pre_fix` and `search_engine` for upstream preparation outputs.
- Improved experimental-design normalization for shared `replicate`, `bio_replicate`, `measurement_batch`, and `plex` handling.
- Added optional GO analysis and blocking controls to the inferential stage configuration.
- Updated result visualization setup to use shared configuration and handle missing GO grouping outputs more gracefully.

### Documentation
- Expanded README configuration guidance with YAML examples, stage-specific parameters, and updated TargetP annotation notes.
- Refreshed Copilot instructions to describe the current four-stage pipeline and configuration-first workflow.
