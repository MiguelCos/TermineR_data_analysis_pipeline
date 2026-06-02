# TermineR Data Analysis Pipeline

This repository contains a streamlined workflow for the analysis of shotgun proteomics data from multiple search engines, specifically designed for the analysis of proteolytic processing events using the TermineR approach.

## Overview

The pipeline consists of three main analysis scripts that support multiple search engine outputs:


1. **Data preparation**: Data loading and preparation, including missing value inputation if necessary.
1. **Exploratory Analysis**: Data loading, quality control, and initial visualization.
2. **Inferential Analysis**: Statistical analysis, protein normalization, and GO enrichment.  
3. **Results Visualization**: Generation of publication-ready plots and summary statistics.

## Supported search engines

The pipeline supports the following search engine outputs through integrated TermineR adapters:

### 1. DIA-NN (Data-Independent Acquisition)
- **File types**: `report.parquet` or `report.tsv`
- **Quantification**: Precursor-level quantification with MAD scaling
- **Modifications**: Automatic N-terminal modification detection if annotated in unimod.

### 2. FragPipe TMT (Tandem Mass Tags)
- **File types**: `psm.tsv` files in mixture directories
- **Quantification**: TMT reporter ion intensities with reference channel normalization
- **Modifications**: N-terminal TMT, Acetyl, Dimethyl detection

### 3. FragPipe Label-Free
- **File types**: `psm.tsv` files with intensity-based quantification
- **Quantification**: Precursor intensity with MAD scaling
- **Modifications**: Acetyl, Dimethyl, TMT, 2PCA, Pyro-Glu

### 4. FragPipe Heavy-Light (SILAC/Dimethyl)
- **File types**: `combined_modified_peptide_label_quant.tsv`
- **Quantification**: Heavy/Light label quantification
- **Modifications**: Dimethyl, Acetyl, Pyro-Glu

### 5. Spectronaut
- **File types**: Spectronaut report TSV files
- **Quantification**: Intensity-based with log2 transformation
- **Modifications**: Acetyl detection (extensible)

## Requirements

### R packages

The following R packages are required:

**Bioconductor packages:**
- `limma`
- `clusterProfiler` 
- `SummarizedExperiment`
- `ComplexHeatmap`
- `dagLogo`
- `org.Hs.eg.db` (or appropriate organism database)

**CRAN packages:**
- `tidyverse`
- `here`
- `yaml`
- `ggpubr`
- `pheatmap`
- `RColorBrewer`
- `mixOmics`
- `naniar`
- `arrow`
- `kableExtra`

**GitHub packages:**
- `TermineR`: `devtools::install_github("MiguelCos/TermineR")`
- `diann`: `devtools::install_github("vdemichev/diann-rpackage")`

## Input data requirements

### 1. DIA-NN report file
- **Location**: `data/report.parquet` or `data/report.tsv`
- **Description**: Output from FragPipe-DIA + DIANN quantitation
- **Format**: Parquet or TSV file containing peptide identifications and quantitative data
- **Note**: The pipeline automatically detects the file format and uses the appropriate loader:
  - `.parquet` files: Uses optimized parquet adapter with MAD scaling
  - `.tsv` files: Uses standard DIA-NN adapter

### 2. Experimental annotation file
- **Location**: `data/experimental_annotation.txt`
- **Description**: Tab-delimited file describing the experimental design and sample metadata
- **Format**: Tab-separated values (.txt or .tsv)

#### Required columns
- **`sample`**: Sample identifier that must exactly match the column names in your search engine output files
- **`sample_name`**: Human-readable sample name or identifier  
- **`condition`**: Experimental condition or treatment group (e.g., "Control", "Treatment", "TU", "NAT", "PC")
- **`bio_replicate`**: Biological replicate identifier
  - **Important**: Use the same `bio_replicate` number when samples represent repeated measures of the same biological entity (e.g., technical replicates, fractionation, or multiple time points from the same individual)
  - Use different `bio_replicate` numbers for independent biological samples
- Additional metadata columns can be included as needed

#### Annotation file guidelines
1. **Header required**: First row must contain column names
2. **Tab-delimited**: Use tab characters to separate columns
3. **No missing values**: All required columns must have values for each sample
4. **Consistent naming**: Sample identifiers must match exactly between annotation and data files
5. **Biological replicates**: 
   - Same `bio_replicate` = repeated measures from same biological source
   - Different `bio_replicate` = independent biological samples
   - This is critical for proper statistical modeling in limma

#### Example annotation files
See `example_annotation.txt` in the repository `data/` for a template with the minimal required columns.

**Basic example:**

```
sample	sample_name	condition	bio_replicate
TF658	T_001	TU	1
TF659	N_002	NAT	2
TF660	T_002	TU	2
TF661	N_003	NAT	3
TF662	C_003	PC	3
TF663	T_003	TU	3
```

**Extended example with additional metadata:**
```
sample	sample_name	condition	bio_replicate	batch	treatment_time
EX001	Control_1	Control	1	A	24h
EX002	Control_2	Control	2	A	24h
EX003	Control_3	Control	3	B	24h
EX004	Treatment_1	Treatment	1	A	24h
EX005	Treatment_2	Treatment	2	A	24h
EX006	Treatment_3	Treatment	3	B	24h
```

### 3. Protein FASTA file
- **Location**: `data/proteome.fasta`
- **Description**: Protein sequence database used for the search
- **Format**: Standard FASTA format

### 4. TargetP processing annotations
- **Description**: TargetP-aware processing annotations are provided by the updated `TermineR::annotate_neo_termini()` for supported organisms.
- **Configuration**: No external TargetP2 output file is required by this pipeline.

## Configuration parameters

All pipeline notebooks now load their parameters from `config.yaml` at the repository root. Update that file once per project instead of editing each notebook separately.

Paths in `config.yaml` can be either:
- relative to the repository root, for example `data/report.tsv`
- absolute paths, for example `Z:/project/data/report.tsv`

Use `null` for optional inputs that are not relevant to your workflow.

### Configuration structure
```yaml
common:
  search_engine: diann
  pre_fix: mc001_project_name
  rds_dir: rds
  results_dir: results

stages:
  data_preparation:
    plot_width: 10
  exploratory:
    plot_width: 10
  inferential:
    defined_contrasts:
      Treatment_vs_Control: Treatment - Control
  results_visualization:
    condition_levels:
      - Control
      - Treatment
```

The `common` section stores parameters shared across the pipeline, while `stages` contains notebook-specific overrides for `data_preparation`, `exploratory`, `inferential`, and `results_visualization`.

### Common parameters

| Parameter | Description | Typical value |
| --- | --- | --- |
| `search_engine` | Search-engine adapter to use. Supported values are `diann`, `fragpipe_tmt`, `fragpipe_lf`, `fragpipe_heavy_light`, and `spectronaut`. | `diann` |
| `diann_report_location` | Path to the DIA-NN `report.tsv` or `report.parquet` file. Used only when `search_engine: diann`. | `data/report.tsv` |
| `fragpipe_parent_dir` | Parent directory containing FragPipe TMT mix folders such as `mix_1/`, `mix_2/`. Used only for TMT workflows. | `data/fragpipe_tmt_search` |
| `fragpipe_lf_parent_dir` | Parent directory for FragPipe label-free outputs. Used only for label-free workflows. | `data/fragpipe_lf_search` |
| `fragpipe_lf_annotation` | Annotation file for FragPipe label-free workflows. | `data/fragpipe_lf_annotation.txt` |
| `fragpipe_hl_file` | Path to the `combined_modified_peptide_label_quant.tsv` file for heavy-light workflows. | `data/combined_modified_peptide_label_quant.tsv` |
| `fragpipe_hl_annotation` | Annotation file for heavy-light workflows. | `data/fragpipe_hl_annotation.txt` |
| `spectronaut_report` | Path to the Spectronaut report file. Used only for Spectronaut workflows. | `data/spectronaut_report.tsv` |
| `location_annotation` | Experimental design file used to define samples, conditions, and biological replicates. | `data/experimental_annotation.txt` |
| `fasta_location` | FASTA file used for TermineR peptide annotation. | `data/proteome.fasta` |
| `uniprot_annotation` | UniProt annotation table used to build protein-to-gene mappings during upstream processing. | `data/uniprot_annotation.tsv` |
| `protein_annotation_path` | UniProt annotation table used downstream in inferential and visualization steps for protein descriptions and gene names. | `data/uniprot_annotation.tsv` |
| `rds_dir` | Directory where cached intermediate `.rds` objects are written. | `rds` |
| `results_dir` | Directory where tabular outputs are written. | `results` |
| `sense_protease` | Cleavage direction expected by TermineR. Use `C` for trypsin-like C-terminal cleavage and `N` for N-terminal cleavage. | `C` |
| `specificity_protease` | Protease specificity pattern passed to TermineR annotation. | `K|R` |
| `organism_annotation` | Organism keyword passed to TermineR annotation functions. | `human` |
| `instrument` | Prefix used to identify quantitative sample columns in search-engine outputs. | `TF` |
| `ref_sample` | Reference channel name for FragPipe TMT normalization. | `norm` |
| `min_purity` | Minimum reporter-ion purity accepted for FragPipe TMT PSMs. | `0.5` |
| `tmt_delta` | TMT delta mass used by the FragPipe TMT adapter. Use `229` for TMT10/11 and `304` for TMT16. | `"229"` |
| `proteotypic_only` | If `true`, keep only proteotypic peptides in DIA-NN and Spectronaut workflows. | `true` |
| `summarization_method` | Precursor summarization used in DIA-NN loading. Supported values are `SUM` and `MAX`. | `SUM` |
| `missing_accepted` | Maximum missing fraction allowed within a condition before a feature is treated as too sparse for direct use. | `0.5` |
| `tune_sigma` | Multiplicative factor applied to the imputation standard deviation for minimum-probability sampling. | `1` |
| `tune_quantile` | Low quantile used to estimate the sample-specific imputation floor. | `0.0000001` |
| `pre_fix` | Prefix used to name cached `.rds` files and result tables across the entire pipeline. | `terminer_analysis_01` |
| `interesting_specificity` | Specificity classes considered “in scope” for PCA and downstream feature filtering. | `[semi_Nterm, semi_Cterm]` |
| `interesting_nterm_modif` | N-terminal modification classes considered "in scope" for PCA and downstream filtering. Quote `"n"` in YAML so it is not parsed as a boolean. | `["n", Acetyl]` |

**Available organisms**: `human`, `mouse`, `arabidopsis`, `medicago_truncatula`, `rhizobium_meliloti`, `pig`, `human_iso`, `ecoli`

### Stage parameters

#### `stages.data_preparation`

| Parameter | Description | Typical value |
| --- | --- | --- |
| `plot_width` | Default figure width for plots generated in the data-preparation notebook. | `10` |
| `plot_height` | Default figure height for plots generated in the data-preparation notebook. | `10` |
| `side_by_side_scale` | Scale factor used when arranging multiple panels side by side. | `2.0` |
| `plot_dpi` | Plot export DPI when high-resolution figures are rendered. | `300` |
| `pca_use_interest_only` | If `true`, PCA uses only features flagged by `interesting_specificity` and `interesting_nterm_modif`. | `false` |
| `do_batch_correction` | If `true`, run batch correction with ComBat where supported in the notebook. | `false` |
| `run_imputation_diagnostics` | If `true`, render cached imputation diagnostic plots after loading saved outputs. | `false` |
| `imputation_tsv_max_mb` | Maximum TSV size allowed for reloading cached imputation summaries when diagnostics are enabled. | `50` |

#### `stages.exploratory`

These parameters have the same meaning as in `stages.data_preparation`, but apply to `terminer_exploratory_analysis.qmd`.

| Parameter | Description | Typical value |
| --- | --- | --- |
| `plot_width` | Default figure width for exploratory plots. | `10` |
| `plot_height` | Default figure height for exploratory plots. | `10` |
| `side_by_side_scale` | Scale factor for side-by-side plot layouts. | `2.0` |
| `plot_dpi` | Plot export DPI. | `300` |
| `pca_use_interest_only` | Restrict PCA to the configured “interesting” feature scope. | `false` |
| `do_batch_correction` | Enable ComBat batch correction. | `false` |
| `run_imputation_diagnostics` | Enable reloading and plotting cached imputation diagnostics. | `false` |
| `imputation_tsv_max_mb` | Maximum TSV size allowed for cached imputation diagnostics. | `50` |

#### `stages.inferential`

| Parameter | Description | Typical value |
| --- | --- | --- |
| `plot_width` | Default figure width for inferential plots. | `10` |
| `plot_height` | Default figure height for inferential plots. | `8` |
| `fc_threshold` | Fold-change threshold used to classify features as upregulated or downregulated. | `1.2` |
| `pval_threshold` | Adjusted p-value threshold used for significance calls. | `0.05` |
| `run_go_analysis` | If `true`, build GO enrichment and GroupGO-derived GO compartment tables from the saved limma results. If `false`, the inferential notebook writes limma outputs only. | `false` |
| `use_blocking` | If `true`, limma runs duplicate-correlation blocking for repeated-measures designs. If `false`, limma runs without a blocking factor. | `false` |
| `blocking_column` | Experimental-design column to use when `use_blocking` is enabled. The loader standardizes `replicate` and `bio_replicate` so either can be supplied in the annotation file. | `replicate` |
| `interesting_processing_type` | Processing types to emphasize in downstream interpretation when applicable. | `[not_canonical]` |
| `defined_contrasts` | Named list of limma contrast formulas. Keys become contrast names in output files. | `Tumor_v_NAT: Tumor - NAT` |
| `datasets_to_analyze` | Which processed datasets should be passed into limma and enrichment analysis. Supported values include `protnorm_pure`, `protnorm_with_lost`, `non_protnorm`, and `protein_abund_norm`. | `[protnorm_with_lost, non_protnorm, protein_abund_norm]` |

#### `stages.results_visualization`

| Parameter | Description | Typical value |
| --- | --- | --- |
| `plot_width` | Default figure width for visualization outputs. | `12` |
| `plot_height` | Default figure height for visualization outputs. | `8` |
| `fc_threshold` | Fold-change threshold used to color volcano plots and summarize regulation status. | `1.2` |
| `pval_threshold` | Adjusted p-value threshold used in volcano plots and result summaries. | `0.05` |
| `dataset_of_interest` | Datasets from the limma results to display in visualization outputs. | `[protnorm_with_lost, non_protnorm, protein_abund_norm]` |
| `condition_levels` | Ordered factor levels for conditions in plots and summaries. | `[NAT, PC, Tumor]` |
| `contrast_levels` | Ordered factor levels for contrasts in plots and summaries. These should match the names in `defined_contrasts`. | `[Tumor_v_NAT, PC_v_NAT, Tumor_v_PC]` |
| `run_sequence_logos` | If `true`, build sequence-logo visualizations with `dagLogo`. | `false` |
| `sequence_logo_contrasts` | Optional subset of contrasts to include in sequence logos. Use `null` to include all configured contrasts. | `null` |
| `sequence_logo_nterm_modifs` | Optional subset of N-terminal modifications to include in sequence logos. Use `null` to include all available modifications. | `null` |

### Search-engine examples

#### DIA-NN configuration
```yaml
common:
  search_engine: diann
  diann_report_location: data/report.parquet
  proteotypic_only: true
  summarization_method: SUM
```

#### FragPipe TMT configuration
```yaml
common:
  search_engine: fragpipe_tmt
  fragpipe_parent_dir: data/fragpipe_search
  ref_sample: norm
  min_purity: 0.5
  tmt_delta: "229"
```

#### FragPipe label-free configuration
```yaml
common:
  search_engine: fragpipe_lf
  fragpipe_lf_parent_dir: data/fragpipe_lf_search
  fragpipe_lf_annotation: data/fragpipe_lf_annotation.txt
```

#### FragPipe heavy-light configuration
```yaml
common:
  search_engine: fragpipe_heavy_light
  fragpipe_hl_file: data/combined_modified_peptide_label_quant.tsv
  fragpipe_hl_annotation: data/fragpipe_hl_annotation.txt
```

#### Spectronaut configuration
```yaml
common:
  search_engine: spectronaut
  spectronaut_report: data/spectronaut_report.tsv
  proteotypic_only: true
```

## Running the analysis

### 1. Setup project structure
Create the following directory structure:
```
project_root/
├── data/
│   ├── report.parquet (or report.tsv)
│   ├── experimental_annotation.txt
│   ├── proteome.fasta
│   └── example_annotation.txt (optional)
├── scr/                    # Helper scripts
│   ├── modif_diann_adapter.R
│   └── helper_functions.R
├── rds/                    # Cached intermediate results
├── results/                # Final output tables
├── terminer_exploratory_analysis.qmd
├── terminer_inferential_analysis.qmd
└── terminer_results_visualization.qmd
```

**Note**: Use `example_annotation.txt` as a template for creating your `data/experimental_annotation.txt` file.

### 2. Data preparation step

1. Open `terminer_data_preparation.qmd`.
2. Update `config.yaml` with your project-specific paths and analysis settings.
3. Execute the script to perform
   - Data loading and formating.
   - Missing value analysis.
   - Missing value imputation if necessary.
   - Store intermediary results for further analysis.

### 2. Exploratory analysis
1. Open `terminer_exploratory_analysis.qmd`
2. Confirm `config.yaml` contains the correct exploratory settings
3. Execute the script to perform:
   - Data loading and quality control
   - Peptide annotation with TermineR
   - Missing value analysis
   - Data imputation
   - Principal component analysis
   - Results caching for downstream analysis

### 3. Inferential analysis
1. Open `terminer_inferential_analysis.qmd`
2. Define your experimental contrasts in `config.yaml` under `stages.inferential.defined_contrasts`
3. Execute the script to perform:
   - Protein-level normalization
   - Differential abundance analysis with limma
  - Optional GO enrichment analysis
   - Statistical results generation

### 4. Results visualization
1. Open `terminer_results_visualization.qmd`
2. Execute the script to generate:
   - Volcano plots
   - Heatmaps of regulated features
   - Summary statistics
   - Publication-ready figures

## Output files

### Results directory (`results/`)
- `*_differential_analysis_results.tsv`: Complete differential analysis results
- `*_go_enrichment_results.tsv`: Optional GO enrichment analysis results
- `*_go_cc_mapping.tsv`: Optional gene-to-GO compartment summary derived from GroupGO
- `*_summary_statistics.tsv`: Summary of analysis results

### RDS cache directory (`rds/`)
- Intermediate results cached for faster re-analysis
- Can be safely deleted to force re-computation

## Key features

### Computational efficiency
- **RDS Caching**: Computationally intensive steps are cached as RDS files
- **Modular Design**: Each analysis step can be run independently
- **Memory Management**: Automatic cleanup of large intermediate objects

### Standardized analysis
- **Consistent Chunk Naming**: All code chunks follow standardized naming conventions
- **Plot Sizing**: Standardized plot dimensions and DPI settings
- **Parameter Configuration**: Centralized parameter definition

### Quality control
- **Missing Value Analysis**: Comprehensive assessment of data completeness
- **PCA Analysis**: Quality control through principal component analysis
- **Summary Statistics**: Detailed reporting of analysis results

## Data loading features

### Automatic file format detection
The pipeline automatically detects whether your DIA-NN report is in `.parquet` or `.tsv` format:

- **Parquet files**: Uses the optimized `diann_adapter_parquet()` function with:
  - Efficient memory usage through Arrow/parquet format
  - MAD (Median Absolute Deviation) scaling for normalization
  - Configurable summarization methods (SUM or MAX)
  - Built-in modification annotation

- **TSV files**: Uses the standard `diann_adapter()` function from the diann R package

### Configuration for parquet files
When using parquet files, you can configure:
```r
# In the data loading section
df_from_diann <- diann_adapter_parquet(
  path_to_file = diann_report_location,
  proteotypic = TRUE,        # Use only proteotypic peptides
  summarization = "SUM"      # Options: "SUM" or "MAX"
)
```

## Customization

### Adding new organisms
To add support for a new organism:
1. Ensure the organism is supported by TermineR
2. Install the appropriate Bioconductor annotation package
3. Update the `organism_annotation` parameter
4. Modify the `library()` call for the organism database in the inferential script

You can always open an issue requesting your organism of interest, if you want to include the Uniprot Processing annotation in your analysis.

### Modifying analysis parameters
Key parameters can be adjusted in the parameter definition sections:
- **Missing value tolerance**: Adjust `missing_accepted`
- **Statistical thresholds**: Modify `fc_threshold` and `pval_threshold`
- **Imputation strategy**: Customize imputation parameters in the exploratory script

### Custom contrasts
Define custom contrasts in the `defined_contrasts` vector using limma syntax:
```r
defined_contrasts <- c(
  "Treatment_vs_Control" = "Treatment - Control",
  "HighDose_vs_LowDose" = "HighDose - LowDose"
)
```

## Troubleshooting

### Common issues
1. **Memory errors**: Increase memory limits or process data in smaller chunks
2. **Missing annotation**: Ensure all required input files are present and properly formatted
3. **Contrast errors**: Verify that condition names in contrasts match those in the annotation file

### Troubleshooting by search engine

#### DIA-NN issues
- **Column naming**: Check that sample names match experimental annotation
- **Memory issues**: Make sure you have enough memory allocated for large datasets

#### FragPipe TMT issues
- **Missing annotation files**: Ensure each mixture directory has annotation.txt
- **Reference channel**: Verify reference channel name matches annotation
- **TMT delta**: Use correct mass delta for your TMT kit

#### FragPipe label-free issues
- **Run mapping**: Ensure run names in annotation match spectrum file names

#### FragPipe heavy-light issues
- **File format**: Ensure combined_modified_peptide_label_quant.tsv exists
- **Label detection**: Verify heavy/light modifications are correctly detected

### Getting help
- Check the TermineR documentation: [GitHub repository](https://github.com/MiguelCos/TermineR)
- Ensure all required packages are installed and up to date
- Verify input file formats match the specifications above
- For specific search engine issues, refer to their respective documentation or forums
- Report an issue on ther TermineR GitHub repository if you encounter problems not covered here.

## Citation

If you use this pipeline in your research, please cite:
- The TermineR paper: https://doi.org/10.1002/pmic.202300491
- The [TermineR repo](https://github.com/MiguelCos/TermineR)
- Relevant Bioconductor packages (limma, clusterProfiler, etc.)
- FragPipe and DIA-NN tools
