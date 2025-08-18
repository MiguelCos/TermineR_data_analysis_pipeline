# TermineR Data Analysis Pipeline

This repository contains a streamlined workflow for the analysis of shotgun proteomics data from multiple search engines, specifically designed for the analysis of proteolytic processing events using the TermineR approach.

## Overview

The pipeline consists of three main analysis scripts that support multiple search engine outputs:

1. **Exploratory Analysis**: Data loading, quality control, and initial visualization
2. **Inferential Analysis**: Statistical analysis, protein normalization, and GO enrichment  
3. **Results Visualization**: Generation of publication-ready plots and summary statistics

## Supported Search Engines

The pipeline supports the following search engine outputs through integrated TermineR adapters:

### 1. DIA-NN (Data-Independent Acquisition)
- **File types**: `report.parquet` or `report.tsv`
- **Quantification**: Precursor-level quantification with MAD scaling
- **Modifications**: Automatic N-terminal modification detection

### 2. FragPipe TMT (Tandem Mass Tags)
- **File types**: `psm.tsv` files in mixture directories
- **Quantification**: TMT reporter ion intensities with reference channel normalization
- **Modifications**: N-terminal TMT, Acetyl, Dimethyl detection

### 3. FragPipe Label-Free
- **File types**: `psm.tsv` files with intensity-based quantification
- **Quantification**: Precursor intensity with MAD scaling
- **Modifications**: N-terminal modification detection

### 4. FragPipe Heavy-Light (SILAC/Dimethyl)
- **File types**: `combined_modified_peptide_label_quant.tsv`
- **Quantification**: Heavy/Light label quantification
- **Modifications**: Dimethyl and Acetyl detection

### 5. Spectronaut
- **File types**: Spectronaut report TSV files
- **Quantification**: Intensity-based with log2 transformation
- **Modifications**: Acetyl detection (extensible)

## Requirements

### R Packages

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

## Input Data Requirements

### 1. DIA-NN Report File
- **Location**: `data/report.parquet` or `data/report.tsv`
- **Description**: Output from FragPipe-DIA + DIANN quantitation
- **Format**: Parquet or TSV file containing peptide identifications and quantitative data
- **Note**: The pipeline automatically detects the file format and uses the appropriate loader:
  - `.parquet` files: Uses optimized parquet adapter with MAD scaling
  - `.tsv` files: Uses standard DIA-NN adapter

### 2. Experimental Annotation File
- **Location**: `data/experimental_annotation.txt`
- **Description**: Tab-delimited file describing the experimental design and sample metadata
- **Format**: Tab-separated values (.txt or .tsv)

#### Required Columns
- **`sample`**: Sample identifier that must exactly match the column names in your search engine output files
- **`sample_name`**: Human-readable sample name or identifier  
- **`condition`**: Experimental condition or treatment group (e.g., "Control", "Treatment", "TU", "NAT", "PC")
- **`bio_replicate`**: Biological replicate identifier
  - **Important**: Use the same `bio_replicate` number when samples represent repeated measures of the same biological entity (e.g., technical replicates, fractionation, or multiple time points from the same individual)
  - Use different `bio_replicate` numbers for independent biological samples
- Additional metadata columns can be included as needed

#### Annotation File Guidelines
1. **Header required**: First row must contain column names
2. **Tab-delimited**: Use tab characters to separate columns
3. **No missing values**: All required columns must have values for each sample
4. **Consistent naming**: Sample identifiers must match exactly between annotation and data files
5. **Biological replicates**: 
   - Same `bio_replicate` = repeated measures from same biological source
   - Different `bio_replicate` = independent biological samples
   - This is critical for proper statistical modeling in limma

#### Example Annotation Files
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

### 3. Protein FASTA File
- **Location**: `data/proteome.fasta`
- **Description**: Protein sequence database used for the search
- **Format**: Standard FASTA format

### 4. TargetP Results (Optional)
- **Location**: `data/targetp_results.targetp2`
- **Description**: TargetP2 prediction results for subcellular localization
- **Format**: TargetP2 output format
- **Note**: Set `targetp_location <- NULL` if not available

## Configuration Parameters

### Search Engine Selection
```r
search_engine <- "diann"  # Options: "diann", "fragpipe_tmt", "fragpipe_lf", "fragpipe_heavy_light", "spectronaut"
```

### Data Locations (Configure based on search engine)

#### DIA-NN Configuration
```r
search_engine <- "diann"
diann_report_location <- here("initial_data/report.parquet")  # or .tsv
proteotypic_only <- TRUE
summarization_method <- "SUM"  # or "MAX"
```

#### FragPipe TMT Configuration
```r
search_engine <- "fragpipe_tmt"
fragpipe_parent_dir <- here("initial_data/fragpipe_search")  # Directory with mix_1, mix_2, etc.
ref_sample <- "norm"  # Reference channel name
min_purity <- 0.5
tmt_delta <- "229"  # "229" for TMT10/11, "304" for TMT16
```

#### FragPipe Label-Free Configuration
```r
search_engine <- "fragpipe_lf"
fragpipe_lf_parent_dir <- here("initial_data/fragpipe_lf_search")
fragpipe_lf_annotation <- here("initial_data/fragpipe_lf_annotation.txt")
```

#### FragPipe Heavy-Light Configuration
```r
search_engine <- "fragpipe_heavy_light"
fragpipe_hl_file <- here("initial_data/combined_modified_peptide_label_quant.tsv")
fragpipe_hl_annotation <- here("initial_data/fragpipe_hl_annotation.txt")
```

#### Spectronaut Configuration
```r
search_engine <- "spectronaut"
spectronaut_report <- here("initial_data/spectronaut_report.tsv")
proteotypic_only <- TRUE
```

### Experimental Parameters
```r
sense_protease <- "C"  # "C" for C-terminal cleavage (trypsin), "N" for N-terminal
specificity_protease <- "K|R"  # Protease specificity (trypsin: "K|R")
organism_annotation <- "human"  # Organism for annotation
instrument <- "EX"  # Instrument prefix in sample names
```

**Available organisms**: `"human"`, `"mouse"`, `"arabidopsis"`, `"medicago_truncatula"`, `"rhizobium_meliloti"`, `"pig"`, `"human_iso"`, `"ecoli"`

### Analysis Parameters
```r
missing_accepted <- 2 / 4  # Maximum missing values per condition (2 out of 4 replicates)
fc_threshold <- 2.5  # Fold-change threshold for significance
pval_threshold <- 0.05  # P-value threshold for significance
pre_fix <- "terminer_analysis_"  # Prefix for output files
```

### Contrast Definition
```r
defined_contrasts <- c(
  "Treatment_vs_Control" = "Treatment - Control",
  "Condition2_vs_Control" = "Condition2 - Control"
  # Add more contrasts as needed
)
```

## Running the Analysis

### 1. Setup Project Structure
Create the following directory structure:
```
project_root/
├── data/
│   ├── report.parquet (or report.tsv)
│   ├── experimental_annotation.txt
│   ├── proteome.fasta
│   ├── targetp_results.targetp2 (optional)
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

### 2. Exploratory Analysis
1. Open `terminer_exploratory_analysis.qmd`
2. Update the parameter section with your specific settings
3. Execute the script to perform:
   - Data loading and quality control
   - Peptide annotation with TermineR
   - Missing value analysis
   - Data imputation
   - Principal component analysis
   - Results caching for downstream analysis

### 3. Inferential Analysis
1. Open `terminer_inferential_analysis.qmd`
2. Define your experimental contrasts
3. Execute the script to perform:
   - Protein-level normalization
   - Differential abundance analysis with limma
   - GO enrichment analysis
   - Statistical results generation

### 4. Results Visualization
1. Open `terminer_results_visualization.qmd`
2. Execute the script to generate:
   - Volcano plots
   - Heatmaps of regulated features
   - Summary statistics
   - Publication-ready figures

## Output Files

### Results Directory (`results/`)
- `*_differential_analysis_results.tsv`: Complete differential analysis results
- `*_go_enrichment_results.tsv`: GO enrichment analysis results
- `*_summary_statistics.tsv`: Summary of analysis results

### RDS Cache Directory (`rds/`)
- Intermediate results cached for faster re-analysis
- Can be safely deleted to force re-computation

## Key Features

### Computational Efficiency
- **RDS Caching**: Computationally intensive steps are cached as RDS files
- **Modular Design**: Each analysis step can be run independently
- **Memory Management**: Automatic cleanup of large intermediate objects

### Standardized Analysis
- **Consistent Chunk Naming**: All code chunks follow standardized naming conventions
- **Plot Sizing**: Standardized plot dimensions and DPI settings
- **Parameter Configuration**: Centralized parameter definition

### Quality Control
- **Missing Value Analysis**: Comprehensive assessment of data completeness
- **PCA Analysis**: Quality control through principal component analysis
- **Summary Statistics**: Detailed reporting of analysis results

## Data Loading Features

### Automatic File Format Detection
The pipeline automatically detects whether your DIA-NN report is in `.parquet` or `.tsv` format:

- **Parquet files**: Uses the optimized `diann_adapter_parquet()` function with:
  - Efficient memory usage through Arrow/parquet format
  - MAD (Median Absolute Deviation) scaling for normalization
  - Configurable summarization methods (SUM or MAX)
  - Built-in modification annotation

- **TSV files**: Uses the standard `diann_adapter()` function from the diann R package

### Configuration for Parquet Files
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

### Adding New Organisms
To add support for a new organism:
1. Ensure the organism is supported by TermineR
2. Install the appropriate Bioconductor annotation package
3. Update the `organism_annotation` parameter
4. Modify the `library()` call for the organism database in the inferential script

### Modifying Analysis Parameters
Key parameters can be adjusted in the parameter definition sections:
- **Missing value tolerance**: Adjust `missing_accepted`
- **Statistical thresholds**: Modify `fc_threshold` and `pval_threshold`
- **Imputation strategy**: Customize imputation parameters in the exploratory script

### Custom Contrasts
Define custom contrasts in the `defined_contrasts` vector using limma syntax:
```r
defined_contrasts <- c(
  "Treatment_vs_Control" = "Treatment - Control",
  "HighDose_vs_LowDose" = "HighDose - LowDose"
)
```

## Troubleshooting

### Common Issues
1. **Memory errors**: Increase memory limits or process data in smaller chunks
2. **Missing annotation**: Ensure all required input files are present and properly formatted
3. **Contrast errors**: Verify that condition names in contrasts match those in the annotation file

### Troubleshooting by Search Engine

#### DIA-NN Issues
- **Column naming**: Check that sample names match experimental annotation
- **Memory issues**: Make sure you have enough memory allocated for large datasets

#### FragPipe TMT Issues
- **Missing annotation files**: Ensure each mixture directory has annotation.txt
- **Reference channel**: Verify reference channel name matches annotation
- **TMT delta**: Use correct mass delta for your TMT kit

#### FragPipe Label-Free Issues
- **Run mapping**: Ensure run names in annotation match spectrum file names

#### FragPipe Heavy-Light Issues
- **File format**: Ensure combined_modified_peptide_label_quant.tsv exists
- **Label detection**: Verify heavy/light modifications are correctly detected

### Getting Help
- Check the TermineR documentation: [GitHub repository](https://github.com/MiguelCos/TermineR)
- Ensure all required packages are installed and up to date
- Verify input file formats match the specifications above
- For specific search engine issues, refer to their respective documentation or forums
- Report an issue on ther TermineR GitHub repository if you encounter problems not covered here.

## Citation

If you use this pipeline in your research, please cite:
- The TermineR package and methodology
- Relevant Bioconductor packages (limma, clusterProfiler, etc.)
- FragPipe and DIA-NN tools
