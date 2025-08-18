diann_adapter_parquet <- function(
    path_to_file,
    proteotypic = TRUE,
    summarization = "SUM" # options: "SUM" or "MAX"
    ) {

  require(diann)
  require(tibble)
  require(dplyr)
  require(tidyr)
  require(stringr)
  require(data.table)
  require(arrow)  # ADD THIS
  require(purrr)  # ADD THIS for map_chr
  
  
  data("unimod_id_to_name_mapping", package = "TermineR")
  
  # the following function modifications are based on the diann R package code hosted on GitHub
  # https://github.com/vdemichev/diann-rpackage/blob/master/R/diann-R.R
  
  # pivot_aggregate function modified to sum up instead of get max value
  pivot_aggregate_2 <- function(
    df, 
    sample.header, 
    id.header, 
    quantity.header) {
    
    x <- melt.data.table(
      df, 
      id.vars = c(
        sample.header, 
        id.header), 
      measure.vars = c(
        quantity.header))
    
    x$value[which(x$value == 0)] <- NA
    
    piv <- as.data.frame(
      dcast.data.table(x, 
                       as.formula(paste0(id.header,'~',sample.header)), 
                       value.var = "value", 
                       fun.aggregate = function(x) sum(x, na.rm=TRUE))) 
    
    rownames(piv) <- piv[[1]]
    
    piv[[1]] <- NULL
    
    piv <- piv[order(rownames(piv)),]
    
    piv = as.matrix(piv)
    
    piv[is.infinite(piv)] <- NA
    
    # check if value is 0
    piv[piv == 0] <- NA
    
    piv
  }

  pivot_2 <- function(df, sample.header, id.header, quantity.header) {
  x <- melt.data.table(df, id.vars = c(sample.header, id.header), measure.vars = c(quantity.header))
  x$value[which(x$value == 0)] <- NA
  piv <- as.data.frame(dcast.data.table(x, as.formula(paste0(id.header,'~',sample.header)), value.var = "value")) 
  rownames(piv) <- piv[[1]]
  piv[[1]] <- NULL
  piv <- piv[order(rownames(piv)),]
  as.matrix(piv)
}
  
  # modify diann_matrix function to sum up instead of calculating max
  diann_matrix_2 <- function(
    x, 
    id.header = "Precursor.Id", 
    quantity.header = "Precursor.Normalised", 
    proteotypic.only = F, 
    q = 0.01, 
    protein.q = 1.0, 
    pg.q = 1.0, 
    gg.q = 1.0) {
    
    df <- as.data.table(x)
    
    if (proteotypic.only) df <- df[which(df[['Proteotypic']] != 0),]
    
    df <- unique(df[which(df[[id.header]] != "" & df[[quantity.header]] > 0 & df[['Q.Value']] <= q & df[['Protein.Q.Value']] <= protein.q & df[['PG.Q.Value']] <= pg.q & df[['GG.Q.Value']] <= gg.q),c("File.Name", id.header, quantity.header),with=FALSE])
    
    is_duplicated = any(duplicated(paste0(df[["File.Name"]],":",df[[id.header]])))
    
    if (is_duplicated) {
      
      warning("Multiple quantities per id: the sum of these will be calculated")
      pivot_aggregate_2(df, "File.Name", id.header, quantity.header)
      
    } else {
      
      pivot_2(df, "File.Name", id.header, quantity.header)
      
    }
  }

# Validate required columns exist
required_cols <- c("Precursor.Id", "Protein.Ids", "Precursor.Normalised", 
                   "Modified.Sequence", "Stripped.Sequence", "Proteotypic",
                   "Q.Value", "Protein.Q.Value", "PG.Q.Value", "GG.Q.Value", "Run")

dataset_schema <- arrow::open_dataset(path_to_file)$schema
  available_cols <- names(dataset_schema)
  
missing_cols <- setdiff(required_cols, available_cols)
if(length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

diann_df_min <- arrow::open_dataset(path_to_file) %>%
  dplyr::select(
    # Add only the columns you actually need for analysis
    Precursor.Id, 
    Protein.Ids,
    Precursor.Normalised,
    Modified.Sequence,
    Stripped.Sequence,
    Proteotypic,
    Q.Value,
    Protein.Q.Value,
    PG.Q.Value,
    GG.Q.Value,
    File.Name = Run
  ) %>%
  collect() %>%
    # Only select columns we actually need
    # extract modifications from modified_sequence check for everything inside '()'
    mutate(
      first_modif = str_extract(Modified.Sequence, "\\((.*?)\\)") %>% 
                                str_remove_all("[\\(\\)]")  # REPLACE map_chr line
    ) %>%
    # substitute the modification, including the brackets, with a Z
    mutate(
      min_first_mod_seq = str_replace(Modified.Sequence, "\\((.*?)\\)", "Z")
    ) %>%
    # get the location of those modifications
    mutate(
      # get the start of the first modification only using str_locate
      first_modif_locat = str_locate(min_first_mod_seq, "Z")[, "start"],
      id_nr = parse_number(first_modif)) %>%
    left_join(.,
              unimod_id_to_name_mapping,
              by = "id_nr") %>%  # ADD by parameter
    mutate(
      nterm_modif = case_when(
        is.na(name) ~ "n",
        first_modif_locat == 1 & !is.na(name) ~ name,
        first_modif_locat != 1 & !is.na(name) ~ "n"
      )
    ) %>%
    mutate(
      nterm_modif = case_when(
        str_detect(nterm_modif, "Dimethyl") ~ "Dimethyl",
        str_detect(nterm_modif, "TMT") ~ "TMT",
        TRUE ~ nterm_modif
      )
    ) %>% 
    mutate(
      nterm_modif_peptide = paste(
        nterm_modif,
        Stripped.Sequence,
        sep = "_"),
    ) %>%
    relocate(
      nterm_modif_peptide,
      nterm_modif,
      .before = Modified.Sequence
    )
  
  peptide2protein <- diann_df_min %>%
    dplyr::select(
      protein = Protein.Ids,
      peptide = Stripped.Sequence,
      nterm_modif_peptide,
      nterm_modif
    ) %>%
    mutate(
      protein = str_remove(
        protein,
        " .*")
    ) %>%
    dplyr::select(
      protein,
      peptide,
      nterm_modif_peptide,
      nterm_modif
    ) %>%
    distinct()
  
  # evaluate if summarization == "SUM" or "MAX" and choose the right diann_amtrix function
  
  if(summarization == "SUM"){
  
    diann_df_mat_min <- diann_matrix_2(
      diann_df_min,
      id.header = "nterm_modif_peptide",
      quantity.header = "Precursor.Normalised",
      proteotypic.only = proteotypic,
      q = 0.01)
  
  } else if(summarization == "MAX"){
  
    diann_df_mat_min <- diann_matrix(
      diann_df_min,
      id.header = "nterm_modif_peptide",
      quantity.header = "Precursor.Normalised",
      proteotypic.only = proteotypic,
      q = 0.01)
  
  }

  colnames(diann_df_mat_min) <- str_remove(
    colnames(diann_df_mat_min),
    ".*\\\\") %>%
    str_remove(
      "\\..*")

  diann_df_mat_min <- diann_df_mat_min %>%
    as.data.frame() %>%
    rownames_to_column("nterm_modif_peptide") %>%
    arrange(nterm_modif_peptide) %>%
    # apply log2 transformation BEFORE pivot_longer
    mutate(
      across(
        where(is.double),
        ~log2(.x)
      )
    ) %>%
    # transform to long format AFTER log2 transformation
    pivot_longer(cols = where(is.double),
                 names_to = "sample",
                 values_to = "log2_ratio_2_ref")
      
# apply MAD scaling
  scaled_diann_df_mat_min <- diann_df_mat_min %>%
    # Rename the column for consistency with downstream code
    dplyr::rename(log2_rat2ref_group = log2_ratio_2_ref) %>%

    # calculate median ratios per sample
    group_by(sample) %>%
    mutate(Mi = median(log2_rat2ref_group,
                       na.rm = TRUE)) %>%
    ungroup() %>%
    
    # global median across all samples
    mutate(M0 = median(Mi,
                       na.rm = TRUE)) %>%
    
    # median centered ratios based on median ratios per sample
    mutate(RCij = log2_rat2ref_group - Mi) %>%
    
    # calculate median absolute deviation per sample
    group_by(sample) %>%
    mutate(MADi = median(abs(RCij),
                         na.rm = TRUE)) %>%
    ungroup() %>%
    
    # calculate global median absolute deviation
    mutate(MAD0 = median(MADi,
                         na.rm = TRUE)) %>%
    ungroup() %>%
    
    # calculate scaled ratios
    mutate(RNij = (RCij / MADi) * MAD0 + M0)
  
  # turn back to wide format using RNij as abundance
  scaled_diann_df_mat_min <- scaled_diann_df_mat_min %>%
    dplyr::select(
      nterm_modif_peptide,
      sample,
      RNij
    ) %>%
    pivot_wider(
      names_from = "sample",
      values_from = "RNij"
    ) %>%
    left_join(
      peptide2protein,
      by = "nterm_modif_peptide",
      multiple = "first"
    ) %>%
    relocate(
      nterm_modif_peptide,
      nterm_modif,
      peptide,
      protein
    ) %>%
    arrange(peptide)

  return(scaled_diann_df_mat_min)
}