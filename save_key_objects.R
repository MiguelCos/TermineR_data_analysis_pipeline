# Script to save key objects for the visualization report
# Run this after your main analysis completes successfully

library(here)

# Create RDS directory
rds_dir <- here("r_project_files/rds")
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

message("Saving key objects for visualization report...")

# Save differential abundance results
if(exists("all_toptables_neo_termini")) {
  write_rds(all_toptables_neo_termini, file.path(rds_dir, "all_toptables_neo_termini.rds"))
  message("✓ Saved all_toptables_neo_termini")
}

if(exists("all_toptables_protnorm_neo_termini")) {
  write_rds(all_toptables_protnorm_neo_termini, file.path(rds_dir, "all_toptables_protnorm_neo_termini.rds"))
  message("✓ Saved all_toptables_protnorm_neo_termini")
}

if(exists("all_toptables_prot_quant")) {
  write_rds(all_toptables_prot_quant, file.path(rds_dir, "all_toptables_prot_quant.rds"))
  message("✓ Saved all_toptables_prot_quant")
}

# Save abundance matrices for PCA
if(exists("abundance_matrix_sum_imp_neo_termini")) {
  write_rds(abundance_matrix_sum_imp_neo_termini, file.path(rds_dir, "abundance_matrix_sum_imp_neo_termini.rds"))
  message("✓ Saved abundance_matrix_sum_imp_neo_termini")
}

if(exists("abundance_matrix_sum_protnorm_neo_termini")) {
  write_rds(abundance_matrix_sum_protnorm_neo_termini, file.path(rds_dir, "abundance_matrix_sum_protnorm_neo_termini.rds"))
  message("✓ Saved abundance_matrix_sum_protnorm_neo_termini")
}

if(exists("mat_prot_quant")) {
  write_rds(mat_prot_quant, file.path(rds_dir, "mat_prot_quant.rds"))
  message("✓ Saved mat_prot_quant")
}

# Save annotation data
if(exists("sample_annotation")) {
  write_rds(sample_annotation, file.path(rds_dir, "sample_annotation.rds"))
  message("✓ Saved sample_annotation")
}

if(exists("interesting_features_table")) {
  write_rds(interesting_features_table, file.path(rds_dir, "interesting_features_table.rds"))
  message("✓ Saved interesting_features_table")
}

# Save proportion data
if(exists("pept_summ_rawpur_semi_3")) {
  write_rds(pept_summ_rawpur_semi_3, file.path(rds_dir, "pept_summ_rawpur_semi_3.rds"))
  message("✓ Saved pept_summ_rawpur_semi_3")
}

# Save additional objects for enhanced visualization
if(exists("fasta_location")) {
  write_rds(fasta_location, file.path(rds_dir, "fasta_location.rds"))
  message("✓ Saved fasta_location")
}

# Save matrices for different PCA analyses
if(exists("abundance_matrix_sum_imp_semi_c")) {
  write_rds(abundance_matrix_sum_imp_semi_c, file.path(rds_dir, "abundance_matrix_sum_imp_semi_c.rds"))
  message("✓ Saved abundance_matrix_sum_imp_semi_c")
}

if(exists("abundance_matrix_sum_protnorm_semi_c")) {
  write_rds(abundance_matrix_sum_protnorm_semi_c, file.path(rds_dir, "abundance_matrix_sum_protnorm_semi_c.rds"))
  message("✓ Saved abundance_matrix_sum_protnorm_semi_c")
}

if(exists("abundance_matrix_sum_imp_semi_n")) {
  write_rds(abundance_matrix_sum_imp_semi_n, file.path(rds_dir, "abundance_matrix_sum_imp_semi_n.rds"))
  message("✓ Saved abundance_matrix_sum_imp_semi_n")
}

if(exists("abundance_matrix_sum_protnorm_semi_n")) {
  write_rds(abundance_matrix_sum_protnorm_semi_n, file.path(rds_dir, "abundance_matrix_sum_protnorm_semi_n.rds"))
  message("✓ Saved abundance_matrix_sum_protnorm_semi_n")
}

# Save volcano data with all annotations for cleavage analysis
if(exists("all_toptables_protnorm_neo_termini_volcano_data")) {
  write_rds(all_toptables_protnorm_neo_termini_volcano_data, file.path(rds_dir, "all_toptables_protnorm_neo_termini_volcano_data.rds"))
  message("✓ Saved all_toptables_protnorm_neo_termini_volcano_data")
}

if(exists("all_toptables_neo_termini_volcano_data")) {
  write_rds(all_toptables_neo_termini_volcano_data, file.path(rds_dir, "all_toptables_neo_termini_volcano_data.rds"))
  message("✓ Saved all_toptables_neo_termini_volcano_data")
}

message("All key objects saved successfully!")
message("You can now render the visualization report: OS001_TermineR_results_visualization_report.qmd")
