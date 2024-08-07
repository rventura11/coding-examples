library(reportGeneration)
library(tidyverse)

# Not sure if we need this for review, but tossing it in
get_ext_patient <- function(accession, SF) {
  
  res <- reportGeneration::combined_get_order(SF, accession)
  Project <- res$project$name
  Patient <- res$patient$external_pid
  Accession <- res$requisition_id
  
  ord <- tibble::tibble(Project, Patient, Accession)
}

run_dir <- "AREX116627"

setwd <- (file.path("/shared_efs/processed_runs", run_dir))

outdir <- "review_csvs"
system(paste0('mkdir ', outdir))

sample_table <- readr::read_csv("radar_qc_loader_files/sample_table.csv",
                                col_types = cols("sampleID" = col_character()))

sample_call <- readr::read_csv("radar_qc_loader_files/sample_call.csv",
                               col_types = cols("sample_name" = col_character()))

sample_qc <- readr::read_csv("radar_qc_loader_files/sample_QC_table.csv",
                             col_types = cols("sampleID" = col_character()))

accs <- dplyr::pull(sample_table, accessionID)

x <- reportGeneration::init_LIS_SF()
x <- reportGeneration::login_LIS_SF(x)

sf_info <- purrr::map_dfr(.x = accs, SF = x, .f = get_ext_patient) %>% 
  unique()

project_summary <- sample_table %>%
dplyr::mutate(type = ifelse(grepl("plasma", type), "plasma", type)) %>%
dplyr::group_by(projectID, type) %>%
dplyr::summarise(total = n()) %>%
tidyr::pivot_wider(names_from = type, values_from = total) %>%
dplyr::mutate(total_samples = rowSums(across(where(is.numeric))))

print(project_summary)
#readr::write_csv(project_summary, file.path(outdir, "project_summary.csv"))

plasma_summary <- sample_table %>% 
  dplyr::mutate(type = ifelse(grepl("plasmaStreck", type), "plasma", type)) %>% 
  dplyr::filter(type == "plasma") %>% 
  dplyr::select(projectID, accessionID, sampleID) %>% 
  dplyr::arrange(projectID, accessionID)

#readr::write_csv(plasma_summary, file.path(outdir, "plasma_summary.csv"))

qc_table <- sample_table %>% 
  dplyr::select(projectID, accessionID, sampleID, panel) %>% 
  dplyr::inner_join(select(sample_qc, sampleID, panel,
                           QC.flag,
                           errorRate.flag, 
                           interbarcodeSNPconsistency.flag, 
                           intersampleSNPconsistency.flag,
                           contamination.flag,gender.flag,
                           depth.flag), by = c("sampleID", "panel"))

readr::write_csv(qc_table, file.path(outdir, "qc_table.csv"))

qc_summary <- qc_table %>% 
  tidyr::pivot_longer(cols = ends_with("flag"),
                      names_to = "QC_Metric",
                      values_to = "QC_Flag") %>% 
  dplyr::group_by(QC_Metric, QC_Flag) %>% 
  dplyr::summarise(Total = n()) %>% 
  dplyr::arrange(desc(QC_Metric), desc(QC_Flag))

readr::write_csv(qc_summary, file.path(outdir, "qc_summary.csv"))

all_pass_summary <- sample_table %>% 
  dplyr::select(projectID, accessionID, sampleID, panel) %>% 
  dplyr::inner_join(select(sample_call, sample_name, supporting_var_n),
                    by = c("sampleID" = "sample_name")) %>% 
  dplyr::arrange(supporting_var_n)

readr::write_csv(all_pass_summary, file.path(outdir, "all_pass_summary.csv"))