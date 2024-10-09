library(reportGeneration)
library(tidyverse)
library(glue)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("-a", "--arex", type = "character", required = TRUE,
                    help = "Analysis to generate project reports for")

args <- parser$parse_args()

exome_run <- args$arex

setwd(file.path("/shared/pipeline-user/run_data", exome_run))

get_exome_run <- function(exome_run) {
  ex_q <- glue('SELECT * FROM exome_analysis WHERE vcf_file_path like "%{exome_run}%";')
  
  q_res <- reportGeneration::query.database(ex_q)
  
}

get_sf_info <- function(accession, SF) {
  
  res <- reportGeneration::combined_get_order(SF, accession)
  Accession <- accession
  Patient <- res$patient$external_pid

  o.df <- tibble::tibble(Accession, Patient)
}

get_panel_info <- function(panel_file) {
  
  ssummary <- readr::read_csv(panel_file) %>% 
    # dplyr::mutate(panel = gsub("_.*", "", amplicon_name)) %>% 
    dplyr::mutate(Panel = gsub('_.*','',basename(panel_file))) %>% 
    dplyr::group_by(Panel) %>% 
    dplyr::summarise(Total = n(),
                     SNVs = sum(grepl("mis", variant_code)),
                     Indels = sum(!grepl("mis", variant_code)),
                     Priority = sum(priority),
                     Priority_SNVs = sum(priority & grepl("mis", variant_code)),
                     Priority_Indels = sum(priority & !grepl("mis", variant_code)),
                     Total_oligos = 4 * Total)
  # Nonpriority_SNVs = sum(!priority & grepl("mis", variant_code)),
  # Nonpriority_Indels = sum(!priority & !grepl("mis", variant_code)))
}

exome_results <- get_exome_run(exome_run)

exome_qc <- exome_results %>% 
  mutate(accession = gsub("_.*$", "", name),
         usable_reads_M = (reads_M * usable_pct)/100) %>% 
  select(accession, tumor_specimen_id, reads_M, usable_reads_M, coverage, usable_pct, 
         dup_pct, tot_var_called, tot_filt_var_called, exome_status)

manifest1 <- readr::read_csv("work/TestIntegrityOfFastqFiles/metadata.csv",
                            col_types = cols("Accession ID" = col_character(),
                                             "Specimen ID" = col_character())) %>% 
  dplyr::filter(!`Accession ID` == "Promega") %>% 
  dplyr::select("project.id", "Accession ID", "Specimen ID", "Panel ID")

manifest <- manifest1

accs <- manifest %>% pull(`Accession ID`) %>% unique()

panels <- pull(manifest, `Panel ID`) %>% unique()

x <- reportGeneration::init_LIS_SF()
x <- reportGeneration::login_LIS_SF(x)

sf_info <- purrr::map_dfr(.x = accs, SF = x, .f = get_sf_info)

summary_files <- list.files(path = "work/RunPrimerDesign/", # can change to "panel_bundle" after panel ordering for updated summary
                            pattern = "selected_summary.csv",
                            full.names = TRUE,
                            recursive=TRUE)

panel_summary <- purrr::map_dfr(.x = summary_files, .f = get_panel_info)

patient_summary <- exome_qc %>% 
  left_join(manifest, by = c("accession" = "Accession ID", "tumor_specimen_id" = "Specimen ID")) %>% 
  left_join(sf_info, by = c("accession" = "Accession")) %>% 
  left_join(panel_summary, by = c("Panel ID" = "Panel")) %>% 
  dplyr::relocate(project.id, Patient, `Panel ID`, .before = accession)  
  
  write_csv(patient_summary, glue("{exome_run}_WES_patient_panel_summary.csv") )
