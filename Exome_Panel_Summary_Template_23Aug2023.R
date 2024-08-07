library(tidyverse)
library(reportGeneration)
library(glue)


# Functions
# Get Salesforce Information
get_orders <- function(accession, SF) {

  res <- reportGeneration::combined_get_order(SF, accession)

  ord <- tibble::tibble(Project = purrr::pluck(res$project, "name", .default = NA),
                        Patient = purrr::pluck(res$patient, "external_pid", .default = NA),
                        Accession = as.character(purrr::pluck(res, "requisition_id", .default = NA)),
#                        Collection_Date = purrr::pluck(res, "collection_date", .default = NA),
#                        Specimen_Type = purrr::pluck(res, "specimen_type", .default = NA),
#                        Order_Date = purrr::pluck(res, "order_date", .default = NA),
#                        Order_Status = purrr::pluck(res, "status", .default = NA),
#                        WES_Complete_Date = purrr::pluck(res$ws, "ws_analysis_date_complete", .default = NA),
#                        Panel = purrr::pluck(res$panel, "panel_id", .default = NA),
#                        Panel_Design_Complete = purrr::pluck(res$panel, "panel_design_complete", .default = NA)
  )
}

# Get panel summary from Panel Design Selected Summary Files
get_panel_info <- function(panel_file) {
  
  ssummary <- readr::read_csv(panel_file) %>% 
    # dplyr::mutate(panel = gsub("_.*", "", amplicon_name)) %>% 
    dplyr::mutate(Panel = gsub('_.*','',basename(panel_file))) %>% 
    dplyr::group_by(Panel) %>% 
    dplyr::summarise(Total = n(),
                     SNVs = sum(grepl("mis", variant_code)),
                     Indels = sum(!grepl("mis", variant_code)),
                     Priority = sum(priority))
                     # Nonpriority_SNVs = sum(!priority & grepl("mis", variant_code)),
                     # Nonpriority_Indels = sum(!priority & !grepl("mis", variant_code)))
}

exome_run <- "AREX12181G"

outdir <- getwd()

# Easiest to simply read in manifest from run directory, no need for s3 download
manifest <- readr::read_csv("work/TestIntegrityOfFastqFiles/metadata.csv",
                            col_types = cols("Accession ID" = col_character(),
                                             "Specimen ID" = col_character())) %>% 
  dplyr::filter(!`Accession ID` == "Promega") %>% 
  dplyr::select("Accession ID", "Specimen ID", "Panel ID")

requests <- unique(dplyr::pull(manifest, 'Accession ID'))

x <- reportGeneration::init_LIS_SF()
x <- reportGeneration::login_LIS_SF(x)

orders <- purrr::map_dfr(.x = requests, SF = x, .f = get_orders)

# Panel Summaries
panels <- manifest %>% pull("Panel ID") %>% unique()

summary_files <- list.files(path = file.path("work/RunPrimerDesign"), 
                            pattern = "selected_summary.csv",
                            full.names = T,
                            recursive = T)

panel_summary <- purrr::map_dfr(.x = summary_files, .f = get_panel_info)

patient_panel_summary <- panel_summary %>% 
  dplyr::inner_join(manifest, by = c("Panel" = "Panel ID")) %>% 
  dplyr::inner_join(orders, by = c("Accession ID" = "Accession")) %>% 
  dplyr::relocate(Project, Patient, `Accession ID`, `Specimen ID`, .before = Panel)

qc <- read_csv("work/DatabaseVarprioSamples/database/Database_loader_1.csv",
               col_types = cols(tumor_specimen_id = col_character())) %>% 
  select(tumor_specimen_id, reads_M, usable_pct, dup_pct, coverage, 
         exome_status) %>% 
  mutate(usable_reads_M = (reads_M * usable_pct)/100) %>% 
  relocate(usable_reads_M, .after = reads_M)

qc_panel <- qc %>% 
  dplyr::inner_join(patient_panel_summary, by = c("tumor_specimen_id" = "Specimen ID")) %>% 
  rename(Specimen = tumor_specimen_id) %>% 
  relocate(Project, Patient, `Accession ID`, Specimen, Panel)

readr::write_csv(qc_panel, file.path(outdir, paste0(exome_run, "_patient_panel_summary.csv")))
