library(tidyverse)
library(reportGeneration)
library(glue)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-a", "--arex_code", type = "character", required = TRUE,
                    help = "Directory of exome data to generate summary for.")
args <- parser$parse_args()

run_dir <- args$arex_code

outdir <- (file.path("review_tables"))
system(paste0('mkdir ', outdir))

## Functions ##

# Get Salesforce Information
get_orders <- function(accession, SF) {

res <- reportGeneration::combined_get_order(SF, accession)
Project <- res$project$name
Patient <- res$patient$external_pid
Accession <- res$requisition_id
#Specimen <- res$specimens[res$specimens$specimen_type == "Extracted Tumor DNA", 2]
Panel <- res$panel$panel_id
Panel_QC1 <- ifelse(is.null(res$panel$panel_first_qc_status), NA_character_, res$panel$panel_first_qc_status)

ord <- tibble::tibble(Project, Patient, Accession, Panel, Panel_QC1)
}

# Get Variant Priority Info
get_variant_priority <- function(panel) {
  query1 <- paste("Select	
    p.name AS panel_name,
		a.name AS amplicon_name,
		tv.priority,
		tv.variant_code,
		tv.AF_T,
		tv.cn,
		tv.effective_AF
From panel AS p
JOIN amplicon AS a ON p.id = a.panel_id
JOIN target_variant AS tv ON a.id = tv.amplicon_id

WHERE p.name = '", panel, "';", sep = "")
  
  res <- reportGeneration::query.database(query1)
  
}

# Get Panel Design Info
get_panel_designs <- function(panel) {
  
  q <- paste("SELECT 
tv.*,
a.id = amplicon_id,
p.id AS panel_id,
p.name AS panel_name
FROM target_variant AS tv
JOIN amplicon AS a ON tv.amplicon_id = a.id
JOIN panel AS p ON a.panel_id = p.id
WHERE p.name = '", panel, "';", sep = "")
  
  q_res <- reportGeneration::query.database(q)
  
}

# Import PanelQC summary files
pnl_qc_in <- readr::read_csv("csvs_for_loader/panel_qc.csv")
amp_qc_in <- readr::read_csv("csvs_for_loader/amplicon_qc.csv")
trgt_info_in <- readr::read_csv("csvs_for_loader/target_variant_info.csv")
sample_in <- readr::read_csv("csvs_for_loader/sample_table.csv")

panels <- dplyr::pull(pnl_qc_in, panel_name)

sample_id <- dplyr::pull(sample_in, sampleID)
accs <- dplyr::pull(sample_in, accessionID)

x <- reportGeneration::init_LIS_SF()
x <- reportGeneration::login_LIS_SF(x)

orders <- purrr::map_dfr(.x = accs, SF = x, .f = get_orders)

project_summary <- orders %>% 
  dplyr::group_by(Project) %>% 
  dplyr::summarise(Panels = n())

print(project_summary)
# readr::write_csv(project_summary, file.path(outdir, "Project_Summary_Table.csv"))

#sample_detail <- orders %>% 
#  dplyr::select(Project, Accession, Panel)

# readr::write_csv(sample_detail, file.path(outdir, "Panel_Sample_Table.csv"))

# QC Detailed Summary Export for JIRA Ticket
pnl_summary <- sample_in %>% 
  dplyr::select(projectID, accessionID, sampleID, panel) %>% 
  dplyr::full_join(select(pnl_qc_in, panel_name, pct_dimers, n_amplicons, n_passed, 
                          amplicon_pass = passed, SNP_passed, NTC_FLAG, 
                          snp_percentage,n_reference_var_pass, 
                          n_tumor_var_pass, tumor_flag), by = c("panel" = "panel_name")) 
  
  
readr::write_csv(pnl_summary, file.path(outdir, glue("{run_dir}_PanelQC_Summary.csv")))

# Priority Query
priority_info <- purrr::map_dfr(.x = panels, .f = get_variant_priority)

variant_summary <- trgt_info_in %>% 
  dplyr::left_join(priority_info, by = c("variant_code", "primer" = "amplicon_name")) %>% 
  dplyr::mutate(Accession = gsub(".{2}$", "", sampleID)) %>% 
  dplyr::group_by(Accession, panel_name) %>% 
  dplyr::summarise(Tumor_Pass = sum(is_present),
                   SNV = sum(is_present & grepl("mis", variant_code)),
                   Indel = sum(is_present & !grepl("mis", variant_code)),
                   Priority = sum(is_present & priority),
                   Priority_SNV = sum(is_present & priority & grepl("mis", variant_code)),
                   Priority_Indel = sum(is_present & priority & !grepl("mis", variant_code)))

patient_var_summary <- orders %>% 
  dplyr::select(Project, Patient, Accession, Panel) %>% 
  dplyr::left_join(variant_summary, by = "Accession") %>% 
  dplyr::arrange(Project, Accession) %>% 
  distinct(Accession, .keep_all = TRUE) %>%
  dplyr::select(-panel_name)

readr::write_csv(patient_var_summary, file.path(outdir, glue("{run_dir}_Patient_Target_Variant_Summary.csv")))

# Get Panel Design Details If Needed
design_res <- purrr::map_dfr(.x = panels, .f = get_panel_designs)

design_summary <- design_res %>%
  dplyr::mutate(Accession = gsub("P0[1-9]$", "", panel_name)) %>%
  dplyr::group_by(panel_name) %>%
  dplyr::summarise(tot = n(),
                   SNVs = sum(grepl("mis", variant_code)),
                   Indels = sum(!grepl("mis", variant_code)),
                   Priority = sum(priority),
                   Priority_SNVs = sum(priority & grepl("mis", variant_code)),
                   Priority_Indels = sum(priority & !grepl("mis", variant_code)),
                   Nonpriority_SNVs = sum(!priority & grepl("mis", variant_code)),
                   Nonpriority_Indels = sum(!priority & !grepl("mis", variant_code)))

readr::write_csv(design_summary, file.path(outdir, glue("{run_dir}_Panel_Design_Summary.csv")))