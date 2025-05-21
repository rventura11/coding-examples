library(tidyverse)
library(reportGeneration)
library(glue)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("-a", "--arex_code", type = "character", required = TRUE,
                    help = "Directory of panel qc data to generate summary for.")
args <- parser$parse_args()

run_dir <- args$arex_code

setwd(file.path("/shared/processed_runs", run_dir))

## Functions ##

# Get LIS Information
get_orders <- function(accession, SF) {
  
  res <- reportGeneration::combined_get_order(SF, accession)
  # Project <- res$project$name
  Patient <- res$patient$external_pid
  Accession <- accession
  # Specimen <- res$specimens[res$specimens$specimen_type == "Extracted Tumor DNA", 4]
  # Panel <- res$panel$panel_id
  Panel_QC1 <- ifelse(is.null(res$panel$panel_first_qc_status), NA_character_, res$panel$panel_first_qc_status)
  
  ord <- tibble::tibble(Patient, Accession, Panel_QC1)
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

# Get Exome QC
get_exome_qc <- function(accession) {

  ex_q <- paste("SELECT
  *
FROM exome_analysis
WHERE name LIKE '", accession, "_%';", sep = "")

  q_res <- reportGeneration::query.database(ex_q)

}


# Import metadata and files
if(file.exists("userInput/metadata.csv.bak")) {
   m <- read_csv("userInput/metadata.csv.bak", skip = 10) %>% 
    filter(specimen.type == "tumor")
    h <- read_lines('userInput/metadata.csv', n_max = 10)
} else if(!file.exists("userInput/metadata.csv.bak")) {
  m <- read_csv("userInput/metadata.csv", skip = 10) %>% 
    filter(specimen.type == "tumor")
    h <- read_lines('userInput/metadata.csv', n_max = 10)
}

pnl_qc_in <- readr::read_csv("csvs_for_loader/panel_qc.csv")
amp_qc_in <- readr::read_csv("csvs_for_loader/amplicon_qc.csv")
trgt_info_in <- readr::read_csv("csvs_for_loader/target_variant_info.csv")
sample_in <- m %>% select(project.id, accession.id, specimen.id, panel)

run_name <- gsub("basespace.run.id,", "", h[5])

panels <- dplyr::pull(m, panel)
accs <- dplyr::pull(m, accession.id)

x <- reportGeneration::init_LIS_SF()
x <- reportGeneration::login_LIS_SF(x)

orders <- purrr::map_dfr(.x = accs, SF = x, .f = get_orders)

sample_lis <- orders %>% 
  left_join(sample_in, by = c("Accession" = "accession.id")) %>%  
  left_join(pnl_qc_in, c("panel" = "panel_name")) %>%
  select(Project = project.id, Patient, Accession, Sample = specimen.id, Panel = panel, "Tumor Depth QC" = "tumor_depth_qc", "QC Flag" = "tumor_flag")

# WES Query
exome_results <- purrr::map_dfr(.x = accs, .f = get_exome_qc)

exome_qc <- exome_results %>% 
  mutate(accession = gsub("_.*$", "", name),
         usable_reads_M = (reads_M * usable_pct)/100) %>% 
  select(accession, usable_reads_M, coverage, dup_pct, exome_status)

# Priority Query
priority_info <- purrr::map_dfr(.x = panels, .f = get_variant_priority)

variant_summary <- trgt_info_in %>% 
  dplyr::left_join(priority_info, by = c("variant_code", "primer" = "amplicon_name")) %>% 
#  dplyr::mutate(Accession = gsub(".{2}$", "", sampleID)) %>% 
  dplyr::group_by(sampleID, panel_name) %>% 
  dplyr::summarise(Tumor_Pass = sum(is_present),
                   SNV = sum(is_present & grepl("mis", variant_code)),
                   Indel = sum(is_present & !grepl("mis", variant_code)),
                   Priority = sum(is_present & priority),
                   Priority_SNV = sum(is_present & priority & grepl("mis", variant_code)),
                   Priority_Indel = sum(is_present & priority & !grepl("mis", variant_code)))

patient_var_summary <- sample_lis %>% 
  left_join(variant_summary, c("Panel" = "panel_name", "Sample" = "sampleID")) %>% 
  left_join(exome_qc, c("Accession" = "accession")) %>%
  arrange(Project, Accession)

readr::write_csv(patient_var_summary, glue("{run_dir}_Patient_Target_Variant_Summary.csv"))

print(orders, n=40, width=4000)
readr::write_csv(orders, glue("{run_dir}_LIS_Panel_Status.csv"))

