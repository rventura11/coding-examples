#!/usr/bin/env Rscript --vanilla
rm(list=ls())

library(reportGeneration)
library(tidyverse)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("-p", "--project", type = "character", required = TRUE,
                    help = "Project to generate summary report for")

parser$add_argument("-pl", "--include_plasma", action = "store_true", required = FALSE,
                    help = "Include plasma pipeline results")

args <- parser$parse_args()

# Functions
get_neo_qc <- function(project) {
  
  s3_cmd <- paste0("aws s3 cp s3://ini-analysis-prod/radar_clinical/neo_qc_sheets/ ./ --recursive  --exclude '*' --include '*.xlsx' --exclude 'archived/*'")
  system(s3_cmd)
  
  neo_file <- list.files(pattern = "Cumulative.xlsx")
  
  neoqc <- readxl::read_excel(neo_file, 
                              sheet = 1, na = "N/A") %>% 
    dplyr::filter(`Sponsor Name` == "Inivata Inc" & Project == project) %>% 
    dplyr::select(8:37)
  
  rm_cmd <- paste0("rm ", neo_file)
  system(rm_cmd)
  
  return(neoqc)
  
}

get_panel_design <- function(project) {
  
  design_q <- paste("SELECT 
pj.name AS project_name,
r.Accession,
p.id AS panel_id,
p.name AS Panel,
a.id AS amplicon_id,
a.name AS amplicon_name,
tv.*
FROM project AS pj
JOIN request AS r on pj.id = r.project_id
JOIN panel AS p on r.id = p.request_id
JOIN amplicon AS a ON p.id = a.panel_id
JOIN target_variant AS tv ON a.id = tv.amplicon_id
WHERE pj.name = '", project, "';", sep = "")
  
  q_res <- reportGeneration::query.database(design_q)
  
}

get_exome <- function(accession) {
  
  ex_q <- paste("SELECT
	*
FROM exome_analysis
WHERE name LIKE '", accession, "_%';", sep = "")
  
  q_res <- reportGeneration::query.database(ex_q)
  
}

get_tumor_order_info <- function(accession, SF) {
  
  res <- reportGeneration::combined_get_order(SF, accession)
  Project <- res$project$name
  Accession <- res$requisition_id
  Patient <- res$patient$external_pid
#  Panel <- res$panel$panel_id > Removing panel due to changes in Salesforce
  
  o.df <- tibble::tibble(Project, Accession, Patient)
}

get_pqc_amps <- function(panel) {
  
  pqc_q <- paste("SELECT
r.accession AS Accession,
p.name AS Panel,
pq.Tumor_Flag,
ar.arex_code,
pr.name AS run_name,
a.name AS amplicon_name,
aq.reads_all_k,
tv.id AS target_id,
tv.priority,
tv.variant_code,
tv.cn,
tv.effective_AF
From panel AS p
JOIN request AS r ON p.request_id = r.id
JOIN panel_qc AS pq ON p.id = pq.panel_id
JOIN analytical_run_experiment AS ar ON pq.analytical_run_experiment_id = ar.id
JOIN physical_run_experiment AS pr ON ar.physical_run_experiment_id = pr.id
JOIN workflow_version AS wv ON pr.workflow_version_id = wv.id
JOIN amplicon AS a ON p.id = a.panel_id
JOIN amplicon_qc AS aq ON a.id = aq.amplicon_id and pq.id = aq.panel_qc_id
JOIN target_variant AS tv ON a.id = tv.amplicon_id
WHERE wv.workflow_id = 4
AND p.name = '", panel, "'
AND pq.updated_at = (SELECT max(pq.updated_at) FROM panel p, panel_qc pq, 
                 analytical_run_experiment ar, physical_run_experiment pr, 
                 workflow_version wv 
                 WHERE p.id = pq.panel_id 
                 AND pq.analytical_run_experiment_id = ar.id 
                 AND ar.physical_run_experiment_id = pr.id 
                 AND pr.workflow_version_id = wv.id 
                 AND wv.workflow_id = 4 
                 AND p.name = '", panel, "');", sep = "")
  
  q_res <- reportGeneration::query.database(pqc_q)
}

get_pqc_targets <- function(panel) {
  
  pqc_q <- paste("SELECT
r.accession AS Accession,
p.name AS Panel,
ex.name AS Tumor,
ex.type,
ar.arex_code,
tvi.target_variant_id,
tvi.is_present,
tvi.AF
From panel AS p
JOIN request AS r ON p.request_id = r.id
JOIN panel_qc AS pq ON p.id = pq.panel_id
JOIN analytical_run_experiment AS ar ON pq.analytical_run_experiment_id = ar.id
JOIN physical_run_experiment AS pr ON ar.physical_run_experiment_id = pr.id
JOIN workflow_version AS wv ON pr.workflow_version_id = wv.id
LEFT JOIN extracted_dna AS ex ON r.id = ex.request_id
LEFT JOIN sample AS s ON ex.id = s.extracted_dna_id AND s.analytical_run_experiment_id = ar.id
LEFT JOIN target_variant_info AS tvi ON s.id = tvi.sample_id
WHERE wv.workflow_id = 4
AND ex.type = 'tumor'
AND p.name = '", panel, "'
AND pq.updated_at = (SELECT max(pq.updated_at) FROM panel p, panel_qc pq, 
                 analytical_run_experiment ar, physical_run_experiment pr, 
                 workflow_version wv 
                 WHERE p.id = pq.panel_id 
                 AND pq.analytical_run_experiment_id = ar.id 
                 AND ar.physical_run_experiment_id = pr.id 
                 AND pr.workflow_version_id = wv.id 
                 AND wv.workflow_id = 4 
                 AND p.name = '", panel, "');", sep = "")
  
  q_res <- reportGeneration::query.database(pqc_q)
}

get_rdr_call_results <- function(project) {
  
  q <- paste("SELECT 
pj.name AS Project,
r.Accession,
s.name AS sample_name,
ex.name AS dna_name,
ex.type AS dna_type,
p.name AS Panel,
ar.arex_code,
pr.name AS run_name,
sqc.extracted_plasma_volume_ml,
sqc.input_amount_copies_short,
sqc.recovery_percent,
sqc.inferred_gender,
sr.qc_flag,
sr.errorrate_flag,
sr.merging_flag,
sr.interbarcodesnpconsistency_flag,
sr.intersamplesnpconsistency_flag,
sr.contamination_flag,
sr.gender_flag,
sr.technical_flag,
sc.mean_VAF,
sc.eVAF,
sc.LR,
sc.supporting_var_n,
sc.number_variants_available,
sc.n_positive_variants,
sc.mutant_molecules,
sc.call
FROM project AS pj
JOIN request AS r ON pj.id = r.project_id
JOIN extracted_dna AS ex ON r.id = ex.request_id
JOIN sample AS s ON ex.id = s.extracted_dna_id
JOIN sample_run_qc AS sqc ON s.id = sqc.sample_id
JOIN panel AS p ON sqc.panel_id = p.id
JOIN sample_reports AS sr ON sqc.id = sr.sample_run_qc_id
JOIN sample_call AS sc ON s.id = sc.sample_id
JOIN analytical_run_experiment AS ar ON s.analytical_run_experiment_id = ar.id
JOIN physical_run_experiment AS pr ON ar.physical_run_experiment_id = pr.id
JOIN workflow_version AS wv ON pr.workflow_version_id = wv.id
WHERE wv.workflow_id = 3
AND ex.type LIKE '%plasma%'
AND sr.is_current_version = 1
AND pj.name = '", project, "';", sep = "")
  
  rdr_res <- reportGeneration::query.database(q)
  
}

get_rdr_orders <- function(accession, SF) {
  
  res <- reportGeneration::combined_get_order(SF, accession)
  Project <- res$project$name
  Accession <- res$requisition_id
  Patient <- res$patient$external_pid
#  Panel <- res$panel$panel_id > Removing panel due to changes in Salesforce
  Specimen <- res$specimens[grepl("Whole Blood|Plasma", res$specimens$specimen_type), 3]
  Date <- as.character(ifelse(is.na(res$specimens[grepl("Whole Blood|Plasma", res$specimens$specimen_type), 2]),
                              res$collection_date, 
                              res$specimens[grepl("Whole Blood|Plasma", res$specimens$specimen_type), 2]))
  Timepoint <- as.character(ifelse(is.na(res$specimens[grepl("Whole Blood|Plasma", res$specimens$specimen_type), 1]),
                                   ifelse(is.null(res$timepoint), NA_character_, 
                                          res$timepoint),
                                   res$specimens[grepl("Whole Blood|Plasma", res$specimens$specimen_type), 1]))
  
  o.df <- tibble::tibble(Project, Accession, Patient, Specimen, Date, Timepoint)
}

project_name <- args$project

report_date <- format(Sys.Date(), format = "%d-%b-%Y")

neoqc <- get_neo_qc(project_name)

panel_design <- purrr::map_dfr(.x = project_name, .f = get_panel_design)

panel_accs <- unique(dplyr::pull(panel_design, Accession))

panels <- unique(dplyr::pull(panel_design, Panel))

exome_results <- purrr::map_dfr(.x = panel_accs, .f = get_exome)

x <- reportGeneration::init_LIS_SF()
x <- reportGeneration::login_LIS_SF(x)

tumor_orders <- purrr::map_dfr(.x = panel_accs, SF = x, .f = get_tumor_order_info)

panel_qc_amps <- purrr::map_dfr(.x = panels, .f = get_pqc_amps)

panel_qc_targets <- purrr::map_dfr(.x = panels, .f = get_pqc_targets)

panel_qc <- panel_qc_amps %>% 
  dplyr::left_join(panel_qc_targets, by = c("Accession", "Panel", "arex_code", "target_id" = "target_variant_id"))

# Exome Summary
exome_summary <- exome_results %>% 
  dplyr::select(Tumor = tumor_specimen_id, `Reads M` = reads_M, `Usable PCT` = usable_pct,
                `Dup PCT` = dup_pct, Coverage = coverage, `Exome QC` = exome_status) %>% 
  dplyr::mutate(Accession = gsub("0[1-9]$", "", Tumor),
                `Usable Reads M` = `Reads M` * (`Usable PCT`/100)) %>%
  dplyr::inner_join(tumor_orders, by = "Accession") %>% 
  dplyr::arrange(Project, Patient, Accession) %>% 
  dplyr::relocate(Project, Patient, Accession, Tumor) %>% 
  dplyr::relocate(`Usable Reads M`, .after = `Usable PCT`)

# Panel Design Summary
design_summary <- panel_design %>% 
  dplyr::group_by(Accession, Panel) %>% 
  dplyr::summarise(`Total Variants` = n(),
                   SNVs = length(grep("mis", variant_code)),
                   Indels = length(grep("ins|del", variant_code)),
                   Priority = length(variant_code[priority == 1]),
                   `Priority SNVs` = length(grep("mis", variant_code[priority == 1])),
                   `Priority Indels` = length(grep("ins|del", variant_code[priority == 1]))) %>% 
  dplyr::mutate(`Panel Design Success` = ifelse(`Total Variants` >= 8, TRUE, FALSE)) %>% 
  dplyr::inner_join(tumor_orders, by = c("Accession")) %>% 
  dplyr::arrange(Project, Patient, Accession) %>%
  dplyr::relocate(Project, Patient, Accession, Panel)

patient_panels <- design_summary %>% 
  dplyr::select(Project, Patient, Accession, Panel, `Total Variants`)

# Panel QC Summary
panel_qc_summary <- panel_qc %>% 
  dplyr::group_by(Panel, Tumor_Flag) %>% 
  dplyr::summarise(# `Total Designed Variants` = n(),
    `Total Passing Variants` = (length(variant_code[is_present ==1 & !is.na(is_present)])),
    `Passing SNVs` = length(grep("mis", variant_code[is_present == 1 & !is.na(is_present)])),
    `Passing Indels` = length(grep("ins|del", variant_code[is_present == 1 & !is.na(is_present)])),
    `Passing Priority` = length(variant_code[priority == 1 & is_present == 1 & !is.na(is_present)]),
    `Passing Priority SNVs` = length(grep("mis", variant_code[priority == 1 & is_present == 1 & !is.na(is_present)])),
    `Passing Priority Indels` = length(grep("ins|del", variant_code[priority == 1 & is_present == 1& !is.na(is_present)]))) %>% 
  dplyr::right_join(patient_panels, by = c("Panel")) %>% 
  dplyr::rename(`Tumor Flag` = Tumor_Flag) %>% 
  dplyr::arrange(Project, Patient, Accession) %>%
  dplyr::relocate(Project, Patient, Accession, Panel, `Total Variants`, `Tumor Flag`)

# RaDaR Summary
if(args$include_plasma) {
  
  radar_results <- purrr::map_dfr(.x = project_name, .f = get_rdr_call_results)
  
  radar_accs <- unique(dplyr::pull(radar_results, Accession))
  
  radar_orders <- purrr::map_dfr(.x = radar_accs, SF = x, .f = get_rdr_orders)
  
    radar_summary <- radar_results %>% 
      dplyr::select(Accession, Panel, `Plasma Sample` = sample_name, 
                    Input = input_amount_copies_short, `QC Flag` = qc_flag, eVAF,
                    Score = LR, `All Pass Variants` = supporting_var_n,
                    `Variants Available` = number_variants_available,
                    `Positive Variants` = n_positive_variants, 
                    `ctDNA Detected` = call) %>% 
      dplyr::inner_join(radar_orders, by = c('Accession', 'Plasma Sample' = 'Specimen')) %>% 
      dplyr::rename(`Collection Date` = Date) %>% 
      dplyr::mutate(`ctDNA Detected` = ifelse(`ctDNA Detected`, TRUE, FALSE),
                    `Sample Input` = as.numeric(Input) * 4) %>% 
      dplyr::relocate(Project, Patient, `Collection Date`, Timepoint, Panel) %>% 
      dplyr::relocate(`Sample Input`, .after = Input) %>% 
      dplyr::arrange(Project, Patient, `Collection Date`)
}

if(args$include_plasma) {
  save(project_name, neoqc, exome_summary, design_summary, panel_qc_summary, radar_summary,
       file = paste0(project_name, "_Summary_Report_", report_date, ".RData"))
} else {
  save(project_name, neoqc, exome_summary, design_summary, panel_qc_summary, 
       file = paste0(project_name, "_Summary_Report_", report_date, ".RData"))
}