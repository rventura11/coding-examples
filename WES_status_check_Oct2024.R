library(reportGeneration)
library(dplyr)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("-a", "--arex", type = "character", required = TRUE,
                    help = "Analysis to generate project reports for")

args <- parser$parse_args()

# Set run directory
rundir <- args$arex 

setwd(file.path("/shared/pipeline-user/run_data", rundir))

get_exome_status <- function(accession, SF) {

res <- reportGeneration::combined_get_order(SF, accession)
Project <- res$project$name
Patient <- res$patient$external_pid
Accession <- res$requisition_id
#Specimen <- res$specimens[res$specimens$specimen_type == "Extracted Tumor DNA", 2]
Panel <- res$panel$panel_id
Status <- res$status

status <- tibble::tibble(Project, Patient, Accession, Panel, Status)
}



meta <- readr::read_csv("metadata.csv", skip=10)
accs <- dplyr::pull(meta, accession.id)



status <- purrr::map_dfr(.x = accs, SF = x, .f = get_exome_status)

print(status, n=50, width=500)