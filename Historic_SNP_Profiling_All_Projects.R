library(tidyverse)
library(reportGeneration)
library(glue)

project <- 'project'

# Functions
get_rdr_orders <- function(accession, SF) {
  
  res <- reportGeneration::combined_get_order(SF, accession)
  Project <- res$project$name
  Accession <- res$requisition_id
  Patient <- res$patient$external_pid
  #  Panel <- res$panel$panel_id > Removing panel due to changes in Salesforce
  Specimen <- res$specimens[grepl("Whole Blood|Plasma|cfDNA", res$specimens$specimen_type), 2]
  Date <- as.character(ifelse(is.na(res$specimens[grepl("Whole Blood|Plasma|cfDNA", res$specimens$specimen_type), 3]),
                              res$collection_date, 
                              res$specimens[grepl("Whole Blood|Plasma|cfDNA", res$specimens$specimen_type), 3]))
  Timepoint <- as.character(ifelse(is.na(res$specimens[grepl("Whole Blood|Plasma|cfDNA", res$specimens$specimen_type), 1]),
                                   ifelse(is.null(res$timepoint), NA_character_, 
                                          res$timepoint),
                                   res$specimens[grepl("Whole Blood|Plasma|cfDNA", res$specimens$specimen_type), 1]))
  
  o.df <- tibble::tibble(Project, Accession, Patient, Specimen, Date, Timepoint)
}

dup_pairs = function(x, y) duplicated(apply(cbind(x, y), 1, function(x) paste(sort(x), collapse=":")))

haplotype_matrix = function(df, snps=NULL) {
  if (is.null(snps)) snps = unique(df$rsid)
  m = matrix(NA, nrow = length(unique(df$sample_id)), ncol=length(snps))
  colnames(m) = snps
  rownames(m) = unique(df$sample_id)
  m[cbind(df$sample_id, df$rsid)] = df$haplotype
  m
}

matching_snps = function(x, y) {
  if (nrow(x) != nrow(y)) stop("must compare equall numbers of SNPs")
  map_int(1:nrow(x), function(i) sum(x[i,] == y[i,], na.rm=TRUE))
}

common_snps = function(x, y) {
  if (nrow(x) != nrow(y)) stop("must compare equall numbers of SNPs")
  map_int(1:nrow(x), function(i) sum(!is.na(x[i,]) & !is.na(y[i,]), na.rm=TRUE))
}

compare_historic_snps <- function(df1, db_snps) {
  
  all_snps = unique(db_snps$rsid)
  
  mat1 = haplotype_matrix(df1, all_snps)
  mat2 = haplotype_matrix(db_snps, all_snps)
  
  # build a data frame of comparisons
  pairs = expand_grid(
    id1 = unique(df1$sample_id),
    id2 = unique(db_snps$sample_id)
  ) %>%
    filter(
      id1 != id2,
      !dup_pairs(id1, id2) # remove reverse comparisons
    )
  
  # calculate common and matching snps
  snp_db_comparisons <- pairs %>%
    mutate(
      common = common_snps(mat1[id1,, drop=F], mat2[id2,, drop=F]),
      matching = matching_snps(mat1[id1,, drop=F], mat2[id2,, drop=F]),
      differences = common - matching
    )
  
  return(snp_db_comparisons)
  
}

# RaDaR SNP Query > Could add panel qc workflow id for tumors too
query <- paste("SELECT 
  pj.name AS project_name,
	r.accession,
	s.name AS sample_id,
	ex.name AS dna_name,
	ex.type AS dna_type,
	ar.arex_code,
	pr.name AS run_name,
	sr.qc_flag,
  sr.errorrate_flag,
  sr.merging_flag,
  sr.interbarcodesnpconsistency_flag,
  sr.intersamplesnpconsistency_flag,
	ssv.value AS haplotype,
	ssv.rsid
FROM project AS pj
JOIN request AS r ON pj.id = r.project_id
JOIN extracted_dna AS ex ON r.id = ex.request_id
JOIN sample AS s ON ex.id = s.extracted_dna_id
JOIN sample_run_qc AS srq ON s.id = srq.sample_id
JOIN sample_reports AS sr ON srq.id = sr.sample_run_qc_id
JOIN sample_snp_variant AS ssv ON srq.id = ssv.sample_run_qc_id 
JOIN analytical_run_experiment AS ar ON s.analytical_run_experiment_id = ar.id
JOIN physical_run_experiment AS pr ON ar.physical_run_experiment_id = pr.id
JOIN workflow_version AS wv ON pr.workflow_version_id = wv.id
WHERE wv.workflow_id = 3;", sep = "")  

db_snps <- reportGeneration::query.database(query)

# Failing samples
fail_snps <- db_snps %>% 
  dplyr::mutate(haplotype = ifelse(haplotype == "-", NA, haplotype)) %>% 
  dplyr::filter(interbarcodesnpconsistency_flag == "FAIL",
                project_name == project)

fail_samps <- fail_snps %>% 
  dplyr::select(project_name, accession, sample_id, dna_name, dna_type, arex_code,
                run_name, qc_flag, interbarcodesnpconsistency_flag, intersamplesnpconsistency_flag) %>% 
  unique()

write_csv(fail_samps, glue("{project}_SNP_Fails.csv"))

snp_compare <- compare_historic_snps(df1 = fail_snps, db_snps = db_snps)

matching_snps <- snp_compare %>% 
  dplyr::filter(common > 16,
                differences < 3)

samps <- c(matching_snps$id1, matching_snps$id2)
accs <- gsub("[0-9][0-9]$", "", gsub("_bc1|_bc|_pl", "", samps)) %>% unique()

# Sf info
x <- reportGeneration::init_LIS_SF()
x <- reportGeneration::login_LIS_SF(x)

radar_orders <- purrr::map_dfr(.x = accs, SF = x, .f = get_rdr_orders)

radar_orders <- radar_orders %>% 
  filter(!is.na(Timepoint))
  
matching_snps_sf <- matching_snps %>% 
  dplyr::left_join(select(radar_orders, Specimen, pat1 = Patient, Timepoint), by = c("id1" = "Specimen")) %>% 
  dplyr::left_join(select(radar_orders, Project, Specimen, pat2 = Patient, Timepoint), by = c("id2" = "Specimen")) %>% 
  dplyr::filter(pat1 != pat2) %>% 
  dplyr::select(failing_sample = id1,failing_patient = pat1, failing_timepoint = Timepoint.x,
                matching_project = Project, matching_sample = id2, matching_patient = pat2, matching_timepoint = Timepoint.y,
                common, matching, differences) %>% 
  dplyr::arrange(failing_patient, failing_sample)

write_csv(matching_snps_sf, glue("{project}_SNP_Alignments_All_Projects.csv"))


