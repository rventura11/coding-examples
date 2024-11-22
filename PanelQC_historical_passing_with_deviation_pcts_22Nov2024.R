rm(list=ls())

library(reportGeneration)
library(tidyverse)


# Starting date from which you would like to gather historical data
start_date <- "2023-01-01"

# Functions
get_pqc_flags_vars <- function(date) {
  
  pqc_info <- paste("SELECT
pj.name as project,
r.accession,
p.name as panel,
pq.passed as amplicon_pass,
pq.snp_passed,pq.ntc_flag,
pq.tumor_flag,pq.n_tumor_var_pass,
ar.arex_code,pr.name as pq_run
FROM project pj,
request r,
panel p,
panel_qc pq,
analytical_run_experiment ar,
physical_run_experiment pr,
workflow_version wv
WHERE pj.id = r.project_id and
r.id = p.request_id and
p.id = pq.panel_id and
pq.analytical_run_experiment_id = ar.id and
ar.physical_run_experiment_id = pr.id and
pr.workflow_version_id = wv.id and
wv.workflow_id = 4 and
ar.report_generated_timestamp > ('", date, "');", sep = "")
  
  q_res <- reportGeneration::query.database(pqc_info)
}

panel_qc_flags_vars <- purrr::map_dfr(.x = start_date, .f = get_pqc_flags_vars)

#removing projects where we didn't receive tumor (causing all panels to pass with deviation)
panel_qc_filtered <- panel_qc_flags_vars %>% filter(!project %in% c("BIO2", "Moderna-mRNA", "NT-125-101"))

#dealing with duplicates...start by splitting data into two dfs - one with all duplicates and one with none
duplicates <- panel_qc_filtered[duplicated(panel_qc_filtered$panel)|duplicated(panel_qc_filtered$panel, fromLast=TRUE),]
nodups <- panel_qc_filtered %>% 
	group_by(panel) %>% 
	filter(n() == 1)

#filtering out only 1 row if n_tumor_var_pass are the same (we want to include these data where duplicate panels passed with deviation (or just passed) in both runs) - first line of function
#also filtering out all duplicate panels where they passed with deviation only once (inconsistent between duplicate runs so we'll toss all of them out) - middle 2 lines of function
#lastly keeping distinct panels where n_tumor_var_pass ≠ 999 (if they passed normally in all duplicates, can include a representative panel in final data) - last line of function
duplicates_filtered <- duplicates %>%
  distinct(panel, n_tumor_var_pass, .keep_all = TRUE) %>%
  group_by(panel) %>%
  filter(!(n() > 1 & sum(n_tumor_var_pass == 999) == 1)) %>%
  distinct(panel, .keep_all = TRUE)

#now bring duplicates back in
panel_qc_final <- rbind(nodups, duplicates_filtered)

write_csv(panel_qc_final, "PanelQC_historical_data_raw.csv")

#creating stats for bar graph
PWD_percent_by_project <- panel_qc_filtered %>% group_by(project) %>% summarize(percentage_passing_with_deviation = mean(n_tumor_var_pass == 999))
add_count1 <- panel_qc_filtered %>% group_by(project) %>% count(project, name = "total_panels")
add_count2 <- panel_qc_filtered %>% group_by(project) %>% summarize(total_PWD = sum(n_tumor_var_pass == 999))

graph_data <- add_count2 %>%left_join(add_count1) %>% left_join(PWD_percent_by_project)

write_csv(graph_data, "PanelQC_historical_PWD_data.csv")

#bar graph creation - done locally

library(readr)
library(ggplot2)
library(scales)

graph_data <- read_csv("PanelQC_historical_PWD_data.csv")

mycolors <- c("darkgreen","yellowgreen","lightyellow","orange","salmon","red")

#Adjust title below as needed (5th line in function, labs(title = "YOUR TITLE"))
bar_graph <- ggplot(graph_data, aes(x=project, y=percentage_passing_with_deviation, fill=percentage_passing_with_deviation)) + 
	geom_bar(stat = "identity") + 
	theme(axis.text.x = element_text(angle = -55, vjust = -0.1, hjust=0.03), plot.title = element_text(hjust = 0.5)) + 
	scale_y_continuous(breaks=seq(0,1, by = .10), labels=scales::percent) + 
	labs(title = "2023-2024 Panel QC - Percent Passing w/Deviation by Project", y = "% Pass w/Deviation", x = "Project") + 
	scale_fill_gradientn(colors=mycolors, space='Lab', labels = label_percent()) + 
	guides(fill=guide_legend(title=""))

ggsave("PanelQC_historical_PWD_graph.png", plot = bar_graph)