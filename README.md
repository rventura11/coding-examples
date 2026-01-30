# Clinical Bioinformatics Code Examples

R and Python scripts developed at NeoGenomics for routine bioinformatics analysis tasks. Demonstrates skills in algorithmic problem-solving, production pipeline integration and automated data processing.

---

## Featured Scripts

### Historic_SNP_Profiling_All_Projects.R

**Purpose:** Quality control tool for comparing historical SNP genotype profiles to investigate possible sample swaps and mix-ups across clinical oncology projects.

**Technical Highlights:**
- Custom algorithmic functions for pairwise haplotype comparison across hundreds of samples
- Matrix-based genotype matching algorithm to identify samples with >95% SNP concordance
- Combinatorial analysis preventing duplicate comparisons in large datasets
- Integration with clinical database and Salesforce API for patient metadata

**Key Techniques:**
```r
# Custom haplotype comparison algorithm
haplotype_matrix()    # Converts SNP calls to matrix format
matching_snps()       # Counts concordant genotypes between sample pairs
compare_historic_snps()  # Orchestrates pairwise comparison workflow
```

**Use Case:** Identifying cross-contamination events where plasma samples from different patients show >16 matching SNPs, flagging potential sample swaps or labeling issues requiring further investigation.

**Technologies:** R, tidyverse, SQL, Salesforce API, custom algorithms

---

### RaDaR_Project_Summary_Report.R

**Purpose:** Automated end-to-end reporting pipeline integrating whole exome sequencing (WES) lab QC/pipeline analysis data, panel design information, panel QC/plasma ctDNA pipeline analysis results, and patient metadata/timepoints for each patient associated with a specified clinical oncology project.

**Technical Highlights:**
- Multi-source data integration: AWS S3 backup data buckets, SQL databases and Salesforce API
- Conditional workflow execution based on argparse flags (tumor-only vs. tumor+plasma analysis)
- Complex SQL queries with 6-8 table joins for extracting nested relational data
- Production-ready error handling and data validation

**Data Sources Integrated:**
```r
get_neo_qc()          # AWS S3: WES QC spreadsheets - from Neo lab backup data bucket
get_panel_design()    # SQL: Panel design information - from analysis results database (Inidata)
get_exome()           # SQL: WES pipeline analysis data - from analysis results database (Inidata)
get_pqc_amps()        # SQL: Panel QC pipeline analysis data - from analysis results database (Inidata)
get_tumor_order_info() # Salesforce: Patient order information
get_rdr_call_results() # SQL: Plasma analysis ctDNA detection results - from analysis results database (Inidata)
get_rdr_orders()      # Salesforce: Patient metadata and timepoint data
```

**Output:** Comprehensive `.RData` object containing whole exome sequencing QC data, panel design statistics, variant detection results, and plasma ctDNA status for entire clinical cohorts that can then be exported to an Excel spreadsheet for stakeholders.

**Technologies:** R, tidyverse, argparse, AWS CLI, SQL, Salesforce API, dplyr

---

### build_script.py

**Purpose:** Automated script generator for re-running variant calling pipelines after identifying and correcting problematic amplicons in panel designs.

**Technical Highlights:**
- Logic to identify failed amplicons by parsing primer design files for specific primer IDs
- Nested data structure (dictionary of lists) mapping panels to problematic variants
- Pandas DataFrame manipulation with lambda functions to update variant presence flags
- Bash command generation for re-running downstream calling and reporting workflows

**Key Techniques:**
```python
#Line 47
def identify_bad_variants(panel_design_folder_path): 
    # Scans design files for primers matching specifed primer IDs
    # Builds mapped panel ID → variant dict for exclusion of affected variants

#Line 74
panel_df['tumour_present'] = panel_df.apply(
    lambda x: False if x['code'] in bad_variants else x['tumour_present'], 
    axis=1
)
    # Updates calling panel CSVs to denote affected variants as 'absent'

#Line 89
def build_calling_command(run_name, sample_name, patient_name, panel_name)
    # Generates calling commands for affected samples using sample data extracted from pipeline analysis metadata
```

**Use Case:** Automated removal of fixed amplicons that are included in sequencing workflows for some clinical oncology studies but not desired for ctDNA detection analysis. Identifies specified primers, updates variant flags, and generates commands for re-analysis.

**Technologies:** Python, pandas, argparse, file I/O, os operations

---

## About

These scripts were developed during clinical bioinformatics work at NeoGenomics (2021-2025), supporting:
- WES panel design, panel QC and ctDNA plasma analysis (RaDaR liquid biopsy/MRD testing suite)
- Clinical quality control and regulatory compliance (CLIA/CAP)
- Multi-center clinical trial data analysis

For production-ready package development demonstrating modular design and API integration, see [vcf-annotator](https://github.com/rventura11/vcf-annotator).

---

## Technical Stack

**Languages:** R, Python  
**Data Processing:** tidyverse (dplyr, purrr, tidyr), pandas  
**APIs & Infrastructure:** Salesforce, AWS S3, SQL databases  
**Visualization:** ggplot2  
**Workflows:** Argparse, functional programming, pipeline automation
