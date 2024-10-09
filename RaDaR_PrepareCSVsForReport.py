# Prepares csvs for calling report only

# Standard library imports
import csv
import os
import argparse
# Related third party imports
import pandas as pd
from luigi import configuration, LocalTarget
# Local application / library specific imports
from inipipeline.external_service_adapters.scripts.luigi_operations import fail_if_outputs_not_created
from inipipeline.external_service_adapters.scripts.global_variables import metadata_header_lines
from inipipeline.external_service_adapters.scripts.utils import (
    get_calling_targets_table_from_metadata,
    get_unique_panels_from_metadata
)
from mrd_pipeline.mrd_analysis.base_tasks.scripts.abstract_inivata_pipeline_task import AbstractInivataPipelineTask
from mrd_pipeline.mrd_analysis.base_tasks.scripts.calling import Calling

parser = argparse.ArgumentParser(prog='PrepareCSVsForReport',
                                 description='prepares csvs for report')

parser.add_argument('-r', '--run')
args = parser.parse_args()


# Update work directory to run script
work_dir = os.path.join('/shared/processed_runs/', args.run)
pos_control_panel = 'SW480'
report_dir = 'reports/calling_report/inputs_from_calling/'


def make_calling_results_paths(working_dir, path_to_metadata):
    # dataframe created explicitly
    df = pd.read_csv(path_to_metadata, keep_default_na=False, skiprows=metadata_header_lines)
    # Keep plasma and cell line specimens, except the positive control (filtered in the 4th line)
    specimen_type_list = df.loc[df['specimen.type'].str.contains('plasma'), 'specimen.type'].unique().tolist()
    specimen_type_list.append('cellLine')
    df = df.loc[df['specimen.type'].isin(specimen_type_list)]
    df = df.loc[df['panel'] != configuration.get_config().get('lims_metadata', 'positive_control_panel')]
    out_path_list = []
    for i in range(len(df.index)):
        out_df = df.iloc[i]
        out_path = os.path.join(working_dir, "calling_output", str(out_df['patient.id']), str(out_df['specimen.id']))
        if not out_path in out_path_list:
            out_path_list.append(out_path)
    return out_path_list


def make_results_csv(calling_result_dirs, outpath):
    if len(calling_result_dirs) == 0:
        results_csv = pd.DataFrame({
            'sample_id': [],
            'input': [],
            'eVAF': [],
            'mutant_molecules': [],
            'mean_VAF': [],
            'LR': [],
            'flag_count': [],
            'all_pass_variants': [],
            'total_variants': [],
            'eVAF_calling': [],
            'n_positive_variants': [],
            'ctDNA_detected': []
        })
        results_csv.to_csv(outpath, index=False)
        return None
    # Read in all of the results.csv files
    results_list = [os.path.join(path, "results.csv") for path in calling_result_dirs]
    results_csv = pd.concat([pd.read_csv(f) for f in results_list])
    # Ensure ctDNA_detected remains a boolean, even if there are NAs
    results_csv['ctDNA_detected'].replace({0: False, 1: True}, inplace=True)
    # Write results to file. Convert NAs to the string 'Fail' only in ctDNA_detected columns
    #LR needs to remain as NA
    results_csv['ctDNA_detected'].fillna(value='Fail', inplace=True)
    results_csv.to_csv(outpath, index=False)


def make_variants_csv(working_dir, target_df, outpath):
    variant_df_list = []
    for i in range(len(target_df.index)):
        temp_df = target_df.iloc[i]
        variant_path = os.path.join(working_dir, "calling_output", str(temp_df['patient.id']),
                                    str(temp_df['specimen.id']), "calling_info.csv")
        variant_df = pd.read_csv(variant_path)
        variant_df['sample_id'] = str(temp_df['specimen.id'])
        variant_df_list.append(variant_df)
    variants_csv = pd.concat(variant_df_list)
    variants_csv.to_csv(outpath, index=False)


path_to_metadata = os.path.join(work_dir, 'userInput', 'metadata.csv')
# Make a list of all the plasma sample directories with results that need to be concatenated
# Note - this list should NOT include positive control
calling_result_dirs = make_calling_results_paths(work_dir, path_to_metadata)
# make the results csv
make_results_csv(calling_result_dirs, os.path.join(work_dir,
                            "reports/calling_report/inputs_from_calling/results.csv"))
# Make variants.csv
target_df = get_calling_targets_table_from_metadata(path_to_metadata)
# Remove positive control
target_df = target_df.loc[target_df['panel'] != pos_control_panel]
make_variants_csv(work_dir, target_df,
                          os.path.join(work_dir, "reports/calling_report/inputs_from_calling/variants.csv"))