import os
# Related third party imports
import numpy as np
import pandas as pd
import argparse
from luigi import configuration, LocalTarget
# Local application / library specific imports
from mrd_pipeline.mrd_analysis.base_tasks.scripts.abstract_inivata_pipeline_task import AbstractInivataPipelineTask
from mrd_pipeline.mrd_analysis.base_tasks.scripts.make_calling_report import MakeCallingReport
from inipipeline.external_service_adapters.scripts.global_variables import metadata_header_lines
from inipipeline.external_service_adapters.scripts.luigi_operations import fail_if_outputs_not_created
from inipipeline.external_service_adapters.scripts.run_scripts import run_script
from inipipeline.external_service_adapters.scripts.utils import convert_nan_inf_to_na, error_handling_in_csv_file

parser = argparse.ArgumentParser(prog='PrepareCSVsForReport',
                                 description='prepares csvs for report')

parser.add_argument('-r', '--run')
args = parser.parse_args()


# Update work directory to run script
work_dir = os.path.join('/shared/processed_runs/', args.run)

########################################################################################################################
# Convenience functions
########################################################################################################################
def make_sample_call(input_path, output_path):
    # Dataframe created explicitly, not referenced
    sample_call = pd.read_csv(input_path)
    sample_call = sample_call.rename(columns={'sample_id': 'sample_name',
                                              'all_pass_variants': 'supporting_var_n',
                                              'total_variants': 'number_variants_available',
                                              'mean_VAF': 'mean_vaf'})
    # Currently we set the LR threshold as 2
    sample_call['call'] = sample_call.apply(
        lambda x: 1 if x['ctDNA_detected'] is True or x['ctDNA_detected'] == "True" else 0, axis=1
    )
    sample_call['call'] = sample_call['call'].astype(int)
    sample_call = sample_call[["sample_name", "eVAF", "mean_vaf", "LR", "supporting_var_n",
                               "number_variants_available", "call", "n_positive_variants", "mutant_molecules"]]
    convert_nan_inf_to_na(sample_call).to_csv(output_path, index=False)
    error_handling_in_csv_file(output_path,
                               ["sample_name", "eVAF", "mean_vaf", "LR", "supporting_var_n",
                                "number_variants_available", "call", "n_positive_variants",
                                "mutant_molecules"])


#make sample_call.csv
make_sample_call(os.path.join(work_dir,
                              "reports/calling_report/inputs_from_calling/results.csv"),
                 os.path.join(work_dir, "radar_qc_loader_files/sample_call.csv"))