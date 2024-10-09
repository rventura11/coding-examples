import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(prog='CallScriptBuilder',
                                 description='What the program does',
                                 epilog='Text at the bottom of help')

parser.add_argument('-r', '--run')
args = parser.parse_args()

run_folder = os.path.join('/shared/processed_runs', args.run)

# Backup Portion
def create_backup_commands(run_folder):
    commands = []
    # calling
    out = os.path.join(run_folder, 'calling_output')
    bak = os.path.join(run_folder, 'calling_output_bak')
    commands.append(f'cp -r {out} {bak}')

    # reports
    out = os.path.join(run_folder, 'reports')
    bak = os.path.join(run_folder, 'reports_bak')
    commands.append(f'cp -r {out} {bak}')

    # qc loader
    out = os.path.join(run_folder, 'radar_qc_loader_files')
    bak = os.path.join(run_folder, 'radar_qc_loader_files_bak')
    commands.append(f'cp -r {out} {bak}')

    # panel
    out = os.path.join(run_folder, 'radar_qc_loader_files')
    bak = os.path.join(run_folder, 'radar_qc_loader_files_bak')
    commands.append(f'cp -r {out} {bak}')

    return commands

commands = []
commands.append('# Backup Files\n')
commands += create_backup_commands()

# Variants Fix
# /shared/processed_runs/AREX12165A/summary_design_files
design_folder_path = os.path.join(run_folder,'summary_design_files')

def identify_bad_variants(panel_design_folder_path):
    files = os.listdir(panel_design_folder_path)
    bad_var_dict = dict()
    for file in files:
        panel = file.split("_")[0]
        # primer,amp_len,orig_code,variant_id
        panel_df = pd.read_csv(file)
        for index, row in panel_df.iterrows():
            primer = row['primer']
            variant = row['orig_code']
            if '_999_C002_003_fg_001_pp' in primer or '_999_ex21_002_se_001_pp' in primer:
                bad_var_dict[panel] = bad_var_dict.get(panel, [])
                bad_var_dict[panel].append(variant)

    return bad_var_dict

bad_variants_dict = identify_bad_variants(design_folder_path)

panels = set(panel in bad_variants_dict.keys())
for panel in bad_variants_dict.keys():
    bad_variants = bad_variants_dict[panel]
    panel_file = os.path.join(run_folder, f'panel/sample/panel/{panel}_calling_panel.csv')
    panel_file_bak = os.path.join(run_folder, f'panel/sample/panel/{panel}_calling_panel_bak.csv')
    # code,transition,vaf_buffy,buffy_present,vaf_tumour,tumour_present
    panel_df = pd.read_csv(panel_file)
    panel_df.to_csv(panel_file_bak, index=False) # create a backup copy

    panel_df['tumour_present'] = panel_df.apply(lambda x: False if x['code'] in bad_variants else  x['tumour_present'], axis=1)
    panel_df.to_csv(panel_file, index=False) # overwrite the old file...


global_panel_path = os.path.join(run_folder, f'panel/global/reports_panel.csv')
global_panel_path_bak = os.path.join(run_folder, f'panel/global/reports_panel_bak.csv')

# code,transition,vaf_buffy,buffy_present,vaf_tumour,tumour_present,panel_id
global_panel_df = pd.read_csv(global_panel_path)
global_panel_df.to_csv(global_panel_path_bak, index=False)
global_panel_df['tumour_present'] = global_panel_df.apply(lambda x: False if x['code'] in bad_variants_dict.get(x['panel_id'],[]) else x['tumour_present'], axis=1)
global_panel_df.to_csv(global_panel_path, index=False) # overwrite the old file...


# Get Calling Commands
def build_calling_command(run_name, sample_name, patient_name, panel_name):
    command = f'Rscript /software/gitRepos/mrd_core/RaDaR_packages/calling_reporting/inst/bin/call_wrapper.R --metadata /shared/processed_runs/{run_name}/panel/sample/metadata/{sample_name}_calling_metadata.csv --panel /shared/processed_runs/{run_name}/panel/sample/panel/{panel_name}_calling_panel.csv  --rates /software/bundle/rates/rates_cpl_novaseq_v1.csv --thresholds /software/bundle/calling_thresholds/LR_thresholds_novaseq_cpl_v1.csv --amplicon /shared/processed_runs/{run_name}/panel_qc/bam2counts_panels/{panel_name}.bed --output /shared/processed_runs/{run_name}/calling_output/{patient_name}/{sample_name} --processed-folder /shared/processed_runs/{run_name}/panel_qc'

def get_calling_commands(panels, run_name):
    metadata_file = os.path.join('/shared/processed_runs', run_name, 'userInput/metadata.csv')
    #project.id,derived.sample.lims.name,derived.sample.lims.id,submitted.sample.lims.id,specimen.id,patient.id,cancer.indication,therapy.name,gender,collection.date,timepoint,specimen.type,extracted.plasma.volume.ml,probe.long,amount.copies.long,probe.short,amount.copies.short,recovery.percent,input.amount.copies.short,extracted.dna.id,replication,panel,barcode,acceptance.discrepancy.code,assay.version,dpcr.qc,tapestation.qc,qubit.qc,library.pooling.qc,accession.id,technical.control,product.region,source.id,product.type,molecularWeightPercent,lane
    metadata_df = pd.read_csv(metadata_file, header=10)
    metadata_df = metadata_df[['patient.id','specimen.id','panels']].drop_duplicates()
    metadata_df = metadata_df[metadata_df['panels'].apply(lambda x: x in panels)]
    commands = []
    for index, row in metadata_df.iterrows():
        patient_name = row['patient.id']
        sample_name = row['specimen.id']
        panel_name = row['panels']
        if panel_name in panels:
            commands.append(build_calling_command(run_name, sample_name, patient_name, panel_name))

    return commands


commands.append('# Calling Commands\n')
commands += get_calling_commands(panels, args.run)

commands.append('# Prepare CSVs Command\n')
commands.append(f'python RaDaR_PrepareCSVsForReport_Manual -r {args.run}')


commands.append('# Make CSVs Command\n')
commands.append(f'python RaDaR_MakeCSVsForLoaders.py -r {args.run}')

commands.append('# Report Command\n')
commands.append(f'Rscript /software/gitRepos/mrd_core/RaDaR_packages/calling_reporting/inst/bin/calling_report_wrapper.R --metadata /shared/processed_runs/{args.run}/panel/global/metadata_for_calling_report.csv --panel /shared/processed_runs/{args.run}/panel/global/reports_panel.csv --results /shared/processed_runs/{args.run}/reports/calling_report/inputs_from_calling/results.csv --variants /shared/processed_runs/{args.run}/reports/calling_report/inputs_from_calling/variants.csv --output-dir /shared/processed_runs/{args.run}/reports/calling_report/output/ --qc-report /shared/processed_runs/{args.run}/reports/QC_report/output/sample_reports.csv --radar-version /shared/processed_runs/{args.run}/versions.txt')


