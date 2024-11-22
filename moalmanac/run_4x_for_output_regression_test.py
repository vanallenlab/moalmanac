import moalmanac
import os
import subprocess
import time

from datetime import date

from reader import Ini

pre_or_post_code_changes = "pre"
#pre_or_post_code_changes = "post"

metadata_dictionary = {
    'patient_id': 'example',
    'reported_tumor_type': 'MEL',
    'stage': 'Metastatic',
    'description': 'Test patient for development runs',
    'purity': 0.85,
    'ploidy': 4.02,
    'WGD': True,
    'microsatellite_status': 'msih'
}

input_dictionary_empty = {
    'snv_handle': '',
    'indel_handle': '',
    'bases_covered_handle': '',
    'called_cn_handle': '',
    'cnv_handle': '',
    'fusion_handle': '',
    'germline_handle': '',
    'validation_handle': '',
    'mutational_signatures_path': '',
    'disable_matchmaking': False
}

input_dictionary = {
    'snv_handle': '../example_data/example_patient.capture.somatic.snvs.maf',
    'indel_handle': '../example_data/example_patient.capture.somatic.indels.maf',
    'bases_covered_handle': '../example_data/example_patient.capture.somatic.coverage.txt',
    'called_cn_handle': '../example_data/example_patient.capture.somatic.called.cna.txt',
    'cnv_handle': '../example_data/example_patient.capture.somatic.seg.annotated',
    'fusion_handle': '../example_data/example_patient.rna.star.fusions.txt',
    'germline_handle': '../example_data/example_patient.capture.germline.maf',
    'validation_handle': '../example_data/example_patient.rna.somatic.snvs.maf',
    'mutational_signatures_path': '../example_data/example_patient.capture.sbs_contributions.txt',
    'disable_matchmaking': False
}

config_ini_path = "config.ini"
config_ini = Ini.read(config_ini_path, extended_interpolation=False, convert_to_dictionary=False)

dbs_ini_path = "annotation-databases.ini"
db_paths = Ini.read(dbs_ini_path, extended_interpolation=True, convert_to_dictionary=True)
db_paths = db_paths['paths']

dbs_preclinical_ini_path = "preclinical-databases.ini"
preclinical_db_paths = Ini.read(dbs_preclinical_ini_path, extended_interpolation=True, convert_to_dictionary=True)
preclinical_db_paths = preclinical_db_paths['paths']


def execute_cmd(command):
    subprocess.call(command, shell=True)


today = date.today().isoformat()

"""
First execution, all example inputs with all settings within the function toggle enabled
"""
for key, value in config_ini.items('function_toggle'):
    config_ini['function_toggle'][key] = 'on'

output_directory = f"{today}-example-outputs-full-all-enabled-{pre_or_post_code_changes}"
if output_directory != "":
    cmd = f"mkdir -p {output_directory}"
    execute_cmd(cmd)
else:
    output_directory = os.getcwd()

start_time = time.time()

moalmanac.main(
    patient=metadata_dictionary,
    inputs=input_dictionary,
    output_folder=output_directory,
    config=config_ini,
    dbs=db_paths,
    dbs_preclinical=preclinical_db_paths
)

end_time = time.time()

time_statement = "Molecular Oncology Almanac runtime: %s seconds" % round((end_time - start_time), 4)
print(time_statement)

"""
Second execution, all empty example inputs with all settings within the function toggle enabled
"""
for key, value in config_ini.items('function_toggle'):
    config_ini['function_toggle'][key] = 'on'

output_directory = f"{today}-example-outputs-empty-all-enabled-{pre_or_post_code_changes}"
if output_directory != "":
    cmd = f"mkdir -p {output_directory}"
    execute_cmd(cmd)
else:
    output_directory = os.getcwd()

start_time = time.time()
moalmanac.main(
    patient=metadata_dictionary,
    inputs=input_dictionary_empty,
    output_folder=output_directory,
    config=config_ini,
    dbs=db_paths,
    dbs_preclinical=preclinical_db_paths
)
end_time = time.time()

time_statement = "Molecular Oncology Almanac runtime: %s seconds" % round((end_time - start_time), 4)
print(time_statement)

"""
Third execution, all example inputs with all settings within the function toggle disabled
"""
for key, value in config_ini.items('function_toggle'):
    config_ini['function_toggle'][key] = 'off'

output_directory = f"{today}-example-outputs-full-all-disabled-{pre_or_post_code_changes}"
if output_directory != "":
    cmd = f"mkdir -p {output_directory}"
    execute_cmd(cmd)
else:
    output_directory = os.getcwd()

start_time = time.time()

moalmanac.main(
    patient=metadata_dictionary,
    inputs=input_dictionary,
    output_folder=output_directory,
    config=config_ini,
    dbs=db_paths,
    dbs_preclinical=preclinical_db_paths
)

end_time = time.time()

time_statement = "Molecular Oncology Almanac runtime: %s seconds" % round((end_time - start_time), 4)
print(time_statement)

"""
Fourth execution, all empty example inputs with all settings within the function toggle disabled
"""
for key, value in config_ini.items('function_toggle'):
    config_ini['function_toggle'][key] = 'off'

output_directory = f"{today}-example-outputs-empty-all-disabled-{pre_or_post_code_changes}"
if output_directory != "":
    cmd = f"mkdir -p {output_directory}"
    execute_cmd(cmd)
else:
    output_directory = os.getcwd()

start_time = time.time()
moalmanac.main(
    patient=metadata_dictionary,
    inputs=input_dictionary_empty,
    output_folder=output_directory,
    config=config_ini,
    dbs=db_paths,
    dbs_preclinical=preclinical_db_paths
)
end_time = time.time()

time_statement = "Molecular Oncology Almanac runtime: %s seconds" % round((end_time - start_time), 4)
print(time_statement)
