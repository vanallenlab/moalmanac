import moalmanac
import os
import time
import subprocess

from reader import Ini

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
dbs_ini_path = "annotation-databases.ini"
dbs_preclinical_ini_path = "preclinical-databases.ini"

config_ini = Ini.read(config_ini_path, extended_interpolation=False, convert_to_dictionary=False)
db_paths = Ini.read(dbs_ini_path, extended_interpolation=False, convert_to_dictionary=True)
preclinical_db_paths = Ini.read(dbs_preclinical_ini_path, extended_interpolation=False, convert_to_dictionary=True)


def execute_cmd(command):
    subprocess.call(command, shell=True)


output_directory = "2024-07-03-config-edits"
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
    dbs=db_paths['databases'],
    dbs_preclinical=preclinical_db_paths['preclinical']
)
end_time = time.time()

time_statement = "Molecular Oncology Almanac runtime: %s seconds" % round((end_time - start_time), 4)
print(time_statement)
