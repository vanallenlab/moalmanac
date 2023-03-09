import moalmanac
import os
import time
import subprocess

patient_dict = {
    'patient_id': 'example',
    'reported_tumor_type': 'SKCM',
    'stage': 'Metastatic',
    'description': 'Test patient for development runs',
    'purity': 0.85,
    'ploidy': 4.02,
    'WGD': True,
    'microsatellite_status': 'msih'
}

empty_dict = {
    'snv_handle': '../example_data/example_patient.capture.somatic.snvs.maf',
    'indel_handle': '',
    'bases_covered_handle': '',
    'called_cn_handle': '',
    'cnv_handle': '',
    'fusion_handle': '',
    'germline_handle': '',
    'validation_handle': '',
    'disable_matchmaking': False
}

example_dict = {
    'snv_handle': '../example_data/example_patient.capture.somatic.snvs.maf',
    'indel_handle': '../example_data/example_patient.capture.somatic.indels.maf',
    'bases_covered_handle': '../example_data/example_patient.capture.somatic.coverage.txt',
    'called_cn_handle': '../example_data/example_patient.capture.somatic.called.cna.txt',
    'cnv_handle': '../example_data/example_patient.capture.somatic.seg.annotated',
    'fusion_handle': '../example_data/example_patient.rna.star.fusions.txt',
    'germline_handle': '../example_data/example_patient.capture.germline.maf',
    'validation_handle': '../example_data/example_patient.rna.somatic.snvs.maf',
    'disable_matchmaking': False
}


def execute_cmd(command):
    subprocess.call(command, shell=True)


output_directory = "example"
if output_directory != "":
    cmd = f"mkdir -p {output_directory}"
    execute_cmd(cmd)
else:
    output_directory = os.getcwd()

start_time = time.time()
moalmanac.main(patient_dict, example_dict, output_directory)
end_time = time.time()

time_statement = "Molecular Oncology Almanac runtime: %s seconds" % round((end_time - start_time), 4)
print(time_statement)
