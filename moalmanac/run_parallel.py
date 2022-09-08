import argparse
import multiprocessing
import numpy as np
import os
import pandas as pd
import subprocess

import annotator
import datasources
import features
import evaluator
import illustrator
import investigator
import matchmaker
import ontologymapper
import reporter
import writer

from config import COLNAMES
from config import CONFIG

oncotree_section = 'oncotree'
ontology = COLNAMES[oncotree_section]['ontology']
code = COLNAMES[oncotree_section]['code']


def format_burden(dictionary, variants):
    df = features.BurdenReader.create_burden_series(dictionary, np.nan)
    df[features.Features.feature_type] = features.BurdenReader.feature_type

    mutations = variants[variants[features.Features.feature_type] == features.BurdenReader.feature_type_mutations].shape[0]
    df[features.BurdenReader.n_nonsyn_mutations] = mutations
    df[features.BurdenReader.mutational_burden] = np.nan
    df[features.BurdenReader.percentile_tcga] = np.nan
    df[features.BurdenReader.percentile_tcga_tissue] = np.nan

    burden_boolean = dictionary['high_tmb']
    df[features.BurdenReader.high_burden_boolean] = burden_boolean
    df[features.Features.feature] = features.BurdenReader.evaluate_high_burden_boolean(burden_boolean)
    df[features.Features.alt] = ''
    return df.to_frame().T


def format_copy_number(df):
    feature_type = CONFIG['feature_types']['cna']
    column_map = features.CalledCNAReader.create_column_map(COLNAMES['called_cn_input'])
    column_map['profile_name'] = 'profile_name'
    df = df.loc[:, column_map.keys()].rename(columns=column_map)
    df = features.Features.preallocate_missing_columns(df)
    if not df.empty:
        df[features.Features.feature_type] = features.Features.annotate_feature_type(feature_type, df.index)
        df[features.Features.feature] = features.CalledCNAReader.format_cn_gene(df[features.Features.feature])

        idx = features.CalledCNAReader.filter_calls(df[features.Features.alt_type])
        seg_accept = df[idx]
        seg_reject = df[~idx]
    else:
        seg_accept = features.Features.create_empty_dataframe()
        seg_reject = features.Features.create_empty_dataframe()
    return seg_accept, seg_reject


def format_variants(df):
    feature_type = CONFIG['feature_types']['mut']
    column_map = features.MutationReader.create_colmap(COLNAMES['snv_input'])
    column_map['profile_name'] = 'profile_name'
    df = df.loc[:, column_map.keys()].rename(columns=column_map)
    df = features.Features.preallocate_missing_columns(df)
    if not df.empty:
        df[features.Features.feature_type] = features.Features.annotate_feature_type(feature_type, df.index)
        df[features.Features.alt_count] = features.CoverageMetrics.format_coverage_col(df[features.Features.alt_count])
        df[features.Features.ref_count] = features.CoverageMetrics.format_coverage_col(df[features.Features.ref_count])
        df[features.Features.coverage] = features.CoverageMetrics.append_coverage(df[features.Features.alt_count], df[features.Features.ref_count])
        df[features.Features.tumor_f] = features.CoverageMetrics.append_tumor_f(df[features.Features.alt_count], df[features.Features.coverage])
        df[features.Features.alt_type] = features.Features.rename_coding_classifications(df[features.Features.alt_type])

        coding_variants = features.Features.return_variants_coding(df)
        non_coding_variants = features.Features.return_variants_non_coding(df)
    else:
        coding_variants = features.Features.create_empty_dataframe()
        non_coding_variants = features.Features.create_empty_dataframe()
    return coding_variants, non_coding_variants


def process_profile(dbs, profile_dictionary, input_dictionary):
    mapped_ontology = ontologymapper.OntologyMapper.map(dbs, profile_dictionary['reported_tumor_type'])
    profile_dictionary[ontology] = mapped_ontology[ontology]
    profile_dictionary[code] = mapped_ontology[code]

    profile_label = profile_dictionary['patient_id']
    ontology_code = profile_dictionary[code]
    ms_status = profile_dictionary['microsatellite_status']

    somatic_accepted = input_dictionary['somatic_accepted']
    somatic_rejected = input_dictionary['somatic_rejected']
    if not somatic_accepted.empty:
        annotated_somatic = annotator.Annotator.annotate_somatic(somatic_accepted, dbs, ontology_code)
        evaluated_somatic = evaluator.Evaluator.evaluate_somatic(annotated_somatic)
    else:
        evaluated_somatic = features.Features.create_empty_dataframe()

    evaluated_burden = format_burden(profile_dictionary, somatic_accepted)
    evaluated_germline = features.Features.create_empty_dataframe()
    evaluated_ms_status = features.Features.create_empty_dataframe()
    evaluated_signatures = features.Features.create_empty_dataframe()
    evaluated_wgd = features.Features.create_empty_dataframe()

    #annotated_ms_status = annotator.Annotator.annotate_almanac(ms_status, dbs, ontology_code)
    evaluated_ms_variants = evaluator.Microsatellite.evaluate_variants(evaluated_somatic, evaluated_germline)
    #evaluated_ms_status = evaluator.Microsatellite.evaluate_status(annotated_ms_status, evaluated_ms_variants)

    actionable = evaluator.Actionable.evaluate(evaluated_somatic, evaluated_germline,
                                               evaluated_ms_variants, evaluated_ms_status,
                                               evaluated_burden, evaluated_signatures, evaluated_wgd)
    report_dictionary = reporter.Reporter.generate_dictionary(evaluated_somatic, profile_dictionary)
    version_dictionary = reporter.Reporter.generate_version_dictionary()

    writer.Actionable.write(actionable, profile_label)
    writer.SomaticScored.write(evaluated_somatic, profile_label)
    writer.SomaticFiltered.write(somatic_rejected, profile_label)

    try:
        os.mkdir(profile_label)
    except OSError as error:
        print(error)

    command = f"mv {profile_label}*.txt {profile_label}/"
    subprocess.call(command, shell=True)

    dbs_preclinical = datasources.Preclinical.import_dbs()
    efficacy_dictionary = investigator.SensitivityDictionary.create(dbs_preclinical, actionable, profile_label)
    efficacy_summary = investigator.SummaryDataFrame.create(efficacy_dictionary, actionable, profile_label)
    matchmaker_results = matchmaker.Matchmaker.compare(dbs, dbs_preclinical, evaluated_somatic, profile_label)

    reporter.Reporter.generate_report(actionable,
                                      report_dictionary,
                                      version_dictionary,
                                      efficacy_dictionary,
                                      efficacy_summary,
                                      False,
                                      matchmaker_results,
                                      dbs_preclinical['dictionary']
                                      )

    command = f"mv build/index.html {profile_label}/{profile_label}.report.html"
    subprocess.call(command, shell=True)


def wrapper(profiles, aggregate_variants=None, aggregate_copy_number_alterations=None, aggregate_fusions=None):
    dbs = datasources.Datasources.generate_db_dict(CONFIG)

    empty_dataframe = features.Features.create_empty_dataframe()
    aggregate_variants = empty_dataframe if aggregate_variants is None else aggregate_variants
    formatted_variants_accepted, formatted_variants_rejected = format_variants(aggregate_variants)

    aggregate_copy_number_alterations = empty_dataframe if aggregate_copy_number_alterations is None else aggregate_copy_number_alterations
    formatted_cn_accepted, formatted_cn_rejected = format_copy_number(aggregate_copy_number_alterations)

    somatic_accepted = pd.concat([formatted_variants_accepted, formatted_cn_accepted])
    somatic_rejected = pd.concat([formatted_variants_rejected, formatted_cn_rejected])

    processes = []
    for label, group in profiles.groupby('profile_name'):
        profile_dictionary = {
            'patient_id': label,
            'reported_tumor_type': 'SKCM',
            'microsatellite_status': 'MSS',
            'high_tmb': True,
            'stage': '',
            'description': '',
            'purity': '',
            'ploidy': '',
            'wgd': ''
        }

        input_dictionary = {
            'somatic_accepted': somatic_accepted[somatic_accepted['profile_name'].eq(label)].reset_index(drop=True),
            'somatic_rejected': somatic_rejected[somatic_rejected['profile_name'].eq(label)].reset_index(drop=True),
        }

        p = multiprocessing.Process(target=process_profile, args=(dbs, profile_dictionary, input_dictionary))
        processes.append(p)
        p.start()

    for process in processes:
        process.join()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='MOAlmanac wrapper', description="Submit patients to different cpus")
    parser.add_argument('--profiles', '-p', help='Aggregate profiles input')
    parser.add_argument('--somatic-variants', '-v', default='', help='Aggregate somatic variants input')
    parser.add_argument('--called-copy-numbers', '-c', default='', help='Aggregate called copy number alterations input')
    parser.add_argument('--fusions', '-f', default='', help='Aggregate fusions input')
    args = parser.parse_args()

    input_profiles = pd.read_csv(args.profiles, sep='\t')
    input_somatic = pd.read_csv(args.somatic_variants, sep='\t')
    input_copy_number = pd.read_csv(args.called_copy_numbers, sep='\t')

    wrapper(
        input_profiles,
        aggregate_variants=input_somatic,
        aggregate_copy_number_alterations=input_copy_number
    )
