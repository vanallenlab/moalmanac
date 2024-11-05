import time
import argparse
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
import moalmanac
import ontologymapper
import reporter
import writer

from config import COLNAMES
from reader import Ini

snv_handle = 'snv_handle'
indel_handle = 'indel_handle'
bases_covered_handle = 'bases_covered_handle'
cnv_handle = 'cnv_handle'
called_cn_handle = 'called_cn_handle'
fusion_handle = 'fusion_handle'
germline_handle = 'germline_handle'
validation_handle = 'validation_handle'
disable_matchmaking = 'disable_matchmaking'

snv_input = 'snv_input'
indel_input = 'indel_input'
seg_input = 'seg_input'
called_cn_input = 'called_cn_input'
fusion_input = 'fusion_input'
germline_input = 'germline_input'

patient_section = 'patient'
patient_id = COLNAMES[patient_section]['patient_id']
tumor_type = COLNAMES[patient_section]['tumor_type']
stage = COLNAMES[patient_section]['stage']
description = COLNAMES[patient_section]['description']
purity = COLNAMES[patient_section]['purity']
ploidy = COLNAMES[patient_section]['ploidy']
wgd = COLNAMES[patient_section]['wgd']
ms_status = COLNAMES[patient_section]['ms_status']

oncotree_section = 'oncotree'
ontology = COLNAMES[oncotree_section]['ontology']
code = COLNAMES[oncotree_section]['code']

feature_type_section = 'feature_types'
feature_type_mut = CONFIG[feature_type_section]['mut']
feature_type_germline = CONFIG[feature_type_section]['germline']
feature_type_cna = CONFIG[feature_type_section]['cna']
feature_type_fusion = CONFIG[feature_type_section]['fusion']
feature_type_burden = CONFIG[feature_type_section]['burden']
feature_type_signature = CONFIG[feature_type_section]['signature']
feature_type_microsatellite = CONFIG[feature_type_section]['microsatellite']
feature_type_aneuploidy = CONFIG[feature_type_section]['aneuploidy']
feature_types = {
    'mutation': feature_type_mut,
    'germline': feature_type_germline,
    'copynumber': feature_type_cna,
    'fusion': feature_type_fusion,
    'burden': feature_type_burden,
    'signature': feature_type_signature,
    'microsatellite': feature_type_microsatellite,
    'aneuploidy': feature_type_aneuploidy
}

generate_illustrations = 'generate_illustrations'
TOGGLE_FEATURES = CONFIG['function_toggle']


def subset_by_feature_type(dataframe):
    feature_type = features.Simple.feature_type
    somatic_types = [feature_type_mut, feature_type_cna, feature_type_fusion]
    germline_types = [feature_type_germline]

    somatic = dataframe[dataframe[feature_type].isin(somatic_types)].reset_index(drop=True)
    germline = dataframe[dataframe[feature_type].isin(germline_types)].reset_index(drop=True)
    return somatic, germline


def main(patient, input_file, output_folder, config, dbs, dbs_preclinical=None):
    metadata_dictionary = moalmanac.create_metadata_dictionary(patient)

    output_folder = moalmanac.format_output_directory(output_folder)
    if output_folder != "":
        moalmanac.execute_cmd(f"mkdir -p {output_folder}")

    string_id = metadata_dictionary[patient_id]

    mapped_ontology = ontologymapper.OntologyMapper.map(dbs, metadata_dictionary[tumor_type])
    metadata_dictionary[ontology] = mapped_ontology[ontology]
    metadata_dictionary[code] = mapped_ontology[code]

    alterations = features.Simple.import_feature(input_file)
    annotated_alterations = annotator.Annotator.annotate_simple(alterations, dbs, patient[code], config=config)
    evaluated_alterations = evaluator.Evaluator.evaluate_somatic(annotated_alterations)

    evaluated_somatic, evaluated_germline = subset_by_feature_type(evaluated_alterations)
    evaluated_somatic = annotator.OverlapSomaticGermline.append_germline_hits(evaluated_somatic, evaluated_germline)
    integrated = evaluator.Integrative.evaluate(evaluated_somatic, evaluated_germline, dbs, feature_types)

    somatic_burden = features.BurdenReader.import_feature(
        handle='',
        patient=metadata_dictionary,
        variants=evaluated_somatic,
        dbs=dbs,
        config=config
    )
    patient_wgd = features.Aneuploidy.summarize(patient[wgd], config=config)
    patient_ms_status = features.MicrosatelliteReader.summarize(patient[ms_status], config=config)
    patient[ms_status] = features.MicrosatelliteReader.map_status(patient[ms_status])

    annotated_burden = annotator.Annotator.annotate_almanac(somatic_burden, dbs, patient[code], config=config)
    annotated_wgd = annotator.Annotator.annotate_almanac(patient_wgd, dbs, patient[code], config=config)
    annotated_ms_status = annotator.Annotator.annotate_almanac(patient_ms_status, dbs, patient[code], config=config)

    evaluated_burden = evaluator.Evaluator.evaluate_almanac(annotated_burden)
    evaluated_wgd = evaluator.Evaluator.evaluate_almanac(annotated_wgd)
    evaluated_ms_variants = evaluator.Microsatellite.evaluate_variants(evaluated_somatic, evaluated_germline)
    evaluated_ms_status = evaluator.Microsatellite.evaluate_status(annotated_ms_status, evaluated_ms_variants)

    actionable = evaluator.Actionable.evaluate(
        somatic=evaluated_somatic,
        germline=evaluated_germline,
        ms_variants=evaluated_ms_variants,
        ms_status=evaluated_ms_status,
        burden=evaluated_burden,
        signatures=features.Features.create_empty_dataframe(),
        wgd=evaluated_wgd,
        config=config
    )

    strategies = evaluator.Strategies.report_therapy_strategies(actionable)

    function_toggle = config['function_toggle']
    efficacy_summary = investigator.SummaryDataFrame.create_empty_dataframe()
    preclinical_efficacy_on = TOGGLE_FEATURES.getboolean('calculate_preclinical_efficacy')

    # The input argument --disable_matchmaking will be removed in the next non-backwards compatible release
    model_similarity_on = TOGGLE_FEATURES.getboolean('calculate_model_similarity')
    similarity_results = matchmaker.Matchmaker.create_empty_output()
    similarity_summary = {}

    if preclinical_efficacy_on or model_similarity_on:
        dbs_preclinical = datasources.Preclinical.import_dbs()
        cell_lines_dictionary = dbs_preclinical['dictionary']
        if preclinical_efficacy_on:
            plot_preclinical = TOGGLE_FEATURES.getboolean('plot_preclinical_efficacy')
            efficacy_results = moalmanac.process_preclinical_efficacy(
                dbs=dbs_preclinical,
                dataframe=actionable,
                folder=output_folder,
                label=string_id,
                config=config,
                plot=plot_preclinical
            )
            efficacy_dictionary = efficacy_results[0]
            efficacy_summary = efficacy_results[1]

            actionable = annotator.PreclinicalEfficacy.annotate(
                actionable,
                efficacy_summary,
                efficacy_dictionary,
                append_lookup=TOGGLE_FEATURES.getboolean('include_preclinical_efficacy_in_actionability_report')
            )
        if model_similarity_on:
            similarity_results = matchmaker.Matchmaker.compare(
                dbs=dbs,
                dbs_preclinical=dbs_preclinical,
                somatic=evaluated_somatic,
                case_sample_id=string_id,
                config=config
            )
            similarity_summary = matchmaker.Report.create_report_dictionary(
                similarity_results,
                cell_lines_dictionary
            )

    writer.Actionable.write(actionable, string_id, output_folder)
    writer.GermlineACMG.write(evaluated_germline, string_id, output_folder)
    writer.GermlineCancer.write(evaluated_germline, string_id, output_folder)
    writer.GermlineHereditary.write(evaluated_germline, string_id, output_folder)
    writer.Integrated.write(integrated, string_id, output_folder)
    writer.MSI.write(evaluated_ms_variants, string_id, output_folder)
    writer.SomaticScored.write(evaluated_somatic, string_id, output_folder)
    writer.SomaticFiltered.write(features.Features.create_empty_dataframe(), string_id, output_folder)
    writer.Strategies.write(strategies, string_id, output_folder)
    writer.PreclinicalEfficacy.write(efficacy_summary, string_id, output_folder)
    writer.PreclinicalMatchmaking.write(similarity_results, string_id, output_folder)

    if function_toggle.getboolean('generate_actionability_report'):
        report_dictionary = reporter.Reporter.generate_dictionary(evaluated_somatic, metadata_dictionary)

        include_similarity = function_toggle.getboolean('include_model_similarity_in_actionability_report')
        reporter.Reporter.generate_actionability_report(
            actionable=actionable,
            report_dictionary=report_dictionary,
            config=config,
            similarity=similarity_summary if include_similarity else None,
            output_directory=output_folder
        )


if __name__ == "__main__":
    start_time = time.time()

    arg_parser = argparse.ArgumentParser(prog='Molecular Oncology Almanac using simplified input',
                                         description='Annotates only using the Molecular Oncology Almanac database')
    arg_parser.add_argument(
        '--patient_id',
        help='patient id label',
        required=True
    )
    arg_parser.add_argument(
        '--stage',
        default='Unknown',
        help='disease stage'
    )
    arg_parser.add_argument(
        '--tumor_type',
        default='Unknown',
        help='reported tumor type'
    )
    arg_parser.add_argument(
        '--input',
        help='Tab delimited file of observed alterations'
    )
    arg_parser.add_argument(
        '--ms_status',
        default='unk',
        choices=['msih', 'msil', 'mss', 'unk'],
        help='microsatellite instability status'
    )
    arg_parser.add_argument(
        '--purity',
        default='Unknown',
        help='Tumor purity'
    )
    arg_parser.add_argument(
        '--ploidy',
        default='Unknown',
        help='Tumor ploidy'
    )
    arg_parser.add_argument(
        '--wgd',
        action='store_true',
        help='Specify the occurrence of whole genome duplication'
    )
    arg_parser.add_argument(
        '--output_directory',
        default=None,
        help='Output directory for generated files'
    )
    arg_parser.add_argument(
        '--config', '-c',
        required=True,
        help='ini file that contains configuration details'
    )
    arg_parser.add_argument(
        '--dbs',
        required=True,
        help='ini file that contains database paths '
    )
    arg_parser.add_argument(
        '--preclinical-dbs',
        required=False,
        help='ini file that contains preclinical file paths'
    )
    args = arg_parser.parse_args()

    patient_dict = {
        patient_id: args.patient_id,
        tumor_type: args.tumor_type,
        stage: "",
        description: "",
        purity: args.purity,
        ploidy: args.ploidy,
        ms_status: args.ms_status,
        wgd: args.wgd
    }

    output_directory = args.output_directory if args.output_directory else os.getcwd()

    config_ini = Ini.read(args.config, extended_interpolation=False, convert_to_dictionary=False)

    db_paths = Ini.read(args.dbs, extended_interpolation=True, convert_to_dictionary=True)
    if args.preclinical_dbs:
        preclinical_db_paths = Ini.read(args.preclinical_dbs, extended_interpolation=True, convert_to_dictionary=True)
    else:
        preclinical_db_paths = None

    main(
        patient=patient_dict,
        input_file=args.input,
        output_folder=output_directory,
        config=config_ini,
        dbs=db_paths['databases'],
        dbs_preclinical=preclinical_db_paths['preclinical'] if preclinical_db_paths else None
    )

    end_time = time.time()
    time_statement = "Molecular Oncology Almanac runtime: %s seconds" % round((end_time - start_time), 4)
    print(time_statement)
