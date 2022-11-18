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
import ontologymapper
import reporter
import writer

from config import COLNAMES
from config import CONFIG

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


def execute_cmd(command):
    subprocess.call(command, shell=True)


def format_output_directory(directory):
    if not directory:
        return os.getcwd()
    # If the directory string is not root and ends with a forward slash, remove it
    elif directory != "/" and directory[-1] == "/":
        return f"{directory[:-1]}"
    else:
        return directory


def load_and_process_mutational_signatures(handle, folder, label, dbs, tumor_type, calculate: bool = False, plot: bool = False):
    if calculate:
        features.CosmicSignatures.calculate_contributions(handle, folder, label)
    contexts = features.CosmicSignatures.import_context(folder, label)
    signatures = features.CosmicSignatures.import_feature(folder, label)
    return process_mutational_signatures(contexts, signatures, folder, label, dbs, tumor_type, plot)


def plot_mutational_signatures(series, folder, label):
    for drawing_function, output_name in [
        (illustrator.Signatures.generate_context_plot, 'counts'),
        (illustrator.Signatures.generate_context_plot_normalized, 'normalized'),
    ]:
        figure = drawing_function(series)
        writer.Illustrations.write(figure, folder, label, f"sigs.tricontext.{output_name}.png")


def plot_preclinical_efficacy(dictionary, folder, label):
    for index_value, nested_dictionary in dictionary.items():
        for therapy_name, therapy_dictionary in nested_dictionary.items():
            figure = therapy_dictionary['figure']
            figure_name = therapy_dictionary['figure_name']
            writer.Illustrations.write(figure, folder, label, f"{figure_name}.png")


def process_mutational_signatures(contexts, signatures, folder, label, dbs, ontology_code, plot: bool = False):
    if contexts is not None and plot:
        plot_mutational_signatures(contexts, folder, label)
    annotated = annotator.Annotator.annotate_almanac(signatures, dbs, ontology_code, label, folder)
    evaluated = evaluator.Evaluator.evaluate_almanac(annotated)
    return evaluated


def process_preclinical_efficacy(dbs, dataframe, folder, label, plot: bool = False):
    efficacy_dictionary = investigator.SensitivityDictionary.create(dbs, dataframe)
    if plot:
        plot_preclinical_efficacy(efficacy_dictionary, folder, label)
    efficacy_summary = investigator.SummaryDataFrame.create(efficacy_dictionary, dataframe, label)
    return efficacy_dictionary, efficacy_summary


def main(patient, inputs, output_folder):
    dbs = datasources.Datasources.generate_db_dict(CONFIG)
    output_folder = format_output_directory(output_folder)
    if output_folder != "":
        execute_cmd(f"mkdir -p {output_folder}")

    string_id = patient[patient_id]

    mapped_ontology = ontologymapper.OntologyMapper.map(dbs, patient[tumor_type])
    patient[ontology] = mapped_ontology[ontology]
    patient[code] = mapped_ontology[code]

    df_snv, df_snv_reject = features.MAFSomatic.import_feature(inputs[snv_handle])
    df_indel, df_indel_reject = features.MAFSomatic.import_feature(inputs[indel_handle])
    df_cnv, df_cnv_reject = features.CopyNumber.import_feature(inputs[called_cn_handle], inputs[cnv_handle])
    df_fusion, df_fusion_reject = features.Fusion.import_feature(inputs[fusion_handle])

    somatic_variants = pd.concat([df_snv, df_indel, df_cnv, df_fusion], ignore_index=True)
    somatic_filtered = pd.concat([df_snv_reject, df_indel_reject, df_cnv_reject, df_fusion_reject], ignore_index=True)

    germline_variants, germline_reject = features.MAFGermline.import_feature(inputs[germline_handle])

    if not somatic_variants.empty:
        annotated_somatic = annotator.Annotator.annotate_somatic(somatic_variants, dbs, patient[code], string_id, output_folder)
        evaluated_somatic = evaluator.Evaluator.evaluate_somatic(annotated_somatic)

        validation_accept, validation_reject = features.MAFValidation.import_feature(inputs[validation_handle])
        if not validation_accept.empty:
            evaluated_somatic = annotator.OverlapValidation.append_validation(evaluated_somatic, validation_accept)
            illustrator.ValidationOverlap.generate_dna_rna_plot(evaluated_somatic, string_id, output_folder)
    else:
        evaluated_somatic = features.Features.create_empty_dataframe()

    if not germline_variants.empty:
        annotated_germline = annotator.Annotator.annotate_germline(germline_variants, dbs, patient[code], string_id, output_folder)
        evaluated_germline = evaluator.Evaluator.evaluate_germline(annotated_germline)
    else:
        evaluated_germline = features.Features.create_empty_dataframe()

    evaluated_somatic = annotator.OverlapSomaticGermline.append_germline_hits(evaluated_somatic, evaluated_germline)
    integrated = evaluator.Integrative.evaluate(evaluated_somatic, evaluated_germline, dbs, feature_types)

    somatic_burden = features.BurdenReader.import_feature(inputs[bases_covered_handle], patient, somatic_variants, dbs)

    patient_wgd = features.Aneuploidy.summarize(patient[wgd])
    patient_ms_status = features.MicrosatelliteReader.summarize(patient[ms_status])
    patient[ms_status] = features.MicrosatelliteReader.map_status(patient[ms_status])

    annotated_burden = annotator.Annotator.annotate_almanac(somatic_burden, dbs, patient[code], string_id, output_folder)
    annotated_wgd = annotator.Annotator.annotate_almanac(patient_wgd, dbs, patient[code], string_id, output_folder)
    annotated_ms_status = annotator.Annotator.annotate_almanac(patient_ms_status, dbs, patient[code], string_id, output_folder)

    evaluated_burden = evaluator.Evaluator.evaluate_almanac(annotated_burden)
    evaluated_wgd = evaluator.Evaluator.evaluate_almanac(annotated_wgd)
    evaluated_ms_variants = evaluator.Microsatellite.evaluate_variants(evaluated_somatic, evaluated_germline)
    evaluated_ms_status = evaluator.Microsatellite.evaluate_status(annotated_ms_status, evaluated_ms_variants)

    plot_signatures = TOGGLE_FEATURES.get('plot_mutational_signatures')
    evaluated_mutational_signatures = load_and_process_mutational_signatures(
        handle=inputs[snv_handle],
        folder=output_folder,
        label=string_id,
        dbs=dbs,
        tumor_type=code,
        plot=plot_signatures
    )

    actionable = evaluator.Actionable.evaluate(
        evaluated_somatic,
        evaluated_germline,
        evaluated_ms_variants,
        evaluated_ms_status,
        evaluated_burden,
        evaluated_mutational_signatures,
        evaluated_wgd
    )

    strategies = evaluator.Strategies.report_therapy_strategies(actionable)

    efficacy_summary = investigator.SummaryDataFrame.create_empty_dataframe()
    efficacy_dictionary = {}
    cell_lines_dictionary = {}
    preclinical_efficacy_on = TOGGLE_FEATURES.getboolean('calculate_preclinical_efficacy')
    matchmaker_results = matchmaker.Matchmaker.create_empty_output()
    model_similarity_on = TOGGLE_FEATURES.getboolean('calculate_model_similarity')
    if preclinical_efficacy_on or model_similarity_on:
        dbs_preclinical = datasources.Preclinical.import_dbs()
        cell_lines_dictionary = dbs_preclinical['dictionary']
        if preclinical_efficacy_on:
            plot_preclinical = TOGGLE_FEATURES.getboolean('plot_preclinical_efficacy')
            efficacy_results = process_preclinical_efficacy(dbs_preclinical, actionable, output_folder, string_id, plot=plot_preclinical)
            efficacy_dictionary = efficacy_results[0]
            efficacy_summary = efficacy_results[1]
            actionable = annotator.PreclinicalEfficacy.annotate(actionable, efficacy_results[1])
        if model_similarity_on:
            matchmaker_results = matchmaker.Matchmaker.compare(dbs, dbs_preclinical, evaluated_somatic, string_id)

    writer.Actionable.write(actionable, string_id, output_folder)
    writer.GermlineACMG.write(evaluated_germline, string_id, output_folder)
    writer.GermlineCancer.write(evaluated_germline, string_id, output_folder)
    writer.GermlineHereditary.write(evaluated_germline, string_id, output_folder)
    writer.Integrated.write(integrated, string_id, output_folder)
    writer.MSI.write(evaluated_ms_variants, string_id, output_folder)
    writer.MutationalBurden.write(evaluated_burden, string_id, output_folder)
    writer.SomaticScored.write(evaluated_somatic, string_id, output_folder)
    writer.SomaticFiltered.write(somatic_filtered, string_id, output_folder)
    writer.Strategies.write(strategies, string_id, output_folder)
    writer.PreclinicalEfficacy.write(efficacy_summary, string_id, output_folder)
    writer.PreclinicalMatchmaking.write(matchmaker_results, string_id, output_folder)

    if TOGGLE_FEATURES.getboolean('generate_actionability_report'):
        report_dictionary = reporter.Reporter.generate_dictionary(evaluated_somatic, patient)
        reporter.generate_report(
            actionable,
            report_dictionary,
            output_directory=output_folder
        )


if __name__ == "__main__":
    start_time = time.time()

    arg_parser = argparse.ArgumentParser(prog='Molecular Oncology Almanac',
                                         description='A clinical interpretation algorithm for cancer genomics.')
    arg_parser.add_argument('--patient_id',
                            help='patient id label',
                            required=True)
    arg_parser.add_argument('--tumor_type',
                            default='Unknown',
                            help='reported tumor type')
    arg_parser.add_argument('--stage',
                            default='Unknown',
                            help='disease stage')
    arg_parser.add_argument('--snv_handle',
                            default='',
                            help='handle for SNV MAF')
    arg_parser.add_argument('--indel_handle',
                            default='',
                            help='handle for InDel MAF')
    arg_parser.add_argument('--bases_covered_handle',
                            default='',
                            help='handle for a text file which contains the numeric number of somatic bases')
    arg_parser.add_argument('--called_cn_handle',
                            default='',
                            help='handle for called copy number alterations file, used over --cnv_handle')
    arg_parser.add_argument('--cnv_handle',
                            default='',
                            help='handle for annotated seg file')
    arg_parser.add_argument('--fusion_handle',
                            default='',
                            help='handle for STAR Fusion output, .final.abridged')
    arg_parser.add_argument('--germline_handle',
                            default='',
                            help='handle for Germline MAF')
    arg_parser.add_argument('--validation_handle',
                            default='',
                            help='handle for SNV MAF called from validation sequencing')
    arg_parser.add_argument('--description',
                            default='',
                            help='description of patient')
    arg_parser.add_argument('--ms_status',
                            default='unk',
                            choices=['msih', 'msil', 'mss', 'unk'],
                            help='microsatellite instability status')
    arg_parser.add_argument('--purity',
                            default='Unknown',
                            help='Tumor purity')
    arg_parser.add_argument('--ploidy',
                            default='Unknown',
                            help='Tumor ploidy')
    arg_parser.add_argument('--wgd',
                            action='store_true',
                            help='Specify the occurrence of whole genome duplication')
    arg_parser.add_argument('--disable_matchmaking',
                            action='store_true',
                            help='Disable matchmaking in report')
    arg_parser.add_argument('--output_directory',
                            default=None,
                            help='Output directory for generated files')
    args = arg_parser.parse_args()

    patient_dict = {
        patient_id: args.patient_id,
        tumor_type: args.tumor_type,
        stage: args.stage,
        description: args.description,
        purity: args.purity,
        ploidy: args.ploidy,
        ms_status: args.ms_status,
        wgd: args.wgd
    }

    inputs_dict = {
        snv_handle: args.snv_handle,
        indel_handle: args.indel_handle,
        bases_covered_handle: args.bases_covered_handle,
        cnv_handle: args.cnv_handle,
        called_cn_handle: args.called_cn_handle,
        fusion_handle: args.fusion_handle,
        germline_handle: args.germline_handle,
        validation_handle: args.validation_handle,
        disable_matchmaking: args.disable_matchmaking
    }

    output_directory = args.output_directory if args.output_directory else os.getcwd()

    main(patient_dict, inputs_dict, output_directory)

    end_time = time.time()
    time_statement = "Molecular Oncology Almanac runtime: %s seconds" % round((end_time - start_time), 4)
    print(time_statement)
