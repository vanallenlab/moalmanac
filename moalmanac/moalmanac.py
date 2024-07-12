import time
import argparse
import os
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
from reader import Ini

snv_handle = 'snv_handle'
indel_handle = 'indel_handle'
bases_covered_handle = 'bases_covered_handle'
cnv_handle = 'cnv_handle'
called_cn_handle = 'called_cn_handle'
fusion_handle = 'fusion_handle'
germline_handle = 'germline_handle'
mutational_signatures_path = 'mutational_signatures_path'
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

generate_illustrations = 'generate_illustrations'


def create_metadata_dictionary(input_dictionary):
    dictionary = {}
    for key, value in input_dictionary.items():
        dictionary[key] = value
    return dictionary


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


def load_and_process_mutational_signatures(input, dbs, tumor_type, config):
    signatures = features.CosmicSignatures.import_feature(input, config)
    annotated = annotator.Annotator.annotate_almanac(signatures, dbs, tumor_type, config)
    evaluated = evaluator.Evaluator.evaluate_almanac(annotated)
    return evaluated


def plot_preclinical_efficacy(dictionary, folder, label):
    for index_value, nested_dictionary in dictionary.items():
        for therapy_name, therapy_dictionary in nested_dictionary.items():
            figure = therapy_dictionary['figure']
            figure_name = therapy_dictionary['figure_name']
            writer.Illustrations.write(figure, folder, label, f"{figure_name}.png")


def process_preclinical_efficacy(dbs, dataframe, folder, label, config, plot: bool = False):
    efficacy_dictionary = investigator.SensitivityDictionary.create(dbs, dataframe, config)
    if plot:
        plot_preclinical_efficacy(efficacy_dictionary, folder, label)
    efficacy_summary = investigator.SummaryDataFrame.create(efficacy_dictionary, dataframe, label)
    return efficacy_dictionary, efficacy_summary


def main(patient, inputs, output_folder, config, dbs, dbs_preclinical=None):
    metadata_dictionary = create_metadata_dictionary(patient)

    output_folder = format_output_directory(output_folder)
    if output_folder != "":
        execute_cmd(f"mkdir -p {output_folder}")

    string_id = metadata_dictionary[patient_id]

    mapped_ontology = ontologymapper.OntologyMapper.map(dbs, metadata_dictionary[tumor_type])
    metadata_dictionary[ontology] = mapped_ontology[ontology]
    metadata_dictionary[code] = mapped_ontology[code]

    df_snv, df_snv_reject = features.MAFSomatic.import_feature(inputs[snv_handle], config)
    df_indel, df_indel_reject = features.MAFSomatic.import_feature(inputs[indel_handle], config)
    df_cnv, df_cnv_reject = features.CopyNumber.import_feature(inputs[called_cn_handle], inputs[cnv_handle], config)
    df_fusion, df_fusion_reject = features.Fusion.import_feature(inputs[fusion_handle], config)

    accepted_variants = [df_snv, df_indel, df_cnv, df_fusion]
    filtered_variants = [df_snv_reject, df_indel_reject, df_cnv_reject, df_fusion_reject]
    somatic_variants = features.Features.concat_list_of_dataframes(accepted_variants)
    somatic_filtered = features.Features.concat_list_of_dataframes(filtered_variants)

    germline_variants, germline_reject = features.MAFGermline.import_feature(inputs[germline_handle], config)

    if not somatic_variants.empty:
        annotated_somatic = annotator.Annotator.annotate_somatic(
            df=somatic_variants,
            dbs=dbs,
            ontology=metadata_dictionary[code],
            config=config
        )
        evaluated_somatic = evaluator.Evaluator.evaluate_somatic(annotated_somatic)

        validation_accept, validation_reject = features.MAFValidation.import_feature(inputs[validation_handle], config)
        if not validation_accept.empty:
            evaluated_somatic = annotator.OverlapValidation.append_validation(
                primary=evaluated_somatic,
                validation=validation_accept,
                biomarker_type=config['feature_types']['mut'])
            illustrator.ValidationOverlap.generate_dna_rna_plot(evaluated_somatic, string_id, output_folder, config)
    else:
        evaluated_somatic = features.Features.create_empty_dataframe()

    if not germline_variants.empty:
        annotated_germline = annotator.Annotator.annotate_germline(
            germline_variants,
            dbs,
            metadata_dictionary[code],
            config=config
        )
        evaluated_germline = evaluator.Evaluator.evaluate_germline(annotated_germline)
    else:
        evaluated_germline = features.Features.create_empty_dataframe()

    evaluated_somatic = annotator.OverlapSomaticGermline.append_germline_hits(evaluated_somatic, evaluated_germline)
    integrated = evaluator.Integrative.evaluate(evaluated_somatic, evaluated_germline, dbs, config)

    somatic_burden = features.BurdenReader.import_feature(
        handle=inputs[bases_covered_handle],
        patient=metadata_dictionary,
        variants=somatic_variants,
        dbs=dbs,
        config=config
    )

    patient_wgd = features.Aneuploidy.summarize(metadata_dictionary[wgd], config)
    patient_ms_status = features.MicrosatelliteReader.summarize(metadata_dictionary[ms_status], config)
    metadata_dictionary[ms_status] = features.MicrosatelliteReader.map_status(metadata_dictionary[ms_status])

    annotated_burden = annotator.Annotator.annotate_almanac(
        df=somatic_burden,
        dbs=dbs,
        ontology=metadata_dictionary[code],
        config=config
    )
    annotated_wgd = annotator.Annotator.annotate_almanac(
        df=patient_wgd,
        dbs=dbs,
        ontology=metadata_dictionary[code],
        config=config
    )
    annotated_ms_status = annotator.Annotator.annotate_almanac(
        df=patient_ms_status,
        dbs=dbs,
        ontology=metadata_dictionary[code],
        config=config
    )

    evaluated_burden = evaluator.Evaluator.evaluate_almanac(annotated_burden)
    evaluated_wgd = evaluator.Evaluator.evaluate_almanac(annotated_wgd)
    evaluated_ms_variants = evaluator.Microsatellite.evaluate_variants(evaluated_somatic, evaluated_germline)
    evaluated_ms_status = evaluator.Microsatellite.evaluate_status(annotated_ms_status, evaluated_ms_variants)

    evaluated_mutational_signatures = load_and_process_mutational_signatures(
        input=inputs[mutational_signatures_path],
        dbs=dbs,
        tumor_type=code,
        config=config
    )

    actionable = evaluator.Actionable.evaluate(
        somatic=evaluated_somatic,
        germline=evaluated_germline,
        ms_variants=evaluated_ms_variants,
        ms_status=evaluated_ms_status,
        burden=evaluated_burden,
        signatures=evaluated_mutational_signatures,
        wgd=evaluated_wgd,
        config=config
    )

    strategies = evaluator.Strategies.report_therapy_strategies(actionable)

    function_toggle = config['function_toggle']

    efficacy_summary = investigator.SummaryDataFrame.create_empty_dataframe()
    preclinical_efficacy_on = function_toggle.getboolean('calculate_preclinical_efficacy')

    # The input argument --disable_matchmaking will be removed in the next non-backwards compatible release
    model_similarity_on = function_toggle.getboolean('calculate_model_similarity') and not inputs[disable_matchmaking]
    similarity_results = matchmaker.Matchmaker.create_empty_output()
    similarity_summary = {}

    if dbs_preclinical is not None:
        if preclinical_efficacy_on or model_similarity_on:
            dbs_preclinical = datasources.Preclinical.import_dbs(dbs_preclinical)
            cell_lines_dictionary = dbs_preclinical['dictionary']
            if preclinical_efficacy_on:
                plot_preclinical = function_toggle.getboolean('plot_preclinical_efficacy')
                efficacy_results = process_preclinical_efficacy(
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
                    append_lookup=function_toggle.getboolean('include_preclinical_efficacy_in_actionability_report')
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
    writer.MutationalBurden.write(evaluated_burden, string_id, output_folder)
    writer.SomaticScored.write(evaluated_somatic, string_id, output_folder)
    writer.SomaticFiltered.write(somatic_filtered, string_id, output_folder)
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

    arg_parser = argparse.ArgumentParser(
        prog='Molecular Oncology Almanac',
        description='A clinical interpretation algorithm for cancer genomics.'
    )
    arg_parser.add_argument(
        '--patient_id',
        help='patient id label',
        required=True
    )
    arg_parser.add_argument(
        '--description',
        default='',
        help='description of patient'
    )
    arg_parser.add_argument(
        '--tumor_type',
        default='Unknown',
        help='reported tumor type'
    )
    arg_parser.add_argument(
        '--stage',
        default='Unknown',
        help='disease stage'
    )
    arg_parser.add_argument(
        '--snv_handle',
        default='',
        help='handle for SNV MAF'
    )
    arg_parser.add_argument(
        '--indel_handle',
        default='',
        help='handle for InDel MAF'
    )
    arg_parser.add_argument(
        '--bases_covered_handle',
        default='',
        help='handle for a text file which contains the numeric number of somatic bases'
    )
    arg_parser.add_argument(
        '--called_cn_handle',
        default='',
        help='handle for called copy number alterations file, used over --cnv_handle'
    )
    arg_parser.add_argument(
        '--cnv_handle',
        default='',
        help='handle for annotated seg file'
    )
    arg_parser.add_argument(
        '--fusion_handle',
        default='',
        help='handle for STAR Fusion output, .final.abridged'
    )
    arg_parser.add_argument(
        '--germline_handle',
        default='',
        help='handle for Germline MAF'
    )
    arg_parser.add_argument(
        '--validation_handle',
        default='',
        help='handle for SNV MAF called from validation sequencing'
    )
    arg_parser.add_argument(
        '--ms_status',
        default='unk',
        choices=['msih', 'msil', 'mss', 'unk'],
        help='microsatellite instability status'
    )
    arg_parser.add_argument(
        '--mutational_signatures',
        default='',
        help='file for SBS signature contributions, version 3.4'
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
        '--disable_matchmaking',
        action='store_true',
        help='Disable matchmaking in report'
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
        help='ini file that contains database paths'
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
        mutational_signatures_path: args.mutational_signatures,
        validation_handle: args.validation_handle,
        disable_matchmaking: args.disable_matchmaking
    }

    output_directory = args.output_directory if args.output_directory else os.getcwd()

    config_ini = Ini.read(args.config, extended_interpolation=False, convert_to_dictionary=False)

    db_paths = Ini.read(args.dbs, extended_interpolation=True, convert_to_dictionary=True)
    if args.preclinical_dbs:
        preclinical_db_paths = Ini.read(args.preclinical_dbs, extended_interpolation=True, convert_to_dictionary=True)
    else:
        preclinical_db_paths = None

    print(db_paths)
    main(
        patient=patient_dict,
        inputs=inputs_dict,
        output_folder=output_directory,
        config=config_ini,
        dbs=db_paths['databases'],
        dbs_preclinical=preclinical_db_paths['preclinical']
    )

    end_time = time.time()
    time_statement = "Molecular Oncology Almanac runtime: %s seconds" % round((end_time - start_time), 4)
    print(time_statement)
