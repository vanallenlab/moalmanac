import time
import argparse
import pandas as pd

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


def main(patient, inputs):
    dbs = datasources.Datasources.generate_db_dict(CONFIG)

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
        annotated_somatic = annotator.Annotator.annotate_somatic(somatic_variants, dbs, patient[code], patient[patient_id])
        evaluated_somatic = evaluator.Evaluator.evaluate_somatic(annotated_somatic)

        validation_accept, validation_reject = features.MAFValidation.import_feature(inputs[validation_handle])
        if not validation_accept.empty:
            evaluated_somatic = annotator.OverlapValidation.append_validation(evaluated_somatic, validation_accept)
            illustrator.ValidationOverlap.generate_dna_rna_plot(evaluated_somatic, patient[patient_id])
    else:
        evaluated_somatic = features.Features.create_empty_dataframe()

    if not germline_variants.empty:
        annotated_germline = annotator.Annotator.annotate_germline(germline_variants, dbs, patient[code], patient[patient_id])
        evaluated_germline = evaluator.Evaluator.evaluate_germline(annotated_germline)
    else:
        evaluated_germline = features.Features.create_empty_dataframe()

    evaluated_somatic = annotator.OverlapSomaticGermline.append_germline_hits(evaluated_somatic, evaluated_germline)
    integrated = evaluator.Integrative.evaluate(evaluated_somatic, evaluated_germline, dbs, feature_types)

    somatic_burden = features.BurdenReader.import_feature(inputs[bases_covered_handle], patient, somatic_variants, dbs)
    somatic_signatures = features.CosmicSignatures.import_feature(inputs[snv_handle], patient)
    patient_wgd = features.Aneuploidy.summarize(patient[wgd])
    patient_ms_status = features.MicrosatelliteReader.summarize(patient[ms_status])
    patient[ms_status] = features.MicrosatelliteReader.map_status(patient[ms_status])

    annotated_burden = annotator.Annotator.annotate_almanac(somatic_burden, dbs, patient[code], patient[patient_id])
    annotated_signatures = annotator.Annotator.annotate_almanac(somatic_signatures, dbs, patient[code], patient[patient_id])
    annotated_wgd = annotator.Annotator.annotate_almanac(patient_wgd, dbs, patient[code], patient[patient_id])
    annotated_ms_status = annotator.Annotator.annotate_almanac(patient_ms_status, dbs, patient[code], patient[patient_id])

    evaluated_burden = evaluator.Evaluator.evaluate_almanac(annotated_burden)
    evaluated_signatures = evaluator.Evaluator.evaluate_almanac(annotated_signatures)
    evaluated_wgd = evaluator.Evaluator.evaluate_almanac(annotated_wgd)
    evaluated_ms_variants = evaluator.Microsatellite.evaluate_variants(evaluated_somatic, evaluated_germline)
    evaluated_ms_status = evaluator.Microsatellite.evaluate_status(annotated_ms_status, evaluated_ms_variants)

    actionable = evaluator.Actionable.evaluate(evaluated_somatic, evaluated_germline,
                                               evaluated_ms_variants, evaluated_ms_status,
                                               evaluated_burden, evaluated_signatures, evaluated_wgd)

    strategies = evaluator.Strategies.report_therapy_strategies(actionable)

    dbs_preclinical = datasources.Preclinical.import_dbs()
    efficacy_dictionary = investigator.SensitivityDictionary.create(dbs_preclinical, actionable, patient[patient_id])
    efficacy_summary = investigator.SummaryDataFrame.create(efficacy_dictionary, actionable, patient[patient_id])
    actionable = annotator.PreclinicalEfficacy.annotate(actionable, efficacy_summary)

    if inputs[disable_matchmaking]:
        matchmaker_results = matchmaker.Matchmaker.create_empty_output()
    else:
        matchmaker_results = matchmaker.Matchmaker.compare(dbs, dbs_preclinical, evaluated_somatic, patient[patient_id])

    report_dictionary = reporter.Reporter.generate_dictionary(evaluated_somatic, patient)
    reporter.Reporter.generate_report(actionable,
                                      report_dictionary,
                                      efficacy_dictionary,
                                      efficacy_summary,
                                      matchmaker_results,
                                      dbs_preclinical['dictionary']
                                      )

    writer.Actionable.write(actionable, patient[patient_id])
    writer.GermlineACMG.write(evaluated_germline, patient[patient_id])
    writer.GermlineCancer.write(evaluated_germline, patient[patient_id])
    writer.GermlineHereditary.write(evaluated_germline, patient[patient_id])
    writer.Integrated.write(integrated, patient[patient_id])
    writer.MSI.write(evaluated_ms_variants, patient[patient_id])
    writer.MutationalBurden.write(evaluated_burden, patient[patient_id])
    writer.SomaticScored.write(evaluated_somatic, patient[patient_id])
    writer.SomaticFiltered.write(somatic_filtered, patient[patient_id])
    writer.Strategies.write(strategies, patient[patient_id])
    writer.PreclinicalEfficacy.write(efficacy_summary, patient[patient_id])
    writer.PreclinicalMatchmaking.write(matchmaker_results, patient[patient_id])


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

    main(patient_dict, inputs_dict)

    end_time = time.time()
    time_statement = "Molecular Oncology Almanac runtime: %s seconds" % round((end_time - start_time), 4)
    print(time_statement)
