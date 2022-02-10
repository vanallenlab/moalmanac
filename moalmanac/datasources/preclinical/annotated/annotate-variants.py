# RUN THIS FROM moalmanac/moalmanac

import argparse
import pandas as pd

import annotator
import datasources


from config import COLNAMES
from config import CONFIG

snv_handle = 'snv_handle'
indel_handle = 'indel_handle'
bases_covered_handle = 'bases_covered_handle'
cnv_handle = 'cnv_handle'
fusion_handle = 'fusion_handle'
germline_handle = 'germline_handle'
validation_handle = 'validation_handle'

snv_input = 'snv_input'
indel_input = 'indel_input'
seg_input = 'seg_input'
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

EVIDENCE_MAP = {
    'FDA-Approved': 5.0, 'Guideline': 4.0, 'Clinical trial': 3.0,
    'Clinical evidence': 2.0, 'Preclinical': 1.0, 'Inferential': 0.0}
INV_MAP = {float(v): k for k, v in EVIDENCE_MAP.items()}


def annotate_match_1(df):
    # Biologically Relevant, no evidence
    idx_match = df['feature'].isin(almanac_genes)
    df.loc[idx_match, 'feature_match_1'] = 1
    return df


def annotate_match_2(df):
    # Gene and Dtype match; feature_match = 2, will have evidence
    match_columns = ['feature']
    df = (df
          .merge(db
                 .loc[:, match_columns + ['evidence_map', 'merged']]
                 .sort_values('evidence_map', ascending=False)
                 .drop_duplicates(match_columns, keep='first'),
                 on=match_columns,
                 how='left'
                 )
          )
    index_match = df['merged'].eq(1)
    df.loc[index_match, 'feature_match_2'] = 1
    df.loc[index_match, 'evidence'] = df.loc[index_match, 'evidence_map']
    df.drop(['merged', 'evidence_map'], axis=1, inplace=True)
    return df


def annotate_match_3(df):
    # Gene, Dtype, AType match; feature_match = 3, will have evidence
    match_columns = ['feature', 'alteration_type']
    df = (df
          .merge(db[~db['alteration_type'].eq('')]
                 .loc[:, match_columns + ['evidence_map', 'merged']]
                 .sort_values('evidence_map', ascending=False)
                 .drop_duplicates(match_columns, keep='first'),
                 on=match_columns,
                 how='left'
                 )
          )
    index_match = df['merged'].eq(1)
    df.loc[index_match, 'feature_match_3'] = 1
    df.loc[index_match, 'evidence'] = df.loc[index_match, 'evidence_map']
    df.drop(['merged', 'evidence_map'], axis=1, inplace=True)
    return df


def annotate_match_4(df):
    # Gene, Dtype, AType match, and alteration match; feature_match = 4, will have evidence
    match_columns = ['feature', 'alteration_type', 'alteration']
    df = (df
          .merge(db[~db['alteration_type'].eq('') & ~db['alteration'].eq('')]
                 .loc[:, match_columns + ['evidence_map', 'merged']]
                 .sort_values('evidence_map', ascending=False)
                 .drop_duplicates(match_columns, keep='first'),
                 on=match_columns,
                 how='left'
                 )
          )
    index_match = df['merged'].eq(1)
    df.loc[index_match, 'feature_match_4'] = 1
    df.loc[index_match, 'evidence'] = df.loc[index_match, 'evidence_map']
    df.drop(['merged', 'evidence_map'], axis=1, inplace=True)
    return df


def annotate_other_datasources(df, dbs):
    df = annotator.CancerHotspots.annotate(df, dbs)
    df = annotator.CancerHotspots3D.annotate(df, dbs)
    df = annotator.CancerGeneCensus.annotate(df, dbs)
    df = annotator.Cosmic.annotate(df, dbs)
    df = annotator.GSEACancerPathways.annotate(df, dbs)
    df = annotator.GSEACancerModules.annotate(df, dbs)
    return df


def load_db(almanac_json, table_label):
    db = pd.DataFrame(almanac_json.table(table_label).all())

    db_columns = ['feature_display', 'gene', 'variant_annotation', 'protein_change', 'predictive_implication']
    db['variant_annotation'].replace({'Oncogenic Mutations': '', 'Activating mutation': ''}, inplace=True)
    db = db.loc[:, db_columns].drop_duplicates()
    db.rename(columns={'gene': 'feature',
                       'variant_annotation': 'alteration_type',
                       'protein_change': 'alteration',
                       'predictive_implication': 'evidence'},
              inplace=True)
    db['evidence_map'] = db['evidence'].replace(EVIDENCE_MAP)
    db.sort_values(['evidence_map', 'feature_display'], ascending=[False, True], inplace=True)
    db['merged'] = 1
    return db


def load_df(handle, feature_type):
    df = pd.read_csv(handle, sep='\t', usecols=['model_id', 'feature', 'alteration_type', 'alteration'])
    df['feature_type'] = feature_type
    df['feature_match_1'] = 0
    df['feature_match_2'] = 0
    df['feature_match_3'] = 0
    df['feature_match_4'] = 0
    df['evidence'] = pd.NA
    return df


def write_df(df, handle, columns):
    df.loc[:, columns].to_csv(handle, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Annotate somatic variants with MOAlmanac')
    parser.add_argument('--directory', type=str, help='output directory')
    parser.add_argument('--input', '-i', type=str, help='input somatic variants')
    parser.add_argument('--output', '-o', type=str, help='output file name')
    args = parser.parse_args()
    root = args.directory

    feature_type = 'Somatic Variant'
    patient = {patient_id: 'TMP', tumor_type: 'NA', stage: 'NA', description: 'NA', purity: 'NA', ploidy: 'NA',
               ms_status: 'mss', wgd: False, code: 'NA', ontology: 'NA'}

    dbs = datasources.Datasources.generate_db_dict(CONFIG)
    almanac_json = datasources.Almanac.import_ds(dbs)
    almanac_genes = almanac_json.table('genes').all()[0]['genes']
    db = load_db(almanac_json, feature_type)
    dataframe = load_df(args.input, feature_type)
    dataframe = annotate_match_1(dataframe)
    dataframe = annotate_match_2(dataframe)
    dataframe = annotate_match_3(dataframe)
    dataframe = annotate_match_4(dataframe)

    dataframe['evidence_map'] = dataframe['evidence'].copy(deep=True)
    dataframe['evidence'] = dataframe['evidence'].fillna('').replace(INV_MAP)
    dataframe = annotate_other_datasources(dataframe, dbs)

    idx_almanac = dataframe[dataframe['feature_match_1'].gt(0)].index
    idx_hotspot = dataframe[dataframe['cancerhotspots_bin'].gt(0)].index
    idx_cgc = dataframe[dataframe['cgc_bin'].gt(0)].index
    #idx_cosmic = dataframe[dataframe['cosmic_bin'].gt(0)].index
    idx = idx_almanac.union(idx_hotspot).union(idx_cgc)#.union(idx_cosmic)
    dataframe = dataframe.loc[idx, :]

    use_columns = ['feature', 'alteration_type', 'alteration', 'model_id',
                   'feature_match_1', 'feature_match_2', 'feature_match_3', 'feature_match_4',
                   'evidence', 'cancerhotspots_bin', 'cancerhotspots3D_bin', 'cgc_bin',
                   'cosmic_bin', 'gsea_pathways_bin', 'gsea_modules_bin']

    output = f"{args.directory}/{args.output}"
    write_df(dataframe, output, use_columns)
