# RUN THIS FROM moalmanac/moalmanac, instead of moalmanac/datasources/gdsc
# Comment line 687 to not add addito

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


def annotate_data(df, db, genes, consider_partner=False):
    # Biologically Relevant, no evidence
    idx_match = df['feature'].isin(genes)
    df.loc[idx_match, 'feature_match_1'] = 1

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

    if consider_partner:
        # Gene, Dtype, AType match, and alteration match; feature_match = 4, will have evidence
        match_columns = ['feature', 'alteration_type', 'partner']
        df = (df
              .merge(db[~db['alteration_type'].eq('') & ~db['partner'].eq('')]
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

    feature_match_columns = ['feature_match_1', 'feature_match_2', 'feature_match_3', 'feature_match_4']
    df['feature_match'] = df.loc[:, feature_match_columns].sum(axis=1)

    inv_map = {float(v): k for k, v in EVIDENCE_MAP.items()}
    df['evidence_map'] = df['evidence'].copy(deep=True)
    df['evidence'] = df['evidence'].fillna('').replace(inv_map)
    return df


def annotate_both_genes(df_gene1, df_gene2, db_gene1, db_gene2, almanac_genes, df, dbs):
    # Gene1 and Gene 2
    group1 = annotate_data(df_gene1, db_gene1, almanac_genes, True)
    group2 = annotate_data(df_gene1, db_gene2, almanac_genes, True)
    group3 = annotate_data(df_gene2, db_gene1, almanac_genes, True)
    group4 = annotate_data(df_gene2, db_gene2, almanac_genes, True)

    values = pd.concat([
        group1.rename(columns={'evidence_map': 'group1'})['group1'],
        group2.rename(columns={'evidence_map': 'group2'})['group2'],
        group3.rename(columns={'evidence_map': 'group3'})['group3'],
        group4.rename(columns={'evidence_map': 'group4'})['group4'],
    ], axis=1)
    values = values.fillna(-1).idxmax(axis=1)

    idx_group1 = values[values.eq('group1')].index
    idx_group2 = values[values.eq('group2')].index
    idx_group3 = values[values.eq('group3')].index
    idx_group4 = values[values.eq('group4')].index

    columns = ['evidence', 'feature_match_1', 'feature_match_2', 'feature_match_3', 'feature_match_4', 'feature_match']
    pairs = [(group1, idx_group1), (group2, idx_group2), (group3, idx_group3), (group4, idx_group4)]

    df['which_match'] = ''
    for index in [idx_group1, idx_group2]:
        df.loc[index, 'which_match'] = df.loc[index, 'feature']
    for index in [idx_group3, idx_group4]:
        df.loc[index, 'which_match'] = df.loc[index, 'partner']

    for group, index in pairs:
        df.loc[index, columns] = group.loc[index, columns]

    columns = ['feature', 'partner']
    for index in df.index:
        df.loc[index, columns] = sorted(df.loc[index, columns].astype(str).tolist())

    df = (df
          .sort_values(['evidence', 'feature_match'], ascending=False)
          .drop_duplicates(['feature', 'partner', 'model_id'], keep='first')
          )

    df['alteration'] = ''
    df = annotate_other_datasources(df, dbs)
    return df


def annotate_gene_1(df_gene1, db_gene1, db_gene2, almanac_genes, dbs):
    # Gene 1
    group1 = annotate_data(df_gene1, db_gene1, almanac_genes, False)
    group2 = annotate_data(df_gene1, db_gene2, almanac_genes, False)

    values = pd.concat([
        group1.rename(columns={'evidence_map': 'group1'})['group1'],
        group2.rename(columns={'evidence_map': 'group2'})['group2'],
    ], axis=1)
    values = values.fillna(-1).idxmax(axis=1)

    idx_group1 = values[values.eq('group1')].index
    idx_group2 = values[values.eq('group2')].index

    columns = ['evidence', 'feature_match']
    group1.loc[idx_group2, columns] = group2.loc[idx_group2, columns]

    group1 = (group1
              .sort_values(['evidence', 'feature_match'], ascending=False)
              .drop_duplicates(['feature', 'model_id'], keep='first')
              )

    group1['alteration'] = ''
    group1 = annotate_other_datasources(group1, dbs)
    return group1


def annotate_gene_2(df_gene2, db_gene1, db_gene2, almanac_genes, dbs):
    # Gene 2
    group3 = annotate_data(df_gene2, db_gene1, almanac_genes, False)
    group4 = annotate_data(df_gene2, db_gene2, almanac_genes, False)

    values = pd.concat([
        group3.rename(columns={'evidence_map': 'group3'})['group3'],
        group4.rename(columns={'evidence_map': 'group4'})['group4'],
    ], axis=1)
    values = values.fillna(-1).idxmax(axis=1)

    idx_group3 = values[values.eq('group3')].index
    idx_group4 = values[values.eq('group4')].index

    columns = ['evidence', 'feature_match']
    group3.loc[idx_group4, columns] = group4.loc[idx_group4, columns]

    group3 = (group3
              .sort_values(['evidence', 'feature_match'], ascending=False)
              .drop_duplicates(['feature', 'model_id'], keep='first')
              )

    group3['alteration'] = ''
    group3 = annotate_other_datasources(group3, dbs)
    return group3


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

    db_columns = ['feature_display', 'gene1', 'gene2', 'rearrangement_type', 'predictive_implication']
    db = db.loc[:, db_columns].drop_duplicates()
    db.rename(columns={
        'rearrangement_type': 'alteration_type',
        'predictive_implication': 'evidence'},
        inplace=True)
    db['evidence_map'] = db['evidence'].replace(EVIDENCE_MAP)
    db.sort_values(['evidence_map', 'feature_display'], ascending=[False, True], inplace=True)
    db['merged'] = 1
    return db


def load_df(handle, feature_type):
    df = pd.read_csv(handle, sep='\t', usecols=['model_id', 'feature', 'partner'])
    df['feature_type'] = feature_type
    df['alteration_type'] = 'Fusion'
    df['feature_match'] = 0
    df['feature_match_1'] = 0
    df['feature_match_2'] = 0
    df['feature_match_3'] = 0
    df['feature_match_4'] = 0
    df['evidence'] = pd.NA
    return df


def rreplace(string, old_value, new_value, occurrence):
    revised = string.rsplit(old_value, occurrence)
    return new_value.join(revised)


def write_df(df, handle, columns):
    df.loc[:, columns].to_csv(handle, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Annotate somatic variants with MOAlmanac')
    parser.add_argument('--directory', type=str, help='output directory')
    parser.add_argument('--input', '-i', type=str, help='input fusions')
    parser.add_argument('--output', '-o', type=str, help='output file name')
    args = parser.parse_args()
    root = args.directory

    feature_type = "Rearrangement"
    patient = {patient_id: 'TMP', tumor_type: 'NA', stage: 'NA', description: 'NA', purity: 'NA', ploidy: 'NA',
               ms_status: 'mss', wgd: False, code: 'NA', ontology: 'NA'}

    dbs = datasources.Datasources.generate_db_dict(CONFIG)
    almanac_json = datasources.Almanac.import_ds(dbs)
    almanac_genes = almanac_json.table('genes').all()[0]['genes']
    db = load_db(almanac_json, feature_type)
    dataframe = load_df(args.input, feature_type)

    df_gene1 = dataframe.copy(deep=True)
    df_gene2 = dataframe.copy(deep=True)

    df_gene1['alteration'] = dataframe['feature'] + '--' + dataframe['partner']
    df_gene2['alteration'] = dataframe['partner'] + '--' + dataframe['feature']
    df_gene2.rename(columns={'partner': 'feature', 'feature': 'partner'}, inplace=True)

    db_gene1 = db.copy(deep=True).rename(columns={'gene1': 'feature', 'gene2': 'partner'}).fillna('')
    db_gene2 = db.copy(deep=True).rename(columns={'gene1': 'partner', 'gene2': 'feature'}).fillna('')

    both_genes = annotate_both_genes(df_gene1, df_gene2, db_gene1, db_gene2, almanac_genes, dataframe, dbs)
    gene1 = annotate_gene_1(df_gene1, db_gene1, db_gene2, almanac_genes, dbs)
    gene2 = annotate_gene_2(df_gene2, db_gene1, db_gene2, almanac_genes, dbs)

    for dataframe in [both_genes, gene1, gene2]:
        idx_almanac = dataframe[dataframe['feature_match_1'].gt(0)].index
        idx_hotspot = dataframe[dataframe['cancerhotspots_bin'].gt(0)].index
        idx_cgc = dataframe[dataframe['cgc_bin'].gt(0)].index
        #idx_cosmic = dataframe[dataframe['cosmic_bin'].gt(0)].index
        idx = idx_almanac.union(idx_hotspot).union(idx_cgc)#.union(idx_cosmic)
        dataframe = dataframe.loc[idx, :]

    use_columns = ['feature', 'partner', 'which_match', 'model_id',
                   'feature_match', 'feature_match_1', 'feature_match_2', 'feature_match_3', 'feature_match_4',
                   'evidence', 'cancerhotspots_bin', 'cancerhotspots3D_bin', 'cgc_bin',
                   'cosmic_bin', 'gsea_pathways_bin', 'gsea_modules_bin']

    output = f"{args.directory}/{args.output}"
    write_df(both_genes, output, use_columns)

    use_columns.remove('which_match')
    write_df(gene1, rreplace(output, '.', '.gene1.', 1), use_columns)
    write_df(gene2, rreplace(output, '.', '.gene2.', 1), use_columns)
