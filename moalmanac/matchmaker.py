import numpy as np
import pandas as pd
import snf
from sklearn import metrics

from annotator import Almanac as AnnotatorAlmanac
from annotator import PreclinicalMatchmaking as AnnotatorPreclinicalMatchmaking
from datasources import Almanac as DatasourceAlmanac
from datasources import CancerGeneCensus as DatasourceCGC
from datasources import Preclinical as DatasourcePreclinical

from config import COLNAMES


class Matchmaker:
    section = 'matchmaking'
    feature_type = COLNAMES[section]['feature_type']
    feature = COLNAMES[section]['feature']
    alt_type = COLNAMES[section]['alteration_type']
    alt = COLNAMES[section]['alteration']
    partner = COLNAMES[section]['partner']
    model_id = COLNAMES[section]['model_id']
    case_profile = 'case-profile'
    summary = 'summary'

    feature_columns = [feature_type, feature, alt_type, alt]

    cgc_bin = 'cgc_bin'
    fusion = 'Fusion'

    @classmethod
    def concat_case_comparisons(cls, somatic, dbs, variant_string, copy_number_string, rearrangement_string):
        somatic[cls.model_id] = cls.case_profile
        case_variants = cls.subset_dataframe_eq(somatic, cls.feature_type, variant_string)
        case_cns = cls.subset_dataframe_eq(somatic, cls.feature_type, copy_number_string)
        case_fusions = cls.subset_dataframe_eq(somatic, cls.feature_type, rearrangement_string)
        if case_fusions.shape[0] > 0:
            case_fusions = cls.format_fusions(case_fusions)

        comparison_variants = dbs[DatasourcePreclinical.variants]
        comparison_cnas = dbs[DatasourcePreclinical.cnas]
        comparison_fusions = dbs[DatasourcePreclinical.fusions]

        columns = [cls.feature, cls.alt_type, cls.alt, cls.model_id]
        variants = cls.concat_dataframes(case_variants, comparison_variants, columns)
        copy_numbers = cls.concat_dataframes(case_cns, comparison_cnas, columns)
        fusion_columns = [cls.feature, cls.partner, cls.model_id]
        fusions = cls.concat_dataframes(case_fusions, comparison_fusions, fusion_columns)

        variants[cls.feature_type] = variant_string
        copy_numbers[cls.feature_type] = copy_number_string
        fusions[cls.feature_type] = rearrangement_string
        fusions[cls.alt_type] = cls.fusion

        return {
            variant_string: variants,
            copy_number_string: copy_numbers,
            rearrangement_string: fusions
        }

    @classmethod
    def compare(cls, dbs, dbs_preclinical, somatic, case_sample_id, config):
        cgc = DatasourceCGC.import_ds(dbs)
        almanac = DatasourceAlmanac.import_ds(dbs)

        somatic_variant_biomarker_type_string = config['feature_types']['mut']
        copy_number_variant_biomarker_type_string = config['feature_types']['cna']
        fusion_biomarker_type_string = config['feature_types']['fusion']

        merged = cls.concat_case_comparisons(
            somatic = somatic,
            dbs = dbs_preclinical,
            variant_string = somatic_variant_biomarker_type_string,
            copy_number_string = copy_number_variant_biomarker_type_string,
            rearrangement_string = fusion_biomarker_type_string
        )
        annotated = AnnotatorPreclinicalMatchmaking.annotate(merged, dbs, config)
        samples_to_use = cls.subset_samples(dbs_preclinical)

        calculated = SNFTypesCGCwithEvidence.calculate(annotated, samples_to_use, cgc, almanac)
        case = cls.subset_dataframe_eq(calculated, 'case', cls.case_profile)
        case = case.reset_index(drop=True)
        case = case.sort_values(by=SNFTypesCGCwithEvidence.label, ascending=True)
        case['case'] = case['case'].replace(cls.case_profile, case_sample_id)
        return case

    @classmethod
    def concat_dataframes(cls, df1, df2, columns):
        df1 = cls.preallocate_missing_columns(df1, columns)
        df2 = cls.preallocate_missing_columns(df2, columns)
        return pd.concat([df1.loc[:, columns], df2.loc[:, columns]], ignore_index=True)

    @classmethod
    def create_empty_output(cls):
        return pd.DataFrame(columns=['case', 'comparison'])

    @classmethod
    def format_fusions(cls, dataframe):
        df = dataframe[cls.alt].fillna('').str.split('--', expand=True).rename(columns={0: 'feature', 1: 'partner'})
        df[cls.model_id] = cls.case_profile
        return df

    @classmethod
    def preallocate_missing_columns(cls, df, columns):
        missing_cols = [col for col in columns if col not in df.columns]
        return pd.concat([df, pd.DataFrame(None, columns=missing_cols, index=df.index)], axis=1)

    @classmethod
    def subset_dataframe_eq(cls, dataframe, column, value):
        return dataframe[dataframe[column].eq(value)]

    @classmethod
    def subset_samples(cls, dbs):
        summary = dbs[cls.summary]
        # this is required for python 3.12 and pandas 2.2.2 to opt into future behavior for type downcasting
        with pd.option_context("future.no_silent_downcasting", True):
            idx_samples_to_use = (
                summary['use']
                .astype(str)
                .replace('True', 1)
                .replace('False', 0)
                .astype(int)
                .eq(1)
            )
        return summary[idx_samples_to_use]['broad'].tolist() + [cls.case_profile]


class Models:
    variants = 'variants'
    cnas = 'copy_number_alterations'
    fusions = 'fusions'
    fusions_gene1 = 'fusions_gene1'
    fusions_gene2 = 'fusions_gene2'
    summary = 'summary'

    feature = 'feature'
    feature_type = 'feature_type'
    model_id = 'model_id'

    case = 'case'
    comparison = 'comparison'

    variant = 'Somatic Variant'
    copy_number = 'Copy Number'
    rearrangement = 'Rearrangement'

    jaccard = metrics.DistanceMetric.get_metric('jaccard')

    @staticmethod
    def aggregate_columns_by_row(dataframe, delimiter):
        return dataframe.agg(delimiter.join, axis=1)

    @staticmethod
    def calculate_distance(dataframe, metric):
        distance = metric(dataframe.values)
        return pd.DataFrame(distance, columns=dataframe.index.tolist(), index=dataframe.index.tolist())

    @staticmethod
    def reset_multi_indexed_dataframe(dataframe, remapped_labels):
        return dataframe.reset_index().rename(columns=remapped_labels)

    @classmethod
    def series_concat(cls, list_of_series, drop_duplicates=True):
        series = pd.concat(list_of_series)
        return cls.series_to_list(series, drop_duplicates)

    @staticmethod
    def series_to_list(series, drop_duplicates=True):
        if drop_duplicates:
            return series.drop_duplicates().dropna().sort_values().tolist()
        else:
            return series.dropna().sort_values().tolist()

    @classmethod
    def stack_distances(cls, dataframe, label):
        # Commented out sorting and dropping of duplicates to leave all N x N comparisons
        stacked = dataframe.stack().reset_index()
        stacked.rename(columns={'level_0': cls.case, 'level_1': cls.comparison, 0: label}, inplace=True)
        stacked.set_index([cls.case, cls.comparison], inplace=True)
        return stacked.loc[stacked.index, label]


class Almanac(Models):
    section = 'matchmaking'
    almanac = COLNAMES[section]['almanac']
    genes = COLNAMES[section]['genes']
    gene = COLNAMES[section]['gene']
    gene1 = COLNAMES[section]['gene1']
    gene2 = COLNAMES[section]['gene2']
    variant_annotation = COLNAMES[section]['variant_annotation']
    protein_change = COLNAMES[section]['protein_change']
    direction = COLNAMES[section]['direction']

    alteration_type = COLNAMES[section]['alteration_type']
    alteration = COLNAMES[section]['alteration']
    partner = COLNAMES[section]['partner']

    missense = 'Missense'
    truncating_types = ['Nonsense', 'Nonstop', 'Frameshift', 'Splice Site']
    truncating = 'Truncating'
    fusion = 'Fusion'

    variant_missense = 'variant_missense'
    variant_truncating = 'variant_truncating'

    feature_match = COLNAMES[section]['feature_match']
    feature_match_1 = COLNAMES[section]['feature_match_1']
    feature_match_2 = COLNAMES[section]['feature_match_2']
    feature_match_3 = COLNAMES[section]['feature_match_3']
    feature_match_4 = COLNAMES[section]['feature_match_4']

    string = 'feature_string'

    @classmethod
    def import_dbs(cls, inputs):
        db = DatasourceAlmanac.import_ds(inputs)
        tables = {}
        for table in [cls.variant, cls.copy_number, cls.rearrangement]:
            tables[table] = AnnotatorAlmanac.subset_records(db['content'], AnnotatorAlmanac.feature_type, table)
        tables[cls.genes] = db['genes']
        return tables

    @classmethod
    def generate_features(cls, dbs):
        missense = cls.generate_features_missense(dbs[cls.variant], cls.gene, cls.variant_annotation, cls.protein_change)
        truncating = cls.generate_features_truncating_aggregated(dbs[cls.variant], cls.gene, cls.variant_annotation)
        nonspecific_variant = cls.generate_features_nonspecific_variant(dbs[cls.variant])
        copy_number = cls.generate_features_copy_number(dbs[cls.copy_number], cls.gene, cls.direction)
        fusions = cls.generate_features_fusions(dbs[cls.rearrangement])
        fusion_partners = cls.generate_features_fusions_partners(dbs[cls.rearrangement])
        series = (pd
                  .Index(missense)
                  .union(pd.Index(truncating))
                  .union(pd.Index(nonspecific_variant))
                  .union(pd.Index(copy_number))
                  .union(pd.Index(fusions))
                  .union(pd.Index(fusion_partners)))
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_copy_number(cls, db, column_gene, column_type):
        series = cls.aggregate_columns_by_row(db.loc[:, [column_gene, column_type]], ' ')
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_nonspecific_variant(cls, db):
        idx = (db[cls.variant_annotation]
               .replace({'Oncogenic Mutations': '', 'Activating mutation': ''})
               .fillna('')
               .eq('')
               )
        series = db.loc[idx, cls.gene] + f' {cls.variant}'
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_missense(cls, db, column_gene, column_type, column_protein):
        idx = db[column_type].eq(cls.missense) & ~db[column_protein].fillna('').eq('')
        series = cls.aggregate_columns_by_row(db.loc[idx, [column_gene, column_protein]], ' ')
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_missense_aggregated(cls, db, column_gene, column_type, column_protein):
        idx = db[column_type].eq(cls.missense) & ~db[column_protein].fillna('').eq('')
        series = db.loc[idx, column_gene] + f' {cls.missense}'
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_truncating(cls, db, column_gene, column_type):
        idx = db[column_type].isin(cls.truncating_types)
        series = cls.aggregate_columns_by_row(db.loc[idx, [column_gene, column_type]], ' ')
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_truncating_aggregated(cls, db, column_gene, column_type):
        idx = db[column_type].isin(cls.truncating_types)
        series = db.loc[idx, column_gene] + f' {cls.truncating}'
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_fusions(cls, db):
        columns = [cls.gene1, cls.gene2]
        paired = cls.sort_fusions(db.loc[:, columns].dropna(), cls.gene1, cls.gene2)
        series = cls.aggregate_columns_by_row(paired.loc[:, columns], '--')
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_fusions_partners(cls, db):
        series1 = db.loc[db[cls.gene2].fillna('').eq(''), cls.gene1] + f' {cls.fusion}'
        series2 = db.loc[db[cls.gene2].fillna('').eq(''), cls.gene2] + f' {cls.fusion}'
        return cls.series_to_list(pd.concat([series1, series2]), drop_duplicates=True)

    @classmethod
    def generate_gene_features(cls, db):
        return cls.series_to_list(db[cls.genes])

    @classmethod
    def generate_gene_features_dtype(cls, dbs):
        variants = dbs[cls.variant][cls.gene].drop_duplicates() + f' {cls.variant}'
        copy_number = dbs[cls.copy_number][cls.gene].drop_duplicates() + f' {cls.copy_number}'
        rearrangements = (pd
                          .concat([dbs[cls.rearrangement][cls.gene1], dbs[cls.rearrangement][cls.gene2]])
                          .drop_duplicates()
                          .dropna()
                          .sort_values()
                          ) + f' {cls.rearrangement}'
        return cls.series_concat([variants, copy_number, rearrangements])

    @classmethod
    def sort_fusions(cls, df, gene1_col, gene2_col):
        for index in df.index:
            gene1 = df.loc[index, gene1_col]
            gene2 = df.loc[index, gene2_col]
            sorted_genes = sorted([gene1, gene2])
            df.loc[index, gene1_col] = sorted_genes[0]
            df.loc[index, gene2_col] = sorted_genes[1]
        return df


class AlmanacEvidence(Almanac):
    label = 'FDA features'
    description = 'We sort by agreement based measure (jaccard) by considering somatic variant, copy number, and ' \
                  'rearrangement molecular features catalogued in the Molecular Oncology Almanac that are associated ' \
                  'either with a FDA approved therapy'

    predictive_implication = 'predictive_implication'

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        almanac = cls.import_dbs(input_dtypes)
        almanac_subset = cls.subset_almanac_by_evidence(almanac, 'FDA-Approved')
        boolean_dataframe = AlmanacFeatures.create_boolean_table(input_dtypes, samples_list, almanac_subset)
        distance_dataframe = AlmanacFeatures.calculate_distance(boolean_dataframe, cls.jaccard.pairwise)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)
        return stacked_dataframe

    @classmethod
    def subset_almanac_by_evidence(cls, db, evidence_tier):
        tables = {}
        for table in [cls.variant, cls.copy_number, cls.rearrangement]:
            db_table = AnnotatorAlmanac.subset_records(db, cls.feature_type, table)
            db_table = AnnotatorAlmanac.subset_records(db_table, cls.predictive_implication, evidence_tier)
            tables[table] = pd.DataFrame(db_table)
        return tables


class AlmanacFeatures(Almanac):
    label = 'almanac_features'
    description = 'We sort by agreement based measure (jaccard) by considering all somatic variant, copy number, and ' \
                  'rearrangement molecular features catalogued in the Molecular Oncology Almanac.'

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        almanac = cls.import_dbs(input_dtypes)
        boolean_dataframe = cls.create_boolean_table(input_dtypes, samples_list, almanac)
        distance_dataframe = cls.calculate_distance(boolean_dataframe, cls.jaccard.pairwise)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)
        return stacked_dataframe

    @classmethod
    def create_bool(cls, dataframe, dtype_column):
        df = dataframe.loc[:, [cls.string, cls.model_id]].drop_duplicates()
        df[dtype_column] = 1
        return df.set_index([cls.string, cls.model_id])

    @classmethod
    def create_boolean_table(cls, inputs, samples, almanac):
        features = cls.generate_features(almanac)
        index = pd.MultiIndex.from_product([features, samples])
        columns = [cls.variant_missense, cls.variant_truncating, 'variant_nonspecific',
                   cls.copy_number,
                   cls.fusions, cls.fusions_gene1, cls.fusions_gene2
                   ]
        df = pd.DataFrame(index=index, columns=columns)

        variants = inputs[cls.variants]
        missense = cls.subset_missense(variants)
        truncating = cls.subset_truncating(variants)
        nonspecific = cls.subset_nonspecific_variant(variants, almanac)
        copy_number_alterations = cls.subset_copy_number(inputs[cls.cnas])
        fusions = cls.subset_fusions(inputs[cls.fusions])
        fusions_gene1 = cls.subset_fusion_members(inputs[cls.fusions_gene1])
        fusions_gene2 = cls.subset_fusion_members(inputs[cls.fusions_gene2])

        missense[cls.string] = cls.aggregate_columns_by_row(missense.loc[:, [cls.feature, cls.alteration]], ' ')
        # missense[cls.feature] + ' Missense'
        truncating[cls.string] = truncating[cls.feature] + f' {cls.truncating}'
        # cls.aggregate_columns_by_row(truncating.loc[:, [cls.feature, cls.alteration_type]], ' ')
        nonspecific[cls.string] = nonspecific[cls.feature] + f' {cls.variant}'
        copy_number_alterations[cls.string] = cls.aggregate_columns_by_row(
            copy_number_alterations.loc[:, [cls.feature, cls.alteration_type]], ' ')
        fusions[cls.string] = cls.aggregate_columns_by_row(fusions.loc[:, [cls.feature, cls.partner]], '--')
        fusions_gene1[cls.string] = fusions_gene1[cls.feature] + f' {cls.fusion}'
        fusions_gene2[cls.string] = fusions_gene2[cls.feature] + f' {cls.fusion}'

        df[cls.variant_missense] = cls.create_bool(missense, cls.variant_missense)
        df[cls.variant_truncating] = cls.create_bool(truncating, cls.variant_truncating)
        df['variant_nonspecific'] = cls.create_bool(nonspecific, 'variant_nonspecific')
        df[cls.copy_number] = cls.create_bool(copy_number_alterations, cls.copy_number)
        df[cls.fusions] = cls.create_bool(fusions, cls.fusions)
        df[cls.fusions_gene1] = cls.create_bool(fusions_gene1, cls.fusions_gene1)
        df[cls.fusions_gene2] = cls.create_bool(fusions_gene2, cls.fusions_gene2)

        df = cls.reset_multi_indexed_dataframe(df.fillna(0), {'level_0': cls.string, 'level_1': cls.model_id})
        df[cls.label] = df.loc[:, columns].sum(axis=1).astype(bool).astype(int)
        return df.pivot_table(index=cls.model_id, columns=cls.string, values=cls.label)

    @classmethod
    def subset_copy_number(cls, df):
        return df[df[cls.feature_match_3].eq(1)].reset_index(drop=True)

    @classmethod
    def subset_fusions(cls, df):
        dataframe = df[df[cls.feature_match_4].eq(1)]
        return cls.sort_fusions(dataframe.reset_index(drop=True), cls.feature, cls.partner)

    @classmethod
    def subset_fusion_members(cls, df):
        return df[df[cls.feature_match_2].eq(1)].reset_index(drop=True)

    @classmethod
    def subset_nonspecific_variant(cls, df, db):
        table = db[cls.variant]
        idx1 = table['variant_annotation'].replace({'Oncogenic Mutations': '', 'Activating mutation': ''}).fillna(
            '').eq('')
        idx2 = table[cls.protein_change].fillna('').eq('')
        relevant_genes = table.loc[(idx1 & idx2), cls.gene]
        return df[df[cls.feature_match_2].eq(1) & df[cls.feature].isin(relevant_genes)].reset_index(drop=True)

    @classmethod
    def subset_missense(cls, df):
        return df[(
                df[cls.feature_match_4].eq(1)
                & df[cls.alteration_type].eq(cls.missense)
                & ~df[cls.alteration_type].fillna('').eq('')
        )].reset_index(drop=True)

    @classmethod
    def subset_truncating(cls, df):
        return df[df[cls.feature_match_3].eq(1) & df[cls.alteration_type].isin(cls.truncating_types)].reset_index(
            drop=True)


class CGC(Models):
    label = 'Jaccard cgc genes'
    input_label = 'cgc'
    description = 'We sort by agreement based measure (jaccard) by considering any variant in a ' \
                  'Cancer Gene Census gene.'

    cgc_bin = 'cgc_bin'
    gene = 'feature'  # 'Gene Symbol'

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        boolean_dataframe = cls.create_boolean_table(input_dtypes, samples_list)
        distance_dataframe = cls.calculate_distance(boolean_dataframe, cls.jaccard.pairwise)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)
        return stacked_dataframe

    @classmethod
    def create_bool(cls, dataframe, value_column, dtype_column):
        df = dataframe[dataframe[value_column].astype(int).eq(1)].loc[:, [cls.feature, cls.model_id]].drop_duplicates()
        df[dtype_column] = 1
        return df.set_index([cls.feature, cls.model_id])

    @classmethod
    def create_boolean_table(cls, inputs, samples):
        cgc = inputs[cls.input_label]
        features = cls.create_features_list(cgc)

        variants = inputs[cls.variants]
        copy_number_alterations = inputs[cls.cnas]
        fusions_gene1 = inputs[cls.fusions_gene1]
        fusions_gene2 = inputs[cls.fusions_gene2]

        index = pd.MultiIndex.from_product([features, samples])
        df = pd.DataFrame(columns=[cls.variants, cls.cnas, cls.fusions_gene1, cls.fusions_gene2], index=index)
        df[cls.variants] = cls.create_bool(variants, cls.cgc_bin, cls.variants)
        df[cls.cnas] = cls.create_bool(copy_number_alterations, cls.cgc_bin, cls.cnas)
        df[cls.fusions_gene1] = cls.create_bool(fusions_gene1, cls.cgc_bin, cls.fusions_gene1)
        df[cls.fusions_gene2] = cls.create_bool(fusions_gene2, cls.cgc_bin, cls.fusions_gene2)
        df = cls.reset_multi_indexed_dataframe(df.fillna(0), {'level_0': cls.feature, 'level_1': cls.model_id})

        sum_columns = [cls.variants, cls.cnas, cls.fusions_gene1, cls.fusions_gene2]
        df[cls.label] = df.loc[:, sum_columns].sum(axis=1).astype(bool).astype(int)
        return df.pivot_table(index=cls.model_id, columns=cls.feature, values=cls.label)

    @classmethod
    def create_boolean_table_single_feature_type(cls, inputs, samples, feature_type, cgc):
        features = cls.create_features_list(cgc)

        if feature_type == cls.fusions:
            columns = [cls.fusions_gene1, cls.fusions_gene2]
        else:
            columns = [feature_type]

        index = pd.MultiIndex.from_product([features, samples])
        df = pd.DataFrame(columns=columns, index=index)
        for column in columns:
            df[column] = CGC.create_bool(inputs[column], CGC.cgc_bin, column)
        df = cls.reset_multi_indexed_dataframe(df.fillna(0), {'level_0': cls.feature, 'level_1': cls.model_id})

        sum_columns = columns
        df[cls.label] = df.loc[:, sum_columns].sum(axis=1).astype(bool).astype(int)
        return df.pivot_table(index=cls.model_id, columns=cls.feature, values=cls.label)

    @classmethod
    def create_features_list(cls, db):
        return cls.series_to_list(db[cls.gene], drop_duplicates=True)


class Report:
    case = Models.case
    comparison = Models.comparison
    case_profile = Matchmaker.case_profile

    @classmethod
    def create_report_dictionary(cls, results, reference, display_count=5):
        dictionary = {}
        if not results.empty:
            results = results[~results[cls.comparison].eq(cls.case_profile)]
            for index in results.index.tolist()[:display_count]:
                comparison = results.loc[index, cls.comparison]
                dictionary[index] = reference[comparison]
        return dictionary


class SNFTypesCGCwithEvidence(Models):
    label = 'SNF: FDA & CGC'
    description = 'We perform similarity network fusion (SNF, Bo Wang in the Goldenberg lab ' \
                  '(https://www.nature.com/articles/nmeth.2810)) using the python implementation by Ross Markello ' \
                  '(https://github.com/rmarkello/snfpy) to fuse networks that contain: ' \
                  '(1) CGC genes that contain a somatic variant, ' \
                  '(2) CGC genes that contain a copy number alteration, ' \
                  '(3) CGC genes that contain a rearrangement, ' \
                  '(4) Almanac features associated with FDA evidence.'

    @classmethod
    def calculate(cls, inputs, samples, cgc, almanac, seed=42):
        boolean_dataframe_variants = CGC.create_boolean_table_single_feature_type(inputs, samples, cls.variants, cgc)
        boolean_dataframe_copy_numbers = CGC.create_boolean_table_single_feature_type(inputs, samples, cls.cnas, cgc)
        boolean_dataframe_fusions = CGC.create_boolean_table_single_feature_type(inputs, samples, cls.fusions, cgc)

        almanac_content = almanac['content']
        almanac_subset = AlmanacEvidence.subset_almanac_by_evidence(almanac_content, 'FDA-Approved')
        boolean_dataframe_1 = AlmanacFeatures.create_boolean_table(inputs, samples, almanac_subset)

        data = [
            boolean_dataframe_variants.loc[samples, :].fillna(0),
            boolean_dataframe_copy_numbers.loc[samples, :].fillna(0),
            boolean_dataframe_fusions.loc[samples, :].fillna(0),
            boolean_dataframe_1.loc[samples, :].fillna(0),
        ]

        np.random.seed(seed)
        affinity_networks = snf.make_affinity(data, metric='jaccard', normalize=False, K=20, mu=0.5)
        fused_network = snf.snf(affinity_networks, K=20)
        fused_dataframe = pd.DataFrame(fused_network, index=samples, columns=samples)
        distance_dataframe = pd.DataFrame(1, index=samples, columns=samples).subtract(fused_dataframe)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)
        return stacked_dataframe.reset_index()
