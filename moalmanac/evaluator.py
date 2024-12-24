import numpy as np
import pandas as pd
import re

import datasources
import features
from config import COLNAMES


class Evaluator:
    """
    Evaluate based on annotated bins
    """

    almanac_bin_map = {
        1: 'Biologically Relevant',
        2: 'Investigate Actionability',  # - Low',
        3: 'Investigate Actionability',  # - High',
        4: 'Putatively Actionable'
    }

    bin_section = 'bin_names'
    score_bin = COLNAMES[bin_section]['score_bin']
    almanac_bin = COLNAMES[bin_section]['almanac']
    sensitive_bin = COLNAMES[bin_section]['sensitive_score_bin']
    resistance_bin = COLNAMES[bin_section]['resistance_score_bin']
    prognostic_bin = COLNAMES[bin_section]['prognostic_score_bin']
    acmg_bin = COLNAMES[bin_section]['acmg']
    cancerhotspots_bin = COLNAMES[bin_section]['cancerhotspots']
    cancerhotspots3d_bin = COLNAMES[bin_section]['cancerhotspots3d']
    cgc_bin = COLNAMES[bin_section]['cgc']
    clinvar_bin = COLNAMES[bin_section]['clinvar']
    cosmic_bin = COLNAMES[bin_section]['cosmic']
    gsea_modules_bin = COLNAMES[bin_section]['gsea_modules']
    gsea_pathways_bin = COLNAMES[bin_section]['gsea_pathways']
    hereditary_bin = COLNAMES[bin_section]['hereditary']
    msi_bin = COLNAMES[bin_section]['msi']
    exac_common_bin = COLNAMES[bin_section]['exac_common']

    datasources_section = 'datasources'
    sensitive_feature_type = COLNAMES[datasources_section]['sensitive_feature_type']
    resistance_feature_type = COLNAMES[datasources_section]['resistance_feature_type']
    prognostic_feature_type = COLNAMES[datasources_section]['prognostic_feature_type']
    sensitive_feature = COLNAMES[datasources_section]['sensitive_feature']
    resistance_feature = COLNAMES[datasources_section]['resistance_feature']
    prognostic_feature = COLNAMES[datasources_section]['prognostic_feature']
    sensitive_alt_type = COLNAMES[datasources_section]['sensitive_alt_type']
    resistance_alt_type = COLNAMES[datasources_section]['resistance_alt_type']
    prognostic_alt_type = COLNAMES[datasources_section]['prognostic_alt_type']
    sensitive_alt = COLNAMES[datasources_section]['sensitive_alt']
    resistance_alt = COLNAMES[datasources_section]['resistance_alt']
    prognostic_alt = COLNAMES[datasources_section]['prognostic_alt']
    sensitive_implication_map = COLNAMES[datasources_section]['sensitive_implication_map']
    resistance_implication_map = COLNAMES[datasources_section]['resistance_implication_map']
    prognostic_implication_map = COLNAMES[datasources_section]['prognostic_implication_map']

    sensitivity = COLNAMES[datasources_section]['sensitivity']
    sensitive_therapy_name = COLNAMES[datasources_section]['sensitive_therapy']
    sensitive_therapy_strategy = COLNAMES[datasources_section]['sensitive_therapy_strategy']
    resistance = COLNAMES[datasources_section]['resistance']
    resistance_therapy_name = COLNAMES[datasources_section]['resistance_therapy']
    resistance_therapy_strategy = COLNAMES[datasources_section]['resistance_therapy_strategy']

    features_section = 'features'
    feature_type = COLNAMES[features_section]['feature_type']
    feature = COLNAMES[features_section]['feature']
    alt_type = COLNAMES[features_section]['alt_type']
    alt = COLNAMES[features_section]['alt']
    tumor_f = COLNAMES[features_section]['tumor_f']
    coverage = COLNAMES[features_section]['coverage']
    feature_str = COLNAMES[features_section]['feature_str']
    feature_display = COLNAMES[features_section]['feature_display']
    left_gene = COLNAMES[features_section]['left_gene']
    right_gene = COLNAMES[features_section]['right_gene']

    burden_section = 'burden'
    high_burden_boolean = COLNAMES[burden_section]['high_burden_boolean']
    mutational_burden = COLNAMES[burden_section]['mutational_burden']

    microsatellite_section = 'microsatellite'
    supporting_variants = COLNAMES[microsatellite_section]['supporting_variants']

    @classmethod
    def assign_bin(cls, df, bin_column, bin_label):
        series_score_bin = df.loc[:, cls.score_bin]
        idx_empty_score_bin = series_score_bin.isnull()
        idx_nonzero_bin = df.loc[:, bin_column].astype(float) != 0.0
        idx = idx_empty_score_bin & idx_nonzero_bin
        df.loc[idx, cls.score_bin] = bin_label
        return df.loc[:, cls.score_bin]

    @classmethod
    def assign_vus(cls, df):
        series_score_bin = df.loc[:, cls.score_bin]
        return series_score_bin.fillna('VUS')

    @classmethod
    def evaluate_almanac(cls, df):
        df[cls.score_bin] = cls.map_almanac_bins(df[cls.almanac_bin])
        for bin_column in [cls.sensitive_bin, cls.resistance_bin, cls.prognostic_bin]:
            df['{}_map'.format(bin_column)] = df[bin_column]
            df[bin_column] = cls.map_almanac_bins(df[bin_column])
        return df

    @classmethod
    def evaluate_somatic(cls, df):
        df = cls.evaluate_almanac(df)
        df[cls.score_bin] = cls.assign_bin(df, cls.cancerhotspots_bin, 'Cancer Hotspot')
        df[cls.score_bin] = cls.assign_bin(df, cls.cancerhotspots3d_bin, 'Cancer Hotspot 3D')
        df[cls.score_bin] = cls.assign_bin(df, cls.cgc_bin, 'Cancer Gene Census')
        df[cls.score_bin] = cls.assign_bin(df, cls.gsea_pathways_bin, 'Cancer Pathway')
        df[cls.score_bin] = cls.assign_bin(df, cls.gsea_modules_bin, 'Cancer Module')
        df[cls.score_bin] = cls.assign_bin(df, cls.cosmic_bin, 'Cosmic')
        df[cls.score_bin] = cls.assign_vus(df)
        return df

    @classmethod
    def evaluate_germline(cls, df):
        df = cls.evaluate_almanac(df)
        df[cls.score_bin] = cls.assign_bin(df, cls.cancerhotspots_bin, 'Cancer Hotspot')
        df[cls.score_bin] = cls.assign_bin(df, cls.cgc_bin, 'Cancer Gene Census')
        return df

    @classmethod
    def map_almanac_bins(cls, series):
        return series.map(cls.almanac_bin_map)

    @staticmethod
    def remap_almanac_bins(series, old_value, new_value):
        return series.astype(int).replace(to_replace=old_value, value=new_value)

    @classmethod
    def remove_low_allele_fraction_variants(cls, df, minimum_allele_fraction=0.05):
        # hard coding somatic and germline variant strings until future refactor
        idx_mut = df[df[cls.feature_type].isin(['Somatic Variant', 'Germline Variant'])].index
        idx_low_quality = df[df[cls.tumor_f].astype(float).lt(float(minimum_allele_fraction))].index
        idx_low_quality_muts = idx_mut.intersection(idx_low_quality)
        idx = df.index.difference(idx_low_quality_muts)
        return df.loc[idx, :]

    @classmethod
    def remove_low_coverage_variants(cls, df, minimum_coverage=15):
        idx_mut = df[df[cls.feature_type].isin(['Somatic Variant', 'Germline Variant'])].index
        idx_low_quality = df[df[cls.coverage].astype(float).le(float(minimum_coverage))].index
        idx_low_quality_muts = idx_mut.intersection(idx_low_quality)
        idx = df.index.difference(idx_low_quality_muts)
        return df.loc[idx, :]

    @classmethod
    def remove_benign_variants(cls, df):
        return df[df[cls.clinvar_bin].astype(float) != 0.0]

    @classmethod
    def remove_common_variants(cls, df):
        return df[df[cls.exac_common_bin].astype(float) != 1.0]

    @staticmethod
    def split_camelcase(string):
        split = re.sub('([a-z])([A-Z])', r'\1 \2', string).split()
        for idx, value in enumerate(split[1:]):
            split[idx + 1] = value.lower()
        return ' '.join(split)

    @classmethod
    def subset_almanac_bin(cls, df):
        return df[df[cls.almanac_bin].astype(float).fillna(0.0) != 0.0]


class Actionable:
    sort_columns = [Evaluator.almanac_bin,
                    Evaluator.sensitive_implication_map,
                    Evaluator.resistance_implication_map,
                    Evaluator.prognostic_implication_map]

    # Map favorable prognosis to favorable and unfavorable

    @staticmethod
    def create_string_list(series):
        return ', '.join(map(str, series.unique()))

    @classmethod
    def display_aneuploidy(cls, df, idx):
        return df.loc[idx, Evaluator.feature]

    @classmethod
    def display_burden(cls, df, idx):
        return df.loc[idx, Evaluator.alt].astype(str)

    @classmethod
    def display_copynumber(cls, df, idx):
        gene = df.loc[idx, Evaluator.feature]
        direction = df.loc[idx, Evaluator.alt_type]
        # Copy Number: CDKN2A Deletion
        return gene + ' ' + direction

    @classmethod
    def display_fusion(cls, df, idx):
        fusion = df.loc[idx, Evaluator.alt]
        # Rearrangement: BCR--ABL1 Fusion
        return fusion + ' Fusion'

    @classmethod
    def display_microsatellite_stability(cls, df, idx):
        return df.loc[idx, Evaluator.feature]

    @classmethod
    def display_microsatellite_variants(cls, df, idx):
        return df.loc[idx, Evaluator.feature] + ': ' + df.loc[idx, Evaluator.alt]

    @classmethod
    def display_signature(cls, df, idx):
        # before_string = "COSMIC Signature"
        # after_string = f"COSMIC Signature (version {version})"
        signature = df.loc[idx, Evaluator.feature]#.str.replace(before_string, after_string)
        contribution = df.loc[idx, Evaluator.alt].astype(float).multiply(100).round(0).astype(int).astype(str)
        # Signature: Cosmic Signature 7 (65%)
        return signature + ' (' + contribution + '%)'

    @classmethod
    def display_variant(cls, df, idx):
        gene = df.loc[idx, Evaluator.feature]
        protein_change = df.loc[idx, Evaluator.alt]
        variant_class = df.loc[idx, Evaluator.alt_type]
        # exon, pathogenic, cDNA, linebreaks
        # Gene p.Foo (c.DNA)
        # Exon 12 Missense
        # Pathogenic
        return gene + ' ' + protein_change + ' (' + variant_class + ')'

    @classmethod
    def evaluate(cls, somatic, germline, ms_variants, ms_status, burden, signatures, wgd, config):
        somatic = cls.format_mutations(somatic)
        germline = cls.format_mutations(germline)

        germline = Evaluator.remove_benign_variants(germline)
        germline = Evaluator.remove_common_variants(germline)

        ms_variants_summary = cls.summarize_ms_variants(ms_variants, config)

        if not burden.loc[0, Evaluator.high_burden_boolean]:
            burden = burden.drop(burden.index[0])

        actionable_list = []
        for dataframe in [somatic, germline, ms_variants_summary, ms_status, burden, signatures, wgd]:
            actionable_list.append(Evaluator.subset_almanac_bin(dataframe))
        df = features.Features.concat_list_of_dataframes(list_of_dataframes=actionable_list).reset_index(drop=True)

        df[Evaluator.feature_display] = cls.format_feature_display(df=df, config=config)
        return df.sort_values(cls.sort_columns, ascending=False)

    @classmethod
    def format_feature_display(cls, df, config):
        display_column = Evaluator.feature_display
        feature_type_column = Evaluator.feature_type
        feature_column = Evaluator.feature

        biomarker_types = config['feature_types']
        idx_somatic = df[feature_type_column].isin([biomarker_types['mut']])
        idx_germline = df[feature_type_column].isin([biomarker_types['germline']])
        idx_cn = df[feature_type_column].isin([biomarker_types['cna']])
        idx_fusion = df[feature_type_column].isin([biomarker_types['fusion']])
        idx_msi = df[feature_type_column].isin([biomarker_types['microsatellite']])
        idx_msi_variants = df[feature_column].isin([Evaluator.supporting_variants])
        idx_msi = idx_msi & ~idx_msi_variants
        idx_burden = df[feature_type_column].isin([biomarker_types['burden']])
        idx_signature = df[feature_type_column].isin([biomarker_types['signature']])
        idx_wgd = df[feature_type_column].isin([biomarker_types['aneuploidy']])

        df.loc[idx_somatic, display_column] = cls.display_variant(df=df, idx=idx_somatic)
        df.loc[idx_germline, display_column] = cls.display_variant(df=df, idx=idx_germline)
        df.loc[idx_cn, display_column] = cls.display_copynumber(df=df, idx=idx_cn)
        df.loc[idx_fusion, display_column] = cls.display_fusion(df=df, idx=idx_fusion)
        df.loc[idx_msi, display_column] = cls.display_microsatellite_stability(df=df, idx=idx_msi)
        df.loc[idx_msi_variants, display_column] = cls.display_microsatellite_variants(df=df, idx=idx_msi_variants)
        df.loc[idx_burden, display_column] = cls.display_burden(df=df, idx=idx_burden)
        df.loc[idx_signature, display_column] = cls.display_signature(df=df, idx=idx_signature)
        df.loc[idx_wgd, display_column] = cls.display_aneuploidy(df=df, idx=idx_wgd)
        return df.loc[:, display_column]

    @classmethod
    def format_mutations(cls, df):
        for column in [Evaluator.alt_type, Evaluator.sensitive_alt_type, Evaluator.resistance_alt_type,
                       Evaluator.prognostic_alt_type]:
            df[column] = cls.format_variant_classification(df[Evaluator.alt_type].fillna(''))
        return df

    @classmethod
    def format_variant_classification(cls, series):
        return series.str.replace('_Mutation', '')

    @classmethod
    def summarize_ms_variants(cls, df, config):
        df = cls.format_mutations(df)
        msi_summary = features.Features.create_empty_dataframe()
        if not df.empty:
            feature = Evaluator.supporting_variants
            feature_displays = cls.format_feature_display(df=df, config=config)
            feature_displays_list = cls.create_string_list(feature_displays)

            msi_summary.loc[0, Evaluator.feature_type] = config['feature_types']['microsatellite']
            msi_summary.loc[0, Evaluator.feature] = feature
            msi_summary.loc[0, Evaluator.alt] = feature_displays_list
            msi_summary.loc[0, Evaluator.almanac_bin] = 1
            msi_summary.loc[0, Evaluator.score_bin] = Evaluator.almanac_bin_map[1]
            msi_summary.loc[0, Evaluator.feature_display] = feature + ': ' + feature_displays_list
        return msi_summary


class Integrative:
    feature = datasources.Datasources.feature
    feature_type = datasources.Datasources.feature_type
    alt_type = datasources.Datasources.alt_type
    alt = datasources.Datasources.alt

    genes = datasources.Almanac.genes

    integrative_section = 'integrative'
    somatic = COLNAMES[integrative_section]['somatic']
    copynumber = COLNAMES[integrative_section]['copynumber']
    fusion = COLNAMES[integrative_section]['fusion']
    germline = COLNAMES[integrative_section]['germline']

    columns = [somatic, copynumber, fusion, germline]

    @classmethod
    def create_integrated_df(cls, genes):
        return pd.DataFrame(None, columns=cls.columns, index=genes)

    @classmethod
    def evaluate(cls, somatic, germline, dbs, config):
        genes = cls.return_datasource_genes(dbs)
        df = cls.create_integrated_df(genes)

        somatic_mutations = cls.extract_feature_type(somatic, config['feature_types']['mut'])
        somatic_copynumbers = cls.extract_feature_type(somatic, config['feature_types']['cna'])
        somatic_fusions = cls.extract_feature_type(somatic, config['feature_types']['fusion'])

        for gene in df.index:
            gene_muts = somatic_mutations[somatic_mutations[cls.feature].astype(str) == gene]
            gene_cnas = somatic_copynumbers[somatic_copynumbers[cls.feature].astype(str) == gene]
            gene_fusion = somatic_fusions[somatic_fusions[cls.feature].astype(str) == gene]
            gene_germline = germline[germline[cls.feature].astype(str) == gene]

            df.loc[gene, cls.copynumber] = cls.join_alterations(gene_cnas, cls.alt_type)
            df.loc[gene, cls.fusion] = cls.join_alterations(gene_fusion, cls.feature)
            df.loc[gene, cls.somatic] = cls.join_alteration_types(gene_muts, [cls.alt_type, cls.alt])
            df.loc[gene, cls.germline] = cls.join_alteration_types(gene_germline, [cls.alt_type, cls.alt])

        df = cls.subset_nonempty_df(df)
        return df

    @classmethod
    def extract_feature_type(cls, df, feature_type):
        return df[df[cls.feature_type] == feature_type]

    @staticmethod
    def join_alteration_types(df, columns):
        return ", ".join(df[columns[0]].fillna("") + " " + df[columns[1]].fillna(""))

    @staticmethod
    def join_alterations(df, column):
        return ', '.join(df[column])

    @classmethod
    def return_datasource_genes(cls, dbs):
        ds_almanac = datasources.Almanac.import_ds(dbs)
        ds_hotspots = datasources.CancerHotspots.import_ds(dbs)
        ds_cgc = datasources.CancerHotspots.import_ds(dbs)

        genes_almanac = cls.return_genes_almanac(ds_almanac)
        genes_hotspots = cls.return_genes(ds_hotspots)
        genes_cgc = cls.return_genes(ds_cgc)

        return sorted(list(set(genes_almanac + genes_hotspots + genes_cgc)))

    @classmethod
    def return_genes(cls, ds):
        return ds[cls.feature].dropna().unique().tolist()

    @classmethod
    def return_genes_almanac(cls, ds):
        return ds[cls.genes]

    @classmethod
    def subset_nonempty_df(cls, df):
        return df[df.astype(bool).sum(axis=1).astype(float) != 0]


class Microsatellite(object):
    microsatellite_section = 'microsatellite'
    msi = COLNAMES[microsatellite_section]['msi']
    msih = COLNAMES[microsatellite_section]['msih']
    msil = COLNAMES[microsatellite_section]['msil']
    mss = COLNAMES[microsatellite_section]['mss']

    @classmethod
    def return_msi_variants(cls, df):
        idx_msi = df[Evaluator.msi_bin].fillna(0.0).astype(float) == float(1.0)
        idx_missense = df[Evaluator.alt_type].fillna('').str.contains('Missense')
        return df[idx_msi & ~idx_missense]

    @classmethod
    def evaluate_status(cls, df, variants):
        columns = [Evaluator.almanac_bin, Evaluator.sensitive_bin, Evaluator.resistance_bin, Evaluator.prognostic_bin]
        if variants.empty:
            for bin_column in columns:
                df[bin_column] = Evaluator.remap_almanac_bins(series=df[bin_column], old_value=3, new_value=2)
        return Evaluator.evaluate_almanac(df)

    @classmethod
    def evaluate_variants(cls, somatic, germline):
        #somatic = Evaluator.remove_low_coverage_variants(somatic)
        #somatic = Evaluator.remove_low_allele_fraction_variants(somatic)

        #germline = Evaluator.remove_low_coverage_variants(germline)
        #germline = Evaluator.remove_low_allele_fraction_variants(germline)
        germline = Evaluator.remove_benign_variants(germline)
        germline = Evaluator.remove_common_variants(germline)

        msi_somatic = cls.return_msi_variants(somatic)
        msi_germline = cls.return_msi_variants(germline)
        msi_variants = [msi_somatic, msi_germline]
        return features.Features.concat_list_of_dataframes(list_of_dataframes=msi_variants)


class Strategies:
    sensitivity = Evaluator.sensitivity
    sensitive_therapy_name = Evaluator.sensitive_therapy_name
    sensitive_therapy_strategy = Evaluator.sensitive_therapy_strategy

    resistance = Evaluator.resistance
    resistance_therapy_name = Evaluator.resistance_therapy_name
    resistance_therapy_strategy = Evaluator.resistance_therapy_strategy

    @classmethod
    def get_union_strategies(cls, sensitive, resistance):
        strategies = sorted(list(set(sensitive + resistance)))
        if '' in strategies:
            strategies.remove('')
        return strategies

    @classmethod
    def iterate_over_strategies(cls, source_df, strategy_column, therapy_column, target_df, target_row):
        for label, group in source_df.groupby(strategy_column):
            if label == '':
                continue
            therapies = group[therapy_column].dropna().drop_duplicates().sort_values()
            therapies_list = cls.list_to_string(therapies.tolist(), ', ')
            target_df.loc[target_row, label] = therapies_list
        return target_df

    @staticmethod
    def list_to_string(list_of_elements, delimiter):
        return f'{delimiter}'.join(list_of_elements)

    @classmethod
    def report_therapy_strategies(cls, dataframe):
        sensitive_strategies = cls.series_to_list(dataframe.loc[:, cls.sensitive_therapy_strategy].dropna())
        resistance_strategies = cls.series_to_list(dataframe.loc[:, cls.resistance_therapy_strategy].dropna())
        union_strategies = cls.get_union_strategies(sensitive_strategies, resistance_strategies)

        df = pd.DataFrame('', index=[cls.sensitivity, cls.resistance], columns=union_strategies)
        df = cls.iterate_over_strategies(dataframe,
                                         cls.sensitive_therapy_strategy, cls.sensitive_therapy_name,
                                         df, cls.sensitivity)
        df = cls.iterate_over_strategies(dataframe,
                                         cls.resistance_therapy_strategy, cls.resistance_therapy_name,
                                         df, cls.resistance)
        return df

    @staticmethod
    def series_to_list(series):
        return series.dropna().tolist()
