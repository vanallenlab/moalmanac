import numpy as np
import pandas as pd
import scipy.stats as stats
import os
import subprocess

from datasources import Lawrence
from reader import Reader
from illustrator import Signatures

from config import COLNAMES
from config import CONFIG

FUSION_CONFIG = CONFIG['fusion']
SEG_CONFIG = CONFIG['seg']
SIGNATURES_CONFIG = CONFIG['signatures']


class Features(object):
    features_section = 'features'
    feature_type = COLNAMES[features_section]['feature_type']
    feature = COLNAMES[features_section]['feature']
    chr = COLNAMES[features_section]['chr']
    start = COLNAMES[features_section]['start']
    end = COLNAMES[features_section]['end']
    ref = COLNAMES[features_section]['ref']
    allele1 = COLNAMES[features_section]['allele1']
    allele2 = COLNAMES[features_section]['allele2']
    alt_type = COLNAMES[features_section]['alt_type']
    alt = COLNAMES[features_section]['alt']
    ref_count = COLNAMES[features_section]['ref_count']
    alt_count = COLNAMES[features_section]['alt_count']
    tumor_f = COLNAMES[features_section]['tumor_f']
    coverage = COLNAMES[features_section]['coverage']
    patient = COLNAMES[features_section]['patient']
    segment_mean = COLNAMES[features_section]['segment_mean']
    spanningfrags = COLNAMES[features_section]['spanningfrags']
    left_gene = COLNAMES[features_section]['left_gene']
    left_chr = COLNAMES[features_section]['left_chr']
    left_start = COLNAMES[features_section]['left_start']
    right_gene = COLNAMES[features_section]['right_gene']
    right_chr = COLNAMES[features_section]['right_chr']
    right_start = COLNAMES[features_section]['right_start']
    tumor = COLNAMES[features_section]['tumor']
    normal = COLNAMES[features_section]['normal']

    annotation_map = {
        'Missense_Mutation': 'Missense',
        'Nonsense_Mutation': 'Nonsense',
        'Nonstop_Mutation': 'Nonstop',
        'Splice_Site': 'Splice Site',
        'Frame_Shift_Ins': 'Frameshift',
        'Frame_Shift_Del': 'Frameshift',
        'In_Frame_Ins': 'Insertion',
        'In_Frame_Del': 'Deletion'
    }

    coding_classifications_raw = []
    coding_classifications_mapped = []
    for key, value in annotation_map.items():
        coding_classifications_raw.append(key)
        coding_classifications_mapped.append(value)

    high_burden = 'High Mutational Burden'
    not_high_burden = 'Not high mutational burden'

    feature_columns = [feature_type, feature, chr, start, end, ref, allele1, allele2, alt_type, alt,
                       ref_count, alt_count, tumor_f, coverage, tumor, normal,
                       spanningfrags, left_gene, right_gene, right_chr, right_start]

    maps_section = 'maps'
    maps_columns = COLNAMES[maps_section].values()

    bin_section = 'bin_names'
    bin_columns = COLNAMES[bin_section].values()

    all_columns = list(set(list(maps_columns) + list(bin_columns)))

    @staticmethod
    def annotate_feature_type(feature_type, idx):
        return pd.Series(feature_type, index=idx)

    @staticmethod
    def calculate_percentile(array, percentile):
        return np.percentile(array.astype(float), float(percentile))

    @classmethod
    def create_empty_dataframe(cls):
        return pd.DataFrame(columns=cls.all_columns)

    @classmethod
    def create_empty_series(cls):
        return pd.Series(index=cls.all_columns, dtype=str)

    @classmethod
    def drop_duplicate_genes(cls, df, sort_column):
        return df.sort_values(sort_column, ascending=False).drop_duplicates([Features.feature], keep='first').index

    @classmethod
    def import_features(cls, handle, importer, column_config, feature_type):
        if os.path.exists(handle):
            df_accept, df_reject = importer.import_feature(handle, column_config, feature_type)
        else:
            df_accept = cls.create_empty_dataframe()
            df_reject = cls.create_empty_dataframe()
        return df_accept, df_reject

    @classmethod
    def preallocate_missing_columns(cls, df):
        missing_cols = [col for col in cls.all_columns if col not in df.columns]
        return pd.concat([df, pd.DataFrame(None, columns=missing_cols, index=df.index)], axis=1)

    @classmethod
    def rename_coding_classifications(cls, series):
        return series.replace(cls.annotation_map)

    @classmethod
    def return_idx_variants_coding(cls, series_alt_type):
        return series_alt_type.isin(cls.coding_classifications_mapped)

    @classmethod
    def return_variants_coding(cls, df):
        return df[cls.return_idx_variants_coding(df[cls.alt_type])]

    @classmethod
    def return_variants_non_coding(cls, df):
        return df[~cls.return_idx_variants_coding(df[cls.alt_type])]


class Aneuploidy(object):
    aneuploidy = 'aneuploidy'

    feature_type_section = 'feature_types'
    feature_type = CONFIG[feature_type_section][aneuploidy]

    aneuploidy_section = aneuploidy
    wgd = COLNAMES[aneuploidy]['wgd']
    wgd_string = COLNAMES[aneuploidy]['wgd_string']

    @classmethod
    def summarize(cls, boolean):
        df = Features.create_empty_dataframe()
        if boolean:
            df.loc[0, Features.feature_type] = cls.feature_type
            df.loc[0, Features.feature] = cls.wgd_string
        return df


class BurdenReader(object):
    feature_type_section = 'feature_types'
    feature_type = CONFIG[feature_type_section]['burden']
    feature_type_mutations = CONFIG[feature_type_section]['mut']

    burden_section = 'burden'
    patient_id = COLNAMES[burden_section]['patient']
    tumor_type = COLNAMES[burden_section]['tumor_type']
    ontology = COLNAMES[burden_section]['ontology']
    code = COLNAMES[burden_section]['code']
    bases_covered = COLNAMES[burden_section]['bases_covered']
    n_nonsyn_mutations = COLNAMES[burden_section]['n_nonsyn_mutations']
    mutational_burden = COLNAMES[burden_section]['mutational_burden']
    percentile_tcga = COLNAMES[burden_section]['percentile_tcga']
    percentile_tcga_tissue = COLNAMES[burden_section]['percentile_tcga_tissue']
    high_burden_boolean = COLNAMES[burden_section]['high_burden_boolean']

    high = 'High Mutational Burden'
    not_high = 'Not high mutational burden'

    @classmethod
    def calculate_burden(cls, mutations, bases_covered):
        if not np.isnan(bases_covered):
            mutations_per_base = float(mutations) / bases_covered
            mutations_per_megabase = mutations_per_base * 10**6
        else:
            mutations_per_megabase = np.nan
        return round(mutations_per_megabase, 3)

    @classmethod
    def calculate_percentile(cls, array, value):
        percentile = stats.percentileofscore(array.astype(float), float(value))
        return round(percentile, 3)

    @classmethod
    def calculate_percentile_tissuetype(cls, value, ds, ontology_code):
        if str.lower(ontology_code) in ds[cls.code].str.lower().unique().tolist():
            ds = ds[ds[cls.code] == ontology_code]
            return cls.calculate_percentile(ds[cls.mutational_burden], value)
        else:
            return np.nan

    @classmethod
    def create_burden_series(cls, patient, bases_covered):
        series = Features.create_empty_series()
        series[cls.patient_id] = patient[cls.patient_id]
        series[cls.tumor_type] = patient[cls.tumor_type]
        series[cls.ontology] = patient[cls.ontology]
        series[cls.code] = patient[cls.code]
        series[cls.bases_covered] = bases_covered
        return series

    @classmethod
    def evaluate_high_burden(cls, series):
        if np.isnan(series[cls.mutational_burden]):
            return False

        if not np.isnan(series[cls.percentile_tcga_tissue]):
            percentile = series[cls.percentile_tcga_tissue]
        else:
            percentile = series[cls.percentile_tcga]

        if (float(percentile) >= 0.80) & (float(series[cls.mutational_burden]) >= 10.0):
            return True
        else:
            return False

    @classmethod
    def evaluate_high_burden_boolean(cls, boolean):
        if boolean:
            return Features.high_burden
        else:
            return Features.not_high_burden

    @classmethod
    def import_feature(cls, handle, patient, variants, dbs):
        if os.path.exists(handle):
            bases_covered = float(Reader.read(handle, '\t', index_col=False).columns.tolist()[0])
        else:
            bases_covered = np.nan
        df = cls.create_burden_series(patient, bases_covered)
        df[Features.feature_type] = cls.feature_type

        mutations = variants[variants[Features.feature_type] == cls.feature_type_mutations].shape[0]
        mutational_burden = cls.calculate_burden(mutations, bases_covered)

        df[cls.n_nonsyn_mutations] = mutations
        df[cls.mutational_burden] = mutational_burden

        code = patient[cls.code]
        lawrence = Lawrence.import_ds(dbs)
        df[cls.percentile_tcga] = cls.calculate_percentile(lawrence[cls.mutational_burden], mutational_burden)
        df[cls.percentile_tcga_tissue] = cls.calculate_percentile_tissuetype(mutational_burden, lawrence, code)

        burden_boolean = cls.evaluate_high_burden(df)
        df[cls.high_burden_boolean] = burden_boolean
        df[Features.feature] = cls.evaluate_high_burden_boolean(burden_boolean)
        df[Features.alt] = f'{round(mutational_burden, 2)} mutations per Mb'
        return df.to_frame().T


class CopyNumber:
    @staticmethod
    def format_cn_gene(series):
        new_series = series.str.split(' ', expand=True).loc[:, 0]
        return new_series


class CalledCNAReader(CopyNumber):
    amplification = SEG_CONFIG['amp']
    deletion = SEG_CONFIG['del']

    @staticmethod
    def create_column_map(config):
        return {
            config['gene']: Features.feature,
            config['call']: Features.alt_type
        }

    @classmethod
    def filter_calls(cls, series):
        return series.fillna('').isin([cls.amplification, cls.deletion])

    @classmethod
    def import_feature(cls, handle, column_names, feature_type):
        column_map = cls.create_column_map(column_names)
        df = Reader.safe_read(handle, '\t', column_map)
        df = Features.preallocate_missing_columns(df)

        if not df.empty:
            df[Features.feature_type] = Features.annotate_feature_type(feature_type, df.index)
            df[Features.feature] = cls.format_cn_gene(df[Features.feature])

            idx = cls.filter_calls(df[Features.alt_type])
            seg_accept = df[idx]
            seg_reject = df[~idx]
        else:
            seg_accept = Features.create_empty_dataframe()
            seg_reject = Features.create_empty_dataframe()
        return seg_accept, seg_reject


class CNVReader(CopyNumber):
    amp_percentile = SEG_CONFIG['amp_percentile']
    del_percentile = SEG_CONFIG['del_percentile']
    amplification = SEG_CONFIG['amp']
    deletion = SEG_CONFIG['del']

    @classmethod
    def annotate_amp_del(cls, idx, idx_amp, idx_del):
        series = pd.Series('', index=idx)
        series[idx_amp] = cls.amplification
        series[idx_del] = cls.deletion
        return series

    @staticmethod
    def create_column_map(config):
        return {
            config['gene']: Features.feature,
            config['contig']: Features.chr,
            config['start']: Features.start,
            config['end']: Features.end,
            config['mean']: Features.segment_mean,
            config['sample']: Features.tumor,
        }

    @classmethod
    def drop_duplicate_genes(cls, df):
        return (
            df
            .sort_values(Features.segment_mean, ascending=False)
            .drop_duplicates([Features.feature], keep='first')
            .index
        )

    @classmethod
    def filter_by_threshold(cls, df, percentile_amp, percentile_del):
        unique_segments = cls.get_unique_segments(df)
        threshold_amp = Features.calculate_percentile(unique_segments, percentile_amp)
        threshold_del = Features.calculate_percentile(unique_segments, percentile_del)

        idx_amp = df[df[Features.segment_mean].astype(float) >= float(threshold_amp)].index
        idx_del = df[df[Features.segment_mean].astype(float) <= float(threshold_del)].index

        df[Features.alt_type] = cls.annotate_amp_del(df.index, idx_amp, idx_del)
        idx_accept = df[df[Features.alt_type] != ''].index
        idx_unique = Features.drop_duplicate_genes(df.loc[idx_accept, :], Features.segment_mean)

        df[Features.segment_mean] = df[Features.segment_mean].astype(float).round(3)
        df_accept = df.loc[idx_unique, :]
        df_reject = df.loc[df.index.difference(df_accept.index), :]

        return df_accept, df_reject

    @staticmethod
    def get_unique_segments(df):
        return df.drop_duplicates([Features.chr, Features.start])[Features.segment_mean]

    @classmethod
    def import_feature(cls, handle, column_names, feature_type):
        column_map = cls.create_column_map(column_names)
        df = Reader.safe_read(handle, '\t', column_map, comment_character='#')
        df = Features.preallocate_missing_columns(df)

        if not df.empty:
            df[Features.feature_type] = Features.annotate_feature_type(feature_type, df.index)
            df[Features.feature] = cls.format_cn_gene(df[Features.feature])

            seg_accept, seg_reject = cls.filter_by_threshold(df, cls.amp_percentile, cls.del_percentile)
        else:
            seg_accept = Features.create_empty_dataframe()
            seg_reject = Features.create_empty_dataframe()

        return seg_accept, seg_reject


class CoverageMetrics(object):
    @staticmethod
    def get_dnp_boolean(series):
        return series.astype(str).str.contains('\|')

    @staticmethod
    def split_counts(series_alt_count):
        return series_alt_count.str.split('\|', expand=True).fillna(np.nan)

    @staticmethod
    def calculate_coverage(series_alt_count, series_ref_count):
        return series_alt_count.astype(float).add(series_ref_count.astype(float))

    @classmethod
    def calculate_coverage_dnp(cls, series_alt_count, series_ref_count):
        df_alt_split = cls.split_counts(series_alt_count)
        added = df_alt_split.astype(float).add(series_ref_count.astype(float), axis=0)
        return added.apply(lambda x: '|'.join(x.dropna().astype(int).astype(str)), axis=1)

    @classmethod
    def append_coverage(cls, series_alt_count, series_ref_count):
        idx = cls.get_dnp_boolean(series_alt_count)
        total = pd.Series(0, index=idx.index)

        idx_snp = series_alt_count[~idx].index
        idx_dnp = series_alt_count[idx].index

        if not idx_snp.empty:
            total[idx_snp] = cls.calculate_coverage(series_alt_count[idx_snp], series_ref_count[idx_snp])
        if not idx_dnp.empty:
            total[idx_dnp] = cls.calculate_coverage_dnp(series_alt_count[idx_dnp], series_ref_count[idx_dnp])

        return total

    @staticmethod
    def calculate_tumor_f(series_alt_count, series_total_coverage):
        return series_alt_count.astype(float).divide(series_total_coverage.astype(float))

    @classmethod
    def calculate_tumor_f_dnp(cls, series_alt_count, series_total_coverage):
        df_alt_split = cls.split_counts(series_alt_count)
        df_cov_split = cls.split_counts(series_total_coverage)
        divided = df_alt_split.astype(float).divide(df_cov_split.astype(float), axis=0)
        return divided.apply(lambda x: '|'.join(x.dropna().astype(int).astype(str)), axis=1)

    @classmethod
    def append_tumor_f(cls, series_alt_count, series_total_coverage):
        idx = cls.get_dnp_boolean(series_alt_count)
        tumor_f = pd.Series(0, index=idx.index)

        idx_snp = series_alt_count[~idx].index
        idx_dnp = series_alt_count[idx].index

        if not idx_snp.empty:
            tumor_f[idx_snp] = cls.calculate_tumor_f(series_alt_count[idx_snp], series_total_coverage[idx_snp])
        if not idx_dnp.empty:
            tumor_f[idx_dnp] = cls.calculate_tumor_f_dnp(series_alt_count[idx_dnp], series_total_coverage[idx_dnp])

        return tumor_f.round(4)

    @staticmethod
    def format_coverage_col(series):
        formatted_series = series.replace('__UNKNOWN__', np.nan)
        formatted_series = formatted_series.fillna(np.nan)
        return formatted_series


class CosmicSignatures(object):
    feature_type_section = 'feature_types'
    feature_type = CONFIG[feature_type_section]['signature']

    signature_section = 'signatures'
    patient_id = COLNAMES[signature_section]['patient']

    min_contribution = SIGNATURES_CONFIG['min_contribution']

    input_section = 'input_data'
    input_sample = COLNAMES[input_section]['tumor']
    input_ref = COLNAMES[input_section]['ref']
    input_alt = COLNAMES[input_section]['allele2']
    input_chr = COLNAMES[input_section]['chr']
    input_pos = COLNAMES[input_section]['start']

    @classmethod
    def calculate_contributions(cls, patient_id, snv_handle):
        if os.path.exists(snv_handle):
            cls.run_deconstructsigs(patient_id, snv_handle)

    @classmethod
    def import_feature(cls, snv_handle, patient):
        patient_id = patient[cls.patient_id]

        cls.calculate_contributions(patient_id, snv_handle)
        context_handle = cls.create_handle(patient_id, 'sigs.context.txt')
        weights_handle = cls.create_handle(patient_id, 'sigs.cosmic.txt')

        if os.path.exists(context_handle):
            contexts = Reader.read(context_handle, '\t')
            cls.illustrate_contexts(patient_id, contexts.loc[0, :])

        if os.path.exists(weights_handle):
            weights = Reader.read(weights_handle, '\t')
            series_weights = cls.format_weights(weights.iloc[0, :30])
        else:
            series_weights = pd.Series()

        series_significant_weights = cls.subset_significant_signatures(series_weights)
        df = cls.create_feature_dataframe(patient, series_significant_weights)
        return df

    @classmethod
    def create_feature_dataframe(cls, patient, series_weights):
        dataframe = Features.create_empty_dataframe()
        dataframe[Features.feature] = series_weights.index
        dataframe[Features.alt] = series_weights.values
        dataframe[Features.alt_type] = 'version 2'
        dataframe.loc[:, Features.feature_type] = cls.feature_type
        dataframe.loc[:, cls.patient_id] = patient[cls.patient_id]
        return dataframe

    @staticmethod
    def create_handle(patient_id, suffix):
        return '.'.join([patient_id, suffix])

    @staticmethod
    def format_weights(weights):
        nonzero = weights[weights.astype(float) != 0.0]
        nonzero.index = nonzero.index.str.replace('weights', 'COSMIC').str.replace('.', ' ').tolist()
        return nonzero.astype(float).round(3)

    @staticmethod
    def illustrate_contexts(patient_id, contexts):
        Signatures.generate_context_plot(patient_id, contexts)
        Signatures.generate_context_plot_normalized(patient_id, contexts)

    @classmethod
    def run_deconstructsigs(cls, patient_id, snv_handle):
        cmd = " ".join(["bash", "wrapper_deconstructsigs.sh", patient_id, snv_handle,
                        cls.input_sample, cls.input_ref, cls.input_alt, cls.input_chr, cls.input_pos])
        cls.run_subprocess(cmd)

    @staticmethod
    def run_subprocess(cmd):
        subprocess.call(cmd, shell=True)

    @classmethod
    def subset_significant_signatures(cls, series):
        return series[series.astype(float) >= float(cls.min_contribution)]


class FusionReader(object):
    leftbreakpoint = FUSION_CONFIG['leftbreakpoint']
    rightbreakpoint = FUSION_CONFIG['rightbreakpoint']
    alt_type = FUSION_CONFIG['alt_type']

    @classmethod
    def create_colmap(cls, config):
        return {
            config['name']: Features.feature,
            config['spanningfrags']: Features.spanningfrags,
            config[cls.leftbreakpoint]: cls.leftbreakpoint,
            config[cls.rightbreakpoint]: cls.rightbreakpoint
        }

    @staticmethod
    def filter_by_spanningfrags(series_spanningfrags, min_spanningfrags):
        return series_spanningfrags[series_spanningfrags.astype(float) >= float(min_spanningfrags)].index

    @classmethod
    def import_feature(cls, handle, colnames, feature_type):
        column_map = cls.create_colmap(colnames)
        df = Reader.safe_read(handle, '\t', column_map, index_col=False)
        df = Features.preallocate_missing_columns(df)
        if not df.empty:
            split_genes = cls.split_genes(df[Features.feature])
            df[Features.left_gene] = split_genes[Features.left_gene]
            df[Features.right_gene] = split_genes[Features.right_gene]

            left = cls.split_breakpoint(df[cls.leftbreakpoint])
            right = cls.split_breakpoint(df[cls.rightbreakpoint])

            df[Features.chr] = left[Features.chr]
            df[Features.start] = left[Features.start]
            df[Features.left_chr] = left[Features.chr]
            df[Features.left_start] = left[Features.start]
            df[Features.right_chr] = right[Features.chr]
            df[Features.right_start] = right[Features.start]
            df.drop([cls.leftbreakpoint, cls.rightbreakpoint], axis=1, inplace=True)

            df[Features.feature_type] = Features.annotate_feature_type(feature_type, df.index)
            df[Features.alt_type] = cls.alt_type
            df[Features.alt] = df[Features.feature]

            idx_min_spanningfrags = cls.filter_by_spanningfrags(df[Features.spanningfrags],
                                                                FUSION_CONFIG['spanningfrags_min'])
            idx_unique = Features.drop_duplicate_genes(df.loc[idx_min_spanningfrags, :], Features.feature)
            fusions_unique = df.loc[idx_unique, :]

            fusions_left = fusions_unique.copy(deep=True)
            fusions_right = fusions_unique.copy(deep=True)
            fusions_left[Features.feature] = fusions_left[Features.left_gene]
            fusions_right[Features.feature] = fusions_right[Features.right_gene]

            fusions_accept = pd.concat([fusions_left, fusions_right], ignore_index=True)
            fusions_reject = df.loc[df.index.difference(fusions_accept.index), :]
        else:
            fusions_accept = Features.create_empty_dataframe()
            fusions_reject = Features.create_empty_dataframe()
        return fusions_accept, fusions_reject

    @staticmethod
    def split_genes(series_gene):
        return series_gene.str.split('--', expand=True).rename(
            columns={0: Features.left_gene, 1: Features.right_gene})

    @staticmethod
    def split_breakpoint(series_breakpoint):
        return series_breakpoint.str.split(':', expand=True).rename(
            columns={0: Features.chr, 1: Features.start}).loc[:, [Features.chr, Features.start]]


class GermlineReader(object):
    @classmethod
    def import_feature(cls, handle, colnames, feature_type):
        column_map = MutationReader.create_colmap(colnames)
        df = Reader.safe_read(handle, '\t', column_map, comment_character='#')
        df = Features.preallocate_missing_columns(df)

        if not df.empty:
            df[Features.feature_type] = Features.annotate_feature_type(feature_type, df.index)
            df[Features.alt_count] = CoverageMetrics.format_coverage_col(df[Features.alt_count])
            df[Features.ref_count] = CoverageMetrics.format_coverage_col(df[Features.ref_count])
            df[Features.coverage] = CoverageMetrics.append_coverage(df[Features.alt_count], df[Features.ref_count])
            df[Features.tumor_f] = CoverageMetrics.append_tumor_f(df[Features.alt_count], df[Features.coverage])
            df[Features.alt_type] = Features.rename_coding_classifications(df[Features.alt_type])

            coding_variants = Features.return_variants_coding(df)
        else:
            coding_variants = Features.create_empty_dataframe()

        return coding_variants, Features.create_empty_dataframe()


class MicrosatelliteReader(object):
    microsatellite = 'microsatellite'
    feature_type_section = 'feature_types'
    feature_type = CONFIG[feature_type_section][microsatellite]

    microsatellite_section = microsatellite
    msih = COLNAMES[microsatellite_section]['msih']
    msil = COLNAMES[microsatellite_section]['msil']
    mss = COLNAMES[microsatellite_section]['mss']
    unk = COLNAMES[microsatellite_section]['unk']

    status_map = {
        'msih': msih,
        'msil': msil,
        'mss': mss,
        'unk': unk
    }

    @classmethod
    def map_status(cls, status):
        return cls.status_map[status]

    @classmethod
    def summarize(cls, status):
        df = Features.create_empty_dataframe()
        df.loc[0, Features.feature_type] = cls.feature_type
        df.loc[0, Features.feature] = cls.map_status(status)
        return df


class MutationReader(object):
    @staticmethod
    def create_colmap(config):
        return {
            config['gene']: Features.feature,
            config['chr']: Features.chr,
            config['start']: Features.start,
            config['end']: Features.end,
            config['ref']: Features.ref,
            config['allele1']: Features.allele1,
            config['allele2']: Features.allele2,
            config['alt_type']: Features.alt_type,
            config['alt']: Features.alt,
            config['ref_count']: Features.ref_count,
            config['alt_count']: Features.alt_count,
            config['tumor']: Features.tumor,
            config['normal']: Features.normal
        }

    @classmethod
    def import_feature(cls, handle, colnames, feature_type, config=''):
        column_map = cls.create_colmap(colnames)
        df = Reader.safe_read(handle, '\t', column_map, comment_character='#')
        df = Features.preallocate_missing_columns(df)
        if not df.empty:
            df[Features.feature_type] = Features.annotate_feature_type(feature_type, df.index)
            df[Features.alt_count] = CoverageMetrics.format_coverage_col(df[Features.alt_count])
            df[Features.ref_count] = CoverageMetrics.format_coverage_col(df[Features.ref_count])
            df[Features.coverage] = CoverageMetrics.append_coverage(df[Features.alt_count], df[Features.ref_count])
            df[Features.tumor_f] = CoverageMetrics.append_tumor_f(df[Features.alt_count], df[Features.coverage])
            df[Features.alt_type] = Features.rename_coding_classifications(df[Features.alt_type])

            coding_variants = Features.return_variants_coding(df)
            non_coding_variants = Features.return_variants_non_coding(df)
        else:
            coding_variants = Features.create_empty_dataframe()
            non_coding_variants = Features.create_empty_dataframe()
        return coding_variants, non_coding_variants


class ValidationMutationReader(Features):
    section = 'validation_sequencing'
    gene = COLNAMES[section]['gene']
    alt_type = COLNAMES[section]['alt_type']
    alt = COLNAMES[section]['alt']
    validation_tumor_f = COLNAMES[section]['validation_tumor_f']
    validation_coverage = COLNAMES[section]['validation_coverage']

    column_map = {
        Features.tumor_f: validation_tumor_f,
        Features.coverage: validation_coverage
    }

    columns = [gene, alt, alt_type, validation_tumor_f, validation_coverage]

    @classmethod
    def import_feature(cls, handle, colnames, feature_type, config=''):
        df, df_reject = MutationReader.import_feature(handle, colnames, feature_type)
        df = df.drop(df.columns[df.columns.str.contains('validation')], axis=1)
        df = df.rename(columns=cls.column_map).loc[:, cls.columns]
        return df, df_reject
