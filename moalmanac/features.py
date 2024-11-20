import numpy as np
import pandas as pd
import scipy.stats as stats
import os
import subprocess
import sys

from datasources import Lawrence
from reader import Reader

from config import COLNAMES


class Features:
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
    def concat_list_of_dataframes(cls, list_of_dataframes, ignore_index=True):
        # Remove empty dataframes
        non_empty_dataframes = [dataframe for dataframe in list_of_dataframes if not dataframe.empty]
        # Remove dataframes with only null values
        dataframes = [dataframe for dataframe in non_empty_dataframes if not dataframe.isnull().all().all()]

        if len(dataframes) > 1:
            # Trim each dataframe to remove columns with all NA values
            trimmed_dataframes = []
            for dataframe in dataframes:
                na_columns = dataframe.columns[dataframe.isnull().all()]
                dataframe = dataframe.drop(columns=na_columns)
                trimmed_dataframes.append(dataframe)

            # Concatenate the trimmed, nonempty, and nonnull dataframes
            # The resulting dataframe will have columns that are an intersection of all concat'd dataframe columns
            # https://pandas.pydata.org/docs/reference/api/pandas.concat.html
            dataframe = pd.concat(trimmed_dataframes, ignore_index=ignore_index)

            # Repopulate removed columns
            return cls.preallocate_missing_columns(dataframe)
        elif len(non_empty_dataframes) == 1:
            return non_empty_dataframes[0]
        else:
            return cls.create_empty_dataframe()

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
    def import_if_path_exists(cls, handle, delimiter, column_map, **kwargs):
        if os.path.exists(handle):
            df = Reader.safe_read(handle, delimiter, column_map, **kwargs)
            return cls.preallocate_missing_columns(df)
        else:
            return cls.create_empty_dataframe()

    @classmethod
    def preallocate_missing_columns(cls, df):
        missing_cols = [col for col in cls.all_columns if col not in df.columns]
        return pd.concat([df, pd.DataFrame(None, columns=missing_cols, index=df.index)], axis=1)


class Aneuploidy:
    aneuploidy = 'aneuploidy'
    aneuploidy_section = aneuploidy
    wgd = COLNAMES[aneuploidy]['wgd']
    wgd_string = COLNAMES[aneuploidy]['wgd_string']

    @classmethod
    def summarize(cls, boolean, config):
        df = Features.create_empty_dataframe()
        feature_type = config['feature_types']['aneuploidy']
        if boolean:
            df.loc[0, Features.feature_type] = feature_type
            df.loc[0, Features.feature] = cls.wgd_string
        return df


class BurdenReader:
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
    def import_feature(cls, handle, patient, variants, dbs, config):
        if os.path.exists(handle):
            bases_covered = float(Reader.read(handle, '\t', index_col=False).columns.tolist()[0])
        else:
            bases_covered = np.nan
        df = cls.create_burden_series(patient, bases_covered)

        biomarker_type = config['feature_types']['burden']
        df[Features.feature_type] = biomarker_type
        mutations = variants[variants[Features.feature_type] == config['feature_types']['mut']].shape[0]
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

    @classmethod
    def import_feature(cls, called_handle, not_called_handle, config):
        if called_handle:
            column_map = CopyNumberCalled.create_column_map()
            handle = called_handle
        else:
            column_map = CopyNumberTotal.create_column_map()
            handle = not_called_handle

        df = Features.import_if_path_exists(handle, '\t', column_map, comment_character="#")

        amplification_string = config['seg']['amp']
        deletion_string = config['seg']['del']
        if not df.empty:
            biomarker_type = config['feature_types']['cna']
            df[Features.feature_type] = Features.annotate_feature_type(biomarker_type, df.index)
            df[Features.feature] = cls.format_cn_gene(df[Features.feature])
            if called_handle:
                seg_accept, seg_reject = CopyNumberCalled.process_calls(df, amplification_string, deletion_string)
            else:
                seg_accept, seg_reject = CopyNumberTotal.process_calls(df, config)
        else:
            seg_accept = Features.create_empty_dataframe()
            seg_reject = Features.create_empty_dataframe()
        return seg_accept, seg_reject


class CopyNumberCalled(CopyNumber):
    @staticmethod
    def create_column_map():
        section = 'called_cn_input'
        column_names = COLNAMES[section]
        return {
            column_names['gene']: Features.feature,
            column_names['call']: Features.alt_type
        }

    @classmethod
    def filter_calls(cls, series, amp_string, del_string):
        return series.fillna('').isin([amp_string, del_string])

    @classmethod
    def process_calls(cls, dataframe, amp_string, del_string):
        idx = cls.filter_calls(dataframe[Features.alt_type], amp_string, del_string)
        return dataframe[idx], dataframe[~idx]


class CopyNumberTotal(CopyNumber):
    @classmethod
    def annotate_amp_del(cls, idx, idx_amp, idx_del, amp_string, del_string):
        series = pd.Series('', index=idx)
        series[idx_amp] = amp_string
        series[idx_del] = del_string
        return series

    @staticmethod
    def create_column_map():
        section = 'seg_input'
        column_names = COLNAMES[section]
        return {
            column_names['gene']: Features.feature,
            column_names['contig']: Features.chr,
            column_names['start']: Features.start,
            column_names['end']: Features.end,
            column_names['mean']: Features.segment_mean,
            column_names['sample']: Features.tumor,
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
    def filter_by_threshold(cls, df, percentile_amp, percentile_del, amp_string, del_string):
        unique_segments = cls.get_unique_segments(df)
        threshold_amp = Features.calculate_percentile(unique_segments, percentile_amp)
        threshold_del = Features.calculate_percentile(unique_segments, percentile_del)

        idx_amp = df[df[Features.segment_mean].astype(float) >= float(threshold_amp)].index
        idx_del = df[df[Features.segment_mean].astype(float) <= float(threshold_del)].index

        df[Features.alt_type] = cls.annotate_amp_del(df.index, idx_amp, idx_del, amp_string, del_string)
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
    def process_calls(cls, dataframe, config):
        amp_percentile = config['seg']['amp_percentile']
        del_percentile = config['seg']['del_percentile']
        amp_string = config['seg']['amp']
        del_string = config['seg']['del']
        return cls.filter_by_threshold(dataframe, amp_percentile, del_percentile, amp_string, del_string)


class CoverageMetrics:
    @classmethod
    def apply_min_coverage_for_onps(cls, series):
        idx = cls.get_onp_boolean(series)
        split_counts = cls.split_counts(series.loc[idx])
        min_counts = cls.get_minimum_values_dataframe(split_counts, axis=1)
        series.loc[idx] = min_counts.loc[idx]
        return series

    @staticmethod
    def calculate_coverage(series_alt_count, series_ref_count):
        series_alt_count = pd.to_numeric(series_alt_count).fillna(0)
        series_ref_count = pd.to_numeric(series_ref_count).fillna(0)
        return series_alt_count.add(series_ref_count).dropna().astype(int)

    @staticmethod
    def calculate_tumor_f(series_alt_count, series_total_coverage):
        series_alt_count = pd.to_numeric(series_alt_count)
        series_total_coverage = pd.to_numeric(series_total_coverage)
        return series_alt_count.divide(series_total_coverage).dropna().astype(float).round(4)

    @classmethod
    def convert_to_pandas_int64(cls, series):
        return (
            series
            .apply(cls.safe_cast)
            .astype("Int64")
        )

    @classmethod
    def format_coverage_col(cls, series):
        # this is required for python 3.12 and pandas 2.2.2 to opt into future behavior for type downcasting
        with pd.option_context("future.no_silent_downcasting", True):
            formatted_series = (
                series
                .replace('__UNKNOWN__', pd.NA)
                .replace('', pd.NA)
                .astype(object)
                .fillna(pd.NA)
            )
        formatted_series = cls.apply_min_coverage_for_onps(formatted_series)
        formatted_series = cls.convert_to_pandas_int64(formatted_series)
        return formatted_series

    @staticmethod
    def get_minimum_values_dataframe(dataframe, axis=1):
        return (
            dataframe
            .replace('', pd.NA)
            .fillna(pd.NA)
            .dropna(how='all')
            .apply(lambda x: x.dropna().min(), axis=axis)
            .astype(int)
        )

    @staticmethod
    def get_onp_boolean(series):
        return series.astype(str).str.contains('|', regex=False)

    @staticmethod
    def safe_cast(value):
        try:
            return int(value)
        except (ValueError, TypeError):
            return pd.NA

    @staticmethod
    def split_counts(series):
        return series.fillna('').astype(str).str.split('|', expand=True, regex=False)


class CosmicSignatures:
    signature_section = 'signatures'
    patient_id = COLNAMES[signature_section]['patient']

    input_section = 'mutational_signature_input'
    input_signature = COLNAMES[input_section]['signature']
    input_contribution = COLNAMES[input_section]['contribution']

    @classmethod
    def create_column_map(cls):
        section = 'mutational_signature_input'
        column_names = COLNAMES[section]
        return {
            column_names[cls.input_signature]: Features.feature,
            column_names[cls.input_contribution]: Features.alt
        }

    @classmethod
    def import_feature(cls, path, config):
        """Loads and formats Cosmic Mutational Signatures based on provided file path."""
        column_map = cls.create_column_map()
        df = Features.import_if_path_exists(path, delimiter='\t', column_map=column_map)
        if not df.empty:
            biomarker_type = config['feature_types']['signature']
            minimum_contribution = config['signatures']['min_contribution']
            df[Features.feature_type] = biomarker_type
            df[Features.alt_type] = 'v3.4'
            df[Features.alt] = cls.round_contributions(df[Features.alt])
            idx = cls.index_for_minimum_contribution(series=df[Features.alt], minimum_value=minimum_contribution)
            return df[idx]
        else:
            return Features.create_empty_dataframe()

    @classmethod
    def index_for_minimum_contribution(cls, series, minimum_value=0.06):
        """Subsets the provided SBS signatures to those that pass the minimum contribution, specified in config.ini"""
        return series.astype(float) >= float(minimum_value)

    @staticmethod
    def round_contributions(series, decimals=3):
        """Rounds a pandas series of float values to the specified number of decimal places, 3 by default."""
        return series.astype(float).round(decimals)


class Fusion:
    @classmethod
    def create_colmap(cls, config):
        section = 'fusion_input'
        column_names = COLNAMES[section]

        leftbreakpoint = config['fusion']['leftbreakpoint']
        rightbreakpoint = config['fusion']['rightbreakpoint']
        return {
            column_names['name']: Features.feature,
            column_names['spanningfrags']: Features.spanningfrags,
            column_names[leftbreakpoint]: leftbreakpoint,
            column_names[rightbreakpoint]: rightbreakpoint
        }

    @staticmethod
    def filter_by_spanning_fragment_count(series, minimum=5.0):
        minimum = int(float(minimum))
        return series[series.astype(int).ge(minimum)].index

    @classmethod
    def import_feature(cls, handle, config):
        column_map = cls.create_colmap(config)
        df = Features.import_if_path_exists(handle, '\t', column_map, index_col=False)
        if not df.empty:
            split_genes = cls.split_genes(df[Features.feature])
            df[Features.left_gene] = split_genes[Features.left_gene]
            df[Features.right_gene] = split_genes[Features.right_gene]

            leftbreakpoint = config['fusion']['leftbreakpoint']
            rightbreakpoint = config['fusion']['rightbreakpoint']
            left = cls.split_breakpoint(df[leftbreakpoint])
            right = cls.split_breakpoint(df[rightbreakpoint])

            df[Features.chr] = left[Features.chr]
            df[Features.start] = left[Features.start]
            df[Features.left_chr] = left[Features.chr]
            df[Features.left_start] = left[Features.start]
            df[Features.right_chr] = right[Features.chr]
            df[Features.right_start] = right[Features.start]
            df.drop([leftbreakpoint, rightbreakpoint], axis=1, inplace=True)

            biomarker_type = config['feature_types']['fusion']
            df[Features.feature_type] = Features.annotate_feature_type(biomarker_type, df.index)
            df[Features.alt_type] = config['fusion']['alt_Type']
            df[Features.alt] = df[Features.feature]

            min_fragments = config['fusion']['spanningfrags_min']
            idx_min_spanning_fragments = cls.filter_by_spanning_fragment_count(
                series=df[Features.spanningfrags],
                minimum=min_fragments
            )
            idx_unique = Features.drop_duplicate_genes(df.loc[idx_min_spanning_fragments, :], Features.feature)
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
        return (
            series_gene
            .str.split('--', expand=True)
            .rename(columns={0: Features.left_gene, 1: Features.right_gene})
        )

    @staticmethod
    def split_breakpoint(series_breakpoint):
        return (
            series_breakpoint
            .str.split(':', expand=True)
            .rename(columns={0: Features.chr, 1: Features.start})
            .loc[:, [Features.chr, Features.start]]
        )


class MicrosatelliteReader:
    microsatellite_section = 'microsatellite'
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
    def summarize(cls, status, config):
        df = Features.create_empty_dataframe()
        biomarker_type = config['feature_types']['microsatellite']
        df.loc[0, Features.feature_type] = biomarker_type
        df.loc[0, Features.feature] = cls.map_status(status)
        return df


class MAF(Features):
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

    @staticmethod
    def check_format(df):
        if "Protein_Change" in df.columns:
            return "tcga_maf_input"
        elif "HGVSp_Short" in df.columns:
            return "gdc_maf_input"
        else:
            sys.exit("Neither 'Protein_Change' nor 'HGVSp_Short' are present columns, cannot map to MAF format")

    @classmethod
    def create_column_map(cls, maf_format):
        column_names = COLNAMES[maf_format]
        return {
            column_names['gene']: cls.feature,
            column_names['chr']: cls.chr,
            column_names['start']: cls.start,
            column_names['end']: cls.end,
            column_names['ref']: cls.ref,
            column_names['allele1']: cls.allele1,
            column_names['allele2']: cls.allele2,
            column_names['alt_type']: cls.alt_type,
            column_names['alt']: cls.alt,
            column_names['ref_count']: cls.ref_count,
            column_names['alt_count']: cls.alt_count,
            column_names['tumor']: cls.tumor,
            column_names['normal']: cls.normal
        }

    @classmethod
    def format_maf(cls, df, feature_type):
        df = Features.preallocate_missing_columns(df)
        df[Features.feature_type] = cls.annotate_feature_type(feature_type, df.index)
        df[Features.alt_count] = CoverageMetrics.format_coverage_col(df[cls.alt_count])
        df[Features.ref_count] = CoverageMetrics.format_coverage_col(df[cls.ref_count])
        df[Features.coverage] = CoverageMetrics.calculate_coverage(df[cls.alt_count], df[cls.ref_count])
        df[Features.tumor_f] = CoverageMetrics.calculate_tumor_f(df[cls.alt_count], df[cls.coverage])
        df[Features.alt_type] = cls.rename_coding_classifications(df[cls.alt_type])
        return df

    @classmethod
    def import_maf(cls, handle):
        if os.path.exists(handle):
            df_mini = Reader.read(handle, delimiter='\t', nrows=2, comment='#')
            maf_format = cls.check_format(df_mini)
            column_map = cls.create_column_map(maf_format)
            df = Reader.safe_read(handle, '\t', column_map, comment_character='#')
            return df
        else:
            return cls.create_empty_dataframe()

    @classmethod
    def rename_coding_classifications(cls, series):
        return series.astype(str).replace(cls.annotation_map)

    @classmethod
    def return_idx_variants_coding(cls, series_alt_type):
        return series_alt_type.isin(cls.coding_classifications_mapped)

    @classmethod
    def return_variants_coding(cls, df):
        return df[cls.return_idx_variants_coding(df[cls.alt_type])]

    @classmethod
    def return_variants_non_coding(cls, df):
        return df[~cls.return_idx_variants_coding(df[cls.alt_type])]


class MAFGermline(MAF):
    @classmethod
    def import_feature(cls, handle, config):
        df = cls.import_maf(handle)
        biomarker_type = config['feature_types']['germline']
        if not df.empty:
            df = cls.format_maf(df, biomarker_type)
            coding_variants = cls.return_variants_coding(df)
        else:
            coding_variants = cls.create_empty_dataframe()
        non_coding_variants = cls.create_empty_dataframe()
        return coding_variants, non_coding_variants


class MAFSomatic(MAF):
    @classmethod
    def import_feature(cls, handle, config):
        df = cls.import_maf(handle)
        biomarker_type = config['feature_types']['mut']
        if not df.empty:
            df = cls.format_maf(df, biomarker_type)
            coding_variants = cls.return_variants_coding(df)
            non_coding_variants = cls.return_variants_non_coding(df)
        else:
            coding_variants = cls.create_empty_dataframe()
            non_coding_variants = cls.create_empty_dataframe()
        return coding_variants, non_coding_variants


class MAFValidation(MAF):
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
    def import_feature(cls, handle, config):
        df, df_reject = MAFSomatic.import_feature(handle, config)
        df = df.drop(df.columns[df.columns.str.contains('validation')], axis=1)
        df = df.rename(columns=cls.column_map).loc[:, cls.columns]
        return df, df_reject


class Simple:
    section = 'simple_input'
    feature_type = COLNAMES[section]['feature_type']
    feature = COLNAMES[section]['feature']
    gene = COLNAMES[section]['gene']
    alt_type = COLNAMES[section]['alt_type']
    alt = COLNAMES[section]['alt']

    column_map_feature = {
        feature: feature,
        feature_type: feature_type,
        alt_type: alt_type,
        alt: alt
    }

    column_map_gene = {
        gene: feature,
        feature_type: feature_type,
        alt_type: alt_type,
        alt: alt
    }

    columns = [feature_type, gene, alt_type, alt]

    @classmethod
    def check_format(cls, df):
        if "feature" in df.columns:
            return cls.column_map_feature
        elif "gene" in df.columns:
            return cls.column_map_gene
        else:
            sys.exit("Neither 'feature' nor 'gene' are present columns, cannot read input file")

    @classmethod
    def import_feature(cls, handle):
        if os.path.exists(handle):
            df_mini = Reader.read(handle, delimiter='\t', nrows=2)
            column_map = cls.check_format(df_mini)
            df = Reader.safe_read(handle, '\t', column_map=column_map)
            return Features.preallocate_missing_columns(df)
        else:
            return Features.create_empty_dataframe()
