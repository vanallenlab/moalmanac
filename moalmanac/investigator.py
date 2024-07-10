import pandas as pd
import math
import numpy as np
import scipy.stats as stats

from datasources import Preclinical
from illustrator import PreclinicalEfficacy

from config import COLNAMES


class Investigator(object):
    summary = Preclinical.summary
    variants = Preclinical.variants
    cnas = Preclinical.cnas
    fusions = Preclinical.fusions
    gdsc = Preclinical.gdsc
    mappings = Preclinical.mappings
    gene = Preclinical.gene
    partner = Preclinical.partner

    preclinical_section = 'preclinical'
    n_wt = COLNAMES[preclinical_section]['n_wt']
    n_mut = COLNAMES[preclinical_section]['n_mut']
    wt_samples = COLNAMES[preclinical_section]['wt_samples']
    mut_samples = COLNAMES[preclinical_section]['mut_samples']

    n_cell_lines_considered = COLNAMES[preclinical_section]['n_cell_lines_considered']
    n_wt_tested = COLNAMES[preclinical_section]['n_wt_tested']
    n_mut_tested = COLNAMES[preclinical_section]['n_mut_tested']
    wt_values = COLNAMES[preclinical_section]['wt_values']
    mut_values = COLNAMES[preclinical_section]['mut_values']
    wt_median = COLNAMES[preclinical_section]['wt_median']
    mut_median = COLNAMES[preclinical_section]['mut_median']
    wt_mean = COLNAMES[preclinical_section]['wt_mean']
    mut_mean = COLNAMES[preclinical_section]['mut_mean']
    wt_std = COLNAMES[preclinical_section]['wt_std']
    mut_std = COLNAMES[preclinical_section]['mut_std']

    statistic = COLNAMES[preclinical_section]['statistic']
    pvalue = COLNAMES[preclinical_section]['pvalue']

    dtype = COLNAMES[preclinical_section]['dtype']
    feature = COLNAMES[preclinical_section]['feature']
    feature_type = COLNAMES[preclinical_section]['feature_type']
    alteration_type = COLNAMES[preclinical_section]['alteration_type']
    alteration = COLNAMES[preclinical_section]['alteration']
    feature_display = COLNAMES[preclinical_section]['feature_display']
    therapy = COLNAMES[preclinical_section]['therapy']
    therapy_mapped = COLNAMES[preclinical_section]['therapy_mapped']
    sensitive_therapy = COLNAMES[preclinical_section]['sensitive_therapy']
    base64 = COLNAMES[preclinical_section]['base64']
    patient_id = COLNAMES[preclinical_section]['patient_id']
    use_column = COLNAMES[preclinical_section]['use_column']
    model_id = COLNAMES[preclinical_section]['model_id']
    ln_ic50 = COLNAMES[preclinical_section]['ln_ic50']
    ic50 = COLNAMES[preclinical_section]['ic50']
    tested_subfeature = COLNAMES[preclinical_section]['tested_subfeature']

    @staticmethod
    def list_feature_combinations(split_feature, feature_length):
        return ['.'.join(split_feature[:i]) for i in range(1, feature_length + 1)]

    @staticmethod
    def split_string(string, delimiter):
        return string.split(delimiter)


class SensitivityDictionary(Investigator):
    rounding_places = 3
    fill_na_value = ''

    @classmethod
    def calculate_series_exp(cls, series):
        return cls.round_value(series.apply(lambda x: math.exp(x)), cls.rounding_places)

    @classmethod
    def calculate_series_mean(cls, series):
        return cls.round_value(series.mean(), cls.rounding_places)

    @classmethod
    def calculate_series_median(cls, series):
        return cls.round_value(series.median(), cls.rounding_places)

    @classmethod
    def calculate_series_std(cls, series):
        return cls.round_value(series.std(), cls.rounding_places)

    @staticmethod
    def calculate_wilcoxon_rank_sum(series1, series2):
        if (len(series1) == 0) | (len(series2) == 0):
            return np.nan, np.nan
        else:
            return stats.ranksums(series1, series2)

    @staticmethod
    def calculate_mann_whitney_u(series1, series2):
        if (len(series1) == 0) | (len(series2) == 0):
            return np.nan, np.nan
        else:
            return stats.mannwhitneyu(series1, series2, alternative='two-sided')

    @classmethod
    def create(cls, dbs, df_actionable, config):
        summary = dbs[cls.summary]
        variants = dbs[cls.variants]
        cnas = dbs[cls.cnas]
        fusions = dbs[cls.fusions]
        gdsc = dbs[cls.gdsc]
        genes = dbs[cls.gene]
        mappings = dbs[cls.mappings]

        input_dtypes = [
            config['feature_types']['mut'],
            config['feature_types']['cna'],
            config['feature_types']['fusion']
        ]

        samples = Preclinical.generate_sample_list(summary, cls.use_column, cls.model_id)
        idx_feature_type = df_actionable[cls.feature_type].isin(input_dtypes)
        idx_sensitive = ~(df_actionable[cls.sensitive_therapy].isnull() | df_actionable[cls.sensitive_therapy].eq(''))

        dictionary = {}
        for index in df_actionable[idx_feature_type & idx_sensitive].index:
            sensitive_therapy = df_actionable.loc[index, cls.sensitive_therapy]
            if sensitive_therapy not in list(mappings.keys()):
                continue
            mapped = list(mappings[sensitive_therapy][cls.gdsc])
            feature_display = df_actionable.loc[index, cls.feature_display]
            index_dict = {}
            if mapped:
                feature_dictionary = cls.split_samples_by_wt_mut(df_actionable.loc[index, :], dbs, samples, config)
                features = list(feature_dictionary)
                for therapy in mapped:
                    therapy_dict = {}
                    for feature in features:
                        feature_dict_copy = feature_dictionary.copy()
                        wt_samples = feature_dictionary[feature]['samples'][0]
                        mut_samples = feature_dictionary[feature]['samples'][1]
                        drug_dict = cls.create_drug_dict(gdsc, therapy, wt_samples, mut_samples)
                        feature_dict_copy.update({'comparison': drug_dict})
                        therapy_dict.update({feature: feature_dict_copy})
                    figure = PreclinicalEfficacy.draw(therapy_dict, therapy, features)
                    figure_base64 = PreclinicalEfficacy.convert_figure_base64(figure)
                    therapy_dict.update({'figure': figure})
                    therapy_dict.update({'figure_name': f"{feature_display}.{therapy.split(' ')[0]}"})
                    therapy_dict.update({'figure_base64': figure_base64})
                    index_dict.update({therapy: therapy_dict})
                dictionary.update({index: index_dict})
        return dictionary

    @classmethod
    def create_drug_dict(cls, gdsc, therapy, wt_samples, mut_samples):
        wt_ln_ic50 = cls.extract_ln_ic50s(gdsc, cls.therapy, therapy, wt_samples, cls.ln_ic50)
        mut_ln_ic50 = cls.extract_ln_ic50s(gdsc, cls.therapy, therapy, mut_samples, cls.ln_ic50)
        wt_series = cls.calculate_series_exp(wt_ln_ic50)
        mut_series = cls.calculate_series_exp(mut_ln_ic50)
        statistic, p_value = cls.calculate_mann_whitney_u(wt_series, mut_series)
        return {
            cls.n_wt: len(wt_samples),
            cls.n_mut: len(mut_samples),
            cls.n_wt_tested: wt_series.shape[0],
            cls.n_mut_tested: mut_series.shape[0],
            cls.wt_values: wt_series.tolist(),
            cls.mut_values: mut_series.tolist(),
            cls.wt_median: cls.calculate_series_median(wt_series),
            cls.mut_median: cls.calculate_series_median(mut_series),
            cls.wt_mean: cls.calculate_series_mean(wt_series),
            cls.mut_mean: cls.calculate_series_mean(mut_series),
            cls.wt_std: cls.calculate_series_std(wt_series),
            cls.mut_std: cls.calculate_series_std(mut_series),
            cls.statistic: cls.round_value(statistic, cls.rounding_places),
            cls.pvalue: cls.round_scientific_notation(p_value, cls.rounding_places)
        }

    @classmethod
    def extract_ln_ic50s(cls, dataframe, subset_column, subset_value, samples, column):
        return (
            dataframe[dataframe[subset_column].eq(subset_value)]
            .set_index(cls.model_id)
            .reindex(samples)
            .loc[:, column]
            .dropna()
            .astype(float)
        )

    @classmethod
    def filter_empty_values(cls, list_of_strings):
        return [string for string in list_of_strings if string is not cls.fill_na_value]

    @classmethod
    def generate_feature_string(cls, list_of_labels):
        list_of_labels = cls.filter_empty_values(list_of_labels)
        return ' '.join(list_of_labels)

    @classmethod
    def generate_feature_strings(cls, list_of_list_of_labels, series):
        feature_strings = []
        for list_of_labels in list_of_list_of_labels:
            feature_strings.append(cls.generate_feature_string(series.loc[list_of_labels].tolist()))
        return feature_strings

    @classmethod
    def populate_feature_dictionary(cls, groups, all_samples):
        dictionary = {}
        for dataframe, condition, feature_string in groups:
            mutated_samples = cls.retrieve_mut_samples(dataframe, condition)
            wt_samples = cls.retrieve_wt_samples(all_samples, mutated_samples)
            mutated_samples = pd.Series(mutated_samples).dropna().tolist()
            wt_samples = pd.Series(wt_samples).dropna().tolist()
            dictionary[feature_string] = {}
            dictionary[feature_string]['samples'] = [wt_samples, mutated_samples]
        return dictionary

    @classmethod
    def select_split_function(cls, feature_type, variant_string, copy_number_string, fusion_string):
        if feature_type == variant_string:
            return cls.split_samples_for_variants
        elif feature_type == copy_number_string:
            return cls.split_samples_for_copy_numbers
        elif feature_type == fusion_string:
            return cls.split_samples_for_fusions
        else:
            return cls.split_exit

    @classmethod
    def split_exit(cls, dbs, series, samples):
        msg = f'Invalid feature type, {series.loc[cls.feature_type]}, passed when evaluating preclinical efficacy.'
        print(msg)
        exit()

    @classmethod
    def split_samples_by_wt_mut(cls, series, dbs, samples, config):
        feature_type = series.loc[cls.feature_type]
        split_function = cls.select_split_function(
            feature_type=feature_type,
            variant_string=config['feature_types']['mut'],
            copy_number_string=config['feature_types']['cna'],
            fusion_string=config['feature_types']['fusion']
        )
        return split_function(
            dbs=dbs,
            series=series,
            all_samples=samples
        )

    @classmethod
    def split_samples_for_copy_numbers(cls, dbs, series, all_samples):
        genes = dbs[cls.gene]
        db = dbs[cls.cnas]

        feature_cond = genes[cls.feature].eq(series.loc[cls.feature])
        feature_type_cond = db[cls.feature].eq(series.loc[cls.feature])
        alteration_type_cond = feature_type_cond & db[cls.alteration_type].eq(series.loc[cls.alteration_type])

        feature_string_labels = [
            [cls.feature],
            [cls.feature, cls.feature_type],
            [cls.feature, cls.feature_type, cls.alteration_type]
        ]
        feature_strings = cls.generate_feature_strings(feature_string_labels, series)

        groups = [
            (genes, feature_cond, feature_strings[0]),
            (db, feature_type_cond, feature_strings[1]),
            (db, alteration_type_cond, feature_strings[2])
        ]
        return cls.populate_feature_dictionary(groups, all_samples)

    @classmethod
    def split_samples_for_fusions(cls, dbs, series, all_samples):
        genes = dbs[cls.gene]
        db = dbs[cls.fusions]

        gene0 = series.loc[cls.alteration].split('--')[0]
        gene1 = series.loc[cls.alteration].split('--')[1]

        gene0_cond = genes[cls.feature].eq(gene0)
        gene1_cond = genes[cls.feature].eq(gene1)
        gene0_feature_type_cond = db[cls.feature].eq(gene0) | db[cls.partner].eq(gene0)
        gene1_feature_type_cond = db[cls.feature].eq(gene1) | db[cls.partner].eq(gene1)
        fusion_cond = gene0_feature_type_cond & gene1_feature_type_cond

        feature_strings = [
            gene0,
            gene1,
            '{gene} Fusions'.format(gene=gene0),
            '{gene} Fusions'.format(gene=gene1),
            '{gene}--{partner}'.format(gene=gene0, partner=gene1)
        ]

        groups = [
            (genes, gene0_cond, feature_strings[0]),
            (genes, gene1_cond, feature_strings[1]),
            (db, gene0_feature_type_cond, feature_strings[2]),
            (db, gene1_feature_type_cond, feature_strings[3]),
            (db, fusion_cond, feature_strings[4])
        ]
        return cls.populate_feature_dictionary(groups, all_samples)

    @classmethod
    def split_samples_for_variants(cls, dbs, series, all_samples):
        genes = dbs[cls.gene]
        db = dbs[cls.variants]

        feature_cond = genes[cls.feature].eq(series.loc[cls.feature])
        feature_type_cond = db[cls.feature].eq(series.loc[cls.feature])
        alteration_type_cond = feature_type_cond & db[cls.alteration_type].eq(series.loc[cls.alteration_type])
        alteration_cond = alteration_type_cond & db[cls.alteration].eq(series.loc[cls.alteration])

        feature_string_labels = [
            [cls.feature],
            [cls.feature, cls.feature_type],
            [cls.feature, cls.feature_type, cls.alteration_type],
            [cls.feature, cls.feature_type, cls.alteration_type, cls.alteration]
        ]
        feature_strings = cls.generate_feature_strings(feature_string_labels, series.fillna(cls.fill_na_value))

        groups = [
            (genes, feature_cond, feature_strings[0]),
            (db, feature_type_cond, feature_strings[1]),
            (db, alteration_type_cond, feature_strings[2]),
            (db, alteration_cond, feature_strings[3])
        ]
        return cls.populate_feature_dictionary(groups, all_samples)

    @classmethod
    def retrieve_mut_samples(cls, dataframe, condition):
        return dataframe[condition][cls.model_id].drop_duplicates().tolist()

    @staticmethod
    def retrieve_wt_samples(all_samples, mut_samples):
        return pd.Index(all_samples).difference(pd.Index(mut_samples)).tolist()

    @staticmethod
    def round_value(value, places):
        return round(value, places)

    @classmethod
    def round_scientific_notation(cls, value, places):
        if value > 0.0001:
            return cls.round_value(value, places)
        else:
            return ''.join(['{:0.', str(places), 'e}']).format(value)


class SummaryDataFrame(Investigator):
    feature_columns = [Investigator.n_mut, Investigator.n_wt]
    therapy_columns = [Investigator.n_mut_tested, Investigator.n_wt_tested,
                       Investigator.mut_median, Investigator.mut_mean, Investigator.mut_std,
                       Investigator.wt_median, Investigator.wt_mean, Investigator.wt_std,
                       Investigator.pvalue, Investigator.statistic]
    columns = [Investigator.patient_id,
               Investigator.feature_display, Investigator.tested_subfeature,
               Investigator.therapy]
    columns += feature_columns + therapy_columns

    @classmethod
    def create(cls, dictionary, dataframe, patient_id):
        list_of_series = []
        for index in dictionary.keys():
            for therapy in dictionary[index].keys():
                for subfeature in list(dictionary[index][therapy].keys()):
                    if 'figure' in subfeature:
                        continue
                    series = pd.Series(dictionary[index][therapy][subfeature]['comparison'])
                    series.loc[Investigator.feature_display] = dataframe.loc[index, Investigator.feature_display]
                    series.loc[Investigator.tested_subfeature] = subfeature
                    series.loc[Investigator.therapy] = therapy
                    series.loc[Investigator.patient_id] = patient_id
                    list_of_series.append(series.loc[cls.columns])
        if list_of_series:
            return pd.concat(list_of_series, axis=1, ignore_index=True).T
        else:
            return cls.create_empty_dataframe()

    @classmethod
    def create_empty_dataframe(cls):
        return pd.DataFrame(columns=cls.columns)
