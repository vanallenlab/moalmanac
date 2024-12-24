import unittest
import pandas as pd
import math

from datasources import Preclinical
from investigator import Investigator, SensitivityDictionary
from reader import Ini

class UnitTestSensitivityDictionary(unittest.TestCase):
    def test_calculate_series_exp(self):
        series = pd.Series([0, 1, 2, 3, 5, 2], name='value')
        observed = SensitivityDictionary.calculate_series_exp(series)
        expected = pd.Series([1.000, 2.718, 7.389, 20.086, 148.413, 7.389], name='value')
        self.assertEqual(observed.eq(expected).sum(), 6)

    def test_calculate_series_mean(self):
        series = pd.Series([0, 1, 2, 3, 5, 2], name='value')
        self.assertEqual(SensitivityDictionary.calculate_series_mean(series), 2.167)

    def test_calculate_series_median(self):
        series = pd.Series([0, 1, 2, 3, 5, 2], name='value')
        self.assertEqual(SensitivityDictionary.calculate_series_median(series), 2.0)

    def test_calculate_series_std(self):
        series = pd.Series([0, 1, 2, 3, 5, 2], name='value')
        self.assertEqual(SensitivityDictionary.calculate_series_std(series), 1.722)

    def test_calculate_wilcoxon_rank_sum(self):
        series1 = pd.Series([0, 1, 2, 3, 5, 2], name='value')
        series2 = pd.Series([1, 2, 3, 4, 5, 6], name='value')
        series3 = pd.Series([], dtype=float, name='value')
        statistic, pvalue = SensitivityDictionary.calculate_wilcoxon_rank_sum(series1, series2)
        statistic = round(statistic, 5)
        pvalue = round(pvalue, 5)
        self.assertEqual(pvalue, 0.22977)
        self.assertEqual(statistic, -1.20096)
        statistic, pvalue = SensitivityDictionary.calculate_wilcoxon_rank_sum(series1, series3)
        self.assertTrue(math.isnan(pvalue))
        self.assertTrue(math.isnan(statistic))
        statistic, pvalue = SensitivityDictionary.calculate_wilcoxon_rank_sum(series2, series3)
        self.assertTrue(math.isnan(pvalue))
        self.assertTrue(math.isnan(statistic))

    def test_calculate_mann_whitney_u(self):
        series1 = pd.Series([0, 1, 2, 3, 5, 2], name='value')
        series2 = pd.Series([1, 2, 3, 4, 5, 6], name='value')
        series3 = pd.Series([], dtype=float, name='value')
        statistic, pvalue = SensitivityDictionary.calculate_mann_whitney_u(series1, series2)
        self.assertEqual(pvalue, 0.25642920468772246)
        self.assertEqual(statistic, 10.5)
        statistic, pvalue = SensitivityDictionary.calculate_mann_whitney_u(series1, series3)
        self.assertTrue(math.isnan(pvalue))
        self.assertTrue(math.isnan(statistic))
        statistic, pvalue = SensitivityDictionary.calculate_mann_whitney_u(series2, series3)
        self.assertTrue(math.isnan(pvalue))
        self.assertTrue(math.isnan(statistic))

    def test_create(self):
        dbs_paths = {
            'almanac_gdsc_mappings': '../datasources/preclinical/formatted/almanac-gdsc-mappings.json',
            'summary': '../datasources/preclinical/formatted/cell-lines.summary.txt',
            'variants': '../datasources/preclinical/annotated/cell-lines.somatic-variants.annotated.txt',
            'copynumbers': '../datasources/preclinical/annotated/cell-lines.copy-numbers.annotated.txt',
            'fusions': '../datasources/preclinical/annotated/cell-lines.fusions.annotated.txt',
            'fusions1': '../datasources/preclinical/annotated/cell-lines.fusions.annotated.gene1.txt',
            'fusions2': '../datasources/preclinical/annotated/cell-lines.fusions.annotated.gene2.txt',
            'gdsc': '../datasources/preclinical/formatted/sanger.gdsc.txt',
            'dictionary': '../datasources/preclinical/cell-lines.pkl'
        }
        dbs = Preclinical.import_dbs(dbs_paths)
        data_dictionary = {
            'feature_type': ['Somatic Variant'],
            'feature': ['BRAF'],
            'alteration_type': ['Missense'],
            'alteration': ['p.V600E'],
            'feature_display': ['BRAF p.V600E'],
            'sensitive_therapy_name': ['Dabrafenib + Trametinib'],
            'preclinical_efficacy_observed': [1]
        }
        actionable = pd.DataFrame(data_dictionary, index=[0])
        expected_dabrafenib = '2.322e-12'
        expected_trametinib = '2.344e-09'
        config = Ini.read('config.ini', extended_interpolation=False, convert_to_dictionary=False)
        result = SensitivityDictionary.create(dbs, actionable, config)
        self.assertEqual(result[0]['Dabrafenib']['BRAF']['comparison']['pvalue_mww'], expected_dabrafenib)
        self.assertEqual(result[0]['Trametinib']['BRAF']['comparison']['pvalue_mww'], expected_trametinib)

    def test_create_drug_dict(self):
        all_samples = ['A', 'B', 'C', 'D']
        test_dataframe = pd.DataFrame({
            SensitivityDictionary.model_id: all_samples,
            'feature': [0, 1, 2, pd.NA]
        })
        groups = [
            (test_dataframe, test_dataframe['feature'].between(1, 2), 'expects B, C'),
            (test_dataframe, test_dataframe['feature'].isnull(), 'expects D'),
            (test_dataframe, test_dataframe['feature'].eq(0), 'expects A'),
            (test_dataframe, test_dataframe['feature'].ge(5), 'expects None')]
        result = SensitivityDictionary.populate_feature_dictionary(groups, all_samples)
        self.assertEqual(result['expects B, C']['samples'][0], ['A', 'D'])
        self.assertEqual(result['expects B, C']['samples'][1], ['B', 'C'])
        self.assertEqual(result['expects D']['samples'][0], ['A', 'B', 'C'])
        self.assertEqual(result['expects D']['samples'][1], ['D'])
        self.assertEqual(result['expects A']['samples'][0], ['B', 'C', 'D'])
        self.assertEqual(result['expects A']['samples'][1], ['A'])
        self.assertEqual(result['expects None']['samples'][0], ['A', 'B',  'C', 'D'])
        self.assertEqual(result['expects None']['samples'][1], [])

    def test_generate_feature_string(self):
        labels = ['BRAF', 'Somatic Variant', 'Missense', 'p.V600E']
        result = SensitivityDictionary.generate_feature_string(labels)
        self.assertEqual(result, 'BRAF Somatic Variant Missense p.V600E')

    def test_generate_feature_strings(self):
        series = pd.Series({'feature_type': 'Copy Number', 'feature': 'CDKN2A', 'alt_type': 'Deletion'})
        labels = [['feature', 'feature_type', 'alt_type'], ['feature', 'feature_type'], ['feature']]
        result = SensitivityDictionary.generate_feature_strings(labels, series)
        self.assertEqual(result, ['CDKN2A Copy Number Deletion', 'CDKN2A Copy Number', 'CDKN2A'])

    def test_select_split_function(self):
        config = Ini.read('config.ini', extended_interpolation=False, convert_to_dictionary=False)
        var_string = config['feature_types']['mut']
        cn_string = config['feature_types']['cna']
        fusion_string = config['feature_types']['fusion']

        var_function = SensitivityDictionary.split_samples_for_variants
        cn_function = SensitivityDictionary.split_samples_for_copy_numbers
        fusion_function = SensitivityDictionary.split_samples_for_fusions

        self.assertEqual(
            SensitivityDictionary.select_split_function(var_string, var_string, cn_string, fusion_string),
            var_function
        )
        self.assertEqual(
            SensitivityDictionary.select_split_function(cn_string, var_string, cn_string, fusion_string),
            cn_function
        )
        self.assertEqual(
            SensitivityDictionary.select_split_function(fusion_string, var_string, cn_string, fusion_string),
            fusion_function
        )

    def test_split_samples_by_wt_mut(self):
        gene = SensitivityDictionary.gene
        variants = SensitivityDictionary.variants
        cnas = SensitivityDictionary.cnas
        fusions = SensitivityDictionary.fusions

        data = pd.DataFrame({
            SensitivityDictionary.feature_type: ['Somatic Variant', 'Copy Number', 'Rearrangement', 'Aneuploidy'],
            SensitivityDictionary.feature: ['BRAF', 'CDKN2A', 'ERG', 'WGD'],
            SensitivityDictionary.alteration_type: ['Missense', 'Deletion', 'Fusion', ''],
            SensitivityDictionary.alteration: ['p.V600E', '', 'TMPRSS2--ERG', '']
        })

        samples = ['A', 'B', 'C', 'D', 'E']
        db_genes = pd.DataFrame({
            SensitivityDictionary.feature: ['BRAF', 'CDKN2A', 'TMPRSS2', 'ERG', 'TMPRSS2', 'BRAF'],
            SensitivityDictionary.model_id: ['A', 'B', 'C', 'C', 'D', 'E']
        })

        db_variants = pd.DataFrame({
            SensitivityDictionary.feature: ['BRAF', 'BRAF', 'BRAF', 'BRAF', 'TP53'],
            SensitivityDictionary.alteration_type: ['Missense', 'Missense', 'Missense', 'Nonsense', 'Missense'],
            SensitivityDictionary.alteration: ['p.V600E', 'p.V600K', 'p.K500A', 'p.V600*', 'p.I100K'],
            SensitivityDictionary.model_id: samples
        })

        db_cnas = pd.DataFrame({
            SensitivityDictionary.feature: ['CDKN2A', 'CDKN2A', 'TP53'],
            SensitivityDictionary.alteration_type: ['Deletion', 'Amplification', 'Deletion'],
            SensitivityDictionary.model_id: ['B', 'C', 'D']
        })

        db_fusions = pd.DataFrame({
            SensitivityDictionary.feature: ['TMPRSS2', 'TMPRSS2', 'BCR'],
            SensitivityDictionary.partner: ['ERG', 'ABL1', 'ABL1'],
            SensitivityDictionary.model_id: ['C', 'D', 'E']
        })

        dbs = {
            gene: db_genes,
            variants: db_variants,
            cnas: db_cnas,
            fusions: db_fusions
        }

        config = Ini.read('config.ini', extended_interpolation=False, convert_to_dictionary=False)

        results_variants = SensitivityDictionary.split_samples_by_wt_mut(data.loc[0, :], dbs, samples, config)
        results_cnas = SensitivityDictionary.split_samples_by_wt_mut(data.loc[1, :], dbs, samples, config)
        results_fusions = SensitivityDictionary.split_samples_by_wt_mut(data.loc[2, :], dbs, samples, config)
        self.assertEqual(results_variants['BRAF']['samples'][0], ['B', 'C', 'D'])
        self.assertEqual(results_variants['BRAF']['samples'][1], ['A', 'E'])
        self.assertEqual(results_variants['BRAF Somatic Variant']['samples'][0], ['E'])
        self.assertEqual(results_variants['BRAF Somatic Variant']['samples'][1], ['A', 'B', 'C', 'D'])
        self.assertEqual(results_variants['BRAF Somatic Variant Missense']['samples'][0], ['D', 'E'])
        self.assertEqual(results_variants['BRAF Somatic Variant Missense']['samples'][1], ['A', 'B', 'C'])
        self.assertEqual(results_variants['BRAF Somatic Variant Missense p.V600E']['samples'][0], ['B', 'C', 'D', 'E'])
        self.assertEqual(results_variants['BRAF Somatic Variant Missense p.V600E']['samples'][1], ['A'])
        self.assertEqual(results_cnas['CDKN2A']['samples'][0], ['A', 'C', 'D', 'E'])
        self.assertEqual(results_cnas['CDKN2A']['samples'][1], ['B'])
        self.assertEqual(results_cnas['CDKN2A Copy Number']['samples'][0], ['A', 'D', 'E'])
        self.assertEqual(results_cnas['CDKN2A Copy Number']['samples'][1], ['B', 'C'])
        self.assertEqual(results_cnas['CDKN2A Copy Number Deletion']['samples'][0], ['A', 'C', 'D', 'E'])
        self.assertEqual(results_cnas['CDKN2A Copy Number Deletion']['samples'][1], ['B'])
        self.assertEqual(results_fusions['TMPRSS2']['samples'][0], ['A', 'B', 'E'])
        self.assertEqual(results_fusions['TMPRSS2']['samples'][1], ['C', 'D'])
        self.assertEqual(results_fusions['ERG']['samples'][0], ['A', 'B', 'D', 'E'])
        self.assertEqual(results_fusions['ERG']['samples'][1], ['C'])
        self.assertEqual(results_fusions['TMPRSS2 Fusions']['samples'][0], ['A', 'B', 'E'])
        self.assertEqual(results_fusions['TMPRSS2 Fusions']['samples'][1], ['C', 'D'])
        self.assertEqual(results_fusions['ERG Fusions']['samples'][0], ['A', 'B', 'D', 'E'])
        self.assertEqual(results_fusions['ERG Fusions']['samples'][1], ['C'])
        self.assertEqual(results_fusions['TMPRSS2--ERG']['samples'][0], ['A', 'B', 'D', 'E'])
        self.assertEqual(results_fusions['TMPRSS2--ERG']['samples'][1], ['C'])

    def test_split_samples_for_copy_numbers(self):
        gene = SensitivityDictionary.gene
        cnas = SensitivityDictionary.cnas

        data = pd.DataFrame({
            SensitivityDictionary.feature_type: ['Somatic Variant', 'Copy Number', 'Rearrangement', 'Aneuploidy'],
            SensitivityDictionary.feature: ['BRAF', 'CDKN2A', 'ERG', 'WGD'],
            SensitivityDictionary.alteration_type: ['Missense', 'Deletion', 'Fusion', ''],
            SensitivityDictionary.alteration: ['p.V600E', '', 'TMPRSS2--ERG', '']
        })

        samples = ['A', 'B', 'C', 'D', 'E']
        db_genes = pd.DataFrame({
            SensitivityDictionary.feature: ['BRAF', 'CDKN2A', 'TMPRSS2', 'ERG', 'TMPRSS2', 'BRAF'],
            SensitivityDictionary.model_id: ['A', 'B', 'C', 'C', 'D', 'E']
        })

        db_cnas = pd.DataFrame({
            SensitivityDictionary.feature: ['CDKN2A', 'CDKN2A', 'TP53'],
            SensitivityDictionary.alteration_type: ['Deletion', 'Amplification', 'Deletion'],
            SensitivityDictionary.model_id: ['B', 'C', 'D']
        })

        dbs = {
            gene: db_genes,
            cnas: db_cnas,
        }

        config = Ini.read('config.ini', extended_interpolation=False, convert_to_dictionary=False)
        results_cnas = SensitivityDictionary.split_samples_by_wt_mut(data.loc[1, :], dbs, samples, config)
        self.assertEqual(results_cnas['CDKN2A']['samples'][0], ['A', 'C', 'D', 'E'])
        self.assertEqual(results_cnas['CDKN2A']['samples'][1], ['B'])
        self.assertEqual(results_cnas['CDKN2A Copy Number']['samples'][0], ['A', 'D', 'E'])
        self.assertEqual(results_cnas['CDKN2A Copy Number']['samples'][1], ['B', 'C'])
        self.assertEqual(results_cnas['CDKN2A Copy Number Deletion']['samples'][0], ['A', 'C', 'D', 'E'])
        self.assertEqual(results_cnas['CDKN2A Copy Number Deletion']['samples'][1], ['B'])

    def test_split_samples_for_fusions(self):
        gene = SensitivityDictionary.gene
        fusions = SensitivityDictionary.fusions

        data = pd.DataFrame({
            SensitivityDictionary.feature_type: ['Somatic Variant', 'Copy Number', 'Rearrangement', 'Aneuploidy'],
            SensitivityDictionary.feature: ['BRAF', 'CDKN2A', 'ERG', 'WGD'],
            SensitivityDictionary.alteration_type: ['Missense', 'Deletion', 'Fusion', ''],
            SensitivityDictionary.alteration: ['p.V600E', '', 'TMPRSS2--ERG', '']
        })

        samples = ['A', 'B', 'C', 'D', 'E']
        db_genes = pd.DataFrame({
            SensitivityDictionary.feature: ['BRAF', 'CDKN2A', 'TMPRSS2', 'ERG', 'TMPRSS2', 'BRAF'],
            SensitivityDictionary.model_id: ['A', 'B', 'C', 'C', 'D', 'E']
        })

        db_fusions = pd.DataFrame({
            SensitivityDictionary.feature: ['TMPRSS2', 'TMPRSS2', 'BCR'],
            SensitivityDictionary.partner: ['ERG', 'ABL1', 'ABL1'],
            SensitivityDictionary.model_id: ['C', 'D', 'E']
        })

        dbs = {
            gene: db_genes,
            fusions: db_fusions
        }

        config = Ini.read('config.ini', extended_interpolation=False, convert_to_dictionary=False)
        results_fusions = SensitivityDictionary.split_samples_by_wt_mut(data.loc[2, :], dbs, samples, config)
        self.assertEqual(results_fusions['TMPRSS2']['samples'][0], ['A', 'B', 'E'])
        self.assertEqual(results_fusions['TMPRSS2']['samples'][1], ['C', 'D'])
        self.assertEqual(results_fusions['ERG']['samples'][0], ['A', 'B', 'D', 'E'])
        self.assertEqual(results_fusions['ERG']['samples'][1], ['C'])
        self.assertEqual(results_fusions['TMPRSS2 Fusions']['samples'][0], ['A', 'B', 'E'])
        self.assertEqual(results_fusions['TMPRSS2 Fusions']['samples'][1], ['C', 'D'])
        self.assertEqual(results_fusions['ERG Fusions']['samples'][0], ['A', 'B', 'D', 'E'])
        self.assertEqual(results_fusions['ERG Fusions']['samples'][1], ['C'])
        self.assertEqual(results_fusions['TMPRSS2--ERG']['samples'][0], ['A', 'B', 'D', 'E'])
        self.assertEqual(results_fusions['TMPRSS2--ERG']['samples'][1], ['C'])

    def test_split_samples_for_variants(self):
        gene = SensitivityDictionary.gene
        variants = SensitivityDictionary.variants

        data = pd.DataFrame({
            SensitivityDictionary.feature_type: ['Somatic Variant', 'Copy Number', 'Rearrangement', 'Aneuploidy'],
            SensitivityDictionary.feature: ['BRAF', 'CDKN2A', 'ERG', 'WGD'],
            SensitivityDictionary.alteration_type: ['Missense', 'Deletion', 'Fusion', ''],
            SensitivityDictionary.alteration: ['p.V600E', '', 'TMPRSS2--ERG', '']
        })

        samples = ['A', 'B', 'C', 'D', 'E']
        db_genes = pd.DataFrame({
            SensitivityDictionary.feature: ['BRAF', 'CDKN2A', 'TMPRSS2', 'ERG', 'TMPRSS2', 'BRAF'],
            SensitivityDictionary.model_id: ['A', 'B', 'C', 'C', 'D', 'E']
        })

        db_variants = pd.DataFrame({
            SensitivityDictionary.feature: ['BRAF', 'BRAF', 'BRAF', 'BRAF', 'TP53'],
            SensitivityDictionary.alteration_type: ['Missense', 'Missense', 'Missense', 'Nonsense', 'Missense'],
            SensitivityDictionary.alteration: ['p.V600E', 'p.V600K', 'p.K500A', 'p.V600*', 'p.I100K'],
            SensitivityDictionary.model_id: samples
        })

        dbs = {
            gene: db_genes,
            variants: db_variants
        }

        config = Ini.read('config.ini', extended_interpolation=False, convert_to_dictionary=False)
        results_variants = SensitivityDictionary.split_samples_by_wt_mut(data.loc[0, :], dbs, samples, config)
        self.assertEqual(results_variants['BRAF']['samples'][0], ['B', 'C', 'D'])
        self.assertEqual(results_variants['BRAF']['samples'][1], ['A', 'E'])
        self.assertEqual(results_variants['BRAF Somatic Variant']['samples'][0], ['E'])
        self.assertEqual(results_variants['BRAF Somatic Variant']['samples'][1], ['A', 'B', 'C', 'D'])
        self.assertEqual(results_variants['BRAF Somatic Variant Missense']['samples'][0], ['D', 'E'])
        self.assertEqual(results_variants['BRAF Somatic Variant Missense']['samples'][1], ['A', 'B', 'C'])
        self.assertEqual(results_variants['BRAF Somatic Variant Missense p.V600E']['samples'][0], ['B', 'C', 'D', 'E'])
        self.assertEqual(results_variants['BRAF Somatic Variant Missense p.V600E']['samples'][1], ['A'])

    def test_retrieve_mut_samples(self):
        samples = ['A', 'B', 'C', 'D', 'E']
        db = pd.DataFrame({
            SensitivityDictionary.feature: ['BRAF', 'BRAF', 'BRAF', 'BRAF', 'TP53'],
            SensitivityDictionary.alteration_type: ['Missense', 'Missense', 'Missense', 'Nonsense', 'Missense'],
            SensitivityDictionary.alteration: ['p.V600E', 'p.V600K', 'p.K500A', 'p.V600*', 'p.I100K'],
            SensitivityDictionary.model_id: samples
        })
        condition = db[SensitivityDictionary.feature].eq('BRAF')
        self.assertEqual(SensitivityDictionary.retrieve_mut_samples(db, condition), ['A', 'B', 'C', 'D'])

    def test_retrieve_wt_samples(self):
        samples = ['A', 'B', 'C', 'D', 'E']
        self.assertEqual(SensitivityDictionary.retrieve_wt_samples(samples, ['A', 'B']), ['C', 'D', 'E'])
        self.assertEqual(SensitivityDictionary.retrieve_wt_samples(samples, []), samples)

    def test_round_value(self):
        self.assertEqual(SensitivityDictionary.round_value(4.2, 0), 4.0)
        self.assertEqual(SensitivityDictionary.round_value(4.75, 1), 4.8)

    def test_round_scientific_notation(self):
        self.assertEqual(SensitivityDictionary.round_scientific_notation(0.0015, 3), 0.002)
        self.assertEqual(SensitivityDictionary.round_scientific_notation(0.00000123, 2), '1.23e-06')
