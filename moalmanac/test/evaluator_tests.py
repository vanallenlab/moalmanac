import unittest
import numpy as np
import pandas as pd

from evaluator import Evaluator, Actionable, Integrative, Microsatellite, Strategies


class UnitTestEvaluator(unittest.TestCase):
    def test_assign_bin(self):
        bin_name = 'test_bin'
        bin_label = 'test_label'
        score_bin = Evaluator.score_bin
        df = pd.DataFrame({score_bin: ['Almanac', 'Cancer Hotspot', np.nan],
                           bin_name: [0, 1, 1]})
        series = Evaluator.assign_bin(df, bin_name, bin_label)
        self.assertEqual(['Almanac', 'Cancer Hotspot', 'test_label'], series.tolist())

    def test_assign_vus(self):
        score_bin = Evaluator.score_bin
        df = pd.DataFrame({score_bin: ['Almanac', 'Cancer Hotspot', np.nan]})
        series = Evaluator.assign_vus(df)
        self.assertEqual(['Almanac', 'Cancer Hotspot', 'VUS'], series.tolist())

    def test_evaluate_almanac(self):
        score_bin = Evaluator.score_bin
        almanac_bin = Evaluator.almanac_bin
        sensitive_bin = Evaluator.sensitive_bin
        resistance_bin = Evaluator.resistance_bin
        prognostic_bin = Evaluator.prognostic_bin

        df = pd.DataFrame({
            almanac_bin: [0, 1, 2, 3, 4],
            sensitive_bin: [0, 1, 2, 3, 4],
            resistance_bin: [0, 1, 2, 3, 4],
            prognostic_bin: [0, 1, 2, 3, 4]})
        df[score_bin] = ''
        annotated = Evaluator.evaluate_almanac(df)
        annotated.replace(0, '', inplace=True)
        for column in [score_bin, sensitive_bin, resistance_bin, prognostic_bin]:
            self.assertEqual([np.nan,
                              'Biologically Relevant',
                              'Investigate Actionability',
                              'Investigate Actionability',
                              'Putatively Actionable'],
                             annotated.loc[:, column].tolist())

    def test_evaluate_somatic(self):
        score_bin = Evaluator.score_bin
        almanac_bin = Evaluator.almanac_bin
        sensitive_bin = Evaluator.sensitive_bin
        resistance_bin = Evaluator.resistance_bin
        prognostic_bin = Evaluator.prognostic_bin
        cancerhotspots_bin = Evaluator.cancerhotspots_bin
        cancerhotspots3d_bin = Evaluator.cancerhotspots3d_bin
        cgc_bin = Evaluator.cgc_bin
        gsea_pathways_bin = Evaluator.gsea_pathways_bin
        gsea_modules_bin = Evaluator.gsea_modules_bin
        cosmic_bin = Evaluator.cosmic_bin

        df = pd.DataFrame({
            almanac_bin: [1, 0, 0, 0, 0, 0, 0, 0],
            sensitive_bin: [0, 0, 0, 0, 0, 0, 0, 0],
            resistance_bin: [0, 0, 0, 0, 0, 0, 0, 0],
            prognostic_bin: [0, 0, 0, 0, 0, 0, 0, 0],
            cancerhotspots_bin: [0, 1, 0, 0, 0, 0, 0, 0],
            cancerhotspots3d_bin: [0, 0, 1, 0, 0, 0, 0, 0],
            cgc_bin: [0, 0, 0, 1, 0, 0, 0, 0],
            gsea_pathways_bin: [0, 0, 0, 0, 1, 0, 0, 0],
            gsea_modules_bin: [0, 0, 0, 0, 0, 1, 0, 0],
            cosmic_bin: [0, 0, 0, 0, 0, 0, 1, 0]})
        annotated = Evaluator.evaluate_somatic(df)
        self.assertEqual(['Biologically Relevant', 'Cancer Hotspot', 'Cancer Hotspot 3D',
                          'Cancer Gene Census', 'Cancer Pathway', 'Cancer Module', 'Cosmic', 'VUS'],
                         annotated[score_bin].tolist())

    def test_evaluate_germline(self):
        score_bin = Evaluator.score_bin
        almanac_bin = Evaluator.almanac_bin
        sensitive_bin = Evaluator.sensitive_bin
        resistance_bin = Evaluator.resistance_bin
        prognostic_bin = Evaluator.prognostic_bin
        cancerhotspots_bin = Evaluator.cancerhotspots_bin
        cgc_bin = Evaluator.cgc_bin

        df = pd.DataFrame({
            almanac_bin: [1, 0, 0],
            sensitive_bin: [0, 0, 0],
            resistance_bin: [0, 0, 0],
            prognostic_bin: [0, 0, 0],
            cancerhotspots_bin: [0, 1, 0],
            cgc_bin: [0, 0, 1]})
        annotated = Evaluator.evaluate_germline(df)
        self.assertEqual(['Biologically Relevant', 'Cancer Hotspot', 'Cancer Gene Census'],
                         annotated[score_bin].tolist())

    def test_map_almanac_bins(self):
        series = pd.Series([0, 1, 2, 3, 4])
        annotated = Evaluator.map_almanac_bins(series)
        expected = [np.nan, 'Biologically Relevant', 'Investigate Actionability',
                    'Investigate Actionability', 'Putatively Actionable']
        self.assertEqual(expected, annotated.tolist())

    def test_remap_almanac_bins(self):
        series = pd.Series([0, 0, 1, 2])
        remapped = Evaluator.remap_almanac_bins(series=series, old_value=0, new_value=4)
        self.assertEqual(4, remapped.loc[0])
        self.assertEqual(4, remapped.loc[1])
        self.assertEqual(1, remapped.loc[2])
        self.assertEqual(2, remapped.loc[3])

    def test_remove_low_allele_fraction_variants(self):
        feature_type = Evaluator.feature_type
        mut_type = 'Somatic Variant'
        germline_type = 'Germline Variant'
        tumor_f = Evaluator.tumor_f
        min_af = 0.05

        low_af = float(min_af) - 0.01
        high_af = float(min_af) + 0.01

        df = pd.DataFrame({feature_type: [mut_type, mut_type, germline_type, germline_type, 'Aneuploidy'],
                           tumor_f: [low_af, high_af, low_af, high_af, np.nan]})
        subsetted = Evaluator.remove_low_allele_fraction_variants(df, minimum_allele_fraction=min_af)
        self.assertEqual([1, 3, 4], subsetted.index.tolist())

    def test_remove_low_coverage_variants(self):
        feature_type = Evaluator.feature_type
        mut_type = 'Somatic Variant'
        germline_type = 'Germline Variant'
        coverage = Evaluator.coverage
        min_coverage = 15

        low_coverage = float(min_coverage) - 10
        high_coverage = float(min_coverage) + 10

        df = pd.DataFrame({feature_type: [mut_type, mut_type, germline_type, germline_type, 'Aneuploidy'],
                           coverage: [low_coverage, high_coverage, low_coverage, high_coverage, np.nan]})
        subsetted = Evaluator.remove_low_coverage_variants(df, minimum_coverage=min_coverage)
        self.assertEqual([1, 3, 4], subsetted.index.tolist())

    def test_remove_benign_variants(self):
        clinvar_bin = Evaluator.clinvar_bin
        df = pd.DataFrame({clinvar_bin: [np.nan, 0, 1, 0],
                           'label': ['a', 'b', 'c', 'd']})
        subsetted = Evaluator.remove_benign_variants(df)
        self.assertEqual(['a', 'c'], subsetted['label'].tolist())

    def test_remove_common_variants(self):
        exac_common_bin = Evaluator.exac_common_bin
        df = pd.DataFrame({exac_common_bin: [np.nan, 0, 1, 0],
                           'label': ['a', 'b', 'c', 'd']})
        subsetted = Evaluator.remove_common_variants(df)
        self.assertEqual(['a', 'b', 'd'], subsetted['label'].tolist())

    def test_subset_almanac_bin(self):
        almanac_bin = Evaluator.almanac_bin
        df = pd.DataFrame({almanac_bin: [np.nan, 0, 1, 0],
                           'label': ['a', 'b', 'c', 'd']})
        subsetted = Evaluator.subset_almanac_bin(df)
        self.assertEqual(['c'], subsetted['label'].tolist())


class UnitTestActionable(unittest.TestCase):
    def test_create_string_list(self):
        series = pd.Series(['the', 'the', 'quick', 'fox'])
        self.assertEqual('the, quick, fox', Actionable.create_string_list(series))

    def test_display_aneuploidy(self):
        feature = Evaluator.feature
        df = pd.DataFrame({feature: ['A', 'B', 'C']})
        idx = [0, 2]
        series = Actionable.display_aneuploidy(df, idx)
        self.assertEqual(['A', 'C'], series.tolist())

    def test_display_burden(self):
        alt = Evaluator.alt
        df = pd.DataFrame({alt: ["10 mutations per Mb", "20", "30 mutations per Mb"]})
        idx = [0, 2]
        series = Actionable.display_burden(df, idx)
        self.assertEqual(['10 mutations per Mb', '30 mutations per Mb'], series.tolist())

    def test_display_copynumber(self):
        feature = Evaluator.feature
        alt_type = Evaluator.alt_type
        df = pd.DataFrame({feature: ['Foo', 'Bar', 'FooBar'],
                           alt_type: ['Amp', 'Amp', 'Del']})
        idx = [0, 2]
        series = Actionable.display_copynumber(df, idx)
        self.assertEqual(['Foo Amp', 'FooBar Del'], series.tolist())

    def test_display_fusion(self):
        alt = Evaluator.alt
        df = pd.DataFrame({alt: ['Foo--Bar', 'Bar--Foo', 'FooBar--Alpha']})
        idx = [0, 2]
        series = Actionable.display_fusion(df, idx)
        self.assertEqual(['Foo--Bar Fusion', 'FooBar--Alpha Fusion'], series.tolist())

    def test_display_microsatellite_stability(self):
        feature = Evaluator.feature
        df = pd.DataFrame({feature: ['A', 'B', 'C']})
        idx = [0, 2]
        series = Actionable.display_microsatellite_stability(df, idx)
        self.assertEqual(['A', 'C'], series.tolist())

    def test_display_microsatellite_variants(self):
        feature = Evaluator.feature
        alt = Evaluator.alt
        df = pd.DataFrame({feature: ['Foo', '', 'Bar'],
                           alt: ['Amp', '', 'Del']})
        idx = [0, 2]
        series = Actionable.display_microsatellite_variants(df, idx)
        self.assertEqual(['Foo: Amp', 'Bar: Del'], series.tolist())

    def test_display_signature(self):
        feature = Evaluator.feature
        alt = Evaluator.alt
        df = pd.DataFrame({feature: ['Signature 1', '', 'Signature 2'], alt: [0.523, '', 0.0145]})
        idx = [0, 2]
        series = Actionable.display_signature(df, idx)
        self.assertEqual(['Signature 1 (52%)', 'Signature 2 (1%)'],
                         series.tolist())

    def test_display_variant(self):
        feature = Evaluator.feature
        alt_type = Evaluator.alt_type
        alt = Evaluator.alt
        df = pd.DataFrame({feature: ['Foo', 'Bar', 'FooBar'],
                           alt_type: ['Missense', 'Nonsense', 'Frameshift'],
                           alt: ['p.V600E', 'p.N500*', 'p.L151fs*']})
        idx = [0, 2]
        series = Actionable.display_variant(df, idx)
        self.assertEqual(['Foo p.V600E (Missense)', 'FooBar p.L151fs* (Frameshift)'], series.tolist())

    def test_format_variant_classification(self):
        series = pd.Series(['Missense_Mutation', 'Frameshift', '', 'Nonsense_Mutation'])
        replaced = Actionable.format_variant_classification(series).tolist()
        self.assertEqual(['Missense', 'Frameshift', '', 'Nonsense'], replaced)


class UnitTestIntegrative(unittest.TestCase):
    def test_create_integrated_df(self):
        genes = ['BRAF', 'TP53', 'CDKN2A']
        dataframe = Integrative.create_integrated_df(genes)
        self.assertEqual(genes, dataframe.index.tolist())
        self.assertEqual(len(Integrative.columns), dataframe.columns.shape[0])

    def test_extract_feature_type(self):
        feature_type = Integrative.feature_type
        df = pd.DataFrame({feature_type: ['A', 'B', 'C', 'A']})
        self.assertEqual([0, 3], Integrative.extract_feature_type(df, 'A').index.tolist())

    def test_join_alteration_types(self):
        column1 = "column1"
        column2 = "column2"
        dataframe = pd.DataFrame({column1: ['A', np.nan, 'C'], column2: ['D', 'E', np.nan]})
        joined = Integrative.join_alteration_types(dataframe, [column1, column2])
        self.assertEqual('A D,  E, C ', joined)

    def test_join_alterations(self):
        dataframe = pd.DataFrame({'A': ['A', 'B', 'C']})
        self.assertEqual('A, B, C', Integrative.join_alterations(dataframe, 'A'))

    def test_return_genes(self):
        feature = Integrative.feature
        dataframe = pd.DataFrame({feature: [np.nan, 'A', 'B', 'A']})
        self.assertEqual(['A', 'B'], Integrative.return_genes(dataframe))

    def test_subset_nonempty_df(self):
        dataframe = pd.DataFrame({'A': [False, False, True], 'B': [True, False, False]})
        self.assertEqual([0, 2], Integrative.subset_nonempty_df(dataframe).index.tolist())


class UnitTestMicrosatellite(unittest.TestCase):
    def test_create_integrated_df(self):
        msi_bin = Evaluator.msi_bin
        alt_type = Evaluator.alt_type
        dataframe = pd.DataFrame({msi_bin: [np.nan, 0, 1, 1], alt_type: ['Missense', 'Nonstop', 'Missense', 'Nonstop']})
        self.assertEqual([3], Microsatellite.return_msi_variants(dataframe).index.tolist())


class UnitTestStrategies(unittest.TestCase):
    def test_get_union_strategies(self):
        list_1 = ['a', 'b']
        list_2 = ['b', 'c']
        union = Strategies.get_union_strategies(list_1, list_2)
        self.assertEqual(['a', 'b', 'c'], union)

    def test_list_to_string(self):
        list_1 = ['a', 'b', 'c']
        expected_result = 'a, b, c'
        returned_result = Strategies.list_to_string(list_1, ', ')
        self.assertEqual(expected_result, returned_result)

    def test_report_therapy_strategies(self):
        columns = [Strategies.sensitive_therapy_name, Strategies.sensitive_therapy_strategy,
                   Strategies.resistance_therapy_name, Strategies.resistance_therapy_strategy]

        df = pd.DataFrame(pd.NA, columns=columns, index=[0, 1, 2])
        df.loc[0, :] = ['AZD3759', 'EGFR inhibition', 'Afatinib', 'EGFR inhibition']
        df.loc[1, :] = [pd.NA, pd.NA, 'Afatinib', 'EGFR inhibition']
        df.loc[2, :] = ['Cetuximab', 'EGFR inhibition', pd.NA, pd.NA]

        expected_result_sensitive = 'AZD3759, Cetuximab'
        expected_result_resistance = 'Afatinib'

        returned_result = Strategies.report_therapy_strategies(df)
        self.assertEqual(returned_result.columns.tolist(), ['EGFR inhibition'])
        self.assertEqual(returned_result.loc[Strategies.sensitivity, 'EGFR inhibition'], expected_result_sensitive)
        self.assertEqual(returned_result.loc[Strategies.resistance, 'EGFR inhibition'], expected_result_resistance)

    def test_series_to_list(self):
        series = pd.Series(['a', 'b', pd.NA, 'c'])
        self.assertEqual(Strategies.series_to_list(series), ['a', 'b', 'c'])
