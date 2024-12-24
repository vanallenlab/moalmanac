import unittest
import pandas as pd

from features import CoverageMetrics
from moalmanac import features
from config import COLNAMES
from reader import Ini

class UnitTestFeatures(unittest.TestCase):
    def test_annotate_feature_type(self):
        feature_type = 'A'
        idx = [0, 1, 2]
        series = features.Features.annotate_feature_type(feature_type, idx)
        self.assertEqual(['A', 'A', 'A'], series.tolist())
        self.assertEqual(idx, series.index.tolist())

    def test_calculate_percentile(self):
        array = pd.Series([1, 2, 3, 4, 5])
        self.assertEqual(4.96, features.Features.calculate_percentile(array, 99))
        self.assertEqual(1.04, features.Features.calculate_percentile(array, 1))
        self.assertAlmostEqual(1.1, features.Features.calculate_percentile(array, 2.5))
        self.assertAlmostEqual(4.9, features.Features.calculate_percentile(array, 97.5))

    def test_create_empty_dataframe(self):
        self.assertEqual(True, features.Features.create_empty_dataframe().dropna().empty)

    def test_create_empty_series(self):
        self.assertEqual(True, features.Features.create_empty_series().dropna().empty)

    def test_drop_duplicate_genes(self):
        feature = features.Features.feature
        dataframe = pd.DataFrame({feature: ['ABL1', 'BRAF', 'BRAF', 'TP53'], 'sort': [0, 0, 1, 1]})
        self.assertEqual([2, 3, 0], features.Features.drop_duplicate_genes(dataframe, 'sort').tolist())

    def test_import_if_path_exists(self):
        handle = "../example_data/example_patient.capture.somatic.snvs.maf"
        column_map = features.MAFSomatic.create_column_map('tcga_maf_input')
        imported_dataframe = features.Features.import_if_path_exists(handle, '\t', column_map)
        empty_dataframe = features.Features.import_if_path_exists("", '\t', column_map)
        empty_dataframe = features.Features.preallocate_missing_columns(empty_dataframe)
        self.assertEqual(sorted(empty_dataframe.columns.tolist()), sorted(imported_dataframe.columns.tolist()))
        self.assertEqual(299, imported_dataframe.shape[0])
        self.assertEqual(0, empty_dataframe.shape[0])

    def test_preallocate_missing_columns(self):
        columns = features.Features.all_columns
        empty_dataframe = pd.DataFrame()
        partial_dataframe = pd.DataFrame(columns=columns[:4])
        self.assertEqual(columns, features.Features.preallocate_missing_columns(empty_dataframe).columns.tolist())
        self.assertEqual(columns, features.Features.preallocate_missing_columns(partial_dataframe).columns.tolist())


class UnitTestAneuploidy(unittest.TestCase):
    feature = features.Features.feature
    feature_type = features.Features.feature_type


class UnitTestCopyNumber(unittest.TestCase):
    feature_type = features.Features.feature_type

    def test_format_cn_gene(self):
        series = pd.Series(['A', 'A B', 'B  C'])
        self.assertEqual(['A', 'A', 'B'], features.CopyNumber.format_cn_gene(series).tolist())


class UnitTestCopyNumberCalled(unittest.TestCase):
    def test_create_column_map(self):
        column_map = features.CopyNumberCalled.create_column_map()
        self.assertEqual(column_map['gene'], features.Features.feature)
        self.assertEqual(column_map['call'], features.Features.alt_type)

    def test_filter_calls(self):
        amp_string = 'Amplification'
        del_string = 'Deletion'
        tmp = pd.Series(['Amplification', 'Deletion', '', 'Deletion'])
        idx = features.CopyNumberCalled.filter_calls(series=tmp, amp_string=amp_string, del_string=del_string)
        self.assertEqual([0, 1, 3], idx[idx].index.tolist())
        self.assertEqual([2], idx[~idx].index.tolist())


class UnitTestCopyNumberTotal(unittest.TestCase):
    def test_annotate_amp_del(self):
        amp_string = 'Amplification'
        del_string = 'Deletion'
        index = pd.Index([0, 1, 2])
        index_amp = pd.Index([0])
        index_del = pd.Index([2])
        expected = [amp_string, '', del_string]
        result = features.CopyNumberTotal.annotate_amp_del(
            idx=index,
            idx_amp=index_amp,
            idx_del=index_del,
            amp_string=amp_string,
            del_string=del_string
        )
        self.assertEqual(expected, result.tolist())

    def test_create_column_map(self):
        column_map = features.CopyNumberTotal.create_column_map()
        column_names = COLNAMES['seg_input']
        self.assertEqual(column_map[column_names['gene']], features.Features.feature)
        self.assertEqual(column_map[column_names['contig']], features.Features.chr)
        self.assertEqual(column_map[column_names['start']], features.Features.start)
        self.assertEqual(column_map[column_names['end']], features.Features.end)
        self.assertEqual(column_map[column_names['mean']], features.Features.segment_mean)
        self.assertEqual(column_map[column_names['sample']], features.Features.tumor)

    def test_drop_duplicate_genes(self):
        genes = ['CDKN2A', 'CDKN2A', 'BRAF', 'BRAF', '', pd.NA, 'TP53']
        segments = [2, 1, 2, 3, 0, -1, 2]
        expected = pd.Index([0, 3, 4, 5, 6])
        df = pd.DataFrame({features.Features.feature: genes, features.Features.segment_mean: segments})
        for idx in expected:
            self.assertEqual(True, idx in features.CopyNumberTotal.drop_duplicate_genes(df))

    def test_filter_by_threshold(self):
        amp_string = 'Amplification'
        del_string = 'Deletion'
        values = pd.Series(range(1, 101))
        df = pd.DataFrame({features.Features.feature: values,
                           features.Features.chr: values,
                           features.Features.start: values,
                           features.Features.segment_mean: values})
        accept, reject = features.CopyNumberTotal.filter_by_threshold(
            df=df,
            percentile_amp=97.5,
            percentile_del=2.5,
            amp_string=amp_string,
            del_string=del_string
        )
        expected = pd.Index([0, 1, 2, 97, 98, 99])
        for idx in expected:
            self.assertEqual(True, idx in accept.index)

    def test_get_unique_segments(self):
        df = pd.DataFrame({features.Features.chr: [1, 1, 2, 2, 3, 4],
                           features.Features.start: [0, 0, 1, 2, 3, 0],
                           features.Features.segment_mean: [0, 1, 2, 3, 4, 5]})
        expected = [0, 2, 3, 4, 5]
        self.assertEqual(expected, features.CopyNumberTotal.get_unique_segments(df).tolist())


class UnitTestCosmicSignatures(unittest.TestCase):
    def test_round_contributions(self):
        contributions = {
            'COSMIC signature 1': 0.0000,
            'COSMIC signature 3': 0.5001,
            'COSMIC signature 7': 0.9999,
            'COSMIC signature 10': 0.2456,
            'COSMIC signature 11': 0.1985,
        }
        series = pd.Series(contributions)
        result = features.CosmicSignatures.round_contributions(series)
        self.assertEqual(5, result.shape[0])
        self.assertEqual(result.loc['COSMIC signature 1'], 0.000)
        self.assertEqual(result.loc['COSMIC signature 3'], 0.500)
        self.assertEqual(result.loc['COSMIC signature 7'], 1.000)
        self.assertEqual(result.loc['COSMIC signature 10'], 0.246)
        self.assertEqual(result.loc['COSMIC signature 11'], 0.198)

    def test_subset_significant_signatures(self):
        contributions = {
            'COSMIC signature 1': 0,
            'COSMIC signature 3': 0.5,
            'COSMIC signature 7': 1.0,
            'COSMIC signature 10': 0.2,
            'COSMIC signature 11': 0.19,
        }
        series = pd.Series(contributions)
        idx = features.CosmicSignatures.index_for_minimum_contribution(series, minimum_value=0.20)
        result = series[idx]
        self.assertEqual(3, result.shape[0])
        self.assertEqual('COSMIC signature 3', result.index[0])
        self.assertEqual('COSMIC signature 7', result.index[1])
        self.assertEqual('COSMIC signature 10', result.index[2])


class UnitTestCoverageMetrics(unittest.TestCase):
    def test_apply_min_coverage_for_onps(self):
        """
        This function is called by features.CoverageMetrics.format_coverage_col
        after it replaces '__UNKNOWN__' and '' with pd.NA
        """
        series = pd.Series([3, pd.NA, '4|5', '4|', '5'])
        result = CoverageMetrics.apply_min_coverage_for_onps(series)
        self.assertEqual(first=3, second=result.loc[0])
        self.assertEqual(first=4, second=result.loc[2])
        self.assertEqual(first=4, second=result.loc[3])
        self.assertEqual(first='5', second=result.loc[4])
        self.assertTrue(isinstance(result.loc[1], type(pd.NA)))

    def test_calculate_coverage(self):
        """
        This function is called after features.CoverageMetrics.format_coverage_col
        so values should be pd.NA or integers
        """
        series_alt_counts = pd.Series([0, 1, pd.NA, 3])
        series_ref_counts = pd.Series([pd.NA, 0, 1, 2])
        result = (
            CoverageMetrics
            .calculate_coverage(
                series_alt_count=series_alt_counts,
                series_ref_count=series_ref_counts
            )
        )
        self.assertEqual(first=0, second=result.loc[0])
        self.assertEqual(first=1, second=result.loc[1])
        self.assertEqual(first=1, second=result.loc[2])
        self.assertEqual(first=5, second=result.loc[3])

    def test_calculate_tumor_f(self):
        series_alt_count = pd.Series([0, 1, pd.NA, 3])
        series_total_coverage = pd.Series([0, 6, 1, 6])
        result = (
            CoverageMetrics
            .calculate_tumor_f(
                series_alt_count=series_alt_count,
                series_total_coverage=series_total_coverage
            )
        )
        self.assertEqual(first=0.1667, second=result.loc[1])
        self.assertEqual(first=0.5000, second=result.loc[3])
        self.assertTrue(0 not in result.index)
        self.assertTrue(2 not in result.index)

    def test_convert_to_pandas_int64(self):
        series = pd.Series([3, None, pd.NA, '4', '5'])
        result = CoverageMetrics.convert_to_pandas_int64(series)
        self.assertEqual(first=3, second=result.loc[0])
        self.assertEqual(first=4, second=result.loc[3])
        self.assertEqual(first=5, second=result.loc[4])
        self.assertTrue(isinstance(result.loc[1], type(pd.NA)))
        self.assertTrue(isinstance(result.loc[2], type(pd.NA)))

    def test_format_coverage_col(self):
        series = pd.Series([3, '', None, pd.NA, '4|5', '4|', '5', '__UNKNOWN__'])
        result = CoverageMetrics.format_coverage_col(series)
        self.assertEqual(first=3, second=result.loc[0])
        self.assertTrue(isinstance(result.loc[1], type(pd.NA)))
        self.assertTrue(isinstance(result.loc[2], type(pd.NA)))
        self.assertTrue(isinstance(result.loc[3], type(pd.NA)))
        self.assertEqual(first=4, second=result.loc[4])
        self.assertEqual(first=4, second=result.loc[5])
        self.assertEqual(first=5, second=result.loc[6])
        self.assertTrue(isinstance(result.loc[7], type(pd.NA)))

    def test_get_minimum_values_dataframe(self):
        """
        This function is run after format_coverage_col, apply_min_coverage_for_onps, get_onp_boolean, split_counts.
        It is also only run on the index values that contain an ONP boolean.
        """
        series = pd.Series([3, '', None, pd.NA, '4|5', '4|', '5'])
        dataframe = features.CoverageMetrics.split_counts(series)
        result = features.CoverageMetrics.get_minimum_values_dataframe(dataframe, axis=1)
        self.assertEqual(first=3, second=result.loc[0])
        self.assertEqual(first=4, second=result.loc[4])
        self.assertEqual(first=4, second=result.loc[5])
        self.assertEqual(first=5, second=result.loc[6])
        self.assertTrue(0 in result.index)
        self.assertTrue(1 not in result.index)
        self.assertTrue(2 not in result.index)
        self.assertTrue(3 not in result.index)
        self.assertTrue(4 in result.index)
        self.assertTrue(5 in result.index)
        self.assertTrue(6 in result.index)

    def test_get_onp_boolean(self):
        series = pd.Series([3, '', None, pd.NA, '4|5', '4|', '|5'])
        result = features.CoverageMetrics.get_onp_boolean(series)
        self.assertEqual(first=False, second=result.loc[0])
        self.assertEqual(first=False, second=result.loc[1])
        self.assertEqual(first=False, second=result.loc[2])
        self.assertEqual(first=False, second=result.loc[3])
        self.assertEqual(first=True, second=result.loc[4])
        self.assertEqual(first=True, second=result.loc[5])
        self.assertEqual(first=True, second=result.loc[6])

    def test_safe_cast(self):
        series = pd.Series([3, '', None, pd.NA, '4|5', '4|', '5', '__UNKNOWN__'])
        result = features.CoverageMetrics.safe_cast(series.loc[0])
        self.assertEqual(first=3, second=result)
        result = features.CoverageMetrics.safe_cast(series.loc[6])
        self.assertEqual(first=5, second=result)
        result = features.CoverageMetrics.safe_cast(series.loc[1])
        self.assertTrue(isinstance(result, type(pd.NA)))
        result = features.CoverageMetrics.safe_cast(series.loc[2])
        self.assertTrue(isinstance(result, type(pd.NA)))
        result = features.CoverageMetrics.safe_cast(series.loc[3])
        self.assertTrue(isinstance(result, type(pd.NA)))
        result = features.CoverageMetrics.safe_cast(series.loc[4])
        self.assertTrue(isinstance(result, type(pd.NA)))
        result = features.CoverageMetrics.safe_cast(series.loc[5])
        self.assertTrue(isinstance(result, type(pd.NA)))
        result = features.CoverageMetrics.safe_cast(series.loc[7])
        self.assertTrue(isinstance(result, type(pd.NA)))

    def test_split_counts(self):
        series = pd.Series([3, '', None, pd.NA, '4|5'])
        result = features.CoverageMetrics.split_counts(series)
        self.assertEqual(first=series.shape[0], second=result.shape[0])
        self.assertEqual(first=2, second=result.shape[1])
        self.assertEqual(first='3', second=result.loc[0, 0])
        self.assertEqual(first='4', second=result.loc[4, 0])
        self.assertEqual(first='5', second=result.loc[4, 1])

        self.assertEqual(first='', second=result.loc[1, 0])
        self.assertEqual(first='', second=result.loc[2, 0])
        self.assertEqual(first='', second=result.loc[3, 0])
        self.assertEqual(first=None, second=result.loc[0, 1])
        self.assertEqual(first=None, second=result.loc[1, 1])
        self.assertEqual(first=None, second=result.loc[2, 1])
        self.assertEqual(first=None, second=result.loc[3, 1])


class UnitTestFusion(unittest.TestCase):
    def test_create_column_map(self):
        config = Ini.read(path='config.ini', extended_interpolation=False, convert_to_dictionary=False)
        column_map = features.Fusion.create_colmap(config)
        leftbreakpoint = 'leftbreakpoint'
        rightbreakpoint = 'rightbreakpoint'

        values = list(column_map.values())
        self.assertEqual(4, len(column_map))
        self.assertEqual(features.Features.feature, values[0])
        self.assertEqual(features.Features.spanningfrags, values[1])
        self.assertEqual(leftbreakpoint, values[2])
        self.assertEqual(rightbreakpoint, values[3])

    def test_filter_by_spanning_fragment_count(self):
        series = pd.Series([4, 5, 6])
        result = features.Fusion.filter_by_spanning_fragment_count(series, 5.0)
        self.assertEqual(False, 0 in result.tolist())
        self.assertEqual(True, 1 in result.tolist())
        self.assertEqual(True, 2 in result.tolist())

    def test_split_genes(self):
        series = pd.Series(['foo--bar', 'foo-bar', 'A--B'])
        result = features.Fusion.split_genes(series)
        self.assertEqual('foo', result.loc[0, features.Features.left_gene])
        self.assertEqual('bar', result.loc[0, features.Features.right_gene])
        self.assertEqual('foo-bar', result.loc[1, features.Features.left_gene])
        self.assertEqual(None, result.loc[1, features.Features.right_gene])
        self.assertEqual('A', result.loc[2, features.Features.left_gene])
        self.assertEqual('B', result.loc[2, features.Features.right_gene])

    def test_split_breakpoint(self):
        series = pd.Series(['9:100:+', 'chr9:100', ''])
        result = features.Fusion.split_breakpoint(series)
        self.assertEqual('9', result.loc[0, features.Features.chr])
        self.assertEqual('100', result.loc[0, features.Features.start])
        self.assertEqual('chr9', result.loc[1, features.Features.chr])
        self.assertEqual('100', result.loc[1, features.Features.start])
        self.assertEqual('', result.loc[2, features.Features.chr])
        self.assertEqual(None, result.loc[2, features.Features.start])


class UnitTestMAF(unittest.TestCase):
    def test_check_format(self):
        df = pd.DataFrame(columns=['Protein_Change', ''])
        self.assertEqual("tcga_maf_input", features.MAF.check_format(df))
        df = pd.DataFrame(columns=['HGVSp_Short', ''])
        self.assertEqual("gdc_maf_input", features.MAF.check_format(df))
        df = pd.DataFrame(columns=[])
        with self.assertRaises(SystemExit) as system_exit:
            features.MAF.check_format(df)
        code = "Neither 'Protein_Change' nor 'HGVSp_Short' are present columns, cannot map to MAF format"
        self.assertEqual(system_exit.exception.code, code)

    def test_create_column_map(self):
        tcga_column_map = features.MAF.create_column_map('tcga_maf_input')
        gdc_column_map = features.MAF.create_column_map('gdc_maf_input')
        tcga_expected_key = "protein_change" in list(tcga_column_map.keys())
        gdc_expected_key = "hgvsp_short" in list(gdc_column_map.keys())
        tcga_unexpected_key = "hgvsp_short" in list(tcga_column_map.keys())
        gdc_unexpected_key = "protein_change" in list(gdc_column_map.keys())
        self.assertEqual(True, tcga_expected_key)
        self.assertEqual(True, gdc_expected_key)
        self.assertEqual(False, tcga_unexpected_key)
        self.assertEqual(False, gdc_unexpected_key)

    def test_format_maf(self):
        data = {
            features.MAF.alt_type: ['Missense_Mutation', 'Nonsense_Mutation', 'A'],
            features.MAF.alt_count: [0, 0, 0],
            features.MAF.ref_count: [1, 1, 1]
        }
        initial = pd.DataFrame(data)
        result = features.MAF.format_maf(initial, 'Somatic Variant')
        self.assertEqual('Missense', result.loc[0, features.MAF.alt_type])
        self.assertEqual('Nonsense', result.loc[1, features.MAF.alt_type])
        self.assertEqual('A', result.loc[2, features.MAF.alt_type])
        self.assertEqual('Somatic Variant', result[features.MAF.feature_type].unique().tolist()[0])
        self.assertEqual(1, result[features.MAF.feature_type].drop_duplicates().shape[0])

    def test_rename_coding_classifications(self):
        alt_types = ['Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site',
                     'Frame_Shift_Ins', 'Frame_Shift_Del', 'In_Frame_Ins', 'In_Frame_Del',
                     '', 'A', 5, pd.NA]
        initial = pd.Series(alt_types)
        result = features.MAF.rename_coding_classifications(initial)
        self.assertEqual('Missense', result.loc[0])
        self.assertEqual('Nonsense', result.loc[1])
        self.assertEqual('Nonstop', result.loc[2])
        self.assertEqual('Splice Site', result.loc[3])
        self.assertEqual('Frameshift', result.loc[4])
        self.assertEqual('Frameshift', result.loc[5])
        self.assertEqual('Insertion', result.loc[6])
        self.assertEqual('Deletion', result.loc[7])
        self.assertEqual('', result.loc[8])
        self.assertEqual('A', result.loc[9])
        self.assertEqual('5', result.loc[10])

    def test_return_idx_variants_coding(self):
        initial = pd.Series(['Missense', 'Silent', 'A', 1, pd.NA])
        result = features.MAF.return_idx_variants_coding(initial)
        self.assertEqual(True, result.loc[0])
        for index in [1, 2, 3]:
            self.assertEqual(False, result.loc[index])

    def test_return_variants_coding(self):
        initial = pd.DataFrame({features.MAF.alt_type: ['Missense', 'Silent', 'A', 1, pd.NA]})
        result = features.MAF.return_variants_coding(initial)
        self.assertEqual([0], result.index.tolist())

    def test_return_variants_non_coding(self):
        initial = pd.DataFrame({features.MAF.alt_type: ['Missense', 'Silent', 'A', 1, pd.NA]})
        result = features.MAF.return_variants_non_coding(initial)
        self.assertEqual([1, 2, 3, 4], result.index.tolist())
