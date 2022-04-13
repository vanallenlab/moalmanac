import unittest
import pandas as pd

from moalmanac import features
from config import CONFIG, COLNAMES


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
        tmp = pd.Series(['Amplification', 'Deletion', '', 'Deletion'])
        idx = features.CopyNumberCalled.filter_calls(tmp)
        self.assertEqual([0, 1, 3], idx[idx].index.tolist())
        self.assertEqual([2], idx[~idx].index.tolist())


class UnitTestCopyNumberTotal(unittest.TestCase):
    def test_annotate_amp_del(self):
        index = pd.Index([0, 1, 2])
        index_amp = pd.Index([0])
        index_del = pd.Index([2])
        expected = [features.CopyNumberTotal.amplification, '', features.CopyNumberTotal.deletion]
        self.assertEqual(expected, features.CopyNumberTotal.annotate_amp_del(index, index_amp, index_del).tolist())

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
        values = pd.Series(range(1, 101))
        df = pd.DataFrame({features.Features.feature: values,
                           features.Features.chr: values,
                           features.Features.start: values,
                           features.Features.segment_mean: values})
        accept, reject = features.CopyNumberTotal.filter_by_threshold(df, 97.5, 2.5)
        expected = pd.Index([0, 1, 2, 97, 98, 99])
        for idx in expected:
            self.assertEqual(True, idx in accept.index)

    def test_get_unique_segments(self):
        df = pd.DataFrame({features.Features.chr: [1, 1, 2, 2, 3, 4],
                           features.Features.start: [0, 0, 1, 2, 3, 0],
                           features.Features.segment_mean: [0, 1, 2, 3, 4, 5]})
        expected = [0, 2, 3, 4, 5]
        self.assertEqual(expected, features.CopyNumberTotal.get_unique_segments(df).tolist())


class UnitTestFusion(unittest.TestCase):
    def test_create_column_map(self):
        column_map = features.Fusion.create_colmap()
        values = list(column_map.values())
        self.assertEqual(4, len(column_map))
        self.assertEqual(features.Features.feature, values[0])
        self.assertEqual(features.Features.spanningfrags, values[1])
        self.assertEqual(features.Fusion.leftbreakpoint, values[2])
        self.assertEqual(features.Fusion.rightbreakpoint, values[3])

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
        initial = pd.DataFrame({features.MAF.alt_type: ['Missense_Mutation', 'Nonsense_Mutation', 'A']})
        result = features.MAF.format_maf(initial, 'Somatic Variant')
        self.assertEqual(249, result.shape[1])
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
