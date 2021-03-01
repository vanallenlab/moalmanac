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


class UnitTestCalledCNAReader(unittest.TestCase):
    feature = features.Features.feature
    feature_type = features.Features.feature_type

    def test_create_column_map(self):
        column_map = features.CalledCNAReader.create_column_map(COLNAMES['called_cn_input'])
        self.assertEqual(column_map['gene'], features.Features.feature)
        self.assertEqual(column_map['call'], features.Features.alt_type)

    def test_filter_calls(self):
        tmp = pd.Series(['Amplification', 'Deletion', '', 'Deletion'])
        idx = features.CalledCNAReader.filter_calls(tmp)
        self.assertEqual([0, 1, 3], idx[idx].index.tolist())
        self.assertEqual([2], idx[~idx].index.tolist())


class UnitTestCNVReader(unittest.TestCase):
    def test_annotate_amp_del(self):
        index = pd.Index([0, 1, 2])
        index_amp = pd.Index([0])
        index_del = pd.Index([2])
        expected = [features.CNVReader.amplification, '', features.CNVReader.deletion]
        self.assertEqual(expected, features.CNVReader.annotate_amp_del(index, index_amp, index_del).tolist())

    def test_create_column_map(self):
        column_map = features.CNVReader.create_column_map(COLNAMES['seg_input'])
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
            self.assertEqual(True, idx in features.CNVReader.drop_duplicate_genes(df))

    def test_filter_by_threshold(self):
        values = pd.Series(range(1, 101))
        df = pd.DataFrame({features.Features.feature: values,
                           features.Features.chr: values,
                           features.Features.start: values,
                           features.Features.segment_mean: values})
        accept, reject = features.CNVReader.filter_by_threshold(df, 97.5, 2.5)
        expected = pd.Index([0, 1, 2, 97, 98, 99])
        for idx in expected:
            self.assertEqual(True, idx in accept.index)

    def test_get_unique_segments(self):
        df = pd.DataFrame({features.Features.chr: [1, 1, 2, 2, 3, 4],
                           features.Features.start: [0, 0, 1, 2, 3, 0],
                           features.Features.segment_mean: [0, 1, 2, 3, 4, 5]})
        expected = [0, 2, 3, 4, 5]
        self.assertEqual(expected, features.CNVReader.get_unique_segments(df).tolist())
