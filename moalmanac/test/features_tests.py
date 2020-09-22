import unittest
import pandas as pd

from features import Features


class UnitTestFeatures(unittest.TestCase):
    def test_annotate_feature_type(self):
        feature_type = 'A'
        idx = [0, 1, 2]
        series = Features.annotate_feature_type(feature_type, idx)
        self.assertEqual(['A', 'A', 'A'], series.tolist())
        self.assertEqual(idx, series.index.tolist())

    def test_calculate_percentile(self):
        array = pd.Series([1, 2, 3, 4, 5])
        self.assertEqual(4.96, Features.calculate_percentile(array, 99))
        self.assertEqual(1.04, Features.calculate_percentile(array, 1))
        self.assertAlmostEqual(1.1, Features.calculate_percentile(array, 2.5))
        self.assertAlmostEqual(4.9, Features.calculate_percentile(array, 97.5))

    def test_create_empty_dataframe(self):
        self.assertEqual(True, Features.create_empty_dataframe().dropna().empty)

    def test_create_empty_series(self):
        self.assertEqual(True, Features.create_empty_series().dropna().empty)

    def test_drop_duplicate_genes(self):
        feature = Features.feature
        dataframe = pd.DataFrame({feature: ['ABL1', 'BRAF', 'BRAF', 'TP53'], 'sort': [0, 0, 1, 1]})
        self.assertEqual([2, 3, 0], Features.drop_duplicate_genes(dataframe, 'sort').tolist())

    def test_preallocate_missing_columns(self):
        columns = Features.all_columns
        empty_dataframe = pd.DataFrame()
        partial_dataframe = pd.DataFrame(columns=columns[:4])
        self.assertEqual(columns, Features.preallocate_missing_columns(empty_dataframe).columns.tolist())
        self.assertEqual(columns, Features.preallocate_missing_columns(partial_dataframe).columns.tolist())


class UnitTestAneuploidy(unittest.TestCase):
    feature = Features.feature
    feature_type = Features.feature_type




