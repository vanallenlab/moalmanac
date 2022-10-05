import unittest
import pandas as pd

from writer import Writer, GermlineCancer


class UnitTestWriter(unittest.TestCase):
    def test_create_output_name(self):
        string1 = '.'
        string2 = 'Foo'
        string3 = 'Bar'
        self.assertEqual('./Foo.Bar', Writer.create_output_name(string1, string2, string3))

    def test_sort_columns(self):
        df = pd.DataFrame({'A': [0, 1, 2], 'B': [1, 1, 2]})
        self.assertEqual([0, 1, 2], list(Writer.sort_columns(df, ['A', 'B'], [True, True]).index))
        self.assertEqual([2, 1, 0], list(Writer.sort_columns(df, ['A', 'B'], [False, True]).index))
        self.assertEqual([1, 0, 2], list(Writer.sort_columns(df, ['B', 'A'], [True, False]).index))
        self.assertEqual([2, 0, 1], list(Writer.sort_columns(df, ['B', 'A'], [False, True]).index))

    def test_return_nonzero_idx(self):
        series = pd.Series([0, 1, 2, 0])
        self.assertEqual([1, 2], list(Writer.return_nonzero_bin_idx(series)))


class UnitTestGermlineCancer(unittest.TestCase):
    def test_get_cancer_idx(self):
        almanac_bin = GermlineCancer.almanac_bin
        hotspots_bin = GermlineCancer.hotspots_bin
        cgc_bin = GermlineCancer.cgc_bin
        df = pd.DataFrame({almanac_bin: [0, 1, 0, 4], hotspots_bin: [0, 1, 0, 0], cgc_bin: [0, 1, 1, 0]})
        self.assertEqual([1, 2, 3], list(GermlineCancer.get_cancer_idx(df)))
