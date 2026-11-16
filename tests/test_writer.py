import pandas as pd
from writer import GermlineCancer
from writer import Writer


class TestWriter:
    def test_create_output_name(self):
        string1 = "."
        string2 = "Foo"
        string3 = "Bar"
        assert Writer.create_output_name(string1, string2, string3) == "Foo.Bar"

    def test_sort_columns(self):
        df = pd.DataFrame({"A": [0, 1, 2], "B": [1, 1, 2]})
        assert list(
            Writer.sort_columns(df, ["A", "B"], [True, True]).index,
        ) == [0, 1, 2]
        assert list(
            Writer.sort_columns(df, ["A", "B"], [False, True]).index,
        ) == [2, 1, 0]
        assert list(
            Writer.sort_columns(df, ["B", "A"], [True, False]).index,
        ) == [1, 0, 2]
        assert list(
            Writer.sort_columns(df, ["B", "A"], [False, True]).index,
        ) == [2, 0, 1]

    def test_return_nonzero_idx(self):
        series = pd.Series([0, 1, 2, 0])
        assert list(Writer.return_nonzero_bin_idx(series)) == [1, 2]


class TestGermlineCancer:
    def test_get_cancer_idx(self):
        almanac_bin = GermlineCancer.almanac_bin
        hotspots_bin = GermlineCancer.hotspots_bin
        cgc_bin = GermlineCancer.cgc_bin
        df = pd.DataFrame(
            {
                almanac_bin: [0, 1, 0, 4],
                hotspots_bin: [0, 1, 0, 0],
                cgc_bin: [0, 1, 1, 0],
            },
        )
        assert list(GermlineCancer.get_cancer_idx(df)) == [1, 2, 3]
