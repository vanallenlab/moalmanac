import features
import pandas as pd
import pytest
from config import COLNAMES
from features import CoverageMetrics
from reader import Ini


class TestFeatures:
    def test_annotate_feature_type(self):
        feature_type = "A"
        idx = [0, 1, 2]
        series = features.Features.annotate_feature_type(feature_type, idx)
        assert series.tolist() == ["A", "A", "A"]
        assert idx == series.index.tolist()

    def test_calculate_percentile(self):
        array = pd.Series([1, 2, 3, 4, 5])
        assert features.Features.calculate_percentile(array, 99) == 4.96
        assert features.Features.calculate_percentile(array, 1) == 1.04
        assert pytest.approx(features.Features.calculate_percentile(array, 2.5)) == 1.1
        assert pytest.approx(features.Features.calculate_percentile(array, 97.5)) == 4.9

    def test_create_empty_dataframe(self):
        assert features.Features.create_empty_dataframe().dropna().empty

    def test_create_empty_series(self):
        assert features.Features.create_empty_series().dropna().empty

    def test_drop_duplicate_genes(self):
        feature = features.Features.feature
        dataframe = pd.DataFrame(
            {feature: ["ABL1", "BRAF", "BRAF", "TP53"], "sort": [0, 0, 1, 1]},
        )
        assert features.Features.drop_duplicate_genes(
            dataframe,
            "sort",
        ).tolist() == [2, 3, 0]

    def test_import_if_path_exists(self):
        handle = "../example_data/example_patient.capture.somatic.snvs.maf"
        column_map = features.MAFSomatic.create_column_map("tcga_maf_input")
        imported_dataframe = features.Features.import_if_path_exists(
            handle,
            "\t",
            column_map,
        )
        empty_dataframe = features.Features.import_if_path_exists("", "\t", column_map)
        empty_dataframe = features.Features.preallocate_missing_columns(empty_dataframe)
        assert sorted(empty_dataframe.columns.tolist()) == sorted(
            imported_dataframe.columns.tolist(),
        )
        assert imported_dataframe.shape[0] == 299
        assert empty_dataframe.shape[0] == 0

    def test_preallocate_missing_columns(self):
        columns = features.Features.all_columns
        empty_dataframe = pd.DataFrame()
        partial_dataframe = pd.DataFrame(columns=columns[:4])
        assert (
            columns
            == features.Features.preallocate_missing_columns(
                empty_dataframe,
            ).columns.tolist()
        )
        assert (
            columns
            == features.Features.preallocate_missing_columns(
                partial_dataframe,
            ).columns.tolist()
        )


class TestAneuploidy:
    feature = features.Features.feature
    feature_type = features.Features.feature_type


class TestCopyNumber:
    feature_type = features.Features.feature_type

    def test_format_cn_gene(self):
        series = pd.Series(["A", "A B", "B  C"])
        assert features.CopyNumber.format_cn_gene(series).tolist() == ["A", "A", "B"]


class TestCopyNumberCalled:
    def test_create_column_map(self):
        column_map = features.CopyNumberCalled.create_column_map()
        assert column_map["gene"] == features.Features.feature
        assert column_map["call"] == features.Features.alt_type

    def test_filter_calls(self):
        amp_string = "Amplification"
        del_string = "Deletion"
        tmp = pd.Series(["Amplification", "Deletion", "", "Deletion"])
        idx = features.CopyNumberCalled.filter_calls(
            series=tmp,
            amp_string=amp_string,
            del_string=del_string,
        )
        assert idx[idx].index.tolist() == [0, 1, 3]
        assert idx[~idx].index.tolist() == [2]


class TestCopyNumberTotal:
    def test_annotate_amp_del(self):
        amp_string = "Amplification"
        del_string = "Deletion"
        index = pd.Index([0, 1, 2])
        index_amp = pd.Index([0])
        index_del = pd.Index([2])
        expected = [amp_string, "", del_string]
        result = features.CopyNumberTotal.annotate_amp_del(
            idx=index,
            idx_amp=index_amp,
            idx_del=index_del,
            amp_string=amp_string,
            del_string=del_string,
        )
        assert expected == result.tolist()

    def test_create_column_map(self):
        column_map = features.CopyNumberTotal.create_column_map()
        column_names = COLNAMES["seg_input"]
        assert column_map[column_names["gene"]] == features.Features.feature
        assert column_map[column_names["contig"]] == features.Features.chr
        assert column_map[column_names["start"]] == features.Features.start
        assert column_map[column_names["end"]] == features.Features.end
        assert column_map[column_names["mean"]] == features.Features.segment_mean
        assert column_map[column_names["sample"]] == features.Features.tumor

    def test_drop_duplicate_genes(self):
        genes = ["CDKN2A", "CDKN2A", "BRAF", "BRAF", "", pd.NA, "TP53"]
        segments = [2, 1, 2, 3, 0, -1, 2]
        expected = pd.Index([0, 3, 4, 5, 6])
        df = pd.DataFrame(
            {
                features.Features.feature: genes,
                features.Features.segment_mean: segments,
            },
        )
        for idx in expected:
            assert idx in features.CopyNumberTotal.drop_duplicate_genes(df)

    def test_filter_by_threshold(self):
        amp_string = "Amplification"
        del_string = "Deletion"
        values = pd.Series(range(1, 101))
        df = pd.DataFrame(
            {
                features.Features.feature: values,
                features.Features.chr: values,
                features.Features.start: values,
                features.Features.segment_mean: values,
            },
        )
        accept, _ = features.CopyNumberTotal.filter_by_threshold(
            df=df,
            percentile_amp=97.5,
            percentile_del=2.5,
            amp_string=amp_string,
            del_string=del_string,
        )
        expected = pd.Index([0, 1, 2, 97, 98, 99])
        for idx in expected:
            assert idx in accept.index

    def test_get_unique_segments(self):
        df = pd.DataFrame(
            {
                features.Features.chr: [1, 1, 2, 2, 3, 4],
                features.Features.start: [0, 0, 1, 2, 3, 0],
                features.Features.segment_mean: [0, 1, 2, 3, 4, 5],
            },
        )
        expected = [0, 2, 3, 4, 5]
        assert expected == features.CopyNumberTotal.get_unique_segments(df).tolist()


class TestCosmicSignatures:
    def test_round_contributions(self):
        contributions = {
            "COSMIC signature 1": 0.0000,
            "COSMIC signature 3": 0.5001,
            "COSMIC signature 7": 0.9999,
            "COSMIC signature 10": 0.2456,
            "COSMIC signature 11": 0.1985,
        }
        series = pd.Series(contributions)
        result = features.CosmicSignatures.round_contributions(series)
        assert result.shape[0] == 5
        assert result.loc["COSMIC signature 1"] == 0.000
        assert result.loc["COSMIC signature 3"] == 0.500
        assert result.loc["COSMIC signature 7"] == 1.000
        assert result.loc["COSMIC signature 10"] == 0.246
        assert result.loc["COSMIC signature 11"] == 0.198

    def test_subset_significant_signatures(self):
        contributions = {
            "COSMIC signature 1": 0,
            "COSMIC signature 3": 0.5,
            "COSMIC signature 7": 1.0,
            "COSMIC signature 10": 0.2,
            "COSMIC signature 11": 0.19,
        }
        series = pd.Series(contributions)
        idx = features.CosmicSignatures.index_for_minimum_contribution(
            series,
            minimum_value=0.20,
        )
        result = series[idx]
        assert result.shape[0] == 3
        assert result.index[0] == "COSMIC signature 3"
        assert result.index[1] == "COSMIC signature 7"
        assert result.index[2] == "COSMIC signature 10"


class TestCoverageMetrics:
    def test_apply_min_coverage_for_onps(self):
        """
        This function is called by features.CoverageMetrics.format_coverage_col
        after it replaces '__UNKNOWN__' and '' with pd.NA
        """
        series = pd.Series([3, pd.NA, "4|5", "4|", "5"])
        result = CoverageMetrics.apply_min_coverage_for_onps(series)
        assert result.loc[0] == 3
        assert result.loc[2] == 4
        assert result.loc[3] == 4
        assert result.loc[4] == "5"
        assert isinstance(result.loc[1], type(pd.NA))

    def test_calculate_coverage(self):
        """
        This function is called after features.CoverageMetrics.format_coverage_col
        so values should be pd.NA or integers
        """
        series_alt_counts = pd.Series([0, 1, pd.NA, 3])
        series_ref_counts = pd.Series([pd.NA, 0, 1, 2])
        result = CoverageMetrics.calculate_coverage(
            series_alt_count=series_alt_counts,
            series_ref_count=series_ref_counts,
        )
        assert result.loc[0] == 0
        assert result.loc[1] == 1
        assert result.loc[2] == 1
        assert result.loc[3] == 5

    def test_calculate_tumor_f(self):
        series_alt_count = pd.Series([0, 1, pd.NA, 3])
        series_total_coverage = pd.Series([0, 6, 1, 6])
        result = CoverageMetrics.calculate_tumor_f(
            series_alt_count=series_alt_count,
            series_total_coverage=series_total_coverage,
        )
        assert result.loc[1] == 0.1667
        assert result.loc[3] == 0.5000
        assert 0 not in result.index
        assert 2 not in result.index

    def test_convert_to_pandas_int64(self):
        series = pd.Series([3, None, pd.NA, "4", "5"])
        result = CoverageMetrics.convert_to_pandas_int64(series)
        assert result.loc[0] == 3
        assert result.loc[3] == 4
        assert result.loc[4] == 5
        assert isinstance(result.loc[1], type(pd.NA))
        assert isinstance(result.loc[2], type(pd.NA))

    def test_format_coverage_col(self):
        series = pd.Series([3, "", None, pd.NA, "4|5", "4|", "5", "__UNKNOWN__"])
        result = CoverageMetrics.format_coverage_col(series)
        assert result.loc[0] == 3
        assert isinstance(result.loc[1], type(pd.NA))
        assert isinstance(result.loc[2], type(pd.NA))
        assert isinstance(result.loc[3], type(pd.NA))
        assert result.loc[4] == 4
        assert result.loc[5] == 4
        assert result.loc[6] == 5
        assert isinstance(result.loc[7], type(pd.NA))

    def test_get_minimum_values_dataframe(self):
        """
        This function is run after format_coverage_col, apply_min_coverage_for_onps,
        get_onp_boolean, split_counts.
        It is also only run on the index values that contain an ONP boolean.
        """
        series = pd.Series([3, "", None, pd.NA, "4|5", "4|", "5"])
        dataframe = features.CoverageMetrics.split_counts(series)
        result = features.CoverageMetrics.get_minimum_values_dataframe(
            dataframe,
            axis=1,
        )
        assert result.loc[0] == 3
        assert result.loc[4] == 4
        assert result.loc[5] == 4
        assert result.loc[6] == 5
        assert 0 in result.index
        assert 1 not in result.index
        assert 2 not in result.index
        assert 3 not in result.index
        assert 4 in result.index
        assert 5 in result.index
        assert 6 in result.index

    def test_get_onp_boolean(self):
        series = pd.Series([3, "", None, pd.NA, "4|5", "4|", "|5"])
        result = features.CoverageMetrics.get_onp_boolean(series)
        assert not result.loc[0]
        assert not result.loc[1]
        assert not result.loc[2]
        assert not result.loc[3]
        assert result.loc[4]
        assert result.loc[5]
        assert result.loc[6]

    def test_safe_cast(self):
        series = pd.Series([3, "", None, pd.NA, "4|5", "4|", "5", "__UNKNOWN__"])
        result = features.CoverageMetrics.safe_cast(series.loc[0])
        assert result == 3
        result = features.CoverageMetrics.safe_cast(series.loc[6])
        assert result == 5
        result = features.CoverageMetrics.safe_cast(series.loc[1])
        assert isinstance(result, type(pd.NA))
        result = features.CoverageMetrics.safe_cast(series.loc[2])
        assert isinstance(result, type(pd.NA))
        result = features.CoverageMetrics.safe_cast(series.loc[3])
        assert isinstance(result, type(pd.NA))
        result = features.CoverageMetrics.safe_cast(series.loc[4])
        assert isinstance(result, type(pd.NA))
        result = features.CoverageMetrics.safe_cast(series.loc[5])
        assert isinstance(result, type(pd.NA))
        result = features.CoverageMetrics.safe_cast(series.loc[7])
        assert isinstance(result, type(pd.NA))

    def test_split_counts(self):
        series = pd.Series([3, "", None, pd.NA, "4|5"])
        result = features.CoverageMetrics.split_counts(series)
        assert series.shape[0] == result.shape[0]
        assert result.shape[1] == 2
        assert result.loc[0, 0] == "3"
        assert result.loc[4, 0] == "4"
        assert result.loc[4, 1] == "5"

        assert result.loc[1, 0] == ""
        assert result.loc[2, 0] == ""
        assert result.loc[3, 0] == ""
        assert pd.isna(result.loc[0, 1])
        assert pd.isna(result.loc[1, 1])
        assert pd.isna(result.loc[2, 1])
        assert pd.isna(result.loc[3, 1])


class TestFusion:
    def test_create_column_map(self):
        config = Ini.read(
            path="config.ini",
            extended_interpolation=False,
            convert_to_dictionary=False,
        )
        column_map = features.Fusion.create_colmap(config)
        leftbreakpoint = "leftbreakpoint"
        rightbreakpoint = "rightbreakpoint"

        values = list(column_map.values())
        assert len(column_map) == 4
        assert features.Features.feature == values[0]
        assert features.Features.spanningfrags == values[1]
        assert leftbreakpoint == values[2]
        assert rightbreakpoint == values[3]

    def test_filter_by_spanning_fragment_count(self):
        series = pd.Series([4, 5, 6])
        result = features.Fusion.filter_by_spanning_fragment_count(series, 5.0)
        assert 0 not in result.tolist()
        assert 1 in result.tolist()
        assert 2 in result.tolist()

    def test_split_genes(self):
        series = pd.Series(["foo--bar", "foo-bar", "A--B"])
        result = features.Fusion.split_genes(series)
        assert result.loc[0, features.Features.left_gene] == "foo"
        assert result.loc[0, features.Features.right_gene] == "bar"
        assert result.loc[1, features.Features.left_gene] == "foo-bar"
        assert pd.isna(result.loc[1, features.Features.right_gene])
        assert result.loc[2, features.Features.left_gene] == "A"
        assert result.loc[2, features.Features.right_gene] == "B"

    def test_split_breakpoint(self):
        series = pd.Series(["9:100:+", "chr9:100", ""])
        result = features.Fusion.split_breakpoint(series)
        assert result.loc[0, features.Features.chr] == "9"
        assert result.loc[0, features.Features.start] == "100"
        assert result.loc[1, features.Features.chr] == "chr9"
        assert result.loc[1, features.Features.start] == "100"
        assert result.loc[2, features.Features.chr] == ""
        assert pd.isna(result.loc[2, features.Features.start])


class TestMAF:
    def test_check_format(self):
        df = pd.DataFrame(columns=["Protein_Change", ""])
        assert features.MAF.check_format(df) == "tcga_maf_input"
        df = pd.DataFrame(columns=["HGVSp_Short", ""])
        assert features.MAF.check_format(df) == "gdc_maf_input"
        df = pd.DataFrame(columns=[])
        with pytest.raises(SystemExit) as system_exit:
            features.MAF.check_format(df)
        code = (
            "Neither 'Protein_Change' nor 'HGVSp_Short' are present columns, "
            "cannot map to MAF format"
        )
        assert system_exit.value.code == code

    def test_create_column_map(self):
        tcga_column_map = features.MAF.create_column_map("tcga_maf_input")
        gdc_column_map = features.MAF.create_column_map("gdc_maf_input")
        tcga_expected_key = "protein_change" in list(tcga_column_map.keys())
        gdc_expected_key = "hgvsp_short" in list(gdc_column_map.keys())
        tcga_unexpected_key = "hgvsp_short" in list(tcga_column_map.keys())
        gdc_unexpected_key = "protein_change" in list(gdc_column_map.keys())
        assert tcga_expected_key
        assert gdc_expected_key
        assert not tcga_unexpected_key
        assert not gdc_unexpected_key

    def test_format_maf(self):
        data = {
            features.MAF.alt_type: ["Missense_Mutation", "Nonsense_Mutation", "A"],
            features.MAF.alt_count: [0, 0, 0],
            features.MAF.ref_count: [1, 1, 1],
        }
        initial = pd.DataFrame(data)
        result = features.MAF.format_maf(initial, "Somatic Variant")
        assert result.loc[0, features.MAF.alt_type] == "Missense"
        assert result.loc[1, features.MAF.alt_type] == "Nonsense"
        assert result.loc[2, features.MAF.alt_type] == "A"
        assert (
            result[features.MAF.feature_type].unique().tolist()[0] == "Somatic Variant"
        )
        assert result[features.MAF.feature_type].drop_duplicates().shape[0] == 1

    def test_rename_coding_classifications(self):
        alt_types = [
            "Missense_Mutation",
            "Nonsense_Mutation",
            "Nonstop_Mutation",
            "Splice_Site",
            "Frame_Shift_Ins",
            "Frame_Shift_Del",
            "In_Frame_Ins",
            "In_Frame_Del",
            "",
            "A",
            5,
            pd.NA,
        ]
        initial = pd.Series(alt_types)
        result = features.MAF.rename_coding_classifications(initial)
        assert result.loc[0] == "Missense"
        assert result.loc[1] == "Nonsense"
        assert result.loc[2] == "Nonstop"
        assert result.loc[3] == "Splice Site"
        assert result.loc[4] == "Frameshift"
        assert result.loc[5] == "Frameshift"
        assert result.loc[6] == "Insertion"
        assert result.loc[7] == "Deletion"
        assert result.loc[8] == ""
        assert result.loc[9] == "A"
        assert result.loc[10] == "5"

    def test_return_idx_variants_coding(self):
        initial = pd.Series(["Missense", "Silent", "A", 1, pd.NA])
        result = features.MAF.return_idx_variants_coding(initial)
        assert result.loc[0]
        for index in [1, 2, 3]:
            assert not (result.loc[index])

    def test_return_variants_coding(self):
        initial = pd.DataFrame(
            {features.MAF.alt_type: ["Missense", "Silent", "A", 1, pd.NA]},
        )
        result = features.MAF.return_variants_coding(initial)
        assert result.index.tolist() == [0]

    def test_return_variants_non_coding(self):
        initial = pd.DataFrame(
            {features.MAF.alt_type: ["Missense", "Silent", "A", 1, pd.NA]},
        )
        result = features.MAF.return_variants_non_coding(initial)
        assert result.index.tolist() == [1, 2, 3, 4]
