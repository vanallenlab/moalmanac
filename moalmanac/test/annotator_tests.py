import unittest
import operator as op
import pandas as pd
import scipy.stats as stats

from annotator import Annotator, ACMG, Almanac, ExAC, OverlapValidation, PreclinicalEfficacy, PreclinicalMatchmaking
from datasources import Datasources
from datasources import Almanac as datasource_Almanac
from datasources import Preclinical as datasources_Preclinical
from features import Features
from investigator import SensitivityDictionary


class UnitTestAnnotator(unittest.TestCase):
    def test_compare_ids_true_index(self):
        df1 = pd.DataFrame({"chromosome": [1, 2, 3], "reference_allele": ["C", "A", "T"]})
        df2 = pd.DataFrame({"chromosome": [1, 2, 3], "reference_allele": ["C", "G", "T"]})
        id1 = Annotator.create_id(df1, ["chromosome", "reference_allele"])
        id2 = Annotator.create_id(df2, ["chromosome", "reference_allele"])
        result = Annotator.compare_ids_true_index(id1, id2)
        self.assertEqual([0, 2], result.tolist())

    def test_create_id(self):
        df = pd.DataFrame({"chromosome": [1, 2, 3], "reference_allele": ["C", "A", "T"]})
        id1 = Annotator.create_id(df, ["chromosome"])
        id2 = Annotator.create_id(df, ["chromosome", "reference_allele"])
        self.assertEqual(['1', '2', '3'], id1.tolist())
        self.assertEqual(['1_C', '2_A', '3_T'], id2.tolist())

    def test_create_id_series(self):
        df = pd.Series({'gene_name': 'TP53', 'chromosome': 1, 'start_position': 100})
        id1 = Annotator.create_id_series(df, ['gene_name', 'chromosome', 'start_position'])
        id2 = Annotator.create_id_series(df, ['gene_name'])
        id3 = Annotator.create_id_series(df, ['chromosome'])
        self.assertEqual('TP53_1_100', id1)
        self.assertEqual('TP53', id2)
        self.assertEqual('1', id3)

    def test_preallocate_bin(self):
        index = [0, 1, 5, 7]
        column_list = 'Column'
        series = pd.Series(0, name=column_list, index=index)
        result = Annotator.preallocate_bin(column_list, index)
        self.assertEqual(True, series.equals(result))

    def test_preallocate_empty_columns(self):
        index = [0, 1, 5, 7]
        columns = ['Column1', 'Column2', 'Column3']
        dataframe1 = pd.DataFrame(None, index=index, columns=columns)
        dataframe2 = pd.DataFrame(None, index=index, columns=columns[0:2])
        result = Annotator.preallocate_empty_columns(dataframe2, [columns[2]])
        self.assertEqual(True, dataframe1.equals(result))


class UnitTestACMG(unittest.TestCase):
    def test_annotate(self):
        gene = ACMG.gene
        bin_name = ACMG.bin_name
        df = pd.DataFrame({gene: ['TP53', 'FOO', 'PMS2', 'TSC1', 'AR']})
        dbs = {
            'acmg_handle': '../datasources/acmg/acmg.secondaryfindings.v3.txt'
        }

        annotated = ACMG.annotate(df, dbs)
        expected_result = pd.Series([1, 0, 1, 1, 0], name=bin_name)
        self.assertEqual(True, annotated[bin_name].equals(expected_result))


class UnitTestAlmanac(unittest.TestCase):
    def test_extract_values_from_list_of_dicts(self):
        dict1 = {'A': 0, 'B': 1, 'C': 2}
        dict2 = {'A': 1, 'B': 2, 'C': 3}
        list_dicts = [dict1, dict2]
        for key, value in [('A', 1), ('B', 3), ('C', 5)]:
            values = Almanac.extract_values_from_list_of_dicts(key, list_dicts)
            self.assertEqual(value, sum(values))

    def test_select_better_fusion_match(self):
        list_of_dictionaries = [{'bin': -1, 'a': 1}, {'bin': 0, 'a': 2}, {'bin': 1, 'a': 3}, {'bin': 1, 'a': 4}]
        better_match = Almanac.select_better_fusion_match(list_of_dictionaries)
        self.assertEqual(1, better_match['bin'])

    def test_sort_and_subset_matches(self):
        implication_map = Almanac.implication_map
        publication_date = Almanac.publication_date
        last_updated = Almanac.last_updated

        table = pd.DataFrame([
            {implication_map: 1, publication_date: 0, 'last_updated': 0, 'label': 'a'},
            {implication_map: 2, publication_date: 0, 'last_updated': 0, 'label': 'b'},
            {implication_map: 2, publication_date: 1, 'last_updated': 0, 'label': 'c'},
            {implication_map: 0, publication_date: 0, 'last_updated': 0, 'label': 'd'},
            {implication_map: 2, publication_date: 1, 'last_updated': 0, 'label': 'e'}
        ])
    
        # key, value pair in each record is if the index matched to the corresponding record in `table`
        records1 = pd.Series({0: True, 1: True, 2: True, 3: False, 4: False})
        records2 = pd.Series({0: False, 1: False, 2: False, 3: True, 4: True})
        records3 = pd.Series({0: False, 1: False, 2: False, 3: False, 4: False})
        self.assertEqual('c', Almanac.sort_and_subset_matches(table, records1, records2)[0]['label'])
        self.assertEqual('e', Almanac.sort_and_subset_matches(table, records3, records2)[0]['label'])

    def test_sort_dictionary_as_dataframe(self):
        records = [{'A': 0, 'B': 3, 'C': 1}, {'A': 1, 'B': 2, 'C': 0}, {'A': 2, 'B': 1, 'C': 2}]
        result1 = Almanac.sort_dictionary_as_dataframe(records, ['A', 'B', 'C'], [True, True, True])
        result2 = Almanac.sort_dictionary_as_dataframe(records, ['A', 'B', 'C'], [False, False, False])
        result3 = Almanac.sort_dictionary_as_dataframe(records, ['C', 'B', 'A'], [True, True, False])
        answer1 = [{'A': 0, 'B': 3, 'C': 1}, {'A': 1, 'B': 2, 'C': 0}, {'A': 2, 'B': 1, 'C': 2}]
        answer2 = [{'A': 2, 'B': 1, 'C': 2}, {'A': 1, 'B': 2, 'C': 0}, {'A': 0, 'B': 3, 'C': 1}]
        answer3 = [{'A': 1, 'B': 2, 'C': 0}, {'A': 0, 'B': 3, 'C': 1}, {'A': 2, 'B': 1, 'C': 2}]
        for answer, result in [(answer1, result1), (answer2, result2), (answer3, result3)]:
            self.assertEqual(answer, result)

    def test_subset_dataframe_by_condition(self):
        df = pd.DataFrame({'A': [-1, 0, 1, 2], 'B': ['cat', 'dog', 'car', 'dare']})
        for answer, condition in [
            ([2, 3], df['A'].gt(0)),
            ([0, 1], df['A'].le(0)),
            ([0, 2, 3], df['B'].str.contains('a')),
            ([0, 2], df['B'].str.contains('ca')),
            ([1, 3], df['B'].str.contains('d')),
            ([1], df['B'].str.match('dog'))
        ]:
            self.assertEqual(answer, list(df[condition].index))

    def test_subset_list_of_dictionaries(self):
        list_of_dicts = [{'bin': -1, 'a': 1}, {'bin': 0, 'a': 2}, {'bin': 1, 'a': 3}, {'bin': 1, 'a': 4}]
        result1 = Almanac.subset_list_of_dictionaries(list_of_dicts, 'bin', op.gt, 0)
        result2 = Almanac.subset_list_of_dictionaries(list_of_dicts, 'bin', op.le, 0)
        result3 = Almanac.subset_list_of_dictionaries(list_of_dicts, 'a', op.ge, 2)
        answer1 = [{'bin': 1, 'a': 3}, {'bin': 1, 'a': 4}]
        answer2 = [{'bin': -1, 'a': 1}, {'bin': 0, 'a': 2}]
        answer3 = [{'bin': 0, 'a': 2}, {'bin': 1, 'a': 3}, {'bin': 1, 'a': 4}]

        for answer, result in [(answer1, result1), (answer2, result2), (answer3, result3)]:
            self.assertEqual(answer, result)

    def test_return_best_evidence_level(self):
        implication_map = Almanac.implication_map
        list_of_dicts1 = [{implication_map: 1}, {implication_map: 3}, {implication_map: 2}, {implication_map: 1}]
        list_of_dicts2 = [{implication_map: 4}, {implication_map: 4}, {implication_map: 3}, {implication_map: 1}]
        list_of_dicts3 = []
        self.assertEqual(3, Almanac.return_best_evidence_level(list_of_dicts1, list_of_dicts2))
        self.assertEqual(4, Almanac.return_best_evidence_level(list_of_dicts2, list_of_dicts1))
        self.assertEqual(3, Almanac.return_best_evidence_level(list_of_dicts1, list_of_dicts3))
        self.assertEqual(3, Almanac.return_best_evidence_level(list_of_dicts3, list_of_dicts1))

    def test_update_series_with_best_match(self):
        matches = [{"gene": "ABL1", "chromosome": "9", "start_position": "133747580", "end_position": "133747580",
                    "reference_allele": "C", "alternate_allele": "T", "cdna_change": "c.944C>T",
                    "protein_change": "p.T315I", "variant_annotation": "Missense", "exon": "5",
                    "rsid": '', "disease": "Chronic Myeloid Leukemia", "context": '',
                    "oncotree_term": "Chronic Myelogenous Leukemia", "oncotree_code": "CML",
                    "therapy_name": "Imatinib",
                    "therapy_strategy": "BCR-ABL inhibition", "therapy_type": "Targeted therapy",
                    "therapy_sensitivity": '',
                    "therapy_resistance": "1", "favorable_prognosis": '',
                    "predictive_implication": "Preclinical",
                    "description": "the quick fox caught the hen", "preferred_assertion": '',
                    "source_type": "Journal",
                    "citation": "peer reviewed", "url": "https://doi.org/10.1126/science.1062538",
                    "doi": "10.1126/science.1062538", "pmid": "11423618", "nct": '',
                    "last_updated": "6/13/19",
                    "feature_display": "ABL1 p.T315I (Missense)", "predictive_implication_map": 1.0}, {}]
        somatic_variant = 'Somatic Variant'
        series = pd.Series(dtype=object)

        for columns in [Almanac.column_map_sensitive, Almanac.column_map_resistance, Almanac.column_map_prognostic]:
            series = Almanac.update_series_with_best_match(matches, columns, series)
            self.assertEqual('ABL1 p.T315I (Missense)', series.loc[columns['feature_display']])
            self.assertEqual('Chronic Myelogenous Leukemia', series.loc[columns['oncotree_term']])
            self.assertEqual('CML', series.loc[columns['oncotree_code']])
            self.assertEqual('', series.loc[columns['context']])
            self.assertEqual('Preclinical', series.loc[columns['predictive_implication']])
            self.assertEqual(1.0, series.loc[columns['predictive_implication_map']])
            self.assertEqual('the quick fox caught the hen', series.loc[columns['description']])
            self.assertEqual('Journal', series.loc[columns['source_type']])
            self.assertEqual('peer reviewed', series.loc[columns['citation']])
            self.assertEqual('https://doi.org/10.1126/science.1062538', series.loc[columns['url']])
            self.assertEqual('10.1126/science.1062538', series.loc[columns['doi']])
            self.assertEqual('11423618', series.loc[columns['pmid']])
            self.assertEqual('', series.loc[columns['nct']])

        for columns in [Almanac.column_map_sensitive, Almanac.column_map_resistance]:
            self.assertEqual('Imatinib', series.loc[columns['therapy_name']])
            self.assertEqual('Targeted therapy', series.loc[columns['therapy_type']])
        self.assertEqual('', series.loc[Almanac.column_map_prognostic[Almanac.prognosis]])


class UnitTestExAC(unittest.TestCase):
    def test_append_exac_af(self):
        chr = ExAC.chr
        start = ExAC.start
        ref = ExAC.ref
        alt = ExAC.alt
        af = ExAC.af
        feature_type = Features.feature_type
        somatic = 'Somatic Variant'
        germline = 'Germline Variant'
        cn = 'Copy Number'

        df = pd.DataFrame({chr: [1, 2, 3, 1],
                           start: [100, 101, 103, 100],
                           ref: ["C", "A", "T", "C"],
                           alt: ["G", "G", "G", "G"],
                           feature_type: [somatic, germline, somatic, cn]
                           })
        exac = pd.DataFrame({chr: [1, 2, 4],
                             start: [100, 99, 103],
                             ref: ["C", "A", "T"],
                             alt: ["G", "G", "G"],
                             af: [1, 0.5, 0.001]})
        biomarker_types = [somatic, germline]
        result = ExAC.append_exac_af(
            df=df,
            ds=exac,
            ds_columns=[chr, start, ref, alt, af],
            variant_biomarker_types=biomarker_types)
        self.assertEqual([1, 0, 0, 0], result[af].tolist())

    def test_annotate_common_af(self):
        exac_common_threshold = 0.001
        series = pd.Series([float(exac_common_threshold) - 0.01, float(exac_common_threshold) + 0.01])
        result = ExAC.annotate_common_af(series, threshold=exac_common_threshold)
        self.assertEqual(0.0, result.loc[0])
        self.assertEqual(1.0, result.loc[1])


class UnitTestValidation(unittest.TestCase):
    columns = ['feature', 'feature_type', 'alteration_type', 'alteration',
               'tumor_f', 'total_coverage',
               'validation_tumor_f', 'validation_total_coverage']

    dataframe1 = pd.DataFrame({
        columns[0]: ['BRAF', 'KRAS', 'NRAS', 'TP53'],
        columns[1]: ['Somatic Variant', 'Copy Number', 'Somatic Variant', 'Somatic Variant'],
        columns[2]: ['Missense', 'Deletion', 'Nonsense', 'Missense'],
        columns[3]: ['p.V600E', '', 'p.Q61*', 'p.f00n'],
        columns[4]: [0.40, pd.NA, 0.30, 0.10],
        columns[5]: [152, pd.NA, 50, 111],
        columns[6]: [pd.NA, pd.NA, pd.NA, pd.NA],
        columns[7]: [pd.NA, pd.NA, pd.NA, pd.NA]
    })

    dataframe2 = pd.DataFrame({
        columns[0]: ['BRAF', 'KRAS', 'NRAS', 'EGFR'],
        columns[2]: ['Missense', 'Missense', 'Nonsense', 'Missense'],
        columns[3]: ['p.V600E', 'p.o00n', 'p.Q61*', 'p.f00n'],
        columns[6]: [0.20, 0.54, 0.66, 0.41],
        columns[7]: [8, 32, 29, 222]
    })

    def test_append_validation(self):
        result = OverlapValidation.append_validation(
            UnitTestValidation.dataframe1,
            UnitTestValidation.dataframe2,
            biomarker_type='Somatic Variant'
        )
        result = result.fillna('')
        self.assertEqual(UnitTestValidation.dataframe1['feature'].tolist(), result['feature'].tolist())
        self.assertEqual([0.20, '', 0.66, 0.0], result['validation_tumor_f'].tolist())
        self.assertEqual([8.0, '', 29.0, 0.0], result['validation_total_coverage'].tolist())
        self.assertEqual([0.4103, '', 0.9715, 0.0], result['validation_detection_power'].tolist())

    def test_calculate_beta_binomial(self):
        self.assertEqual(0.7006851149078945, stats.betabinom.sf(k=3, n=5, a=2.3, b=0.63))
        self.assertEqual(0.7007, round(stats.betabinom.sf(k=3, n=5, a=2.3, b=0.63), 4))
        series_n = pd.Series([10, 50, 100])
        series_a = pd.Series([5, 25, 50])
        series_b = pd.Series([15, 40, 100])
        solutions = list(stats.betabinom.sf(k=3, n=series_n, a=series_a, b=series_b, loc=0))
        self.assertEqual(0.25543941316055174, solutions[0])
        self.assertEqual(0.9999749071547652, solutions[1])
        self.assertEqual(0.9999999996130999, solutions[2])

    def test_calculate_validation_detection_power(self):
        tumor_f = pd.Series([0, 0.33, 0.5, 1])
        coverage = pd.Series([5, 100, 50, 66])
        validation_coverage = pd.Series([12, 70, 40, 20])
        result = OverlapValidation.calculate_validation_detection_power(tumor_f, coverage, validation_coverage)
        self.assertEqual(0.16176470588235325, result[0])
        self.assertEqual(0.999998846495356, result[1])
        self.assertEqual(0.9999945730050155, result[2])
        self.assertEqual(0.9999999999999977, result[3])

    def test_drop_validation_columns(self):
        columns = ['A', 'B', OverlapValidation.validation_tumor_f, OverlapValidation.validation_coverage, 'C']
        solution = pd.Index(['A', 'B', 'C'])
        dataframe = pd.DataFrame(columns=columns)
        result = OverlapValidation.drop_validation_columns(dataframe)
        self.assertEqual(True, solution.equals(result.columns))

    def test_get_mutation_index(self):
        dataframe = pd.DataFrame(['Somatic Variant', 'bar', 'foo'], columns=[OverlapValidation.feature_type])
        solution = ['Somatic Variant']
        solution_index = pd.Index([0])
        result = OverlapValidation.get_mutation_index(dataframe, biomarker_type='Somatic Variant')
        self.assertEqual(solution[0], dataframe.loc[result[0], OverlapValidation.feature_type])
        self.assertEqual(solution_index[0], result[0])

    def test_merge_data_frames(self):
        dataframe1 = pd.DataFrame(columns=['A', 'B', 'C', 'D'])
        dataframe2 = pd.DataFrame(columns=['A', 'B', 'E'])
        dataframe1['A'] = ['a', 'b', 'c']
        dataframe1['B'] = [1, 2, 3]
        dataframe1['C'] = ['d', 'e', 'f']
        dataframe1['D'] = [4, 5, 6]
        dataframe2['A'] = ['a', 'b', 'g']
        dataframe2['B'] = [1, 4, 5]
        dataframe2['E'] = ['hello', 'goodbye', 'good day']
        result = OverlapValidation.merge_data_frames(dataframe1, dataframe2, ['A', 'B']).fillna(pd.NA)
        self.assertEqual(['hello', pd.NA, pd.NA], result['E'].tolist())
        self.assertEqual(['A', 'B', 'C', 'D', 'E'], result.columns.tolist())

    def test_round_series(self):
        series = pd.Series([0.000000, 0.1111111, 0.66888])
        result = OverlapValidation.round_series(series, 2)
        self.assertEqual(0.00, result.loc[0])
        self.assertEqual(0.11, result.loc[1])
        self.assertEqual(0.67, result.loc[2])


class UnitTestPreclinicalEfficacy(unittest.TestCase):
    data_dictionary = {
        'feature_type': ['Somatic Variant', 'Rearrangement'],
        'feature': ['BRAF', 'COL1A1'],
        'alteration_type': ['Missense', 'Fusion'],
        'alteration': ['p.V600E', 'COL1A1--CITED4'],
        'feature_display': ['BRAF p.V600E', 'COL1A1--CITED4 Fusion'],
        'sensitive_therapy_name': ['Dabrafenib + Trametinib', 'Imatinib'],
        'preclinical_efficacy_observed': [1, 0]
    }
    df1 = pd.DataFrame(data_dictionary, index=[0, 1])
    data_dictionary = {
        'feature_display': ['BRAF p.V600E', 'BRAF p.V600E', 'COL1A1--CITED4'],
        'tested_subfeatre': ['BRAF', 'BRAF Somatic Variant', 'COL1A1'],
        'pvalue_mww': [2.322E-12, 7.627E-17, 0.835]
    }
    df2 = pd.DataFrame(data_dictionary, index=[0, 1, 2])
    dbs_dictionary = {
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
    config = {
        'feature_types': {
            'mut': 'Somatic Variant',
            'cna': 'Copy Number',
            'fusion': 'Rearrangement'
        }
    }
    dbs_preclinical = datasources_Preclinical.import_dbs(dbs_dictionary)
    efficacy_dictionary = SensitivityDictionary.create(dbs_preclinical, df1, config=config)

    def test_annotate(self):
        result = PreclinicalEfficacy.annotate(
            actionable=self.df1,
            efficacy=self.df2,
            dictionary=self.efficacy_dictionary,
        )

        result_braf = result.loc[0, 'preclinical_efficacy_lookup'][0]
        result_pvalue_index_0 = float(result_braf['Dabrafenib']['BRAF']['comparison']['pvalue_mww'])
        result_pvalue_index_1 = float(result_braf['Dabrafenib']['BRAF Somatic Variant']['comparison']['pvalue_mww'])
        self.assertEqual(result_pvalue_index_0, float(self.df2.loc[0, 'pvalue_mww']))
        self.assertEqual(result_pvalue_index_1, float(self.df2.loc[1, 'pvalue_mww']))

    def test_series_for_significance(self):
        series1 = pd.Series([0, 0.06, pd.NA])
        series2 = pd.Series([0.65, 0.062, pd.NA])
        self.assertEqual(PreclinicalEfficacy.search_for_significance(series1), 1)
        self.assertEqual(PreclinicalEfficacy.search_for_significance(series2), 0)


class UnitTestPreclinicalMatchmaking(unittest.TestCase):
    def test_annotate_copy_numbers(self):
        dbs = {
            'almanac_handle': '../datasources/moalmanac/molecular-oncology-almanac.json',
            'cancerhotspots_handle': '../datasources/cancerhotspots/hotspots_v2.txt',
            '3dcancerhotspots_handle': '../datasources/cancerhotspots/hotspots3d.txt',
            'cgc_handle': '../datasources/cancergenecensus/cancer_gene_census_v97.genes.tsv',
            'cosmic_handle': '../datasources/cosmic/CosmicMutantExport_v97.lite.txt',
            'gsea_pathways_handle': '../datasources/gsea_gene_sets/GSEA_cancer_gene_sets.txt',
            'gsea_modules_handle': '../datasources/gsea_gene_sets/c4.cm.v6.0.symbols.txt',
            'exac_handle': '../datasources/exac/exac.expanded.r1.txt',
            'acmg_handle': '../datasources/acmg/acmg.secondaryfindings.v3.txt',
            'clinvar_handle': '../datasources/clinvar/variant_summary.lite.txt',
            'hereditary_handle': '../datasources/hereditary/hereditary.txt',
            'oncotree_handle': '../datasources/oncotree/oncotree.2023-03-09.txt',
            'lawrence_handle': '../datasources/lawrence/lawrence_mapped_ontology.txt'
        }
        feature = PreclinicalMatchmaking.feature
        feature_type = PreclinicalMatchmaking.feature_type
        alteration_type = PreclinicalMatchmaking.alteration_type
        alteration = PreclinicalMatchmaking.alteration

        copy_number = 'Copy Number'

        df = pd.DataFrame({
            feature: ['CDKN2A', 'CDKN2A', 'KRAS'],
            alteration_type: ['Deletion', 'Amplification', 'Amplification']
        })
        df[feature_type] = copy_number
        df[alteration] = pd.NA
        result = PreclinicalMatchmaking.annotate_copy_numbers(df, dbs, biomarker_type_string=copy_number)

        expected_cdkn2a_del = {
            'feature_match_1': 1,
            'feature_match_2': 1,
            'feature_match_3': 1,
            'feature_match_4': 0,
            'feature_match': 0,
            'evidence': 2,
            'cancerhotspots_bin': 1,
            'cancerhotspots3D_bin': 0,
            'cgc_bin': 1,
            'cosmic_bin': 1,
            'gsea_pathways_bin': 1,
            'gsea_modules_bin': 1
        }

        expected_cdkn2a_amp = {
            'feature_match_1': 1,
            'feature_match_2': 1,
            'feature_match_3': 0,
            'feature_match_4': 0,
            'feature_match': 0,
            'evidence': 2,
            'cancerhotspots_bin': 1,
            'cancerhotspots3D_bin': 0,
            'cgc_bin': 1,
            'cosmic_bin': 1,
            'gsea_pathways_bin': 1,
            'gsea_modules_bin': 1
        }

        expected_kras_amp = {
            'feature_match_1': 1,
            'feature_match_2': 0,
            'feature_match_3': 0,
            'feature_match_4': 0,
            'feature_match': 0,
            'cancerhotspots_bin': 1,
            'cancerhotspots3D_bin': 1,
            'cgc_bin': 1,
            'cosmic_bin': 1,
            'gsea_pathways_bin': 1,
            'gsea_modules_bin': 1
        }

        for key, value in expected_cdkn2a_del.items():
            self.assertEqual(result.loc[0, key], expected_cdkn2a_del[key])
        for key, value in expected_cdkn2a_amp.items():
            self.assertEqual(result.loc[1, key], expected_cdkn2a_amp[key])
        for key, value in expected_kras_amp.items():
            self.assertEqual(result.loc[2, key], expected_kras_amp[key])

    def test_annotate_fusions(self):
        dbs = {
            'almanac_handle': '../datasources/moalmanac/molecular-oncology-almanac.json',
            'cancerhotspots_handle': '../datasources/cancerhotspots/hotspots_v2.txt',
            '3dcancerhotspots_handle': '../datasources/cancerhotspots/hotspots3d.txt',
            'cgc_handle': '../datasources/cancergenecensus/cancer_gene_census_v97.genes.tsv',
            'cosmic_handle': '../datasources/cosmic/CosmicMutantExport_v97.lite.txt',
            'gsea_pathways_handle': '../datasources/gsea_gene_sets/GSEA_cancer_gene_sets.txt',
            'gsea_modules_handle': '../datasources/gsea_gene_sets/c4.cm.v6.0.symbols.txt',
            'exac_handle': '../datasources/exac/exac.expanded.r1.txt',
            'acmg_handle': '../datasources/acmg/acmg.secondaryfindings.v3.txt',
            'clinvar_handle': '../datasources/clinvar/variant_summary.lite.txt',
            'hereditary_handle': '../datasources/hereditary/hereditary.txt',
            'oncotree_handle': '../datasources/oncotree/oncotree.2023-03-09.txt',
            'lawrence_handle': '../datasources/lawrence/lawrence_mapped_ontology.txt'
        }
        feature = PreclinicalMatchmaking.feature
        feature_type = PreclinicalMatchmaking.feature_type
        alteration_type = PreclinicalMatchmaking.alteration_type
        partner = PreclinicalMatchmaking.partner
        fusion = 'Rearrangement'
        model_id = PreclinicalMatchmaking.model_id

        df = pd.DataFrame({
            feature: ['BCR', 'BCR', 'NTRK1', 'CDKN2A'],
            partner: ['ABL1', 'foo', 'ABL1', 'bar']})
        df[feature_type] = fusion
        df[alteration_type] = 'Fusion'
        df[model_id] = 'case'

        result, group1, group2 = PreclinicalMatchmaking.annotate_fusions(df, dbs, biomarker_type_string=fusion)

        expected_index_0 = {
            'feature_match_1': 1,
            'feature_match_2': 1,
            'feature_match_3': 1,
            'feature_match_4': 1,
            'feature_match': 4,
            'evidence': 'FDA-Approved',
            'cancerhotspots_bin': 0,
            'cancerhotspots3D_bin': 0,
            'cgc_bin': 1,
            'cosmic_bin': 1,
            'gsea_pathways_bin': 1,
            'gsea_modules_bin': 0
        }
        self.assertEqual(result.loc[0, 'feature_match_1'], expected_index_0['feature_match_1'])
        self.assertEqual(result.loc[0, 'feature_match_2'], expected_index_0['feature_match_2'])
        self.assertEqual(result.loc[0, 'feature_match_3'], expected_index_0['feature_match_3'])
        self.assertEqual(result.loc[0, 'feature_match_4'], expected_index_0['feature_match_4'])
        self.assertEqual(result.loc[0, 'evidence'], expected_index_0['evidence'])
        self.assertEqual(result.loc[0, 'cancerhotspots_bin'], expected_index_0['cancerhotspots_bin'])
        self.assertEqual(result.loc[0, 'cancerhotspots3D_bin'], expected_index_0['cancerhotspots3D_bin'])
        self.assertEqual(result.loc[0, 'cgc_bin'], expected_index_0['cgc_bin'])
        self.assertEqual(result.loc[0, 'cosmic_bin'], expected_index_0['cosmic_bin'])
        self.assertEqual(result.loc[0, 'gsea_pathways_bin'], expected_index_0['gsea_pathways_bin'])
        self.assertEqual(result.loc[0, 'gsea_modules_bin'], expected_index_0['gsea_modules_bin'])

        expected_index_1 = {
            'feature_match_1': 1,
            'feature_match_2': 1,
            'feature_match_3': 1,
            'feature_match_4': 0,
            'feature_match': 3,
            'evidence': 'FDA-Approved',
            'cancerhotspots_bin': 0,
            'cancerhotspots3D_bin': 0,
            'cgc_bin': 1,
            'cosmic_bin': 1,
            'gsea_pathways_bin': 1,
            'gsea_modules_bin': 1
        }
        self.assertEqual(result.loc[1, 'feature_match_1'], expected_index_1['feature_match_1'])
        self.assertEqual(result.loc[1, 'feature_match_2'], expected_index_1['feature_match_2'])
        self.assertEqual(result.loc[1, 'feature_match_3'], expected_index_1['feature_match_3'])
        self.assertEqual(result.loc[1, 'feature_match_4'], expected_index_1['feature_match_4'])
        self.assertEqual(result.loc[1, 'evidence'], expected_index_1['evidence'])
        self.assertEqual(result.loc[1, 'cancerhotspots_bin'], expected_index_1['cancerhotspots_bin'])
        self.assertEqual(result.loc[1, 'cancerhotspots3D_bin'], expected_index_1['cancerhotspots3D_bin'])
        self.assertEqual(result.loc[1, 'cgc_bin'], expected_index_1['cgc_bin'])
        self.assertEqual(result.loc[1, 'cosmic_bin'], expected_index_1['cosmic_bin'])
        self.assertEqual(result.loc[1, 'gsea_pathways_bin'], expected_index_1['gsea_pathways_bin'])
        self.assertEqual(result.loc[1, 'gsea_modules_bin'], expected_index_1['gsea_modules_bin'])

        expected_index_2 = {
            'feature_match_1': 1,
            'feature_match_2': 1,
            'feature_match_3': 1,
            'feature_match_4': 0,
            'feature_match': 2,
            'evidence': 'FDA-Approved',
            'cancerhotspots_bin': 0,
            'cancerhotspots3D_bin': 0,
            'cgc_bin': 1,
            'cosmic_bin': 1,
            'gsea_pathways_bin': 1,
            'gsea_modules_bin': 0
        }
        self.assertEqual(result.loc[2, 'feature_match_1'], expected_index_2['feature_match_1'])
        self.assertEqual(result.loc[2, 'feature_match_2'], expected_index_2['feature_match_2'])
        self.assertEqual(result.loc[2, 'feature_match_3'], expected_index_2['feature_match_3'])
        self.assertEqual(result.loc[2, 'feature_match_4'], expected_index_2['feature_match_4'])
        self.assertEqual(result.loc[2, 'evidence'], expected_index_2['evidence'])
        self.assertEqual(result.loc[2, 'cancerhotspots_bin'], expected_index_2['cancerhotspots_bin'])
        self.assertEqual(result.loc[2, 'cancerhotspots3D_bin'], expected_index_2['cancerhotspots3D_bin'])
        self.assertEqual(result.loc[2, 'cgc_bin'], expected_index_2['cgc_bin'])
        self.assertEqual(result.loc[2, 'cosmic_bin'], expected_index_2['cosmic_bin'])
        self.assertEqual(result.loc[2, 'gsea_pathways_bin'], expected_index_2['gsea_pathways_bin'])
        self.assertEqual(result.loc[2, 'gsea_modules_bin'], expected_index_2['gsea_modules_bin'])

        expected_index_3 = {
            'feature_match_1': 1,
            'feature_match_2': 0,
            'feature_match_3': 0,
            'feature_match_4': 0,
            'feature_match': 1,
            'evidence': '',
            'cancerhotspots_bin': 1,
            'cancerhotspots3D_bin': 0,
            'cgc_bin': 1,
            'cosmic_bin': 1,
            'gsea_pathways_bin': 1,
            'gsea_modules_bin': 1
        }
        self.assertEqual(result.loc[3, 'feature_match_1'], expected_index_3['feature_match_1'])
        self.assertEqual(result.loc[3, 'feature_match_2'], expected_index_3['feature_match_2'])
        self.assertEqual(result.loc[3, 'feature_match_3'], expected_index_3['feature_match_3'])
        self.assertEqual(result.loc[3, 'feature_match_4'], expected_index_3['feature_match_4'])
        self.assertEqual(result.loc[3, 'evidence'], expected_index_3['evidence'])
        self.assertEqual(result.loc[3, 'cancerhotspots_bin'], expected_index_3['cancerhotspots_bin'])
        self.assertEqual(result.loc[3, 'cancerhotspots3D_bin'], expected_index_3['cancerhotspots3D_bin'])
        self.assertEqual(result.loc[3, 'cgc_bin'], expected_index_3['cgc_bin'])
        self.assertEqual(result.loc[3, 'cosmic_bin'], expected_index_3['cosmic_bin'])
        self.assertEqual(result.loc[3, 'gsea_pathways_bin'], expected_index_3['gsea_pathways_bin'])
        self.assertEqual(result.loc[3, 'gsea_modules_bin'], expected_index_3['gsea_modules_bin'])

        expected_index_0_group1 = {
            'feature_match_1': 1,
            'feature_match_2': 1,
            'feature_match_3': 1,
            'feature_match_4': 0,
            'feature_match': 3,
            'evidence': 'FDA-Approved',
            'cancerhotspots_bin': 0,
            'cancerhotspots3D_bin': 0,
            'cgc_bin': 1,
            'cosmic_bin': 1,
            'gsea_pathways_bin': 1,
            'gsea_modules_bin': 1
        }
        self.assertEqual(group1.loc[0, 'feature_match_1'], expected_index_0_group1['feature_match_1'])
        self.assertEqual(group1.loc[0, 'feature_match_2'], expected_index_0_group1['feature_match_2'])
        self.assertEqual(group1.loc[0, 'feature_match_3'], expected_index_0_group1['feature_match_3'])
        self.assertEqual(group1.loc[0, 'feature_match_4'], expected_index_0_group1['feature_match_4'])
        self.assertEqual(group1.loc[0, 'evidence'], expected_index_0_group1['evidence'])
        self.assertEqual(group1.loc[0, 'cancerhotspots_bin'], expected_index_0_group1['cancerhotspots_bin'])
        self.assertEqual(group1.loc[0, 'cancerhotspots3D_bin'], expected_index_0_group1['cancerhotspots3D_bin'])
        self.assertEqual(group1.loc[0, 'cgc_bin'], expected_index_0_group1['cgc_bin'])
        self.assertEqual(group1.loc[0, 'cosmic_bin'], expected_index_0_group1['cosmic_bin'])
        self.assertEqual(group1.loc[0, 'gsea_pathways_bin'], expected_index_0_group1['gsea_pathways_bin'])
        self.assertEqual(group1.loc[0, 'gsea_modules_bin'], expected_index_0_group1['gsea_modules_bin'])

        expected_index_2_group1 = {
            'feature_match_1': 1,
            'feature_match_2': 1,
            'feature_match_3': 1,
            'feature_match_4': 0,
            'feature_match': 2,
            'evidence': 'FDA-Approved',
            'cancerhotspots_bin': 0,
            'cancerhotspots3D_bin': 0,
            'cgc_bin': 1,
            'cosmic_bin': 1,
            'gsea_pathways_bin': 1,
            'gsea_modules_bin': 1
        }
        self.assertEqual(group1.loc[2, 'feature_match_1'], expected_index_2_group1['feature_match_1'])
        self.assertEqual(group1.loc[2, 'feature_match_2'], expected_index_2_group1['feature_match_2'])
        self.assertEqual(group1.loc[2, 'feature_match_3'], expected_index_2_group1['feature_match_3'])
        self.assertEqual(group1.loc[2, 'feature_match_4'], expected_index_2_group1['feature_match_4'])
        self.assertEqual(group1.loc[2, 'evidence'], expected_index_2_group1['evidence'])
        self.assertEqual(group1.loc[2, 'cancerhotspots_bin'], expected_index_2_group1['cancerhotspots_bin'])
        self.assertEqual(group1.loc[2, 'cancerhotspots3D_bin'], expected_index_2_group1['cancerhotspots3D_bin'])
        self.assertEqual(group1.loc[2, 'cgc_bin'], expected_index_2_group1['cgc_bin'])
        self.assertEqual(group1.loc[2, 'cosmic_bin'], expected_index_2_group1['cosmic_bin'])
        self.assertEqual(group1.loc[2, 'gsea_pathways_bin'], expected_index_2_group1['gsea_pathways_bin'])
        self.assertEqual(group1.loc[2, 'gsea_modules_bin'], expected_index_2_group1['gsea_modules_bin'])

        expected_index_3_group1 = {
            'feature_match_1': 1,
            'feature_match_2': 0,
            'feature_match_3': 0,
            'feature_match_4': 0,
            'feature_match': 1,
            'cancerhotspots_bin': 1,
            'cancerhotspots3D_bin': 0,
            'cgc_bin': 1,
            'cosmic_bin': 1,
            'gsea_pathways_bin': 1,
            'gsea_modules_bin': 1
        }
        self.assertEqual(group1.loc[3, 'feature_match_1'], expected_index_3_group1['feature_match_1'])
        self.assertEqual(group1.loc[3, 'feature_match_2'], expected_index_3_group1['feature_match_2'])
        self.assertEqual(group1.loc[3, 'feature_match_3'], expected_index_3_group1['feature_match_3'])
        self.assertEqual(group1.loc[3, 'feature_match_4'], expected_index_3_group1['feature_match_4'])
        self.assertEqual(group1.loc[3, 'cancerhotspots_bin'], expected_index_3_group1['cancerhotspots_bin'])
        self.assertEqual(group1.loc[3, 'cancerhotspots3D_bin'], expected_index_3_group1['cancerhotspots3D_bin'])
        self.assertEqual(group1.loc[3, 'cgc_bin'], expected_index_3_group1['cgc_bin'])
        self.assertEqual(group1.loc[3, 'cosmic_bin'], expected_index_3_group1['cosmic_bin'])
        self.assertEqual(group1.loc[3, 'gsea_pathways_bin'], expected_index_3_group1['gsea_pathways_bin'])
        self.assertEqual(group1.loc[3, 'gsea_modules_bin'], expected_index_3_group1['gsea_modules_bin'])

        expected_index_0_group2 = {
            'feature_match_1': 1,
            'feature_match_2': 0,
            'feature_match_3': 0,
            'feature_match_4': 0,
            'feature_match': 3,
            'cancerhotspots_bin': 0,
            'cancerhotspots3D_bin': 0,
            'cgc_bin': 1,
            'cosmic_bin': 1,
            'gsea_pathways_bin': 1,
            'gsea_modules_bin': 0
        }
        self.assertEqual(group2.loc[0, 'feature_match_1'], expected_index_0_group2['feature_match_1'])
        self.assertEqual(group2.loc[0, 'feature_match_2'], expected_index_0_group2['feature_match_2'])
        self.assertEqual(group2.loc[0, 'feature_match_3'], expected_index_0_group2['feature_match_3'])
        self.assertEqual(group2.loc[0, 'feature_match_4'], expected_index_0_group2['feature_match_4'])
        self.assertEqual(group2.loc[0, 'cancerhotspots_bin'], expected_index_0_group2['cancerhotspots_bin'])
        self.assertEqual(group2.loc[0, 'cancerhotspots3D_bin'], expected_index_0_group2['cancerhotspots3D_bin'])
        self.assertEqual(group2.loc[0, 'cgc_bin'], expected_index_0_group2['cgc_bin'])
        self.assertEqual(group2.loc[0, 'cosmic_bin'], expected_index_0_group2['cosmic_bin'])
        self.assertEqual(group2.loc[0, 'gsea_pathways_bin'], expected_index_0_group2['gsea_pathways_bin'])
        self.assertEqual(group2.loc[0, 'gsea_modules_bin'], expected_index_0_group2['gsea_modules_bin'])

        expected_index_1_group2 = {
            'feature_match_1': 0,
            'feature_match_2': 0,
            'feature_match_3': 0,
            'feature_match_4': 0,
            'feature_match': 0,
            'evidence': '',
            'cancerhotspots_bin': 0,
            'cancerhotspots3D_bin': 0,
            'cgc_bin': 0,
            'cosmic_bin': 0,
            'gsea_pathways_bin': 0,
            'gsea_modules_bin': 0
        }
        self.assertEqual(group2.loc[1, 'feature_match_1'], expected_index_1_group2['feature_match_1'])
        self.assertEqual(group2.loc[1, 'feature_match_2'], expected_index_1_group2['feature_match_2'])
        self.assertEqual(group2.loc[1, 'feature_match_3'], expected_index_1_group2['feature_match_3'])
        self.assertEqual(group2.loc[1, 'feature_match_4'], expected_index_1_group2['feature_match_4'])
        self.assertEqual(group2.loc[1, 'evidence'], expected_index_1_group2['evidence'])
        self.assertEqual(group2.loc[1, 'cancerhotspots_bin'], expected_index_1_group2['cancerhotspots_bin'])
        self.assertEqual(group2.loc[1, 'cancerhotspots3D_bin'], expected_index_1_group2['cancerhotspots3D_bin'])
        self.assertEqual(group2.loc[1, 'cgc_bin'], expected_index_1_group2['cgc_bin'])
        self.assertEqual(group2.loc[1, 'cosmic_bin'], expected_index_1_group2['cosmic_bin'])
        self.assertEqual(group2.loc[1, 'gsea_pathways_bin'], expected_index_1_group2['gsea_pathways_bin'])
        self.assertEqual(group2.loc[1, 'gsea_modules_bin'], expected_index_1_group2['gsea_modules_bin'])

    def test_annotate_fusions_matching(self):
        dbs = {
            'almanac_handle': '../datasources/moalmanac/molecular-oncology-almanac.json',
            'cancerhotspots_handle': '../datasources/cancerhotspots/hotspots_v2.txt',
            '3dcancerhotspots_handle': '../datasources/cancerhotspots/hotspots3d.txt',
            'cgc_handle': '../datasources/cancergenecensus/cancer_gene_census_v97.genes.tsv',
            'cosmic_handle': '../datasources/cosmic/CosmicMutantExport_v97.lite.txt',
            'gsea_pathways_handle': '../datasources/gsea_gene_sets/GSEA_cancer_gene_sets.txt',
            'gsea_modules_handle': '../datasources/gsea_gene_sets/c4.cm.v6.0.symbols.txt',
            'exac_handle': '../datasources/exac/exac.expanded.r1.txt',
            'acmg_handle': '../datasources/acmg/acmg.secondaryfindings.v3.txt',
            'clinvar_handle': '../datasources/clinvar/variant_summary.lite.txt',
            'hereditary_handle': '../datasources/hereditary/hereditary.txt',
            'oncotree_handle': '../datasources/oncotree/oncotree.2023-03-09.txt',
            'lawrence_handle': '../datasources/lawrence/lawrence_mapped_ontology.txt'
        }

        feature = PreclinicalMatchmaking.feature
        feature_type = PreclinicalMatchmaking.feature_type
        alteration_type = PreclinicalMatchmaking.alteration_type
        partner = PreclinicalMatchmaking.partner
        fusion = 'Rearrangement'
        model_id = PreclinicalMatchmaking.model_id
        evidence_map_str = PreclinicalMatchmaking.evidence_map_str
        merged = PreclinicalMatchmaking.merged

        gene1 = PreclinicalMatchmaking.gene1
        gene2 = PreclinicalMatchmaking.gene2

        match = PreclinicalMatchmaking.feature_match
        match_1 = PreclinicalMatchmaking.match_1
        match_2 = PreclinicalMatchmaking.match_2
        match_3 = PreclinicalMatchmaking.match_3
        match_4 = PreclinicalMatchmaking.match_4

        df = pd.DataFrame({
            feature: ['BCR', 'BCR', 'NTRK1', 'CDKN2A', ''],
            partner: ['ABL1', '', 'ABL1', '', '']})
        df[feature_type] = fusion
        df[alteration_type] = 'Fusion'
        df[model_id] = 'case'

        almanac = datasource_Almanac.import_ds(dbs)
        db = Almanac.subset_records(almanac['content'], Almanac.feature_type, fusion)
        db = pd.DataFrame(db)
        db = db.rename(columns={gene1: feature, gene2: partner})
        db[merged] = 1
        db[evidence_map_str] = 5
        db[alteration_type] = 'Fusion'
        db_genes = ['ABL1', 'BCR', 'NTRK1', 'CDKN2A']
        result = PreclinicalMatchmaking.annotate_fusions_matching(df, db, db_genes, consider_partner=True).fillna(0)

        expected_0 = {match_1: 1, match_2: 1, match_3: 1, match_4: 1, match: 4}
        expected_1 = {match_1: 1, match_2: 1, match_3: 1, match_4: 0, match: 3}
        expected_2 = {match_1: 1, match_2: 1, match_3: 1, match_4: 0, match: 3}
        expected_3 = {match_1: 1, match_2: 0, match_3: 0, match_4: 0, match: 1}
        expected_4 = {match_1: 0, match_2: 0, match_3: 0, match_4: 0, match: 0}

        for index, dictionary in [(0, expected_0), (1, expected_1), (2, expected_2), (3, expected_3), (4, expected_4)]:
            for key, value in dictionary.items():
                self.assertEqual(result.loc[index, key], value)

    def test_annotate_somatic_variants(self):
        dbs = {
            'almanac_handle': '../datasources/moalmanac/molecular-oncology-almanac.json',
            'cancerhotspots_handle': '../datasources/cancerhotspots/hotspots_v2.txt',
            '3dcancerhotspots_handle': '../datasources/cancerhotspots/hotspots3d.txt',
            'cgc_handle': '../datasources/cancergenecensus/cancer_gene_census_v97.genes.tsv',
            'cosmic_handle': '../datasources/cosmic/CosmicMutantExport_v97.lite.txt',
            'gsea_pathways_handle': '../datasources/gsea_gene_sets/GSEA_cancer_gene_sets.txt',
            'gsea_modules_handle': '../datasources/gsea_gene_sets/c4.cm.v6.0.symbols.txt',
            'exac_handle': '../datasources/exac/exac.expanded.r1.txt',
            'acmg_handle': '../datasources/acmg/acmg.secondaryfindings.v3.txt',
            'clinvar_handle': '../datasources/clinvar/variant_summary.lite.txt',
            'hereditary_handle': '../datasources/hereditary/hereditary.txt',
            'oncotree_handle': '../datasources/oncotree/oncotree.2023-03-09.txt',
            'lawrence_handle': '../datasources/lawrence/lawrence_mapped_ontology.txt'
        }
        feature = PreclinicalMatchmaking.feature
        feature_type = PreclinicalMatchmaking.feature_type
        alteration_type = PreclinicalMatchmaking.alteration_type
        alteration = PreclinicalMatchmaking.alteration

        somatic_variant = 'Somatic Variant'

        df = pd.DataFrame({
            feature: ['BRAF', 'BRAF', 'IDH1', 'CDKN2A'],
            alteration_type: ['Missense', 'Missense', 'Nonsense', 'Nonsense'],
            alteration: ['p.V600E', 'p.F00N', 'p.F00*', 'p.F42*']
        })
        df[feature_type] = somatic_variant

        result = PreclinicalMatchmaking.annotate_somatic_variants(df, dbs, biomarker_type_string=somatic_variant)

        expected_braf_1 = {
            'feature_match_1': 1,
            'feature_match_2': 1,
            'feature_match_3': 1,
            'feature_match_4': 1,
            'feature_match': 0,
            'evidence': 5,
            'cancerhotspots_bin': 2,
            'cancerhotspots3D_bin': 2,
            'cgc_bin': 1,
            'cosmic_bin': 2,
            'gsea_pathways_bin': 1,
            'gsea_modules_bin': 1
        }

        expected_braf_2 = {
            'feature_match_1': 1,
            'feature_match_2': 1,
            'feature_match_3': 1,
            'feature_match_4': 0,
            'feature_match': 0,
            'evidence': 5,
            'cancerhotspots_bin': 1,
            'cancerhotspots3D_bin': 1,
            'cgc_bin': 1,
            'cosmic_bin': 1,
            'gsea_pathways_bin': 1,
            'gsea_modules_bin': 1
        }

        expected_idh1 = {
            'feature_match_1': 1,
            'feature_match_2': 1,
            'feature_match_3': 0,
            'feature_match_4': 0,
            'feature_match': 0,
            'evidence': 5,
            'cancerhotspots_bin': 1,
            'cancerhotspots3D_bin': 0,
            'cgc_bin': 1,
            'cosmic_bin': 1,
            'gsea_pathways_bin': 0,
            'gsea_modules_bin': 1
        }

        expected_cdkn2a = {
            'feature_match_1': 1,
            'feature_match_2': 0,
            'feature_match_3': 0,
            'feature_match_4': 0,
            'feature_match': 0,
            'cancerhotspots_bin': 1,
            'cancerhotspots3D_bin': 0,
            'cgc_bin': 1,
            'cosmic_bin': 1,
            'gsea_pathways_bin': 1,
            'gsea_modules_bin': 1
        }

        for key, value in expected_braf_1.items():
            self.assertEqual(result.loc[0, key], expected_braf_1[key])
        for key, value in expected_braf_2.items():
            self.assertEqual(result.loc[1, key], expected_braf_2[key])
        for key, value in expected_idh1.items():
            self.assertEqual(result.loc[2, key], expected_idh1[key])
        for key, value in expected_cdkn2a.items():
            self.assertEqual(result.loc[3, key], expected_cdkn2a[key])

    def test_annotate_match_1(self):
        feature = PreclinicalMatchmaking.feature
        match_1 = PreclinicalMatchmaking.match_1

        df = pd.DataFrame({
            feature: ['Gene1', 'Gene2', '', ''],
            'expectation': [1, 0, 1, 1]
        })
        genes = ['Gene1', '']
        result = PreclinicalMatchmaking.annotate_match_1(df, genes).fillna(0)
        for index in result.index:
            self.assertEqual(result.loc[index, match_1], result.loc[index, 'expectation'])

    def test_annotate_match_2(self):
        match_2 = PreclinicalMatchmaking.match_2

        dbs = {
            'almanac_handle': '../datasources/moalmanac/molecular-oncology-almanac.json',
            'cancerhotspots_handle': '../datasources/cancerhotspots/hotspots_v2.txt',
            '3dcancerhotspots_handle': '../datasources/cancerhotspots/hotspots3d.txt',
            'cgc_handle': '../datasources/cancergenecensus/cancer_gene_census_v97.genes.tsv',
            'cosmic_handle': '../datasources/cosmic/CosmicMutantExport_v97.lite.txt',
            'gsea_pathways_handle': '../datasources/gsea_gene_sets/GSEA_cancer_gene_sets.txt',
            'gsea_modules_handle': '../datasources/gsea_gene_sets/c4.cm.v6.0.symbols.txt',
            'exac_handle': '../datasources/exac/exac.expanded.r1.txt',
            'acmg_handle': '../datasources/acmg/acmg.secondaryfindings.v3.txt',
            'clinvar_handle': '../datasources/clinvar/variant_summary.lite.txt',
            'hereditary_handle': '../datasources/hereditary/hereditary.txt',
            'oncotree_handle': '../datasources/oncotree/oncotree.2023-03-09.txt',
            'lawrence_handle': '../datasources/lawrence/lawrence_mapped_ontology.txt'
        }
        almanac = datasource_Almanac.import_ds(dbs)
        copy_number = 'Copy Number'
        fusion = 'Rearrangement'
        somatic_variant = 'Somatic Variant'
        
        feature = PreclinicalMatchmaking.feature
        alteration_type = PreclinicalMatchmaking.alteration_type
        alteration = PreclinicalMatchmaking.alteration
        gene = PreclinicalMatchmaking.gene
        variant_annotation = PreclinicalMatchmaking.variant_annotation
        protein_change = PreclinicalMatchmaking.protein_change
        direction = PreclinicalMatchmaking.direction
        rearrangement_type = PreclinicalMatchmaking.rearrangement_type
        gene1 = PreclinicalMatchmaking.gene1
        gene2 = PreclinicalMatchmaking.gene2

        df = pd.DataFrame({feature: ['BRAF', 'CDKN2A', ''], 'expectation': [1, 0, 0]})
        column_map = {gene: feature,
                      variant_annotation: alteration_type,
                      protein_change: alteration}
        db = Almanac.subset_records(almanac['content'], Almanac.feature_type, somatic_variant)
        db = pd.DataFrame(db)
        db = PreclinicalMatchmaking.format_db(list(column_map.keys()), column_map, db)

        result = PreclinicalMatchmaking.annotate_match_2(df, db, feature=feature).fillna(0)
        self.assertEqual(result.loc[0, match_2], result.loc[0, 'expectation'])
        self.assertEqual(result.loc[1, match_2], result.loc[1, 'expectation'])
        self.assertEqual(result.loc[2, match_2], result.loc[2, 'expectation'])

        df = pd.DataFrame({feature: ['ERBB2', 'KRAS', ''], 'expectation': [1, 0, 0]})
        column_map = {gene: feature, direction: alteration_type}
        db = Almanac.subset_records(almanac['content'], Almanac.feature_type, copy_number)
        db = pd.DataFrame(db)
        db = PreclinicalMatchmaking.format_db(list(column_map.keys()), column_map, db)

        result = PreclinicalMatchmaking.annotate_match_2(df, db, feature=feature).fillna(0)
        self.assertEqual(result.loc[0, match_2], result.loc[0, 'expectation'])
        self.assertEqual(result.loc[1, match_2], result.loc[1, 'expectation'])
        self.assertEqual(result.loc[2, match_2], result.loc[2, 'expectation'])

        df = pd.DataFrame({'gene1': ['ALK', 'ABL1', 'BCR', ''],
                           'gene2': ['ALK', 'ABL1', 'BCR', ''],
                           'expectation1': [1, 0, 1, 0],
                           'expectation2': [1, 1, 0, 0]})
        column_map = {rearrangement_type: alteration_type}
        db = Almanac.subset_records(almanac['content'], Almanac.feature_type, fusion)
        db = pd.DataFrame(db)
        db = PreclinicalMatchmaking.format_db([gene1, gene2, rearrangement_type], column_map, db)
        result = PreclinicalMatchmaking.annotate_match_2(df, db, feature='gene1').fillna(0)
        self.assertEqual(result.loc[0, match_2], result.loc[0, 'expectation1'])
        self.assertEqual(result.loc[1, match_2], result.loc[1, 'expectation1'])
        self.assertEqual(result.loc[2, match_2], result.loc[2, 'expectation1'])
        self.assertEqual(result.loc[3, match_2], result.loc[3, 'expectation1'])

        result = PreclinicalMatchmaking.annotate_match_2(df, db, feature='gene2').fillna(0)
        self.assertEqual(result.loc[0, match_2], result.loc[0, 'expectation2'])
        self.assertEqual(result.loc[1, match_2], result.loc[1, 'expectation2'])
        self.assertEqual(result.loc[2, match_2], result.loc[2, 'expectation2'])
        self.assertEqual(result.loc[3, match_2], result.loc[3, 'expectation2'])

    def test_annotate_match_3(self):
        match_3 = PreclinicalMatchmaking.match_3

        dbs = {
            'almanac_handle': '../datasources/moalmanac/molecular-oncology-almanac.json',
            'cancerhotspots_handle': '../datasources/cancerhotspots/hotspots_v2.txt',
            '3dcancerhotspots_handle': '../datasources/cancerhotspots/hotspots3d.txt',
            'cgc_handle': '../datasources/cancergenecensus/cancer_gene_census_v97.genes.tsv',
            'cosmic_handle': '../datasources/cosmic/CosmicMutantExport_v97.lite.txt',
            'gsea_pathways_handle': '../datasources/gsea_gene_sets/GSEA_cancer_gene_sets.txt',
            'gsea_modules_handle': '../datasources/gsea_gene_sets/c4.cm.v6.0.symbols.txt',
            'exac_handle': '../datasources/exac/exac.expanded.r1.txt',
            'acmg_handle': '../datasources/acmg/acmg.secondaryfindings.v3.txt',
            'clinvar_handle': '../datasources/clinvar/variant_summary.lite.txt',
            'hereditary_handle': '../datasources/hereditary/hereditary.txt',
            'oncotree_handle': '../datasources/oncotree/oncotree.2023-03-09.txt',
            'lawrence_handle': '../datasources/lawrence/lawrence_mapped_ontology.txt'
        }
        almanac = datasource_Almanac.import_ds(dbs)
        copy_number = 'Copy Number'
        fusion = 'Rearrangement'
        somatic_variant = 'Somatic Variant'

        feature = PreclinicalMatchmaking.feature
        alteration_type = PreclinicalMatchmaking.alteration_type
        alteration = PreclinicalMatchmaking.alteration
        gene = PreclinicalMatchmaking.gene
        variant_annotation = PreclinicalMatchmaking.variant_annotation
        protein_change = PreclinicalMatchmaking.protein_change
        direction = PreclinicalMatchmaking.direction
        rearrangement_type = PreclinicalMatchmaking.rearrangement_type
        gene1 = PreclinicalMatchmaking.gene1
        gene2 = PreclinicalMatchmaking.gene2

        df = pd.DataFrame({feature: ['BRAF', 'CDKN2A', 'BRAF'],
                           alteration_type: ['Missense', 'Missense', 'Nonsense'],
                           'expectation': [1, 0, 0]})
        column_map = {gene: feature,
                      variant_annotation: alteration_type,
                      protein_change: alteration}
        db = Almanac.subset_records(almanac['content'], Almanac.feature_type, somatic_variant)
        db = pd.DataFrame(db)
        db = PreclinicalMatchmaking.format_db(list(column_map.keys()), column_map, db)

        result = PreclinicalMatchmaking.annotate_match_3(df, db,
                                                         feature=feature,
                                                         alteration_type=alteration_type).fillna(0)
        for index in result.index:
            self.assertEqual(result.loc[index, match_3], result.loc[index, 'expectation'])

        df = pd.DataFrame({feature: ['ERBB2', 'KRAS', 'CDKN2A'],
                           alteration_type: ['Amplification', 'Deletion', 'Deletion'],
                           'expectation': [1, 0, 1]})
        column_map = {gene: feature, direction: alteration_type}
        db = Almanac.subset_records(almanac['content'], Almanac.feature_type, copy_number)
        db = pd.DataFrame(db)
        db = PreclinicalMatchmaking.format_db(list(column_map.keys()), column_map, db)

        result = PreclinicalMatchmaking.annotate_match_3(df, db,
                                                         feature=feature,
                                                         alteration_type=alteration_type).fillna(0)
        for index in result.index:
            self.assertEqual(result.loc[index, match_3], result.loc[index, 'expectation'])

        df = pd.DataFrame({gene1: ['ALK', 'ABL1', 'BCR', ''],
                           gene2: ['ALK', 'ABL1', 'BCR', ''],
                           alteration_type: ['Fusion', 'Fusion', 'Fusion', 'Fusion'],
                           'expectation1': [1, 0, 1, 0],
                           'expectation2': [1, 1, 0, 0]})
        column_map = {rearrangement_type: alteration_type}
        db = Almanac.subset_records(almanac['content'], Almanac.feature_type, fusion)
        db = pd.DataFrame(db)
        db = PreclinicalMatchmaking.format_db([gene1, gene2, rearrangement_type], column_map, db)
        result = PreclinicalMatchmaking.annotate_match_3(df, db,
                                                         feature=gene1,
                                                         alteration_type=alteration_type).fillna(0)
        self.assertEqual(result.loc[0, match_3], result.loc[0, 'expectation1'])
        self.assertEqual(result.loc[1, match_3], result.loc[1, 'expectation1'])
        self.assertEqual(result.loc[2, match_3], result.loc[2, 'expectation1'])

        result = PreclinicalMatchmaking.annotate_match_3(df, db,
                                                         feature=gene2,
                                                         alteration_type=alteration_type).fillna(0)
        self.assertEqual(result.loc[0, match_3], result.loc[0, 'expectation2'])
        self.assertEqual(result.loc[1, match_3], result.loc[1, 'expectation2'])
        self.assertEqual(result.loc[2, match_3], result.loc[2, 'expectation2'])

    def test_annotate_match_4(self):
        match_4 = PreclinicalMatchmaking.match_4

        dbs = {
            'almanac_handle': '../datasources/moalmanac/molecular-oncology-almanac.json',
            'cancerhotspots_handle': '../datasources/cancerhotspots/hotspots_v2.txt',
            '3dcancerhotspots_handle': '../datasources/cancerhotspots/hotspots3d.txt',
            'cgc_handle': '../datasources/cancergenecensus/cancer_gene_census_v97.genes.tsv',
            'cosmic_handle': '../datasources/cosmic/CosmicMutantExport_v97.lite.txt',
            'gsea_pathways_handle': '../datasources/gsea_gene_sets/GSEA_cancer_gene_sets.txt',
            'gsea_modules_handle': '../datasources/gsea_gene_sets/c4.cm.v6.0.symbols.txt',
            'exac_handle': '../datasources/exac/exac.expanded.r1.txt',
            'acmg_handle': '../datasources/acmg/acmg.secondaryfindings.v3.txt',
            'clinvar_handle': '../datasources/clinvar/variant_summary.lite.txt',
            'hereditary_handle': '../datasources/hereditary/hereditary.txt',
            'oncotree_handle': '../datasources/oncotree/oncotree.2023-03-09.txt',
            'lawrence_handle': '../datasources/lawrence/lawrence_mapped_ontology.txt'
        }
        almanac = datasource_Almanac.import_ds(dbs)
        copy_number = 'Copy Number'
        fusion = 'Rearrangement'
        somatic_variant = 'Somatic Variant'

        feature = PreclinicalMatchmaking.feature
        alteration_type = PreclinicalMatchmaking.alteration_type
        alteration = PreclinicalMatchmaking.alteration
        gene = PreclinicalMatchmaking.gene
        variant_annotation = PreclinicalMatchmaking.variant_annotation
        protein_change = PreclinicalMatchmaking.protein_change
        rearrangement_type = PreclinicalMatchmaking.rearrangement_type
        gene1 = PreclinicalMatchmaking.gene1
        gene2 = PreclinicalMatchmaking.gene2

        df = pd.DataFrame({feature: ['BRAF', 'CDKN2A', 'BRAF'],
                           alteration_type: ['Missense', 'Missense', 'Nonsense'],
                           alteration: ['p.V600E', 'p.F00N', 'p.F00*'],
                           'expectation': [1, 0, 0]})
        column_map = {gene: feature,
                      variant_annotation: alteration_type,
                      protein_change: alteration}
        db = Almanac.subset_records(almanac['content'], Almanac.feature_type, somatic_variant)
        db = pd.DataFrame(db)
        db = PreclinicalMatchmaking.format_db(list(column_map.keys()), column_map, db)

        result = PreclinicalMatchmaking.annotate_match_4(df, db).fillna(0)
        self.assertEqual(result.loc[0, match_4], result.loc[0, 'expectation'])
        self.assertEqual(result.loc[1, match_4], result.loc[1, 'expectation'])
        self.assertEqual(result.loc[2, match_4], result.loc[2, 'expectation'])

        df = pd.DataFrame({gene1: ['BCR', 'EML4', 'ALK', ''],
                           gene2: ['ABL1', 'ALK', 'BCR', ''],
                           alteration_type: ['Fusion', 'Fusion', 'Fusion', 'Fusion'],
                           'expectation1': [1, 1, 0, 0]})
        column_map = {rearrangement_type: alteration_type}
        db = Almanac.subset_records(almanac['content'], Almanac.feature_type, fusion)
        db = pd.DataFrame(db)
        db = PreclinicalMatchmaking.format_db([gene1, gene2, rearrangement_type], column_map, db)
        result = PreclinicalMatchmaking.annotate_match_4(df, db,
                                                         feature=gene1,
                                                         alteration=gene2).fillna(0)

        self.assertEqual(result.loc[0, match_4], result.loc[0, 'expectation1'])
        self.assertEqual(result.loc[1, match_4], result.loc[1, 'expectation1'])
        self.assertEqual(result.loc[2, match_4], result.loc[2, 'expectation1'])
        self.assertEqual(result.loc[3, match_4], result.loc[3, 'expectation1'])

    def test_format_db(self):
        dbs = {
            'almanac_handle': '../datasources/moalmanac/molecular-oncology-almanac.json',
            'cancerhotspots_handle': '../datasources/cancerhotspots/hotspots_v2.txt',
            '3dcancerhotspots_handle': '../datasources/cancerhotspots/hotspots3d.txt',
            'cgc_handle': '../datasources/cancergenecensus/cancer_gene_census_v97.genes.tsv',
            'cosmic_handle': '../datasources/cosmic/CosmicMutantExport_v97.lite.txt',
            'gsea_pathways_handle': '../datasources/gsea_gene_sets/GSEA_cancer_gene_sets.txt',
            'gsea_modules_handle': '../datasources/gsea_gene_sets/c4.cm.v6.0.symbols.txt',
            'exac_handle': '../datasources/exac/exac.expanded.r1.txt',
            'acmg_handle': '../datasources/acmg/acmg.secondaryfindings.v3.txt',
            'clinvar_handle': '../datasources/clinvar/variant_summary.lite.txt',
            'hereditary_handle': '../datasources/hereditary/hereditary.txt',
            'oncotree_handle': '../datasources/oncotree/oncotree.2023-03-09.txt',
            'lawrence_handle': '../datasources/lawrence/lawrence_mapped_ontology.txt'
        }
        almanac = datasource_Almanac.import_ds(dbs)
        copy_number = 'Copy Number'
        fusion = 'Rearrangement'
        somatic_variant = 'Somatic Variant'

        feature = PreclinicalMatchmaking.feature
        alteration_type = PreclinicalMatchmaking.alteration_type
        alteration = PreclinicalMatchmaking.alteration
        gene = PreclinicalMatchmaking.gene
        variant_annotation = PreclinicalMatchmaking.variant_annotation
        protein_change = PreclinicalMatchmaking.protein_change
        direction = PreclinicalMatchmaking.direction
        rearrangement_type = PreclinicalMatchmaking.rearrangement_type
        gene1 = PreclinicalMatchmaking.gene1
        gene2 = PreclinicalMatchmaking.gene2

        evidence_str = PreclinicalMatchmaking.evidence
        evidence_map_str = PreclinicalMatchmaking.evidence_map_str
        evidence_map = PreclinicalMatchmaking.evidence_map

        column_map = {gene: feature,
                      variant_annotation: alteration_type,
                      protein_change: alteration}
        db = Almanac.subset_records(almanac['content'], Almanac.feature_type, somatic_variant)
        db = pd.DataFrame(db)
        result = PreclinicalMatchmaking.format_db(list(column_map.keys()), column_map, db)
        for index in result.index:
            self.assertEqual(result.loc[index, evidence_map_str], evidence_map[result.loc[index, evidence_str]])

        column_map = {gene: feature, direction: alteration_type}
        db = Almanac.subset_records(almanac['content'], Almanac.feature_type, copy_number)
        db = pd.DataFrame(db)
        result = PreclinicalMatchmaking.format_db(list(column_map.keys()), column_map, db)
        for index in result.index:
            self.assertEqual(result.loc[index, evidence_map_str], evidence_map[result.loc[index, evidence_str]])

        column_map = {rearrangement_type: alteration_type}
        db = Almanac.subset_records(almanac['content'], Almanac.feature_type, fusion)
        db = pd.DataFrame(db)
        result = PreclinicalMatchmaking.format_db([gene1, gene2, rearrangement_type], column_map, db).fillna('')
        for index in result.index:
            self.assertEqual(result.loc[index, evidence_map_str], evidence_map[result.loc[index, evidence_str]])

    def test_map_evidence(self):
        evidence_str = PreclinicalMatchmaking.evidence
        evidence_map_str = PreclinicalMatchmaking.evidence_map_str
        evidence_map = {'FDA-Approved': 5.0, 'Guideline': 4.0, 'Clinical trial': 3.0,
                        'Clinical evidence': 2.0, 'Preclinical': 1.0, 'Inferential': 0.0}

        dataframe = pd.Series(list(evidence_map.values()), name=evidence_str).to_frame()
        result = PreclinicalMatchmaking.map_evidence(dataframe).set_index(evidence_map_str)
        for label, value in evidence_map.items():
            self.assertEqual(result.loc[value, evidence_str], label)

    def test_preallocate_almanac_matches(self):
        empty = pd.DataFrame(index=[0])
        result = PreclinicalMatchmaking.preallocate_almanac_matches(empty)
        columns = [PreclinicalMatchmaking.match_1, PreclinicalMatchmaking.match_2, PreclinicalMatchmaking.match_3,
                   PreclinicalMatchmaking.match_4, PreclinicalMatchmaking.feature_match]
        for column in columns:
            self.assertEqual(result.loc[0, column], 0)


if __name__ == '__main__':
    unittest.main()
