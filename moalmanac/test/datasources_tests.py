import unittest
import pandas as pd

from datasources import CancerHotspots, Preclinical


class UnitTestCancerHotspots(unittest.TestCase):
    def test_format_cancerhotspots(self):
        alt = CancerHotspots.alt
        aa_pos = CancerHotspots.aa_pos
        aa_ref = CancerHotspots.aa_ref
        aa_var = CancerHotspots.aa_var

        initial_df = pd.DataFrame(
            {aa_ref: ['Q:100', 'V:101', 'S:5'],
             aa_pos: ['1', '2', '3'],
             aa_var: ['A:50', 'T:1', 'L:123'],
             alt: ['', '', '']})
        df = CancerHotspots.format_cancerhotspots(initial_df)
        self.assertEqual(['p.Q1A', 'p.V2T', 'p.S3L'], df.loc[:, alt].tolist())


class UnitTestPreclinical(unittest.TestCase):
    def test_import_dbs(self):
        summary = Preclinical.summary
        variants = Preclinical.variants
        cnas = Preclinical.cnas
        fusions = Preclinical.fusions
        gdsc = Preclinical.gdsc
        mappings = Preclinical.mappings
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
        dbs = Preclinical.import_dbs(dbs_dictionary)

        for label in [summary, variants, cnas, fusions, gdsc]:
            self.assertEqual(type(dbs[label]), type(pd.DataFrame()))
        for dictionary in [dbs[mappings], dbs]:
            self.assertEqual(type(dictionary), type({}))

    def test_subset_dataframe(self):
        dataframe = pd.DataFrame({
            'A': [1, 2, 3],
            'B': [3, 4, 5]
        })
        result = Preclinical.subset_dataframe(dataframe, 'B', [3, 4])
        self.assertEqual(result['A'].tolist(), [1, 2])

    def test_record_gene_hits(self):
        variants = pd.DataFrame({'feature': ['A', 'B'], 'model_id': [1, 2]})
        cnas = pd.DataFrame({'feature': ['B', 'C'], 'model_id': [2, 3]})
        fusions = pd.DataFrame({'feature': ['C', 'D'], 'model_id': [1, 3], 'partner': ['A', 'B']})

        results = Preclinical.record_gene_hits(variants, cnas, fusions)
        self.assertEqual(results['feature'].tolist(), ['A', 'B', 'C', 'C', 'D', 'B'])
        self.assertEqual(results['model_id'].tolist(), [1, 2, 3, 1, 3, 3])
