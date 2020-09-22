import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--clinvarVariantSummary', help='')
args = parser.parse_args()

_columns = ['GeneSymbol', 'Chromosome', 'Start', 'Stop', 'ReferenceAllele', 'AlternateAllele',
            'ClinicalSignificance', 'ClinSigSimple']

df = pd.read_csv(args.clinvarVariantSummary, sep = '\t', usecols = _columns)
output_name = str(args.clinvarVariantSummary)[:-4] + '.lite.txt'
df.to_csv(output_name, sep = '\t', index = False)
