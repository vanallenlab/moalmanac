import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--cosmicMutantExport', help='CosmicMutantExport_vXY.tsv')
args = parser.parse_args()

columns_ = ['Gene name', 'Mutation AA']
df = pd.read_csv(args.cosmicMutantExport, sep = '\t', usecols = columns_)
output_name = str(args.cosmicMutantExport)[:-4] + '.lite.txt'
df.to_csv(output_name, sep = '\t', index = False)
