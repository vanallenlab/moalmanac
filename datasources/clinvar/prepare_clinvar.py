import pandas
import argparse

COLUMNS = [
    'GeneSymbol',
    'Chromosome',
    'Start',
    'Stop',
    'ReferenceAllele',
    'AlternateAllele',
    'ClinicalSignificance',
    'ClinSigSimple'
]


def create_output_filename(date):
    return f"variant_summary.reduced.{date}.txt"


def read_file(file, relevant_columns):
    return pandas.read_csv(file, sep='\t', usecols=relevant_columns, low_memory=False)


def write_file(dataframe, date):
    output_name = create_output_filename(date)
    dataframe.to_csv(output_name, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='prepare clinvar', description='Prepare ClinVar for use with MOAlmanac')
    parser.add_argument('--input', '-i', help='input file, CosmicMutantExport.tsv', required=True)
    parser.add_argument('--date', '-d', help='date of access; e.g. 2023-03-09', required=True)
    args = parser.parse_args()

    df = read_file(args.input, COLUMNS)
    df.drop_duplicates(inplace=True)
    write_file(df, args.date)

    gene_count = df['GeneSymbol'].drop_duplicates().shape[0]
    total_count = df.shape[0]
    print(f"As of {args.date}, ClinVar contains {gene_count} genes and {total_count} variants.")
