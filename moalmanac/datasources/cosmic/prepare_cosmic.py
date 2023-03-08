import pandas
import argparse


def create_output_filename(version):
	return f"CosmicMutantExport_{version}.lite.txt"


def read_file(file, relevant_columns):
	return pandas.read_csv(file, sep='\t', usecols=relevant_columns)


def write_file(dataframe, version):
	output_name = f"CosmicMutantExport_{version}.lite.txt"
	dataframe.to_csv(output_name, sep='\t', index=False)


if __name__ == "__main__":
	arg_parser = argparse.ArgumentParser(prog='extract genes', description='Extract genes from tab separated values.')
	arg_parser.add_argument('--input', '-i', help='input file, CosmicMutantExport.tsv', required=True)
	arg_parser.add_argument('--version', '-v', help='input file version; e.g., v97', required=True)
	arg_parser.add_argument('--gene_column_name', '-g', help='column that contains gene names', default="Gene name")
	arg_parser.add_argument('--protein_column_name', '-p', help='column that contains protein changes', default="Mutation AA")
	args = arg_parser.parse_args()

	df = read_file(args.input, [args.gene_column_name, args.protein_column_name])
	df.drop_duplicates(inplace=True)
	write_file(df, args.version)

	gene_count = df[args.gene_column_name].drop_duplicates().shape[0]
	total_count = df.shape[0]
	print(f"COSMIC {args.version} contains {gene_count} genes and {total_count} protein changes")
