import argparse
import pandas


def read_file(file, column_name):
	return pandas.read_csv(file, sep='\t', usecols=[column_name])


def write_file(dataframe, output_name):
	dataframe.to_csv(output_name, sep='\t', index=False)


if __name__ == "__main__":
	arg_parser = argparse.ArgumentParser(prog='extract genes', description='Extract genes from tab separated values.')
	arg_parser.add_argument('--input', '-i', help='input file', required=True)
	arg_parser.add_argument('--output', '-o', help='output file', required=True)
	arg_parser.add_argument('--gene_column_name', '-g', help='column which contains gene names', default="Gene Symbol")
	args = arg_parser.parse_args()

	df = read_file(args.input, args.gene_column_name)
	df.to_csv(args.output, sep='\t', index=False)
