import pandas
import argparse


def read_file(file, relevant_columns):
	return pandas.read_csv(file, sep='\t', usecols=relevant_columns)


def write_file(dataframe, version):
	output_name = f"CosmicMutantExport_{version}.lite.txt"
	dataframe.to_csv(output_name, sep='\t', index=False)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='compare outputs', description='Compare changed columns relative to unchanged columns')
	parser.add_argument('--file1', '-f1', help='input file 1', required=True)
	parser.add_argument('--file2', '-f2', help='input file 2', required=True)
	parser.add_argument('--output', '-o', help='output file', required=True)
	parser.add_argument('--column_unchanged', '-u', action="append", help='columns that did not change', required=True)
	parser.add_argument('--column_changed', '-c', help='column that did change', required=True)
	args = parser.parse_args()
	print(args)

	all_columns = args.column_unchanged + [args.column_changed]
	df1 = read_file(args.file1, all_columns)
	df2 = read_file(args.file2, all_columns)
	merged = df1.merge(df2, on=args.column_unchanged, how='outer')
	print(merged.head())

	idx_difference = merged[f"{args.column_changed}_x"].ne(merged[f"{args.column_changed}_y"])
	result = (
		merged
		.loc[idx_difference, :]
		.sort_values(by=[f"{args.column_changed}_x", f"{args.column_changed}_y"], axis=0, ascending=True)
	)
	result.to_csv(args.output, sep='\t', index=False)
