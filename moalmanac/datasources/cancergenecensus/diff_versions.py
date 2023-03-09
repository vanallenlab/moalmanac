import argparse
import pandas


def get_set_difference(case_series, comparison_series, column_name):
	new_members = case_series.difference(comparison_series)
	return pandas.Series(new_members, name=column_name)


def read_file(file, column_name):
	dataframe = pandas.read_csv(file, sep='\t', usecols=[column_name])
	return dataframe.set_index(column_name).index


def write_file(dataframe, output_name):
	dataframe.to_csv(output_name, sep='\t', index=False)


if __name__ == "__main__":
	description = "Identify genes added and removed between versions of Cancer Gene Census"
	arg_parser = argparse.ArgumentParser(prog='diff versions', description=description)
	arg_parser.add_argument('--old_version', '-o', help='input file, old version of datasource', required=True)
	arg_parser.add_argument('--new_version', '-n', help='input file, new version of datasource', required=True)
	arg_parser.add_argument('--gene_column_name', '-g', help='column which contains gene names', default="Gene Symbol")
	args = arg_parser.parse_args()

	old = read_file(args.old_version, args.gene_column_name)
	new = read_file(args.new_version, args.gene_column_name)

	removals = get_set_difference(old, new, args.gene_column_name)
	additions = get_set_difference(new, old, args.gene_column_name)

	print(f"{len(removals)} genes have been removed between {args.old_version} and {args.new_version}")
	print(f"{', '.join(removals.tolist())}")
	print('')

	print(f"{len(additions)} new genes appear in {args.new_version} that were not present in {args.old_version}")
	print(f"{', '.join(additions.tolist())}")
	print('')
