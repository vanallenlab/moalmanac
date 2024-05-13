import argparse
import pandas as pd

FILES = [
    "example.actionable.txt",
    "example.germline.acmg.txt",
    "example.germline.cancer_related.txt",
    "example.germline.hereditary_cancers.txt",
    "example.integrated.summary.txt",
    "example.matchmaker.txt",
    "example.msi_variants.txt",
    "example.mutational_burden.txt",
    "example.preclinical.efficacy.txt",
    "example.somatic.filtered.txt",
    "example.somatic.scored.txt",
    "example.therapeutic_strategies.txt"
]


def main(output_folder_1, output_folder_2):
    print(f"Comparing moalmanac outputs from {output_folder_1} and {output_folder_2}")
    print("")

    if output_folder_1[-1] == "/":
        output_folder_1 = output_folder_1[:-1]
    if output_folder_2[-1] == "/":
        output_folder_2 = output_folder_2[:-1]

    for file in FILES:
        file1 = f"{output_folder_1}/{file}"
        file2 = f"{output_folder_2}/{file}"
        print("***")
        print(f'{file}')
        print('')

        df1 = pd.read_csv(file1, sep='\t')
        df2 = pd.read_csv(file2, sep='\t')
        df1_row_hashes = pd.util.hash_pandas_object(df1, index=True)
        df2_row_hashes = pd.util.hash_pandas_object(df2, index=True)

        if df1.shape != df2.shape:
            print(f"Shape differs: {df1.shape} vs. {df2.shape}")
            print("")

        row_hashes_match = df1_row_hashes.eq(df2_row_hashes)
        if False not in row_hashes_match.tolist():
            print(f"matches")
            print('')
            continue

        non_matching_rows = df1_row_hashes[~row_hashes_match]
        for row_index in non_matching_rows.index:
            series_1_hash = pd.util.hash_pandas_object(df1.loc[row_index, :], index=True)
            series_2_hash = pd.util.hash_pandas_object(df2.loc[row_index, :], index=True)

            column_hashes_match = series_1_hash.eq(series_2_hash)
            if False not in column_hashes_match.tolist():
                print(f"{file1} index {row_index} matches {file2} {row_index}")
                print('')
                continue

            non_matching_columns = series_1_hash[~column_hashes_match]
            for column_name in non_matching_columns.index:
                df1_value = df1.loc[row_index, column_name]
                df2_value = df2.loc[row_index, column_name]
                print(f"row {row_index} {column_name} do not match.")
                print(f"{df1_value} vs {df2_value}")
                print('')


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(prog='Molecular Oncology Almanac output regression test',
                                         description='Hash and compare output files from two different folders')
    arg_parser.add_argument('--folder1', '-f1',
                            type=str,
                            help='First folder of output files',
                            required=True)
    arg_parser.add_argument('--folder2', '-f2',
                            type=str,
                            help='Second folder of output files',
                            required=True)
    args = arg_parser.parse_args()

    main(args.folder1, args.folder2)
