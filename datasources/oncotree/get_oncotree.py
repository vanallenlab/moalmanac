import argparse
import pandas
import requests
import sys


def convert_dictionary_to_dataframe(dictionary):
    return pandas.DataFrame.from_dict(dictionary, orient='index').T


def convert_response_to_df(json):
    df = convert_dictionary_to_dataframe(json[0])
    for entry in json[1:]:
        tmp = convert_dictionary_to_dataframe(entry)
        df = pandas.concat([df, tmp], ignore_index=True, axis=0)
    return df


def create_output_filename(date):
    return f"oncotree.{date}.txt"


def request_get(request_uri):
    return requests.get(request_uri)


def read_file(file, relevant_columns):
    return pandas.read_csv(file, sep='\t', usecols=relevant_columns, low_memory=False)


def subset_oncotree(dataframe):
    column_map = {
        'code': 'code',
        'mainType': 'main_type',
        'name': 'name',
        'tissue': 'tissue'
    }

    return dataframe.loc[:, column_map.keys()].sort_values('code').rename(columns=column_map)


def write_file(dataframe, date):
    output_name = create_output_filename(date)
    dataframe.to_csv(output_name, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='get oncotree', description='Download Oncotree for use with MOAlmanac')
    parser.add_argument('--date', '-d', help='date of access; e.g. 2023-03-09', required=True)
    args = parser.parse_args()

    endpoint = "http://oncotree.mskcc.org/api/tumorTypes"
    r = request_get(endpoint)
    if r.status_code == 200:
        print(f"Request status code: {r.status_code}, Successful request from {endpoint}")
    else:
        print(f"Request status code: {r.status_code}, Unsuccessful request from {endpoint}")
        print(r.content)
        sys.exit()

    response_json = r.json()
    response_dataframe = convert_response_to_df(response_json)
    result = subset_oncotree(response_dataframe)
    write_file(result, args.date)
