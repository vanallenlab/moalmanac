import pandas as pd
import requests

outname = 'oncotree_2018-06-01.txt'

def dict_to_df(dictionary):
    return pd.DataFrame.from_dict(dictionary, orient='index').T

def json_to_df(json):
    df = dict_to_df(json[0])
    for entry in json[1:]:
        tmp = dict_to_df(entry)
        df = pd.concat([df, tmp], ignore_index=True, axis=0)
    return df

def subset_df(df):
    column_map = {
        'code': 'code',
        'mainType': 'main_type',
        'name': 'name',
        'tissue': 'tissue'
    }

    return df.loc[:, column_map.keys()].sort_values('code').rename(columns=column_map)

def main():
	request = 'http://oncotree.mskcc.org/api/tumorTypes'
	r = requests.get(request)
	df = json_to_df(r.json())
	df = subset_df(df)
	df.to_csv(outname, sep='\t', index=False)

if __name__ == '__main__':
    main()

