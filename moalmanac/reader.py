import json
import pandas as pd
import pickle
import tinydb


class Reader(object):
    @staticmethod
    def check_comment_rows(handle, comment_character):
        skip_rows = 0
        with open(handle, 'r') as f:
            for line in f:
                if line.startswith(comment_character):
                    skip_rows += 1
                else:
                    break
        return skip_rows

    @staticmethod
    def check_column_names(df, columns_map):
        for column_name in columns_map.keys():
            assert column_name in df.columns, \
                'Expected column %s not found among %s' % (column_name, df.columns)

    @staticmethod
    def read(handle, delimiter, **kwargs):
        return pd.read_csv(handle, sep=delimiter, dtype='object', **kwargs)

    @staticmethod
    def read_json(handle):
        with open(handle, 'r') as json_handle:
            return json.load(json_handle)

    @staticmethod
    def read_pickle(handle):
        return pickle.load(open(handle, "rb"))

    @staticmethod
    def read_tinydb(handle):
        return tinydb.TinyDB(handle)

    @classmethod
    def safe_read(cls, handle, delimiter, column_map, comment_character='', **kwargs):
        if comment_character != '':
            n_comment_rows = cls.check_comment_rows(handle, comment_character)
            cls.check_column_names(cls.read(handle, delimiter, nrows=3, header=n_comment_rows, **kwargs), column_map)
            return cls.read(handle, delimiter, header=n_comment_rows,
                            usecols=column_map.keys()).rename(columns=column_map)
        else:
            cls.check_column_names(cls.read(handle, delimiter, nrows=3, **kwargs), column_map)
            return cls.read(handle, delimiter, usecols=column_map.keys(), **kwargs).rename(columns=column_map)
