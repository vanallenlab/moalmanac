import configparser
import json
import pandas as pd
import pathlib
import pickle


class Ini:
    @staticmethod
    def convert_ini_to_dictionary(ini):
        dictionary = {}
        for section in ini.sections():
            dictionary[section] = {}
            for key, value in ini.items(section):
                dictionary[section][key] = value
        return dictionary

    @staticmethod
    def load(path, extended_interpolation=False):
        if extended_interpolation:
            config = configparser.ConfigParser(
                interpolation=configparser.ExtendedInterpolation()
            )
        else:
            config = configparser.ConfigParser()
        config.read(path)
        return config

    @classmethod
    def read(
        cls,
        path,
        extended_interpolation=False,
        convert_to_dictionary=False,
        resolve_paths=False,
    ):
        ini = cls.load(path, extended_interpolation=extended_interpolation)
        if convert_to_dictionary:
            data = cls.convert_ini_to_dictionary(ini)
            if resolve_paths:
                base = pathlib.Path(path).resolve().parent.parent
                data = cls.resolve_paths_dict(data=data, base_directory=base)
            return data
        else:
            return ini

    @staticmethod
    def resolve_paths_dict(data: dict, base_directory: pathlib.Path) -> dict:
        """
        Resolve all values under any 'paths' section in Ini to absolute paths.
        """
        resolved = dict(data)

        if "paths" not in resolved:
            return resolved

        out = {}
        for key, value in resolved["paths"].items():
            if value is None:
                out[key] = value
                continue

            p = pathlib.Path(value)
            if not p.is_absolute():
                revised_path = (base_directory / p).resolve()
            out[key] = str(revised_path)

        resolved["paths"] = out
        return resolved


class Reader:
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
            assert (
                str.lower(column_name) in df.columns.str.lower()
            ), 'Expected column %s not found among %s' % (
                str.lower(column_name),
                df.columns.str.lower(),
            )

    @staticmethod
    def read(handle, delimiter, **kwargs):
        return pd.read_csv(handle, sep=delimiter, dtype='object', **kwargs)

    @staticmethod
    def read_json(handle):
        with open(handle, 'r') as json_handle:
            return json.load(json_handle)

    @staticmethod
    def read_pickle(handle):
        with open(handle, 'rb') as pickle_handle:
            return pickle.load(pickle_handle)

    @staticmethod
    def return_columns_as_lowercase(dataframe):
        dataframe.columns = dataframe.columns.str.lower()
        return dataframe

    @staticmethod
    def return_keys_as_lowercase(dictionary):
        new_dictionary = {}
        for key, value in dictionary.items():
            new_dictionary[key.lower()] = value
        return new_dictionary

    @classmethod
    def safe_read(cls, handle, delimiter, column_map, comment_character='', **kwargs):
        lowercase_column_map = cls.return_keys_as_lowercase(column_map)
        if comment_character != '':
            n_comment_rows = cls.check_comment_rows(handle, comment_character)
            cls.check_column_names(
                cls.read(handle, delimiter, nrows=3, header=n_comment_rows, **kwargs),
                column_map,
            )
            df = cls.read(
                handle,
                delimiter,
                header=n_comment_rows,
                usecols=(lambda x: str.lower(str(x)) in lowercase_column_map.keys()),
                **kwargs,
            )
        else:
            cls.check_column_names(
                cls.read(handle, delimiter, nrows=3, **kwargs), column_map
            )
            df = cls.read(
                handle,
                delimiter,
                usecols=(lambda x: str.lower(str(x)) in lowercase_column_map.keys()),
                **kwargs,
            )
        df = cls.return_columns_as_lowercase(df)
        df = df.rename(columns=lowercase_column_map)
        return df
