import configparser

default_colnames = 'colnames.ini'


def create_colnames_dict(config):
    dictionary = {}
    for section in config.sections():
        dictionary[section] = {}
        for (key, value) in config.items(section):
            dictionary[section][key] = value
    return dictionary


def create_colnames():
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(default_colnames)
    return create_colnames_dict(config)


COLNAMES = create_colnames()
