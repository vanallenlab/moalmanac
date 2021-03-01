import argparse
import configparser
import tinydb
import pandas as pd


def create_config():
    config = configparser.ConfigParser()
    config.read('../../config.ini')
    return config


CONFIG = create_config()


class Reader(object):
    feature_type_section = 'feature_types'
    feature_type_somatic = CONFIG[feature_type_section]['mut']
    feature_type_germline = CONFIG[feature_type_section]['germline']
    feature_type_cna = CONFIG[feature_type_section]['cna']
    feature_type_fusion = CONFIG[feature_type_section]['fusion']
    feature_type_burden = CONFIG[feature_type_section]['burden']
    feature_type_signature = CONFIG[feature_type_section]['signature']
    feature_type_microsatellite = CONFIG[feature_type_section]['microsatellite']
    feature_type_aneuploidy = CONFIG[feature_type_section]['aneuploidy']
    feature_types = {
        'somatic': feature_type_somatic,
        'germline': feature_type_germline,
        'copynumber': feature_type_cna,
        'rearrangement': feature_type_fusion,
        'burden': feature_type_burden,
        'signature': feature_type_signature,
        'microsatellite': feature_type_microsatellite,
        'aneuploidy': feature_type_aneuploidy
    }

    source_columns = {
        'source_type': str, 'citation': str, 'url': str, 'doi': str, 'pmid': object, 'nct': str
    }

    assertion_columns = {
        'disease': str, 'context': str, 'oncotree_term': str, 'oncotree_code': str,
        'therapy_name': str, 'therapy_strategy': str, 'therapy_type': str,
        'therapy_sensitivity': object, 'therapy_resistance': object, 'favorable_prognosis': object,
        'predictive_implication': str, 'description': str, 'last_updated': str, 'preferred_assertion': object
    }

    variant_columns = {
        'gene': str, 'chromosome': object, 'start_position': object, 'end_position': object, 'exon': object,
        'reference_allele': str, 'alternate_allele': str, 'cdna_change': str, 'protein_change': str,
        'variant_annotation': str, 'rsid': str
    }

    somatic_columns = {
        **variant_columns,
        **assertion_columns, **source_columns}
    germline_columns = {
        **variant_columns, 'pathogenic': str,
        **assertion_columns, **source_columns}
    copynumber_columns = {
        'gene': str, 'direction': str, 'cytoband': str,
        **assertion_columns, **source_columns}
    rearrangement_columns = {
        'gene1': str, 'gene2': str, 'rearrangement_type': str, 'locus': str,
        **assertion_columns, **source_columns}
    aneuploidy_columns = {
        'event': str,
        **assertion_columns, **source_columns}
    microsatellite_stability_columns = {
        'status': str,
        **assertion_columns, **source_columns}
    mutational_burden_columns = {
        'classification': str, 'minimum_mutations': object, 'mutations_per_mb': float,
        **assertion_columns, **source_columns}
    mutational_signature_columns = {
        'cosmic_signature_number': object,
        **assertion_columns, **source_columns}

    feature_type_dictionary = {
        'somatic': {'handle': 'somatic_variant',
                    'columns': somatic_columns},
        'germline': {'handle': 'germline_variant',
                     'columns': germline_columns},
        'copynumber': {'handle': 'copy_number',
                       'columns': copynumber_columns},
        'rearrangement': {'handle': 'rearrangement',
                          'columns': rearrangement_columns},
        'aneuploidy': {'handle': 'aneuploidy',
                       'columns': aneuploidy_columns},
        'microsatellite': {'handle': 'microsatellite_stability',
                           'columns': microsatellite_stability_columns},
        'burden': {'handle': 'mutational_burden',
                   'columns': mutational_burden_columns},
        'signature': {'handle': 'mutational_signature',
                      'columns': mutational_signature_columns}
    }

    for feature_type in feature_type_dictionary.keys():
        feature_type_dictionary[feature_type]['label'] = feature_types[feature_type]

    @classmethod
    def read_csv(cls, label, dictionary, root):
        column_types = dictionary[label]['columns']
        columns = column_types.keys()
        handle_feature_type = dictionary[label]['handle']
        handle = '{}{}.tsv'.format(root, handle_feature_type)
        return pd.read_csv(handle, sep='\t', usecols=columns, dtype=column_types)


#    predictive_implication_map = {
#        'FDA-Approved': 5.0, 'Guideline': 4.0, 'Clinical trial': 3.0,
#        'Clinical evidence': 2.0, 'Preclinical': 1.0, 'Inferential': 0.0}
#    @classmethod
#    def map_predictive_implication(cls, series):
#        return series.replace(cls.predictive_implication_map)


class DBCreator(object):
    basename = 'moalmanac'

    @classmethod
    def initialize(cls):
        handle = f'{cls.basename}.json'
        print('Documents will be added to {}'.format(handle))
        return tinydb.TinyDB(handle)

    @classmethod
    def create_table(cls, db, table_name):
        return db.table(table_name)

    @classmethod
    def insert_documents(cls, table, documents):
        table.insert_multiple(documents)

    @classmethod
    def create_gene_list(cls, db):
        gene_list = []
        for table_name in db.tables():
            table = db.table(table_name)
            for gene_column in ['gene', 'gene1', 'gene2']:
                query = tinydb.Query()[gene_column].exists()
                results = table.search(query)
                table_genes = [document[gene_column] for document in results if str(document[gene_column]) != 'nan']
                gene_list.extend(table_genes)
        return sorted(list((set(gene_list))))

    @classmethod
    def create_gene_table(cls, db):
        gene_list = cls.create_gene_list(db)
        genes_dictionary = {'genes': gene_list}

        table = db.table('genes')
        table.insert(genes_dictionary)

    @classmethod
    def create_release_table(cls, db, version):
        release_dictionary = {'release': version}
        table = db.table('Release')
        table.insert(release_dictionary)


class Formatter(object):
    feature_types = Reader.feature_types
    feature_display = 'feature_display'

    @classmethod
    def format_aneuploidy(cls, df):
        df.loc[:, cls.feature_display] = df.loc[:, 'event']
        return df

    @classmethod
    def format_mutational_burden(cls, df):
        df.loc[:, cls.feature_display] = df.loc[:, 'classification'] + ' mutational burden'
        return df

    @classmethod
    def format_copy_number_variants(cls, df):
        df.loc[:, cls.feature_display] = df.loc[:, 'gene'] + ' ' + df.loc[:, 'direction']
        return df

    @classmethod
    def format_germline_variants(cls, df):
        df = cls.format_somatic_variants(df)
        df.loc[:, 'pathogenic'] = df.loc[:, 'pathogenic'].fillna('')

        for idx in df.index:
            pathogenic = df.loc[idx, 'pathogenic']
            if pathogenic == 1:
                display = df.loc[idx, cls.feature_display] + ' (Pathogenic)'
            else:
                continue
            df.loc[idx, cls.feature_display] = display
        return df

    @classmethod
    def format_microsatellite_stability(cls, df):
        df.loc[:, cls.feature_display] = df.loc[:, 'status']
        return df

    @classmethod
    def format_rearrangements(cls, df):
        df.loc[:, cls.feature_display] = df.loc[:, 'gene1']

        idx_gene2_not_null = df[~df['gene2'].isnull()].index
        df.loc[idx_gene2_not_null, cls.feature_display] = (
                df.loc[idx_gene2_not_null, cls.feature_display] + '--' + df.loc[idx_gene2_not_null, 'gene2'])

        idx_type_not_null = df[~df['rearrangement_type'].isnull()].index
        df.loc[idx_type_not_null, cls.feature_display] = (
                df.loc[idx_type_not_null, cls.feature_display] + ' ' + df.loc[idx_type_not_null, 'rearrangement_type'])
        return df

    @classmethod
    def format_signatures(cls, df):
        df.loc[:, cls.feature_display] = 'Cosmic signature ' + df.loc[:, 'cosmic_signature_number']
        return df

    @classmethod
    def format_somatic_variants(cls, df):
        df.loc[:, cls.feature_display] = ''

        columns = ['gene', 'protein_change', 'exon', 'variant_annotation']
        df.loc[:, columns] = df.loc[:, columns].fillna('')
        for idx in df.index:
            gene = df.loc[idx, 'gene']
            protein_change = df.loc[idx, 'protein_change']
            exon = df.loc[idx, 'exon']
            annotation = df.loc[idx, 'variant_annotation']

            if protein_change != '':
                display = '{} {} ({})'.format(gene, protein_change, annotation)
            elif exon != '' and annotation != '':
                display = '{} Exon {} ({})'.format(gene, exon, annotation)
            elif exon != '' and annotation == '':
                display = '{} Exon {}'.format(gene, exon)
            elif annotation != '':
                display = '{} ({})'.format(gene, annotation)
            else:
                display = '{}'.format(gene)
            df.loc[idx, cls.feature_display] = display
        return df

    @classmethod
    def format_feature_display(cls, feature_type, dataframe):
        if feature_type == cls.feature_types['somatic']:
            return cls.format_somatic_variants(dataframe)
        elif feature_type == cls.feature_types['germline']:
            return cls.format_germline_variants(dataframe)
        elif feature_type == cls.feature_types['copynumber']:
            return cls.format_copy_number_variants(dataframe)
        elif feature_type == cls.feature_types['rearrangement']:
            return cls.format_rearrangements(dataframe)
        elif feature_type == cls.feature_types['aneuploidy']:
            return cls.format_aneuploidy(dataframe)
        elif feature_type == cls.feature_types['burden']:
            return cls.format_mutational_burden(dataframe)
        elif feature_type == cls.feature_types['microsatellite']:
            return cls.format_microsatellite_stability(dataframe)
        elif feature_type == cls.feature_types['signature']:
            return cls.format_signatures(dataframe)
        else:
            print('ERROR: {} does not have a feature display format function!'.format(feature_type))
            return dataframe


def main(directory, version):
    database = DBCreator.initialize()

    for feature_type in Reader.feature_type_dictionary.keys():
        table = DBCreator.create_table(database, Reader.feature_type_dictionary[feature_type]['label'])
        tmp = Reader.read_csv(feature_type, Reader.feature_type_dictionary, directory)
        tmp = Formatter.format_feature_display(Reader.feature_types[feature_type], tmp)
        # tmp['predictive_implication_map'] = Reader.map_predictive_implication(tmp.loc[:, 'predictive_implication'])
        documents = tmp.to_dict(orient='records')
        DBCreator.insert_documents(table, documents)

    DBCreator.create_gene_table(database)
    DBCreator.create_release_table(database, version)

    print(database.tables)
    print('creation, complete!')


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(prog='Create Molecular Oncology Almanac datasource',
                                         description='Compiles molecular oncology almanac into TinyDB for the method.')
    arg_parser.add_argument('--directory', '-d',
                            help='Root path to database flat files of feature types',
                            default='/Users/brendan/Github/moalmanac-db/content/')
    arg_parser.add_argument('--version', '-v',
                            help='Version for output database',
                            default='')
    args = arg_parser.parse_args()

    main(args.directory, args.version)
