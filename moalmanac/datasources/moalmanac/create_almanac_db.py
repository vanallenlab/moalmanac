import argparse
import configparser
import json
import sys


def create_config():
    config = configparser.ConfigParser()
    config.read('../../config.ini')
    return config


CONFIG = create_config()
feature_type_section = 'feature_types'
FEATURE_TYPE_SOMATIC = CONFIG[feature_type_section]['mut']
FEATURE_TYPE_GERMLINE = CONFIG[feature_type_section]['germline']
FEATURE_TYPE_COPY_NUMBER = CONFIG[feature_type_section]['cna']
FEATURE_TYPE_FUSION = CONFIG[feature_type_section]['fusion']
FEATURE_TYPE_MUTATIONAL_BURDEN = CONFIG[feature_type_section]['burden']
FEATURE_TYPE_MUTATIONAL_SIGNATURE = CONFIG[feature_type_section]['signature']
FEATURE_TYPE_MICROSATELLITE_STABILITY = CONFIG[feature_type_section]['microsatellite']
FEATURE_TYPE_ANEUPLOIDY = CONFIG[feature_type_section]['aneuploidy']
FEATURE_TYPE_KNOCKDOWN = CONFIG[feature_type_section]['knockdown']

ASSERTION_FIELDS = [
    'disease', 'context', 'oncotree_term', 'oncotree_code',
    'therapy_name', 'therapy_strategy', 'therapy_type',
    'therapy_sensitivity', 'therapy_resistance', 'favorable_prognosis',
    'predictive_implication', 'description', 'last_updated'
]

SOURCE_FIELDS = [
    'source_type', 'citation', 'url', 'doi', 'pmid', 'nct', 'publication_date', 'last_updated'
]

VARIANT_FIELDS = [
    'gene', 'chromosome', 'start_position', 'end_position', 'exon',
    'reference_allele', 'alternate_allele', 'cdna_change', 'protein_change',
    'variant_annotation', 'rsid'
]

GERMLINE_FIELDS = [
    'gene', 'chromosome', 'start_position', 'end_position', 'exon',
    'reference_allele', 'alternate_allele', 'cdna_change', 'protein_change',
    'variant_annotation', 'rsid', 'pathogenic'
]

COPY_NUMBER_FIELDS = [
    'gene', 'cytoband', 'direction'

]

REARRANGEMENT_FIELDS = [
    'gene1', 'gene2', 'rearrangement_type', 'locus'
]

MUTATIONAL_BURDEN_FIELDS = [
    'classification', 'minimum_mutations', 'mutations_per_mb'
]

MUTATIONAL_SIGNATURE_FIELDS = [
    'cosmic_signature'
]

MICROSATELLITE_FIELDS = [
    'status'
]

ANEUPLOIDY_FIELDS = [
    'event'
]

KNOCKDOWN_FIELDS = [
    'gene', 'technique'
]


class DisplayAlteration:
    @classmethod
    def aneuploidy(cls, record):
        return f"{record['event']}"

    @classmethod
    def copy_number(cls, record):
        return f"{record['gene']} {record['direction']}"

    @classmethod
    def cosmic_mutational_signature(cls, record):
        return f"Cosmic signature {record['cosmic_signature']} (version 3.4)"

    @classmethod
    def germline_variant(cls, record):
        string = cls.somatic_variant(record)
        if record['pathogenic'] == 1:
            return f"{string} (Pathogenic)"
        else:
            return string

    @classmethod
    def knockdown(cls, record):
        return f"{record['gene']} knockdown ({record['technique']})"

    @classmethod
    def microsatellite_stability(cls, record):
        return f"{record['status']}"

    @classmethod
    def mutational_burden(cls, record):
        return f"{record['classification']} mutational burden"

    @classmethod
    def rearrangements(cls, record):
        string = f"{record['gene1']}"
        if record['gene2'] != '':
            string = f"{string}--{record['gene2']}"
        if record['rearrangement_type'] != '':
            string = f"{string} {record['rearrangement_type']}"
        return string

    @classmethod
    def somatic_variant(cls, record):
        gene = record['gene']
        protein_change = record['protein_change']
        exon = record['exon']
        annotation = record['variant_annotation']
        if protein_change != '':
            return f"{gene} {protein_change} ({annotation})"
        elif exon != "" and annotation != "":
            return f"{gene} Exon {exon} ({annotation})"
        elif exon != "" and annotation == "":
            return f"{gene} Exon {exon}"
        elif annotation != "":
            return f"{gene} ({annotation})"
        else:
            return f"{gene}"

    @classmethod
    def generate(cls, record):
        feature_type = record['feature_type']
        if feature_type == FEATURE_TYPE_SOMATIC:
            return cls.somatic_variant(record)
        elif feature_type == FEATURE_TYPE_GERMLINE:
            return cls.germline_variant(record)
        elif feature_type == FEATURE_TYPE_COPY_NUMBER:
            return cls.copy_number(record)
        elif feature_type == FEATURE_TYPE_FUSION:
            return cls.rearrangements(record)
        elif feature_type == FEATURE_TYPE_ANEUPLOIDY:
            return cls.aneuploidy(record)
        elif feature_type == FEATURE_TYPE_MUTATIONAL_BURDEN:
            return cls.mutational_burden(record)
        elif feature_type == FEATURE_TYPE_MICROSATELLITE_STABILITY:
            return cls.microsatellite_stability(record)
        elif feature_type == FEATURE_TYPE_MUTATIONAL_SIGNATURE:
            return cls.cosmic_mutational_signature(record)
        elif feature_type == FEATURE_TYPE_KNOCKDOWN:
            return cls.knockdown(record)
        else:
            print(f'ERROR: {feature_type} does not have a feature display format function.')


def check_fields(fields, record):
    for field in fields:
        if field not in record.keys():
            sys.exit(f'{field} not present in record.\n{record}')


def check_format(record):
    check_fields(SOURCE_FIELDS, record)
    check_fields(ASSERTION_FIELDS, record)

    feature_type = record['feature_type']
    feature_type_fields = get_feature_type_fields(feature_type)
    check_fields(feature_type_fields, record)


def extract_genes(records):
    genes = []
    for record in records:
        for string in ['gene', 'gene1', 'gene2']:
            if string in record.keys():
                if record[string] != '':
                    genes.append(record[string])
    genes_unique = set(genes)
    genes_sorted = sorted(genes_unique)
    return genes_sorted


def get_feature_type_fields(feature_type):
    if feature_type == FEATURE_TYPE_SOMATIC:
        return VARIANT_FIELDS
    elif feature_type == FEATURE_TYPE_GERMLINE:
        return GERMLINE_FIELDS
    elif feature_type == FEATURE_TYPE_COPY_NUMBER:
        return COPY_NUMBER_FIELDS
    elif feature_type == FEATURE_TYPE_FUSION:
        return REARRANGEMENT_FIELDS
    elif feature_type == FEATURE_TYPE_MUTATIONAL_BURDEN:
        return MUTATIONAL_BURDEN_FIELDS
    elif feature_type == FEATURE_TYPE_MUTATIONAL_SIGNATURE:
        return MUTATIONAL_SIGNATURE_FIELDS
    elif feature_type == FEATURE_TYPE_MICROSATELLITE_STABILITY:
        return MICROSATELLITE_FIELDS
    elif feature_type == FEATURE_TYPE_ANEUPLOIDY:
        return ANEUPLOIDY_FIELDS
    elif feature_type == FEATURE_TYPE_KNOCKDOWN:
        return KNOCKDOWN_FIELDS
    else:
        sys.exit(f'feature type {feature_type} present in database but not accounted for.')


def initialize():
    return {
        'release': '',
        'version': '',
        'genes': [],
        'content': []
    }


def load_json(json_file):
    with open(json_file, 'r') as f:
        json_data = json.load(f)
    return json_data


def write_json(file, data):
    with open(file, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)


def main(version, release, content):
    db = initialize()
    db['release'] = release
    db['version'] = version

    db_genes = extract_genes(content)
    db['genes'] = db_genes

    for record in content:
        check_format(record)
        record['feature_display'] = DisplayAlteration.generate(record)
    db['content'] = content
    write_json('molecular-oncology-almanac.json', db)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(prog='Create Molecular Oncology Almanac datasource',
                                         description='Compiles MOAlmanac db for use with algorithm')
    arg_parser.add_argument('--file', '-f',
                            help='Molecular Oncology Almanac db file from https://github.com/vanallenlab/moalmanac-db')
    arg_parser.add_argument('--version', '-v',
                            help='Database version; e.g. 1.0.0')
    arg_parser.add_argument('--release', '-r',
                            help='Database content release; e.g. v.2022-12-01')
    args = arg_parser.parse_args()
    print(args)

    moalmanac_json = load_json(args.file)
    main(args.version, args.release, moalmanac_json)
