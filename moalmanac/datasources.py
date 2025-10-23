import pandas as pd

from reader import Reader
from config import COLNAMES


class Datasources:
    """
    Loads datasources
    """
    datasources_section = 'datasources'
    feature_type = COLNAMES[datasources_section]['feature_type']
    feature = COLNAMES[datasources_section]['feature']
    alt_type = COLNAMES[datasources_section]['alt_type']
    alt = COLNAMES[datasources_section]['alt']

    chr = COLNAMES[datasources_section]['chr']
    start = COLNAMES[datasources_section]['start']
    end = COLNAMES[datasources_section]['end']
    ref = COLNAMES[datasources_section]['ref']
    allele2 = COLNAMES[datasources_section]['allele2']

    disease = COLNAMES[datasources_section]['disease']
    ontology = COLNAMES[datasources_section]['ontology']
    code = COLNAMES[datasources_section]['code']
    context = COLNAMES[datasources_section]['context']
    mutational_burden = COLNAMES[datasources_section]['mutational_burden']
    therapy = COLNAMES[datasources_section]['therapy']
    therapy_strategy = COLNAMES[datasources_section]['therapy_strategy']
    therapy_type = COLNAMES[datasources_section]['therapy_type']

    sensitivity = COLNAMES[datasources_section]['sensitivity']
    resistance = COLNAMES[datasources_section]['resistance']
    prognosis = COLNAMES[datasources_section]['prognosis']
    implication = COLNAMES[datasources_section]['implication']
    implication_map = COLNAMES[datasources_section]['implication_map']
    connections = COLNAMES[datasources_section]['connections']
    description = COLNAMES[datasources_section]['description']

    source_type = COLNAMES[datasources_section]['source_type']
    citation = COLNAMES[datasources_section]['citation']
    url = COLNAMES[datasources_section]['url']
    doi = COLNAMES[datasources_section]['doi']
    pmid = COLNAMES[datasources_section]['pmid']
    nct = COLNAMES[datasources_section]['nct']
    publication_date = COLNAMES[datasources_section]['publication_date']
    last_updated = COLNAMES[datasources_section]['last_updated']

    sensitivity_matches = COLNAMES[datasources_section]['sensitivity_matches']
    resistance_matches = COLNAMES[datasources_section]['resistance_matches']
    prognostic_matches = COLNAMES[datasources_section]['prognostic_matches']

    columns = COLNAMES[datasources_section]['columns']
    query = COLNAMES[datasources_section]['query']
    genes = COLNAMES[datasources_section]['genes']

    af = COLNAMES[datasources_section]['exac_af']
    ac = COLNAMES[datasources_section]['exac_ac']
    an = COLNAMES[datasources_section]['exac_an']
    ac_afr = COLNAMES[datasources_section]['exac_afr_ac']
    ac_amr = COLNAMES[datasources_section]['exac_amr_ac']
    ac_eas = COLNAMES[datasources_section]['exac_eas_ac']
    ac_fin = COLNAMES[datasources_section]['exac_fin_ac']
    ac_nfe = COLNAMES[datasources_section]['exac_nfe_ac']
    ac_sas = COLNAMES[datasources_section]['exac_sas_ac']
    ac_oth = COLNAMES[datasources_section]['exac_oth_ac']
    an_afr = COLNAMES[datasources_section]['exac_afr_an']
    an_amr = COLNAMES[datasources_section]['exac_amr_an']
    an_eas = COLNAMES[datasources_section]['exac_eas_an']
    an_fin = COLNAMES[datasources_section]['exac_fin_an']
    an_nfe = COLNAMES[datasources_section]['exac_nfe_an']
    an_sas = COLNAMES[datasources_section]['exac_sas_an']
    an_oth = COLNAMES[datasources_section]['exac_oth_an']


class ACMG:
    gene = Datasources.feature

    column_map = {
        'Gene': gene
    }

    @classmethod
    def import_ds(cls, dbs):
        return Reader.safe_read(dbs['acmg_handle'], '\t', cls.column_map, comment_character='#')


class Almanac:
    feature_type = Datasources.feature_type
    feature = Datasources.feature
    alt_type = Datasources.alt_type
    alt = Datasources.alt

    disease = Datasources.disease
    ontology = Datasources.ontology
    code = Datasources.code
    context = Datasources.context
    therapy = Datasources.therapy
    therapy_strategy = Datasources.therapy_strategy
    therapy_type = Datasources.therapy_type
    sensitivity = Datasources.sensitivity
    resistance = Datasources.resistance
    prognosis = Datasources.prognosis
    implication = Datasources.implication
    implication_map = Datasources.implication_map
    description = Datasources.description

    source_type = Datasources.source_type
    citation = Datasources.citation
    url = Datasources.url
    doi = Datasources.doi
    pmid = Datasources.pmid
    nct = Datasources.nct
    publication_date = Datasources.publication_date
    last_updated = Datasources.last_updated

    sensitivity_matches = Datasources.sensitivity_matches
    resistance_matches = Datasources.resistance_matches
    prognostic_matches = Datasources.prognostic_matches

    columns = Datasources.columns
    query = Datasources.query
    genes = Datasources.genes

    predictive_implication_map = {
        'FDA-Approved': 5.0,
        'Guideline': 4.0,
        'Clinical trial': 3.0,
        'Clinical evidence': 2.0,
        'Preclinical': 1.0,
        'Inferential': 0.0
    }

    @staticmethod
    def import_ds(dbs):
        ds = Reader.read_json(dbs['almanac_handle'])
        return ds

    @classmethod
    def import_genes(cls, dbs):
        ds = cls.import_ds(dbs)
        return ds['genes']


class CancerGeneCensus:
    gene = Datasources.feature

    column_map = {
        'Gene Symbol': gene
    }

    @classmethod
    def import_ds(cls, dbs):
        return Reader.safe_read(dbs['cgc_handle'], '\t', cls.column_map)


class CancerHotspots:
    gene = Datasources.feature
    alt = Datasources.alt
    aa_pos = 'aa_pos'
    aa_ref = 'aa_ref'
    aa_var = 'aa_var'

    column_map = {
        'Hugo_Symbol': gene,
        'Amino_Acid_Position': aa_pos,
        'Reference_Amino_Acid': aa_ref,
        'Variant_Amino_Acid': aa_var
    }

    @classmethod
    def format_cancerhotspots(cls, df):
        df[cls.alt] = 'p.' + df.loc[:, cls.aa_ref].str.split(':', expand=True).iloc[:, 0].astype(str)
        df[cls.alt] += df[cls.aa_pos].astype(str)
        df[cls.alt] += df[cls.aa_var].str.split(':', expand=True).iloc[:, 0].astype(str)
        return df

    @classmethod
    def import_ds(cls, dbs):
        df = Reader.safe_read(dbs['cancerhotspots_handle'], '\t', cls.column_map)
        return cls.format_cancerhotspots(df)


class CancerHotspots3D:
    gene = Datasources.feature
    alt = Datasources.alt

    column_map = {
        'Gene': gene,
        'alteration': alt
    }

    @classmethod
    def import_ds(cls, dbs):
        return Reader.safe_read(dbs['3dcancerhotspots_handle'], '\t', cls.column_map)


class ClinVar:
    gene = Datasources.feature
    chr = Datasources.chr
    start = Datasources.start
    end = Datasources.end
    ref = Datasources.ref
    allele2 = Datasources.allele2

    bin_section = 'bin_names'
    bin_name = COLNAMES[bin_section]['clinvar']
    bin_map = bin_name.split('_')[0]

    column_map = {
        'GeneSymbol': gene,
        'Chromosome': chr,
        'Start': start,
        'Stop': end,
        'ReferenceAllele': ref,
        'AlternateAllele': allele2,
        'ClinicalSignificance': bin_map,
        'ClinSigSimple': bin_name
    }

    @classmethod
    def import_ds(cls, dbs):
        return Reader.safe_read(dbs['clinvar_handle'], '\t', cls.column_map)


class Cosmic:
    gene = Datasources.feature
    alt = Datasources.alt

    column_map = {
        'Gene name': gene,
        'Mutation AA': alt
    }

    @classmethod
    def import_ds(cls, dbs):
        return Reader.safe_read(dbs['cosmic_handle'], '\t', cls.column_map)


class ExAC:
    chr = Datasources.chr
    start = Datasources.start
    ref = Datasources.ref
    alt = Datasources.allele2
    af = Datasources.af

    column_map = {
        'CHROM': chr,
        'POS': start,
        'REF': ref,
        'ALT': alt,
        'AF': af
    }

    @classmethod
    def import_ds(cls, dbs):
        return Reader.safe_read(dbs['exac_handle'], '\t', cls.column_map)


class ExACExtended:
    chr = Datasources.chr
    start = Datasources.start
    ref = Datasources.ref
    alt = Datasources.allele2
    af = Datasources.af
    ac = Datasources.ac
    an = Datasources.an
    ac_afr = Datasources.ac_afr
    ac_amr = Datasources.ac_amr
    ac_eas = Datasources.ac_eas
    ac_fin = Datasources.ac_fin
    ac_nfe = Datasources.ac_nfe
    ac_sas = Datasources.ac_sas
    ac_oth = Datasources.ac_oth
    an_afr = Datasources.an_afr
    an_amr = Datasources.an_amr
    an_eas = Datasources.an_eas
    an_fin = Datasources.an_fin
    an_nfe = Datasources.an_nfe
    an_sas = Datasources.an_sas
    an_oth = Datasources.an_oth

    column_map = {
        'CHROM': chr,
        'POS': start,
        'REF': ref,
        'ALT': alt,
        'AF': af,
        'AC': ac,
        'AN': an,
        'AC_AFR': ac_afr,
        'AC_AMR': ac_amr,
        'AC_EAS': ac_eas,
        'AC_FIN': ac_fin,
        'AC_NFE': ac_nfe,
        'AC_SAS': ac_sas,
        'AC_OTH': ac_oth,
        'AN_AFR': an_afr,
        'AN_AMR': an_amr,
        'AN_EAS': an_eas,
        'AN_FIN': an_fin,
        'AN_NFE': an_nfe,
        'AN_SAS': an_sas,
        'AN_OTH': an_oth,
    }

    @classmethod
    def import_ds(cls, dbs):
        return Reader.safe_read(dbs['exac_handle'], '\t', cls.column_map)


class GSEACancerPathways:
    gene = Datasources.feature

    @classmethod
    def format_pathways(cls, df):
        pathway_genes = []
        for idx in df.index:
            pathway_genes += df.iloc[idx, 2:].dropna().tolist()
        return pd.DataFrame(pathway_genes, columns=[cls.gene])

    @classmethod
    def import_ds(cls, dbs):
        df = Reader.read(dbs['gsea_pathways_handle'], '\t', header=None)
        return cls.format_pathways(df)


class GSEACancerModules:
    gene = Datasources.feature

    @classmethod
    def format_modules(cls, df):
        df = df.loc[:, 0].str.split('\t')

        module_genes = []
        for module in df.tolist():
            module_genes += module[2:]

        module_genes_unique = sorted(list(set(module_genes)))
        return pd.DataFrame(module_genes_unique, columns=[cls.gene])

    @classmethod
    def import_ds(cls, dbs):
        df = Reader.read(dbs['gsea_modules_handle'], ',', header=None)
        return cls.format_modules(df)


class Hereditary:
    gene = Datasources.feature

    column_map = {
        'gene_name': gene,
    }

    @classmethod
    def import_ds(cls, dbs):
        return Reader.safe_read(dbs['hereditary_handle'], '\t', cls.column_map)


class Lawrence:
    code = Datasources.code
    ontology = Datasources.ontology
    mutational_burden = Datasources.mutational_burden

    column_map = {
        'code': code,
        'ontology': ontology,
        'non-silent per Mb': mutational_burden
    }

    @classmethod
    def import_ds(cls, dbs):
        return Reader.safe_read(dbs['lawrence_handle'], '\t', cls.column_map)


class Oncotree:
    ontology = Datasources.ontology
    code = Datasources.code

    column_map = {
        'name': ontology,
        'code': code
    }

    @classmethod
    def import_ds(cls, dbs):
        return Reader.safe_read(dbs['oncotree_handle'], '\t', cls.column_map)


class Preclinical:
    section = 'preclinical'
    feature = COLNAMES[section]['feature']
    partner = COLNAMES[section]['partner']
    gene = COLNAMES[section]['gene']
    use_column = COLNAMES[section]['use_column']
    model_id = COLNAMES[section]['model_id']
    broad = COLNAMES[section]['broad']
    ccle_name = COLNAMES[section]['ccle_name']
    sanger = COLNAMES[section]['sanger']

    summary = 'summary'
    variants = 'variants'
    cnas = 'copy_number_alterations'
    fusions = 'fusions'
    fusions_gene1 = 'fusions_gene1'
    fusions_gene2 = 'fusions_gene2'
    gdsc = 'gdsc'
    mappings = 'mappings'
    dictionary = 'dictionary'

    @classmethod
    def create_convert_names_dict(cls, dataframe, map_from, map_to):
        return dataframe.loc[:, [map_from, map_to]].dropna().set_index(map_from)[map_to].to_dict()

    @staticmethod
    def generate_sample_list(dataframe, use_column, sample_column):
        return dataframe[dataframe[use_column].astype(bool).astype(int).eq(1)][sample_column].sort_values().tolist()

    @classmethod
    def import_dbs(cls, paths_dictionary):
        summary = Reader.read(paths_dictionary['summary'], delimiter='\t')
        variants = Reader.read(paths_dictionary['variants'], delimiter='\t', low_memory=False)
        cnas = Reader.read(paths_dictionary['copynumbers'], delimiter='\t', low_memory=False)
        fusions = Reader.read(paths_dictionary['fusions'], delimiter='\t', low_memory=False)
        fusions1 = Reader.read(paths_dictionary['fusions1'], delimiter='\t', low_memory=False)
        fusions2 = Reader.read(paths_dictionary['fusions2'], delimiter='\t', low_memory=False)
        gdsc = Reader.read(paths_dictionary['gdsc'], delimiter='\t', low_memory=False)
        mappings = Reader.read_json(paths_dictionary['almanac_gdsc_mappings'])
        dictionary = Reader.read_pickle(paths_dictionary['dictionary'])

        ccle_map = cls.create_convert_names_dict(summary, cls.ccle_name, cls.broad)
        sanger_map = cls.create_convert_names_dict(summary, cls.sanger, cls.broad)

        summary[cls.model_id] = summary[cls.broad]
        variants[cls.model_id] = variants[cls.model_id].replace(ccle_map)
        cnas[cls.model_id] = cnas[cls.model_id].replace(ccle_map)
        fusions[cls.model_id] = fusions[cls.model_id].replace(sanger_map)
        fusions1[cls.model_id] = fusions1[cls.model_id].replace(sanger_map)
        fusions2[cls.model_id] = fusions2[cls.model_id].replace(sanger_map)
        gdsc[cls.model_id] = gdsc[cls.model_id].replace(sanger_map)
        samples = cls.generate_sample_list(summary, cls.use_column, cls.model_id)

        dbs = {
            cls.summary: cls.subset_dataframe(summary, cls.model_id, samples),
            cls.variants: cls.subset_dataframe(variants, cls.model_id, samples),
            cls.cnas: cls.subset_dataframe(cnas, cls.model_id, samples),
            cls.fusions: cls.subset_dataframe(fusions, cls.model_id, samples),
            cls.gdsc: cls.subset_dataframe(gdsc, cls.model_id, samples),
            cls.mappings: mappings,
            cls.dictionary: dictionary
        }
        dbs[cls.gene] = cls.record_gene_hits(dbs[cls.variants], dbs[cls.cnas], dbs[cls.fusions])
        return dbs

    @staticmethod
    def subset_dataframe(dataframe, column, subset_list):
        return dataframe[dataframe[column].isin(subset_list)].reset_index(drop=True)

    @classmethod
    def record_gene_hits(cls, variants, copy_numbers, fusions):
        df = pd.concat([
            variants.loc[:, [cls.feature, cls.model_id]],
            copy_numbers.loc[:, [cls.feature, cls.model_id]],
            fusions.loc[:, [cls.feature, cls.model_id]],
            fusions.loc[:, [cls.partner, cls.model_id]].rename(columns={cls.partner: cls.feature}),
        ]).drop_duplicates()
        return df
