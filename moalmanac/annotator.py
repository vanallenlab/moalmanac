import pandas as pd
import numpy as np
import scipy.stats as stats
import operator as op
import tinydb
import copy

import datasources
from features import Features

from config import COLNAMES
from config import CONFIG

EXAC_CONFIG = CONFIG['exac']


class Annotator(object):
    """
    Annotates variants using datasources
    """
    bin_section = 'bin_names'
    score_bin = COLNAMES[bin_section]['score_bin']
    almanac_bin = COLNAMES[bin_section]['almanac']
    sensitive_bin = COLNAMES[bin_section]['sensitive']
    resistance_bin = COLNAMES[bin_section]['resistance']
    prognostic_bin = COLNAMES[bin_section]['prognostic']
    sensitive_score_bin = COLNAMES[bin_section]['sensitive_score_bin']
    resistance_score_bin = COLNAMES[bin_section]['resistance_score_bin']
    prognostic_score_bin = COLNAMES[bin_section]['prognostic_score_bin']
    acmg_bin = COLNAMES[bin_section]['acmg']
    cancerhotspots_bin = COLNAMES[bin_section]['cancerhotspots']
    cancerhotspots3d_bin = COLNAMES[bin_section]['cancerhotspots3d']
    cgc_bin = COLNAMES[bin_section]['cgc']
    clinvar_bin = COLNAMES[bin_section]['clinvar']
    cosmic_bin = COLNAMES[bin_section]['cosmic']
    gsea_module_bin = COLNAMES[bin_section]['gsea_modules']
    gsea_pathway_bin = COLNAMES[bin_section]['gsea_pathways']
    hereditary_bin = COLNAMES[bin_section]['hereditary']
    msi_bin = COLNAMES[bin_section]['msi']
    exac_common_bin = COLNAMES[bin_section]['exac_common']

    @classmethod
    def annotate(cls, df, dbs, importer, bin_name, comparison_columns):
        ds = importer.import_ds(dbs)
        df[bin_name] = cls.match_ds(df, ds, bin_name, comparison_columns)
        return df

    @classmethod
    def annotate_almanac(cls, df, dbs, ontology):
        df[cls.score_bin] = cls.preallocate_bin(cls.score_bin, df.index)
        df = Almanac.annotate(df, dbs, ontology)
        return df

    @classmethod
    def annotate_germline(cls, df, dbs, ontology):
        df[cls.score_bin] = cls.preallocate_bin(cls.score_bin, df.index)
        df = Almanac.annotate(df, dbs, ontology)
        df = CancerHotspots.annotate(df, dbs)
        df = CancerGeneCensus.annotate(df, dbs)
        df = ACMG.annotate(df, dbs)
        df = ClinVar.annotate(df, dbs)
        df = Hereditary.annotate(df, dbs)
        df = ExACExtended.annotate(df, dbs)
        df = MSI.annotate(df)
        return df

    @classmethod
    def annotate_somatic(cls, df, dbs, ontology):
        df[cls.score_bin] = cls.preallocate_bin(cls.score_bin, df.index)
        df = Almanac.annotate(df, dbs, ontology)
        df = CancerHotspots.annotate(df, dbs)
        df = CancerHotspots3D.annotate(df, dbs)
        df = CancerGeneCensus.annotate(df, dbs)
        df = Cosmic.annotate(df, dbs)
        df = GSEACancerPathways.annotate(df, dbs)
        df = GSEACancerModules.annotate(df, dbs)
        df = ExAC.annotate(df, dbs)
        df = MSI.annotate(df)
        return df

    @staticmethod
    def compare_ids_true_index(id1, id2):
        return id1[id1.isin(id2)].index

    @staticmethod
    def create_id(df, columns):
        idx = df.loc[:, columns[0]].astype(str)
        if len(columns) > 1:
            for col in columns[1:]:
                idx += '_' + df.loc[:, col].astype(str)
        idx = idx.dropna()
        return idx

    @staticmethod
    def create_id_series(series, columns):
        idx = str(series[columns[0]])
        if len(columns) > 1:
            for col in columns[1:]:
                idx += '_' + str(series[col])
        return idx

    @classmethod
    def match_ds(cls, df, ds, bin_column, compare_columns):
        df[bin_column] = cls.preallocate_bin(bin_column, df.index)

        for i in range(len(compare_columns)):
            cols = compare_columns[:i + 1]
            idx_match = cls.compare_ids_true_index(cls.create_id(df, cols), cls.create_id(ds, cols))
            df.loc[idx_match, bin_column] = int(i + 1)
        return df[bin_column]

    @staticmethod
    def preallocate_bin(bin_name, idx):
        return pd.Series(0, name=bin_name, index=idx)

    @staticmethod
    def preallocate_empty_columns(df, columns):
        for col in columns:
            df.loc[:, col] = None
        return df


class ACMG(object):
    gene = datasources.ACMG.gene

    bin_name = Annotator.acmg_bin
    comparison_columns = [gene]

    @classmethod
    def annotate(cls, df, dbs):
        return Annotator.annotate(df, dbs, datasources.ACMG, cls.bin_name, cls.comparison_columns)


class Almanac(object):
    feature_type = datasources.Almanac.feature_type
    feature = datasources.Almanac.feature
    alt_type = datasources.Almanac.alt_type
    alt = datasources.Almanac.alt

    disease = datasources.Almanac.disease
    oncotree_term = datasources.Almanac.ontology
    oncotree_code = datasources.Almanac.code
    context = datasources.Almanac.context
    therapy = datasources.Almanac.therapy
    therapy_type = datasources.Almanac.therapy_type
    sensitivity = datasources.Almanac.sensitivity
    resistance = datasources.Almanac.resistance
    prognosis = datasources.Almanac.prognosis
    implication = datasources.Almanac.implication
    implication_map = datasources.Almanac.implication_map
    description = datasources.Almanac.description
    preferred_assertion = datasources.Almanac.preferred_assertion

    source_type = datasources.Almanac.source_type
    citation = datasources.Almanac.citation
    url = datasources.Almanac.url
    doi = datasources.Almanac.doi
    pmid = datasources.Almanac.pmid
    nct = datasources.Almanac.nct

    predictive_implication_map = datasources.Almanac.predictive_implication_map

    datasources_section = 'datasources'
    datasources_colnames = COLNAMES[datasources_section]
    sensitive_oncotree_code = COLNAMES[datasources_section]['sensitive_oncotree_code']
    sensitive_oncotree_term = COLNAMES[datasources_section]['sensitive_oncotree_term']
    sensitive_context = COLNAMES[datasources_section]['sensitive_context']
    sensitive_feature_type = COLNAMES[datasources_section]['sensitive_feature_type']
    sensitive_feature = COLNAMES[datasources_section]['sensitive_feature']
    sensitive_alt_type = COLNAMES[datasources_section]['sensitive_alt_type']
    sensitive_alt = COLNAMES[datasources_section]['sensitive_alt']
    sensitive_feature_display = COLNAMES[datasources_section]['sensitive_feature_display']
    sensitive_implication = COLNAMES[datasources_section]['sensitive_implication']
    sensitive_implication_map = COLNAMES[datasources_section]['sensitive_implication_map']
    sensitive_therapy = COLNAMES[datasources_section]['sensitive_therapy']
    sensitive_therapy_type = COLNAMES[datasources_section]['sensitive_therapy_type']
    sensitive_description = COLNAMES[datasources_section]['sensitive_description']
    sensitive_source_type = COLNAMES[datasources_section]['sensitive_source_type']
    sensitive_citation = COLNAMES[datasources_section]['sensitive_citation']
    sensitive_url = COLNAMES[datasources_section]['sensitive_url']
    sensitive_doi = COLNAMES[datasources_section]['sensitive_doi']
    sensitive_pmid = COLNAMES[datasources_section]['sensitive_pmid']
    sensitive_nct = COLNAMES[datasources_section]['sensitive_nct']

    resistance_oncotree_code = COLNAMES[datasources_section]['resistance_oncotree_code']
    resistance_oncotree_term = COLNAMES[datasources_section]['resistance_oncotree_term']
    resistance_context = COLNAMES[datasources_section]['resistance_context']
    resistance_feature_type = COLNAMES[datasources_section]['resistance_feature_type']
    resistance_feature = COLNAMES[datasources_section]['resistance_feature']
    resistance_alt_type = COLNAMES[datasources_section]['resistance_alt_type']
    resistance_alt = COLNAMES[datasources_section]['resistance_alt']
    resistance_feature_display = COLNAMES[datasources_section]['resistance_feature_display']
    resistance_implication = COLNAMES[datasources_section]['resistance_implication']
    resistance_implication_map = COLNAMES[datasources_section]['resistance_implication_map']
    resistance_therapy = COLNAMES[datasources_section]['resistance_therapy']
    resistance_therapy_type = COLNAMES[datasources_section]['resistance_therapy_type']
    resistance_description = COLNAMES[datasources_section]['resistance_description']
    resistance_source_type = COLNAMES[datasources_section]['resistance_source_type']
    resistance_citation = COLNAMES[datasources_section]['resistance_citation']
    resistance_url = COLNAMES[datasources_section]['resistance_url']
    resistance_doi = COLNAMES[datasources_section]['resistance_doi']
    resistance_pmid = COLNAMES[datasources_section]['resistance_pmid']
    resistance_nct = COLNAMES[datasources_section]['resistance_nct']

    prognostic_oncotree_code = COLNAMES[datasources_section]['prognostic_oncotree_code']
    prognostic_oncotree_term = COLNAMES[datasources_section]['prognostic_oncotree_term']
    prognostic_context = COLNAMES[datasources_section]['prognostic_context']
    prognostic_feature_type = COLNAMES[datasources_section]['prognostic_feature_type']
    prognostic_feature = COLNAMES[datasources_section]['prognostic_feature']
    prognostic_alt_type = COLNAMES[datasources_section]['prognostic_alt_type']
    prognostic_alt = COLNAMES[datasources_section]['prognostic_alt']
    prognostic_feature_display = COLNAMES[datasources_section]['prognostic_feature_display']
    prognostic_implication = COLNAMES[datasources_section]['prognostic_implication']
    prognostic_implication_map = COLNAMES[datasources_section]['prognostic_implication_map']
    prognostic_description = COLNAMES[datasources_section]['prognostic_description']
    prognostic_source_type = COLNAMES[datasources_section]['prognostic_source_type']
    prognostic_citation = COLNAMES[datasources_section]['prognostic_citation']
    prognostic_url = COLNAMES[datasources_section]['prognostic_url']
    prognostic_doi = COLNAMES[datasources_section]['prognostic_doi']
    prognostic_pmid = COLNAMES[datasources_section]['prognostic_pmid']
    prognostic_nct = COLNAMES[datasources_section]['prognostic_nct']

    feature_display = COLNAMES[datasources_section]['feature_display']
    gene = COLNAMES[datasources_section]['gene']
    variant_annotation = COLNAMES[datasources_section]['variant_annotation']
    protein_change = COLNAMES[datasources_section]['protein_change']
    direction = COLNAMES[datasources_section]['direction']
    gene1 = COLNAMES[datasources_section]['gene1']
    gene2 = COLNAMES[datasources_section]['gene2']
    rearrangement_type = COLNAMES[datasources_section]['rearrangement_type']
    classification = COLNAMES[datasources_section]['classification']
    event = COLNAMES[datasources_section]['event']
    cosmic_signature_number = COLNAMES[datasources_section]['cosmic_signature_number']
    status = COLNAMES[datasources_section]['status']

    column_map_sensitive = {
        feature_display: sensitive_feature_display,
        oncotree_term: sensitive_oncotree_term,
        oncotree_code: sensitive_oncotree_code,
        context: sensitive_context,
        therapy: sensitive_therapy,
        therapy_type: sensitive_therapy_type,
        implication: sensitive_implication,
        implication_map: sensitive_implication_map,
        description: sensitive_description,
        source_type: sensitive_source_type,
        citation: sensitive_citation,
        url: sensitive_url,
        doi: sensitive_doi,
        pmid: sensitive_pmid,
        nct: sensitive_nct
    }

    column_map_resistance = {
        feature_display: resistance_feature_display,
        oncotree_term: resistance_oncotree_term,
        oncotree_code: resistance_oncotree_code,
        context: resistance_context,
        therapy: resistance_therapy,
        therapy_type: resistance_therapy_type,
        implication: resistance_implication,
        implication_map: resistance_implication_map,
        description: resistance_description,
        source_type: resistance_source_type,
        citation: resistance_citation,
        url: resistance_url,
        doi: resistance_doi,
        pmid: resistance_pmid,
        nct: resistance_nct
    }

    column_map_prognostic = {
        feature_display: prognostic_feature_display,
        oncotree_term: prognostic_oncotree_term,
        oncotree_code: prognostic_oncotree_code,
        context: prognostic_context,
        prognosis: prognosis,
        implication: prognostic_implication,
        implication_map: prognostic_implication_map,
        description: prognostic_description,
        source_type: prognostic_source_type,
        citation: prognostic_citation,
        url: prognostic_url,
        doi: prognostic_doi,
        pmid: prognostic_pmid,
        nct: prognostic_nct
    }

    score_bin = Annotator.score_bin
    bin_name = Annotator.almanac_bin
    sensitive_bin = Annotator.sensitive_bin
    resistance_bin = Annotator.resistance_bin
    prognostic_bin = Annotator.prognostic_bin
    sensitive_score_bin = Annotator.sensitive_score_bin
    resistance_score_bin = Annotator.resistance_score_bin
    prognostic_score_bin = Annotator.prognostic_score_bin

    columns = datasources.Almanac.columns
    query = datasources.Almanac.query
    genes = datasources.Almanac.genes

    matches = 'matches'
    sensitivity_matches_tablename = datasources.Almanac.sensitivity_matches_tablename
    resistance_matches_tablename = datasources.Almanac.resistance_matches_tablename
    prognostic_matches_tablename = datasources.Almanac.prognostic_matches_tablename

    Query = tinydb.Query()

    assertion_types_dict = {
        sensitivity: {
            query: (Query[sensitivity] == '1'),
            score_bin: sensitive_score_bin,
            columns: column_map_sensitive,
            feature: sensitive_feature,
            feature_type: sensitive_feature_type,
            alt_type: sensitive_alt_type,
            alt: sensitive_alt,
            matches: sensitivity_matches_tablename},
        resistance: {
            query: (Query[resistance] == '1'),
            score_bin: resistance_score_bin,
            columns: column_map_resistance,
            feature: resistance_feature,
            feature_type: resistance_feature_type,
            alt_type: resistance_alt_type,
            alt: resistance_alt,
            matches: resistance_matches_tablename},
        prognosis: {
            query: (Query[prognosis].one_of(['0', '1'])),
            score_bin: prognostic_score_bin,
            columns: column_map_prognostic,
            feature: prognostic_feature,
            feature_type: prognostic_feature_type,
            alt_type: prognostic_alt_type,
            alt: prognostic_alt,
            matches: prognostic_matches_tablename}
    }

    feature_types_section = 'feature_types'
    feature_types_config = CONFIG[feature_types_section]
    aneuploidy = feature_types_config['aneuploidy']
    burden = feature_types_config['burden']
    copynumber_variant = feature_types_config['cna']
    fusion = feature_types_config['fusion']
    germline_variant = feature_types_config['germline']
    microsatellite_status = feature_types_config['microsatellite']
    signature = feature_types_config['signature']
    somatic_variant = feature_types_config['mut']

    @classmethod
    def annotate(cls, df, dbs, ontology):
        ds = datasources.Almanac.import_ds(dbs)
        list_genes = ds.table(cls.genes).all()[0][cls.genes]

        annotation_function_dict = {
            cls.aneuploidy: cls.annotate_aneuploidy,
            cls.burden: cls.annotate_burden,
            cls.copynumber_variant: cls.annotate_copy_number,
            cls.fusion: cls.annotate_fusion,
            cls.germline_variant: cls.annotate_variants,
            cls.microsatellite_status: cls.annotate_microsatellite_stability,
            cls.signature: cls.annotate_signatures,
            cls.somatic_variant: cls.annotate_variants
        }

        for feature_type, group in df.groupby(cls.feature_type):
            table = ds.table(feature_type)
            for key, value in cls.predictive_implication_map.items():
                table.update({cls.implication_map: value}, (cls.Query[cls.implication] == key))

            if feature_type in [cls.somatic_variant, cls.germline_variant, cls.copynumber_variant, cls.fusion]:
                idx = group[cls.feature].isin(list_genes)
                df.loc[group[~idx].index, cls.bin_name] = 0
                group = group[group[cls.feature].isin(list_genes)]

            for index in group.index:
                df.loc[index, :] = annotation_function_dict[feature_type](df.loc[index, :], ontology, table)
        return df

    @classmethod
    def annotate_aneuploidy(cls, sliced_series, ontology, table):
        series = sliced_series.copy(deep=True)
        feature = series.loc[cls.feature]
        feature_type = series.loc[cls.feature_type]

        query_same_ontology = (cls.Query[cls.oncotree_code] == ontology)
        query_diff_ontology = (cls.Query[cls.oncotree_code] != ontology)
        query_feature = (cls.Query[cls.event] == feature)

        match_bins = []

        for assertion_type in cls.assertion_types_dict.keys():
            query = cls.assertion_types_dict[assertion_type][cls.query]
            score_bin = cls.assertion_types_dict[assertion_type][cls.score_bin]
            columns = copy.deepcopy(cls.assertion_types_dict[assertion_type][cls.columns])
            columns[cls.feature_type] = cls.assertion_types_dict[assertion_type][cls.feature_type]
            columns[cls.event] = cls.assertion_types_dict[assertion_type][cls.feature]

            if table.contains(query_feature & query):
                results_same_ontology = table.search(query_same_ontology & query_feature & query)
                results_diff_ontology = table.search(query_diff_ontology & query_feature & query)
                feature_match_to_assertion_bin = 3

            else:
                feature_match_to_assertion_bin = 0
                match_bins.append(feature_match_to_assertion_bin)
                continue

            matches = cls.sort_and_subset_matches(results_same_ontology, results_diff_ontology)
            series = cls.update_series_with_best_match(matches, feature_type, columns, series)
            series.loc[score_bin] = feature_match_to_assertion_bin
            match_bins.append(feature_match_to_assertion_bin)

            series = cls.insert_additional_matches(feature_type, sliced_series.name, assertion_type, matches, series)

        series.loc[cls.bin_name] = max(match_bins)
        return series

    @classmethod
    def annotate_burden(cls, sliced_series, ontology, table):
        series = sliced_series.copy(deep=True)
        feature = series.loc[cls.feature].split(' ')[0]
        feature_type = series.loc[cls.feature_type]

        query_same_ontology = (cls.Query[cls.oncotree_code] == ontology)
        query_diff_ontology = (cls.Query[cls.oncotree_code] != ontology)
        query_feature = (cls.Query[cls.classification] == feature)

        match_bins = []

        for assertion_type in cls.assertion_types_dict.keys():
            query = cls.assertion_types_dict[assertion_type][cls.query]
            score_bin = cls.assertion_types_dict[assertion_type][cls.score_bin]
            columns = copy.deepcopy(cls.assertion_types_dict[assertion_type][cls.columns])
            columns[cls.feature_type] = cls.assertion_types_dict[assertion_type][cls.feature_type]
            columns[cls.classification] = cls.assertion_types_dict[assertion_type][cls.feature]

            if table.contains(query_feature & query):
                results_same_ontology = table.search(query_same_ontology & query_feature & query)
                results_diff_ontology = table.search(query_diff_ontology & query_feature & query)
                feature_match_to_assertion_bin = 3

            else:
                feature_match_to_assertion_bin = 0
                match_bins.append(feature_match_to_assertion_bin)
                continue

            matches = cls.sort_and_subset_matches(results_same_ontology, results_diff_ontology)
            series = cls.update_series_with_best_match(matches, feature_type, columns, series)
            series.loc[score_bin] = feature_match_to_assertion_bin
            match_bins.append(feature_match_to_assertion_bin)

            series = cls.insert_additional_matches(feature_type, sliced_series.name, assertion_type, matches, series)

        series.loc[cls.bin_name] = max(match_bins)
        return series

    @classmethod
    def annotate_copy_number(cls, sliced_series, ontology, table):
        series = sliced_series.copy(deep=True)
        feature = series.loc[cls.feature]
        feature_type = series.loc[cls.feature_type]
        alt_type = series.loc[cls.alt_type]
        alt = series.loc[Features.segment_mean]

        query_same_ontology = (cls.Query[cls.oncotree_code] == ontology)
        query_diff_ontology = (cls.Query[cls.oncotree_code] != ontology)
        query_feature = (cls.Query[cls.gene] == feature)
        query_alt_type = (cls.Query[cls.direction] == alt_type)

        query_to_alt_type = query_feature & query_alt_type

        match_bins = []

        for assertion_type in cls.assertion_types_dict.keys():
            query = cls.assertion_types_dict[assertion_type][cls.query]
            score_bin = cls.assertion_types_dict[assertion_type][cls.score_bin]
            columns = copy.deepcopy(cls.assertion_types_dict[assertion_type][cls.columns])
            columns[cls.feature_type] = cls.assertion_types_dict[assertion_type][cls.feature_type]
            columns[cls.gene] = cls.assertion_types_dict[assertion_type][cls.feature]
            columns[cls.direction] = cls.assertion_types_dict[assertion_type][cls.alt_type]

            if table.contains(query_to_alt_type & query):
                results_same_ontology = table.search(query_same_ontology & query_to_alt_type & query)
                results_diff_ontology = table.search(query_diff_ontology & query_to_alt_type & query)
                if abs(float(alt)) >= 1.0:
                    feature_match_to_assertion_bin = 4
                else:
                    feature_match_to_assertion_bin = 3

            elif table.contains(query_feature & query):
                results_same_ontology = table.search(query_same_ontology & query_feature & query)
                results_diff_ontology = table.search(query_diff_ontology & query_feature & query)
                feature_match_to_assertion_bin = 2

            else:
                feature_match_to_assertion_bin = 1
                match_bins.append(feature_match_to_assertion_bin)
                continue

            matches = cls.sort_and_subset_matches(results_same_ontology, results_diff_ontology)
            series = cls.update_series_with_best_match(matches, feature_type, columns, series)
            series.loc[score_bin] = feature_match_to_assertion_bin
            match_bins.append(feature_match_to_assertion_bin)

            series = cls.insert_additional_matches(feature_type, sliced_series.name, assertion_type, matches, series)

        series.loc[cls.bin_name] = max(match_bins)
        return series

    @classmethod
    def annotate_fusion(cls, sliced_series, ontology, table):
        series = sliced_series.fillna('').copy(deep=True)
        feature = series.loc[cls.feature]
        feature_type = series.loc[cls.feature_type]
        alt_type = series.loc[cls.alt_type]
        alt = series.loc[cls.alt]

        if '--' in str(alt):
            split_alt = alt.split('--')
            gene1 = split_alt[0]
            gene2 = split_alt[1]
            # gene1 = feature
        else:
            gene1 = feature
            gene2 = ''

        query_same_ontology = (cls.Query[cls.oncotree_code] == ontology)
        query_diff_ontology = (cls.Query[cls.oncotree_code] != ontology)
        query_alt_type = (cls.Query[cls.rearrangement_type] == alt_type)

        match_bins = []

        for assertion_type in cls.assertion_types_dict.keys():
            query = cls.assertion_types_dict[assertion_type][cls.query]
            score_bin = cls.assertion_types_dict[assertion_type][cls.score_bin]
            columns = copy.deepcopy(cls.assertion_types_dict[assertion_type][cls.columns])
            columns[cls.feature_type] = cls.assertion_types_dict[assertion_type][cls.feature_type]
            columns[cls.rearrangement_type] = cls.assertion_types_dict[assertion_type][cls.alt_type]
            alt_column = cls.assertion_types_dict[assertion_type][cls.alt]

            fusion_matches = []
            for first_gene, second_gene, feature_col in [(gene1, gene2, cls.gene1), (gene2, gene1, cls.gene2)]:
                columns[feature_col] = cls.assertion_types_dict[assertion_type][cls.feature]

                query_first_gene = (cls.Query[cls.gene1] == first_gene)
                query_second_gene = (cls.Query[cls.gene2] == second_gene)
                query_first_gene_reverse = (cls.Query[cls.gene1] == second_gene)
                query_second_gene_reverse = (cls.Query[cls.gene2] == first_gene)

                query_to_alt = query_first_gene & query_second_gene & query_alt_type
                query_to_alt_reverse = query_first_gene_reverse & query_second_gene_reverse & query_alt_type

                query_to_alt_type = query_first_gene & query_alt_type
                query_to_alt_type_reverse = query_first_gene_reverse & query_alt_type

                if table.contains(query_to_alt & query):
                    results_same_ontology = table.search(query_same_ontology & query_to_alt & query)
                    results_diff_ontology = table.search(query_diff_ontology & query_to_alt & query)
                    match_bin = 4

                elif table.contains(query_to_alt_reverse & query):
                    results_same_ontology = table.search(query_same_ontology & query_to_alt_reverse & query)
                    results_diff_ontology = table.search(query_diff_ontology & query_to_alt_reverse & query)
                    match_bin = 4

                elif table.contains(query_to_alt_type & query):
                    results_same_ontology = table.search(query_same_ontology & query_to_alt_type & query)
                    results_diff_ontology = table.search(query_diff_ontology & query_to_alt_type & query)
                    match_bin = 3

                elif table.contains(query_to_alt_type_reverse & query):
                    results_same_ontology = table.search(query_same_ontology & query_to_alt_type_reverse & query)
                    results_diff_ontology = table.search(query_diff_ontology & query_to_alt_type_reverse & query)
                    match_bin = 3

                elif table.contains(query_first_gene & query):
                    results_same_ontology = table.search(query_same_ontology & query_first_gene & query)
                    results_diff_ontology = table.search(query_diff_ontology & query_first_gene & query)
                    match_bin = 2

                elif table.contains(query_second_gene & query):
                    results_same_ontology = table.search(query_same_ontology & query_second_gene & query)
                    results_diff_ontology = table.search(query_diff_ontology & query_second_gene & query)
                    match_bin = 2

                else:
                    results_same_ontology = []
                    results_diff_ontology = []
                    match_bin = 1

                tmp_dict = {'bin': match_bin,
                            'results_same_ontology': results_same_ontology,
                            'results_diff_ontology': results_diff_ontology,
                            'alt': '{}--{}'.format(gene1, gene2),
                            'feature_col': feature_col}
                fusion_matches.append(tmp_dict)

            better_match = cls.select_better_fusion_match(fusion_matches)
            feature_match_to_assertion_bin = better_match['bin']
            if feature_match_to_assertion_bin == 1:
                match_bins.append(feature_match_to_assertion_bin)
                continue

            results_same_ontology = better_match['results_same_ontology']
            results_diff_ontology = better_match['results_diff_ontology']
            matches = cls.sort_and_subset_matches(results_same_ontology, results_diff_ontology)

            feature_col = better_match['feature_col']
            columns[feature_col] = cls.assertion_types_dict[assertion_type][cls.feature]
            series = cls.update_series_with_best_match(matches, feature_type, columns, series)
            series.loc[alt_column] = better_match['alt']
            series.loc[score_bin] = feature_match_to_assertion_bin
            match_bins.append(feature_match_to_assertion_bin)

            series = cls.insert_additional_matches(feature_type, sliced_series.name, assertion_type, matches, series)
        series.loc[cls.bin_name] = max(match_bins)
        return series

    @classmethod
    def annotate_microsatellite_stability(cls, sliced_series, ontology, table):
        series = sliced_series.copy(deep=True)
        feature = series.loc[cls.feature].split(' ')[-1]
        feature_type = series.loc[cls.feature_type]

        query_same_ontology = (cls.Query[cls.oncotree_code] == ontology)
        query_diff_ontology = (cls.Query[cls.oncotree_code] != ontology)
        query_feature = (cls.Query[cls.status] == feature)

        match_bins = []

        for assertion_type in cls.assertion_types_dict.keys():
            query = cls.assertion_types_dict[assertion_type][cls.query]
            score_bin = cls.assertion_types_dict[assertion_type][cls.score_bin]
            columns = copy.deepcopy(cls.assertion_types_dict[assertion_type][cls.columns])
            columns[cls.feature_type] = cls.assertion_types_dict[assertion_type][cls.feature_type]
            columns[cls.status] = cls.assertion_types_dict[assertion_type][cls.feature]

            if table.contains(query_feature & query):
                results_same_ontology = table.search(query_same_ontology & query_feature & query)
                results_diff_ontology = table.search(query_diff_ontology & query_feature & query)
                feature_match_to_assertion_bin = 3

            else:
                feature_match_to_assertion_bin = 0
                match_bins.append(feature_match_to_assertion_bin)
                continue

            matches = cls.sort_and_subset_matches(results_same_ontology, results_diff_ontology)
            series = cls.update_series_with_best_match(matches, feature_type, columns, series)
            series.loc[score_bin] = feature_match_to_assertion_bin
            match_bins.append(feature_match_to_assertion_bin)

            series = cls.insert_additional_matches(feature_type, sliced_series.name, assertion_type, matches, series)

        series.loc[cls.bin_name] = max(match_bins)
        return series

    @classmethod
    def annotate_signatures(cls, sliced_series, ontology, table):
        series = sliced_series.copy(deep=True)
        feature = series.loc[cls.feature].split(' ')[-1]
        feature_type = series.loc[cls.feature_type]

        query_same_ontology = (cls.Query[cls.oncotree_code] == ontology)
        query_diff_ontology = (cls.Query[cls.oncotree_code] != ontology)
        query_feature = (cls.Query[cls.cosmic_signature_number] == feature)

        match_bins = []

        for assertion_type in cls.assertion_types_dict.keys():
            query = cls.assertion_types_dict[assertion_type][cls.query]
            score_bin = cls.assertion_types_dict[assertion_type][cls.score_bin]
            columns = copy.deepcopy(cls.assertion_types_dict[assertion_type][cls.columns])
            columns[cls.feature_type] = cls.assertion_types_dict[assertion_type][cls.feature_type]
            columns[cls.cosmic_signature_number] = cls.assertion_types_dict[assertion_type][cls.feature]

            if table.contains(query_feature & query):
                results_same_ontology = table.search(query_same_ontology & query_feature & query)
                results_diff_ontology = table.search(query_diff_ontology & query_feature & query)
                feature_match_to_assertion_bin = 3

            else:
                feature_match_to_assertion_bin = 1
                match_bins.append(feature_match_to_assertion_bin)
                continue

            matches = cls.sort_and_subset_matches(results_same_ontology, results_diff_ontology)
            series = cls.update_series_with_best_match(matches, feature_type, columns, series)
            series.loc[score_bin] = feature_match_to_assertion_bin
            match_bins.append(feature_match_to_assertion_bin)

            series = cls.insert_additional_matches(feature_type, sliced_series.name, assertion_type, matches, series)

        series.loc[cls.bin_name] = max(match_bins)
        return series

    @classmethod
    def annotate_variants(cls, sliced_series, ontology, table):
        series = sliced_series.copy(deep=True)
        feature = series.loc[cls.feature]
        feature_type = series.loc[cls.feature_type]
        alt_type = series.loc[cls.alt_type]
        alt = series.loc[cls.alt]

        query_same_ontology = (cls.Query[cls.oncotree_code] == ontology)
        query_diff_ontology = (cls.Query[cls.oncotree_code] != ontology)
        query_feature = (cls.Query[cls.gene] == feature)
        query_alt_type = (cls.Query[cls.variant_annotation] == alt_type)
        query_alt = (cls.Query[cls.protein_change] == alt)

        query_to_alt = query_feature & query_alt_type & query_alt
        query_to_alt_type = query_feature & query_alt_type

        match_bins = []

        for assertion_type in cls.assertion_types_dict.keys():
            query = cls.assertion_types_dict[assertion_type][cls.query]
            score_bin = cls.assertion_types_dict[assertion_type][cls.score_bin]
            columns = copy.deepcopy(cls.assertion_types_dict[assertion_type][cls.columns])
            columns[cls.feature_type] = cls.assertion_types_dict[assertion_type][cls.feature_type]
            columns[cls.gene] = cls.assertion_types_dict[assertion_type][cls.feature]
            columns[cls.variant_annotation] = cls.assertion_types_dict[assertion_type][cls.alt_type]
            columns[cls.protein_change] = cls.assertion_types_dict[assertion_type][cls.alt]

            if table.contains(query_to_alt & query):
                results_same_ontology = table.search(query_same_ontology & query_to_alt & query)
                results_diff_ontology = table.search(query_diff_ontology & query_to_alt & query)
                feature_match_to_assertion_bin = 4

            elif table.contains(query_to_alt_type & query):
                results_same_ontology = table.search(query_same_ontology & query_to_alt_type & query)
                results_diff_ontology = table.search(query_diff_ontology & query_to_alt_type & query)
                feature_match_to_assertion_bin = 3

            elif table.contains(query_feature & query):
                results_same_ontology = table.search(query_same_ontology & query_feature & query)
                results_diff_ontology = table.search(query_diff_ontology & query_feature & query)
                feature_match_to_assertion_bin = 2

            else:
                feature_match_to_assertion_bin = 1
                match_bins.append(feature_match_to_assertion_bin)
                continue

            matches = cls.sort_and_subset_matches(results_same_ontology, results_diff_ontology)
            series = cls.update_series_with_best_match(matches, feature_type, columns, series)
            series.loc[score_bin] = feature_match_to_assertion_bin
            match_bins.append(feature_match_to_assertion_bin)

            series = cls.insert_additional_matches(feature_type, sliced_series.name, assertion_type, matches, series)

        series.loc[cls.bin_name] = max(match_bins)
        return series

    @staticmethod
    def extract_values_from_list_of_dicts(key, list_of_dicts):
        return [x[key] for x in list_of_dicts]

    @staticmethod
    def initialize_matches_db(table_name):
        handle = CONFIG['databases']['additional_matches']
        db = tinydb.TinyDB(handle)
        return db.table(table_name)

    @classmethod
    def insert_additional_matches(cls, feature_type, index, assertion_type, matches, series):
        table_name = '{}-{}-{}'.format(feature_type, index, assertion_type)
        assertion_table = cls.initialize_matches_db(table_name)
        if len(matches) > 1:
            cls.insert_documents(assertion_table, matches[1:])
        else:
            cls.insert_documents(assertion_table, [])
        column = cls.assertion_types_dict[assertion_type][cls.matches]
        series.loc[column] = table_name
        return series

    @staticmethod
    def insert_documents(table, documents):
        table.insert_multiple(documents)

    @staticmethod
    def query_by_search(db, query):
        return db.search(query)

    @classmethod
    def select_better_fusion_match(cls, list_of_matches):
        max_bin = max(match['bin'] for match in list_of_matches)
        better_match = [match for match in list_of_matches if (match['bin'] == max_bin)]
        return better_match[0]

    @classmethod
    def sort_and_subset_matches(cls, results_same_ontology, results_diff_ontology):
        if results_same_ontology:
            sorted_results_same_ontology = cls.sort_dictionary_as_dataframe(results_same_ontology,
                                                                            [cls.implication_map,
                                                                             cls.preferred_assertion],
                                                                            [False, False])
        else:
            sorted_results_same_ontology = []

        if results_diff_ontology:
            sorted_results_diff_ontology = cls.sort_dictionary_as_dataframe(results_diff_ontology,
                                                                            [cls.implication_map,
                                                                             cls.preferred_assertion],
                                                                            [False, False])
        else:
            sorted_results_diff_ontology = []

        best_evidence = cls.return_best_evidence_level(results_same_ontology, results_diff_ontology)

        matches = []
        for results in [sorted_results_same_ontology, sorted_results_diff_ontology]:
            subsetted_results = cls.subset_list_of_dictionaries(results, cls.implication_map, op.ge, best_evidence)
            matches.extend(subsetted_results)
        return matches

    @staticmethod
    def sort_dictionary_as_dataframe(dictionary, sort_columns, ascending_boolean):
        return (pd.DataFrame(dictionary)
                .sort_values(sort_columns, ascending=ascending_boolean)
                ).to_dict(orient='records')

    @staticmethod
    def subset_dataframe_by_condition(dataframe, condition):
        return dataframe[condition]

    @staticmethod
    def subset_list_of_dictionaries(list_of_dicts, key, operator, subset_value):
        return [dictionary for dictionary in list_of_dicts if operator(dictionary[key], subset_value)]

    @classmethod
    def return_best_evidence_level(cls, results_same_ontology, results_diff_ontology):
        evidence_same_ontology = cls.extract_values_from_list_of_dicts(cls.implication_map, results_same_ontology)
        evidence_diff_ontology = cls.extract_values_from_list_of_dicts(cls.implication_map, results_diff_ontology)

        if evidence_same_ontology:
            best_evidence = max(evidence_same_ontology)
        else:
            best_evidence = max(evidence_diff_ontology)
        return best_evidence

    @classmethod
    def update_series_with_best_match(cls, matches, feature_type, columns, series):
        best_match = matches[0]
        best_match[cls.feature_type] = feature_type
        for column, assertion_column in columns.items():
            series.loc[assertion_column] = best_match[column]
        return series


class CancerHotspots(object):
    gene = datasources.CancerHotspots.gene
    alteration = datasources.CancerHotspots.alt

    bin_name = Annotator.cancerhotspots_bin
    comparison_columns = [gene, alteration]

    @classmethod
    def annotate(cls, df, dbs):
        return Annotator.annotate(df, dbs, datasources.CancerHotspots, cls.bin_name, cls.comparison_columns)


class CancerHotspots3D(object):
    gene = datasources.CancerHotspots3D.gene
    alteration = datasources.CancerHotspots3D.alt

    bin_name = Annotator.cancerhotspots3d_bin
    comparison_columns = [gene, alteration]

    @classmethod
    def annotate(cls, df, dbs):
        return Annotator.annotate(df, dbs, datasources.CancerHotspots3D, cls.bin_name, cls.comparison_columns)


class CancerGeneCensus(object):
    gene = datasources.CancerGeneCensus.gene

    bin_name = Annotator.cgc_bin
    comparison_columns = [gene]

    @classmethod
    def annotate(cls, df, dbs):
        return Annotator.annotate(df, dbs, datasources.CancerGeneCensus, cls.bin_name, cls.comparison_columns)


class ClinVar(object):
    chr = datasources.ClinVar.chr
    start = datasources.ClinVar.start
    end = datasources.ClinVar.end

    bin_name = Annotator.clinvar_bin

    @classmethod
    def append_clinvar(cls, df, ds):
        position_columns = [cls.start, cls.end]
        for col in position_columns:
            df.loc[:, col] = pd.to_numeric(df.loc[:, col])
            ds.loc[:, col] = pd.to_numeric(ds.loc[:, col])

        df[cls.chr] = df[cls.chr].astype(str)
        ds[cls.chr] = ds[cls.chr].astype(str)

        return df.merge(ds, how='left')

    @classmethod
    def annotate(cls, df, dbs):
        df.drop(df.columns[df.columns.str.contains('clinvar')], axis=1, inplace=True)
        ds = datasources.ClinVar.import_ds(dbs)
        df = cls.append_clinvar(df, ds)
        return Features.preallocate_missing_columns(df)


class Cosmic(object):
    gene = datasources.Cosmic.gene
    alteration = datasources.Cosmic.alt

    bin_name = Annotator.cosmic_bin
    comparison_columns = [gene, alteration]

    @classmethod
    def annotate(cls, df, dbs):
        return Annotator.annotate(df, dbs, datasources.Cosmic, cls.bin_name, cls.comparison_columns)


class ExAC:
    chr = datasources.ExAC.chr
    start = datasources.ExAC.start
    ref = datasources.ExAC.ref
    alt = datasources.ExAC.alt
    af = datasources.ExAC.af

    bin_name = Annotator.exac_common_bin
    exac_common_threshold = EXAC_CONFIG['exac_common_af_threshold']

    str_columns = [chr, ref, alt]
    int_columns = [start]

    somatic = CONFIG['feature_types']['mut']
    germline = CONFIG['feature_types']['germline']
    feature_type = Features.feature_type

    @classmethod
    def append_exac_af(cls, df, ds):
        variants, not_variants = cls.subset_for_variants(df)
        ds = ds.loc[:, [cls.chr, cls.start, cls.ref, cls.alt, cls.af]]

        for column, data_type in [(cls.str_columns, str), (cls.int_columns, int)]:
            variants.loc[variants.index, column] = cls.format_columns(variants, column, data_type)
            ds.loc[ds.index, column] = cls.format_columns(ds, column, data_type)

        merged = variants.merge(ds, how='left')
        merged.loc[merged.index, cls.af] = cls.fill_na(merged, cls.af, 0.0, float, 6)
        not_variants.loc[not_variants.index, cls.af] = cls.fill_na(not_variants, cls.af, 0.0, float, 6)
        return pd.concat([merged, not_variants]).sort_index()

    @classmethod
    def annotate(cls, df, dbs):
        df_dropped = cls.drop_existing_columns(df)
        ds = datasources.ExAC.import_ds(dbs)
        df_annotated = cls.append_exac_af(df_dropped, ds)
        df_annotated[cls.bin_name] = cls.annotate_common_af(df_annotated[cls.af])
        return Features.preallocate_missing_columns(df_annotated)

    @classmethod
    def annotate_common_af(cls, series_exac_af):
        if not series_exac_af.empty:
            series = pd.Series(float(0.0), index=series_exac_af.index.tolist())
            condition = (series_exac_af
                         .astype(str).str.split(',', expand=True)
                         .fillna(value=np.nan)
                         .astype(float).mean(axis=1)
                         .fillna(0.0)
                         )
            idx = condition.astype(float) >= float(cls.exac_common_threshold)
            series[idx] = float(1.0)
            return series
        else:
            return pd.Series()

    @classmethod
    def drop_existing_columns(cls, dataframe):
        return dataframe.drop(dataframe.columns[dataframe.columns.str.contains('exac')], axis=1)

    @classmethod
    def fill_na(cls, dataframe, column, fill_value, fill_data_type, round_places):
        if column not in dataframe.columns:
            dataframe.loc[dataframe.index, column] = pd.NA
        return dataframe.loc[dataframe.index, column].fillna(fill_value).astype(fill_data_type).round(round_places)

    @classmethod
    def format_columns(cls, dataframe, column, data_type):
        return dataframe.loc[dataframe.index, column].astype(data_type)

    @classmethod
    def subset_for_variants(cls, dataframe):
        idx = dataframe[cls.feature_type].isin([cls.somatic, cls.germline])
        return dataframe[idx].copy(), dataframe[~idx].copy()


class ExACExtended:
    @classmethod
    def annotate(cls, df, dbs):
        df_dropped = ExAC.drop_existing_columns(df)
        ds = datasources.ExACExtended.import_ds(dbs)
        df_annotated = ExAC.append_exac_af(df_dropped, ds)
        df_annotated[ExAC.bin_name] = ExAC.annotate_common_af(df_annotated[ExAC.af])
        return Features.preallocate_missing_columns(df_annotated)


class GSEACancerModules(object):
    gene = datasources.GSEACancerModules.gene

    bin_name = Annotator.gsea_module_bin
    comparison_columns = [gene]

    @classmethod
    def annotate(cls, df, dbs):
        return Annotator.annotate(df, dbs, datasources.GSEACancerModules, cls.bin_name, cls.comparison_columns)


class GSEACancerPathways(object):
    gene = datasources.GSEACancerPathways.gene

    bin_name = Annotator.gsea_pathway_bin
    comparison_columns = [gene]

    @classmethod
    def annotate(cls, df, dbs):
        return Annotator.annotate(df, dbs, datasources.GSEACancerPathways, cls.bin_name, cls.comparison_columns)


class Hereditary(object):
    gene = datasources.Hereditary.gene

    bin_name = Annotator.hereditary_bin
    comparison_columns = [gene]

    @classmethod
    def annotate(cls, df, dbs):
        return Annotator.annotate(df, dbs, datasources.Hereditary, cls.bin_name, cls.comparison_columns)


class MSI(object):
    gene = datasources.Datasources.feature

    bin_name = Annotator.msi_bin
    comparison_columns = [gene]

    msi_genes = ['MSH2', 'MSH6', 'MLH1', 'PMS2', 'POLE', 'POLE2',
                 'ACVR2A', 'RNF43', 'JAK1', 'MSH3',
                 'ESRP1', 'PRDM2', 'DOCK3']

    @classmethod
    def create_msi_df(cls):
        return pd.DataFrame(cls.msi_genes, columns=[cls.gene])

    @classmethod
    def annotate(cls, df):
        df[cls.bin_name] = Annotator.match_ds(df, cls.create_msi_df(), cls.bin_name, cls.comparison_columns)
        return df


class OverlapValidation(object):
    section = 'validation_sequencing'
    gene = COLNAMES[section]['gene']
    feature_type = COLNAMES[section]['feature_type']
    alt_type = COLNAMES[section]['alt_type']
    alt = COLNAMES[section]['alt']
    tumor_f = COLNAMES[section]['tumor_f']
    coverage = COLNAMES[section]['coverage']
    validation_tumor_f = COLNAMES[section]['validation_tumor_f']
    validation_coverage = COLNAMES[section]['validation_coverage']
    validation_detection_power = COLNAMES[section]['validation_detection_power']

    somatic_variants = CONFIG['feature_types']['mut']

    merge_cols = [gene, alt_type, alt]
    fill_cols = [tumor_f, validation_tumor_f, validation_coverage]

    @classmethod
    def append_validation(cls, primary, validation):
        df = cls.drop_validation_columns(primary)
        df = cls.merge_data_frames(df, validation, cls.merge_cols)
        idx = cls.get_mutation_index(df)
        df.loc[idx, cls.fill_cols] = df.loc[idx, cls.fill_cols].fillna(0.0)

        df.loc[idx, cls.validation_detection_power] = cls.calculate_validation_detection_power(
            df.loc[idx, cls.tumor_f].astype(float),
            df.loc[idx, cls.coverage].astype(float),
            df.loc[idx, cls.validation_coverage].astype(float)
        )

        df[cls.validation_detection_power] = cls.round_series(df[cls.validation_detection_power], 4)
        return df

    @staticmethod
    def calculate_beta_binomial(k, n, a, b):
        return stats.betabinom.sf(k=k, n=n, a=a, b=b, loc=0)

    @classmethod
    def calculate_validation_detection_power(cls, primary_tumor_f, primary_coverage, validation_coverage):
        alt_count = primary_tumor_f.multiply(primary_coverage)
        ref_count = primary_coverage.subtract(alt_count)
        return cls.calculate_beta_binomial(k=3, n=validation_coverage, a=alt_count.add(1), b=ref_count.add(1))

    @classmethod
    def drop_validation_columns(cls, df):
        return df.drop([cls.validation_tumor_f, cls.validation_coverage], axis=1)

    @classmethod
    def get_mutation_index(cls, df):
        return df[df[cls.feature_type].eq(cls.somatic_variants)].index

    @classmethod
    def merge_data_frames(cls, df1, df2, columns):
        return df1.merge(df2, on=columns, how='left')

    @staticmethod
    def round_series(series, n):
        series = pd.to_numeric(series, errors='coerce')
        return series.round(n)


class OverlapSomaticGermline(object):
    section = 'overlap_somatic_germline'
    gene = COLNAMES[section]['gene']
    alt_type = COLNAMES[section]['alt_type']
    number_germline_mutations = COLNAMES[section]['number_germline_mutations']

    @classmethod
    def count_germline_hits(cls, df):
        variants = df[df[cls.alt_type].isin(Features.coding_classifications_mapped)]
        value_counts = variants[cls.gene].value_counts()
        return cls.value_counts_to_df(value_counts)

    @classmethod
    def merge_germline_hits(cls, somatic, germline_counts):
        return pd.merge(somatic.drop([cls.number_germline_mutations], axis=1), germline_counts,
                        on=[cls.gene], how='left')

    @classmethod
    def append_germline_hits(cls, somatic, germline):
        germline_counts = cls.count_germline_hits(germline)
        merged_df = cls.merge_germline_hits(somatic, germline_counts)
        return merged_df

    @classmethod
    def value_counts_to_df(cls, value_counts):
        return pd.DataFrame({cls.gene: value_counts.index, cls.number_germline_mutations: value_counts.values})


class PreclinicalEfficacy:
    section = 'preclinical'
    feature_display = COLNAMES[section]['feature_display']
    pvalue = COLNAMES[section]['pvalue']
    efficacy = COLNAMES[section]['efficacy_obs']

    @classmethod
    def annotate(cls, actionable, efficacy):
        series_all_features = actionable[cls.feature_display]
        series_features = series_all_features[series_all_features.isin(efficacy[cls.feature_display])]
        for index in series_features.index:
            feature = series_features.loc[index]
            dataframe = efficacy[efficacy[cls.feature_display].eq(feature)]
            efficacy_observed = cls.search_for_significance(dataframe[cls.pvalue])
            actionable.loc[index, cls.efficacy] = efficacy_observed
        actionable[cls.efficacy].fillna(pd.NA, inplace=True)
        return actionable

    @classmethod
    def search_for_significance(cls, series):
        test = series.dropna().astype(float).le(0.05).tolist()
        if True in test:
            return 1
        else:
            return 0


class PreclinicalMatchmaking:
    section = 'matchmaking'
    feature_type = COLNAMES[section]['feature_type']
    feature = COLNAMES[section]['feature']
    alteration_type = COLNAMES[section]['alteration_type']
    alteration = COLNAMES[section]['alteration']
    partner = COLNAMES[section]['partner']
    gene = COLNAMES[section]['gene']
    gene1 = COLNAMES[section]['gene1']
    gene2 = COLNAMES[section]['gene2']
    genes = COLNAMES[section]['genes']
    direction = COLNAMES[section]['direction']
    rearrangement_type = COLNAMES[section]['rearrangement_type']
    evidence = COLNAMES[section]['evidence']
    evidence_map_str = COLNAMES[section]['evidence_map']
    group1 = COLNAMES[section]['group1']
    group2 = COLNAMES[section]['group2']
    group3 = COLNAMES[section]['group3']
    group4 = COLNAMES[section]['group4']
    feature_match = COLNAMES[section]['feature_match']
    match_1 = COLNAMES[section]['feature_match_1']
    match_2 = COLNAMES[section]['feature_match_2']
    match_3 = COLNAMES[section]['feature_match_3']
    match_4 = COLNAMES[section]['feature_match_4']
    model_id = COLNAMES[section]['model_id']
    variant_annotation = COLNAMES[section]['variant_annotation']
    protein_change = COLNAMES[section]['protein_change']
    merged = COLNAMES[section]['merged']
    feature_display = COLNAMES[section]['feature_display']
    predictive_implication = COLNAMES[section]['predictive_implication']

    feature_types_section = 'feature_types'
    feature_types_config = CONFIG[feature_types_section]
    copy_number = feature_types_config['cna']
    fusion = feature_types_config['fusion']
    somatic_variant = feature_types_config['mut']

    evidence_map = {
        'FDA-Approved': 5.0, 'Guideline': 4.0, 'Clinical trial': 3.0,
        'Clinical evidence': 2.0, 'Preclinical': 1.0, 'Inferential': 0.0}
    inv_evidence_map = {float(v): k for k, v in evidence_map.items()}

    variants = 'variants'
    cnas = 'copy_number_alterations'
    fusions = 'fusions'
    fusions_gene1 = 'fusions_gene1'
    fusions_gene2 = 'fusions_gene2'

    @classmethod
    def annotate(cls, input_dict, dbs):
        input_variants = input_dict[cls.somatic_variant]
        input_copy_number_alterations = input_dict[cls.copy_number]
        input_fusions = input_dict[cls.fusion]

        variants = cls.annotate_somatic_variants(input_variants, dbs)
        copy_number_alterations = cls.annotate_copy_numbers(input_copy_number_alterations, dbs)
        fusions, fusions_gene1, fusions_gene2 = cls.annotate_fusions(input_fusions, dbs)
        return {
            cls.variants: variants,
            cls.cnas: copy_number_alterations,
            cls.fusions: fusions,
            cls.fusions_gene1: fusions_gene1,
            cls.fusions_gene2: fusions_gene2
        }

    @classmethod
    def annotate_copy_numbers(cls, df, dbs):
        almanac = datasources.Almanac.import_ds(dbs)
        almanac_genes = almanac.table(cls.genes).all()[0][cls.genes]

        df = df[df[cls.feature_type].eq(cls.copy_number)]
        db = pd.DataFrame(almanac.table(cls.copy_number).all())

        column_map = {cls.gene: cls.feature, cls.direction: cls.alteration_type}
        db = cls.format_db(list(column_map.keys()), column_map, db)
        df = cls.preallocate_almanac_matches(df)
        df = cls.annotate_match_1(df, almanac_genes)
        df = cls.annotate_match_2(df, db)
        df = cls.annotate_match_3(df, db)
        df = cls.annotate_other_datasources(df, dbs)
        return df

    @classmethod
    def annotate_fusions(cls, df, dbs):
        almanac = datasources.Almanac.import_ds(dbs)
        almanac_genes = almanac.table(cls.genes).all()[0][cls.genes]

        df = df[df[cls.feature_type].eq(cls.fusion)]
        db = pd.DataFrame(almanac.table(cls.fusion).all())

        column_map = {cls.rearrangement_type: cls.alteration_type}
        db = cls.format_db([cls.gene1, cls.gene2, cls.rearrangement_type], column_map, db)
        df = cls.preallocate_almanac_matches(df)

        df_gene1 = df.copy(deep=True)
        df_gene2 = df.copy(deep=True)
        df_gene1[cls.alteration] = df[cls.feature] + '--' + df[cls.partner]
        df_gene2[cls.alteration] = df[cls.partner] + '--' + df[cls.feature]
        df_gene2.rename(columns={cls.partner: cls.feature, cls.feature: cls.partner}, inplace=True)
        db_gene1 = db.copy(deep=True).rename(columns={cls.gene1: cls.feature, cls.gene2: cls.partner}).fillna('')
        db_gene2 = db.copy(deep=True).rename(columns={cls.gene1: cls.partner, cls.gene2: cls.feature}).fillna('')

        # Consider both genes
        group1 = cls.annotate_fusions_matching(df_gene1, db_gene1, almanac_genes, consider_partner=True)
        group2 = cls.annotate_fusions_matching(df_gene1, db_gene2, almanac_genes, consider_partner=True)
        group3 = cls.annotate_fusions_matching(df_gene2, db_gene1, almanac_genes, consider_partner=True)
        group4 = cls.annotate_fusions_matching(df_gene2, db_gene2, almanac_genes, consider_partner=True)

        values = pd.concat([
            group1.rename(columns={cls.evidence_map_str: cls.group1})[cls.group1],
            group2.rename(columns={cls.evidence_map_str: cls.group2})[cls.group2],
            group3.rename(columns={cls.evidence_map_str: cls.group3})[cls.group3],
            group4.rename(columns={cls.evidence_map_str: cls.group4})[cls.group4],
        ], axis=1)
        values = values.fillna(-1).idxmax(axis=1)

        idx_group1 = values[values.eq(cls.group1)].index
        idx_group2 = values[values.eq(cls.group2)].index
        idx_group3 = values[values.eq(cls.group3)].index
        idx_group4 = values[values.eq(cls.group4)].index

        columns = [cls.evidence, cls.match_1, cls.match_2, cls.match_3, cls.match_4, cls.feature_match]
        pairs = [(group1, idx_group1), (group2, idx_group2), (group3, idx_group3), (group4, idx_group4)]
        for group, index in pairs:
            df.loc[index, columns] = group.loc[index, columns]

        columns = [cls.feature, cls.partner]
        for index in df.index:
            df.loc[index, columns] = sorted(df.loc[index, columns].astype(str).tolist())

        df = (df
              .sort_values([cls.evidence, cls.feature_match], ascending=False)
              .drop_duplicates([cls.feature, cls.partner, cls.model_id], keep='first')
              )
        df[cls.alteration] = ''

        # Consider Gene1
        group1 = cls.annotate_fusions_matching(df_gene1, db_gene1, almanac_genes, False)
        group2 = cls.annotate_fusions_matching(df_gene1, db_gene2, almanac_genes, False)

        values = pd.concat([
            group1.rename(columns={cls.evidence_map_str: cls.group1})[cls.group1],
            group2.rename(columns={cls.evidence_map_str: cls.group2})[cls.group2],
        ], axis=1)
        values = values.fillna(-1).idxmax(axis=1)

        idx_group1 = values[values.eq(cls.group1)].index
        idx_group2 = values[values.eq(cls.group2)].index

        columns = [cls.evidence, cls.feature_match]
        group1.loc[idx_group2, columns] = group2.loc[idx_group2, columns]

        group1 = (group1
                  .sort_values([cls.evidence, cls.feature_match], ascending=False)
                  .drop_duplicates([cls.feature, cls.model_id], keep='first')
                  )
        group1[cls.alteration] = ''

        # Consider Gene2
        group3 = cls.annotate_fusions_matching(df_gene2, db_gene1, almanac_genes, False)
        group4 = cls.annotate_fusions_matching(df_gene2, db_gene2, almanac_genes, False)

        values = pd.concat([
            group3.rename(columns={cls.evidence_map_str: cls.group3})[cls.group3],
            group4.rename(columns={cls.evidence_map_str: cls.group4})[cls.group4],
        ], axis=1)
        values = values.fillna(-1).idxmax(axis=1)

        idx_group3 = values[values.eq(cls.group3)].index
        idx_group4 = values[values.eq(cls.group4)].index

        columns = [cls.evidence, cls.feature_match]
        group3.loc[idx_group4, columns] = group4.loc[idx_group4, columns]

        group3 = (group3
                  .sort_values([cls.evidence, cls.feature_match], ascending=False)
                  .drop_duplicates([cls.feature, cls.model_id], keep='first')
                  )

        group3[cls.alteration] = ''

        # Annotate with other datasources and return
        df = cls.annotate_other_datasources(df, dbs)
        group1 = cls.annotate_other_datasources(group1, dbs)
        group3 = cls.annotate_other_datasources(group3, dbs)
        return df, group1, group3

    @classmethod
    def annotate_fusions_matching(cls, df, db, db_genes, consider_partner=False):
        df = cls.annotate_match_1(df, db_genes)
        df = cls.annotate_match_2(df, db)
        df = cls.annotate_match_3(df, db)
        if consider_partner:
            df = cls.annotate_match_4(df, db, alteration=cls.partner)
        feature_match_columns = [cls.match_1, cls.match_2, cls.match_3, cls.match_4]
        df[cls.feature_match] = df.loc[:, feature_match_columns].sum(axis=1)
        df[cls.evidence_map_str] = df[cls.evidence].copy(deep=True)
        df[cls.evidence] = df[cls.evidence].fillna('').replace(cls.inv_evidence_map)
        return df

    @classmethod
    def annotate_somatic_variants(cls, df, dbs):
        almanac = datasources.Almanac.import_ds(dbs)
        almanac_genes = almanac.table(cls.genes).all()[0][cls.genes]

        df = df[df[cls.feature_type].eq(cls.somatic_variant)]
        db = pd.DataFrame(almanac.table(cls.somatic_variant).all())
        db[cls.variant_annotation].replace({'Oncogenic Mutations': '', 'Activating mutation': ''}, inplace=True)

        column_map = {cls.gene: cls.feature,
                      cls.variant_annotation: cls.alteration_type,
                      cls.protein_change: cls.alteration}
        db = cls.format_db(list(column_map.keys()), column_map, db)
        df = cls.preallocate_almanac_matches(df)
        df = cls.annotate_match_1(df, almanac_genes)
        df = cls.annotate_match_2(df, db)
        df = cls.annotate_match_3(df, db)
        df = cls.annotate_match_4(df, db)
        df = cls.annotate_other_datasources(df, dbs)
        return df

    @classmethod
    def annotate_match_1(cls, df, genes):
        idx_match = df[cls.feature].isin(genes)
        df.loc[idx_match, cls.match_1] = 1
        return df

    @classmethod
    def annotate_match_2(cls, df, db, feature=feature):
        match_columns = [feature]
        df = (df
              .merge(db
                     .loc[:, match_columns + [cls.evidence_map_str, cls.merged]]
                     .sort_values(cls.evidence_map_str, ascending=False)
                     .drop_duplicates(match_columns, keep='first'),
                     on=match_columns,
                     how='left'
                     )
              )
        index_match = df[cls.merged].eq(1)
        df.loc[index_match, cls.match_2] = 1
        df.loc[index_match, cls.evidence] = df.loc[index_match, cls.evidence_map_str]
        df.drop([cls.merged, cls.evidence_map_str], axis=1, inplace=True)
        return df

    @classmethod
    def annotate_match_3(cls, df, db, feature=feature, alteration_type=alteration_type):
        match_columns = [feature, alteration_type]
        df = (df
              .merge(db[~db[alteration_type].eq('')]
                     .loc[:, match_columns + [cls.evidence_map_str, cls.merged]]
                     .sort_values(cls.evidence_map_str, ascending=False)
                     .drop_duplicates(match_columns, keep='first'),
                     on=match_columns,
                     how='left'
                     )
              )
        index_match = df[cls.merged].eq(1)
        df.loc[index_match, cls.match_3] = 1
        df.loc[index_match, cls.evidence] = df.loc[index_match, cls.evidence_map_str]
        df.drop([cls.merged, cls.evidence_map_str], axis=1, inplace=True)
        return df

    @classmethod
    def annotate_match_4(cls, df, db, feature=feature, alteration_type=alteration_type, alteration=alteration):
        match_columns = [feature, alteration_type, alteration]
        df = (df
              .merge(db[~db[alteration_type].eq('') & ~db[alteration].eq('')]
                     .loc[:, match_columns + [cls.evidence_map_str, cls.merged]]
                     .sort_values(cls.evidence_map_str, ascending=False)
                     .drop_duplicates(match_columns, keep='first'),
                     on=match_columns,
                     how='left'
                     )
              )
        index_match = df[cls.merged].eq(1)
        df.loc[index_match, cls.match_4] = 1
        df.loc[index_match, cls.evidence] = df.loc[index_match, cls.evidence_map_str]
        df.drop([cls.merged, cls.evidence_map_str], axis=1, inplace=True)
        return df

    @classmethod
    def annotate_other_datasources(cls, df, dbs):
        df = CancerHotspots.annotate(df, dbs)
        df = CancerHotspots3D.annotate(df, dbs)
        df = CancerGeneCensus.annotate(df, dbs)
        df = Cosmic.annotate(df, dbs)
        df = GSEACancerPathways.annotate(df, dbs)
        df = GSEACancerModules.annotate(df, dbs)
        return df

    @classmethod
    def format_db(cls, feature_columns, column_map, db):
        columns = feature_columns + [cls.feature_display, cls.predictive_implication]
        column_map[cls.predictive_implication] = cls.evidence
        db = db.loc[:, columns].drop_duplicates()
        db.rename(columns=column_map, inplace=True)
        db[cls.evidence_map_str] = db[cls.evidence].replace(cls.evidence_map)
        db.sort_values([cls.evidence_map_str, cls.feature_display], ascending=[False, True], inplace=True)
        db[cls.merged] = 1
        return db

    @classmethod
    def map_evidence(cls, df):
        df[cls.evidence_map_str] = df[cls.evidence].copy(deep=True)
        df[cls.evidence] = df[cls.evidence].fillna('').replace(cls.inv_evidence_map)
        return df

    @classmethod
    def preallocate_almanac_matches(cls, df):
        df[cls.match_1] = 0
        df[cls.match_2] = 0
        df[cls.match_3] = 0
        df[cls.match_4] = 0
        df[cls.feature_match] = 0
        df[cls.evidence] = pd.NA
        return df
