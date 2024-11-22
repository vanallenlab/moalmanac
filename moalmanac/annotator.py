import pandas as pd
import numpy as np
import scipy.stats as stats
import operator as op
import copy

import datasources
import features
import logger

from config import COLNAMES


class Annotator:
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
        logger.Messages.dataframe_size(label="...datasource size", dataframe=ds)

        df[bin_name] = cls.match_ds(df, ds, bin_name, comparison_columns)
        message = f"...annotation complete, {df[bin_name].notnull().shape[0]} records annotated"
        logger.Messages.general(message=message)
        for value, group in df.groupby(bin_name):
            message = f"...{bin_name} == {value} for {group.shape[0]} records"
            logger.Messages.general(message=message)
        logger.Messages.general("")
        return df

    @classmethod
    def annotate_almanac(cls, df, dbs, ontology, config):
        df[cls.score_bin] = cls.preallocate_bin(cls.score_bin, df.index)
        df = Almanac.annotate(df, dbs, ontology, config)
        return df

    @classmethod
    def annotate_germline(cls, df, dbs, ontology, config):
        df[cls.score_bin] = cls.preallocate_bin(cls.score_bin, df.index)
        df = Almanac.annotate(df, dbs, ontology, config)
        df = CancerHotspots.annotate(df, dbs)
        df = CancerGeneCensus.annotate(df, dbs)
        df = ACMG.annotate(df, dbs)
        df = ClinVar.annotate(df, dbs)
        df = Hereditary.annotate(df, dbs)
        df = ExACExtended.annotate(df, dbs, config)
        df = MSI.annotate(df)
        return df

    @classmethod
    def annotate_simple(cls, df, dbs, ontology, config):
        df[cls.score_bin] = cls.preallocate_bin(cls.score_bin, df.index)
        df = Almanac.annotate(df, dbs, ontology, config)
        df = CancerHotspots.annotate(df, dbs)
        df = CancerHotspots3D.annotate(df, dbs)
        df = CancerGeneCensus.annotate(df, dbs)
        df = Cosmic.annotate(df, dbs)
        df = GSEACancerPathways.annotate(df, dbs)
        df = GSEACancerModules.annotate(df, dbs)
        df = ACMG.annotate(df, dbs)
        df = Hereditary.annotate(df, dbs)
        df = MSI.annotate(df)
        return df

    @classmethod
    def annotate_somatic(cls, df, dbs, ontology, config):
        df[cls.score_bin] = cls.preallocate_bin(cls.score_bin, df.index)
        df = Almanac.annotate(df, dbs, ontology, config)
        df = CancerHotspots.annotate(df, dbs)
        df = CancerHotspots3D.annotate(df, dbs)
        df = CancerGeneCensus.annotate(df, dbs)
        df = Cosmic.annotate(df, dbs)
        df = GSEACancerPathways.annotate(df, dbs)
        df = GSEACancerModules.annotate(df, dbs)
        df = ExAC.annotate(df, dbs, config)
        df = MSI.annotate(df)
        return df

    @classmethod
    def annotate_somatic_no_exac(cls, df, dbs, ontology, config):
        df[cls.score_bin] = cls.preallocate_bin(cls.score_bin, df.index)
        df = Almanac.annotate(df, dbs, ontology, config)
        df = CancerHotspots.annotate(df, dbs)
        df = CancerHotspots3D.annotate(df, dbs)
        df = CancerGeneCensus.annotate(df, dbs)
        df = Cosmic.annotate(df, dbs)
        df = GSEACancerPathways.annotate(df, dbs)
        df = GSEACancerModules.annotate(df, dbs)
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

    @staticmethod
    def dataframe_drop_na(dataframe, columns, how='any', axis=0):
        return dataframe.loc[:, columns].dropna(how=how, axis=axis)

    @classmethod
    def fill_na(cls, dataframe, column, fill_value, fill_data_type):
        if column not in dataframe.columns:
            dataframe = cls.preallocate_empty_columns(dataframe, [column])
        return dataframe.loc[dataframe.index, column].astype(fill_data_type).fillna(fill_value)

    @classmethod
    def get_idx_without_na(cls, dataframe, columns, datasource_label):
        subset = cls.dataframe_drop_na(dataframe, columns, how='any', axis=0)
        idx = subset.index
        removed_rows = dataframe.index.difference(idx).shape[0]
        message = (f"...{removed_rows} of {dataframe.shape[0]} rows have null values in "
                   f"one or more required columns: {', '.join(columns)}")
        logger.Messages.general(message=message)
        if removed_rows != 0:
            message = f"...these rows will not be annotated with {datasource_label}"
            logger.Messages.general(message=message)
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
        for column in columns:
            df[column] = None
        return df


class ACMG:
    gene = datasources.ACMG.gene

    bin_name = Annotator.acmg_bin
    comparison_columns = [gene]

    @classmethod
    def annotate(cls, df, dbs):
        return Annotator.annotate(df, dbs, datasources.ACMG, cls.bin_name, cls.comparison_columns)


class Almanac:
    feature_type = datasources.Almanac.feature_type
    feature = datasources.Almanac.feature
    alt_type = datasources.Almanac.alt_type
    alt = datasources.Almanac.alt

    disease = datasources.Almanac.disease
    oncotree_term = datasources.Almanac.ontology
    oncotree_code = datasources.Almanac.code
    context = datasources.Almanac.context
    therapy = datasources.Almanac.therapy
    therapy_strategy = datasources.Almanac.therapy_strategy
    therapy_type = datasources.Almanac.therapy_type
    sensitivity = datasources.Almanac.sensitivity
    resistance = datasources.Almanac.resistance
    prognosis = datasources.Almanac.prognosis
    implication = datasources.Almanac.implication
    implication_map = datasources.Almanac.implication_map
    description = datasources.Almanac.description

    source_type = datasources.Almanac.source_type
    citation = datasources.Almanac.citation
    url = datasources.Almanac.url
    doi = datasources.Almanac.doi
    pmid = datasources.Almanac.pmid
    nct = datasources.Almanac.nct
    publication_date = datasources.Almanac.publication_date
    last_updated = datasources.Almanac.last_updated

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
    sensitive_therapy_strategy = COLNAMES[datasources_section]['sensitive_therapy_strategy']
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
    resistance_therapy_strategy = COLNAMES[datasources_section]['resistance_therapy_strategy']
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
        therapy_strategy: sensitive_therapy_strategy,
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
        therapy_strategy: resistance_therapy_strategy,
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
    genes = datasources.Almanac.genes

    query_key = 'query'
    query_values = 'query_values'

    matches = 'matches'
    sensitivity_matches = datasources.Almanac.sensitivity_matches
    resistance_matches = datasources.Almanac.resistance_matches
    prognostic_matches = datasources.Almanac.prognostic_matches

    assertion_types_dict = {
        sensitivity: {
            query_key: sensitivity,
            query_values: [1],
            score_bin: sensitive_score_bin,
            columns: column_map_sensitive,
            feature: sensitive_feature,
            feature_type: sensitive_feature_type,
            alt_type: sensitive_alt_type,
            alt: sensitive_alt,
            matches: sensitivity_matches
        },
        resistance: {
            query_key: resistance,
            query_values: [1],
            score_bin: resistance_score_bin,
            columns: column_map_resistance,
            feature: resistance_feature,
            feature_type: resistance_feature_type,
            alt_type: resistance_alt_type,
            alt: resistance_alt,
            matches: resistance_matches
        },
        prognosis: {
            query_key: prognosis,
            query_values: [0, 1],
            score_bin: prognostic_score_bin,
            columns: column_map_prognostic,
            feature: prognostic_feature,
            feature_type: prognostic_feature_type,
            alt_type: prognostic_alt_type,
            alt: prognostic_alt,
            matches: prognostic_matches
        }
    }

    @classmethod
    def annotate(cls, df, dbs, ontology, config):
        logger.Messages.general(message="...with MOAlmanac's database")
        db = datasources.Almanac.import_ds(dbs)
        ds = db['content']
        list_genes = db['genes']

        annotation_function_dict = {
            config['feature_types']['aneuploidy']: cls.annotate_aneuploidy,
            config['feature_types']['burden']: cls.annotate_burden,
            config['feature_types']['cna']: cls.annotate_copy_number,
            config['feature_types']['fusion']: cls.annotate_fusion,
            config['feature_types']['germline']: cls.annotate_variants,
            config['feature_types']['microsatellite']: cls.annotate_microsatellite_stability,
            config['feature_types']['signature']: cls.annotate_signatures,
            config['feature_types']['mut']: cls.annotate_variants
        }

        for feature_type, group in df.groupby(cls.feature_type):
            logger.Messages.general(message=f"...annotating inputs of type {feature_type} with MOAlmanac's database")
            feature_type_records = cls.subset_records(ds, cls.feature_type, feature_type)
            table = pd.DataFrame(feature_type_records)
            logger.Messages.general(message=f"...records of {feature_type} in the database: {table.shape[0]}")

            # this is required for python 3.12 and pandas 2.2.2 to opt into future behavior for type downcasting
            with pd.option_context("future.no_silent_downcasting", True):
                table[cls.implication_map] = (
                    table[cls.implication]
                    .astype(str)
                    .replace(cls.predictive_implication_map)
                    .astype(float)
                )

            simple_biomarkers = [
                config['feature_types']['mut'],
                config['feature_types']['germline'],
                config['feature_types']['cna'],
                config['feature_types']['fusion']
            ]
            if feature_type in simple_biomarkers:
                idx = group[cls.feature].isin(list_genes)
                df.loc[group[~idx].index, cls.bin_name] = 0
                group = group[group[cls.feature].isin(list_genes)]

            logger.Messages.general(message=f"...records of {feature_type} provided: {group.shape[0]}")
            for index in group.index:
                annotation_function = annotation_function_dict[feature_type]
                new_series = annotation_function(sliced_series=df.loc[index, :], ontology=ontology, table=table)
                df.loc[index, new_series.index] = new_series
            message = f"...annotating input {feature_type}s with MOAlmanac's database complete"
            logger.Messages.general(message=message, add_line_break=True)

        return df

    @classmethod
    def annotate_aneuploidy(cls, sliced_series, ontology, table):
        series = sliced_series.copy(deep=True)
        feature = series.loc[cls.feature]

        query_same_ontology = (table[cls.oncotree_code] == ontology)
        query_diff_ontology = (table[cls.oncotree_code] != ontology)
        query_feature = (table[cls.event] == feature)

        match_bins = []

        for assertion_type in cls.assertion_types_dict.keys():
            message = f"...annotating {feature} {assertion_type}..."
            logger.Messages.general(message=message)

            query_key = cls.assertion_types_dict[assertion_type][cls.query_key]
            query_values = cls.assertion_types_dict[assertion_type][cls.query_values]
            query = table[query_key].isin(query_values)
            score_bin = cls.assertion_types_dict[assertion_type][cls.score_bin]
            columns = copy.deepcopy(cls.assertion_types_dict[assertion_type][cls.columns])
            columns[cls.feature_type] = cls.assertion_types_dict[assertion_type][cls.feature_type]
            columns[cls.event] = cls.assertion_types_dict[assertion_type][cls.feature]

            if (query_feature & query).any():
                results_same_ontology = (query_same_ontology & query_feature & query)
                results_diff_ontology = (query_diff_ontology & query_feature & query)
                feature_match_to_assertion_bin = 3

                message = "......was able to match on biomarker"
                logger.Messages.general(message=message)

            else:
                feature_match_to_assertion_bin = 0
                match_bins.append(feature_match_to_assertion_bin)

                message = "......was unable to match on biomarker, no specific database records matched"
                logger.Messages.general(message=message)

                continue

            matches = cls.sort_and_subset_matches(table, results_same_ontology, results_diff_ontology)
            series = cls.update_series_with_best_match(matches, columns, series)
            series.loc[score_bin] = feature_match_to_assertion_bin
            match_bins.append(feature_match_to_assertion_bin)

            column = cls.assertion_types_dict[assertion_type][cls.matches]
            series.loc[column] = matches

            message = f"......{len(matches)} matches for {assertion_type} within the database using this criteria"
            logger.Messages.general(message=message)

        series.loc[cls.bin_name] = max(match_bins)
        return series

    @classmethod
    def annotate_burden(cls, sliced_series, ontology, table):
        series = sliced_series.copy(deep=True)
        feature = series.loc[cls.feature].split(' ')[0]

        query_same_ontology = (table[cls.oncotree_code] == ontology)
        query_diff_ontology = (table[cls.oncotree_code] != ontology)
        query_feature = (table[cls.classification] == feature)

        match_bins = []

        for assertion_type in cls.assertion_types_dict.keys():
            message = f"...annotating {feature} {assertion_type}..."
            logger.Messages.general(message=message)

            query_key = cls.assertion_types_dict[assertion_type][cls.query_key]
            query_values = cls.assertion_types_dict[assertion_type][cls.query_values]
            query = table[query_key].isin(query_values)
            score_bin = cls.assertion_types_dict[assertion_type][cls.score_bin]
            columns = copy.deepcopy(cls.assertion_types_dict[assertion_type][cls.columns])
            columns[cls.feature_type] = cls.assertion_types_dict[assertion_type][cls.feature_type]
            columns[cls.classification] = cls.assertion_types_dict[assertion_type][cls.feature]

            if (query_feature & query).any():
                results_same_ontology = (query_same_ontology & query_feature & query)
                results_diff_ontology = (query_diff_ontology & query_feature & query)
                feature_match_to_assertion_bin = 3

                message = "......was able to match on biomarker"
                logger.Messages.general(message=message)

            else:
                feature_match_to_assertion_bin = 0
                match_bins.append(feature_match_to_assertion_bin)

                message = "......was unable to match on biomarker, no specific database records matched"
                logger.Messages.general(message=message)

                continue

            matches = cls.sort_and_subset_matches(table, results_same_ontology, results_diff_ontology)
            series = cls.update_series_with_best_match(matches, columns, series)
            series.loc[score_bin] = feature_match_to_assertion_bin
            match_bins.append(feature_match_to_assertion_bin)

            column = cls.assertion_types_dict[assertion_type][cls.matches]
            series.loc[column] = matches

            message = f"......{len(matches)} matches for {assertion_type} within the database using this criteria"
            logger.Messages.general(message=message)

        series.loc[cls.bin_name] = max(match_bins)
        return series

    @classmethod
    def annotate_copy_number(cls, sliced_series, ontology, table):
        series = sliced_series.copy(deep=True)
        feature = series.loc[cls.feature]
        alt_type = series.loc[cls.alt_type]

        query_same_ontology = (table[cls.oncotree_code] == ontology)
        query_diff_ontology = (table[cls.oncotree_code] != ontology)
        query_feature = (table[cls.gene] == feature)
        query_alt_type = (table[cls.direction] == alt_type)

        query_to_alt_type = query_feature & query_alt_type

        match_bins = []

        for assertion_type in cls.assertion_types_dict.keys():
            message = f"...annotating {feature} {alt_type} for {assertion_type}..."
            logger.Messages.general(message=message)

            query_key = cls.assertion_types_dict[assertion_type][cls.query_key]
            query_values = cls.assertion_types_dict[assertion_type][cls.query_values]
            query = table[query_key].isin(query_values)
            score_bin = cls.assertion_types_dict[assertion_type][cls.score_bin]
            columns = copy.deepcopy(cls.assertion_types_dict[assertion_type][cls.columns])
            columns[cls.feature_type] = cls.assertion_types_dict[assertion_type][cls.feature_type]
            columns[cls.gene] = cls.assertion_types_dict[assertion_type][cls.feature]
            columns[cls.direction] = cls.assertion_types_dict[assertion_type][cls.alt_type]

            if (query_to_alt_type & query).any():
                results_same_ontology = (query_same_ontology & query_to_alt_type & query)
                results_diff_ontology = (query_diff_ontology & query_to_alt_type & query)
                feature_match_to_assertion_bin = 4

                message = "......was able to match on gene, biomarker type, and direction"
                logger.Messages.general(message=message)

            elif (query_feature & query).any():
                results_same_ontology = (query_same_ontology & query_feature & query)
                results_diff_ontology = (query_diff_ontology & query_feature & query)
                feature_match_to_assertion_bin = 2

                message = "......was unable to match on gene, biomarker type, and direction"
                logger.Messages.general(message=message)

                message = "......was able to match on gene and biomarker type"
                logger.Messages.general(message=message)

            else:
                feature_match_to_assertion_bin = 1
                match_bins.append(feature_match_to_assertion_bin)

                message = "......was unable to match on gene, biomarker type, and direction"
                logger.Messages.general(message=message)

                message = "......was unable to match on gene and biomarker type"
                logger.Messages.general(message=message)

                message = "......was only able to match on gene, no specific database records matched"
                logger.Messages.general(message=message)

                continue

            matches = cls.sort_and_subset_matches(table, results_same_ontology, results_diff_ontology)
            series = cls.update_series_with_best_match(matches, columns, series)
            series.loc[score_bin] = feature_match_to_assertion_bin
            match_bins.append(feature_match_to_assertion_bin)

            column = cls.assertion_types_dict[assertion_type][cls.matches]
            series.loc[column] = matches

            message = f"......{len(matches)} matches for {assertion_type} within the database using this criteria"
            logger.Messages.general(message=message)

        series.loc[cls.bin_name] = max(match_bins)
        return series

    @classmethod
    def annotate_fusion(cls, sliced_series, ontology, table):
        series = sliced_series.fillna(pd.NA).copy(deep=True)
        feature = series.loc[cls.feature]
        alt_type = series.loc[cls.alt_type]
        alt = series.loc[cls.alt]

        if '--' in str(alt):
            split_alt = alt.split('--')
            gene1 = split_alt[0]
            gene2 = split_alt[1]
        else:
            gene1 = feature
            gene2 = ''

        query_same_ontology = (table[cls.oncotree_code] == ontology)
        query_diff_ontology = (table[cls.oncotree_code] != ontology)
        query_alt_type = (table[cls.rearrangement_type] == alt_type)

        match_bins = []

        for assertion_type in cls.assertion_types_dict.keys():
            message = f"...annotating {feature} {alt_type} {alt} for {assertion_type}..."
            logger.Messages.general(message=message)

            query_key = cls.assertion_types_dict[assertion_type][cls.query_key]
            query_values = cls.assertion_types_dict[assertion_type][cls.query_values]
            query = table[query_key].isin(query_values)
            score_bin = cls.assertion_types_dict[assertion_type][cls.score_bin]
            columns = copy.deepcopy(cls.assertion_types_dict[assertion_type][cls.columns])
            columns[cls.feature_type] = cls.assertion_types_dict[assertion_type][cls.feature_type]
            columns[cls.rearrangement_type] = cls.assertion_types_dict[assertion_type][cls.alt_type]
            alt_column = cls.assertion_types_dict[assertion_type][cls.alt]

            fusion_matches = []
            for first_gene, second_gene, feature_col in [(gene1, gene2, cls.gene1), (gene2, gene1, cls.gene2)]:
                columns[feature_col] = cls.assertion_types_dict[assertion_type][cls.feature]

                message = f"......processing as {first_gene}::{second_gene}..."
                logger.Messages.general(message=message)

                query_first_gene = (table[cls.gene1] == first_gene)
                query_second_gene = (table[cls.gene2] == second_gene)
                query_first_gene_reverse = (table[cls.gene1] == second_gene)
                query_second_gene_reverse = (table[cls.gene2] == first_gene)

                query_to_alt = query_first_gene & query_second_gene & query_alt_type
                query_to_alt_reverse = query_first_gene_reverse & query_second_gene_reverse & query_alt_type

                query_to_alt_type = query_first_gene & query_alt_type
                query_to_alt_type_reverse = query_first_gene_reverse & query_alt_type

                if (query_to_alt & query).any():
                    results_same_ontology = (query_same_ontology & query_to_alt & query)
                    results_diff_ontology = (query_diff_ontology & query_to_alt & query)
                    match_bin = 4

                    message = "......was able to match on first gene, second gene, biomarker type, and rearrangement type"
                    logger.Messages.general(message=message)

                elif (query_to_alt_reverse & query).any():
                    results_same_ontology = (query_same_ontology & query_to_alt_reverse & query)
                    results_diff_ontology = (query_diff_ontology & query_to_alt_reverse & query)
                    match_bin = 4

                    message = "......was able to match on second gene, first gene, biomarker type, and rearrangement type"
                    logger.Messages.general(message=message)

                elif (query_to_alt_type & query).any():
                    results_same_ontology = (query_same_ontology & query_to_alt_type & query)
                    results_diff_ontology = (query_diff_ontology & query_to_alt_type & query)
                    match_bin = 3

                    message = "......was unable to match on first gene, second gene, biomarker type, and rearrangement type"
                    logger.Messages.general(message=message)

                    message = "......was able to match on first gene, biomarker type, and rearrangement type"
                    logger.Messages.general(message=message)

                elif (query_to_alt_type_reverse & query).any():
                    results_same_ontology = (query_same_ontology & query_to_alt_type_reverse & query)
                    results_diff_ontology = (query_diff_ontology & query_to_alt_type_reverse & query)
                    match_bin = 3

                    message = "......was unable to match on second gene, first gene, biomarker type, and rearrangement type"
                    logger.Messages.general(message=message)

                    message = "......was able to match on second gene, biomarker type, and rearrangement type"
                    logger.Messages.general(message=message)

                elif (query_first_gene & query).any():
                    results_same_ontology = (query_same_ontology & query_first_gene & query)
                    results_diff_ontology = (query_diff_ontology & query_first_gene & query)
                    match_bin = 2

                    message = "......was unable to match on first gene, second gene, biomarker type, and rearrangement type"
                    logger.Messages.general(message=message)

                    message = "......was unable to match on first gene, biomarker type, and rearrangement type"
                    logger.Messages.general(message=message)

                    message = "......was able to match on first gene and biomarker type"
                    logger.Messages.general(message=message)

                elif (query_second_gene & query).any():
                    results_same_ontology = (query_same_ontology & query_second_gene & query)
                    results_diff_ontology = (query_diff_ontology & query_second_gene & query)
                    match_bin = 2

                    message = "......was unable to match on second gene, first gene, biomarker type, and rearrangement type"
                    logger.Messages.general(message=message)

                    message = "......was unable to match on second gene, biomarker type, and rearrangement type"
                    logger.Messages.general(message=message)

                    message = "......was able to match on second gene and biomarker type"
                    logger.Messages.general(message=message)

                else:
                    results_same_ontology = []
                    results_diff_ontology = []
                    match_bin = 1

                    message = "......was only able to match on gene, no specific database records matched"
                    logger.Messages.general(message=message)


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
            matches = cls.sort_and_subset_matches(table, results_same_ontology, results_diff_ontology)

            feature_col = better_match['feature_col']
            columns[feature_col] = cls.assertion_types_dict[assertion_type][cls.feature]
            series = cls.update_series_with_best_match(matches, columns, series)
            series.loc[alt_column] = better_match['alt']
            series.loc[score_bin] = feature_match_to_assertion_bin
            match_bins.append(feature_match_to_assertion_bin)

            column = cls.assertion_types_dict[assertion_type][cls.matches]
            series.loc[column] = matches

            message = f"......{len(matches)} matches for {assertion_type} within the database using this criteria"
            logger.Messages.general(message=message)

        series.loc[cls.bin_name] = max(match_bins)
        return series

    @classmethod
    def annotate_microsatellite_stability(cls, sliced_series, ontology, table):
        series = sliced_series.copy(deep=True)
        feature = series.loc[cls.feature].split(' ')[-1]

        query_same_ontology = (table[cls.oncotree_code] == ontology)
        query_diff_ontology = (table[cls.oncotree_code] != ontology)
        query_feature = (table[cls.status] == feature)

        match_bins = []

        for assertion_type in cls.assertion_types_dict.keys():
            message = f"...annotating {feature} {assertion_type}..."
            logger.Messages.general(message=message)

            query_key = cls.assertion_types_dict[assertion_type][cls.query_key]
            query_values = cls.assertion_types_dict[assertion_type][cls.query_values]
            query = table[query_key].isin(query_values)
            score_bin = cls.assertion_types_dict[assertion_type][cls.score_bin]
            columns = copy.deepcopy(cls.assertion_types_dict[assertion_type][cls.columns])
            columns[cls.feature_type] = cls.assertion_types_dict[assertion_type][cls.feature_type]
            columns[cls.status] = cls.assertion_types_dict[assertion_type][cls.feature]

            if (query_feature & query).any():
                results_same_ontology = (query_same_ontology & query_feature & query)
                results_diff_ontology = (query_diff_ontology & query_feature & query)
                feature_match_to_assertion_bin = 3

                message = "......was able to match on biomarker"
                logger.Messages.general(message=message)

            else:
                feature_match_to_assertion_bin = 0
                match_bins.append(feature_match_to_assertion_bin)

                message = "......was unable to match on biomarker, no specific database records matched"
                logger.Messages.general(message=message)

                continue

            matches = cls.sort_and_subset_matches(table, results_same_ontology, results_diff_ontology)
            series = cls.update_series_with_best_match(matches, columns, series)
            series.loc[score_bin] = feature_match_to_assertion_bin
            match_bins.append(feature_match_to_assertion_bin)

            column = cls.assertion_types_dict[assertion_type][cls.matches]
            series.loc[column] = matches

            message = f"......{len(matches)} matches for {assertion_type} within the database using this criteria"
            logger.Messages.general(message=message)

        series.loc[cls.bin_name] = max(match_bins)
        return series

    @classmethod
    def annotate_signatures(cls, sliced_series, ontology, table):
        series = sliced_series.copy(deep=True)
        feature = series.loc[cls.feature]

        query_same_ontology = (table[cls.oncotree_code] == ontology)
        query_diff_ontology = (table[cls.oncotree_code] != ontology)
        query_feature = (table[cls.cosmic_signature_number].astype(str) == str(feature))

        match_bins = []

        for assertion_type in cls.assertion_types_dict.keys():
            message = f"...annotating {feature} {assertion_type}..."
            logger.Messages.general(message=message)

            query_key = cls.assertion_types_dict[assertion_type][cls.query_key]
            query_values = cls.assertion_types_dict[assertion_type][cls.query_values]
            query = table[query_key].isin(query_values)
            score_bin = cls.assertion_types_dict[assertion_type][cls.score_bin]
            columns = copy.deepcopy(cls.assertion_types_dict[assertion_type][cls.columns])
            columns[cls.feature_type] = cls.assertion_types_dict[assertion_type][cls.feature_type]
            columns[cls.cosmic_signature_number] = cls.assertion_types_dict[assertion_type][cls.feature]

            if (query_feature & query).any():
                results_same_ontology = (query_same_ontology & query_feature & query)
                results_diff_ontology = (query_diff_ontology & query_feature & query)
                feature_match_to_assertion_bin = 3

                message = "......was able to match on biomarker"
                logger.Messages.general(message=message)

            else:
                feature_match_to_assertion_bin = 1
                match_bins.append(feature_match_to_assertion_bin)

                message = "......was unable to match on biomarker, no specific database records matched"
                logger.Messages.general(message=message)

                continue

            matches = cls.sort_and_subset_matches(table, results_same_ontology, results_diff_ontology)
            series = cls.update_series_with_best_match(matches, columns, series)
            series.loc[score_bin] = feature_match_to_assertion_bin
            match_bins.append(feature_match_to_assertion_bin)

            column = cls.assertion_types_dict[assertion_type][cls.matches]
            series.loc[column] = matches

            message = f"......{len(matches)} matches for {assertion_type} within the database using this criteria"
            logger.Messages.general(message=message)

        series.loc[cls.bin_name] = max(match_bins)
        return series

    @classmethod
    def annotate_variants(cls, sliced_series, ontology, table):
        series = sliced_series.copy(deep=True)
        feature = series.loc[cls.feature]
        alt_type = series.loc[cls.alt_type]
        alt = series.loc[cls.alt]

        query_same_ontology = (table[cls.oncotree_code] == ontology)
        query_diff_ontology = (table[cls.oncotree_code] != ontology)
        query_feature = (table[cls.gene] == feature)
        query_alt_type = (table[cls.variant_annotation] == alt_type)
        query_alt = (table[cls.protein_change] == alt)

        query_to_alt = query_feature & query_alt_type & query_alt
        query_to_alt_type = query_feature & query_alt_type

        match_bins = []

        for assertion_type in cls.assertion_types_dict.keys():
            message = f"...annotating {feature} {alt_type} {alt} for {assertion_type}..."
            logger.Messages.general(message=message)

            query_key = cls.assertion_types_dict[assertion_type][cls.query_key]
            query_values = cls.assertion_types_dict[assertion_type][cls.query_values]
            query = table[query_key].isin(query_values)
            score_bin = cls.assertion_types_dict[assertion_type][cls.score_bin]
            columns = copy.deepcopy(cls.assertion_types_dict[assertion_type][cls.columns])
            columns[cls.feature_type] = cls.assertion_types_dict[assertion_type][cls.feature_type]
            columns[cls.gene] = cls.assertion_types_dict[assertion_type][cls.feature]
            columns[cls.variant_annotation] = cls.assertion_types_dict[assertion_type][cls.alt_type]
            columns[cls.protein_change] = cls.assertion_types_dict[assertion_type][cls.alt]

            if (query_to_alt & query).any():
                results_same_ontology = (query_same_ontology & query_to_alt & query)
                results_diff_ontology = (query_diff_ontology & query_to_alt & query)
                feature_match_to_assertion_bin = 4

                message = "......was able to match on gene, biomarker type, variant classification, and protein change"
                logger.Messages.general(message=message)

            elif (query_to_alt_type & query).any():
                results_same_ontology = (query_same_ontology & query_to_alt_type & query)
                results_diff_ontology = (query_diff_ontology & query_to_alt_type & query)
                feature_match_to_assertion_bin = 3

                message = "......was unable to match on gene, biomarker type, variant classification, and protein change"
                logger.Messages.general(message=message)

                message = "......was able to match on gene, biomarker type, and variant classification"
                logger.Messages.general(message=message)

            elif (query_feature & query).any():
                results_same_ontology = (query_same_ontology & query_feature & query)
                results_diff_ontology = (query_diff_ontology & query_feature & query)
                feature_match_to_assertion_bin = 2

                message = "......was unable to match on gene, biomarker type, variant classification, and protein change"
                logger.Messages.general(message=message)

                message = "......was unable to match on gene, biomarker type, and variant classification"
                logger.Messages.general(message=message)

                message = "......was able to match on gene and biomarker type"
                logger.Messages.general(message=message)

            else:
                feature_match_to_assertion_bin = 1
                match_bins.append(feature_match_to_assertion_bin)

                message = "......was unable to match on gene, biomarker type, variant classification, and protein change"
                logger.Messages.general(message=message)

                message = "......was unable to match on gene, biomarker type, and variant classification"
                logger.Messages.general(message=message)

                message = "......was unable to match on gene and biomarker type"
                logger.Messages.general(message=message)

                message = "......was only able to match on gene, no specific database records matched"
                logger.Messages.general(message=message)

                continue

            matches = cls.sort_and_subset_matches(table, results_same_ontology, results_diff_ontology)
            series = cls.update_series_with_best_match(matches, columns, series)
            series.loc[score_bin] = feature_match_to_assertion_bin
            match_bins.append(feature_match_to_assertion_bin)

            column = cls.assertion_types_dict[assertion_type][cls.matches]
            series.loc[column] = matches

            message = f"......{len(matches)} matches for {assertion_type} within the database using this criteria"
            logger.Messages.general(message=message)

        series.loc[cls.bin_name] = max(match_bins)
        return series

    @staticmethod
    def convert_dataframe_to_records(df):
        return df.to_dict(orient='records')

    @staticmethod
    def extract_values_from_list_of_dicts(key, list_of_dicts):
        return [x[key] for x in list_of_dicts]

    @classmethod
    def select_better_fusion_match(cls, list_of_matches):
        max_bin = max(match['bin'] for match in list_of_matches)
        better_match = [match for match in list_of_matches if (match['bin'] == max_bin)]
        return better_match[0]

    @classmethod
    def sort_and_subset_matches(cls, table, results_same_ontology, results_diff_ontology):
        sort_columns = [cls.implication_map, cls.publication_date, cls.last_updated]
        sort_ascending = [False, False, False]
        if results_same_ontology.any():
            results_same_ontology_sorted = cls.sort_dataframe(table[results_same_ontology], sort_columns, sort_ascending)
            results_same_ontology_records = cls.convert_dataframe_to_records(results_same_ontology_sorted)
        else:
            results_same_ontology_records = []

        if results_diff_ontology.any():
            results_diff_ontology_sorted = cls.sort_dataframe(table[results_diff_ontology], sort_columns, sort_ascending)
            results_diff_ontology_records = cls.convert_dataframe_to_records(results_diff_ontology_sorted)
        else:
            results_diff_ontology_records = []

        best_evidence = cls.return_best_evidence_level(results_same_ontology_records, results_diff_ontology_records)

        matches = []
        for results in [results_same_ontology_records, results_diff_ontology_records]:
            results_returned = cls.subset_list_of_dictionaries(results, cls.implication_map, op.ge, best_evidence)
            matches.extend(results_returned)
        return matches

    @staticmethod
    def sort_dataframe(df, sort_columns, ascending_boolean):
        return df.sort_values(sort_columns, ascending=ascending_boolean)

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

    @staticmethod
    def subset_records(records, key, value):
        return [record for record in records if record[key] == value]

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
    def update_series_with_best_match(cls, matches, columns, series):
        best_match = matches[0]
        for column, assertion_column in columns.items():
            series.loc[assertion_column] = best_match[column]
        return series


class CancerHotspots:
    gene = datasources.CancerHotspots.gene
    alteration = datasources.CancerHotspots.alt

    bin_name = Annotator.cancerhotspots_bin
    comparison_columns = [gene, alteration]

    @classmethod
    def annotate(cls, df, dbs):
        logger.Messages.general(message="...with Cancer Hotspots")
        return Annotator.annotate(df, dbs, datasources.CancerHotspots, cls.bin_name, cls.comparison_columns)


class CancerHotspots3D:
    gene = datasources.CancerHotspots3D.gene
    alteration = datasources.CancerHotspots3D.alt

    bin_name = Annotator.cancerhotspots3d_bin
    comparison_columns = [gene, alteration]

    @classmethod
    def annotate(cls, df, dbs):
        logger.Messages.general(message="...with Cancer Hotspots 3D")
        return Annotator.annotate(df, dbs, datasources.CancerHotspots3D, cls.bin_name, cls.comparison_columns)


class CancerGeneCensus:
    gene = datasources.CancerGeneCensus.gene

    bin_name = Annotator.cgc_bin
    comparison_columns = [gene]

    @classmethod
    def annotate(cls, df, dbs):
        logger.Messages.general(message="...with Cancer Gene Census")
        return Annotator.annotate(df, dbs, datasources.CancerGeneCensus, cls.bin_name, cls.comparison_columns)


class ClinVar:
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

        return df.reset_index().merge(ds, how='left').set_index('index')

    @classmethod
    def annotate(cls, df, dbs):
        logger.Messages.general(message="...with ClinVar")
        columns = [cls.chr, cls.start, cls.end]
        idx = Annotator.get_idx_without_na(dataframe=df, columns=columns, datasource_label='ClinVar')
        idx_result = df.index.copy(deep=True)

        if df.index.difference(idx).empty:
            remaining_variants = features.Features.create_empty_dataframe()
        else:
            remaining_variants = df.loc[df.index.difference(idx), :].copy(deep=True)
            df = df.loc[idx, :].copy(deep=True)

        subset = df.drop(df.columns[df.columns.str.contains('clinvar')], axis=1)
        subset = subset.loc[idx, :]
        ds = datasources.ClinVar.import_ds(dbs)
        logger.Messages.dataframe_size(label="...datasource size", dataframe=ds)

        subset = cls.append_clinvar(subset, ds)
        for value, group in subset.groupby(cls.bin_name):
            message = f"...{cls.bin_name} == {value} for {group.shape[0]} records"
            logger.Messages.general(message=message)
        logger.Messages.general("")

        if remaining_variants.empty:
            result = subset.loc[idx_result]
        else:
            subset = features.Features.preallocate_missing_columns(subset)
            remaining_variants = features.Features.preallocate_missing_columns(remaining_variants)
            list_dataframes = [subset, remaining_variants]
            result = features.Features.concat_list_of_dataframes(list_of_dataframes=list_dataframes, ignore_index=False)
            result = result.loc[idx_result]
        return features.Features.preallocate_missing_columns(result)


class Cosmic:
    gene = datasources.Cosmic.gene
    alteration = datasources.Cosmic.alt

    bin_name = Annotator.cosmic_bin
    comparison_columns = [gene, alteration]

    @classmethod
    def annotate(cls, df, dbs):
        logger.Messages.general(message="...with COSMIC")
        return Annotator.annotate(df, dbs, datasources.Cosmic, cls.bin_name, cls.comparison_columns)


class ExAC:
    chr = datasources.ExAC.chr
    start = datasources.ExAC.start
    ref = datasources.ExAC.ref
    alt = datasources.ExAC.alt
    af = datasources.ExAC.af

    bin_name = Annotator.exac_common_bin

    str_columns = [chr, ref, alt]
    int_columns = [start]

    feature_type = features.Features.feature_type

    @classmethod
    def append_exac_af(cls, df, ds, ds_columns, variant_biomarker_types):
        variants, not_variants = cls.subset_for_variants(df, variant_biomarker_types)
        ds = ds.loc[:, ds_columns]

        for column, data_type in [(cls.str_columns, str), (cls.int_columns, float), (cls.int_columns, int)]:
            variants[column] = variants[column].astype(data_type)
            ds[column] = ds[column].astype(data_type)

        merged = variants.reset_index().merge(ds, how='left').set_index('index')
        merged.loc[merged.index, cls.af] = Annotator.fill_na(
            dataframe=merged,
            column=cls.af,
            fill_value=0.0,
            fill_data_type=float
        )
        not_variants.loc[not_variants.index, cls.af] = Annotator.fill_na(
            dataframe=not_variants,
            column=cls.af,
            fill_value=0.0,
            fill_data_type=float
        )
        result = pd.concat([merged, not_variants])
        result[cls.af] = result[cls.af].astype(float).round(6)
        return result

    @classmethod
    def annotate(cls, df, dbs, config):
        logger.Messages.general(message="...with ExAC")
        columns = [cls.chr, cls.start, cls.ref, cls.alt]
        idx = Annotator.get_idx_without_na(dataframe=df, columns=columns, datasource_label='ExAC')
        idx_result = df.index.copy(deep=True)

        if df.index.difference(idx).empty:
            remaining_variants = features.Features.create_empty_dataframe()
        else:
            remaining_variants = df.loc[df.index.difference(idx), :].copy(deep=True)
            df = df.loc[idx, :].copy(deep=True)

        subset = cls.drop_existing_columns(df)
        subset = subset.loc[idx, :]
        ds = datasources.ExAC.import_ds(dbs)
        logger.Messages.dataframe_size(label="...datasource size", dataframe=ds)

        subset = cls.append_exac_af(
            df=subset,
            ds=ds,
            ds_columns=[cls.chr, cls.start, cls.ref, cls.alt, cls.af],
            variant_biomarker_types=[config['feature_types']['mut'], config['feature_types']['germline']]
        )
        common_allele_frequency_threshold = config['exac']['exac_common_af_threshold']
        subset[cls.bin_name] = cls.annotate_common_af(
            series_exac_af=subset[cls.af],
            threshold=common_allele_frequency_threshold
        )

        for value, group in subset.groupby(cls.bin_name):
            message = f"...{cls.bin_name} == {value} for {group.shape[0]} records"
            logger.Messages.general(message=message)
        logger.Messages.general("")

        if remaining_variants.empty:
            result = subset.loc[idx_result]
        else:
            subset = features.Features.preallocate_missing_columns(subset)
            remaining_variants = features.Features.preallocate_missing_columns(remaining_variants)
            list_dataframes = [subset, remaining_variants]
            result = features.Features.concat_list_of_dataframes(list_of_dataframes=list_dataframes, ignore_index=False)
            result = result.loc[idx_result]
        return features.Features.preallocate_missing_columns(result)

    @classmethod
    def annotate_common_af(cls, series_exac_af, threshold):
        if not series_exac_af.empty:
            series = pd.Series(float(0.0), index=series_exac_af.index.tolist())
            condition = (series_exac_af
                         .astype(str).str.split(',', expand=True)
                         .fillna(value=np.nan)
                         .astype(float).mean(axis=1)
                         .fillna(0.0)
                         )
            idx = condition.astype(float) >= float(threshold)
            series[idx] = float(1.0)
            return series
        else:
            return pd.Series()

    @classmethod
    def drop_existing_columns(cls, dataframe):
        return dataframe.drop(dataframe.columns[dataframe.columns.str.contains('exac')], axis=1)

    @classmethod
    def format_columns(cls, dataframe, column, data_type):
        return dataframe.loc[dataframe.index, column].astype(data_type)

    @classmethod
    def subset_for_variants(cls, dataframe, variant_biomarker_types):
        idx = dataframe[cls.feature_type].isin(variant_biomarker_types)
        return dataframe[idx].copy(), dataframe[~idx].copy()


class ExACExtended:
    chr = datasources.ExACExtended.chr
    start = datasources.ExACExtended.start
    ref = datasources.ExACExtended.ref
    alt = datasources.ExACExtended.alt
    af = datasources.ExACExtended.af
    ac = datasources.ExACExtended.ac
    an = datasources.ExACExtended.an
    ac_afr = datasources.ExACExtended.ac_afr
    ac_amr = datasources.ExACExtended.ac_amr
    ac_eas = datasources.ExACExtended.ac_eas
    ac_fin = datasources.ExACExtended.ac_fin
    ac_nfe = datasources.ExACExtended.ac_nfe
    ac_sas = datasources.ExACExtended.ac_sas
    ac_oth = datasources.ExACExtended.ac_oth
    an_afr = datasources.ExACExtended.an_afr
    an_amr = datasources.ExACExtended.an_amr
    an_eas = datasources.ExACExtended.an_eas
    an_fin = datasources.ExACExtended.an_fin
    an_nfe = datasources.ExACExtended.an_nfe
    an_sas = datasources.ExACExtended.an_sas
    an_oth = datasources.ExACExtended.an_oth

    ds_columns = [chr, start, ref, alt,
                  af, ac, an,
                  ac_afr, ac_amr, ac_eas, ac_fin, ac_nfe, ac_sas, ac_oth,
                  an_afr, an_amr, an_eas, an_fin, an_nfe, an_sas, an_oth]

    @classmethod
    def annotate(cls, df, dbs, config):
        logger.Messages.general(message="...with ExAC, Extended")
        columns = [cls.chr, cls.start, cls.ref, cls.alt]
        idx = Annotator.get_idx_without_na(dataframe=df, columns=columns, datasource_label='ExAC')
        idx_result = df.index.copy(deep=True)

        if df.index.difference(idx).empty:
            remaining_variants = features.Features.create_empty_dataframe()
        else:
            remaining_variants = df.loc[df.index.difference(idx), :].copy(deep=True)
            df = df.loc[idx, :].copy(deep=True)

        subset = ExAC.drop_existing_columns(df)
        subset = subset.loc[idx, :]
        ds = datasources.ExACExtended.import_ds(dbs)
        logger.Messages.dataframe_size(label="...datasource size", dataframe=ds)

        subset = ExAC.append_exac_af(
            df=subset,
            ds=ds,
            ds_columns=cls.ds_columns,
            variant_biomarker_types=[config['feature_types']['mut'], config['feature_types']['germline']]
        )

        common_allele_frequency_threshold = config['exac']['exac_common_af_threshold']
        subset[ExAC.bin_name] = ExAC.annotate_common_af(
            series_exac_af=subset[ExAC.af],
            threshold=common_allele_frequency_threshold
        )

        for value, group in subset.groupby(ExAC.bin_name):
            message = f"...{ExAC.bin_name} == {value} for {group.shape[0]} records"
            logger.Messages.general(message=message)
        logger.Messages.general("")

        if remaining_variants.empty:
            result = subset.loc[idx_result]
        else:
            subset = features.Features.preallocate_missing_columns(subset)
            remaining_variants = features.Features.preallocate_missing_columns(remaining_variants)
            list_dataframes = [subset, remaining_variants]
            result = features.Features.concat_list_of_dataframes(list_of_dataframes=list_dataframes, ignore_index=False)
            result = result.loc[idx_result]
            #result = pd.concat([subset, remaining_variants]).loc[idx_result, :]
        return features.Features.preallocate_missing_columns(result)


class GSEACancerModules:
    gene = datasources.GSEACancerModules.gene

    bin_name = Annotator.gsea_module_bin
    comparison_columns = [gene]

    @classmethod
    def annotate(cls, df, dbs):
        logger.Messages.general(message="...with GSEA Cancer Modules")
        return Annotator.annotate(df, dbs, datasources.GSEACancerModules, cls.bin_name, cls.comparison_columns)


class GSEACancerPathways:
    gene = datasources.GSEACancerPathways.gene

    bin_name = Annotator.gsea_pathway_bin
    comparison_columns = [gene]

    @classmethod
    def annotate(cls, df, dbs):
        logger.Messages.general(message="...with GSEA Cancer Pathways")
        return Annotator.annotate(df, dbs, datasources.GSEACancerPathways, cls.bin_name, cls.comparison_columns)


class Hereditary:
    gene = datasources.Hereditary.gene

    bin_name = Annotator.hereditary_bin
    comparison_columns = [gene]

    @classmethod
    def annotate(cls, df, dbs):
        logger.Messages.general(message="...with gene list containing genes related to hereditary cancers")
        return Annotator.annotate(df, dbs, datasources.Hereditary, cls.bin_name, cls.comparison_columns)


class MSI:
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
        logger.Messages.general(message="...with gene list with containing genes related to microsatellite instability")
        logger.Messages.general(message=f"...genes: {','.join(cls.msi_genes)}")
        df[cls.bin_name] = Annotator.match_ds(df, cls.create_msi_df(), cls.bin_name, cls.comparison_columns)
        for value, group in df.groupby(cls.bin_name):
            message = f"...{cls.bin_name} == {value} for {group.shape[0]} records"
            logger.Messages.general(message=message)

        return df


class OverlapValidation:
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

    merge_cols = [gene, alt_type, alt]
    fill_cols = [tumor_f, validation_tumor_f, validation_coverage]

    @classmethod
    def append_validation(cls, primary, validation, biomarker_type):
        df = cls.drop_validation_columns(primary)

        logger.Messages.general("...merging dataframes")
        df = cls.merge_data_frames(df, validation, cls.merge_cols)
        idx = cls.get_mutation_index(df, biomarker_type)
        for column in cls.fill_cols:
            df.loc[idx, column] = Annotator.fill_na(
                dataframe=df.loc[idx, :],
                column=column,
                fill_value=0.0,
                fill_data_type=float
            )
        count_match = df.loc[df[cls.validation_coverage].gt(0.0), :].shape[0]
        count_total = df.loc[idx, :].shape[0]
        message = f"...{count_match} of {count_total} variants were observed in validation sequencing"
        logger.Messages.general(message=message)

        logger.Messages.general(message="...calculating power to detect variants in validation sequencing")
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
    def get_mutation_index(cls, df, biomarker_type):
        return df[df[cls.feature_type].eq(biomarker_type)].index

    @classmethod
    def merge_data_frames(cls, df1, df2, columns):
        return df1.merge(df2, on=columns, how='left')

    @staticmethod
    def round_series(series, n):
        series = pd.to_numeric(series, errors='coerce')
        return series.round(n)


class OverlapSomaticGermline:
    section = 'overlap_somatic_germline'
    gene = COLNAMES[section]['gene']
    alt_type = COLNAMES[section]['alt_type']
    number_germline_mutations = COLNAMES[section]['number_germline_mutations']

    @classmethod
    def count_germline_hits(cls, df):
        variants = df[df[cls.alt_type].isin(features.MAF.coding_classifications_mapped)]
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

        count_match = merged_df.loc[merged_df[cls.number_germline_mutations].notnull(), :].shape[0]
        count_total = merged_df.shape[0]
        message = f"...{count_match} of {count_total} somatic variants had germline variants in the same gene"
        logger.Messages.general(message=message)

        return merged_df

    @classmethod
    def value_counts_to_df(cls, value_counts):
        return pd.DataFrame({cls.gene: value_counts.index, cls.number_germline_mutations: value_counts.values})


class PreclinicalEfficacy:
    section = 'preclinical'
    feature_display = COLNAMES[section]['feature_display']
    pvalue = COLNAMES[section]['pvalue']
    efficacy = COLNAMES[section]['efficacy_obs']
    lookup = COLNAMES[section]['efficacy_lookup']

    @classmethod
    def annotate(cls, actionable, efficacy, dictionary, append_lookup=True):
        series_all_features = actionable[cls.feature_display]
        series_features = series_all_features[series_all_features.isin(efficacy[cls.feature_display])]
        for index in series_features.index:
            feature = series_features.loc[index]
            dataframe = efficacy[efficacy[cls.feature_display].eq(feature)]
            efficacy_observed = cls.search_for_significance(dataframe[cls.pvalue])
            actionable.loc[index, cls.efficacy] = efficacy_observed
        actionable[cls.efficacy] = actionable[cls.efficacy].fillna(pd.NA)
        idx = actionable.index
        if append_lookup:
            actionable.loc[idx, cls.lookup] = cls.create_lookup(idx, series_features.index, dictionary)
        else:
            actionable.loc[idx, cls.lookup] = ''
        return actionable

    @classmethod
    def create_lookup(cls, all_index_values, relevant_index_values, dictionary):
        series = pd.Series('', index=all_index_values, name=cls.lookup)
        for index in relevant_index_values:
            series.loc[index] = [dictionary[index]]
        return series

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

    evidence_map = {
        'FDA-Approved': 5, 'Guideline': 4, 'Clinical trial': 3,
        'Clinical evidence': 2, 'Preclinical': 1, 'Inferential': 0}
    inv_evidence_map = {float(v): k for k, v in evidence_map.items()}

    variants = 'variants'
    cnas = 'copy_number_alterations'
    fusions = 'fusions'
    fusions_gene1 = 'fusions_gene1'
    fusions_gene2 = 'fusions_gene2'

    @classmethod
    def annotate(cls, input_dict, dbs, config):
        copy_number = config['feature_types']['cna']
        fusion = config['feature_types']['fusion']
        somatic_variant = config['feature_types']['mut']

        input_variants = input_dict[somatic_variant]
        input_copy_number_alterations = input_dict[copy_number]
        input_fusions = input_dict[fusion]

        variants = cls.annotate_somatic_variants(input_variants, dbs, somatic_variant)
        copy_number_alterations = cls.annotate_copy_numbers(input_copy_number_alterations, dbs, copy_number)
        fusions, fusions_gene1, fusions_gene2 = cls.annotate_fusions(input_fusions, dbs, fusion)
        return {
            cls.variants: variants,
            cls.cnas: copy_number_alterations,
            cls.fusions: fusions,
            cls.fusions_gene1: fusions_gene1,
            cls.fusions_gene2: fusions_gene2
        }

    @classmethod
    def annotate_copy_numbers(cls, df, dbs, biomarker_type_string):
        almanac = datasources.Almanac.import_ds(dbs)
        almanac_genes = datasources.Almanac.import_genes(dbs)

        df = df[df[cls.feature_type].eq(biomarker_type_string)]
        db = Almanac.subset_records(almanac['content'], cls.feature_type, biomarker_type_string)
        db = pd.DataFrame(db)

        column_map = {cls.gene: cls.feature, cls.direction: cls.alteration_type}
        db = cls.format_db(list(column_map.keys()), column_map, db)
        df = cls.preallocate_almanac_matches(df)
        df = cls.annotate_match_1(df, almanac_genes)
        df = cls.annotate_match_2(df, db)
        df = cls.annotate_match_3(df, db)
        df = cls.annotate_other_datasources(df, dbs)
        return df

    @classmethod
    def annotate_fusions(cls, df, dbs, biomarker_type_string):
        almanac = datasources.Almanac.import_ds(dbs)
        almanac_genes = datasources.Almanac.import_genes(dbs)

        df = df[df[cls.feature_type].eq(biomarker_type_string)]
        db = Almanac.subset_records(almanac['content'], cls.feature_type, biomarker_type_string)
        db = pd.DataFrame(db)

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

        # this is required for python 3.12 and pandas 2.2.2 to opt into future behavior for type downcasting
        with pd.option_context("future.no_silent_downcasting", True):
            values = (
                values
                .fillna(-1.0)
                .idxmax(axis=1)
            )

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

        # this is required for python 3.12 and pandas 2.2.2 to opt into future behavior for type downcasting
        with pd.option_context("future.no_silent_downcasting", True):
            values = (
                values
                .fillna(-1.0)
                .idxmax(axis=1)
            )

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

        # this is required for python 3.12 and pandas 2.2.2 to opt into future behavior for type downcasting
        with pd.option_context("future.no_silent_downcasting", True):
            values = (
                values
                .fillna(-1.0)
                .idxmax(axis=1)
            )

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
    def annotate_somatic_variants(cls, df, dbs, biomarker_type_string):
        almanac = datasources.Almanac.import_ds(dbs)
        almanac_genes = datasources.Almanac.import_genes(dbs)

        df = df[df[cls.feature_type].eq(biomarker_type_string)]
        db = Almanac.subset_records(almanac['content'], cls.feature_type, biomarker_type_string)
        db = pd.DataFrame(db)

        replacement_dictionary = {'Oncogenic Mutations': '', 'Activating mutation': ''}
        db[cls.variant_annotation] = db[cls.variant_annotation].replace(replacement_dictionary)

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
                     .loc[(~db[feature].eq('')), :]
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
              .merge(db
                     .loc[(~db[feature].eq('') & ~db[alteration_type].eq('')), :]
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
              .merge(db
                     .loc[(~db[feature].eq('') & ~db[alteration_type].eq('') & ~db[alteration].eq('')), :]
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

        # this is required for python 3.12 and pandas 2.2.2 to opt into future behavior for type downcasting
        with pd.option_context("future.no_silent_downcasting", True):
            db[cls.evidence_map_str] = (
                db[cls.evidence]
                .astype(str)
                .replace(cls.evidence_map)
                .astype(int)
            )
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
