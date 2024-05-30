import json


class Writer:
    section = "outputs"

    def __init__(self, strings):
        self.strings = strings

    """
    Defining properties for each datasource bin
    """

    @property
    def score_bin(self):
        return self.strings[self.section]['score_bin']

    @property
    def almanac_bin(self):
        return self.strings[self.section]['almanac_bin']

    @property
    def cancer_hotspots_bin(self):
        return self.strings[self.section]['cancer_hotspots_bin']

    @property
    def cancer_hotspots_3d_bin(self):
        return self.strings[self.section]['cancer_hotspots_3d_bin']

    @property
    def cgc_bin(self):
        return self.strings[self.section]['cgc_bin']

    @property
    def gsea_pathways_bin(self):
        return self.strings[self.section]['gsea_pathways_bin']

    @property
    def gsea_cm_bin(self):
        return self.strings[self.section]['gsea_cm_bin']

    @property
    def cosmic_bin(self):
        return self.strings[self.section]['cosmic_bin']

    @property
    def clinvar_bin(self):
        return self.strings[self.section]['clinvar_bin']

    @property
    def acmg_bin(self):
        return self.strings[self.section]['acmg_bin']

    @property
    def hereditary_bin(self):
        return self.strings[self.section]['hereditary_bin']

    @property
    def msi_bin(self):
        return self.strings[self.section]['msi_bin']

    """
    Defining properties for moalmanac specific annotations
    """

    @property
    def sensitive_bin(self):
        return self.strings[self.section]['sensitive_bin']

    @property
    def resistance_bin(self):
        return self.strings[self.section]['resistance_bin']

    @property
    def prognostic_bin(self):
        return self.strings[self.section]['prognostic_bin']

    @property
    def sensitive_implication(self):
        return self.strings[self.section]['sensitive_implication']

    @property
    def resistance_implication(self):
        return self.strings[self.section]['resistance_implication']

    @property
    def prognostic_implication(self):
        return self.strings[self.section]['prognostic_implication']

    @property
    def sensitive_map(self):
        return self.strings[self.section]['sensitive_implication_map']

    @property
    def resistance_map(self):
        return self.strings[self.section]['resistance_implication_map']

    @property
    def prognostic_map(self):
        return self.strings[self.section]['prognostic_implication_map']

    @property
    def sensitive_therapy(self):
        return self.strings[self.section]['sensitive_therapy']

    @property
    def resistance_therapy(self):
        return self.strings[self.section]['resistance_therapy']

    @property
    def sensitive_therapy_strategy(self):
        return self.strings[self.section]['sensitive_therapy_strategy']

    @property
    def resistance_therapy_strategy(self):
        return self.strings[self.section]['resistance_therapy_strategy']

    @property
    def sensitive_therapy_type(self):
        return self.strings[self.section]['sensitive_therapy_type']

    @property
    def resistance_therapy_type(self):
        return self.strings[self.section]['resistance_therapy_type']

    @property
    def favorable_prognosis(self):
        return self.strings[self.section]['favorable_prognosis']

    @property
    def sensitive_oncotree_code(self):
        return self.strings[self.section]['sensitive_oncotree_code']

    @property
    def resistance_oncotree_code(self):
        return self.strings[self.section]['resistance_oncotree_code']

    @property
    def prognostic_oncotree_code(self):
        return self.strings[self.section]['prognostic_oncotree_code']

    @property
    def sensitive_description(self):
        return self.strings[self.section]['sensitive_description']

    @property
    def resistance_description(self):
        return self.strings[self.section]['resistance_description']

    @property
    def prognostic_description(self):
        return self.strings[self.section]['prognostic_description']

    @property
    def sensitive_url(self):
        return self.strings[self.section]['sensitive_url']

    @property
    def resistance_url(self):
        return self.strings[self.section]['resistance_url']

    @property
    def prognostic_url(self):
        return self.strings[self.section]['prognostic_url']

    @property
    def sensitive_citation(self):
        return self.strings[self.section]['sensitive_citation']

    @property
    def resistance_citation(self):
        return self.strings[self.section]['resistance_citation']

    @property
    def prognostic_citation(self):
        return self.strings[self.section]['prognostic_citation']

    """
    Defining properties for describing biomarkers
    """

    @property
    def feature_type(self):
        return self.strings[self.section]['feature_type']

    @property
    def feature(self):
        return self.strings[self.section]['feature']

    @property
    def alt_type(self):
        return self.strings[self.section]['alt_type']

    @property
    def alt(self):
        return self.strings[self.section]['alt']

    @property
    def chr(self):
        return self.strings[self.section]['chr']

    @property
    def start(self):
        return self.strings[self.section]['start']

    @property
    def end(self):
        return self.strings[self.section]['end']

    @property
    def ref(self):
        return self.strings[self.section]['ref']

    @property
    def allele1(self):
        return self.strings[self.section]['allele1']

    @property
    def allele2(self):
        return self.strings[self.section]['allele2']

    @property
    def tumor_f(self):
        return self.strings[self.section]['tumor_f']

    @property
    def coverage(self):
        return self.strings[self.section]['coverage']

    @property
    def clinvar(self):
        return self.strings[self.section]['clinvar']

    @property
    def number_germline_mutations(self):
        return self.strings[self.section]['number_germline_mutations']

    """
    Defining properties for ExAC annotations
    """

    @property
    def exac_common(self):
        return self.strings[self.section]['exac_common']

    @property
    def exac_af(self):
        return self.strings[self.section]['exac_af']

    @property
    def exac_ac(self):
        return self.strings[self.section]['exac_ac']

    @property
    def exac_an(self):
        return self.strings[self.section]['exac_an']

    @property
    def exac_afr_ac(self):
        return self.strings[self.section]['exac_afr_ac']

    @property
    def exac_amr_ac(self):
        return self.strings[self.section]['exac_eas_ac']

    @property
    def exac_eas_ac(self):
        return self.strings[self.section]['exac_eas_ac']

    @property
    def exac_fin_ac(self):
        return self.strings[self.section]['exac_nfe_ac']

    @property
    def exac_nfe_ac(self):
        return self.strings[self.section]['exac_nfe_ac']

    @property
    def exac_sas_ac(self):
        return self.strings[self.section]['exac_sas_ac']

    @property
    def exac_oth_ac(self):
        return self.strings[self.section]['exac_oth_ac']

    @property
    def exac_afr_an(self):
        return self.strings[self.section]['exac_afr_an']

    @property
    def exac_amr_an(self):
        return self.strings[self.section]['exac_amr_an']

    @property
    def exac_eas_an(self):
        return self.strings[self.section]['exac_eas_an']

    @property
    def exac_fin_an(self):
        return self.strings[self.section]['exac_fin_an']

    @property
    def exac_nfe_an(self):
        return self.strings[self.section]['exac_nfe_an']

    @property
    def exac_sas_an(self):
        return self.strings[self.section]['exac_sas_an']

    @property
    def exac_oth_an(self):
        return self.strings[self.section]['exac_oth_an']

    """
    Defining properties to describe fusions
    """

    @property
    def spanningfrags(self):
        return self.strings[self.section]['spanningfrags']

    @property
    def left_gene(self):
        return self.strings[self.section]['left_gene']

    @property
    def left_chr(self):
        return self.strings[self.section]['left_chr']

    @property
    def left_start(self):
        return self.strings[self.section]['left_start']

    @property
    def right_gene(self):
        return self.strings[self.section]['right_gene']

    @property
    def right_chr(self):
        return self.strings[self.section]['right_chr']

    @property
    def right_start(self):
        return self.strings[self.section]['right_start']

    """
    Define properties for annotations to annotate somatic variants with those observed in validation sequencing 
    """

    @property
    def validation_tumor_f(self):
        return self.strings[self.section]['validation_tumor_f']

    @property
    def validation_coverage(self):
        return self.strings[self.section]['validation_coverage']

    @property
    def validation_detection_power(self):
        return self.strings[self.section]['validation_detection_power']

    """
    Defining remaining properties
    """

    @property
    def feature_display(self):
        return self.strings[self.section]['feature_display']

    @property
    def preclinical_efficacy(self):
        return self.strings[self.section]['preclinical_efficacy']

    @property
    def connections(self):
        return self.strings[self.section]['connections']

    @property
    def rationale(self):
        return self.strings[self.section]['rationale']

    @property
    def patient_id(self):
        return self.strings[self.section]['patient_id']

    @property
    def tumor(self):
        return self.strings[self.section]['tumor']

    @property
    def normal(self):
        return self.strings[self.section]['normal']

    """
    Writer class functions
    """

    @staticmethod
    def create_output_name(folder, patient_id, output_suffix):
        return f'{folder}/{patient_id}.{output_suffix}'

    @staticmethod
    def export_series(series, output_name):
        series.to_csv(output_name, sep='\t', index=True)

    @staticmethod
    def export_dataframe(df, output_name):
        df.to_csv(output_name, sep='\t', index=False)

    @staticmethod
    def export_dataframe_indexed(df, output_name, index_label):
        df.to_csv(output_name, sep='\t', index_label=index_label)

    @staticmethod
    def save_figure(figure, output_name):
        figure.savefig(output_name, bbox_inches='tight')

    @staticmethod
    def sort_columns(df, columns, ascending_boolean):
        return df.sort_values(columns, ascending=ascending_boolean)

    @staticmethod
    def return_nonzero_bin_idx(series):
        return series[series.astype(float) != 0].index


class Actionable:
    output_suffix = "actionable.txt"

    def __init__(self, strings):
        self.writer = Writer(strings)
        self.sort_columns = [
            self.writer.almanac_bin,
            self.writer.sensitive_map,
            self.writer.resistance_map,
            self.writer.prognostic_map
        ]
        self.output_columns = [
            self.writer.score_bin,
            self.writer.sensitive_implication,
            self.writer.resistance_implication,
            self.writer.prognostic_implication,
            self.writer.feature_type,
            self.writer.feature,
            self.writer.alt_type,
            self.writer.alt,
            self.writer.tumor_f,
            self.writer.coverage,
            self.writer.exac_af,
            self.writer.exac_common,
            self.writer.clinvar,
            self.writer.sensitive_bin,
            self.writer.sensitive_therapy,
            self.writer.sensitive_therapy_strategy,
            self.writer.sensitive_therapy_type,
            self.writer.sensitive_oncotree_code,
            self.writer.sensitive_description,
            self.writer.sensitive_citation,
            self.writer.sensitive_url,
            self.writer.resistance_bin,
            self.writer.resistance_therapy,
            self.writer.resistance_therapy_strategy,
            self.writer.resistance_therapy_type,
            self.writer.resistance_oncotree_code,
            self.writer.resistance_description,
            self.writer.resistance_citation,
            self.writer.resistance_url,
            self.writer.prognostic_bin,
            self.writer.favorable_prognosis,
            self.writer.prognostic_oncotree_code,
            self.writer.prognostic_description,
            self.writer.prognostic_citation,
            self.writer.prognostic_url,
            self.writer.number_germline_mutations,
            self.writer.validation_coverage,
            self.writer.validation_tumor_f,
            self.writer.validation_detection_power,
            self.writer.feature_display,
            self.writer.preclinical_efficacy,
            self.writer.patient_id,
            self.writer.tumor,
            self.writer.normal
        ]

    def write(self, dataframe, patient_label, folder):
        dataframe[self.writer.patient_id] = patient_label
        dataframe_sorted = Writer.sort_columns(
            df=dataframe,
            columns=self.sort_columns,
            ascending_boolean=False
        )
        output_name = Writer.create_output_name(
            folder=folder,
            patient_id=patient_label,
            output_suffix=self.__class__.output_suffix
        )
        Writer.export_dataframe(
            df=dataframe_sorted.loc[:, self.output_columns],
            output_name=output_name
        )
        return dataframe_sorted


class GermlineACMG:
    output_suffix = "germline.acmg.txt"

    def __init__(self, strings):
        self.writer = Writer(strings)
        self.bin = self.writer.acmg_bin
        self.sort_columns = [
            self.writer.clinvar_bin,
            self.writer.feature,
            self.writer.feature_type
        ]
        self.output_columns = [
            self.writer.feature_type,
            self.writer.feature,
            self.writer.alt_type,
            self.writer.alt,
            self.writer.chr,
            self.writer.start,
            self.writer.end,
            self.writer.ref,
            self.writer.allele1,
            self.writer.allele2,
            self.writer.tumor_f,
            self.writer.coverage,
            self.writer.clinvar,
            self.writer.exac_common,
            self.writer.exac_af,
            self.writer.exac_ac,
            self.writer.exac_an,
            self.writer.exac_afr_ac,
            self.writer.exac_amr_ac,
            self.writer.exac_eas_ac,
            self.writer.exac_fin_ac,
            self.writer.exac_nfe_ac,
            self.writer.exac_sas_ac,
            self.writer.exac_oth_ac,
            self.writer.exac_afr_an,
            self.writer.exac_amr_an,
            self.writer.exac_eas_an,
            self.writer.exac_fin_an,
            self.writer.exac_nfe_an,
            self.writer.exac_sas_ac,
            self.writer.exac_oth_an,
            self.writer.patient_id,
            self.writer.tumor,
            self.writer.normal
        ]

    def write(self, dataframe, patient_label, folder):
        dataframe[self.writer.patient_id] = patient_label
        dataframe_sorted = Writer.sort_columns(dataframe, self.sort_columns, ascending_boolean=False)
        idx = Writer.return_nonzero_bin_idx(dataframe.loc[:, self.bin])
        output_name = Writer.create_output_name(folder, patient_label, self.__class__.output_suffix)
        Writer.export_dataframe(dataframe_sorted.loc[idx, self.output_columns], output_name)
        return dataframe_sorted


class GermlineCancer:
    output_suffix = 'germline.cancer_related.txt'

    def __init__(self, strings):
        self.writer = Writer(strings)
        self.almanac_bin = self.writer.almanac_bin
        self.hotspots_bin = self.writer.cancer_hotspots_bin
        self.cgc_bin = self.writer.cgc_bin
        self.sort_columns = [
            self.writer.almanac_bin,
            self.writer.cancer_hotspots_bin,
            self.writer.cancer_hotspots_3d_bin,
            self.writer.cgc_bin,
            self.writer.gsea_pathways_bin,
            self.writer.gsea_cm_bin,
            self.writer.cosmic_bin,
            self.writer.exac_common,
            self.writer.exac_af
        ]
        self.ascending = [
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            True,
            True
        ]
        self.output_columns = [
            self.writer.score_bin,
            self.writer.sensitive_bin,
            self.writer.resistance_bin,
            self.writer.prognostic_bin,
            self.writer.feature_type,
            self.writer.feature,
            self.writer.alt_type,
            self.writer.alt,
            self.writer.chr,
            self.writer.start,
            self.writer.end,
            self.writer.ref,
            self.writer.allele1,
            self.writer.allele2,
            self.writer.tumor_f,
            self.writer.coverage,
            self.writer.clinvar,
            self.writer.exac_common,
            self.writer.exac_af,
            self.writer.exac_ac,
            self.writer.exac_an,
            self.writer.exac_afr_ac,
            self.writer.exac_amr_ac,
            self.writer.exac_eas_ac,
            self.writer.exac_fin_ac,
            self.writer.exac_nfe_ac,
            self.writer.exac_sas_ac,
            self.writer.exac_oth_ac,
            self.writer.exac_afr_an,
            self.writer.exac_amr_an,
            self.writer.exac_eas_an,
            self.writer.exac_fin_an,
            self.writer.exac_nfe_an,
            self.writer.exac_sas_ac,
            self.writer.exac_oth_an,
            self.writer.patient_id,
            self.writer.tumor,
            self.writer.normal
        ]

    def get_cancer_idx(self, dataframe):
        idx_almanac = Writer.return_nonzero_bin_idx(dataframe.loc[:, self.almanac_bin])
        idx_hotspot = Writer.return_nonzero_bin_idx(dataframe.loc[:, self.hotspots_bin])
        idx_cgc = Writer.return_nonzero_bin_idx(dataframe.loc[:, self.cgc_bin])
        return idx_almanac.union(idx_hotspot).union(idx_cgc)

    def write(self, dataframe, patient_label, folder):
        dataframe[self.writer.patient_id] = patient_label
        dataframe_sorted = Writer.sort_columns(
            df=dataframe,
            columns=self.sort_columns,
            ascending_boolean=self.ascending
        )
        idx = self.__class__.get_cancer_idx(self=self, dataframe=dataframe_sorted)
        output_name = Writer.create_output_name(
            folder=folder,
            patient_id=patient_label,
            output_suffix=self.__class__.output_suffix
        )
        Writer.export_dataframe(
            df=dataframe_sorted.loc[idx, self.output_columns],
            output_name=output_name
        )


class GermlineHereditary:
    output_suffix = 'germline.hereditary_cancers.txt'

    def __init__(self, strings):
        self.writer = Writer(strings)
        self.bin = self.writer.hereditary_bin
        self.sort_columns = [
            self.writer.clinvar_bin,
            self.writer.feature,
            self.writer.feature_type
        ]
        self.output_columns = [
            self.writer.feature_type,
            self.writer.feature,
            self.writer.alt_type,
            self.writer.alt,
            self.writer.chr,
            self.writer.start,
            self.writer.end,
            self.writer.ref,
            self.writer.allele1,
            self.writer.allele2,
            self.writer.tumor_f,
            self.writer.coverage,
            self.writer.clinvar,
            self.writer.exac_common,
            self.writer.exac_af,
            self.writer.exac_ac,
            self.writer.exac_an,
            self.writer.exac_afr_ac,
            self.writer.exac_amr_ac,
            self.writer.exac_eas_ac,
            self.writer.exac_fin_ac,
            self.writer.exac_nfe_ac,
            self.writer.exac_sas_ac,
            self.writer.exac_oth_ac,
            self.writer.exac_afr_an,
            self.writer.exac_amr_an,
            self.writer.exac_eas_an,
            self.writer.exac_fin_an,
            self.writer.exac_nfe_an,
            self.writer.exac_sas_ac,
            self.writer.exac_oth_an,
            self.writer.patient_id,
            self.writer.tumor,
            self.writer.normal
        ]

    def write(self, dataframe, patient_label, folder):
        dataframe[self.writer.patient_id] = patient_label
        dataframe_sorted = Writer.sort_columns(
            df=dataframe,
            columns=self.sort_columns,
            ascending_boolean=False
        )
        idx = Writer.return_nonzero_bin_idx(series=dataframe.loc[:, self.bin])
        output_name = Writer.create_output_name(
            folder=folder,
            patient_id=patient_label,
            output_suffix=self.__class__.output_suffix
        )
        Writer.export_dataframe(
            df=dataframe_sorted.loc[idx, self.output_columns],
            output_name=output_name
        )


class Illustrations:
    @staticmethod
    def write(figure, folder, profile_label, output_suffix):
        output_name = Writer.create_output_name(folder, profile_label, output_suffix)
        Writer.save_figure(figure=figure, output_name=output_name)


class Integrated:
    output_suffix = 'integrated.summary.txt'
    section = 'integrative'

    def __init__(self, strings):
        self.strings = strings
        self.writer = Writer(strings)
        self.output_columns = [self.somatic, self.copynumber, self.fusion, self.germline]

    @property
    def somatic(self):
        return self.strings[self.section]['somatic']

    @property
    def copynumber(self):
        return self.strings[self.section]['copynumber']

    @property
    def fusion(self):
        return self.strings[self.section]['fusion']

    @property
    def germline(self):
        return self.strings[self.section]['germline']

    def write(self, dataframe, patient_label, folder):
        dataframe_sorted = dataframe.sort_index()
        output_name = Writer.create_output_name(
            folder=folder,
            patient_id=patient_label,
            output_suffix=self.__class__.output_suffix
        )
        Writer.export_dataframe_indexed(
            df=dataframe_sorted.loc[:, self.output_columns],
            output_name=output_name,
            index_label=Writer.feature
        )


class MSI:
    output_suffix = 'msi_variants.txt'

    def __init__(self, strings):
        self.writer = Writer(strings)
        self.bin = self.writer.msi_bin
        self.sort_columns = [
            self.writer.almanac_bin,
            self.writer.cancer_hotspots_bin,
            self.writer.cancer_hotspots_3d_bin,
            self.writer.cgc_bin,
            self.writer.gsea_pathways_bin,
            self.writer.gsea_cm_bin,
            self.writer.cosmic_bin,
            self.writer.exac_common,
            self.writer.exac_af
        ],
        self.ascending=[
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            True,
            True
        ]
        self.output_columns = [
            self.writer.score_bin,
            self.writer.sensitive_bin,
            self.writer.resistance_bin,
            self.writer.prognostic_bin,
            self.writer.feature_type,
            self.writer.feature,
            self.writer.alt_type,
            self.writer.alt,
            self.writer.chr,
            self.writer.start,
            self.writer.end,
            self.writer.ref,
            self.writer.allele1,
            self.writer.allele2,
            self.writer.tumor_f,
            self.writer.coverage,
            self.writer.exac_af,
            self.writer.exac_common,
            self.writer.clinvar,
            self.writer.number_germline_mutations,
            self.writer.spanningfrags,
            self.writer.left_gene,
            self.writer.left_chr,
            self.writer.left_start,
            self.writer.right_gene,
            self.writer.right_chr,
            self.writer.right_start,
            self.writer.rationale,
            self.writer.patient_id,
            self.writer.tumor,
            self.writer.normal,
            self.writer.almanac_bin,
            self.writer.cancer_hotspots_bin,
            self.writer.cancer_hotspots_3d_bin,
            self.writer.cgc_bin,
            self.writer.gsea_pathways_bin,
            self.writer.gsea_cm_bin,
            self.writer.cosmic_bin
        ]

    def write(self, dataframe, patient_label, folder):
        dataframe[Writer.patient_id] = patient_label
        dataframe_sorted = Writer.sort_columns(
            df=dataframe,
            columns=self.sort_columns,
            ascending_boolean=False
        )
        output_name = Writer.create_output_name(
            folder=folder,
            patient_id=patient_label,
            output_suffix=self.__class__.output_suffix
        )
        Writer.export_dataframe(
            df=dataframe_sorted.loc[:, self.output_columns],
            output_name=output_name
        )


class MutationalBurden:
    output_suffix = 'mutational_burden.txt'
    section = 'burden'

    def __init__(self, strings):
        self.strings = strings
        self.writer = Writer(strings)
        self.output_columns = [
            self.patient,
            self.tumor_type,
            self.ontology,
            self.code,
            self.bases_covered,
            self.n_nonsyn_mutations,
            self.mutational_burden,
            self.percentile_tcga,
            self.percentile_tcga_tissue,
            self.high_burden_boolean
        ]

    @property
    def patient(self):
        return self.strings[self.section]['patient']

    @property
    def tumor_type(self):
        return self.strings[self.section]['tumor_type']

    @property
    def ontology(self):
        return self.strings[self.section]['ontology']

    @property
    def code(self):
        return self.strings[self.section]['code']

    @property
    def bases_covered(self):
        return self.strings[self.section]['bases_covered']

    @property
    def n_nonsyn_mutations(self):
        return self.strings[self.section]['n_nonsyn_mutations']

    @property
    def mutational_burden(self):
        return self.strings[self.section]['mutational_burden']

    @property
    def percentile_tcga(self):
        return self.strings[self.section]['percentile_tcga']

    @property
    def percentile_tcga_tissue(self):
        return self.strings[self.section]['percentile_tcga_tissue']

    @property
    def high_burden_boolean(self):
        return self.strings[self.section]['high_burden_boolean']

    def write(self, dataframe, patient_label, folder):
        dataframe[self.writer.patient_id] = patient_label
        output_name = Writer.create_output_name(
            folder=folder,
            patient_id=patient_label,
            output_suffix=self.__class__.output_suffix
        )
        Writer.export_dataframe(
            df=dataframe,
            output_name=output_name
        )


class PreclinicalEfficacy:
    output_suffix = 'preclinical.efficacy.txt'

    def __init__(self, strings):
        self.writer = Writer(strings)

    def write(self, dataframe, patient_label, folder):
        output_name = Writer.create_output_name(
            folder=folder,
            patient_id=patient_label,
            output_suffix=self.__class__.output_suffix
        )
        Writer.export_dataframe(
            df=dataframe,
            output_name=output_name
        )


class PreclinicalMatchmaking:
    output_suffix = 'matchmaker.txt'

    def __init__(self, strings):
        self.writer = Writer(strings)

    def write(self, dataframe, patient_label, folder):
        output_name = Writer.create_output_name(
            folder=folder,
            patient_id=patient_label,
            output_suffix=self.__class__.output_suffix
        )
        Writer.export_dataframe(
            df=dataframe,
            output_name=output_name
        )


class SomaticFiltered:
    output_suffix = 'somatic.filtered.txt'

    def __init__(self, strings):
        self.writer = Writer(strings)
        self.sort_columns = [
            self.writer.feature,
            self.writer.feature_type
        ]
        self.output_columns = [
            self.writer.feature_type,
            self.writer.feature,
            self.writer.alt_type,
            self.writer.alt,
            self.writer.chr,
            self.writer.start,
            self.writer.end,
            self.writer.ref,
            self.writer.allele1,
            self.writer.allele2,
            self.writer.tumor_f,
            self.writer.coverage,
            self.writer.spanningfrags,
            self.writer.left_gene,
            self.writer.left_chr,
            self.writer.left_start,
            self.writer.right_gene,
            self.writer.right_chr,
            self.writer.right_start,
            self.writer.rationale,
            self.writer.patient_id,
            self.writer.tumor,
            self.writer.normal
        ]

    def write(self, dataframe, patient_label, folder):
        dataframe[Writer.patient_id] = patient_label
        dataframe_sorted = Writer.sort_columns(
            df=dataframe,
            columns=self.sort_columns,
            ascending_boolean=False
        )
        output_name = Writer.create_output_name(
            folder=folder,
            patient_id=patient_label,
            output_suffix=self.__class__.output_suffix
        )
        Writer.export_dataframe(
            df=dataframe_sorted.loc[:, self.output_columns],
            output_name=output_name
        )


class SomaticScored:
    output_suffix = 'somatic.scored.txt'

    def __init__(self, strings):
        self.writer = Writer(strings)
        self.sort_columns = [
            self.writer.almanac_bin,
            self.writer.cancer_hotspots_bin,
            self.writer.cancer_hotspots_3d_bin,
            self.writer.cgc_bin,
            self.writer.gsea_pathways_bin,
            self.writer.gsea_cm_bin,
            self.writer.cosmic_bin,
            self.writer.validation_detection_power,
            self.writer.validation_coverage,
            self.writer.number_germline_mutations,
            self.writer.exac_common,
            self.writer.exac_af
        ]
        self.ascending = [
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            True,
            True
        ]
        self.output_columns = [
            self.writer.score_bin,
            self.writer.sensitive_bin,
            self.writer.resistance_bin,
            self.writer.prognostic_bin,
            self.writer.feature_type,
            self.writer.feature,
            self.writer.alt_type,
            self.writer.alt,
            self.writer.chr,
            self.writer.start,
            self.writer.end,
            self.writer.ref,
            self.writer.allele1,
            self.writer.allele2,
            self.writer.tumor_f,
            self.writer.coverage,
            self.writer.exac_af,
            self.writer.exac_common,
            self.writer.clinvar,
            self.writer.number_germline_mutations,
            self.writer.spanningfrags,
            self.writer.left_gene,
            self.writer.left_chr,
            self.writer.left_start,
            self.writer.right_gene,
            self.writer.right_chr,
            self.writer.right_start,
            self.writer.validation_coverage,
            self.writer.validation_tumor_f,
            self.writer.validation_detection_power,
            self.writer.rationale,
            self.writer.patient_id,
            self.writer.tumor,
            self.writer.normal,
            self.writer.almanac_bin,
            self.writer.cancer_hotspots_bin,
            self.writer.cancer_hotspots_3d_bin,
            self.writer.cgc_bin,
            self.writer.gsea_pathways_bin,
            self.writer.gsea_cm_bin,
            self.writer.cosmic_bin
        ]

    def write(self, dataframe, patient_label, folder):
        dataframe[self.writer.patient_id] = patient_label
        dataframe_sorted = Writer.sort_columns(
            df=dataframe,
            columns=self.sort_columns,
            ascending_boolean=self.ascending
        )
        output_name = Writer.create_output_name(
            folder=folder,
            patient_id=patient_label,
            output_suffix=self.__class__.output_suffix)
        Writer.export_dataframe(
            df=dataframe_sorted.loc[:, self.output_columns],
            output_name=output_name
        )


class Strategies:
    output_suffix = 'therapeutic_strategies.txt'

    def __init__(self, strings):
        self.writer = Writer(strings)

    def write(self, dataframe, patient_label, folder):
        output_name = Writer.create_output_name(
            folder=folder,
            patient_id=patient_label,
            output_suffix=self.__class__.output_suffix
        )
        Writer.export_dataframe_indexed(
            df=dataframe,
            output_name=output_name,
            index_label='Assertion / Strategy'
        )


class Json:
    @staticmethod
    def write(handle, dictionary):
        with open(handle, 'w') as json_handle:
            json.dump(dictionary, json_handle, sort_keys=True, indent=4)
